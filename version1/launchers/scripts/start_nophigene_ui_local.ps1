[CmdletBinding()]
param(
    [int]$Port = 8766,
    [switch]$UseDocker,
    [string]$ImageName = "nophigene:latest",
    [string]$ContainerName = "nophigene-ui",
    [switch]$SkipBuild,
    [switch]$EnableLocalExtraction,
    [switch]$NoOpenBrowser,
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = (Resolve-Path (Join-Path $scriptDir "..")).Path
$dockerStartScript = Join-Path $scriptDir "start_nophigene_ui.ps1"
$pythonPath = Join-Path $repoRoot ".venv\Scripts\python.exe"
$dataDir = Join-Path $repoRoot "data"
$referenceDir = Join-Path $dataDir "reference\hg38"
$extractedDir = Join-Path $dataDir "extracted"
$resultsDir = Join-Path $repoRoot "results"
$pidFile = Join-Path $repoRoot ".nophigene-ui.pid"
$portFile = Join-Path $repoRoot ".nophigene-ui.port"
$stdoutLog = Join-Path $repoRoot ".nophigene-ui.log"
$stderrLog = Join-Path $repoRoot ".nophigene-ui.err.log"
$portWasExplicit = $PSBoundParameters.ContainsKey("Port")

if ($UseDocker) {
    Write-Host ""
    Write-Host "NophiGene local launcher is delegating to Docker for Extraction support." -ForegroundColor Green
    Write-Host "Docker image will include samtools and bcftools when rebuilt from the current Dockerfile."
    Write-Host ""

    $dockerParams = @{
        ImageName = $ImageName
        ContainerName = $ContainerName
    }
    if ($portWasExplicit) { $dockerParams.Port = $Port }
    if ($SkipBuild) { $dockerParams.SkipBuild = $true }
    if ($NoOpenBrowser) { $dockerParams.NoOpenBrowser = $true }
    if ($DryRun) { $dockerParams.DryRun = $true }

    & $dockerStartScript @dockerParams
    exit $LASTEXITCODE
}

function Write-Step {
    param([string]$Message)
    Write-Host "==> $Message" -ForegroundColor Cyan
}

function Ensure-Directories {
    foreach ($dir in @($dataDir, $referenceDir, $extractedDir, $resultsDir)) {
        if (-not (Test-Path $dir)) {
            Write-Step "Creating $dir"
            if (-not $DryRun) {
                New-Item -ItemType Directory -Path $dir -Force | Out-Null
            }
        }
    }
}

function Test-CommandAvailable {
    param([string]$Name)
    return $null -ne (Get-Command $Name -ErrorAction SilentlyContinue)
}

function Configure-LocalExtraction {
    $envOverrideEnabled = $env:NOPHIGENE_ENABLE_LOCAL_EXTRACTION -eq "1"
    if (-not $EnableLocalExtraction -and -not $envOverrideEnabled) {
        Write-Host "Extraction  : disabled in local mode; use Docker for BAM extraction." -ForegroundColor Yellow
        return
    }

    $missingTools = @()
    foreach ($toolName in @("samtools", "bcftools")) {
        if (-not (Test-CommandAvailable -Name $toolName)) {
            $missingTools += $toolName
        }
    }

    if ($missingTools.Count -gt 0) {
        throw "Local extraction was requested, but these tools were not found on PATH: $($missingTools -join ', '). Use the Docker launcher or install the tools locally."
    }

    $env:NOPHIGENE_ENABLE_LOCAL_EXTRACTION = "1"
    Write-Host "Extraction  : enabled for this local process via NOPHIGENE_ENABLE_LOCAL_EXTRACTION=1." -ForegroundColor Green
}

function Test-WebReady {
    param([int]$WebPort)

    if ($DryRun) {
        return $true
    }

    try {
        $response = Invoke-WebRequest -Uri "http://127.0.0.1:$WebPort/api/v1/health" -UseBasicParsing -TimeoutSec 3
        if ($response.StatusCode -ne 200) {
            return $false
        }
        $payload = $response.Content | ConvertFrom-Json
        return $payload.status -in @("ok", "degraded")
    }
    catch {
        return $false
    }
}

function Test-PortAvailable {
    param([int]$WebPort)

    if ($DryRun) {
        return $true
    }

    $listener = $null
    try {
        $listener = New-Object System.Net.Sockets.TcpListener([System.Net.IPAddress]::Loopback, $WebPort)
        $listener.Start()
        return $true
    }
    catch {
        return $false
    }
    finally {
        if ($null -ne $listener) {
            $listener.Stop()
        }
    }
}

function Resolve-WebPort {
    param(
        [int]$RequestedPort,
        [bool]$Explicit
    )

    if ($RequestedPort -lt 1 -or $RequestedPort -gt 65535) {
        throw "Port must be between 1 and 65535."
    }
    if (Test-PortAvailable -WebPort $RequestedPort) {
        return $RequestedPort
    }
    if ($Explicit) {
        throw "Port $RequestedPort is already in use or cannot be bound. Choose another port with -Port, for example -Port 9000."
    }

    $lastCandidate = [Math]::Min($RequestedPort + 50, 65535)
    if ($RequestedPort -lt $lastCandidate) {
        foreach ($candidate in (($RequestedPort + 1)..$lastCandidate)) {
            if (Test-PortAvailable -WebPort $candidate) {
                Write-Host "Port $RequestedPort is unavailable; using port $candidate instead." -ForegroundColor Yellow
                return $candidate
            }
        }
    }
    throw "No available local port was found between $RequestedPort and $lastCandidate."
}

function Wait-ForWebApp {
    param(
        [int]$WebPort,
        [int]$ProcessId
    )

    if ($DryRun) {
        Write-Step "Would wait for http://127.0.0.1:$WebPort to respond"
        return
    }

    Write-Step "Waiting for the local web UI to respond on port $WebPort"
    $deadline = (Get-Date).AddMinutes(2)
    while ((Get-Date) -lt $deadline) {
        if (Test-WebReady -WebPort $WebPort) {
            return
        }
        if ($ProcessId -gt 0 -and -not (Get-Process -Id $ProcessId -ErrorAction SilentlyContinue)) {
            $errorTail = ""
            if (Test-Path $stderrLog) {
                $errorTail = (Get-Content $stderrLog -Tail 20 -ErrorAction SilentlyContinue) -join [Environment]::NewLine
            }
            Remove-Item -LiteralPath $pidFile -Force -ErrorAction SilentlyContinue
            Remove-Item -LiteralPath $portFile -Force -ErrorAction SilentlyContinue
            throw "The local UI process exited before becoming ready on port $WebPort.`n$errorTail"
        }
        Start-Sleep -Seconds 2
    }

    Remove-Item -LiteralPath $pidFile -Force -ErrorAction SilentlyContinue
    Remove-Item -LiteralPath $portFile -Force -ErrorAction SilentlyContinue
    throw "The local UI did not respond on http://127.0.0.1:$WebPort within 2 minutes. Check $stdoutLog and $stderrLog."
}

function Get-TrackedProcess {
    if (-not (Test-Path $pidFile)) {
        return $null
    }

    $rawPid = Get-Content $pidFile -ErrorAction SilentlyContinue | Select-Object -First 1
    if (-not $rawPid) {
        return $null
    }

    try {
        return Get-Process -Id ([int]$rawPid) -ErrorAction Stop
    }
    catch {
        return $null
    }
}

function Get-ListeningProcessId {
    param([int]$WebPort)

    if ($DryRun) {
        return $null
    }

    $pattern = "^\s*TCP\s+\S+:$WebPort\s+\S+\s+\S+\s+(?<ProcessId>\d+)\s*$"
    foreach ($line in (netstat -ano -p tcp)) {
        if ($line -match $pattern) {
            return [int]$Matches.ProcessId
        }
    }
    return $null
}

function Get-TrackedPort {
    if (-not (Test-Path $portFile)) {
        return $null
    }
    $rawPort = Get-Content $portFile -ErrorAction SilentlyContinue | Select-Object -First 1
    if (-not $rawPort) {
        return $null
    }
    try {
        return [int]$rawPort
    }
    catch {
        return $null
    }
}

Write-Host ""
Write-Host "NophiGene local launcher" -ForegroundColor Green
Write-Host "Repository : $repoRoot"
Write-Host "Python     : $pythonPath"
Write-Host "Requested  : $Port"
Write-Host "Reference  : $referenceDir"
Write-Host "Extracted  : $extractedDir"
Write-Host ""

if (-not (Test-Path $pythonPath)) {
    throw "No local Python environment was found at $pythonPath. Recreate .venv first."
}

Ensure-Directories
Configure-LocalExtraction

$existing = Get-TrackedProcess
if ($existing) {
    $trackedPort = Get-TrackedPort
    if ($trackedPort -and (Test-WebReady -WebPort $trackedPort)) {
        Write-Step "A tracked UI process is already running with PID $($existing.Id)"
        if (-not $NoOpenBrowser) {
            $url = "http://127.0.0.1:$trackedPort"
            Write-Step "Opening $url"
            if (-not $DryRun) {
                Start-Process $url | Out-Null
            }
        }
        exit 0
    }

    Write-Step "A stale tracked UI process was found. Restarting it."
    if (-not $DryRun) {
        Stop-Process -Id $existing.Id -Force -ErrorAction SilentlyContinue
    }
}
elseif (Test-Path $pidFile) {
    Remove-Item -LiteralPath $pidFile -Force
}
Remove-Item -LiteralPath $portFile -Force -ErrorAction SilentlyContinue

$Port = Resolve-WebPort -RequestedPort $Port -Explicit $portWasExplicit
Write-Host "Selected   : $Port"

Write-Step "Checking local app dependencies"
if (-not $DryRun) {
    & $pythonPath -c "import flask, allel, pandas, methylprep"
    if ($LASTEXITCODE -ne 0) {
        throw "The local environment is missing app dependencies. Install requirements-app.txt into .venv first."
    }
}

Write-Step "Starting the local web UI"
if (-not $DryRun) {
    if (Test-Path $stdoutLog) { Remove-Item -LiteralPath $stdoutLog -Force }
    if (Test-Path $stderrLog) { Remove-Item -LiteralPath $stderrLog -Force }

    $process = Start-Process `
        -FilePath $pythonPath `
        -ArgumentList @("src/app.py", "web", "--host", "127.0.0.1", "--port", "$Port") `
        -WorkingDirectory $repoRoot `
        -RedirectStandardOutput $stdoutLog `
        -RedirectStandardError $stderrLog `
        -PassThru

    Set-Content -LiteralPath $pidFile -Value $process.Id -NoNewline
    Set-Content -LiteralPath $portFile -Value $Port -NoNewline
}

$startedProcessId = if ($DryRun) { 0 } else { $process.Id }
Wait-ForWebApp -WebPort $Port -ProcessId $startedProcessId
if (-not $DryRun) {
    $processIds = @($process.Id)
    $listenerProcessId = Get-ListeningProcessId -WebPort $Port
    if ($listenerProcessId -and $listenerProcessId -notin $processIds) {
        $processIds += $listenerProcessId
    }
    Set-Content -LiteralPath $pidFile -Value $processIds
}

if (-not $NoOpenBrowser) {
    $url = "http://127.0.0.1:$Port"
    Write-Step "Opening $url"
    if (-not $DryRun) {
        Start-Process $url | Out-Null
    }
}

Write-Host ""
Write-Host "Local UI is ready." -ForegroundColor Green
Write-Host "Open http://127.0.0.1:$Port if the browser did not appear automatically."
Write-Host "BAM extraction remains Docker-first unless local samtools/bcftools are explicitly enabled."
Write-Host "Use 'Stop NophiGene UI.cmd' to stop the local server later."
