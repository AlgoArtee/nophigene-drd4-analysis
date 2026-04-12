[CmdletBinding()]
param(
    [int]$Port = 8000,
    [switch]$NoOpenBrowser,
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = (Resolve-Path (Join-Path $scriptDir "..")).Path
$pythonPath = Join-Path $repoRoot ".venv\Scripts\python.exe"
$dataDir = Join-Path $repoRoot "data"
$resultsDir = Join-Path $repoRoot "results"
$pidFile = Join-Path $repoRoot ".nophigene-ui.pid"
$stdoutLog = Join-Path $repoRoot ".nophigene-ui.log"
$stderrLog = Join-Path $repoRoot ".nophigene-ui.err.log"

function Write-Step {
    param([string]$Message)
    Write-Host "==> $Message" -ForegroundColor Cyan
}

function Ensure-Directories {
    foreach ($dir in @($dataDir, $resultsDir)) {
        if (-not (Test-Path $dir)) {
            Write-Step "Creating $dir"
            if (-not $DryRun) {
                New-Item -ItemType Directory -Path $dir -Force | Out-Null
            }
        }
    }
}

function Test-WebReady {
    param([int]$WebPort)

    if ($DryRun) {
        return $true
    }

    try {
        $response = Invoke-WebRequest -Uri "http://127.0.0.1:$WebPort" -UseBasicParsing -TimeoutSec 3
        return $response.StatusCode -ge 200 -and $response.StatusCode -lt 500
    }
    catch {
        return $false
    }
}

function Wait-ForWebApp {
    param([int]$WebPort)

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
        Start-Sleep -Seconds 2
    }

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

Write-Host ""
Write-Host "NophiGene local launcher" -ForegroundColor Green
Write-Host "Repository : $repoRoot"
Write-Host "Python     : $pythonPath"
Write-Host "Port       : $Port"
Write-Host ""

if (-not (Test-Path $pythonPath)) {
    throw "No local Python environment was found at $pythonPath. Recreate .venv first."
}

Ensure-Directories

$existing = Get-TrackedProcess
if ($existing) {
    if (Test-WebReady -WebPort $Port) {
        Write-Step "A tracked UI process is already running with PID $($existing.Id)"
        if (-not $NoOpenBrowser) {
            $url = "http://127.0.0.1:$Port"
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
}

Wait-ForWebApp -WebPort $Port

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
Write-Host "Use 'Stop NophiGene UI.cmd' to stop the local server later."
