[CmdletBinding()]
param(
    [int]$Port = 8766,
    [string]$ImageName = "nophigene:latest",
    [string]$ContainerName = "nophigene-ui",
    [switch]$SkipBuild,
    [switch]$NoOpenBrowser,
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = (Resolve-Path (Join-Path $scriptDir "..")).Path
$dockerConfigDir = Join-Path $repoRoot ".docker-local"
$dataDir = Join-Path $repoRoot "data"
$referenceDir = Join-Path $dataDir "reference\hg38"
$extractedDir = Join-Path $dataDir "extracted"
$resultsDir = Join-Path $repoRoot "results"
$portWasExplicit = $PSBoundParameters.ContainsKey("Port")
$dockerDesktopCandidates = @(
    "C:\Program Files\Docker\Docker\Docker Desktop.exe",
    (Join-Path $env:LocalAppData "Programs\Docker\Docker\Docker Desktop.exe")
)

function Write-Step {
    param([string]$Message)
    Write-Host "==> $Message" -ForegroundColor Cyan
}

function Invoke-Checked {
    param(
        [scriptblock]$Command,
        [string]$Description
    )

    Write-Step $Description
    if ($DryRun) {
        return
    }

    & $Command
    if ($LASTEXITCODE -ne 0) {
        throw "Failed: $Description"
    }
}

function Test-DockerReady {
    if ($DryRun) {
        return $true
    }

    try {
        docker info *> $null
        return $true
    }
    catch {
        return $false
    }
}

function Ensure-DockerDesktop {
    if (Test-DockerReady) {
        return
    }

    $dockerDesktopPath = $dockerDesktopCandidates | Where-Object { Test-Path $_ } | Select-Object -First 1
    if (-not $dockerDesktopPath) {
        throw "Docker Desktop was not found. Install Docker Desktop or update scripts/start_nophigene_ui.ps1 with the correct path."
    }

    Write-Step "Starting Docker Desktop"
    if (-not $DryRun) {
        Start-Process -FilePath $dockerDesktopPath | Out-Null
        $deadline = (Get-Date).AddMinutes(3)
        while ((Get-Date) -lt $deadline) {
            if (Test-DockerReady) {
                return
            }
            Start-Sleep -Seconds 3
        }
        throw "Docker Desktop started, but the Docker engine did not become ready within 3 minutes."
    }
}

function Ensure-Directories {
    foreach ($dir in @($dockerConfigDir, $dataDir, $referenceDir, $extractedDir, $resultsDir)) {
        if (-not (Test-Path $dir)) {
            Write-Step "Creating $dir"
            if (-not $DryRun) {
                New-Item -ItemType Directory -Path $dir -Force | Out-Null
            }
        }
    }
}

function Test-DockerImageHasExtractionTools {
    if ($DryRun) {
        Write-Step "Would verify samtools and bcftools inside $ImageName"
        return $true
    }

    docker image inspect $ImageName *> $null
    if ($LASTEXITCODE -ne 0) {
        return $false
    }

    docker run --rm --entrypoint sh $ImageName -c "command -v samtools >/dev/null 2>&1 && command -v bcftools >/dev/null 2>&1" *> $null
    return $LASTEXITCODE -eq 0
}

function Remove-ExistingContainer {
    if ($DryRun) {
        Write-Step "Would remove existing container named $ContainerName if present"
        return
    }

    $existing = docker ps -a --filter "name=^/${ContainerName}$" --format "{{.Names}}"
    if ($LASTEXITCODE -ne 0) {
        throw "Failed while checking for an existing container."
    }

    if ($existing) {
        Invoke-Checked -Description "Removing existing container $ContainerName" -Command {
            docker rm -f $ContainerName | Out-Null
        }
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
        [string]$RunningContainer
    )

    if ($DryRun) {
        Write-Step "Would wait for http://127.0.0.1:$WebPort to respond"
        return
    }

    Write-Step "Waiting for the web UI to respond on port $WebPort"
    $deadline = (Get-Date).AddMinutes(2)
    while ((Get-Date) -lt $deadline) {
        try {
            $response = Invoke-WebRequest -Uri "http://127.0.0.1:$WebPort/api/v1/health" -UseBasicParsing -TimeoutSec 5
            $payload = $response.Content | ConvertFrom-Json
            if ($response.StatusCode -eq 200 -and $payload.status -in @("ok", "degraded")) {
                return
            }
        }
        catch {
            $running = docker inspect -f "{{.State.Running}}" $RunningContainer 2>$null
            if ($LASTEXITCODE -ne 0 -or $running -ne "true") {
                $containerLogs = docker logs --tail 40 $RunningContainer 2>&1
                throw "The container exited before the UI became ready on port $WebPort.`n$($containerLogs -join [Environment]::NewLine)"
            }
        }
        Start-Sleep -Seconds 2
    }

    throw "The container started, but the web UI did not respond on http://127.0.0.1:$WebPort within 2 minutes."
}

$env:DOCKER_CONFIG = $dockerConfigDir

Write-Host ""
Write-Host "NophiGene launcher" -ForegroundColor Green
Write-Host "Repository : $repoRoot"
Write-Host "Image      : $ImageName"
Write-Host "Container  : $ContainerName"
Write-Host "Requested  : $Port"
Write-Host "Reference  : $referenceDir"
Write-Host "Extracted  : $extractedDir"
Write-Host ""

Ensure-Directories
Ensure-DockerDesktop

if (-not $SkipBuild) {
    Invoke-Checked -Description "Building Docker image $ImageName" -Command {
        docker build -t $ImageName $repoRoot
    }
}
else {
    Write-Step "Skipping image build"
}

if (-not (Test-DockerImageHasExtractionTools)) {
    throw "Docker image $ImageName does not include samtools and bcftools. Re-run this launcher without -SkipBuild so the updated Extraction-capable image is built."
}

Remove-ExistingContainer
$Port = Resolve-WebPort -RequestedPort $Port -Explicit $portWasExplicit
Write-Host "Selected   : $Port"

Invoke-Checked -Description "Starting container $ContainerName" -Command {
    docker run -d `
        --name $ContainerName `
        -p "${Port}:8766" `
        -e "NOPHIGENE_IN_DOCKER=1" `
        -v "${dataDir}:/home/appuser/app/data" `
        -v "${resultsDir}:/home/appuser/app/results" `
        $ImageName | Out-Null
}

Wait-ForWebApp -WebPort $Port -RunningContainer $ContainerName

if (-not $NoOpenBrowser) {
    $url = "http://127.0.0.1:$Port"
    Write-Step "Opening $url"
    if (-not $DryRun) {
        Start-Process $url | Out-Null
    }
}

Write-Host ""
Write-Host "UI is ready." -ForegroundColor Green
Write-Host "Open http://127.0.0.1:$Port if the browser did not appear automatically."
Write-Host "The Extraction tab is enabled in Docker and writes reference/VCF files under data\reference\hg38 and data\extracted."
Write-Host "Use 'Stop NophiGene UI (Docker).cmd' to stop the container later."
