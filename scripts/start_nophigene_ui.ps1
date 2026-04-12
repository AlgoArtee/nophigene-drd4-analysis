[CmdletBinding()]
param(
    [int]$Port = 8000,
    [string]$ImageName = "nophigene-drd4-analysis:latest",
    [string]$ContainerName = "nophigene-drd4-analysis-ui",
    [switch]$SkipBuild,
    [switch]$NoOpenBrowser,
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = (Resolve-Path (Join-Path $scriptDir "..")).Path
$dockerConfigDir = Join-Path $repoRoot ".docker-local"
$dataDir = Join-Path $repoRoot "data"
$resultsDir = Join-Path $repoRoot "results"
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
    foreach ($dir in @($dockerConfigDir, $dataDir, $resultsDir)) {
        if (-not (Test-Path $dir)) {
            Write-Step "Creating $dir"
            if (-not $DryRun) {
                New-Item -ItemType Directory -Path $dir -Force | Out-Null
            }
        }
    }
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

function Wait-ForWebApp {
    param([int]$WebPort)

    if ($DryRun) {
        Write-Step "Would wait for http://127.0.0.1:$WebPort to respond"
        return
    }

    Write-Step "Waiting for the web UI to respond on port $WebPort"
    $deadline = (Get-Date).AddMinutes(2)
    while ((Get-Date) -lt $deadline) {
        try {
            $response = Invoke-WebRequest -Uri "http://127.0.0.1:$WebPort" -UseBasicParsing -TimeoutSec 5
            if ($response.StatusCode -ge 200 -and $response.StatusCode -lt 500) {
                return
            }
        }
        catch {
            Start-Sleep -Seconds 2
        }
    }

    throw "The container started, but the web UI did not respond on http://127.0.0.1:$WebPort within 2 minutes."
}

$env:DOCKER_CONFIG = $dockerConfigDir

Write-Host ""
Write-Host "NophiGene launcher" -ForegroundColor Green
Write-Host "Repository : $repoRoot"
Write-Host "Image      : $ImageName"
Write-Host "Container  : $ContainerName"
Write-Host "Port       : $Port"
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

Remove-ExistingContainer

Invoke-Checked -Description "Starting container $ContainerName" -Command {
    docker run -d `
        --name $ContainerName `
        -p "${Port}:8000" `
        -v "${dataDir}:/home/appuser/app/data" `
        -v "${resultsDir}:/home/appuser/app/results" `
        $ImageName | Out-Null
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
Write-Host "UI is ready." -ForegroundColor Green
Write-Host "Open http://127.0.0.1:$Port if the browser did not appear automatically."
Write-Host "Use 'Stop NophiGene UI.cmd' to stop the container later."
