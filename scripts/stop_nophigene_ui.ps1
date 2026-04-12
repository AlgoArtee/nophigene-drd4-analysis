[CmdletBinding()]
param(
    [string]$ContainerName = "nophigene-drd4-analysis-ui"
)

$ErrorActionPreference = "Stop"
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = (Resolve-Path (Join-Path $scriptDir "..")).Path
$dockerConfigDir = Join-Path $repoRoot ".docker-local"
$env:DOCKER_CONFIG = $dockerConfigDir

function Test-DockerReady {
    try {
        docker info *> $null
        return $true
    }
    catch {
        return $false
    }
}

Write-Host ""
Write-Host "Stopping NophiGene UI" -ForegroundColor Cyan
Write-Host "Container: $ContainerName"
Write-Host ""

if (-not (Test-DockerReady)) {
    Write-Host "Docker is not running. If no container is active, there is nothing to stop." -ForegroundColor Yellow
    exit 0
}

$existing = docker ps -a --filter "name=^/${ContainerName}$" --format "{{.Names}}"
if ($LASTEXITCODE -ne 0) {
    throw "Failed while checking for an existing container."
}

if (-not $existing) {
    Write-Host "No matching container was found. Nothing to stop." -ForegroundColor Yellow
    exit 0
}

docker rm -f $ContainerName | Out-Null
if ($LASTEXITCODE -ne 0) {
    throw "Failed to stop and remove container $ContainerName."
}

Write-Host "Container stopped and removed." -ForegroundColor Green
