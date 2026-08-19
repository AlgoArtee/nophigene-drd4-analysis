[CmdletBinding()]
param(
    [ValidateRange(1, 300)]
    [int]$ShutdownTimeoutSeconds = 30,
    [switch]$KeepRuntimeSecrets,
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"
$startedAt = Get-Date
$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$composePath = Join-Path $projectRoot "docker-compose.yml"
$dockerConfigDir = Join-Path $projectRoot ".docker-local"
$secretRoot = Join-Path $projectRoot ".runtime-secrets"
$databaseSecretPath = Join-Path $secretRoot "database-key"
$sessionSecretPath = Join-Path $secretRoot "session-token"
$dandelionSecretPath = Join-Path $secretRoot "dandelion-runner-key"
$dandelionArtifactSecretPath = Join-Path $secretRoot "dandelion-artifact-key"
$modelRunnerSecretPath = Join-Path $secretRoot "model-runner-key"
$modelCredentialSecretPath = Join-Path $secretRoot "model-credential-key"
$stageCount = 5

function Write-Stage {
    param([int]$Number, [string]$Message)
    $timestamp = Get-Date -Format "HH:mm:ss"
    Write-Host "[$timestamp] [$Number/$stageCount] $Message" -ForegroundColor Cyan
}

function Write-Detail {
    param([string]$Label, [string]$Value)
    Write-Host ("  {0,-20} {1}" -f ($Label + ":"), $Value)
}

function Write-Success {
    param([string]$Message)
    Write-Host "  OK  $Message" -ForegroundColor Green
}

function Write-Notice {
    param([string]$Message)
    Write-Host "  --  $Message" -ForegroundColor DarkGray
}

function Remove-SecretFile {
    param([string]$Path, [string]$Label)
    if (-not (Test-Path -LiteralPath $Path)) {
        Write-Notice "$Label was already absent."
        return
    }
    if ($DryRun) {
        Write-Notice "Dry run: would remove $Label at $Path."
        return
    }
    Remove-Item -LiteralPath $Path -Force
    Write-Success "$Label removed."
}

Write-Host ""
Write-Host "============================================================" -ForegroundColor Green
Write-Host " NophiGene Version 2 - verbose shutdown" -ForegroundColor Green
Write-Host "============================================================" -ForegroundColor Green
Write-Detail "Started" $startedAt.ToString("yyyy-MM-dd HH:mm:ss K")
Write-Detail "Project" $projectRoot
Write-Detail "Compose file" $composePath
Write-Detail "Docker config" $dockerConfigDir
Write-Detail "Timeout" "$ShutdownTimeoutSeconds seconds"
Write-Detail "Keep secrets" $KeepRuntimeSecrets.IsPresent
Write-Detail "Dry run" $DryRun.IsPresent
Write-Host ""

Write-Stage 1 "Validating the shutdown target"
if (-not (Test-Path -LiteralPath $composePath -PathType Leaf)) {
    throw "Compose file not found: $composePath"
}
$env:DOCKER_CONFIG = $dockerConfigDir
Write-Success "Resolved the Version 2 Compose project."

Write-Stage 2 "Checking Docker availability and current service state"
$dockerAvailable = $null -ne (Get-Command docker -ErrorAction SilentlyContinue)
$dockerReady = $false
if ($dockerAvailable -and -not $DryRun) {
    & docker info *> $null
    $dockerReady = $LASTEXITCODE -eq 0
}
elseif ($DryRun) {
    $dockerReady = $dockerAvailable
}

if (-not $dockerAvailable) {
    Write-Host "  WARN Docker CLI is not installed or is not on PATH; container state cannot be inspected." -ForegroundColor Yellow
}
elseif (-not $dockerReady) {
    Write-Host "  WARN Docker engine is not running; no running container can be stopped right now." -ForegroundColor Yellow
}
else {
    Write-Notice "Current service state:"
    if ($DryRun) {
        Write-Notice "Dry run: would run docker compose ps."
    }
    else {
        & docker compose --project-directory $projectRoot ps 2>&1 | ForEach-Object { Write-Host "    $_" }
        if ($LASTEXITCODE -ne 0) { throw "Could not inspect the Compose project." }
    }
}

Write-Stage 3 "Stopping and removing Version 2 containers and networks"
if (-not $dockerAvailable -or -not $dockerReady) {
    Write-Notice "Compose shutdown skipped because the Docker engine is unavailable."
}
elseif ($DryRun) {
    Write-Notice "Dry run: would run docker compose down --timeout $ShutdownTimeoutSeconds --remove-orphans."
}
else {
    & docker compose --project-directory $projectRoot down --timeout $ShutdownTimeoutSeconds --remove-orphans
    if ($LASTEXITCODE -ne 0) {
        throw "docker compose down failed with exit code $LASTEXITCODE. Runtime secrets were left in place for diagnosis."
    }
    Write-Success "Containers and the Compose network were stopped and removed."
}

Write-Stage 4 "Removing ephemeral runtime secret files"
if ($KeepRuntimeSecrets) {
    Write-Notice "Secret cleanup skipped by -KeepRuntimeSecrets."
    Write-Detail "Secret directory" $secretRoot
}
else {
    Remove-SecretFile -Path $databaseSecretPath -Label "Database-key runtime file"
    Remove-SecretFile -Path $sessionSecretPath -Label "Browser-session runtime file"
    Remove-SecretFile -Path $dandelionSecretPath -Label "DANDELION runner-key runtime file"
    Remove-SecretFile -Path $dandelionArtifactSecretPath -Label "DANDELION managed-artifact runtime file"
    Remove-SecretFile -Path $modelRunnerSecretPath -Label "Model-runner signing-key runtime file"
    Remove-SecretFile -Path $modelCredentialSecretPath -Label "Model-credential encryption-key runtime file"
    if ((Test-Path -LiteralPath $secretRoot) -and -not (Get-ChildItem -LiteralPath $secretRoot -Force)) {
        if ($DryRun) {
            Write-Notice "Dry run: would remove the empty secret directory."
        }
        else {
            Remove-Item -LiteralPath $secretRoot -Force
            Write-Success "Empty runtime-secret directory removed."
        }
    }
}

Write-Stage 5 "Confirming retained data and shutdown result"
Write-Detail "Results retained" (Join-Path $projectRoot "results")
Write-Detail "Input data retained" (Join-Path $projectRoot "data")
Write-Detail "Credential key" "retained in Windows Credential Manager"
if ($dockerReady -and -not $DryRun) {
    $remaining = & docker compose --project-directory $projectRoot ps --status running --quiet 2>$null
    if ($LASTEXITCODE -eq 0 -and -not $remaining) {
        Write-Success "No Version 2 Compose services remain running."
    }
    elseif ($remaining) {
        Write-Host "  WARN One or more Compose services still report as running." -ForegroundColor Yellow
    }
}

$elapsed = (Get-Date) - $startedAt
Write-Host ""
Write-Host "NophiGene Version 2 shutdown completed." -ForegroundColor Green
Write-Detail "Elapsed" ("{0:n1} seconds" -f $elapsed.TotalSeconds)
Write-Host ""
