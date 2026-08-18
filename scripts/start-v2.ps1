[CmdletBinding()]
param(
    [ValidateRange(1, 65535)]
    [int]$Port = 8766,
    [switch]$NoBrowser,
    [switch]$SkipBuild,
    [ValidateRange(10, 600)]
    [int]$StartupTimeoutSeconds = 180,
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"
$ProgressPreference = "SilentlyContinue"
$startedAt = Get-Date
$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$composePath = Join-Path $projectRoot "docker-compose.yml"
$dockerConfigDir = Join-Path $projectRoot ".docker-local"
$secretRoot = Join-Path $projectRoot ".runtime-secrets"
$databaseSecretPath = Join-Path $secretRoot "database-key"
$sessionSecretPath = Join-Path $secretRoot "session-token"
$dandelionSecretPath = Join-Path $secretRoot "dandelion-runner-key"
$dandelionArtifactSecretPath = Join-Path $secretRoot "dandelion-artifact-key"
$vaultResource = "NophiGene-v2"
$vaultUser = "database-key"
$dandelionVaultUser = "dandelion-runner-key"
$dandelionArtifactVaultUser = "dandelion-artifact-key"
$healthUrl = "http://127.0.0.1:$Port/api/v2/health"
$publicUrl = "http://127.0.0.1:$Port/"
$stageCount = 7

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

function Invoke-DockerCommand {
    param(
        [string[]]$Arguments,
        [string]$Description,
        [switch]$AllowFailure,
        [switch]$RunDuringDryRun
    )
    Write-Notice "$Description"
    Write-Verbose ("docker " + ($Arguments -join " "))
    if ($DryRun -and -not $RunDuringDryRun) {
        Write-Notice "Dry run: command was not executed."
        return 0
    }
    & docker @Arguments
    $exitCode = $LASTEXITCODE
    if ($exitCode -ne 0 -and -not $AllowFailure) {
        throw "$Description failed with exit code $exitCode."
    }
    return $exitCode
}

function Test-DockerEngine {
    if ($DryRun) { return $true }
    & docker info *> $null
    return $LASTEXITCODE -eq 0
}

function Remove-RuntimeSecretFiles {
    foreach ($path in @($databaseSecretPath, $sessionSecretPath, $dandelionSecretPath, $dandelionArtifactSecretPath)) {
        if (Test-Path -LiteralPath $path) {
            Remove-Item -LiteralPath $path -Force -ErrorAction SilentlyContinue
        }
    }
    if ((Test-Path -LiteralPath $secretRoot) -and -not (Get-ChildItem -LiteralPath $secretRoot -Force)) {
        Remove-Item -LiteralPath $secretRoot -Force -ErrorAction SilentlyContinue
    }
}

function Show-ComposeDiagnostics {
    if ($DryRun) { return }
    Write-Host ""
    Write-Host "Docker service state:" -ForegroundColor Yellow
    & docker compose --project-directory $projectRoot ps 2>&1 | ForEach-Object { Write-Host "  $_" }
    Write-Host ""
    Write-Host "Recent app logs (secrets are never printed by the launcher):" -ForegroundColor Yellow
    & docker compose --project-directory $projectRoot logs --tail 80 app 2>&1 | ForEach-Object { Write-Host "  $_" }
}

Write-Host ""
Write-Host "============================================================" -ForegroundColor Green
Write-Host " NophiGene Version 2 - secure local launcher" -ForegroundColor Green
Write-Host "============================================================" -ForegroundColor Green
Write-Detail "Started" $startedAt.ToString("yyyy-MM-dd HH:mm:ss K")
Write-Detail "Project" $projectRoot
Write-Detail "Compose file" $composePath
Write-Detail "Docker config" $dockerConfigDir
Write-Detail "Host binding" "127.0.0.1:$Port"
Write-Detail "Build policy" $(if ($SkipBuild) { "reuse existing image" } else { "build before start" })
Write-Detail "Browser" $(if ($NoBrowser) { "do not open" } else { "open authenticated session" })
Write-Detail "Dry run" $DryRun.IsPresent
Write-Host ""

Write-Stage 1 "Validating launcher inputs and repository files"
if (-not (Test-Path -LiteralPath $composePath -PathType Leaf)) {
    throw "Compose file not found: $composePath"
}
if (-not (Test-Path -LiteralPath (Join-Path $projectRoot "Dockerfile") -PathType Leaf)) {
    throw "Dockerfile not found under $projectRoot."
}
foreach ($runtimeDirectory in @(
    (Join-Path $projectRoot "data\dandelion"),
    (Join-Path $projectRoot "results\dandelion")
)) {
    if (-not (Test-Path -LiteralPath $runtimeDirectory)) {
        if ($DryRun) {
            Write-Notice "Dry run: would create runtime directory $runtimeDirectory."
        }
        else {
            New-Item -ItemType Directory -Path $runtimeDirectory -Force | Out-Null
        }
    }
}
if (-not (Test-Path -LiteralPath $dockerConfigDir)) {
    if ($DryRun) {
        Write-Notice "Dry run: would create isolated Docker config directory $dockerConfigDir."
    }
    else {
        New-Item -ItemType Directory -Path $dockerConfigDir -Force | Out-Null
    }
}
$env:DOCKER_CONFIG = $dockerConfigDir
Write-Success "Repository and Compose files are present."
Write-Detail "DANDELION imports" (Join-Path $projectRoot "data\dandelion")

Write-Stage 2 "Checking Docker CLI, Compose, and engine readiness"
if (-not (Get-Command docker -ErrorAction SilentlyContinue)) {
    throw "Docker CLI was not found on PATH. Install or start Docker Desktop, then retry."
}
Invoke-DockerCommand -Arguments @("compose", "version") -Description "Reading Docker Compose version" -RunDuringDryRun | Out-Null
if ($DryRun) {
    Write-Notice "Dry run: Docker engine readiness probe skipped."
}
elseif (-not (Test-DockerEngine)) {
    throw "Docker Desktop is installed but its Linux engine is not ready. Start Docker Desktop and wait for it to report that the engine is running."
}
else {
    Write-Success "Docker engine is ready."
}

Write-Stage 3 "Validating the resolved Compose configuration"
$env:NOPHIGENE_PORT = [string]$Port
Invoke-DockerCommand -Arguments @("compose", "--project-directory", $projectRoot, "config", "--quiet") -Description "Validating docker-compose.yml" -RunDuringDryRun | Out-Null
Write-Success "Compose configuration is valid."

Write-Stage 4 "Preparing access-restricted runtime secrets"
if ($DryRun) {
    Write-Notice "Dry run: Credential Manager and runtime secret files were not changed."
    $databaseKey = ""
    $sessionToken = ""
}
else {
    Add-Type -AssemblyName System.Runtime.WindowsRuntime
    $vaultType = [Windows.Security.Credentials.PasswordVault,Windows.Security.Credentials,ContentType=WindowsRuntime]
    $credentialType = [Windows.Security.Credentials.PasswordCredential,Windows.Security.Credentials,ContentType=WindowsRuntime]
    $vault = [Activator]::CreateInstance($vaultType)
    try {
        $credential = $vault.Retrieve($vaultResource, $vaultUser)
        $credential.RetrievePassword()
        $databaseKey = $credential.Password
        Write-Notice "Reused the database key stored in Windows Credential Manager."
    }
    catch {
        $databaseBytes = New-Object byte[] 48
        [System.Security.Cryptography.RandomNumberGenerator]::Fill($databaseBytes)
        $databaseKey = [Convert]::ToBase64String($databaseBytes)
        $credential = [Activator]::CreateInstance($credentialType, @($vaultResource, $vaultUser, $databaseKey))
        $vault.Add($credential)
        Write-Notice "Created a new database key in Windows Credential Manager."
    }

    $sessionBytes = New-Object byte[] 32
    [System.Security.Cryptography.RandomNumberGenerator]::Fill($sessionBytes)
    $sessionToken = [Convert]::ToBase64String($sessionBytes).Replace("+", "-").Replace("/", "_").TrimEnd("=")
    try {
        $dandelionCredential = $vault.Retrieve($vaultResource, $dandelionVaultUser)
        $dandelionCredential.RetrievePassword()
        $dandelionRunnerKey = $dandelionCredential.Password
        Write-Notice "Reused the offline-runner signing key stored in Windows Credential Manager."
    }
    catch {
        $dandelionBytes = New-Object byte[] 48
        [System.Security.Cryptography.RandomNumberGenerator]::Fill($dandelionBytes)
        $dandelionRunnerKey = [Convert]::ToBase64String($dandelionBytes)
        $dandelionCredential = [Activator]::CreateInstance($credentialType, @($vaultResource, $dandelionVaultUser, $dandelionRunnerKey))
        $vault.Add($dandelionCredential)
        Write-Notice "Created a new offline-runner signing key in Windows Credential Manager."
    }
    try {
        $artifactCredential = $vault.Retrieve($vaultResource, $dandelionArtifactVaultUser)
        $artifactCredential.RetrievePassword()
        $dandelionArtifactKey = $artifactCredential.Password
        Write-Notice "Reused the managed-artifact encryption key stored in Windows Credential Manager."
    }
    catch {
        $artifactBytes = New-Object byte[] 48
        [System.Security.Cryptography.RandomNumberGenerator]::Fill($artifactBytes)
        $dandelionArtifactKey = [Convert]::ToBase64String($artifactBytes)
        $artifactCredential = [Activator]::CreateInstance($credentialType, @($vaultResource, $dandelionArtifactVaultUser, $dandelionArtifactKey))
        $vault.Add($artifactCredential)
        Write-Notice "Created a new managed-artifact encryption key in Windows Credential Manager."
    }
    New-Item -ItemType Directory -Path $secretRoot -Force | Out-Null
    [System.IO.File]::WriteAllText($databaseSecretPath, $databaseKey, [System.Text.UTF8Encoding]::new($false))
    [System.IO.File]::WriteAllText($sessionSecretPath, $sessionToken, [System.Text.UTF8Encoding]::new($false))
    [System.IO.File]::WriteAllText($dandelionSecretPath, $dandelionRunnerKey, [System.Text.UTF8Encoding]::new($false))
    [System.IO.File]::WriteAllText($dandelionArtifactSecretPath, $dandelionArtifactKey, [System.Text.UTF8Encoding]::new($false))

    $currentIdentity = [System.Security.Principal.WindowsIdentity]::GetCurrent().Name
    foreach ($secretPath in @($databaseSecretPath, $sessionSecretPath, $dandelionSecretPath, $dandelionArtifactSecretPath)) {
        $acl = New-Object System.Security.AccessControl.FileSecurity
        $acl.SetAccessRuleProtection($true, $false)
        $rule = New-Object System.Security.AccessControl.FileSystemAccessRule(
            $currentIdentity,
            [System.Security.AccessControl.FileSystemRights]::FullControl,
            [System.Security.AccessControl.AccessControlType]::Allow
        )
        $acl.AddAccessRule($rule)
        Set-Acl -LiteralPath $secretPath -AclObject $acl
    }
    Write-Success "Database, browser-session, offline-runner, and managed-artifact secrets were written with a restricted ACL."
    Write-Detail "Secret directory" $secretRoot
    Write-Detail "Secret values" "redacted"
}

Write-Stage 5 "Building and starting the Version 2 services"
$composeArguments = @("compose", "--project-directory", $projectRoot, "up", "-d", "--remove-orphans")
if (-not $SkipBuild) { $composeArguments += "--build" }
try {
    Invoke-DockerCommand -Arguments $composeArguments -Description "Starting the app service" | Out-Null
}
catch {
    Show-ComposeDiagnostics
    Remove-RuntimeSecretFiles
    throw
}
if ($DryRun) {
    Write-Stage 6 "Skipping health polling in dry-run mode"
    Write-Notice "No service was started, so no health request was sent."
    Write-Stage 7 "Skipping browser launch in dry-run mode"
    Write-Notice "No authenticated URL was created or opened."
    Write-Host ""
    Write-Host "Dry run completed successfully; no container, secret, or browser state was changed." -ForegroundColor Green
    exit 0
}
Write-Success "Compose accepted the start request."

Write-Stage 6 "Waiting for the authenticated web application to become healthy"
$deadline = (Get-Date).AddSeconds($StartupTimeoutSeconds)
$attempt = 0
$healthPayload = $null
while ((Get-Date) -lt $deadline) {
    $attempt++
    try {
        $response = Invoke-WebRequest -Uri $healthUrl -UseBasicParsing -TimeoutSec 5
        if ($response.StatusCode -eq 200) {
            $healthPayload = $response.Content | ConvertFrom-Json
            if ($healthPayload.status -in @("ok", "degraded")) { break }
        }
    }
    catch {
        if (($attempt % 5) -eq 0) {
            Write-Notice "Still waiting (attempt $attempt); checking container state."
            & docker compose --project-directory $projectRoot ps --status running 2>&1 | ForEach-Object { Write-Host "    $_" }
        }
    }
    Start-Sleep -Seconds 2
}
if (-not $healthPayload -or $healthPayload.status -notin @("ok", "degraded")) {
    Show-ComposeDiagnostics
    throw "The app did not become healthy at $healthUrl within $StartupTimeoutSeconds seconds."
}
Write-Success "Health endpoint returned '$($healthPayload.status)'."
Write-Detail "Worker" $(if ($healthPayload.worker.alive) { "running" } else { "degraded" })
Write-Detail "Queue depth" ([string]$healthPayload.worker.queue_depth)
Write-Detail "Database schema" ([string]$healthPayload.database.schema)
Write-Detail "GPU preflight" ([string]$healthPayload.gpu.status)
Write-Detail "DANDELION worker" ([string]$healthPayload.dandelion.worker.status)
Write-Detail "DANDELION queue" ([string]$healthPayload.dandelion.queue_depth)
Write-Notice "Current service state:"
& docker compose --project-directory $projectRoot ps 2>&1 | ForEach-Object { Write-Host "    $_" }

Write-Stage 7 "Opening the one-time authenticated browser session"
if ($NoBrowser) {
    Write-Notice "Browser launch skipped by -NoBrowser."
    Write-Notice "To authenticate manually, rerun without -NoBrowser; the session token is intentionally not printed."
}
else {
    $authenticatedUrl = "$publicUrl`?access_token=$([Uri]::EscapeDataString($sessionToken))"
    Write-Notice "Opening the local URL with a redacted one-time token."
    Start-Process $authenticatedUrl | Out-Null
    Write-Success "Browser launch requested."
}

$elapsed = (Get-Date) - $startedAt
Write-Host ""
Write-Host "NophiGene Version 2 is ready." -ForegroundColor Green
Write-Detail "URL" $publicUrl
Write-Detail "Elapsed" ("{0:n1} seconds" -f $elapsed.TotalSeconds)
Write-Detail "Stop command" ".\scripts\stop-v2.ps1"
Write-Detail "Persistent data" (Join-Path $projectRoot "results")
Write-Notice "Runtime secret files remain access-restricted until the stop script removes them."
Write-Host ""
