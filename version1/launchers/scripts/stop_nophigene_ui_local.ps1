[CmdletBinding()]
param()

$ErrorActionPreference = "Stop"
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = (Resolve-Path (Join-Path $scriptDir "..")).Path
$pidFile = Join-Path $repoRoot ".nophigene-ui.pid"
$portFile = Join-Path $repoRoot ".nophigene-ui.port"
$dataDir = Join-Path $repoRoot "data"
$referenceDir = Join-Path $dataDir "reference\hg38"
$extractedDir = Join-Path $dataDir "extracted"

Write-Host ""
Write-Host "Stopping NophiGene local UI" -ForegroundColor Cyan
Write-Host "Persistent data: $dataDir"
Write-Host ""

function Test-NophiGeneReady {
    param([int]$WebPort)
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

function Get-ListeningProcessId {
    param([int]$WebPort)
    $pattern = "^\s*TCP\s+\S+:$WebPort\s+\S+\s+\S+\s+(?<ProcessId>\d+)\s*$"
    foreach ($line in (netstat -ano -p tcp)) {
        if ($line -match $pattern) {
            return [int]$Matches.ProcessId
        }
    }
    return $null
}

if (-not (Test-Path $pidFile)) {
    Remove-Item -LiteralPath $portFile -Force -ErrorAction SilentlyContinue
    Write-Host "No PID file was found. Nothing to stop." -ForegroundColor Yellow
    exit 0
}

$rawProcessIds = @(
    Get-Content $pidFile -ErrorAction SilentlyContinue |
        ForEach-Object {
            try { [int]$_ } catch { $null }
        } |
        Where-Object { $null -ne $_ }
)
$trackedPort = $null
if (Test-Path $portFile) {
    try {
        $trackedPort = [int](Get-Content $portFile -ErrorAction Stop | Select-Object -First 1)
    }
    catch {
        $trackedPort = $null
    }
}
if ($trackedPort -and (Test-NophiGeneReady -WebPort $trackedPort)) {
    $listenerProcessId = Get-ListeningProcessId -WebPort $trackedPort
    if ($listenerProcessId -and $listenerProcessId -notin $rawProcessIds) {
        $rawProcessIds += $listenerProcessId
    }
}
if ($rawProcessIds.Count -eq 0) {
    Remove-Item -LiteralPath $pidFile -Force -ErrorAction SilentlyContinue
    Remove-Item -LiteralPath $portFile -Force -ErrorAction SilentlyContinue
    Write-Host "The PID file was empty, so it was removed." -ForegroundColor Yellow
    exit 0
}

$stoppedProcessIds = @()
foreach ($processIdValue in ($rawProcessIds | Select-Object -Unique)) {
    try {
        $process = Get-Process -Id $processIdValue -ErrorAction Stop
        Stop-Process -Id $process.Id -Force
        $stoppedProcessIds += $process.Id
    }
    catch {
        continue
    }
}

if ($stoppedProcessIds.Count -gt 0) {
    Write-Host "Stopped local UI process(es): $($stoppedProcessIds -join ', ')." -ForegroundColor Green
}
else {
    Write-Host "The tracked processes were not running anymore. Cleaning up the launcher files." -ForegroundColor Yellow
}

Remove-Item -LiteralPath $pidFile -Force -ErrorAction SilentlyContinue
Remove-Item -LiteralPath $portFile -Force -ErrorAction SilentlyContinue
Write-Host "Reference files remain in $referenceDir."
Write-Host "Extracted VCFs remain in $extractedDir."
