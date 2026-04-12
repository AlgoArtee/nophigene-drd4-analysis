[CmdletBinding()]
param()

$ErrorActionPreference = "Stop"
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = (Resolve-Path (Join-Path $scriptDir "..")).Path
$pidFile = Join-Path $repoRoot ".nophigene-ui.pid"

Write-Host ""
Write-Host "Stopping NophiGene local UI" -ForegroundColor Cyan
Write-Host ""

if (-not (Test-Path $pidFile)) {
    Write-Host "No PID file was found. Nothing to stop." -ForegroundColor Yellow
    exit 0
}

$rawPid = Get-Content $pidFile -ErrorAction SilentlyContinue | Select-Object -First 1
if (-not $rawPid) {
    Remove-Item -LiteralPath $pidFile -Force -ErrorAction SilentlyContinue
    Write-Host "The PID file was empty, so it was removed." -ForegroundColor Yellow
    exit 0
}

try {
    $process = Get-Process -Id ([int]$rawPid) -ErrorAction Stop
    Stop-Process -Id $process.Id -Force
    Write-Host "Stopped local UI process $($process.Id)." -ForegroundColor Green
}
catch {
    Write-Host "The tracked process was not running anymore. Cleaning up the PID file." -ForegroundColor Yellow
}

Remove-Item -LiteralPath $pidFile -Force -ErrorAction SilentlyContinue
