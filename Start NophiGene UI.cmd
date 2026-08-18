@echo off
setlocal EnableExtensions
set "PROJECT_ROOT=%~dp0"
set "TARGET_SCRIPT=%PROJECT_ROOT%scripts\start-v2.ps1"

echo.
echo ============================================================
echo  NophiGene Version 2 - verbose start
echo ============================================================
echo  Timestamp    : %DATE% %TIME%
echo  Project      : %PROJECT_ROOT%
echo  PowerShell   : powershell.exe
echo  Script       : %TARGET_SCRIPT%
echo  Runtime      : secured Docker Compose + SQLCipher
echo  Arguments    : forwarded without displaying secret values
echo.

where powershell.exe >nul 2>&1
if errorlevel 1 (
  echo ERROR PowerShell was not found on PATH.
  echo Install or enable Windows PowerShell and try again.
  pause
  exit /b 9009
)

if not exist "%TARGET_SCRIPT%" (
  echo ERROR Start script was not found: %TARGET_SCRIPT%
  pause
  exit /b 2
)

echo Launching Version 2 now...
powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%TARGET_SCRIPT%" %*
set "EXIT_CODE=%ERRORLEVEL%"

echo.
echo ------------------------------------------------------------
echo  Start command finished at %DATE% %TIME%
echo  Exit code: %EXIT_CODE%
echo ------------------------------------------------------------
if not "%EXIT_CODE%"=="0" (
  echo ERROR NophiGene did not start successfully. Review the diagnostics above.
  pause
) else (
  echo NophiGene Version 2 start completed successfully.
)
exit /b %EXIT_CODE%
