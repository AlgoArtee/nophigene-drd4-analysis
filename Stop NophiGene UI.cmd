@echo off
setlocal EnableExtensions
set "PROJECT_ROOT=%~dp0"
set "TARGET_SCRIPT=%PROJECT_ROOT%scripts\stop-v2.ps1"

echo.
echo ============================================================
echo  NophiGene Version 2 - verbose stop
echo ============================================================
echo  Timestamp    : %DATE% %TIME%
echo  Project      : %PROJECT_ROOT%
echo  Script       : %TARGET_SCRIPT%
echo  Data policy  : input data and results are retained
echo.

where powershell.exe >nul 2>&1
if errorlevel 1 (
  echo ERROR PowerShell was not found on PATH.
  pause
  exit /b 9009
)
if not exist "%TARGET_SCRIPT%" (
  echo ERROR Stop script was not found: %TARGET_SCRIPT%
  pause
  exit /b 2
)

echo Stopping Version 2 services and cleaning runtime secrets...
powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%TARGET_SCRIPT%" %*
set "EXIT_CODE=%ERRORLEVEL%"

echo.
echo ------------------------------------------------------------
echo  Stop command finished at %DATE% %TIME%
echo  Exit code: %EXIT_CODE%
echo ------------------------------------------------------------
if not "%EXIT_CODE%"=="0" (
  echo ERROR NophiGene did not stop cleanly. Review the diagnostics above.
  pause
) else (
  echo NophiGene Version 2 shutdown completed successfully.
)
exit /b %EXIT_CODE%
