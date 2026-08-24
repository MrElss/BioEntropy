@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"
python bootstrap_build.py
set "EC=%ERRORLEVEL%"
if not "%EC%"=="0" (
  echo Build failed with exit code %EC%.
  pause
  exit /b %EC%
)
echo Build complete. EXE should be in dist\BioEntropy_Desktop.exe
pause
