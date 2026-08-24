@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"
if not exist "bootstrap_run.py" (
  echo Missing bootstrap_run.py
  pause
  exit /b 1
)
python bootstrap_run.py
if errorlevel 1 (
  echo.
  echo Run failed. Please check the message above.
  pause
)
