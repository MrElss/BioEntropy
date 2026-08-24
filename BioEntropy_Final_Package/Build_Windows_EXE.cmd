@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"
title BioEntropy - Build Windows EXE

echo ============================================================
echo  BioEntropy - One-click Windows EXE builder
echo  Folder: %~dp0
echo ============================================================
echo.

REM --- Locate a Python 3 interpreter (prefer the py launcher) ---
set "PYCMD="
py -3 --version >nul 2>&1 && set "PYCMD=py -3"
if not defined PYCMD (
  python --version >nul 2>&1 && set "PYCMD=python"
)
if not defined PYCMD (
  echo [ERROR] Python 3 was not found on this computer.
  echo.
  echo   1^) Install Python 3.10 or newer from:
  echo        https://www.python.org/downloads/
  echo   2^) During installation, tick "Add python.exe to PATH".
  echo   3^) Then run this script again.
  echo.
  pause
  exit /b 1
)

echo Using Python:
%PYCMD% --version
echo.

REM --- Build (bootstrap_build.py creates .venv_build, installs deps, runs PyInstaller) ---
echo The first build downloads PyInstaller and the dependencies into
echo .venv_build and may take several minutes. An internet connection
echo is required the first time.
echo.

%PYCMD% bootstrap_build.py
set "EC=%ERRORLEVEL%"
echo.

if not "%EC%"=="0" (
  echo ============================================================
  echo  BUILD FAILED  ^(exit code %EC%^)
  echo ------------------------------------------------------------
  echo  Scroll up to read the error. Common causes:
  echo    - No internet for the first-time dependency download
  echo    - Antivirus blocking PyInstaller
  echo    - Python too old ^(needs 3.10+^)
  echo ============================================================
  pause
  exit /b %EC%
)

if exist "%~dp0dist\BioEntropy_Desktop.exe" (
  echo ============================================================
  echo  BUILD COMPLETE
  echo ------------------------------------------------------------
  echo  EXE: "%~dp0dist\BioEntropy_Desktop.exe"
  echo  Double-click that file to launch BioEntropy.
  echo ============================================================
) else (
  echo [WARN] Build reported success but dist\BioEntropy_Desktop.exe
  echo        was not found. Scroll up for details.
)

pause
exit /b %EC%
