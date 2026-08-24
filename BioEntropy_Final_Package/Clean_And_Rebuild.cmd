@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"
title BioEntropy Clean And Rebuild
echo ============================================
echo BioEntropy - clean old cache and rebuild EXE
echo Folder: %~dp0
echo ============================================
if exist build rmdir /s /q build
if exist dist rmdir /s /q dist
if exist __pycache__ rmdir /s /q __pycache__
for %%f in (*.spec) do del /f /q "%%f"
for /d /r %%d in (__pycache__) do if exist "%%d" rmdir /s /q "%%d"
echo.
echo Rebuilding...
call "%~dp0Rebuild_EXE.cmd"
set "EC=%ERRORLEVEL%"
if "%EC%"=="0" (
  echo.
  echo Rebuild finished successfully.
  if exist "%~dp0dist\BioEntropy_Desktop.exe" echo EXE: "%~dp0dist\BioEntropy_Desktop.exe"
) else (
  echo.
  echo Rebuild failed with exit code %EC%.
)
pause
exit /b %EC%
