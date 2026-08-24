@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"
if exist "dist\BioEntropy_Desktop.exe" (
  certutil -hashfile "dist\BioEntropy_Desktop.exe" SHA256
) else (
  echo dist\BioEntropy_Desktop.exe not found. Please rebuild first.
)
pause
