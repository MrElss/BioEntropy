@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"
if exist "dist\BioEntropy_V1_0.exe" (
  certutil -hashfile "dist\BioEntropy_V1_0.exe" SHA256
) else (
  echo dist\BioEntropy_V1_0.exe not found. Please rebuild first.
)
pause
