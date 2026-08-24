@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"
set "LOG=%~dp0rebuild_log.txt"
echo BioEntropy rebuild log > "%LOG%"
echo Time: %DATE% %TIME% >> "%LOG%"
echo Folder: %~dp0 >> "%LOG%"
echo. >> "%LOG%"
call "%~dp0Clean_And_Rebuild.cmd" >> "%LOG%" 2>&1
echo.
echo Log saved to: "%LOG%"
pause
