@echo off
setlocal
cd /d "%~dp0"
set STAGE=%~1
if "%STAGE%"=="" set STAGE=derive
if "%STAGE%"=="setup" (Rscript code\00_setup.R & exit /b)
if "%STAGE%"=="prepare" (Rscript code\01_prepare.R & exit /b)
if "%STAGE%"=="derive" (Rscript code\02_derive.R & exit /b)
if "%STAGE%"=="test" (Rscript tests\run_tests.R & exit /b)
if "%STAGE%"=="model" (
  if "%~2"=="" (echo Usage: run.bat model path\to\consensus_mapping.csv & exit /b 1)
  Rscript code\03_model.R "%~2"
  exit /b
)
echo Usage: run.bat [setup^|prepare^|derive^|model MAPPING^|test]
exit /b 1
