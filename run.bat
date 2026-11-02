@echo off
setlocal
cd /d "%~dp0"
set STAGE=%~1
if "%STAGE%"=="" (echo See README.md for setup, stage1 and coordinator-issued stage2 commands. & exit /b 0)
if "%STAGE%"=="setup" (Rscript --vanilla scripts\setup.R & exit /b)
if "%STAGE%"=="prepare" (Rscript code\01_prepare.R & exit /b)
if "%STAGE%"=="derive" (Rscript code\02_derive.R & exit /b)
if "%STAGE%"=="supplement" (Rscript code\04_supplement.R & exit /b)
if "%STAGE%"=="test" (
  Rscript tests\run_tests.R
  if errorlevel 1 exit /b 1
  Rscript tests\apache_adapted.R
  if errorlevel 1 exit /b 1
  Rscript tests\sofa_recovery.R
  if errorlevel 1 exit /b 1
  Rscript tests\site_supplement.R
  if errorlevel 1 exit /b 1
  Rscript tests\site_return_collection.R
  if errorlevel 1 exit /b 1
  Rscript tests\publication_descriptives.R
  if errorlevel 1 exit /b 1
  Rscript tests\competing_event_curves.R
  if errorlevel 1 exit /b 1
  Rscript tests\survival_tdauc.R
  if errorlevel 1 exit /b 1
  Rscript tests\exact_aggregates.R
  if errorlevel 1 exit /b 1
  Rscript tests\pooled_survival.R
  if errorlevel 1 exit /b 1
  Rscript tests\pooled_roc.R
  if errorlevel 1 exit /b 1
  Rscript tests\four_endpoint_cox.R
  if errorlevel 1 exit /b 1
  Rscript tests\nearest_assignment.R
  if errorlevel 1 exit /b 1
  Rscript tests\model_supplement_diagnostics.R
  exit /b
)
if "%STAGE%"=="model" (
  if "%~2"=="" (echo Usage: run.bat model path\to\consensus_mapping.csv & exit /b 1)
  Rscript code\03_model.R "%~2"
  exit /b
)
echo Usage: run.bat [setup^|prepare^|derive^|supplement^|model MAPPING^|test]
exit /b 1
