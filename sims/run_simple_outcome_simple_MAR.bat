:: run all estimators, sequentially, for the base case scenario with 12% outcome proportion, 40% missing-data proportion, simple outcome model and simple MAR missingness
@echo off
setlocal enabledelayedexpansion

:: grab R version for the current machine
For /F "Delims=" %%0 In ('where /r "C:\Program Files\R" Rscript.exe') do set scriptpath="%%~0"
echo Using R executable !scriptpath!

:: pull command-line arguments for the Y and M scenario
set Y=1.1
set M=1.1

:: set up arguments to pass
set nreps=2500
set n=10000
set data_only=0
set use_cached_datasets=0
set S=1
set X=1
:: send output to directory
set outdir="%cd%\rout"

if not exist %outdir% mkdir %outdir%

:: run the simulation
:: outermost loop is over M
:: loop over estimators
for %%E in ("cc_oracle" "cc_population" "cc_noW" "cc" "ipw" "gr" "mice" "rf" "xgb" "ipcw-tmle_m" "ipcw-tmle_mto" "ipcw-a-tmle_m" "ipcw-a-tmle_mto") do (
  set str="%%E"
  set nice_E=!str:;=_%!

  echo Running Y scenario !Y!, M scenario !M!, X scenario !X!, Seed !S!, estimator !nice_E!

  set this_outfile=!outdir!\output_m!M!_y!Y!_x!X!_n!n!_s!S!_est_!nice_E!.out

  !scriptpath! 04_main.R --nreps-total !nreps! --n !n! --yscenario !Y! --mscenario !M! --xscenario !X! --estimator %%E --data-only !data_only! --seed !S! --use-cached-datasets !use_cached_datasets! 1>!this_outfile! 2>&1
)


