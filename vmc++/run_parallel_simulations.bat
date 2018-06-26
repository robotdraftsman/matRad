@echo off
setlocal
set "lock=%temp%\wait%random%.lock"

start "" 9>"%lock%1" /B .\bin\vmc_Windows.exe MCpencilbeam_temp_1
start "" 9>"%lock%2" /B .\bin\vmc_Windows.exe MCpencilbeam_temp_2
start "" 9>"%lock%3" /B .\bin\vmc_Windows.exe MCpencilbeam_temp_3
start "" 9>"%lock%4" /B .\bin\vmc_Windows.exe MCpencilbeam_temp_4

:Wait for all processes to finish (wait until lock files are no longer locked)
1>nul 2>nul ping /n 2 ::1
for %%N in (1 2 3 4) do (
  (call ) 9>"%lock%%%N" || goto :Wait
) 2>nul

del "%lock%*"

