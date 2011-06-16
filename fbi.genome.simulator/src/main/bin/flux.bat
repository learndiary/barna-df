@echo off
setlocal ENABLEDELAYEDEXPANSION
::CONFIGURATION
IF "%FLUX_MEM%" == "" set FLUX_MEM="1G"

:: test for a java installation in the path
for %%X in (java.exe) do (set FOUND=%%~$PATH:X)
if not defined FOUND (
   echo No Java installation found, please install Java ^>= 1.6 and make sure it is in your PATH.
   exit /B 1
)

set parent=%~dp0%..\lib
FOR /R %parent% %%G IN (*.jar) DO set CLASSPATH=!CLASSPATH!;%%G
java -Xmx%FLUX_MEM% -cp %CLASSPATH% fbi.commons.flux.Flux %*