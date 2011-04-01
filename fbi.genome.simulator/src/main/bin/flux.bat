@echo off
java  -Xmx1G -XX:-UseBiasedLocking -DwrapperDir="'%~dp0'" -jar "%~dp0..\lib\FluxSimulator.jar" %*
