@echo off
setlocal
cd /d "%~dp0"
if defined VSCMD_VER goto compile
set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist "%VSWHERE%" exit /b 1
for /f "usebackq delims=" %%I in (`"%VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "VCROOT=%%I"
if not defined VCROOT exit /b 1
call "%VCROOT%\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 exit /b 1
:compile
if not exist build mkdir build
pushd build
cl /nologo /O2 /DNOMINMAX /std:c++17 /EHsc ..\src\compare.cpp /Fe:blockdrift.exe /link Psapi.lib
if errorlevel 1 exit /b 1
popd
