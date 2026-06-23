@echo off
REM ============================================================
REM  DEVELOPER BUILD SCRIPT — NOT AN INSTALLER
REM  This compiles Oligolia from source. End users should
REM  download the pre-built Setup.exe from:
REM  https://github.com/moonsoup/oligolia/releases/latest
REM ============================================================
REM Build Oligolia Windows installer
REM Run from project root: build\build_win.bat
REM Requirements: Python 3.11+, pip
REM Output: dist\Oligolia-Setup.exe

setlocal
cd /d "%~dp0\.."

set VENV=backend\.venv
set PYTHON=%VENV%\Scripts\python.exe
set PYI=%VENV%\Scripts\pyinstaller.exe

echo === Oligolia Windows Build ===

REM 1. Create venv if needed
if not exist "%PYTHON%" (
    echo Creating virtual environment...
    python -m venv %VENV%
    %VENV%\Scripts\pip install --quiet -r backend\requirements.txt PyInstaller PyQt6 Pillow
)

REM 2. Generate icons
echo Generating icons...
%PYTHON% assets\make_icon.py

REM 3. PyInstaller
echo Building .exe (this takes 3-8 minutes)...
%PYI% build\oligolia.spec ^
    --distpath dist ^
    --workpath build\work ^
    --noconfirm ^
    --clean

if not exist "dist\Oligolia\Oligolia.exe" (
    echo ERROR: Build failed. Check output above.
    exit /b 1
)
echo Built: dist\Oligolia\Oligolia.exe

REM 4. Create installer with Inno Setup
where iscc >nul 2>&1
if %errorlevel% == 0 (
    echo Creating installer with Inno Setup...
    iscc build\inno_setup.iss
    echo.
    echo Done! Installer: dist\Oligolia-0.1.0-Setup.exe
) else (
    echo.
    echo NOTE: Inno Setup not found. Standalone exe is at dist\Oligolia\Oligolia.exe
    echo Download Inno Setup from https://jrsoftware.org to create a proper installer.
    echo.
    echo Packaging as zip instead...
    powershell -Command "Compress-Archive -Path 'dist\Oligolia\*' -DestinationPath 'dist\Oligolia-0.1.0-win.zip' -Force"
    echo Done! dist\Oligolia-0.1.0-win.zip
)

echo.
echo User instructions:
echo   Option A (installer): Run Oligolia-0.1.0-Setup.exe, follow prompts
echo   Option B (portable):  Unzip Oligolia-0.1.0-win.zip, run Oligolia.exe
