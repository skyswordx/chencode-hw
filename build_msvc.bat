@echo off
REM ============================================
REM Channel Coding Simulation - MSVC Build Script
REM Usage: 
REM   1. Open "Developer Command Prompt for VS"
REM   2. cd to project root directory
REM   3. Run build_msvc.bat
REM ============================================

echo.
echo ============================================
echo  Building Channel Coding Simulator (MSVC)
echo ============================================
echo.

REM Create output directory for CSV files
if not exist output mkdir output

REM Clean up any previous build artifacts first
echo [1/3] Cleaning previous build artifacts...
del /Q *.obj 2>nul
REM Skip exe deletion to avoid conflict with running simulations

REM Set source files
set SOURCES=src\main.c src\convolutional_code.c src\trellis.c src\viterbi.c src\bcjr.c src\log_map_decoder.c src\turbo_code.c src\ccsds_turbo.c src\sim_runner.c src\csv_export.c

REM Set compiler flags
REM /nologo  - Suppress copyright message
REM /O2      - Optimize for speed
REM /W3      - Warning level 3
REM /wd4819  - Suppress Unicode warning (Chinese comments)
REM /I       - Include directory
set CFLAGS=/nologo /O2 /W3 /wd4819 /I include

REM Compile
echo [2/3] Compiling source files...
cl %CFLAGS% %SOURCES% /Fe:chencode_r12.exe

if %ERRORLEVEL% EQU 0 (
    echo.
    echo [3/3] Cleaning intermediate files...
    del /Q *.obj 2>nul
    
    echo.
    echo ============================================
    echo  BUILD SUCCESSFUL!
    echo  Output: chencode_r12.exe
    echo ============================================
    echo.
    echo To run: .\chencode_r12.exe
    echo.
) else (
    echo.
    echo ============================================
    echo  BUILD FAILED - Check errors above
    echo ============================================
    
    REM Clean up partial build artifacts
    del /Q *.obj 2>nul
)
