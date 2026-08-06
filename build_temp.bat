@echo off
REM ============================================
REM Temporary Build Script (Avoids file lock)
REM ============================================

echo.
echo ============================================
echo  Building Channel Coding Simulator (TEMP)
echo ============================================
echo.

if not exist output mkdir output

echo [1/3] Cleaning previous build artifacts...
del /Q *.obj 2>nul
del /Q chencode_temp.exe 2>nul

set SOURCES=src\main.c src\convolutional_code.c src\trellis.c src\viterbi.c src\bcjr.c src\log_map_decoder.c src\turbo_code.c src\ccsds_turbo.c src\sim_runner.c src\csv_export.c
set CFLAGS=/nologo /O2 /W3 /wd4819 /I include

echo [2/3] Compiling source files...
cl %CFLAGS% %SOURCES% /Fe:chencode_temp.exe

if %ERRORLEVEL% EQU 0 (
    echo.
    echo [3/3] Cleaning intermediate files...
    del /Q *.obj 2>nul
    
    echo.
    echo ============================================
    echo  BUILD SUCCESSFUL!
    echo  Output: chencode_temp.exe
    echo ============================================
) else (
    echo.
    echo  BUILD FAILED!
)
