@echo off
REM ========================================
REM DPPUv2 v3.0 All modes and all topologies batch execution script
REM ========================================
REM 
REM usage:
REM   run_all.bat
REM   run_all.bat quick
REM 
REM ========================================

setlocal enabledelayedexpansion

set "QUICK_MODE=%1"

REM Timestamp Generation
for /f "tokens=1-6 delims=/: " %%a in ("%date% %time%") do (
    set "timestamp=%%a%%b%%c_%%d%%e%%f"
)
set "timestamp=%timestamp: =0%"

REM OUTPUT DIR
set "OUTPUT_DIR=run_results_%timestamp%"
if not exist "%OUTPUT_DIR%" mkdir "%OUTPUT_DIR%"

set "LOG_FILE=%OUTPUT_DIR%\summary_%timestamp%.log"
set "CSV_FILE=%OUTPUT_DIR%\results_%timestamp%.csv"

echo ================================================================================ > "%LOG_FILE%"
echo DPPUv2 v3.0 All modes and all topologies batch execution >> "%LOG_FILE%"
echo ================================================================================ >> "%LOG_FILE%"
echo start timestamp: %date% %time% >> "%LOG_FILE%"
echo output dir: %OUTPUT_DIR% >> "%LOG_FILE%"
echo ================================================================================ >> "%LOG_FILE%"
echo. >> "%LOG_FILE%"

echo ================================================================================
echo DPPUv2 v3.0 All modes and all topologies batch execution
echo ================================================================================
echo start timestamp: %date% %time%
echo output dir: %OUTPUT_DIR%
echo ================================================================================
echo.

REM CSV headers
echo Topology,Mode,Variant,Status,LogFile > "%CSV_FILE%"

REM counters
set "TOTAL_COUNT=0"
set "PASSED_COUNT=0"
set "FAILED_COUNT=0"

REM Topology and Mode Definitions
set "TOPOLOGIES=S3S1 T3S1 Nil3S1"

if /i "%QUICK_MODE%"=="quick" (
    set "MODES=MX"
    set "VARIANTS=FULL"
    echo Quick mode: 1 modes x 1 variants = 1 cases per topology
    echo.
) else (
    set "MODES=AX VT MX"
    set "VARIANTS=TT REE FULL"
    echo Full mode: 3 modes x 3 variants = 9 cases per topology
    echo.
)

REM main
for %%T in (%TOPOLOGIES%) do (
    echo [%%T] Starting tests...
    echo [%%T] Starting tests... >> "%LOG_FILE%"
    
    for %%M in (%MODES%) do (
        for %%V in (%VARIANTS%) do (
            set /a TOTAL_COUNT+=1
            set "RUNNER=DPPUv2_runner_%%T_v3.py"
            set "LOG_NAME=%%T_%%M_%%V_%timestamp%.log"
            set "LOG_PATH=%OUTPUT_DIR%\!LOG_NAME!"
            
            echo   [!TOTAL_COUNT!] Running: %%T %%M %%V ...
            
            REM Python execution (capture both standard output and error output)
            set "TEMP_OUTPUT=%OUTPUT_DIR%\temp_%%T_%%M_%%V.txt"
            python !RUNNER! --mode %%M --ny-variant %%V --log-file "!LOG_PATH!" > "!TEMP_OUTPUT!" 2>&1
            
            if !errorlevel! equ 0 (
                REM Check SUCCESS from standard output
                findstr /C:"SUCCESS: Computation completed" "!TEMP_OUTPUT!" > nul
                if !errorlevel! equ 0 (
                    set "STATUS=PASS"
                    set /a PASSED_COUNT+=1
                    echo     PASS
                    echo   [!TOTAL_COUNT!] %%T %%M %%V : PASS >> "%LOG_FILE%"
                    echo %%T,%%M,%%V,PASS,!LOG_NAME! >> "%CSV_FILE%"
                ) else (
                    set "STATUS=FAIL"
                    set /a FAILED_COUNT+=1
                    echo     FAIL ^(No SUCCESS message^)
                    echo   [!TOTAL_COUNT!] %%T %%M %%V : FAIL ^(No SUCCESS message^) >> "%LOG_FILE%"
                    echo %%T,%%M,%%V,FAIL,!LOG_NAME! >> "%CSV_FILE%"
                )
            ) else (
                set "STATUS=FAIL"
                set /a FAILED_COUNT+=1
                echo     FAIL ^(Exit code: !errorlevel!^)
                echo   [!TOTAL_COUNT!] %%T %%M %%V : FAIL ^(Exit code: !errorlevel!^) >> "%LOG_FILE%"
                echo %%T,%%M,%%V,FAIL,!LOG_NAME! >> "%CSV_FILE%"
            )
            
            REM Delete temporary files
            if exist "!TEMP_OUTPUT!" del "!TEMP_OUTPUT!"
        )
    )
    echo.
)

REM summary
echo ================================================================================ >> "%LOG_FILE%"
echo Execution completion summary >> "%LOG_FILE%"
echo ================================================================================ >> "%LOG_FILE%"
echo end timestamp: %date% %time% >> "%LOG_FILE%"
echo Total number of cases: %TOTAL_COUNT% >> "%LOG_FILE%"
echo PASSED: %PASSED_COUNT% >> "%LOG_FILE%"
echo FALIED: %FAILED_COUNT% >> "%LOG_FILE%"
echo ================================================================================ >> "%LOG_FILE%"
echo. >> "%LOG_FILE%"

echo ================================================================================
echo Execution completion summary
echo ================================================================================
echo end timestamp: %date% %time%
echo Total number of cases: %TOTAL_COUNT%
echo PASSED: %PASSED_COUNT%
echo FALIED: %FAILED_COUNT%
echo ================================================================================
echo.
echo result file:
echo   - SUMMARY LOG: %LOG_FILE%
echo   - CSV RESULT : %CSV_FILE%
echo   - Logs for each case: %OUTPUT_DIR%\*.log
echo.

if %FAILED_COUNT% equ 0 (
    echo ALL TEST PASS！
    exit /b 0
) else (
    echo %FAILED_COUNT% Tests failed.
    exit /b 1
)
