#!/bin/bash
# ========================================
# DPPUv2 v3.0 All modes and all topologies batch execution script
# ========================================
# 
# usage:
#   ./run_all.sh
#   ./run_all.sh quick
# 
# ========================================

set -u

QUICK_MODE="${1:-}"

# Timestamp Generation
timestamp=$(date +%Y%m%d_%H%M%S)

# OUTPUT DIR
OUTPUT_DIR="run_results_${timestamp}"
mkdir -p "${OUTPUT_DIR}"

LOG_FILE="${OUTPUT_DIR}/summary_${timestamp}.log"
CSV_FILE="${OUTPUT_DIR}/results_${timestamp}.csv"

{
    echo "================================================================================"
    echo "DPPUv2 v3.0 All modes and all topologies batch execution"
    echo "================================================================================"
    echo "start timestamp: $(date)"
    echo "output dir: ${OUTPUT_DIR}"
    echo "================================================================================"
    echo ""
} > "${LOG_FILE}"

echo "================================================================================"
echo "DPPUv2 v3.0 All modes and all topologies batch execution"
echo "================================================================================"
echo "start timestamp: $(date)"
echo "output dir: ${OUTPUT_DIR}"
echo "================================================================================"
echo ""

# CSV headers
echo "Topology,Mode,Variant,Status,LogFile" > "${CSV_FILE}"

# counters
TOTAL_COUNT=0
PASSED_COUNT=0
FAILED_COUNT=0

# Topology and Mode Definitions
TOPOLOGIES="S3S1 T3S1 Nil3S1"

if [[ "${QUICK_MODE}" == "quick" ]]; then
    MODES="MX"
    VARIANTS="FULL"
    echo "Quick mode: 1 modes x 1 variants = 1 cases per topology"
    echo ""
else
    MODES="AX VT MX"
    VARIANTS="TT REE FULL"
    echo "Full mode: 3 modes x 3 variants = 9 cases per topology"
    echo ""
fi

# main
for T in ${TOPOLOGIES}; do
    echo "[${T}] Starting tests..."
    echo "[${T}] Starting tests..." >> "${LOG_FILE}"
    
    for M in ${MODES}; do
        for V in ${VARIANTS}; do
            ((TOTAL_COUNT++))
            RUNNER="DPPUv2_runner_${T}_v3.py"
            LOG_NAME="${T}_${M}_${V}_${timestamp}.log"
            LOG_PATH="${OUTPUT_DIR}/${LOG_NAME}"
            
            echo "  [${TOTAL_COUNT}] Running: ${T} ${M} ${V} ..."
            
            # Python execution (capture both standard output and error output)
            TEMP_OUTPUT="${OUTPUT_DIR}/temp_${T}_${M}_${V}.txt"
            if python "${RUNNER}" --mode "${M}" --ny-variant "${V}" --log-file "${LOG_PATH}" > "${TEMP_OUTPUT}" 2>&1; then
                # Check SUCCESS from standard output
                if grep -q "SUCCESS: Computation completed" "${TEMP_OUTPUT}"; then
                    STATUS="PASS"
                    ((PASSED_COUNT++))
                    echo "    PASS"
                    echo "  [${TOTAL_COUNT}] ${T} ${M} ${V} : PASS" >> "${LOG_FILE}"
                    echo "${T},${M},${V},PASS,${LOG_NAME}" >> "${CSV_FILE}"
                else
                    STATUS="FAIL"
                    ((FAILED_COUNT++))
                    echo "    FAIL (No SUCCESS message)"
                    echo "  [${TOTAL_COUNT}] ${T} ${M} ${V} : FAIL (No SUCCESS message)" >> "${LOG_FILE}"
                    echo "${T},${M},${V},FAIL,${LOG_NAME}" >> "${CSV_FILE}"
                fi
            else
                EXIT_CODE=$?
                STATUS="FAIL"
                ((FAILED_COUNT++))
                echo "    FAIL (Exit code: ${EXIT_CODE})"
                echo "  [${TOTAL_COUNT}] ${T} ${M} ${V} : FAIL (Exit code: ${EXIT_CODE})" >> "${LOG_FILE}"
                echo "${T},${M},${V},FAIL,${LOG_NAME}" >> "${CSV_FILE}"
            fi
            
            # Delete temporary files
            rm -f "${TEMP_OUTPUT}"
        done
    done
    echo ""
done

# summary
{
    echo "================================================================================"
    echo "Execution completion summary"
    echo "================================================================================"
    echo "end timestamp: $(date)"
    echo "Total number of cases: ${TOTAL_COUNT}"
    echo "PASSED: ${PASSED_COUNT}"
    echo "FALIED: ${FAILED_COUNT}"
    echo "================================================================================"
    echo ""
} >> "${LOG_FILE}"

echo "================================================================================"
echo "Execution completion summary"
echo "================================================================================"
echo "end timestamp: $(date)"
echo "Total number of cases: ${TOTAL_COUNT}"
echo "PASSED: ${PASSED_COUNT}"
echo "FALIED: ${FAILED_COUNT}"
echo "================================================================================"
echo ""
echo "result file:"
echo "  - SUMMARY LOG: ${LOG_FILE}"
echo "  - CSV RESULT : ${CSV_FILE}"
echo "  - Logs for each case: ${OUTPUT_DIR}/*.log"
echo ""

if [[ ${FAILED_COUNT} -eq 0 ]]; then
    echo "ALL TEST PASS！"
    exit 0
else
    echo "${FAILED_COUNT} Tests failed."
    exit 1
fi
