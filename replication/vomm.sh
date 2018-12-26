#!/bin/sh
# Script to create an index with the software described 
# in the submitted paper.
#
P=$1  # Variable $\tau_1$ in Figure 1 of the submitted paper.
CODE_DIR=$2   # The directory that contains the executables

ROOT_DIR=$(pwd)
DATA_DIR="${ROOT_DIR}/data"
TEXT_FILE="${DATA_DIR}/text.txt"
RESULTS_DIR="${DATA_DIR}/vomm-P-${P}"
LOG_FILE="${RESULTS_DIR}/log"
INDEX_COMMAND="./build_model_optimized"
OTHER_CONTEXTS="0.001 0.952 1.050"

rm -rf ${RESULTS_DIR}
mkdir ${RESULTS_DIR}

cd ${CODE_DIR}
/usr/bin/time -v -o ${LOG_FILE}.time ${INDEX_COMMAND} --reference-fasta ${TEXT_FILE} --outputdir ${RESULTS_DIR} --four-thresholds ${P} ${OTHER_CONTEXTS} > ${LOG_FILE}.txt 2>&1
cd ${ROOT_DIR}
