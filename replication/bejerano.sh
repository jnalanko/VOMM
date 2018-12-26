#!/bin/bash
# Script to create an index with the following software:
#
# http://bejerano.stanford.edu/pst_v201.zip
#
P=$1  # Variable $\tau_1$ in Figure 1 of the submitted paper.
CODE_DIR=$2  # The directory that contains the executables

ROOT_DIR=$(pwd)
DATA_DIR="${ROOT_DIR}/data"
ALPHABET_FILE="${DATA_DIR}/alphabet.txt"
TEXT_FILE="${DATA_DIR}/text.txt"
RESULTS_DIR="${DATA_DIR}/bejerano-P-${P}"
LOG_FILE="${RESULTS_DIR}/log"
INDEX_COMMAND="./train"
OTHER_CONTEXTS="0 0.001 1.050"
MAX_LENGTH="1000000"

rm -rf ${RESULTS_DIR}
mkdir ${RESULTS_DIR}

cd ${CODE_DIR}
/usr/bin/time -v -o ${LOG_FILE}.time ${INDEX_COMMAND} ${ALPHABET_FILE} ${TEXT_FILE} @ @ ${TEXT_FILE}.pst ${TEXT_FILE}.stats @ 1 ${P} ${OTHER_CONTEXTS} ${MAX_LENGTH} > ${LOG_FILE}.txt 2>&1
cd ${ROOT_DIR}
mv ${DATA_DIR}/*.pst ${RESULTS_DIR}
mv ${DATA_DIR}/*.stats ${RESULTS_DIR}
