#!/bin/bash
# Script to create an index with the following software:
#
# http://www.csee.wvu.edu/~adjeroh/projects/PSA/PSA_software.tar
#
CODE_DIR=$1   # The directory that contains the executables

ROOT_DIR=$(pwd)
DATA_DIR="${ROOT_DIR}/data"
TEXT_FILE="${DATA_DIR}/text-noHeaders.txt"
RESULTS_DIR="${DATA_DIR}/psa"
LOG_FILE="${RESULTS_DIR}/log"
INDEX_COMMAND="./psa.construction"
LENGTH=$(wc -c < ${TEXT_FILE})

rm -rf ${RESULTS_DIR}
mkdir ${RESULTS_DIR}

cd ${CODE_DIR}
/usr/bin/time -v -o ${LOG_FILE}.time ${INDEX_COMMAND} ${TEXT_FILE} ${LENGTH} > ${LOG_FILE}.txt 2>&1
cd ${ROOT_DIR}
mv ${DATA_DIR}/*.psa ${RESULTS_DIR}
