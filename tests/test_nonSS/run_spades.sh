#!/bin/bash

export NXF_ANSI_LOG='false'
GIT_DIR=$(dirname $(readlink -f ./transXpress-trinity.nf))"/.git"
GIT_HASH=$(git --git-dir=${GIT_DIR} log --pretty=format:'%H' -n 1)
OUTFILE="transXpress-rnaspades.stdout.log"
ERRFILE="transXpress-rnaspades.stderr.log"
echo "$(date)" | tee -a $OUTFILE
echo "transXpress now running. git hash: "${GIT_HASH} | tee -a $OUTFILE
echo "Logs are being written to $OUTFILE and $ERRFILE in the current directory" | tee -a $OUTFILE
echo "Try 'lsof transXpress.stdout.log' if you need to get the process id of the nextflow manager" | tee -a $OUTFILE
echo "transXpress dropping to background on host "$HOSTNAME"..." | tee -a $OUTFILE
/lab/solexa_weng/testtube/nextflow/nextflow-19.02.0-edge-all run transXpress-rnaspades.nf --samples 'samples.txt' -profile test_nonSS -resume 1>>$OUTFILE 2>$ERRFILE &
disown
