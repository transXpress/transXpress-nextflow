#!/bin/bash

ASSEMBLER=$1
export NXF_ANSI_LOG='false'
GIT_DIR=$(dirname $(readlink -f ./transXpress.nf))"/.git"
GIT_HASH=$(git --git-dir=${GIT_DIR} log --pretty=format:'%H' -n 1)
OUTFILE="transXpress-$ASSEMBLER.stdout.log"
ERRFILE="transXpress-$ASSEMBLER.stderr.log"
echo "$(date)" | tee -a $OUTFILE
echo "transXpress now running. git hash: "${GIT_HASH} | tee -a $OUTFILE
echo "Logs are being written to $OUTFILE and $ERRFILE in the current directory" | tee -a $OUTFILE
echo "Try 'lsof $OUTFILE' if you need to get the process id of the nextflow manager" | tee -a $OUTFILE
echo "transXpress dropping to background on host "$HOSTNAME"..." | tee -a $OUTFILE
/lab/solexa_weng/testtube/nextflow/nextflow-19.02.0-edge-all run transXpress-$ASSEMBLER.nf -w work-$ASSEMBLER --samples 'samples.txt' --species 'species.txt' -resume 1>>$OUTFILE 2>$ERRFILE &
disown


