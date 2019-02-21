#!/bin/bash

export NXF_ANSI_LOG='false'
GIT_DIR=$(dirname $(readlink -f ./transXpress-trinity.nf))"/.git"
GIT_HASH=$(git --git-dir=${GIT_DIR} log --pretty=format:'%H' -n 1)
echo "$(date)" | tee -a transXpress.stdout.log
echo "transXpress now running. git hash: "${GIT_HASH} | tee -a transXpress.stdout.log
echo "Logs are being written to transXpress.stdout.log and transXpress.stderr.log in the current directory" | tee -a transXpress.stdout.log
echo "Try 'lsof transXpress.stdout.log' if you need to get the process id of the nextflow manager" | tee -a transXpress.stdout.log
echo "transXpress dropping to background on host "$HOSTNAME"..." | tee -a transXpress.stdout.log
/lab/solexa_weng/testtube/nextflow/nextflow-19.02.0-edge-all run transXpress-trinity.nf -profile test_nonSS -resume 1>>transXpress-trinity.stdout.log 2>transXpress-trinity.stderr.log &
disown
