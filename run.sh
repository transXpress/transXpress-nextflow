#!/bin/bash

GIT_DIR=$(dirname $(readlink -f ./transXpress.nf))"/.git"
GIT_HASH=$(git --git-dir=${GIT_DIR} log --pretty=format:'%H' -n 1)
echo "$(date)" | tee -a transXpress.stdout.log
echo "transXpress now running. git hash: "${GIT_HASH} | tee -a transXpress.stdout.log
/lab/weng_scratch/Tomas/nextflow/nextflow/launch.sh run transXpress.nf -resume 1>>transXpress.stdout.log 2>>transXpress.stderr.log &
disown
echo "transXpress dropping to background on host "$HOSTNAME"..." | tee -a transXpress.stdout.log
echo "Logs are being written to transXpress.stdout.log and transXpress.stderr.log in the current directory" | tee -a transXpress.stdout.log
echo "Try 'lsof transXpress.stdout.log' if you need to get the process id of the nextflow manager" | tee -a transXpress.stdout.log


