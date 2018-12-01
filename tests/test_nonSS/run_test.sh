#!/bin/bash

rm -rf ./work ./*.html
/lab/weng_scratch/Tomas/nextflow/nextflow/launch.sh run transXpress.nf -resume -profile test_nonSS


