#! /bin/bash

##Little bit of help for when running internally in Weng lab
if [ -d "/lab/solexa_weng/testtube/miniconda3/bin" ]; then
echo "Setting up Weng lab miniconda environment..."
source /lab/solexa_weng/testtube/miniconda3/bin/activate
fi

../../run.sh "trinity" "test_notStrandSpecific_local"


