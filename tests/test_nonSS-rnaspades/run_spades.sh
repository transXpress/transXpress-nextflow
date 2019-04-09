#! /bin/bash

##Little bit of help for when running internally in Weng lab
if [ -d "/lab/solexa_weng/testtube/miniconda3/bin" ]; then
echo "Setting up Weng lab miniconda environment..."
source /lab/solexa_weng/testtube/miniconda3/bin/activate
conda activate transxpress-nextflow
else 
echo "Weng lab conda environment not found. Maybe you are running transXpress on non"
echo "Weng lab hardware."
echo "If so, Please make sure you've installed all"
echo "the dependencies & sourced the right conda environment!"
fi

../../run.sh "rnaspades" "test_notStrandSpecific_local"
