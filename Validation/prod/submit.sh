#!/bin/bash

for DS in $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_*.txt; do
  NAME=`basename $DS`
  NAME=${NAME/.txt/}
  NAME=${NAME/dataset_/}
  create-batch --jobName $NAME --nJobs 1 --fileList $DS --cfg run_validation_cfg.py -n
done
