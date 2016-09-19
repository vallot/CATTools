#!/usr/bin/env python

import os
datasetList=  open("dataset.txt").readlines()

for filename in datasetList :
  filename = filename.strip()
  if ( filename == "") : continue
  print filename
  dataset = filename.replace("dataset_","").replace(".txt","")
  print dataset
  cmd = "create-batch --jobName %s --fileList ../../../data/dataset/%s --maxFiles 20 --cfg run_CATJetCharge_cfg.py"%(dataset,filename) 
  if ( dataset.find("TT") != -1 ) :
    cmd += " --args \"isTT=True\""
  print cmd
  os.system(cmd)
