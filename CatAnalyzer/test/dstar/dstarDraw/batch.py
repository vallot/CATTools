#!/usr/bin/env python

import os,sys

if ( len(sys.argv) != 2) : 
  print "Wrong argument"
  os.exit(-1)

if ( not os.path.isfile(sys.argv[1])) :
  print "Wrong filename"
  os.exit(-1)
  
datasetList =  open(sys.argv[1]).readlines()


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
