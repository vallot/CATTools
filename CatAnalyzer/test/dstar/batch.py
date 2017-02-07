#!/usr/bin/env python

import os,sys

#strDest = "/store/user/quark2930/dilepton_mass_v801_16102701"
strFilenameRootdir = "rootdirPath.txt"

fRootDir = open(strFilenameRootdir, "r")
strDest = fRootDir.readline().split("\r")[0].split("\n")[0]
fRootDir.close()

if ( len(sys.argv) != 2) : 
  print "Wrong argument"
  sys.exit(-1)

if ( not os.path.isfile(sys.argv[1])) :
  print "Wrong filename"
  sys.exit(-1)
  
datasetList =  open(sys.argv[1]).readlines()


for filename in datasetList :
  filename = filename.strip()
  if ( filename == "") : continue
  print filename
  dataset = filename.replace("dataset_","").replace(".txt","")
  print dataset
  cmd = "create-batch --jobName samples_%s --fileList ../../data/dataset/%s --maxFiles 20 --cfg run_CATDstar_cfg.py --transferDest %s/%s"%(dataset,filename,strDest,dataset) 
  if ( dataset.find("TT") != -1 ) :
    cmd += " --args \"isTT=True\""
  print cmd
  os.system(cmd)
