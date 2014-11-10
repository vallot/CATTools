#!/usr/bin/env python

import os,sys

if len(sys.argv) < 2 :
  print "Usage : ./createMulticrab.py ttbar_mc.txt ttbar_rd.txt ..."
  exit(-1)

output = open("multicrab.cfg","w")

datasets =[]
for file in sys.argv[1:] :
  lines = open(file)
  datasets += lines.readlines()

#MC or Data?
isMC = True
type_label = "MC"

dataset = datasets[0]
datatype = dataset.strip().split("/")[-1]
if datatype == "AOD" :
  isMC = False
  type_label = "RD"



#print datasets
#Init multicrab.cfg
header = """[MULTICRAB]\ncfg=crab%s.cfg\n[COMMON]\n"""%(type_label)

output.write(header)

for dataset in datasets :
  dataset = dataset.strip()
  if isMC :
    label = dataset.split("/")[1]
  else :
    label = dataset.split("/")[1]+"_"+dataset.split("/")[2]
  output.write("[%s]\n"%(label))
  output.write("CMSSW.datasetpath = %s\n"%(dataset))
