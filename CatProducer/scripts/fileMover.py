#!/usr/bin/env python
# Author: Tae.Jeong.Kim@cern.ch
# Filename: analysis.py
# How to run: ./copy.py status mc ElEl
import os
import re
import sys
import time
import commands
#os.system("voms-proxy-init -voms cms")
pfnBaseURLs = {}
pfnBaseURLs['T2_BE_IIHE'] = 'srm://maite.iihe.ac.be:8443/srm/managerv2?SFN=/pnfs/iihe/cms/ph/sc4'
pfnBaseURLs['T2_UK_London_IC'] = 'srm://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms'
pfnBaseURLs['T2_IT_Roma'] = 'srm://cmsrm-se01.roma1.infn.it/pnfs/roma1.infn.it/data/cms'
pfnBaseURLs['T2_US_Florida'] = 'srm://srmb.ihepa.ufl.edu/cms/data'
pfnBaseURLs['T2_DE_DESY'] = 'srm://dcache-se-cms.desy.de:8443/srm/managerv2?SFN=/pnfs/desy.de/cms/tier2'
pfnBaseURLs['T2_CH_CERN'] = 'srm://srm-eoscms.cern.ch:8443/srm/v2/server?SFN=/eos/cms/' 


location = pfnBaseURLs['T2_CH_CERN']
dir = '/store/relval/CMSSW_7_0_7/RelValTTbar_13/GEN-SIM-RECO/PU25ns_PLS170_V7AN1-v1/00000/'

dest="/cms/home/tjkim/"+dir
os.system("mkdir -p "+dest)

cmd="lcg-ls --verbose -D  srmv2 --nobdii "+location+dir

fileNames = commands.getoutput(cmd).split("\n")

nfile = 1
i=0

for fileName in fileNames:
    if i==nfile:
      break

    if fileName.endswith("root") > 0:
      print "copying..." + fileName
      files = fileName.split("/")
      destfile = files[len(files)-1]
      os.system("lcg-cp -b -D srmv2 --vo cms "+location+dir+destfile+" "+dest+"/"+destfile)
      i=i+1
  





