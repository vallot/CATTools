#!/usr/bin/env python
# Author: Tae.Jeong.Kim@cern.ch
import os
import re
import sys
import time
import commands

os.system("voms-proxy-init -voms cms")

location = "root://xrootd.unl.edu/" 
fileNames = [
'/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/003B199E-0F81-E411-8E76-0025905A60B0.root',
'/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/004ED5E6-4E7F-E411-88E7-0025905A606A.root',
'/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/00B13EC1-4D7F-E411-B9E5-0026189437E8.root',
'/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/02131EBA-627E-E411-B166-0025905A605E.root',
'/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/02A66D37-C07F-E411-9291-0025905A609A.root',
'/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/02C100E7-4E7F-E411-A6F1-0025905A60A6.root',
'/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/04809019-F37F-E411-87D1-0025905AA9F0.root',
'/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/048FB111-F37F-E411-A2D5-00261894380B.root',
'/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/0684E1E4-4E7F-E411-8909-00261894392F.root',
'/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/069C95E6-4E7F-E411-AA3F-002618943869.root'
]


i = 0
for file in fileNames:
  prefix = file.split("/")
  outname = prefix[8]
  destination = "/eos/cms/store/caf/user/tjkim/" + prefix[2] + "/" + prefix[3] + "/" + prefix[4] + "/" + prefix[5] + "/" + prefix[6] + "/" + prefix[7] + "/"
  #if i == 0:
    #currently eos command does not work in python. so the directory needs to be created manually
    #os.system("eos mkdir -p "+destination) 
    #print "created data path" 
  print "copying... " + outname
  print "to... " + destination
 
  cmd = "xrdcp "+location+file+" root://eoscms/" + destination
  os.system(cmd)
  i = i + 1;
