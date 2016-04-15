#!/usr/bin/env python
import sys,os
os.system("hostname")
os.system("source /cvmfs/cms.cern.ch/cmsset_default.sh")
os.system("scramv1 project CMSSW CMSSW_7_6_3_patch2")
os.system("cd CMSSW_7_6_3_patch2")
os.system("eval `scramv1 runtime -sh`")
os.system("cd ..")

cmd = "root -l -b -q \"opti_d0.C(%f,%f)\" "%(float(sys.argv[1]),float(sys.argv[2]))
print cmd
os.system(cmd)
