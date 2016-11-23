#!/usr/bin/env python
import sys, os
from getopt import getopt

year = 16

baseAFS = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification"
baseURL = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification"
baseCAT = "%s/src/CATTools/CatProducer/data/LumiMask" % os.environ["CMSSW_BASE"]

configs = {
    16:{
        "lumiJSON":"Collisions16/13TeV/ReReco",
        "pileupJSON":"Collisions16/13TeV/PileUp/pileup_latest.txt",
        "minBiasXsec":69200.,
        "minBiasXsecUnc":0.046,
    },
    15:{
        "lumiJSON":"Collisions15/13TeV/Reprocessing",
        "pileupJSON":"Collisions15/13TeV/PileUp/pileup_latest.txt",
        "minBiasXsec":69000.,
        "minBiasXsecUnc":0.05,
    },
}

config = configs[year]

## Check JSON file in the configuration
import CATTools.CatProducer.catDefinitions_cfi as cat
jsonFiles = []
for p in dir(cat):
    if 'lumiJSON' not in p: continue
    jsonFiles.append(getattr(cat, p)+".txt")

## Check JSON file in the CATTools, copy them if it is new one
for jsonFile in jsonFiles:
    if os.path.exists("%s/%s" % (baseCAT, jsonFile)): continue
    subPath = "%s/%s" % (config["lumiJSON"], jsonFile)

    print "JSON file not in the CATTools. Checking the central JSON file repository."
    if os.path.exists(baseAFS):
        print "Copying JSON file to CATTools..."
        print "cp %s/%s %s/" % (baseAFS, subPath, baseCAT)
        os.system("cp %s/%s %s/" % (baseAFS, subPath, baseCAT))
    else:
        print "Downloading JSON file to CATTools..."
        print "wget %s/%s -O %s/%s" % (baseURL, subPath, baseCAT, jsonFile)
        os.system("wget --no-check-certificate %s/%s -O %s/%s" % (baseURL, subPath, baseCAT, jsonFile))

    if not os.path.exists("%s/%s" % (baseCAT, jsonFile)):
        print "!!! Failed to get JSON file ", jsonFile
        sys.exit(255)

puJSON = config["pileupJSON"]
puSubPath = os.path.dirname(puJSON)
puJSON = os.path.basename(puJSON)
if not os.path.exists("%s/%s" % (baseCAT, puJSON)):
    print "Pileup JSON not in the CATTools. Checking the central JSON file repository."
    if os.path.exists(baseAFS):
        print "Copying pileup JSON to CATTools..."
        print "cp %s/%s %s/" % (baseAFS, puSubPath, baseCAT)
        os.system("cp %s/%s %s/" % (baseAFS, puSubPath, baseCAT))
    else:
        print "Downloading pileup JSON to CATTools..."
        print "wget %s/%s/%s -O %s/%s" % (baseURL, puSubPath, puJSON, baseCAT, puJSON)
        os.system("wget --no-check-certificate %s/%s/%s -O %s/%s" % (baseURL, puSubPath, puJSON, baseCAT, puJSON))

    if not os.path.exists("%s/%s" % (baseCAT, puJSON)):
        print "!!! Failed to get Pileup JSON", puJSON
        sys.exit(255)

## Write pileup weight config file
#outfile = open('%s/src/CATTools/CatProducer/python/pileupWeight/pileupWeight20%d_cff.py' % (os.environ["CMSSW_BASE"], year), 'a')
outfile = open('pileupWeight.py', 'w')

import ROOT
for jsonFile in jsonFiles:
    puname = jsonFile[:-4]
    syst = ['', '_Up', '_Dn']
    for f in syst:
        minBiasXsec = config["minBiasXsec"]
        if   f == '_Up': minBiasXsec *= 1+config["minBiasXsecUnc"]
        elif f == '_Dn': minBiasXsec *= 1-config["minBiasXsecUnc"]

        fPU = "PileUpData%s.root" % f

        command  = "pileupCalc.py -i %s/%s " % (baseCAT, jsonFile)
        command += "--inputLumiJSON %s/%s " % (baseCAT, puJSON)
        command += "--minBiasXsec %i " % int(minBiasXsec)
        command += "--calcMode true --maxPileupBin 100 --numPileupBins 100 "
        command += fPU

        print command
        os.system(command)
        tt = ROOT.TFile(fPU)
        histo = tt.Get("pileup")
        outfile.write('#%s \n'%command.replace(os.environ["CMSSW_BASE"]+"/src",""))
        outfile.write('pileupWeightMap["%s%s"] = cms.vdouble(\n'%(puname,f))
        for b in range(1, histo.GetNbinsX()+1):
            outfile.write("%e,"%histo.GetBinContent(b))
            if b%5 == 0:
                outfile.write('\n')
        outfile.write(')\n')

outfile.close()

