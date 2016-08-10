#!/usr/bin/env python
import sys, os
from getopt import getopt

"""try:
    opts, args = getopt(sys.argv[1:],"hl:y:",["lumiMask",'year'])
except getopt.GetoptError:
    print 'Usage : getPileUpData.py -l <lumiMask>'
    sys.exit(2)
"""

year = 16

configs = {
    16:{
        "lumiJSON":"/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV",
        "pileupJSON":"/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt",
        "minBiasXsec":69200.,
        "minBiasXsecUnc":0.046,
    },
    15:{
        "lumiJSON":"/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Reprocessing",
        "pileupJSON":"/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt",
        "minBiasXsec":69000.,
        "minBiasXsecUnc":0.05,
    },
}

config = configs[year]
jsonBasedir = config["lumiJSON"]
jsonBasedirCAT = "%s/src/CATTools/CatProducer/data/LumiMask" % os.environ["CMSSW_BASE"]

## Check JSON file in the configuration
import CATTools.CatProducer.catDefinitions_cfi as catDefinitions
jsonFiles = []
for p in dir(catDefinitions):
    if 'lumiJSON' not in p: continue
    jsonFiles.append(getattr(catDefinitions, p)+".txt")

## Check JSON file in the CATTools, copy them if it is new one
for jsonFile in jsonFiles:
    if os.path.exists("%s/%s" % (jsonBasedirCAT, jsonFile)): continue

    print "JSON file not in the CATTools. Checking the central JSON file repository."
    if not os.path.exists("%s/%s" % (jsonBasedir, jsonFile)):
        print "!!! This file is not in the central JSON file repository, neither !!!"
        sys.exit(1)

    print "Copying JSON file to CATTools..."
    print "cp %s/%s %s/" % (jsonBasedir, jsonFile, jsonBasedirCAT)
    os.system("cp %s/%s %s/" % (jsonBasedir, jsonFile, jsonBasedirCAT))

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

        command  = "pileupCalc.py -i %s/%s " % (jsonBasedirCAT, jsonFile)
        command += "--inputLumiJSON %s " % config["pileupJSON"]
        command += "--minBiasXsec %i " % int(minBiasXsec)
        command += "--calcMode true --maxPileupBin 50 --numPileupBins 50 "
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

