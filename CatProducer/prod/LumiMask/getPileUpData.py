#!/usr/bin/env python

certURL = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV"
puname = 'Run2015_25nsSilver'
certJSON = 'Cert_246908-259891_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt'
#puname = 'Run2015_25ns'
#certJSON = 'Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
minBiasXsec = 69000.

import os
from urllib import urlretrieve
if not os.path.exists(certJSON):
    print "Downloading Lumi JSON file..."
    print "  "+certURL+"/"+certJSON
    urlretrieve(certURL+"/"+certJSON, certJSON)

if not os.path.exists("pileup_latest.txt"):
    print "Downloading Pileup JSON file..."
    print "  "+certURL+"/PileUp/pileup_latest.txt"
    urlretrieve(certURL+"/PileUp/pileup_latest.txt", 'pileup_latest.txt')

import ROOT
syst = ['', 'Up', 'Dn']
for i, f in enumerate(syst):
    PileUpData = 'PileUpData%s.root'%(f)
    if i == 1: minBiasXsec = minBiasXsec*1.05
    if i == 2: minBiasXsec = minBiasXsec*0.95
    command = 'pileupCalc.py -i %s --inputLumiJSON pileup_latest.txt --calcMode true --minBiasXsec %i --maxPileupBin 50 --numPileupBins 50 %s'%(certJSON,minBiasXsec,PileUpData)
    os.system(command)
    tt = ROOT.TFile(PileUpData)
    histo = tt.Get("pileup")
    print '#%s'%command
    print '"%s%s":cms.vdouble('%(puname,f)
    for b in range(1, histo.GetNbinsX()+1):
        print ("%e,"%histo.GetBinContent(b)),
        if b%5 == 0:
            print
    print '),'

