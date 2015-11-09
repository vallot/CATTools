#!/usr/bin/env python
import ROOT,os

puname = 'Run2015_25nsSilver'
certJSON = 'Cert_246908-259891_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt'
#puname = 'Run2015_25ns'
#certJSON = 'Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
minBiasXsec = 69000.

syst = ['', 'Up', 'Dn']
for i, f in enumerate(syst):
    PileUpData = 'PileUpData%s.root'%(f)
    if i == 1: minBiasXsec = minBiasXsec*1.05
    if i == 2: minBiasXsec = minBiasXsec*0.95
    command = 'pileupCalc.py -i %s --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec %i --maxPileupBin 50 --numPileupBins 50 %s'%(certJSON,minBiasXsec,PileUpData)
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

