import ROOT,os
"""
to get the RD pileup histo
pileupCalc.py -i Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt --calcMode true  --minBiasXsec 80000 --maxPileupBin 80 --numPileupBins 80 PileUpData.root
pileupCalc.py -i Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt --calcMode true  --minBiasXsec 88000 --maxPileupBin 80 --numPileupBins 80 PileUpData_up.root
pileupCalc.py -i Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt --calcMode true  --minBiasXsec 72000 --maxPileupBin 80 --numPileupBins 80 PileUpData_down.root
"""
pileupFiles = ["PileUpData.root", "PileUpData_up.root", "PileUpData_down.root"]
for i, f in enumerate(pileupFiles):
    tt = ROOT.TFile(f)
    histo = tt.Get("pileup")
    print f
    for b in range(histo.GetNbinsX()):
        if b == 0:
            continue        
        print ("%e,"%histo.GetBinContent(b)),
        if b%5 == 0:
            print
    print 

