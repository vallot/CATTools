#!/usr/bin/env python
import ROOT,os,getopt,sys

lumiMask = None
minBiasXsec = 69000.

try:
    opts, args = getopt.getopt(sys.argv[1:],"hl:c:",["lumiMask",'minBiasXsec'])
except getopt.GetoptError:          
    print 'Usage : getPileUpData.py -l <lumiMask>'
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        print 'Usage : getPileUpData.py -l <lumiMask>'
        sys.exit()
    elif opt in ("-l", "--lumiMask"):
        lumiMask = arg
    elif opt in ("-c", "--minBiasXsec"):
        minBiasXsec = eval(arg)

puname = lumiMask[lumiMask.find('Cert_'):-4]
print puname
outfile = open('pileup.py', 'w')
outfile.write('import FWCore.ParameterSet.Config as cms\n')
outfile.write('pileupMap = {\n')
syst = ['', '_Up', '_Dn']
for i, f in enumerate(syst):
    PileUpData = 'PileUpData%s.root'%(f)
    if i == 1: minBiasXsec = minBiasXsec*1.05
    if i == 2: minBiasXsec = minBiasXsec*0.95
    command = 'pileupCalc.py -i %s --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec %i --maxPileupBin 50 --numPileupBins 50 %s'%(lumiMask,minBiasXsec,PileUpData)
    os.system(command)
    tt = ROOT.TFile(PileUpData)
    histo = tt.Get("pileup")
    outfile.write('#%s \n'%command)
    outfile.write('"%s%s":cms.vdouble('%(puname,f))
    for b in range(1, histo.GetNbinsX()+1):
        outfile.write("%e,"%histo.GetBinContent(b))
        if b%5 == 0:
            outfile.write('\n')
    outfile.write('),\n')      
outfile.write('}\n')
