#!/usr/bin/env python

from ROOT import *
import sys, os
import json
from multiprocessing import Pool, cpu_count

sys.path.append("%s/src/CATTools/CatAnalyzer/python" % os.environ['CMSSW_BASE'])
from HistMaker import *

basePath = 'root://cms-xrdr.sdfarm.kr///xrd/store/user/youngjo/Cattools/v736'
outPath = "hist"

def makeHist(dset):
    timer = TStopwatch()
    timer.Start()

    name = dset['name']
    if 'nevt' in dset: nevt = dset['nevt']
    else: nevt = 0

    an = HistMaker(
        inFileNames = ("%s/ntuple_%s.root" % (basePath, name)),
        outFileName = ("%s/hist_%s.root" % (outPath, name)),
        modName = "ntuple", treeName = "event",
        #eventCount = "hNEvent",
        eventCount = nevt,
    )
    if not an.isOK: return

    an.addH1("nvtx", "nGoodPV", "Vertex multiplicity;Vertex multiplicity;Events", 100, 0, 100)

    an.addH1("l1_pt", "muons_pt[0]", "Leading lepton p_{T};p_{T}^{l, 1st lead} (GeV);Events/20 GeV", 100, 0, 200)
    an.addH1("l2_pt", "muons_pt[1]", "Trailing lepton p_{T};p_{T}^{l, 2nd lead} (GeV);Events/20 GeV", 100, 0, 200)
    an.addH1("l1_eta", "muons_eta[0]", "Leading lepton #eta;#eta^{l, 1st lead};Events/0.1", 60, -3, 3)
    an.addH1("l2_eta", "muons_eta[1]", "Trailing lepton #eta;#eta^{l, 2nd lead};Events/0.1", 60, -3, 3)

    an.addH1("j1_pt", "jets_pt[0]", "Leading jet p_{T};p_{T}^{j, 1st} (GeV);Events/20 GeV", 100, 0, 200)
    an.addH1("j2_pt", "jets_pt[1]", "Trailing jet p_{T};p_{T}^{j, 2nd} (GeV);Events/20 GeV", 100, 0, 200)
    an.addH1("j3_pt", "jets_pt[2]", "3rd jet p_{T};p_{T}^{j, 3rd} (GeV);Events/20 GeV", 100, 0, 200)
    an.addH1("j4_pt", "jets_pt[3]", "4th jet p_{T};p_{T}^{j, 4th} (GeV);Events/20 GeV", 100, 0, 200)
    an.addH1("j1_eta", "jets_eta[0]", "Leading jet #eta;#eta^{j, 1st};Events/0.1", 60, -3, 3)
    an.addH1("j2_eta", "jets_eta[1]", "Trailing jet #eta;#eta^{j, 2nd};Events/0.1", 60, -3, 3)
    an.addH1("j3_eta", "jets_eta[2]", "3rd jet #eta;#eta^{j, 3rd};Events/0.1", 60, -3, 3)
    an.addH1("j4_eta", "jets_eta[3]", "4th jet #eta;#eta^{j, 4th};Events/0.1", 60, -3, 3)

    an.addCutStep("S1", "", "l1_pt,l2_pt,l1_eta,l2_eta,nvtx")
    #an.addCutStep("S2", "abs(mLL-91.2) > 15", "l1_pt,l2_pt,l1_eta,l2_eta")
    an.addCutStep("S3", "@jets_pt.size() >= 2", "j1_pt,j2_pt,j1_eta,j2_eta")

    an.process()

    timer.Stop()
    return (name, timer.CpuTime(), timer.RealTime())

if __name__ == '__main__':
    if not os.path.exists(outPath): os.makedirs(outPath)
    datasets = json.load(open("%s/src/CatAnalyzer/data/dataset/dataset.json" % os.environ['CMSSW_BASE']))

    nCPU = cpu_count()
    p = Pool(nCPU)
    print "="*40
    print "Start processing using %d CPUs" % nCPU
    print "="*40
    timer = TStopwatch()
    timer.Start()
    cpuTimes = p.map(makeHist, datasets)
    timer.Stop()
    print "="*40
    print "End of processing"
    print "."*40
    print "Total real time = %.1f s" % timer.RealTime()
    print "Total CPU time  = %.1f s" % sum(x[1] for x in cpuTimes)
    print "-"*40
    for l in cpuTimes:
        print "%s: CPU=%.1f s, Real=%.1f s" % (l[0], l[1], l[2])
    print "="*40

