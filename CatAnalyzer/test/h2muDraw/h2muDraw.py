import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os,copy
from CATTools.CatAnalyzer.histoHelper import *
ROOT.gROOT.SetBatch(True)

datalumi = 1300

mcfilelist = ['TT_powheg','DYJets']
rdfilelist = ['SingleMuon_Run2015']
rootfileDir = "/afs/cern.ch/work/j/jlee/h2muAnalyzer_"
datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset.json" % os.environ['CMSSW_BASE']))

mchistList = []
binning = [30, 0, 30]
#plotvar = 'll_m'
plotvar = 'nvertex'
cut = '(step >1)'
x_name = 'mass [GeV]'
y_name = 'events'
CMS_lumi.lumi_sqrtS = "13 TeV, 1.3fb"

for i, mcname in enumerate(mcfilelist):
    rootfilename = rootfileDir + mcname +".root"
    tt = ROOT.TFile(rootfilename)
    tree = tt.h2mu.Get("tree")
    mccut = cut+'*(puweight)'
    scale = 1.
    for data in datasets:
        if data["name"] == mcname:
            scale = datalumi*data["xsec"]
    scale = scale/tree.GetEntries()
    mchist = getTH1(mcname, binning, tree, plotvar, mccut)    
    mchist.Scale(scale)
    mchist.SetFillColor(i+2)
    mchistList.append(mchist)


rootfilename = rootfileDir + rdfilelist[0] +".root"
tt = ROOT.TFile(rootfilename)
tree = tt.h2mu.Get("tree")
rdhist = getTH1('data', binning, tree, plotvar, cut)
rdhist.SetLineColor(1)

drawTH1withRatio(plotvar, CMS_lumi, mchistList, rdhist, x_name, y_name)
