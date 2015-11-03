import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os,copy
from CATTools.CatAnalyzer.histoHelper import *
ROOT.gROOT.SetBatch(True)

datalumi = 1280.23

mcfilelist = ['TT_powheg','DYJets']
rdfilelist = ['MuonEG','DoubleEG','DoubleMuon']
rootfileDir = "/cms/scratch/tt8888tt/cattools_v744/src/CATTools/CatAnalyzer/test/result3/"
datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset.json" % os.environ['CMSSW_BASE']))

mchistList = []
#plotvar = 'nvertex'
#binning = [30, 0, 30]
plotvar = 'll_m'
binning = [50, 0, 200]
channel = 2
cut = '(step>=1 && channel == %i && filtered == 1)'%(channel)
print "TCut =",cut
x_name = 'mass [GeV]'
y_name = 'events'
CMS_lumi.lumi_sqrtS = "%.0f pb^{-1}, #sqrt{s} = 13 TeV 25ns "%(datalumi)

for i, mcname in enumerate(mcfilelist):
    rootfilename = rootfileDir + mcname +".root"
    tt = ROOT.TFile(rootfilename)
    tree = tt.ttll.Get("tree")
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


rootfilename = rootfileDir + rdfilelist[channel-1] +".root"
tt = ROOT.TFile(rootfilename)
tree = tt.ttll.Get("tree")
rdhist = getTH1('data', binning, tree, plotvar, cut)
rdhist.SetLineColor(1)

drawTH1(plotvar+".png", CMS_lumi, mchistList, rdhist, x_name, y_name)
