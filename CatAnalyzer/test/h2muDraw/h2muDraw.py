import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os,copy
from CATTools.CatAnalyzer.histoHelper import *
ROOT.gROOT.SetBatch(True)

datalumi = 1280.23

mcfilelist = ['TT_powheg','DYJets','DYJets_10to50']
rdfilelist = ['SingleMuon_Run2015']
rootfileDir = "/cms/scratch/jlee/v7-4-4/h2muAnalyzer_"
datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset.json" % os.environ['CMSSW_BASE']))

mchistList = []
binning = [30, 0, 30]
#plotvar = 'll_m'
plotvar = 'nvertex'
cut = '(step >1)*puweight'
x_name = 'mass [GeV]'
y_name = 'events'
CMS_lumi.lumi_sqrtS = "%.0f pb^{-1}, #sqrt{s} = 13 TeV 25ns "%(datalumi)
tname = "cattree/nom"

for i, mcname in enumerate(mcfilelist):
    data = findDataSet(mcname, datasets)
    scale = datalumi*data["xsec"]
    colour = data["colour"]
    title = data["title"]
    
    rootfilename = rootfileDir + mcname +".root"
    mchist = makeTH1(rootfilename, tname, title, binning, plotvar, cut, scale)    
    mchist.SetFillColor(colour)
    mchist.SetLineColor(colour)
    mchistList.append(mchist)

rootfilename = rootfileDir + rdfilelist[0] +".root"
rdhist = makeTH1(rootfilename, tname, 'data', binning, plotvar, cut)

drawTH1(plotvar+".png", CMS_lumi, mchistList, rdhist, x_name, y_name)
