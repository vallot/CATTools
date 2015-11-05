import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os,copy
from CATTools.CatAnalyzer.histoHelper import *
ROOT.gROOT.SetBatch(True)

datalumi = 1280.23

mcfilelist = ['GG_HToMuMu','VBF_HToMuMu','WW','WZ','ZZ','TT_powheg','DYJets','DYJets_10to50']
rdfilelist = ['SingleMuon_Run2015']
rootfileDir = "/cms/scratch/jlee/v7-4-4/h2muAnalyzer_"
datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset.json" % os.environ['CMSSW_BASE']))

mchistList = []
cut = '(step>3&&isTight==1)*puweight'
y_name = 'events'
dolog = False
plot = 1
if plot == 1:
    plotvar = 'll_m'
    x_name = 'mass [GeV]'
    binning = [200, 0, 200]
    #binning = [ 0, 10,20,40,50]
    dolog = True
if plot == 2:
    plotvar = 'nvertex'
    x_name = 'no. vertex'
    binning = [30, 0, 30]

CMS_lumi.lumi_sqrtS = "%.0f pb^{-1}, #sqrt{s} = 13 TeV 25ns "%(datalumi)
tname = "cattree/nom"

for mcname in mcfilelist:
    data = findDataSet(mcname, datasets)
    scale = datalumi*data["xsec"]
    colour = data["colour"]
    title = data["title"]
    if "HToMuMu" in mcname:
        scale = scale*30.
        title = title+" #times 30"
    rootfilename = rootfileDir + mcname +".root"
    mchist = makeTH1(rootfilename, tname, title, binning, plotvar, cut, scale)    
    mchist.SetFillColor(colour)
    mchist.SetLineColor(colour)
    mchistList.append(mchist)

rootfilename = rootfileDir + rdfilelist[0] +".root"
rdhist = makeTH1(rootfilename, tname, 'data', binning, plotvar, cut)
if plot == 1:# blind data around higgs mass
    for i in range(11):
        rdhist.SetBinContent(120+i,-1)

drawTH1(plotvar+cut+".png", CMS_lumi, mchistList, rdhist, x_name, y_name,dolog)
