#!/usr/bin/env python
import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os, getopt, sys
from CATTools.CatAnalyzer.histoHelper import *
ROOT.gROOT.SetBatch(True)
'''
h2muDraw.py -c 'll_m>50&&step>=5&&isTight==1&&filtered==1' -b [40,0,40] -p nvertex -x 'no. vertex' &
h2muDraw.py -c 'll_m>50&&step>=5&&isTight==1&&filtered==1' -b [200,0,200] -p ll_m -x 'mass [GeV]' &
h2muDraw.py -c 'll_m>50&&step>=5&&isTight==1&&filtered==1' -b [100,0,100] -p ll_pt -x 'diMuon p_{T}  [GeV]' &
h2muDraw.py -c 'll_m>50&&step>=5&&isTight==1&&filtered==1' -b [100,0,100] -p lep1_pt -x 'leading muon p_{T}  [GeV]' &
h2muDraw.py -c 'll_m>50&&step>=5&&isTight==1&&filtered==1' -b [100,0,100] -p lep2_pt -x 'sub-leading muon p_{T}  [GeV]' &
h2muDraw.py -c 'll_m>50&&step>=5&&isTight==1&&filtered==1' -b [100,0,100] -p lep1_pt,lep2_pt -x 'muon p_{T}  [GeV]' &
h2muDraw.py -c 'll_m>50&&step>=5&&isTight==1&&filtered==1' -b [100,0,100] -p met -x 'met  [GeV]' &
h2muDraw.py -c 'll_m>50&&step>=5&&isTight==1&&filtered==1' -b [100,-3,3] -p lep1_eta,lep2_eta -x '#eta' &
'''

json_used = 'Golden'
datalumi = 2110
rootfileDir = "/cms/scratch/jlee/v7-4-6/h2muAnalyzer_"

CMS_lumi.lumi_sqrtS = "%.0f pb^{-1}, #sqrt{s} = 13 TeV 25ns "%(datalumi)
mcfilelist = ['VBF_HToMuMu','DYJets','ZZTo4L_powheg','ZZTo2L2Q','ZZTo2L2Nu_powheg','WWTo2L2Nu_powheg','WZTo2L2Q','WZTo3LNu_powheg','GluGluToZZTo2mu2tau','GluGluToZZTo2e2mu','GluGluToZZTo4mu','TTJets_aMC','ttZToLLNuNu']#ref : https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToMuMu
#mcfilelist = ['VBF_HToMuMu','WW','WZ','ZZ','TT_powheg','DYJets','DYJets_10to50']#,'WJets']
rdfilelist = ['SingleMuon_Run2015']

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/dataset.json" % os.environ['CMSSW_BASE']))

cut = 'll_m>50&&step>=5&&isTight==1&&filtered==1'
weight = 'weight'
plotvar = 'll_m'
binning = [200, 0, 200]
x_name = 'mass [GeV]'
y_name = 'events'
dolog = False
f_name = plotvar

try:
    opts, args = getopt.getopt(sys.argv[1:],"hdc:w:b:p:x:y:f:j:",["cut","weight","binning","plotvar","x_name","y_name","f_name","json_used","dolog"])
except getopt.GetoptError:          
    print 'Usage : ./h2muDraw.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -f <f_name> -j <json_used> -d <dolog>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./h2muDraw.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -f <f_name> -j <json_used> -d <dolog>'
        sys.exit()
    elif opt in ("-c", "--cut"):
        cut = arg
    elif opt in ("-w", "--weight"):
        weight = arg
    elif opt in ("-b", "--binning"):
        binning = eval(arg)
    elif opt in ("-p", "--plotvar"):
        plotvar = arg
    elif opt in ("-x", "--x_name"):
        x_name = arg
    elif opt in ("-y", "--y_name"):
        y_name = arg
    elif opt in ("-f", "--f_name"):
        f_name = arg
    elif opt in ("-j", "--json_used"):
        json_used = arg
    elif opt in ("-d", "--dolog"):
        dolog = True
print plotvar, x_name, f_name

if json_used=='Silver':
    datalumi = 2460
    rootfileDir = "/cms/scratch/jlee/v7-4-6/h2muAnalyzerSilver_"

    CMS_lumi.lumi_sqrtS = "%.0f pb^{-1}, #sqrt{s} = 13 TeV 25ns "%(datalumi)
    mcfilelist = ['GG_HToMuMu','DYJets','ZZTo4L_powheg','ZZTo2L2Q','ZZTo2L2Nu_powheg','WWTo2L2Nu_powheg','WZTo2L2Q','WZTo3LNu_powheg','GluGluToZZTo2mu2tau','GluGluToZZTo2e2mu','GluGluToZZTo4mu','TTJets_aMC','ttZToLLNuNu']#ref : https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToMuMu
    #mcfilelist = ['GG_HToMuMu','WW','WZ','ZZ','TT_powheg','DYJets','DYJets_10to50']#,'WJets'] silver
    rdfilelist = ['SingleMuon_Run2015']
    f_name = "Silver_"+f_name 

tname = "cattree/nom"
mchistList = []
tcut = '(%s)*%s'%(cut,weight)

for mcname in mcfilelist:
    data = findDataSet(mcname, datasets)
    scale = datalumi*data["xsec"]
    colour = data["colour"]
    title = data["title"]
    if "HToMuMu" in mcname:
        scale = scale*30.
        title = title+" #times 30"
    rfname = rootfileDir + mcname +".root"

    wentries = getWeightedEntries(rfname, tname, "tri",'weight')
    scale = scale/wentries
    
    mchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)    
    mchist.SetFillColor(colour)
    mchist.SetLineColor(colour)
    mchistList.append(mchist)

rfname = rootfileDir + rdfilelist[0] +".root"
rdhist = makeTH1(rfname, tname, 'data', binning, plotvar, tcut)
if plotvar == 'll_m':# blind data around higgs mass
    for i in range(11):
        rdhist.SetBinContent(120+i,0)

drawTH1(f_name+".png", CMS_lumi, mchistList, rdhist, x_name, y_name,dolog)
