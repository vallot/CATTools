#!/usr/bin/env python
import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os, getopt, sys
from CATTools.CatAnalyzer.histoHelper import *
from ROOT import TLorentzVector
#import DYestimation
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
datalumi = 2260
version = os.environ['CMSSW_VERSION']

rootfileDir = "/xrootd/store/user/pseudotop/ntuples/results_merged/%s/h2muAnalyzer_"%version
#rootfileDir = "/xrootd/store/user/pseudotop/ntuples/results_merged/v7-6-6/h2muAnalyzer_"
#rootfileDir = "root:///cms-xrdr.sdfarm.kr:1094//xrd/store/user/pseudotop/ntuples/results_merged/v7-6-6/h2muAnalyzer_"
#rootfileDir = "%s/src/CATTools/CatAnalyzer/test/results_merged/h2muAnalyzer_" % os.environ['CMSSW_BASE']
#rootfileDir = "%s/cattuples/20160324_163101/results_merged/h2muAnalyzer_" % os.environ['HOME_SCRATCH']

CMS_lumi.lumi_sqrtS = "%.0f pb^{-1}, #sqrt{s} = 13 TeV 25ns "%(datalumi)
mcfilelist = [
             # 'GG_HToMuMu',
             # 'GluGluToZZTo2mu2tau',
             # 'GluGluToZZTo2e2mu',
             # 'GluGluToZZTo4mu',
             # 'ttZToLLNuNu',
             # 'VBF_HToMuMu',
              'ZZTo4L_powheg',
              'ZZTo2L2Q',
              'ZZTo2L2Nu_powheg',
              'WWTo2L2Nu_powheg',
              'WZTo2L2Q',
              'WZTo3LNu_powheg',
              'TTJets_aMC',
              'DYJets',
              'DYJets_10to50',
             ]#ref : https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToMuMu
#mcfilelist = ['VBF_HToMuMu','WW','WZ','ZZ','TT_powheg','DYJets','DYJets_10to50']#,'WJets']
rdfilelist = [
              'MuonEG_Run2015',#emu
              'SingleElectron_Run2015',#ee
              'SingleMuon_Run2015',#mumu
              #'DoubleEG_Run2015', #compare emu to mc
              #'SingleMuon_Run2015C',
              #'SingleMuon_Run2015D'
             ]

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/dataset.json" % os.environ['CMSSW_BASE']))

#cut_step = "(step>=5)"
cut = 'dilep.M()>20&&step>=5&&filtered==1'
#cut = 'filtered==1&&%s&&%s'%(cut_step,emu_pid)
#cut = 'channel==2'
print cut
#weight = 'genweight*puweight*mueffweight*eleffweight*tri'
weight = 'weight'
plotvar = 'dilep.M()'
binning = [300, 0, 300]
x_name = 'mass [GeV]'
y_name = 'events'
dolog = False
f_name = 'll_m'

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


#tname = "cattree/nom"
ltname = ["/nom","/mu_u","/mu_d","/jes_u","/jes_d","/jer_u","/jer_d"]
ltcut = ["weight","(genweight)*(puweight_up)","(genweight)*(puweight_dn)"]
lhsum = []

dolog = True
#tcut = '(%s)*%s'%(cut,weight)
rdfname = rootfileDir + rdfilelist[2] +".root"

sig=[0,0,0,0,0,0]
bg=[0,0,0,0,0,0]
#lumilist= [datalumi,300*1000,900*1000,3000*1000] 

for itcut,tcut in enumerate(ltcut):
    tcut = '(%s)*%s'%(cut,tcut)
    print tcut
    for itname,tname in enumerate(ltname):
        mchistList = []
        tname = "cattree"+tname
        for imc,mcname in enumerate(mcfilelist):
            print tname
            data = findDataSet(mcname, datasets)
            scale = datalumi*data["xsec"]
            colour = data["colour"]
            title = data["title"]
            #if 'DYJets' in mcname: 
            #scale = scale*dyratio[channel][step] 
            #    scale = scale*dyratio[channel][1] 
            #if "HToMuMu" in mcname:
                #scale = scale*30.
                #title = title+" #times 30"
            rfname = rootfileDir + mcname +".root"
            print rfname

            #tfile = ROOT.TFile(rfname,"NET")
            tfile = ROOT.TFile(rfname,"READ")
            wentries = tfile.Get("cattree/nevents").Integral()
            scale = scale/wentries
 
            mchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)    
            mchist.SetFillColor(colour)
            mchist.SetLineColor(colour)
            mchistList.append(mchist)
        tname=""
        hsum=setLastHist(mchistList)
        hsum.SetLineColor(itname+1)
        lhsum.append(hsum)

k=0
rfile = ROOT.TFile("histo_background_%s.root"%(f_name),"RECREATE")
for i in range(len(ltcut)):
  for j in range(len(ltname)):
    lhsum[k].SetName(ltname[j]+"_"+ltcut[i])
    lhsum[k].Write()
    k+=1
rfile.Write()
rfile.Close()

