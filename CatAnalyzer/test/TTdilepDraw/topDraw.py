#!/usr/bin/env python
import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os, getopt, sys
from CATTools.CatAnalyzer.histoHelper import *
ROOT.gROOT.SetBatch(True)
'''
topDraw.py -a 1 -s 1 -c 'tri==1&&filtered==1' -b [40,0,40] -p nvertex -x 'no. vertex' &
topDraw.py -a 1 -s 1 -b [100,-3,3] -p lep1_eta,lep2_eta -x '#eta' &
'''
datalumi = 1.56
CMS_lumi.lumi_sqrtS = "%.2f fb^{-1}, #sqrt{s} = 13 TeV "%(datalumi)
datalumi = datalumi*1000 # due to fb

mcfilelist = ['WW','WZ','ZZ','TT_powheg','DYJets','DYJets_10to50']#,'WJets']
rdfilelist = ['MuonEG_Run2015','DoubleEG_Run2015','DoubleMuon_Run2015']
rootfileDir = "/cms/scratch/jlee/v7-4-5/TtbarDiLeptonAnalyzer_"

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset.json" % os.environ['CMSSW_BASE']))

cut = 'tri==1&&filtered==1'
weight = 'weight'
plotvar = 'll_m'
binning = [200, 0, 200]
x_name = 'mass [GeV]'
y_name = 'events'
dolog = False
channel = 1
step = 1
try:
    opts, args = getopt.getopt(sys.argv[1:],"hdc:w:b:p:x:y:a:s:",["cut","weight","binning","plotvar","x_name","y_name","dolog","channel","step"])
except getopt.GetoptError:          
    print 'Usage : ./topDraw.py.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./topDraw.py.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog>'
        sys.exit()
    elif opt in ("-c", "--cut"):
        cut = arg
    elif opt in ("-a", "--channel"):
        channel = int(arg)
    elif opt in ("-s", "--step"):
        step = int(arg)
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
    elif opt in ("-d", "--dolog"):
        dolog = True

tname = "cattree/nom"
mchistList = []

if channel == 1: ttother_tcut = "!(parton_channel==2 && ((parton_mode1==1 && parton_mode2==2) || (parton_mode1==2 && parton_mode2==1)))"
elif channel == 2: ttother_tcut = "!(parton_channel==2 && (parton_mode1==2 && parton_mode2==2))"
elif channel == 3: ttother_tcut = "!(parton_channel==2 && (parton_mode1==1 && parton_mode2==1))"

stepch_tcut =  'step>=%i&&channel==%i'%(step,channel)
tcut = '(%s&&%s)*%s'%(stepch_tcut,cut,weight)
ttother_tcut = '(%s&&%s&&%s)*%s'%(stepch_tcut,cut,ttother_tcut,weight)
print "TCut =",tcut

#DY estimation
dyratio = [[0 for x in range(6)] for x in range(4)]
dyratio[1][step] = 1.
if channel !=1:
    scale = 1.
    dycut = ""
    if step == 1: dycut = "(step1==1)*"
    if step == 2: dycut = "(step1==1)*"
    if step == 3: dycut = "(step1==1)*(step3==1)*"
    if step == 4: dycut = "(step1==1)*(step3==1)*(step4==1)*"
    if step == 5: dycut = "(step1==1)*(step3==1)*(step4==1)*(step5==1)*"

    rfname = rootfileDir + 'DYJets' +".root"
    data = findDataSet('DYJets', datasets)
    scale = datalumi*data["xsec"]
    wentries = getWeightedEntries(rfname, tname, "tri",weight)
    scale = scale/wentries
    
    mc_ee_in = makeTH1(rfname,tname,"mc_ee_in", binning, plotvar, dycut+'(%s && channel==2 && step2==0)*(%s)'%(cut,weight), scale)
    mc_mm_in = makeTH1(rfname,tname,"mc_mm_in", binning, plotvar, dycut+'(%s && channel==3 && step2==0)*(%s)'%(cut,weight), scale)
    mc_ee_out = makeTH1(rfname,tname,"mc_ee_in", binning, plotvar, dycut+'(%s && channel==2 && step2==1)*(%s)'%(cut,weight), scale)
    mc_mm_out = makeTH1(rfname,tname,"mc_mm_in", binning, plotvar, dycut+'(%s && channel==3 && step2==1)*(%s)'%(cut,weight), scale)

    rfname = rootfileDir + 'DYJets_10to50' +".root"
    data = findDataSet('DYJets_10to50', datasets)
    scale = datalumi*data["xsec"]
    wentries = getWeightedEntries(rfname, tname, "tri",weight)
    scale = scale/wentries
    mc_ee_in.Add(makeTH1(rfname,tname,"mc_ee_in", binning, plotvar, dycut+'(%s && channel==2 && step2==0)*(%s)'%(cut,weight), scale))
    mc_mm_in.Add(makeTH1(rfname,tname,"mc_mm_in", binning, plotvar, dycut+'(%s && channel==3 && step2==0)*(%s)'%(cut,weight), scale))
    mc_ee_out.Add(makeTH1(rfname,tname,"mc_ee_in", binning, plotvar, dycut+'(%s && channel==2 && step2==1)*(%s)'%(cut,weight), scale))
    mc_mm_out.Add(makeTH1(rfname,tname,"mc_mm_in", binning, plotvar, dycut+'(%s && channel==3 && step2==1)*(%s)'%(cut,weight), scale))
    
    rfname = rootfileDir+rdfilelist[1-1]+".root"
    rd_em_in = makeTH1(rfname, tname,'rd_em_in', binning, plotvar, dycut+'(%s && channel==1 && ((ll_m > 76) && (ll_m < 106)))'%(cut))
    rfname = rootfileDir + rdfilelist[2-1] +".root"
    rd_ee_in = makeTH1(rfname, tname,'rd_ee_in', binning, plotvar, dycut+'(%s && channel==2 && step2 ==0)'%(cut))
    rfname = rootfileDir + rdfilelist[3-1] +".root"
    rd_mm_in = makeTH1(rfname, tname,'rd_mm_in', binning, plotvar, dycut+'(%s && channel==3 && step2 ==0)'%(cut))

    dyest = drellYanEstimation(mc_ee_in.Integral(), mc_ee_out.Integral(), mc_mm_in.Integral(), mc_mm_out.Integral(),
                               rd_ee_in.Integral(), rd_mm_in.Integral(), rd_em_in.Integral())
    print "DY estimation for", step, "ee =",dyest[0], "mm =",dyest[1]   
    dyratio[2][step] = dyest[0]
    dyratio[3][step] = dyest[1]


for i, mcname in enumerate(mcfilelist):
    data = findDataSet(mcname, datasets)
    scale = datalumi*data["xsec"]
    colour = data["colour"]
    title = data["title"]
    if 'DYJets' in mcname:
        scale = scale*dyratio[channel][step]

    rfname = rootfileDir + mcname +".root"
    wentries = getWeightedEntries(rfname, tname, "tri",weight)
    scale = scale/wentries
    
    mchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)
    mchist.SetLineColor(colour)
    mchist.SetFillColor(colour)
    mchistList.append(mchist)
    if 'TT_powheg' == mcname:
        ttothers = makeTH1(rfname, tname, title+' others', binning, plotvar, ttother_tcut, scale)
        ttothers.SetLineColor(906)
        ttothers.SetFillColor(906)
        mchistList.append(ttothers)
        mchist.Add(ttothers, -1)

rfname = rootfileDir + rdfilelist[channel-1] +".root"
rdhist = makeTH1(rfname, tname, 'data', binning, plotvar, tcut)
rdhist.SetLineColor(1)

drawTH1(plotvar+tcut+".png", CMS_lumi, mchistList, rdhist, x_name, y_name, True)
