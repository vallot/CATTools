#!/usr/bin/env python
import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os, getopt, sys
from CATTools.CatAnalyzer.histoHelper import *
import DYestimation
ROOT.gROOT.SetBatch(True)
'''
topDraw.py -a 1 -s 1 -c 'tri==1&&filtered==1' -b [40,0,40] -p nvertex -x 'no. vertex' &
topDraw.py -a 1 -s 1 -b [100,-3,3] -p lep1_eta,lep2_eta -x '#eta' &
'''
datalumi = 2.11
CMS_lumi.lumi_sqrtS = "%.2f fb^{-1}, #sqrt{s} = 13 TeV"%(datalumi)
datalumi = datalumi*1000 # due to fb

mcfilelist = ['TT_powheg', 'WJets', 'SingleTbar_tW', 'SingleTop_tW', 'ZZ', 'WW', 'WZ', 'DYJets', 'DYJets_10to50']
rdfilelist = ['MuonEG_Run2015','DoubleEG_Run2015','DoubleMuon_Run2015']
rootfileDir = "file:/xrootd/store/user/tt8888tt/v7-4-6/"
channel_name = ['MuEl', 'ElEl', 'MuMu']

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset.json" % os.environ['CMSSW_BASE']))

#defalts
step = 1
channel = 1
cut = 'tri==1&&filtered==1'
weight = 'genweight*puweight*lepweight'
binning = [60, 20, 320]
plotvar = 'll_m'
x_name = 'mass [GeV]'
y_name = 'Events'
dolog = False

#get input
try:
    opts, args = getopt.getopt(sys.argv[1:],"hdc:w:b:p:x:y:a:s:",["cut","weight","binning","plotvar","x_name","y_name","dolog","channel","step"])
except getopt.GetoptError:          
    print 'Usage : ./topDraw.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./topDraw.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog>'
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

#cut define
if channel == 1: ttother_tcut = "!(parton_channel==2 && ((parton_mode1==1 && parton_mode2==2) || (parton_mode1==2 && parton_mode2==1)))"
elif channel == 2: ttother_tcut = "!(parton_channel==2 && (parton_mode1==2 && parton_mode2==2))"
elif channel == 3: ttother_tcut = "!(parton_channel==2 && (parton_mode1==1 && parton_mode2==1))"
stepch_tcut =  'step>=%i&&channel==%i'%(step,channel)

tcut = '(%s&&%s)*(%s)'%(stepch_tcut,cut,weight)
ttother_tcut = '(%s&&%s&&%s)*(%s)'%(stepch_tcut,cut,ttother_tcut,weight)
print "TCut =",tcut

#namming
x_name = channel_name[channel-1]+" "+x_name
if len(binning) <= 3:
    num = (binning[2]-binning[1])/float(binning[0])
    if num != 1:
        if x_name.endswith(']'):
            unit = "["+x_name.split('[')[1]
        else: unit = ""
        y_name = y_name + "/%g%s"%(num,unit)

#DYestimation
if not os.path.exists('./DYFactor.json'):
	DYestimation.printDYFactor(rootfileDir, tname, datasets, datalumi, cut, weight, rdfilelist)# <------ This will create 'DYFactor.json' on same dir.
dyratio=json.load(open('./DYFactor.json'))

#saving histos
mchistList = []
for i, mcname in enumerate(mcfilelist):
	data = findDataSet(mcname, datasets)
	scale = datalumi*data["xsec"]
	colour = data["colour"]
	title = data["title"]
	if 'DYJets' in mcname:
		scale = scale*dyratio[channel][step]

	rfname = rootfileDir + mcname +".root"

	wentries = getWeightedEntries(rfname, tname, "tri", weight)
	scale = scale/wentries
		
	mchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)
	mchist.SetLineColor(colour)
	mchist.SetFillColor(colour)
	mchistList.append(mchist)
	if 'TT' in mcname:
		ttothers = makeTH1(rfname, tname, title+' others', binning, plotvar, ttother_tcut, scale)
		ttothers.SetLineColor(906)
		ttothers.SetFillColor(906)
		mchistList.append(ttothers)
		mchist.Add(ttothers, -1)

rdtcut = '%s&&%s'%(stepch_tcut,cut)
rfname = rootfileDir + rdfilelist[channel-1] +".root"
rdhist = makeTH1(rfname, tname, 'data', binning, plotvar, rdtcut)
rdhist.SetLineColor(1)

#bin normalize
for hist in mchistList:
	for i in range(len(binning)):
		hist.SetBinContent(i, hist.GetBinContent(i)/hist.GetBinWidth(i))
		hist.SetBinError(i, hist.GetBinError(i)/hist.GetBinWidth(i))
for i in range(len(binning)):
	rdhist.SetBinContent(i, rdhist.GetBinContent(i)/rdhist.GetBinWidth(i))
	rdhist.SetBinError(i, rdhist.GetBinError(i)/rdhist.GetBinWidth(i))

#Drawing plots on canvas
var = plotvar.split(',')[0]
var = ''.join(i for i in var if not i.isdigit())
outfile = "%s_s%d_%s.png"%(channel_name[channel-1],step,var)
drawTH1(outfile, CMS_lumi, mchistList, rdhist, x_name, y_name, dolog)
print outfile

