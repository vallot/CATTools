#!/usr/bin/env python
import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os, getopt, sys
from CATTools.CatAnalyzer.histoHelper import *
import DYestimation
ROOT.gROOT.SetBatch(True)
'''
dstarDraw.py -a 1 -s 1 -c 'tri==1&&filtered==1' -b [40,0,40] -p nvertex -x 'no. vertex' &
dstarDraw.py -a 1 -s 1 -b [100,-3,3] -p lep1_eta,lep2_eta -x '#eta' &
'''
datalumi = 15.92 # Run2016 B & C & D & E, v8-0-1
CMS_lumi.lumi_sqrtS = "%.1f fb^{-1}, #sqrt{s} = 13 TeV"%(datalumi)
datalumi = datalumi*1000 # due to fb

mcfilelist = ['TT_powheg', 'WJets', 'SingleTbar_tW', 'SingleTop_tW', 'ZZ', 'WW', 'WZ', 'DYJets', 'DYJets_10to50']
#rdfilelist = ['MuonEG_Run2016','DoubleEG_Run2016','DoubleMuon_Run2016']
#rootfileDir = "/cms/scratch/geonmo/for2016KPS_Ana/src/CATTools/CatAnalyzer/test/cattools/cattree_"
rootfileDir = "/xrootd/store/user/quark2930/dilepton_mass_v801_16092901/cattree_"
#channel_name = ['Combined', 'MuEl', 'ElEl', 'MuMu']

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/dataset.json" % os.environ['CMSSW_BASE']))

#defalts
step = 1
channel = 3
cut = 'tri!=0&&filtered==1&&'
acut = 'tri!=0&&filtered==1&&'
weight = 'genweight*puweight*mueffweight*eleffweight*tri'
binning = [60, 20, 320]
plotvar = 'll_m'
x_name = 'mass [GeV]'
y_name = 'Events'
dolog = False
overflow = False
binNormalize = False
suffix = ''
#get input
try:
    opts, args = getopt.getopt(sys.argv[1:],"hdnoc:b:p:x:y:a:s:f:",["binNormalize","overflow","cut","binning","plotvar","x_name","y_name","dolog","acut","step","suffix"])
except getopt.GetoptError:          
    print 'Usage : ./compare.py -c <cut> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog> -f <suffix>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./topDraw.py -c <cut> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog> -f <suffix>'
        sys.exit()
    elif opt in ("-c", "--cut"):
        cut += arg
    elif opt in ("-a", "--acut"):
        acut += arg
    elif opt in ("-x", "--x_name"):
        x_name = arg
    elif opt in ("-y", "--y_name"):
        y_name = arg
    elif opt in ("-d", "--dolog"):
        dolog = True
    elif opt in ("-p", "--plotvar"):
        plotvar = arg
    elif opt in ("-o", "--overflow"):
        overflow = True
    elif opt in ("-s", "--step"):
        step = int(arg)
    elif opt in ("-n", "--binNormalize"):
        binNormalize = True
    elif opt in ("-f", "--suffix"):
        suffix = "_"+arg

tname = "cattree/nom"

#cut define
if   channel == 1: ttother_tcut = "!(gen_partonChannel==2 && ((gen_partonMode1==1 && gen_partonMode2==2) || (gen_partonMode1==2 && gen_partonMode2==1)))"
elif channel == 2: ttother_tcut = "!(gen_partonChannel==2 && (gen_partonMode1==2 && gen_partonMode2==2))"
elif channel == 3: ttother_tcut = "!(gen_partonChannel==2 && (gen_partonMode1==1 && gen_partonMode2==1))"
stepch_tcut =  'step>=%i'%(step)
tcut = '(%s&&%s)*(%s)'%(stepch_tcut,cut,weight)
atcut = '(%s&&%s)*(%s)'%(stepch_tcut,acut,weight)
ttother_tcut = '(%s&&%s&&%s)*(%s)'%(stepch_tcut,cut,ttother_tcut,weight)
rd_tcut = '%s&&%s'%(stepch_tcut,cut)
print "TCut =",tcut
print "anti-TCut =",atcut

#namming
x_name = "Dilepton channel "+x_name
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

#saving mc histos
mchistList = []
amchistList = []
for i, mcname in enumerate(mcfilelist):
	data = findDataSet(mcname, datasets)
	scale = datalumi*data["xsec"]
	colour = data["colour"]
	title = data["title"]
	if 'DYJets' in mcname:
		scale = scale*dyratio[channel][step]

	rfname = rootfileDir + mcname +".root"
	tfile = ROOT.TFile(rfname)
	wentries = tfile.Get("cattree/nevents").Integral()
	scale = scale/wentries
		
	mchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)
	amchist = makeTH1(rfname, tname, title, binning, plotvar,atcut, scale)
	mchist.SetLineColor(colour)
	mchist.SetFillColor(colour)
	amchist.SetLineColor(colour)
	amchist.SetFillColor(colour)

	mchistList.append(mchist)
	amchistList.append(amchist)

	if 'TT' in mcname:
		if len(binning) == 3:
			ttothershist = ROOT.TH1D("name", title+' others', binning[0], binning[1], binning[2])
			attothershist = ROOT.TH1D("name", title+' others', binning[0], binning[1], binning[2])
		else:
			ttothershist = ROOT.TH1D("name", title+' others', len(binning)-1, array.array('f', binning))
			attothershist = ROOT.TH1D("name", title+' others', len(binning)-1, array.array('f', binning))
		for channel in range(1,4):
			if channel == 1: ttother_tcut = "!(gen_partonChannel==2 && ((gen_partonMode1==1 && gen_partonMode2==2) || (gen_partonMode1==2 && gen_partonMode2==1)))"
			elif channel == 2: ttother_tcut = "!(gen_partonChannel==2 && (gen_partonMode1==2 && gen_partonMode2==2))"
			elif channel == 3: ttother_tcut = "!(gen_partonChannel==2 && (gen_partonMode1==1 && gen_partonMode2==1))"
			ttother_tcut = '(channel==%d&&%s&&%s&&%s)*%s'%(channel,stepch_tcut,cut,ttother_tcut,weight)
			ttother_atcut = '(channel==%d&&%s&&%s&&%s)*%s'%(channel,stepch_tcut,acut,ttother_tcut,weight)
			ttothers = makeTH1(rfname, tname, title+' others', binning, plotvar, ttother_tcut, scale)
			attothers = makeTH1(rfname, tname, title+' others', binning, plotvar, ttother_atcut, scale)
			ttothershist.Add(ttothers)
			attothershist.Add(attothers)
		ttothershist.SetLineColor(906)
		ttothershist.SetFillColor(906)
		mchistList.append(ttothershist)
		mchist.Add(ttothershist, -1)

		attothershist.SetLineColor(906)
		attothershist.SetFillColor(906)
		amchistList.append(attothershist)
		amchist.Add(attothershist, -1)


#overflow
if overflow:
	if len(binning) == 3 : nbin = binning[0]
	else : nbin = len(binnin)-1
	for hist in mchistList+amchistList:
		hist.SetBinContent(nbin, hist.GetBinContent(nbin+1))

#bin normalize
if binNormalize and len(binning)!=3:
	for hist in mchistList+amchistList:
		for i in range(len(binning)):
			hist.SetBinContent(i, hist.GetBinContent(i)/hist.GetBinWidth(i))
			hist.SetBinError(i, hist.GetBinError(i)/hist.GetBinWidth(i))
	y_name = y_name + "/%s"%(unit)

#Drawing plots on canvas
var = plotvar.split(',')[0]
#var = ''.join(i for i in var if not i.isdigit())
var = ''.join(i for i in var )
outfile = "Dilepton_s%d_%s%s.png"%(step,var,suffix)


#drawTH1(outfile, CMS_lumi, mchistList, rdhist, x_name, y_name, dolog)

drawTHS_1 = ROOT.THStack("h1","h1")
drawTHS_2 = ROOT.THStack("h1","h1")
for hist in mchistList :
  drawTHS_1.Add(hist)
for hist in amchistList :
  drawTHS_2.Add(hist)

c1 = ROOT.TCanvas("c1","c1",600,600)
h1 = drawTHS_1.GetStack().Last()
h1.SetLineColor(ROOT.kBlue)
h1.SetMarkerColor(ROOT.kBlue)
h2 = drawTHS_2.GetStack().Last()
h2.SetLineColor(ROOT.kRed)
h2.SetMarkerColor(ROOT.kRed)
h2.Draw()
h1.Draw("same")
c1.SaveAs(outfile)


#print outfile

