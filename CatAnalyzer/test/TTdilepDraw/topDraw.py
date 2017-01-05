#!/usr/bin/env python
import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os, getopt, sys, math
from CATTools.CatAnalyzer.histoHelper import *
import DYestimation
ROOT.gROOT.SetBatch(True)
'''
topDraw.py -a 1 -s 1 -c 'tri==1&&filtered==1' -b [40,0,40] -p nvertex -x 'no. vertex' &
topDraw.py -a 1 -s 1 -b [100,-3,3] -p lep1_eta,lep2_eta -x '#eta' &
'''
datalumi = 2.17
CMS_lumi.lumi_sqrtS = "%.1f fb^{-1}, #sqrt{s} = 13 TeV"%(datalumi)
datalumi = datalumi*1000 # due to fb

mcfilelist = ['TT_powheg', 'WJets', 'SingleTbar_tW', 'SingleTop_tW', 'ZZ', 'WW', 'WZ', 'DYJets', 'DYJets_10to50']
rdfilelist = ['MuonEG_Run2015','DoubleEG_Run2015','DoubleMuon_Run2015']
rootfileDir = "/xrootd/store/user/tt8888tt/v765/TtbarDiLeptonAnalyzer_"
channel_name = ['MuEl', 'ElEl', 'MuMu']

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/dataset.json" % os.environ['CMSSW_BASE']))

#defalts
step = 1
channel = 1 #combined: channel = 0
cut = 'tri!=0&&filtered==1&&is3lep==2'
weight = 'genweight*puweight*mueffweight*eleffweight*tri'
binning = [60, 20, 320]
plotvar = 'dilep.M()'
x_name = 'mass [GeV]'
y_name = 'Events'
dolog = False
overflow = False
binNormalize = False

#get input
try:
    opts, args = getopt.getopt(sys.argv[1:],"hdnoc:w:b:p:x:y:a:s:",["binNormalize","overflow","cut","weight","binning","plotvar","x_name","y_name","dolog","channel","step"])
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
    elif opt in ("-o", "--overflow"):
        overflow = True
    elif opt in ("-n", "--binNormalize"):
        binNormalize = True

tname = "cattree/nom"

#cut define
ttother_tcut = '!(gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel && gen_partonMode==%d)'%(channel)
stepch_tcut =  'step>=%i&&channel==%i'%(step,channel)
if channel == 0:
    ttother_tcut = '!(gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel && gen_partonMode==channel && gen_partonMode!=0)'
    stepch_tcut =  'step>=%i'%(step)
if step == 6: stepch_tcut = stepch_tcut+'&&step6'#desy smeared

tcut = '(%s&&%s)*(%s)'%(stepch_tcut,cut,weight)
ttother_tcut = '(%s&&%s&&%s)*(%s)'%(stepch_tcut,cut,ttother_tcut,weight)
rd_tcut = '%s&&%s'%(stepch_tcut,cut)
print "TCut =",tcut

#namming
unit = ""
if len(binning) == 3:
    num = (binning[2]-binning[1])/float(binning[0])
    if num != 1:
        if x_name.endswith(']'):
            unit = x_name.split('[')[1]
            unit = unit.split(']')[0]
        y_name = y_name + " / %.2g %s"%(num,unit)
if len(binning) == 3 : nbins = binning[0]
else : nbins = len(binning)-1

#DYestimation
if not os.path.exists('./DYFactor.json'):
    DYestimation.printDYFactor(rootfileDir, tname, datasets, datalumi, cut, weight, rdfilelist)#This will create 'DYFactor.json' on same dir.
dyratio=json.load(open('./DYFactor.json'))

title_l = ["t#bar{t}", "W+jets", "Single top", "Single top", "Dibosons", "Dibosons", "Dibosons", "Z/#gamma^{*}#rightarrow#font[12]{l#lower[-0.4]{+}l#lower[-0.4]{#font[122]{\55}}}", "Z/#gamma^{*}#rightarrow#font[12]{l#lower[-0.4]{+}l#lower[-0.4]{#font[122]{\55}}}"]

#saving mc histos
errList = []
mchistList = []
for i, mcname in enumerate(mcfilelist):
    data = findDataSet(mcname, datasets)
    scale = datalumi*data["xsec"]
    colour = data["colour"]
    #title = data["title"]
    title = title_l[i]
    if ('DYJets' in mcname) and channel!=0:
        scale = scale*dyratio[channel][step]

    rfname = rootfileDir + mcname +".root"
    tfile = ROOT.TFile(rfname)
    wentries = tfile.Get("cattree/nevents").Integral()
    scale = scale/wentries
		
    N = getWeightedEntries(rfname, tname, 'tri', tcut)
    err = math.sqrt(abs(N))*scale

    mchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)
    mchist.SetLineColor(colour)
    mchist.SetFillColor(colour)
    mchistList.append(mchist)
    if 'TT' in mcname:
        mchist.SetTitle(mchist.GetTitle()+"-signal (visible)")
        N = getWeightedEntries(rfname, tname, 'tri', ttother_tcut)
        ttothers_err = math.sqrt(abs(N))*scale
        errList.append(adderrs(err, ttothers_err, -1))
        err = ttothers_err

        tthist = mchist.Clone()
        ttothers = makeTH1(rfname, tname, title+'-others', binning, plotvar, ttother_tcut, scale)
        ttothers.SetLineColor(906)
        ttothers.SetFillColor(906)
        mchistList.append(ttothers)
        mchist.Add(ttothers, -1)
        ttsignal = mchist.Clone()

    errList.append(err)

#data histo
if channel != 0:
    rfname = rootfileDir + rdfilelist[channel-1] +".root"
    rdhist = makeTH1(rfname, tname, 'data', binning, plotvar, rd_tcut)
else:
    rdhist = mchistList[0].Clone()
    rdhist.Reset()
    for i, rdfile in enumerate(rdfilelist):
        rfname = rootfileDir + rdfile +".root"
        rdtcut = 'channel==%d&&%s&&%s'%((i+1),stepch_tcut,cut)
        rdhist.Add(makeTH1(rfname, tname, 'data', binning, plotvar, rdtcut))
rdhist.SetLineColor(1)

#get event yeild table
num, err = table(mchistList, errList, mchistList[0], errList[0])
for k in num.keys():
    print '%s  ~&~ $%8d \\pm %6.2f$'%(k, max(0,num[k]), err[k])
print 'data  ~&~ $%8d           $ \n'%rdhist.Integral(0,nbins+1)

#overflow
if overflow:
    for hist in mchistList:
        hist.SetBinContent(nbins, hist.GetBinContent(nbins)+hist.GetBinContent(nbins+1))
        hist.SetBinContent(1, hist.GetBinContent(1)+hist.GetBinContent(0))
    rdhist.SetBinContent(nbins, rdhist.GetBinContent(nbins)+rdhist.GetBinContent(nbins+1))
    rdhist.SetBinContent(1, rdhist.GetBinContent(1)+rdhist.GetBinContent(0))

#bin normalize
if binNormalize and len(binning)!=3:
    for hist in mchistList:
        for i in range(len(binning)):
            hist.SetBinContent(i, hist.GetBinContent(i)/hist.GetBinWidth(i))
            hist.SetBinError(i, hist.GetBinError(i)/hist.GetBinWidth(i))
    for i in range(len(binning)):
        rdhist.SetBinContent(i, rdhist.GetBinContent(i)/rdhist.GetBinWidth(i))
        rdhist.SetBinError(i, rdhist.GetBinError(i)/rdhist.GetBinWidth(i))
    y_name = y_name + "/%s"%(unit)

#Drawing plots on canvas
var = plotvar.split(',')[0]
var = ''.join(i for i in var if not i.isdigit())
var = var.replace('.','_').lower()
var = var.replace('()','')
outfile = "%s_s%d_%s"%(channel_name[channel-1],step,var)
if channel == 0: outfile = "Dilepton_s%d_%s"%(step,var)
drawTH1(outfile, CMS_lumi, mchistList, rdhist, x_name, y_name, dolog, True, 0.5)
print outfile

