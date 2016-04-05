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
#rootfileDir = "/cms/scratch/jlee/v763desy/TtbarDiLeptonAnalyzer_"
#rootfileDir = "/xrootd/store/user/tt8888tt/v763/TtbarDiLeptonAnalyzer_"
rootfileDir = "/xrootd/store/user/tt8888tt/TtbarDiLeptonAnalyzer_"
channel_name = ['MuEl', 'ElEl', 'MuMu']

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/dataset.json" % os.environ['CMSSW_BASE']))

#defalts
step = 1
channel = 1
cut = 'tri!=0&&filtered==1'
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
if   channel == 1: ttother_tcut = "!(gen_partonChannel==2 && ((gen_partonMode1==1 && gen_partonMode2==2) || (gen_partonMode1==2 && gen_partonMode2==1)))"
elif channel == 2: ttother_tcut = "!(gen_partonChannel==2 && (gen_partonMode1==2 && gen_partonMode2==2))"
elif channel == 3: ttother_tcut = "!(gen_partonChannel==2 && (gen_partonMode1==1 && gen_partonMode2==1))"
#stepch_tcut =  'step7&&step>=5&&channel==%i'%(channel)
stepch_tcut =  'step>=%i&&channel==%i'%(step,channel)

tcut = '(%s&&%s)*(%s)'%(stepch_tcut,cut,weight)
gentop_tcut = '(%s&&%s)*(%s)'%(cut,'!'+ttother_tcut,weight)
signal_tcut = '(%s&&%s&&%s)*(%s)'%(stepch_tcut,cut,'!'+ttother_tcut,weight)
ttother_tcut = '(%s&&%s&&%s)*(%s)'%(stepch_tcut,cut,ttother_tcut,weight)
rd_tcut = '%s&&%s'%(stepch_tcut,cut)
print "TCut =",tcut

#namming
x_name = channel_name[channel-1]+" "+x_name
unit = ""
if len(binning) == 3:
    num = (binning[2]-binning[1])/float(binning[0])
    if num != 1:
        if x_name.endswith(']'):
            unit = "["+x_name.split('[')[1]
        y_name = y_name + "/%g%s"%(num,unit)

#DYestimation
if not os.path.exists('./DYFactor.json'):
    DYestimation.printDYFactor(rootfileDir, tname, datasets, datalumi, cut, weight, rdfilelist)# <------ This will create 'DYFactor.json' on same dir.
dyratio=json.load(open('./DYFactor.json'))

#saving mc histos
errList = []
mchistList = []
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
		
    N = getWeightedEntries(rfname, tname, 'tri', tcut)
    err = math.sqrt(abs(N))*scale

    mchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)
    mchist.SetLineColor(colour)
    mchist.SetFillColor(colour)
    mchistList.append(mchist)
    if 'TT' in mcname:
        N = getWeightedEntries(rfname, tname, 'tri', ttother_tcut)
        ttothers_err = math.sqrt(abs(N))*scale
        errList.append(adderrs(err, ttothers_err, -1))
        err = ttothers_err

        tthist = mchist.Clone()
        ttothers = makeTH1(rfname, tname, title+' others', binning, plotvar, ttother_tcut, scale)
        ttothers.SetLineColor(906)
        ttothers.SetFillColor(906)
        mchistList.append(ttothers)
        mchist.Add(ttothers, -1)
        ttsignal = mchist.Clone()

        if step >= 6:
            parton_plotvar = ','.join('gen_%s'%var for var in plotvar.split(','))
            genhist = makeTH1(rfname, tname, title+' gen', binning, parton_plotvar, gentop_tcut, scale)
            gensignal = makeTH1(rfname, tname, title+' gensignal', binning, parton_plotvar, signal_tcut, scale)

            eff_hist = ROOT.TH1D('title','name',100,-2,2)
            residual = ['(gen_%s-%s)/gen_%s'%(var,var,var) for var in plotvar.split(',')]
            for resi in residual:
                eff_hist.Add(makeTH1(rfname, tname, title+' efficiency', [100,-2,2], resi, signal_tcut))

    errList.append(err)

#data histo
rfname = rootfileDir + rdfilelist[channel-1] +".root"
rdhist = makeTH1(rfname, tname, 'data', binning, plotvar, rd_tcut)
rdhist.SetLineColor(1)

#get event yeild table
num, err = table(mchistList, errList, mchistList[0], errList[0])
for k in num.keys():
    print '%s  ~&~ $%8d \\pm %6.2f$'%(k, max(0,num[k]), err[k])
print 'data  ~&~ $%8d           $ \n'%rdhist.Integral(0,binning[0]+1)

#overflow
if overflow:
    if len(binning) == 3 : nbin = binning[0]
    else : nbin = len(binnin)-1
    for hist in mchistList:
        hist.SetBinContent(nbin, hist.GetBinContent(nbin+1))
    rdhist.SetBinContent(nbin, rdhist.GetBinContent(nbin+1))

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
outfile = "%s_s%d_%s.png"%(channel_name[channel-1],step,var)
drawTH1(outfile, CMS_lumi, mchistList, rdhist, x_name, y_name, dolog)
print outfile

if step >= 6:
    #purity
    purity = ttsignal.Clone()
    stability = gensignal.Clone()
    purity.Divide(tthist)
    stability.Divide(genhist)

    purity.GetXaxis().SetTitle(x_name)
    purity.GetYaxis().SetTitle(y_name)
    purity.SetMaximum(1.1)
    purity.SetMinimum(0)
    purity.SetLineColor(4)
    purity.SetLineWidth(2)
    purity.SetFillColor(0)
    stability.SetLineColor(2)
    stability.SetLineWidth(2)
    stability.SetFillColor(0)

    cn=ROOT.TCanvas()
    purity.Draw('hist')
    stability.Draw('histsame')
    purity.Draw('histsame')
    outfile = "%s_s%d_%s_purity.png"%(channel_name[channel-1],step,var)
    cn.Print(outfile)

    #efficiency
    cn=ROOT.TCanvas()
    eff_hist.SetLineColor(1)
    eff_hist.SetFillColor(0)
    eff_hist.Draw('hist')
    outfile = "%s_s%d_%s_eff.png"%(channel_name[channel-1],step,var)
    cn.Print(outfile)
