#!/usr/bin/env python
import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os, getopt, sys
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
sgfilelist = ['TTbarDMJets_scalar_Mchi-10_Mphi-10']
sglist = [
    {'name':'TTbarDMJets_scalar_Mchi-10_Mphi-10',    'title':'M_{\chi}=10 GeV, M_{\phi}=10 GeV', 
        'color':4, 'style':1}, 
    {'name':'TTbarDMJets_scalar_Mchi-10_Mphi-100',   'title':'M_{\chi}=10 GeV, M_{\phi}=100 GeV', 
        'color':4, 'style':7}, 
    {'name':'TTbarDMJets_scalar_Mchi-1_Mphi-20',     'title':'M_{\chi}=1 GeV, M_{\phi}=20 GeV', 
        'color':8, 'style':1}, 
    {'name':'TTbarDMJets_scalar_Mchi-1_Mphi-500',    'title':'M_{\chi}=1 GeV, M_{\phi}=500 GeV', 
        'color':8, 'style':7}, 
]
rdfilelist = ['MuonEG_Run2015','DoubleEG_Run2015','DoubleMuon_Run2015']
#rootfileDir = "/cms/scratch/jlee/v763desy/TtbarDiLeptonAnalyzer_"
rootfileDir = "/xrootd/store/user/tt8888tt/v763/TtbarDiLeptonAnalyzer_"
channel_name = ['MuEl', 'ElEl', 'MuMu', 'Dilepton']

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset.json" % os.environ['CMSSW_BASE']))

'''
for entrydata in datasets:
    if 'TTbarDMJets_scalar_' in entrydata["name"] and entrydata["name"] != sgfilelist[0]:
        #print entrydata["name"]
        sgfilelist.append(entrydata["name"])
'''

#defalts
step = 1
channel = 1
cut = 'tri!=0&&filtered==1'
#genweight: , puweight: pileup weight, mueffweignt: , 
#elefffweight: , tri:
weight = 'genweight*puweight*mueffweight*eleffweight*tri'
binning = [60, 20, 320]
plotvar = 'll_m'
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
if   channel == 1: ttother_tcut = "!(parton_channel==2 && ((parton_mode1==1 && parton_mode2==2) || (parton_mode1==2 && parton_mode2==1)))"
elif channel == 2: ttother_tcut = "!(parton_channel==2 && (parton_mode1==2 && parton_mode2==2))"
elif channel == 3: ttother_tcut = "!(parton_channel==2 && (parton_mode1==1 && parton_mode2==1))"

if channel < 4:
    stepch_tcut =  'step>=%i&&channel==%i'%(step,channel)
else:
    stepch_tcut =  'step>=%i'%(step)
'''
tcut = '(%s&&%s)*(%s)'%(stepch_tcut,cut,weight)
ttother_tcut = '(%s&&%s&&%s)*(%s)'%(stepch_tcut,cut,ttother_tcut,weight)
rd_tcut = '%s&&%s'%(stepch_tcut,cut)
'''
'''
rd_tcut = '%s&&(jet1.Pt()+jet2.Pt())<400&&acos(cos(lep1.Phi()-lep2.Phi()))<2&&%s'%(stepch_tcut,cut) #met
#rd_tcut = '%s&&(jet1.Pt()+jet2.Pt())<400&&met>320&&acos(cos(lep1.Phi()-lep2.Phi()))<2&&%s'%(stepch_tcut,cut) #lpt
#rd_tcut = '%s&&met>320&&acos(cos(lep1.Phi()-lep2.Phi()))<2&&%s'%(stepch_tcut,cut) #jetpt
#rd_tcut = '%s&&(jet1.Pt()+jet2.Pt())<400&&met>320&&%s'%(stepch_tcut,cut) #dphi
'''
rd_tcut = '%s'%(stepch_tcut)
if 'met' != plotvar: rd_tcut = '%s&&%s'%(rd_tcut, 'met>320')
if not ('Pt' in plotvar and 'lep' in plotvar): rd_tcut = '%s&&%s'%(rd_tcut, '(lep1.Pt()+lep2.Pt())>120')
if not ('Pt' in plotvar and 'jet' in plotvar): rd_tcut = '%s&&%s'%(rd_tcut, '(jet1.Pt()+jet2.Pt())<400')
if 'Phi' not in plotvar: rd_tcut = '%s&&%s'%(rd_tcut, 'acos(cos(lep1.Phi()-lep2.Phi()))<2')
rd_tcut = '%s&&%s'%(rd_tcut, cut)
tcut = '(%s)*weight'%rd_tcut
print "TCut =",tcut

#namming
x_name = channel_name[channel-1]+" "+x_name
if len(binning) == 3:
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
for i, mcname in enumerate(mcfilelist):
    print "Analyzing ", mcname
    data = findDataSet(mcname, datasets)
    scale = datalumi*data["xsec"]
    colour = data["colour"]
    title = data["title"]
    if channel != 4 and 'DYJets' in mcname:
        scale = scale*dyratio[channel][step]

    rfname = rootfileDir + mcname +".root"
    tfile = ROOT.TFile(rfname)
    wentries = tfile.Get("cattree/nevents").Integral()
    scale = scale/wentries
    
    mchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)
    mchist.SetLineColor(colour)
    mchist.SetFillColor(colour)
    mchistList.append(mchist)
    
'''
    if 'TT_' in mcname:
        ttothers = makeTH1(rfname, tname, title+' others', binning, plotvar, ttother_tcut, scale)
        ttothers.SetLineColor(906)
        ttothers.SetFillColor(906)
        mchistList.append(ttothers)
        mchist.Add(ttothers, -1)
'''

#data histo
print "Analyzing data"
if channel < 4:
    print "In Single channel"
    rfname = rootfileDir + rdfilelist[channel-1] +".root"
    rdhist = makeTH1(rfname, tname, 'data', binning, plotvar, rd_tcut)
    rdhist.SetLineColor(1)
else:
    print "(In combined channel)"
    if len(binning) == 3:
        rdhist = ROOT.TH1D("data", "data", binning[0], binning[1], binning[2])
    else:
        rdhist = ROOT.TH1D("data", "data", len(binning)-1, array.array('f', binning))
    for i, rdfile in enumerate(rdfilelist):
        rfname = rootfileDir + rdfile +".root"
        rd_tcut_tmp = 'channel==%d&&%s&&%s'%((i+1),stepch_tcut,rd_tcut)
        rdhist_tmp = makeTH1(rfname, tname, 'data_tmp', binning, plotvar, rd_tcut_tmp)
        #rdhist_tmp = makeTH1(rfname, tname, 'data_tmp', binning, plotvar, tcut)
        rdhist.SetLineColor(1)
        rdhist.Add(rdhist_tmp)

#signal histo
print "Analyzing signal"

rshistList = []
for i, sig in enumerate(sglist):
    data = findDataSet(sig["name"], datasets)
    scale = datalumi*data["xsec"]
    colour = data["colour"]
    #title = data["title"]
    title = sig["title"]

    rfname = rootfileDir + sig["name"] + ".root"
    tfile = ROOT.TFile(rfname)
    wentries = tfile.Get("cattree/nevents").Integral()
    scale = scale/wentries

    rshist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)
    #rshist.SetLineColor(colour)
    rshist.SetLineColor(sig["color"])
    rshist.SetLineWidth(3)
    rshist.SetLineStyle(sig["style"])
    #rshist.SetLineColor(15 + i)
    #rshist.SetFillColor(0)
    rshistList.append(rshist)

#overflow
if overflow:
    if len(binning) == 3 : nbin = binning[0]
    else : nbin = len(binnin)-1
    
    for hist in mchistList:
        hist.SetBinContent(nbin, hist.GetBinContent(nbin+1))
    
    rdhist.SetBinContent(nbin, rdhist.GetBinContent(nbin+1))
    
    for hist in rshistList:
        hist.SetBinContent(nbin, hist.GetBinContent(nbin+1))

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
outfile = "%s_s%d_%s.png"%(channel_name[channel-1],step,var)
drawTH1(outfile, CMS_lumi, mchistList, rdhist, x_name, y_name, dolog, True, 0.45, rshistList, -0.1)
print outfile


