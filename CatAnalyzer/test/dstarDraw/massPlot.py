#!/usr/bin/env python
import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os, getopt, sys
from CATTools.CatAnalyzer.histoHelper import *
import DYestimation
ROOT.gROOT.SetBatch(True)


datalumi = 2.17
CMS_lumi.lumi_sqrtS = "%.1f fb^{-1}, #sqrt{s} = 13 TeV"%(datalumi)
datalumi = datalumi*1000 # due to fb

topMassList = ['TT_powheg_mtop1665','TT_powheg_mtop1695','TT_powheg_mtop1715','TT_powheg_mtop1735','TT_powheg_mtop1755','TT_powheg_mtop1785']
mcfilelist = ['TT_powheg', 'WJets', 'SingleTbar_tW', 'SingleTop_tW', 'ZZ', 'WW', 'WZ', 'DYJets', 'DYJets_10to50']
rdfilelist = ['MuonEG_Run2015','DoubleEG_Run2015','DoubleMuon_Run2015']
#rootfileDir = "/xrootd/store/user/tt8888tt/v763_desy/TtbarDiLeptonAnalyzer_"
rootfileDir = "/cms/scratch/geonmo/for2016KPS_Ana/src/CATTools/CatAnalyzer/test/cattools/20160414/cattree_"
channel_name = ['Combined', 'MuEl', 'ElEl', 'MuMu']

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/dataset.json" % os.environ['CMSSW_BASE']))

#defalts
step = 1
channel = 1
cut = 'tri!=0&&filtered==1'
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
    opts, args = getopt.getopt(sys.argv[1:],"hdnoc:w:b:p:x:y:a:s:f:",["binNormalize","overflow","cut","weight","binning","plotvar","x_name","y_name","dolog","channel","step","suffix"])
except getopt.GetoptError:          
    print 'Usage : ./.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog> -f <suffix>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./topDraw.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog> -f <suffix>'
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
    elif opt in ("-f", "--suffix"):
        suffix = "_"+arg

tname = "cattree/nom"

#cut define
if   channel == 1: ttother_tcut = "!(gen_partonChannel==2 && ((gen_partonMode1==1 && gen_partonMode2==2) || (gen_partonMode1==2 && gen_partonMode2==1)))"
elif channel == 2: ttother_tcut = "!(gen_partonChannel==2 && (gen_partonMode1==2 && gen_partonMode2==2))"
elif channel == 3: ttother_tcut = "!(gen_partonChannel==2 && (gen_partonMode1==1 && gen_partonMode2==1))"
stepch_tcut =  'step>=%i'%(step)
tcut = '(%s&&%s)*(%s)'%(stepch_tcut,cut,weight)
ttother_tcut = '(%s&&%s&&%s)*(%s)'%(stepch_tcut,cut,ttother_tcut,weight)
rd_tcut = '%s&&%s'%(stepch_tcut,cut)
print "TCut =",tcut

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
for topMass in topMassList :
  mcfilelist[0]=topMass
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
    mchist.SetLineColor(colour)
    mchist.SetFillColor(colour)
    mchistList.append(mchist)
  #overflow
  if overflow:
    if len(binning) == 3 : nbin = binning[0]
    else : nbin = len(binnin)-1
    for hist in mchistList:
      hist.SetBinContent(nbin, hist.GetBinContent(nbin+1))

  #bin normalize
  if binNormalize and len(binning)!=3:
    for hist in mchistList:
      for i in range(len(binning)):
        hist.SetBinContent(i, hist.GetBinContent(i)/hist.GetBinWidth(i))
        hist.SetBinError(i, hist.GetBinError(i)/hist.GetBinWidth(i))
  massValue = int(topMass.split("mtop")[-1])
  output=ROOT.TFile("invMass_"+str(massValue)+"_%s.root"%(plotvar),"RECREATE")
  hs = ROOT.THStack("hs","hs")
  for hist in mchistList :
    hs.Add( hist)
  hs.Draw()
  h_mc = hs.GetStack().Last()
  h_mc.SetName("invMass")
  h_mc.SetTitle("Invariant mass; (l+D*)_{M};Entries/%f"%(h_mc.GetBinWidth(i)))
  h_mc.Write()  
  output.Close()

output = ROOT.TFile.Open("data_%s.root"%(plotvar),"RECREATE")
rdhist=None
mchist=None
if len(binning) == 3:
  rdhist = ROOT.TH1D("Run2015", "RealData in 2015", binning[0], binning[1], binning[2])
else:
  rdhist = ROOT.TH1D("Run2015", "RealData in 2015", len(binning)-1, array.array('f', binning))
for i, rdfile in enumerate(rdfilelist):
  rfname = rootfileDir + rdfile +".root"
  rdtcut = 'channel==%d&&%s&&%s'%((i+1),stepch_tcut,cut)
  rdhist_tmp = makeTH1(rfname, tname, 'data', binning, plotvar, rdtcut)
  rdhist.SetLineColor(1)
  rdhist.Add(rdhist_tmp)
#overflow
if overflow:
  if len(binning) == 3 : nbin = binning[0]
  else : nbin = len(binnin)-1
  rdhist.SetBinContent(nbin, rdhist.GetBinContent(nbin+1))
if binNormalize and len(binning)!=3:
  for i in range(len(binning)):
    rdhist.SetBinContent(i, rdhist.GetBinContent(i)/rdhist.GetBinWidth(i))
    rdhist.SetBinError(i, rdhist.GetBinError(i)/rdhist.GetBinWidth(i))
output.cd()
output.Write()
output.Close()
