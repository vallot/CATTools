#!/usr/bin/env python
import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os, getopt, sys
from CATTools.CatAnalyzer.histoHelper import *
import DYestimation
ROOT.gROOT.SetBatch(True)

#rootfileDir = "/xrootd/store/user/tt8888tt/v763_desy/TtbarDiLeptonAnalyzer_"
#rootfileDir = "/cms/scratch/geonmo/for2016KPS_Ana/src/CATTools/CatAnalyzer/test/cattools/cattree_"
#rootfileDir = "/xrootd/store/user/quark2930/dilepton_mass_v801_16092901/cattree_"
strFilenameRootdir = "../rootdirPath.txt"

fRootDir = open(strFilenameRootdir, "r")
rootfileDir = "/xrootd" + fRootDir.readline().split("\r")[0].split("\n")[0] + "/cattree_"
fRootDir.close()


'''
dstarDraw.py -a 1 -s 1 -c 'tri==1&&filtered==1' -b [40,0,40] -p nvertex -x 'no. vertex' &
dstarDraw.py -a 1 -s 1 -b [100,-3,3] -p lep1_eta,lep2_eta -x '#eta' &
'''
datalumi = 27.64 # Run2016 B & C & D & E & F & G, v8-0-2 (H is not yet)
#datalumi = 16.94 # Run2016 B & C & D & E, v8-0-2
#datalumi = 15.92 # Run2016 B & C & D & E, v8-0-1
CMS_lumi.lumi_sqrtS = "%.1f fb^{-1}, #sqrt{s} = 13 TeV"%(datalumi)
datalumi = datalumi*1000 # due to fb
CMS_lumi.extraText   = "Private work"
mcfilelist = ['TT_powheg', 'WJets', 'SingleTbar_tW', 'SingleTop_tW', 'ZZ', 'WW', 'WZ', 'DYJets', 'DYJets_10to50']
rdfilelist = ['MuonEG_Run2016','DoubleEG_Run2016','DoubleMuon_Run2016']
#rdfilelist = ['MuonEG_Run2016_16fb','DoubleEG_Run2016_16fb','DoubleMuon_Run2016_16fb']
channel_name = ['Combined', 'MuEl', 'ElEl', 'MuMu']

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/dataset.json" % os.environ['CMSSW_BASE']))

###############################################################################
## 
## One day, if the problem about trigger in ee and em channel is solved, 
## you MUST return to the original code in DYestimation.py.
## 
###############################################################################

#defalts
#treename = 'nom'
treename = 'nom'
step = 1
channel = 3
#cut = 'tri!=0&&filtered==1'
cut = ''
#weight = 'tri*genweight*puweight*mueffweight*eleffweight*topPtWeight'
weight = 'genweight*puweight*mueffweight*eleffweight*topPtWeight'
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
    opts, args = getopt.getopt(sys.argv[1:],"hdnot:c:w:b:p:x:y:a:s:f:",["binNormalize","overflow","treename","cut","weight","binning","plotvar","x_name","y_name","dolog","channel","step","suffix"])
except getopt.GetoptError:          
    print 'Usage : ./topDraw.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog> -f <suffix>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./topDraw.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog> -f <suffix>'
        sys.exit()
    elif opt in ("-t", "--treename"):
        treename = arg
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

tname = "cattree/%s"%(treename)

#cut define
stepch_tcut = 'step>=%i%s'%(step, "&&channel==%i"%(channel) if channel != 0 else "")
tcutonly = '%s%s'%(stepch_tcut, "&&" + cut if cut != "" else "")
tcut = '(%s)*(%s)'%(tcutonly, weight)

#if   channel == 1: 
#  ttother_tcut = "!(gen_partonChannel==2 && ((gen_partonMode1==1 && gen_partonMode2==2) || (gen_partonMode1==2 && gen_partonMode2==1)))"
#elif channel == 2: ttother_tcut = "!(gen_partonChannel==2 && (gen_partonMode1==2 && gen_partonMode2==2))"
#elif channel == 3: ttother_tcut = "!(gen_partonChannel==2 && (gen_partonMode1==1 && gen_partonMode2==1))"

#ttother_tcut = '(%s&&%s&&%s)*(%s)'%(stepch_tcut,cut,ttother_tcut,weight)
#rd_tcut = '%s&&%s'%(stepch_tcut,cut)

print "TCut =",tcut

#namming
x_name = x_name
if len(binning) <= 3:
  num = (binning[2]-binning[1])/float(binning[0])
  if num != 1:
    unit = "["+x_name.split('[')[1] if x_name.endswith(']') else ""
    y_name = y_name + "/%g%s"%(num,unit)

#DYestimation
if not os.path.exists('./DYFactor.json'):
  DYestimation.printDYFactor(rootfileDir, tname, datasets, datalumi, cut, weight, rdfilelist)# <------ This will create 'DYFactor.json' on same dir.
dyratio=json.load(open('./DYFactor.json'))

#saving mc histos
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

  mchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)
  mchist.SetLineColor(colour)
  mchist.SetFillColor(colour)
  mchistList.append(mchist)
  
  print "%s : %f, %f, %f, %f"%(mcname, scale, scale * wentries, wentries, mchist.Integral())
 
  if 'TT' in mcname:
    if len(binning) == 3:
      ttothershist = ROOT.TH1D("name_others", title+' others', binning[0], binning[1], binning[2])
    else:
      ttothershist = ROOT.TH1D("name_others", title+' others', len(binning)-1, array.array('f', binning))
    dstar_tcut = "((%s)&&abs(dstar_relPtTrue)<0.1&&abs(dstar_dRTrue)<0.1&&abs(d0_dRTrue)<0.1&&abs(d0_relPtTrue)<0.1)*(%s)"%(tcutonly,weight)
    d0_tcut = "((%s)&&(abs(d0_relPtTrue)<0.1&&abs(d0_dRTrue)<0.1)&&!(abs(dstar_relPtTrue)<0.1&&abs(dstar_dRTrue)<0.1))*(%s)"%(tcutonly,weight)
    d0_true_hist = makeTH1(rfname, tname, title+' D0 Signal', binning, plotvar, d0_tcut, scale)
    dstar_true_hist = makeTH1(rfname, tname, title+' D* Signal', binning, plotvar, dstar_tcut, scale)
    #ttddothershist.Add(ttothers)
    d0_true_hist.SetLineColor(906)
    d0_true_hist.SetFillColor(906)
    dstar_true_hist.SetLineColor(806)
    dstar_true_hist.SetFillColor(806)

    mchistList.append(d0_true_hist)
    mchistList.append(dstar_true_hist)
    mchist.Add(d0_true_hist, -1)
    mchist.Add(dstar_true_hist, -1)

#data histo
if len(binning) == 3:
  rdhist = ROOT.TH1D("name_data", title, binning[0], binning[1], binning[2])
else:
  rdhist = ROOT.TH1D("name_data", title, len(binning)-1, array.array('f', binning))

testweight = [1.0 / 0.75, 1.0 / 0.6, 1.0 / 1.0]
for i, rdfile in enumerate(rdfilelist):
  if channel != 0 and i + 1 != channel : continue
  rfname = rootfileDir + rdfile +".root"
  rdtcut = '(channel==%d&&%s%s)*(%f)'%((i+1), stepch_tcut, "&&" + cut if cut != "" else "", testweight[ i ])
  rdhist_tmp = makeTH1(rfname, tname, 'data', binning, plotvar, rdtcut)
  rdhist.SetLineColor(1)
  rdhist.Add(rdhist_tmp)

#overflow
if overflow:
  nbin = binning[0] if len(binning) == 3 else len(binning)-1
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
#var = ''.join(i for i in var if not i.isdigit())
var = ''.join(i for i in var )
outfile = "Dilepton_%s_ch_%i_%s_s%d%s"%(treename,channel,var,step,suffix)
drawTH1(outfile, CMS_lumi, mchistList, rdhist, x_name, y_name, dolog)
print outfile

