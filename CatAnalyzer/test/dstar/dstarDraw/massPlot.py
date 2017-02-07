#!/usr/bin/env python
import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os, getopt, sys, copy
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


datalumi = 27.633 # Run2016 B & C & D & E & F & G, v8-0-2 (H is not yet)
#datalumi = 15.92 # Run2016 B & C & D & E, v8-0-1 (F and latters cannot be used; ICHEP)
CMS_lumi.lumi_sqrtS = "%.1f fb^{-1}, #sqrt{s} = 13 TeV"%(datalumi)
datalumi = datalumi*1000 # due to fb

#topMassList = ['TT_powheg_mtop1665','TT_powheg_mtop1695','TT_powheg_mtop1715','TT_powheg_mtop1735','TT_powheg_mtop1755','TT_powheg_mtop1785','TT_powheg']
topMassList = ['TT_powheg_mtop1665','TT_powheg_mtop1695','TT_powheg_mtop1755','TT_powheg_mtop1785','TT_powheg'] # it will be filled fully
#topMassList = ['TT_powheg_mtop1695','TT_powheg_mtop1755','TT_powheg'] # it will be filled fully
mcfilelist = ['WJets', 'SingleTbar_tW', 'SingleTop_tW', 'ZZ', 'WW', 'WZ', 'DYJets', 'DYJets_10to50']
rdfilelist = ['MuonEG_Run2016','DoubleEG_Run2016','DoubleMuon_Run2016']
channel_name = ['Combined', 'MuEl', 'ElEl', 'MuMu']

dMassNomial = 172.50

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/dataset.json" % os.environ['CMSSW_BASE']))

################################################################
## defalts
################################################################
step = 1
channel = 0
cut = 'tri!=0&&filtered==1'
# In DYJet, genweight yields negative value in histogram(!!!)
weight = 'genweight*puweight*mueffweight*eleffweight*tri'
weightTopPT = '*topPtWeight'
binning = [60, 20, 320]
plotvar = 'll_m'
strType = ""
x_name = 'mass [GeV]'
y_name = 'Events'
dolog = False
overflow = False
binNormalize = False
suffix = ''
strTypeSuffix = ""

################################################################
## get input
################################################################
try:
    opts, args = getopt.getopt(sys.argv[1:],"hdnoa:p:s:c:w:b:t:x:y:f:",["channel","plotvar","step","cut","weight","binning","type","binNormalize","overflow","x_name","y_name","dolog","suffix"])
except getopt.GetoptError:          
    print 'Usage : ./massPlot.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog> -f <suffix>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./topDraw.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog> -f <suffix>'
        sys.exit()
    elif opt in ("-a", "--channel"):
        channel = int(arg)
    elif opt in ("-p", "--plotvar"):
        plotvar = arg
    elif opt in ("-s", "--step"):
        step = int(arg)
    elif opt in ("-c", "--cut"):
        #cut = arg
        cut = "%s&&%s"%(cut,arg)
    elif opt in ("-w", "--weight"):
        weight = arg
    elif opt in ("-b", "--binning"):
        binning = eval(arg)
    elif opt in ("-t", "--type"):
        strType = arg
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
tnameData = tname

################################################################
## Read type
################################################################
if strType != "" :
  strType = "".join(ch for ch in strType if not ch == ' ')

  for strEntry in strType.split(",") :
    listOptArg = strEntry.split("=")
    
    if listOptArg[ 0 ] == "noTopPtW" : 
      weightTopPT = ""
      strTypeSuffix = strTypeSuffix + "_" + listOptArg[ 0 ]
    elif listOptArg[ 0 ] == "TTnominal" : 
      for i in range(len(topMassList)) :
        if topMassList[ i ] == "TT_powheg" : 
          topMassList[ i ] = listOptArg[ 1 ]
          break
      
      strTypeSuffix = strTypeSuffix + "_" + listOptArg[ 0 ] + "_" + listOptArg[ 1 ]
    elif listOptArg[ 0 ] == "JES_Up" : 
      tname = "cattree/jes_u"
      strTypeSuffix = strTypeSuffix + "_" + listOptArg[ 0 ]
    elif listOptArg[ 0 ] == "JES_Down" : 
      tname = "cattree/jes_d"
      strTypeSuffix = strTypeSuffix + "_" + listOptArg[ 0 ]
    elif listOptArg[ 0 ] == "JER_Up" : 
      tname = "cattree/jer_u"
      strTypeSuffix = strTypeSuffix + "_" + listOptArg[ 0 ]
    elif listOptArg[ 0 ] == "JER_Down" : 
      tname = "cattree/jer_d"
      strTypeSuffix = strTypeSuffix + "_" + listOptArg[ 0 ]
    elif listOptArg[ 0 ] == "El_Up" : 
      tname = "cattree/el_u"
      strTypeSuffix = strTypeSuffix + "_" + listOptArg[ 0 ]
    elif listOptArg[ 0 ] == "El_Down" : 
      tname = "cattree/el_d"
      strTypeSuffix = strTypeSuffix + "_" + listOptArg[ 0 ]
    elif listOptArg[ 0 ] == "Mu_Up" : 
      tname = "cattree/mu_u"
      strTypeSuffix = strTypeSuffix + "_" + listOptArg[ 0 ]
    elif listOptArg[ 0 ] == "Mu_Down" : 
      tname = "cattree/mu_d"
      strTypeSuffix = strTypeSuffix + "_" + listOptArg[ 0 ]
      

################################################################
## cut define
################################################################
#if   channel == 1: ttother_tcut = "!(gen_partonChannel==2 && ((gen_partonMode1==1 && gen_partonMode2==2) || (gen_partonMode1==2 && gen_partonMode2==1)))"
#elif channel == 2: ttother_tcut = "!(gen_partonChannel==2 && (gen_partonMode1==2 && gen_partonMode2==2))"
#elif channel == 3: ttother_tcut = "!(gen_partonChannel==2 && (gen_partonMode1==1 && gen_partonMode2==1))"
stepch_tcut = 'step>=%i%s'%(step, "&&channel==%i"%(channel) if channel != 0 else "")
tcutonly = '%s%s'%(stepch_tcut, "&&" + cut if cut != "" else "")
tcut = '(%s)*(%s)'%(tcutonly, weight + weightTopPT)
#ttother_tcut = '(%s&&%s&&%s)*(%s)'%(stepch_tcut,cut,ttother_tcut,weight)
rd_tcut = '%s&&%s'%(stepch_tcut, cut)
print "TCut =",tcut

################################################################
## namming
################################################################
if len(binning) <= 3:
  num = (binning[2]-binning[1])/float(binning[0])
  if num != 1:
    unit = "["+x_name.split('[')[1] if x_name.endswith(']') else ""
    y_name = y_name + "/%g%s"%(num,unit)

################################################################
## DYestimation
################################################################
if not os.path.exists('./DYFactor.json'):
  DYestimation.printDYFactor(rootfileDir, tname, datasets, datalumi, cut, weight, rdfilelist)# <------ This will create 'DYFactor.json' on same dir.
dyratio=json.load(open('./DYFactor.json'))

################################################################
## Initializing the result root file
################################################################
strFilename = "invMass_%s%s%s"%(plotvar,suffix,strTypeSuffix)
outMassHist = ROOT.TFile.Open("invmass/" + strFilename + ".root","RECREATE")
dicListHist = {"rootfilename":strFilename + ".root", 
  "x_name":x_name, "y_name":y_name, "binning":binning, "Gen":0}

################################################################
## Saving MC histograms for backgrounds
################################################################
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
  print "Bkg scales : ", mcname, scale
    
  mchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)
  mchist.SetLineColor(colour)
  mchist.SetFillColor(colour)
  mchistList.append(mchist)
  
  mchist.SetName("hist_bkg_" + mcname)
  outMassHist.cd()
  mchist.Write()

  dicListHist[mchist.GetName()] = {"type":"bkg_part"}
  
#overflow
if overflow:
  nbin = binning[0] if len(binning) == 3 else len(binning)-1
  for hist in mchistList:
    hist.SetBinContent(nbin, hist.GetBinContent(nbin+1))

#bin normalize
if binNormalize and len(binning)!=3:
  for hist in mchistList:
    for i in range(len(binning)):
      hist.SetBinContent(i, hist.GetBinContent(i)/hist.GetBinWidth(i))
      hist.SetBinError(i, hist.GetBinError(i)/hist.GetBinWidth(i))

hs_bkg = ROOT.THStack("bkg_hs","bkg_hs")
for hist in mchistList :
  hs_bkg.Add( hist)

hs_bkg.Draw()
bkgs = hs_bkg.GetStack().Last()

bkgs.SetName("hist_bkg")
bkgs.Draw()

outMassHist.cd()
bkgs.Write()

dicListHist[bkgs.GetName()] = {"type":"bkg"}

print "bkg entries: ",bkgs.GetEntries()

################################################################
##  Saving TT samples with side-band
################################################################
"""
arrInfoSideBandWin = [
  {"name":"center", "mv":0.0, }, 
]

for dMoveWin in [-1.0, 0.0, 1.0]:
  dicInfoWin = 
"""

for topMass in topMassList :
  massValue = topMass.split("mtop")[-1] if ( topMass.find("mtop") != -1 ) else "nominal"
  sum_hs =  hs_bkg.Clone()
  data = findDataSet(topMass, datasets)
  scale = datalumi*data["xsec"]
  colour = data["colour"]
  title = data["title"]

  dMassCurr = int(massValue) * 0.1 if massValue != "nominal" else dMassNomial

  rfname = rootfileDir + topMass + ".root"
  tfile = ROOT.TFile(rfname)
  wentries = tfile.Get("cattree/nevents").Integral()
  scale = scale/wentries
  print topMass, scale, wentries, colour, title
  
  #for dicInfoWin in arrInfoSideBandWin:
  mchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)
  mchist.SetLineColor(colour)
  mchist.SetFillColor(colour)
  print "topmass hsit : ",mchist.Integral()
  # -- Overflow
  if overflow:
    nbin = binning[0] if len(binning) == 3 else len(binning)-1
    mchist.SetBinContent(nbin, mchist.GetBinContent(nbin+1))

  # -- Bin normalize
  if binNormalize and len(binning) != 3:
    for i in range(len(binning)):
      mchist.SetBinContent(i, mchist.GetBinContent(i) / mchist.GetBinWidth(i))
      mchist.SetBinError(i, mchist.GetBinError(i) / mchist.GetBinWidth(i))

  # -- Getting the plot 
  sum_hs.Add( mchist )
  masshist = sum_hs.GetStack().Last()
 
  # -- Saving results into the root file
  masshist.Draw()
  print masshist.GetEntries()
  
  outMassHist.cd()
  mchist.SetName("hist_TT_onlytt_%s"%(massValue))
  mchist.Write()
  
  dicListHist[mchist.GetName()] = {"type":"TT_onlytt", "mass":dMassCurr}
  
  outMassHist.cd()
  masshist.SetName("hist_TT_withbkg_%s"%(massValue))
  masshist.Write()
  
  dicListHist[masshist.GetName()] = {"type":"TT_withbkg", "mass":dMassCurr}

if "correctM" in plotvar:
  dicListHist[ "Gen" ] = 1

################################################################
##  Saving data samples
################################################################
#output = ROOT.TFile.Open("data_%s.root"%(plotvar),"RECREATE")
plotvarData = plotvar
rdtcutPre = '%s%s'%(stepch_tcut, "&&" + cut if cut != "" else "")

if "correctM" in plotvar: 
  if "d0" in plotvar: 
    plotvarData = "d0_lepSV_lowM"
    rdtcutPre = "%s&&%s"%(rdtcutPre, "d0_L3D>0.2&&d0_LXY>0.1&&abs(d0.M()-1.8648)<0.040")
  else:
    plotvarData = "dstar_lepSV_lowM"
    rdtcutPre = "%s&&%s"%(rdtcutPre, "d0_L3D>0.2&&d0_LXY>0.1&&dstar_L3D>0.2&&dstar_LXY>0.1&&abs(dstar_diffMass-0.145)<0.01")

if len(binning) == 3:
  rdhist = ROOT.TH1D("hist_data", "RealData in 2016", binning[0], binning[1], binning[2])
else:
  rdhist = ROOT.TH1D("hist_data", "RealData in 2016", len(binning)-1, array.array('f', binning))
for i, rdfile in enumerate(rdfilelist):
  rfname = rootfileDir + rdfile +".root"
  #rdtcut = 'channel==%d&&%s&&%s'%((i+1),stepch_tcut,cut)
  rdtcut = 'channel==%d&&%s'%((i+1),rdtcutPre)
  rdhist_tmp = makeTH1(rfname, tnameData, 'data', binning, plotvarData, rdtcut)
  rdhist.SetLineColor(1)
  rdhist.Add(rdhist_tmp)
#overflow
if overflow:
  nbin = binning[0] if len(binning) == 3 else len(binning)-1
  rdhist.SetBinContent(nbin, rdhist.GetBinContent(nbin+1))
if binNormalize and len(binning)!=3:
  for i in range(len(binning)):
    rdhist.SetBinContent(i, rdhist.GetBinContent(i)/rdhist.GetBinWidth(i))
    rdhist.SetBinError(i, rdhist.GetBinError(i)/rdhist.GetBinWidth(i))

outMassHist.cd()
rdhist.Write()

dicListHist[rdhist.GetName()] = {"type":"data"}

outMassHist.Write()
outMassHist.Close()

################################################################
##  Saving informations about histograms
################################################################

fileDicHist = open("invmass/" + strFilename + ".json", "w")
fileDicHist.write(json.dumps(dicListHist))
fileDicHist.close()

print ''
print ''
