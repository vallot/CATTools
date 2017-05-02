#!/usr/bin/env python
import ROOT,json, os, getopt, sys, math
from CATTools.CatAnalyzer.histoHelper import *
import DYEstimation
ROOT.gROOT.SetBatch(True)
'''
topDraw.py -a 1 -s 1 -c 'tri==1&&filtered==1' -b [40,0,40] -p nvertex -x 'no. vertex' &
topDraw.py -a 1 -s 1 -b [100,-3,3] -p lep1_eta,lep2_eta -x '#eta' &
'''
datalumi = 37.06 #Run2016 v806
#datalumi = 2.17 #Run2015 v765
CMS_lumi.lumi_sqrtS = "%.1f fb^{-1}, #sqrt{s} = 13 TeV"%(datalumi)
datalumi = datalumi*1000 # due to fb
CMS_lumi.writeExtraText = False

#mcfilelist = ['TT_powheg', 'WJets', 'SingleTbar_tW', 'SingleTop_tW', 'ZZ', 'WW', 'WZ', 'DYJets', 'DYJets_10to50']
mcfilelist = ['TT_powheg', 'WJets', 'ZZ', 'WW', 'WZ', 'DYJets', 'DYJets_10to50']
rdfilelist = ['MuonEG_Run2016','DoubleEG_Run2016','DoubleMuon_Run2016']
#rootfileDir = "/xrootd/store/user/tt8888tt/v765/TtbarDiLeptonAnalyzer_"
rootfileDir = "/xrootd/store/user/king11kr/ntuples_TtbarDstar_v806/CMSSW_8_0_26_patch1/"
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

title_l = ["t#bar{t}", "W+jets", "Single top", "Single top", "Diboson", "Diboson", "Diboson", "Z/#gamma^{*}#rightarrow#font[12]{l#lower[-0.4]{+}l#lower[-0.4]{#font[122]{\55}}}", "Z/#gamma^{*}#rightarrow#font[12]{l#lower[-0.4]{+}l#lower[-0.4]{#font[122]{\55}}}"]
sysNameList = ['jer','jes','mu','el','puweight','mueffweight','eleffweight','btagweight']

#DYEstimation
if not os.path.exists('./DYFactor.json'):
    DYEstimation.printDYFactor(rootfileDir, tname, datasets, datalumi, cut, weight, rdfilelist)#This will create 'DYFactor.json' on same dir.
dyratio=json.load(open('./DYFactor.json'))

#saving mc histos
errList = []
mchistList = []
sysErr_up = []
sysErr_dn = []
for sysname in sysNameList:
    sysErr_up.append(defTH1(sysname+'_up', sysname+'_up', binning))
    sysErr_dn.append(defTH1(sysname+'_dn', sysname+'_dn', binning))

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
    for i, sysname in enumerate(sysNameList):
        if 'weight' in sysname:
            sysErr_up[i].Add(makeTH1(rfname, tname, title, binning, plotvar, tcut.replace(sysname,sysname+'_up'), scale))
            sysErr_dn[i].Add(makeTH1(rfname, tname, title, binning, plotvar, tcut.replace(sysname,sysname+'_dn'), scale))
        else:
            sysErr_up[i].Add(makeTH1(rfname, "cattree/%s_u"%sysname, title, binning, plotvar, tcut, scale))
            sysErr_dn[i].Add(makeTH1(rfname, "cattree/%s_d"%sysname, title, binning, plotvar, tcut, scale))

    mchist.SetLineColor(colour)
    mchist.SetFillColor(colour)
    if overflow: overFlow(mchist)
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
        if overflow: overFlow(ttothers)
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
if overflow: overFlow(rdhist)
rdhist.SetLineColor(1)
nbins = rdhist.GetNbinsX()

#namming
unit = ""
if len(binning) == 3 and rdhist.GetBinWidth(1) != 1:
    if x_name.endswith(']'):
        unit = x_name.split('[')[1]
        unit = unit.split(']')[0]
    y_name = y_name + " / %g %s"%(rdhist.GetBinWidth(1),unit)

#error band(sys)
errorBand = copy.deepcopy(rdhist)
errorBand.SetFillColor(14)
errorBand.SetFillStyle(3001)
errorBand.SetMarkerStyle(0)

h_nom = defTH1("nom", "nom", binning)
for h in mchistList:
    h_nom.Add(h)

for i in range(len(sysNameList)):
    overFlow(sysErr_up[i])
    overFlow(sysErr_dn[i])
    sysErr_up[i].Add(h_nom, -1)
    sysErr_dn[i].Add(h_nom, -1)
    for j in range(1,nbins+1):
        maxErr = max(abs(sysErr_up[i].GetBinContent(j)), abs(sysErr_dn[i].GetBinContent(j)))
        sumErr = math.sqrt(errorBand.GetBinError(j)**2+maxErr**2) 
        errorBand.SetBinError(j, sumErr)

#get event yeild table
num, err = table(mchistList, errList, mchistList[0], errList[0])
for k in num.keys():
    print '%s  ~&~ $%8d \\pm %6.2f$'%(k, max(0,num[k]), err[k])
print 'data  ~&~ $%8d           $ \n'%rdhist.Integral(0,nbins+1)

#diff
"""
mcdiff = defTH1('MC (Simulation)', "mctotal", binning)
for h in mchistList:
    mcdiff.Add(h)
mcdiff.Scale(1/float(mcdiff.Integral()))
rdhist.Scale(1/float(rdhist.Integral()))
errorBand.Scale(1/float(errorBand.Integral()))
mcdiff.SetLineColor(2)
mcdiff.SetFillColor(0)
mchistList = [mcdiff]
"""

#bin normalize
if binNormalize:
    for hist in mchistList:
        hist.Scale(1,"width")
    rdhist.Scale(1,"width")
    errorBand.Scale(1,"width")

#Drawing plots on canvas
var = plotvar.split(',')[0]
var = ''.join(i for i in var if not i.isdigit())
var = var.replace('.','_').lower()
var = var.replace('()','')
outfile = "%s_s%d_%s.png"%(channel_name[channel-1],step,var)
if channel == 0: outfile = "Dilepton_s%d_%s.png"%(step,var)
canv = drawTH1(outfile, CMS_lumi, mchistList, rdhist, x_name, y_name, dolog, True, 0.5)

## Additional draw
## you can get list of objects on canvavs via 'canv.GetListOfPrimitives()'
## in case you know the name, 'obj = canv.GetPrivmitive(name)'
mainPad = canv.GetPrimitive("mainPad")
mainPad.cd()
errorBand.Draw("e2same")
mainPad.GetPrimitive("data").Draw("esamex0")
extraText(canv, [0.3,0.85], outfile.split("_")[0])
canv.Update()

ratioPad = canv.GetPrimitive("ratioPad")
ratioPad.cd()
sysErrRatio = errorBand.Clone()
sysErrRatio.Divide(rdhist)
sysErrRatio.Draw("e2same")
ratioPad.GetPrimitive("hratio").Draw("esame")
canv.Update()

canv.SaveAs(outfile)
#canv.SaveAs(outfile.replace(".png",".pdf"))
print outfile

