import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os,copy
from CATTools.CatAnalyzer.histoHelper import *
ROOT.gROOT.SetBatch(True)

datalumi = 1280.23

mcfilelist = ['TT_powheg','DYJets','DYJets_10to50']
rdfilelist = ['MuonEG_Run2015','DoubleEG_Run2015','DoubleMuon_Run2015']
rootfileDir = "/cms/scratch/jlee/v7-4-4/TtbarDiLeptonAnalyzer_"
datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset.json" % os.environ['CMSSW_BASE']))

mchistList = []
#plotvar = 'nvertex'
#binning = [30, 0, 30]
plotvar = 'll_m'
binning = [50, 0, 200]
channel = 2
step = 4
x_name = 'mass [GeV]'
y_name = 'events'
CMS_lumi.lumi_sqrtS = "%.0f pb^{-1}, #sqrt{s} = 13 TeV 25ns "%(datalumi)
tname = "cattree/nom"

weight = 'weight'
if channel == 1: ttother_tcut = "!(parton_channel==2 && ((parton_mode1==1 && parton_mode2==2) || (parton_mode1==2 && parton_mode2==1)))"
elif channel == 2: ttother_tcut = "!(parton_channel==2 && (parton_mode1==2 && parton_mode2==2))"
elif channel == 3: ttother_tcut = "!(parton_channel==2 && (parton_mode1==1 && parton_mode2==1))"
cut = '(step>=%i&&channel==%i&&filtered==1)*%s'%(step,channel,weight)
ttother_tcut = '(step>=%i && channel == %i && filtered == 1 && %s)*%s'%(step,channel,ttother_tcut,weight)
print "TCut =",cut

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
    
    mc_ee_in = makeTH1(rfname,tname,"mc_ee_in", binning, plotvar, dycut+'(filtered==1 && channel==2 && step2==0)*(%s)'%(weight), scale)
    mc_mm_in = makeTH1(rfname,tname,"mc_mm_in", binning, plotvar, dycut+'(filtered==1 && channel==3 && step2==0)*(%s)'%(weight), scale)
    mc_ee_out = makeTH1(rfname,tname,"mc_ee_in", binning, plotvar, dycut+'(filtered==1 && channel==2 && step2==1)*(%s)'%(weight), scale)
    mc_mm_out = makeTH1(rfname,tname,"mc_mm_in", binning, plotvar, dycut+'(filtered==1 && channel==3 && step2==1)*(%s)'%(weight), scale)

    rfname = rootfileDir + 'DYJets_10to50' +".root"
    data = findDataSet('DYJets_10to50', datasets)
    scale = datalumi*data["xsec"]
    wentries = getWeightedEntries(rfname, tname, "tri",weight)
    scale = scale/wentries
    mc_ee_in.Add(makeTH1(rfname,tname,"mc_ee_in", binning, plotvar, dycut+'(filtered==1 && channel==2 && step2==0)*(%s)'%(weight), scale))
    mc_mm_in.Add(makeTH1(rfname,tname,"mc_mm_in", binning, plotvar, dycut+'(filtered==1 && channel==3 && step2==0)*(%s)'%(weight), scale))
    mc_ee_out.Add(makeTH1(rfname,tname,"mc_ee_in", binning, plotvar, dycut+'(filtered==1 && channel==2 && step2==1)*(%s)'%(weight), scale))
    mc_mm_out.Add(makeTH1(rfname,tname,"mc_mm_in", binning, plotvar, dycut+'(filtered==1 && channel==3 && step2==1)*(%s)'%(weight), scale))
    
    rfname = rootfileDir+rdfilelist[1-1]+".root"
    rd_em_in = makeTH1(rfname, tname,'rd_ee_in', binning, plotvar, dycut+'(filtered==1 && channel==1 && ((ll_m > 76) && (ll_m < 106)))')
    rfname = rootfileDir + rdfilelist[2-1] +".root"
    rd_ee_in = makeTH1(rfname, tname,'rd_ee_in', binning, plotvar, dycut+'(filtered==1 && channel==2 && step2 ==0)')
    rfname = rootfileDir + rdfilelist[3-1] +".root"
    rd_mm_in = makeTH1(rfname, tname,'rd_ee_in', binning, plotvar, dycut+'(filtered==1 && channel==3 && step2 ==0)')

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
    
    mchist = makeTH1(rfname, tname, title, binning, plotvar, cut, scale)
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
rdhist = makeTH1(rfname, tname, 'data', binning, plotvar, cut)
rdhist.SetLineColor(1)

drawTH1(plotvar+cut+".png", CMS_lumi, mchistList, rdhist, x_name, y_name, True)
