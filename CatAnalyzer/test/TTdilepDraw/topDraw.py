import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os,copy
from CATTools.CatAnalyzer.histoHelper import *
ROOT.gROOT.SetBatch(True)

datalumi = 1280.23

mcfilelist = ['TT_powheg','DYJets']
rdfilelist = ['MuonEG','DoubleEG','DoubleMuon']
rootfileDir = "/cms/scratch/tt8888tt/cattools_v744/src/CATTools/CatAnalyzer/test/result3/"
datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset.json" % os.environ['CMSSW_BASE']))

mchistList = []
#plotvar = 'nvertex'
#binning = [30, 0, 30]
plotvar = 'll_m'
binning = [50, 0, 200]
channel = 2
step = 4
cut = '(step>=%i && channel == %i && filtered == 1)'%(step,channel)
print "TCut =",cut
x_name = 'mass [GeV]'
y_name = 'events'
CMS_lumi.lumi_sqrtS = "%.0f pb^{-1}, #sqrt{s} = 13 TeV 25ns "%(datalumi)
tname = tname

#DY estimation
dyratio = [[0 for x in range(6)] for x in range(4)]
dyratio[1][step] = 1.
if channel !=1:
    rfname = rootfileDir + 'DYJets' +".root"
    scale = 1.
    for data in datasets:
        if data["name"] == 'DYJets':
            scale = datalumi*data["xsec"]

    dycut = ""
    if step == 1:
        dycut = "(step1==1)*"
    if step == 2:
        dycut = "(step1==1)*"
    if step == 3:
        dycut = "(step1==1)*(step3==1)*"
    if step == 4:
        dycut = "(step1==1)*(step3==1)*(step4==1)*"
    if step == 5:
        dycut = "(step1==1)*(step3==1)*(step4==1)*(step5==1)*"

    mc_ee_in = makeTH1(rfname,tname,"mc_ee_in", binning, plotvar, dycut+'(filtered==1 && channel==2 && step2==0)*(puweight)', scale)
    mc_mm_in = makeTH1(rfname,tname,"mc_mm_in", binning, plotvar, dycut+'(filtered==1 && channel==3 && step2==0)*(puweight)', scale)
    mc_ee_out = makeTH1(rfname,tname,"mc_ee_in", binning, plotvar, dycut+'(filtered==1 && channel==2 && step2==1)*(puweight)', scale)
    mc_mm_out = makeTH1(rfname,tname,"mc_mm_in", binning, plotvar, dycut+'(filtered==1 && channel==3 && step2==1)*(puweight)', scale)

    rfname = rootfileDir+rdfilelist[1-1]+".root"
    rd_em_in = makeTH1(rfname, tname,'rd_ee_in', binning, plotvar, dycut+'(filtered==1 && channel==1 && ((ll_m > 76) && (ll_m < 106)))')
    rfname = rootfileDir + rdfilelist[2-1] +".root"
    rd_ee_in = makeTH1(rfname, tname,'rd_ee_in', binning, plotvar, dycut+'(filtered==1 && channel==2 && step2 ==0)')
    rfname = rootfileDir + rdfilelist[3-1] +".root"
    rd_mm_in = makeTH1(rfname, tname,'rd_ee_in', binning, plotvar, dycut+'(filtered==1 && channel==3 && step2 ==0)')

    print mc_ee_in.Integral(), mc_ee_out.Integral(), mc_mm_in.Integral(), mc_mm_out.Integral(), rd_ee_in.Integral(), rd_mm_in.Integral(), rd_em_in.Integral()
    dyest = drellYanEstimation(mc_ee_in.Integral(), mc_ee_out.Integral(), mc_mm_in.Integral(), mc_mm_out.Integral(),
                               rd_ee_in.Integral(), rd_mm_in.Integral(), rd_em_in.Integral())
    print step, "ee =",dyest[0], "mm =",dyest[1]   
    dyratio[2][step] = dyest[0]
    dyratio[3][step] = dyest[1]


for i, mcname in enumerate(mcfilelist):
    rfname = rootfileDir + mcname +".root"
    mccut = cut+'*(puweight)'
    scale = 1.
    for data in datasets:
        if data["name"] == mcname:
            scale = datalumi*data["xsec"]
            if mcname == 'DYJets':
                scale = scale*dyratio[channel][step]
                print dyratio[channel][step]
    mchist = makeTH1(rfname, tname, mcname, binning, plotvar, mccut,scale)    
    
    mchist.SetFillColor(i+2)
    mchistList.append(mchist)


rfname = rootfileDir + rdfilelist[channel-1] +".root"
rdhist = makeTH1(rfname, tname,'data', binning, plotvar, cut)
rdhist.SetLineColor(1)

drawTH1(plotvar+".png", CMS_lumi, mchistList, rdhist, x_name, y_name, True)
