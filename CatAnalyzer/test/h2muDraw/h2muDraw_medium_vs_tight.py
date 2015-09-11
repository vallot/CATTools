import os
import ROOT,copy,array
ROOT.gROOT.SetBatch(True)

def hist_maker(name, title, bin_set, x_name, y_name, tr, br, cut):
    if (len(bin_set)==3):
        hist = ROOT.TH1F(name, title, bin_set[0], bin_set[1], bin_set[2])
    elif (len(bin_set)==2):
        hist = ROOT.TH1D(name, title, bin_set[0], bin_set[1])
    hist.GetXaxis().SetTitle(x_name)
    hist.GetYaxis().SetTitle(y_name)
    #hist.SetLineColor(color)
    hist.Sumw2()
    #hist.SetStats(0)
    tr.Project(name, br, cut)
    return hist

datalumi = 49.502
currentdir = os.getcwd()
saveddir = '/PLOTS'
if not os.path.isdir(currentdir+saveddir):
    os.mkdir(currentdir+saveddir)
filelist1 = os.listdir("./results_merged")
resultdir = "./results_merged/"
filelist2 = []
for i in filelist1:
    filelist2.append(i.split(".")[0])
filenames = [0,0,0,0,0,0,0]
date = "20150817"
for i in filelist2:#replace 'filelist2' replace to 'filelist1'
    if (('DYJets' in i) and ( date in i)):filenames[0]=i
    if (('TTJets' in i) and ( date in i)):filenames[1]=i
    if (('ZZ' in i) and ( date in i)):filenames[2]=i
    if (('WW' in i) and ( date in i)):filenames[3]=i
    if (('WZ' in i) and ( date in i)):filenames[4]=i
    if (('Double' in i) and ( date in i)):filenames[5]=i
    if (('Single' in i) and ( date in i)):filenames[6]=i
'''
mcfilelist = {'mc_DYJets_merged':6025.2 ,
              'mc_TTJets_mad_merged':831.8,
              'mc_ZZ_merged':31.8,
              'mc_WW_merged':65.9,
              'mc_WZ_merged':118.7}
rdfilelist = [#'DoubleMuon',
              'data_MuonEG_merged']
'''
mcfilelist = {filenames[0]:6025.2 ,
              filenames[1]:831.8,
              filenames[2]:31.8,
              filenames[3]:65.9,
              filenames[4]:118.7}
rdfilelist = [filenames[5],
              filenames[6]]


## Jet Category cut
jet0_tight = "(jetcat_f_hier == 1)"
jet0_loose = "(jetcat_f_hier == 2)"
jet1_tight = "(jetcat_f_hier == 3)"
jet1_loose = "(jetcat_f_hier == 4)"
jet2_vbf = "(jetcat_f_hier == 5 || jetcat_f_hier == 7)"
jet2_ggf = "(jetcat_f_hier == 6)"
jet2_loose = "(jetcat_f_hier == 8)"

jetcat = ["0jet_tight","0jet_loose","1jet_tight","1jet_loose","2jet_VBF_tight","2jet_ggF_tight","2jet_loose"]
jetcat_cut = [jet0_tight, jet0_loose, jet1_tight, jet1_loose, jet2_vbf, jet2_ggf, jet2_loose]

##jetcat, in case 0,1jet
BB = "(jetcat_GC == 2)"
BO = "(jetcat_GC == 11)"
BE = "(jetcat_GC == 101)"
OO = "(jetcat_GC == 20)"
OE = "(jetcat_GC == 110)"
EE = "(jetcat_GC == 200)"

jetcat_GC = ["BB","BO","BE","OO","OE","EE"]
jetcat_GC_cut = [BB,BO,BE,OO,OE,EE]

## initial cut
init_cut = ["(step == 2 && isTight)","(step == 2 && isMedium)"]
init_config = ["-Tight","-Medium"]

os.chdir(resultdir)
plotvar = "ll_m"
title = plotvar
x_name = "M [GeV]"
y_name = "efficiency [Tight/Medium]"
canvas = ROOT.TCanvas(title,title)
bins_N = 20
# 0 to 80, 10 bins; 80 to 150, 7 bins; 150 to 200, 2bins;
Edges = [0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 90, 100, 110, 120, 130, 140, 150, 175, 200, 250]
Edges_bin = array.array('d', Edges)
bin_set = [bins_N, Edges_bin]
h_eff = []
for i in range(len(filenames)):
    h_eff.append(ROOT.TH1D("pT_efficiency_of_TightMuon_per_MediumMuon%d"%i, "pT-efficiency of TightMuon/MediumMuon; %s; %s"%(x_name,y_name), bins_N, Edges_bin))
h_eff_tmp = ROOT.TH1D(plotvar,plotvar, bin_set[0], bin_set[1])
h_tight = []
h_medium = []
leg = ROOT.TLegend(0.7,0.7,0.9,0.9)
leg1 = ROOT.TLegend(0.5,0.35,0.7,0.55)
logscale = False
for ii,init in enumerate(init_cut):
    hs = ROOT.THStack(plotvar,plotvar)
    tcut=init
    j=1
    for i in mcfilelist:
        j=j+1
        rootfilename = i+".root"
        print rootfilename
        samplename = i.strip().split("_")[0]+init_config[ii]
        tt = ROOT.TFile(rootfilename)
        tree = tt.h2mu.Get("tree")
        # untill better way to get nentries
        tempdraw = plotvar +" >> temp" +samplename
        tree.Draw(tempdraw)
        temphist = ROOT.gDirectory.Get("temp" +samplename)
        # untill better way to get nentries
        scale = mcfilelist[i]*datalumi / temphist.GetEntries()
        print mcfilelist[i]
        histo = copy.deepcopy(hist_maker(samplename, title, bin_set, x_name, y_name, tree, plotvar, tcut))
        histo.SetLineColor(j)
        histo.Scale(scale)
        if ii == 0:
            h_tight.append(histo)
        if ii == 1:
            h_medium.append(histo)
            legendname = i.strip().split("_")[0]
            leg1.AddEntry(h_eff[j-2], legendname,"f")
        print histo
        tt.Close()
    for num,i in enumerate(rdfilelist):
        rootfilename = i+".root"
        print rootfilename
        samplename = i.strip().split("_")[0]+init_config[ii]
        tt = ROOT.TFile(rootfilename)
        tree = tt.h2mu.Get("tree")
        histo = copy.deepcopy(hist_maker(samplename, title, bin_set, x_name, y_name, tree, plotvar, tcut))
        if ii == 0:
            h_tight.append(histo)
        if ii == 1:
            h_medium.append(histo)
            legendname = i.strip().split("_")[0]
            leg1.AddEntry(h_eff[-2+num],"Data-"+legendname,"lep")
        print histo
        tt.Close()
#### plot style ####
st_marker = [20,25]
st_color = [2,6,7,8,9]
####
h_eff[0].SetStats(0)
h_eff[0].Divide(h_tight[0],h_medium[0], 1.0, 1.0, "B")
h_eff[0].GetYaxis().SetRangeUser(0.3,1.1)
h_eff[0].Draw("E1")
for i in range(len(mcfilelist)-1):
    h_eff[i+1].SetStats(0)
    h_eff[i+1].Divide(h_tight[i+1],h_medium[i+1], 1.0, 1.0, "B")
    h_eff[i+1].SetLineColor(st_color[i])
    h_eff[i+1].GetYaxis().SetRangeUser(0.3,1.1)
    h_eff[i+1].Draw("E1same")
k=0
for i in range(len(mcfilelist),len(filenames)):
    h_eff[i].SetStats(0)
    h_eff[i].Divide(h_tight[i],h_medium[i], 1.0, 1.0, "B")
    h_eff[i].SetMarkerStyle(st_marker[k])
    h_eff[i].SetMarkerSize(0.6)
    k += 1
    h_eff[i].GetYaxis().SetRangeUser(0.3,1.1)
    h_eff[i].Draw("same")
leg1.Draw("same")

canvas.SaveAs(currentdir+saveddir+"/"+title+"_vs"+".root")
canvas.SaveAs(currentdir+saveddir+"/"+title+"_vs"+".eps")
canvas.SaveAs(currentdir+saveddir+"/"+title+"_vs"+".png")

del leg

