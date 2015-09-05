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
def ini_hist_maker(name, title, bin_set, x_name, y_name):
    if (len(bin_set)==3):
        hist = ROOT.TH1F(name, title, bin_set[0], bin_set[1], bin_set[2])
    elif (len(bin_set)==2):
        hist = ROOT.TH1D(name, title, bin_set[0], bin_set[1])
    hist.GetXaxis().SetTitle(x_name)
    hist.GetYaxis().SetTitle(y_name)
    return hist
def TypeOfMuons(type1, type2, h_gen_t, h_reco_t, h_gen_m, h_reco_m, h_gen_l, h_reco_l,histo):
    ## Usage #####################
    ##############################
    # num      0       1      2  #
    # type1   gen    reco        #
    # type2  tight  medium  loose#
    ##############################
    if (type1 == 0 and type2 == 0):
        h_gen_t.append(histo)
        return h_gen_t
    if (type1 == 1 and type2 == 0):
        h_reco_t.append(histo)
        return h_reco_t
    if (type1 == 0 and type2 == 1):
        h_gen_m.append(histo)
        return h_gen_m
    if (type1 == 1 and type2 == 1):
        h_reco_m.append(histo)    
        return h_reco_m
    if (type1 == 0 and type2 == 2):
        h_gen_l.append(histo)
        return h_gen_l
    if (type1 == 1 and type2 == 2):
        h_reco_l.append(histo)    
        return h_reco_l

datalumi = 49.502
currentdir = os.getcwd()
saveddir = '/PLOTS'
if not os.path.isdir(currentdir+saveddir):
    os.mkdir(currentdir+saveddir)
filelist1 = os.listdir("./results_merged")
#filelist1 = os.listdir("./")
resultdir = "./results_merged/"
#resultdir = "./"
filelist2 = []
for i in filelist1:
    if not '.root' in i:continue
    filelist2.append(i.split(".")[0])
filenames = [0,0,0,0,0,0,0]
date = "20150905"
#date = ''
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

## jetcat, in case 0,1jet
BB = "(jetcat_GC == 2)"
BO = "(jetcat_GC == 11)"
BE = "(jetcat_GC == 101)"
OO = "(jetcat_GC == 20)"
OE = "(jetcat_GC == 110)"
EE = "(jetcat_GC == 200)"

jetcat_GC = ["BB","BO","BE","OO","OE","EE"]
jetcat_GC_cut = [BB,BO,BE,OO,OE,EE]


##### initial cut
init_cut = ["(lep_isTight)","(lep_isMedium)","(lep_isLoose)"]
init_config = ["-Tight","-Medium","-Loose"]

##### eta cut
eta_cut = ["(gen_lep_eta < 2.4)","(gen_lep_eta < 0.8)","(gen_lep_eta > 0.8)*(gen_lep_eta < 1.5)","(gen_lep_eta > 1.5)*(gen_lep_eta < 2.4)"]
eta_config = ["all","barrel","overlap","endcap"]

##### variables
plotvar_pt = ["gen_lep_pt","reco_lep_pt"]
plotvar_eta = ["gen_lep_eta","reco_lep_eta"]
plotvar_rsl = ["resolution"]

##### orders : pt, eta, resolution
name_l = ["pT_efficiency_of_recoMuon_per_genMuon","eta_efficiency_of_recoMuon_per_genMuon","resolution"]
title_l = ["pT-efficiency of recoMuon/genMuon","eta-efficiency of recoMuon/genMuon","Muon resolution"]
x_name_l = ["M [GeV]","#eta","p_{T}^{gen}-p_{T}^{reco}/p_{T}^{gen}"]
y_name_l = ["Efficiency [mu_{reco}/mu_{gen}]","Efficiency [mu_{reco}/mu_{gen}]","number of entries"]
bins_N = 20
Edges = [0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 90, 100, 110, 120, 130, 140, 150, 175, 200, 250]
Edges_bin = array.array('d', Edges)
bin_set_l = [[bins_N, Edges_bin], [35, -3.5, 3.5], [100, -0.2, 0.2]]

h_eff = []
for i in range(3):
    h_eff_tmp = []
    for j in range(len(filenames)):
        h_eff_tmp.append(ini_hist_maker(name_l[i]+"%d"%j, title_l[i], bin_set_l[i], x_name_l[i], y_name_l[i]))
    h_eff.append(h_eff_tmp)
#####

h_gen_t = []
h_reco_t = []
h_gen_m = []
h_reco_m = []
h_gen_l = []
h_reco_l = []
os.chdir(resultdir)

## pt-eff
j=1
leg = ROOT.TLegend(0.7,0.7,0.9,0.9)
leg1 = ROOT.TLegend(0.5,0.35,0.7,0.55)
for init_i, init_loop in enumerate(init_cut):
    for pt_i, pt_loop in enumerate(plotvar_pt):
        title = pt_loop
        canvas = ROOT.TCanvas(title,title)
        logscale = False
        if pt_i==0:
            tcut = ''
        else:
            tcut = init_loop
        print tcut
        for i in mcfilelist:
            if not 'DY' in i:
                continue
            j=j+1
            rootfilename = i+".root"
            print rootfilename
            samplename = i.strip().split("_")[0]+init_config[init_i]
            tt = ROOT.TFile(rootfilename)
            tree = tt.h2mu.Get("tree2")
            # untill better way to get nentries
            tempdraw = pt_loop +" >> temp" +samplename
            tree.Draw(tempdraw)
            temphist = ROOT.gDirectory.Get("temp" +samplename)
            # untill better way to get nentries
            scale = mcfilelist[i]*datalumi / temphist.GetEntries()
            print mcfilelist[i]
            histo = copy.deepcopy(hist_maker(samplename, title, bin_set_l[0], x_name_l[0], y_name_l[0], tree, pt_loop, tcut))
            histo.SetLineColor(j)
            histo.Scale(scale)
            TypeOfMuons(pt_i, init_i, h_gen_t, h_reco_t, h_gen_m, h_reco_m, h_gen_l, h_reco_l, histo)
            print histo
            tt.Close()
#### plot style ####
st_marker = [20,25]
st_color = [2,6,7,8,9]
st_gen_tml = [h_gen_t,h_gen_m,h_gen_l]
st_reco_tml = [h_reco_t,h_reco_m,h_reco_l]
st_draw = ["","same","same"]
samplename = filenames[0].strip().split("_")[0]
####
for i in range(3):
    h_eff[0][i].SetStats(0)
    h_eff[0][i].Divide(st_reco_tml[i][0],st_gen_tml[i][0], 1.0, 1.0, "B")
    h_eff[0][i].SetLineColor(st_color[i])
    h_eff[0][i].Draw("E1%s"%st_draw[i])
    leg1.AddEntry(h_eff[0][i],samplename+init_config[i],"f")
leg1.Draw("same")

canvas.SaveAs(currentdir+saveddir+"/"+"pt_vs2"+".root")
canvas.SaveAs(currentdir+saveddir+"/"+"pt_vs2"+".eps")
canvas.SaveAs(currentdir+saveddir+"/"+"pt_vs2"+".png")
del leg1

del h_gen_t[:]
del h_reco_t[:]
del h_gen_m[:]
del h_reco_m[:]
del h_gen_l[:]
del h_reco_l[:]

## eta-eff
j=1
leg = ROOT.TLegend(0.7,0.7,0.9,0.9)
leg1 = ROOT.TLegend(0.5,0.35,0.7,0.55)
for init_i, init_loop in enumerate(init_cut):
    for eta_i, eta_loop in enumerate(plotvar_eta):
        title = eta_loop
        canvas = ROOT.TCanvas(title,title)
        logscale = False
        if eta_i==0:
            tcut = ''
        else:
            tcut = init_loop
        for i in mcfilelist:
            if not 'DY' in i:
                continue
            j=j+1
            rootfilename = i+".root"
            print rootfilename
            samplename = i.strip().split("_")[0]+init_config[init_i]
            tt = ROOT.TFile(rootfilename)
            tree = tt.h2mu.Get("tree2")
            # untill better way to get nentries
            tempdraw = eta_loop +" >> temp" +samplename
            tree.Draw(tempdraw)
            temphist = ROOT.gDirectory.Get("temp" +samplename)
            # untill better way to get nentries
            scale = mcfilelist[i]*datalumi / temphist.GetEntries()
            print mcfilelist[i]
            histo = copy.deepcopy(hist_maker(samplename, title, bin_set_l[1], x_name_l[1], y_name_l[1], tree, eta_loop, tcut))
            histo.SetLineColor(j)
            histo.Scale(scale)
            TypeOfMuons(eta_i, init_i, h_gen_t, h_reco_t, h_gen_m, h_reco_m, h_gen_l,h_reco_l,histo)
            print histo
            tt.Close()
#### plot style ####
st_marker = [20,25]
st_color = [2,6,7,8,9]
st_gen_tml = [h_gen_t,h_gen_m,h_gen_l]
st_reco_tml = [h_reco_t,h_reco_m,h_reco_l]
st_draw = ["","same","same"]
samplename = filenames[0].strip().split("_")[0]
####
for i in range(3):
    h_eff[1][i].SetStats(0)
    h_eff[1][i].Divide(st_reco_tml[i][0],st_gen_tml[i][0], 1.0, 1.0, "B")
    h_eff[1][i].SetLineColor(st_color[i])
    h_eff[1][i].Draw("E1%s"%st_draw[i])
    leg1.AddEntry(h_eff[1][i],samplename+init_config[i],"f")
leg1.Draw("same")
canvas.SaveAs(currentdir+saveddir+"/"+"eta_vs2"+".root")
canvas.SaveAs(currentdir+saveddir+"/"+"eta_vs2"+".eps")
canvas.SaveAs(currentdir+saveddir+"/"+"eta_vs2"+".png")
del leg1

## resoultion

for hist_i,i in enumerate(mcfilelist):
    for rsl_i, rsl_loop in enumerate(plotvar_rsl):
        title = rsl_loop
        h_rsl = [] 
        leg = ROOT.TLegend(0.6,0.6,0.9,0.9)
        leg1 = ROOT.TLegend(0.5,0.35,0.7,0.55)
        logscale = False
        j=1
        for init_i, init_loop in enumerate(init_cut):
            tcut=init_loop
            j=j+1
            rootfilename = i+".root"
            print rootfilename
            samplename = i.strip().split("_")[0]+init_config[init_i]
            tt = ROOT.TFile(rootfilename)
            tree = tt.h2mu.Get("tree2")
            # untill better way to get nentries
            tempdraw = rsl_loop +" >> temp" +samplename
            tree.Draw(tempdraw)
            temphist = ROOT.gDirectory.Get("temp" +samplename)
            # untill better way to get nentries
            scale = mcfilelist[i]*datalumi / temphist.GetEntries()
            print mcfilelist[i]
            histo = copy.deepcopy(hist_maker(samplename, title, bin_set_l[2], x_name_l[2], y_name_l[2], tree, rsl_loop, tcut))
            histo.SetLineColor(j)
            histo.Scale(scale)
            histo.Fit("gaus","0","",-0.05,0.05)
            fitresult = ROOT.TVirtualFitter.GetFitter()
            sig = fitresult.GetParameter(2)
            print sig
            h_rsl.append(histo)
            leg.AddEntry(histo, samplename+", %4.4f%%"%(sig*100),"f")
            print histo
            tt.Close()
        #h_eff[2][hist_i].SetLabelOffset(0.5,"X")
        canvas = ROOT.TCanvas(title,title)
        h_eff[2][hist_i].Add(h_rsl[0])
        h_eff[2][hist_i].SetLineColor(2)
        h_eff[2][hist_i].Sumw2(False)
        h_eff[2][hist_i].SetStats(0)
        h_eff[2][hist_i].Draw()
        h_rsl[1].Sumw2(False)
        h_rsl[1].Draw("same")
        h_rsl[2].Sumw2(False)
        h_rsl[2].Draw("same")
        leg.Draw("same")
        canvas.SaveAs(currentdir+saveddir+"/"+title+"_"+i.strip().split("_")[0]+"_test3"+".root")
        canvas.SaveAs(currentdir+saveddir+"/"+title+"_"+i.strip().split("_")[0]+"_test3"+".eps")
        canvas.SaveAs(currentdir+saveddir+"/"+title+"_"+i.strip().split("_")[0]+"_test3"+".png")
        del leg
         
    #j=j+1

