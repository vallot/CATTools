import os
import ROOT,copy
ROOT.gROOT.SetBatch(True)

def hist_maker(name, title, bin_set, x_name, y_name, tr, br, cut):
    hist = ROOT.TH1F(name, title, bin_set[0], bin_set[1], bin_set[2])
    #hist.GetXaxis().SetTitle(x_name)
    #hist.GetYaxis().SetTitle(y_name)
    #hist.SetLineColor(color)
    #hist.Sumw2()
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
for i in filelist2:#replace 'filelist2' replace to 'filelist1'
    if 'DYJets' in i:filenames[0]=i
    if 'TTJets' in i:filenames[1]=i
    if 'ZZ' in i:filenames[2]=i
    if 'WW' in i:filenames[3]=i
    if 'WZ' in i:filenames[4]=i
    if 'Double' in i:filenames[5]=i
    if 'Single' in i:filenames[6]=i
'''
mcfilelist = {'mc_DYJets_merged':6025.2 ,
              'mc_TTJets_mad_merged':831.8,
              'mc_ZZ_merged':31.8,
              'mc_WW_merged':65.9,
              'mc_WZ_merged':118.7}
rdfilelist = [#'DoubleMuon',
              'data_MuonEG_merged']
'''
mcfilelist = {filenames[0]:6025.2,
              filenames[1]:831.8,
              filenames[2]:31.8,
              filenames[3]:65.9,
              filenames[4]:118.7}
rdfilelist = [filenames[5],
              filenames[6]]
mcfilelist_order = sorted(mcfilelist.keys(), key=lambda key_value: key_value[0])
print mcfilelist
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
init_cuts = ["(step == 2 && isTight)","(step == 2 && isMedium)"]
whatiscut = ["_tight","_medium"]

os.chdir(resultdir)
for ps,init_cut in enumerate(init_cuts):
    for plot in range(4):
        plotvar = "ll_m"
        title = plotvar
        tcut = init_cut
        bin_set = [250, 0, 250]
        x_name = "a.u."
        y_name = "M [GeV]"
        if plot == 1:
          title = plotvar+"_"+jetcat[4]
          tcut = init_cut+"*"+jetcat_cut[4]
        if plot == 2:
          title = plotvar+"_"+jetcat[5]
          tcut = init_cut+"*"+jetcat_cut[5]
        if plot == 3:
          title = plotvar+"_"+jetcat[6]
          tcut = init_cut+"*"+jetcat_cut[6]
        hs = ROOT.THStack(plotvar,plotvar)
        h_rd = ROOT.TH1F(plotvar, plotvar, bin_set[0], bin_set[1], bin_set[2])
        h_rd.SetMarkerStyle(20)
        leg = ROOT.TLegend(0.7,0.7,0.9,0.9)
        leg.AddEntry(h_rd,"Data","p")
        logscale = False
        j=1
        for i in mcfilelist:
            j=j+1
            rootfilename = i+".root"
            print rootfilename
            samplename = i.strip().split("_")[0]
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
            histo.SetFillColor(j)
            histo.Scale(scale)
            if histo.GetMaximum()>2000:
                logscale = True
            hs.Add(histo)
            leg.AddEntry(histo,samplename,"f")
            tt.Close()
        for i in rdfilelist:
            rootfilename = i+".root"
            print rootfilename
            samplename = i.strip().split("_")[0]
            tt = ROOT.TFile(rootfilename)
            tree = tt.h2mu.Get("tree")
            histo = copy.deepcopy(hist_maker(samplename, title, bin_set, x_name, y_name, tree, plotvar, tcut))
            if histo.GetMaximum()>2000:
                logscale = True
            h_rd.Add(histo)
            tt.Close()
            

        canvas = ROOT.TCanvas(title,title)
        hs.Draw()
        h_rd.Draw("psame")
        leg.Draw("same")
        if logscale:
            canvas.SetLogy()
        canvas.SaveAs(currentdir+saveddir+"/"+title+whatiscut[ps]+".root")
        canvas.SaveAs(currentdir+saveddir+"/"+title+whatiscut[ps]+".png")

        del leg

    for jet01 in range(4):
      for gc in range(6):
        title = plotvar+"_"+jetcat[jet01]+"_"+jetcat_GC[gc]
        tcut = init_cut+"*"+jetcat_cut[jet01]+"*"+jetcat_GC_cut[gc]
        hs = ROOT.THStack(plotvar,plotvar)
        h_rd = ROOT.TH1F(plotvar, plotvar, bin_set[0], bin_set[1], bin_set[2])
        h_rd.SetMarkerStyle(20)
        leg = ROOT.TLegend(0.7,0.7,0.9,0.9)
        leg.AddEntry(h_rd,"Data","p")
        logscale=False    

        j=1
        for i in mcfilelist:
            j=j+1
            rootfilename = i+".root"
            print rootfilename
            samplename = i.strip().split("_")[0]
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
            histo.SetFillColor(j)
            histo.Scale(scale)
            if histo.GetMaximum()>2000:
                logscale = True
            hs.Add(histo)
            leg.AddEntry(histo,samplename,"f")
            tt.Close()

        for i in rdfilelist:
            rootfilename = i+".root"
            print rootfilename
            samplename = i.strip().split("_")[0]
            tt = ROOT.TFile(rootfilename)
            tree = tt.h2mu.Get("tree")
            histo = copy.deepcopy(hist_maker(samplename, title, bin_set, x_name, y_name, tree, plotvar, tcut))
            if histo.GetMaximum()>2000:
                logscale = True
            h_rd.Add(histo)
            tt.Close()        

        canvas = ROOT.TCanvas(title,title)
        hs.Draw()
        h_rd.Draw("psame")
        leg.Draw("same")
        if logscale:
            canvas.SetLogy()
        canvas.SaveAs(currentdir+saveddir+"/"+title+whatiscut[ps]+".root")
        canvas.SaveAs(currentdir+saveddir+"/"+title+whatiscut[ps]+".png")

        del leg

#dilepmass = ROOT.gDirectory.Get("dilepmass")
