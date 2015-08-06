import ROOT,copy

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

mcfilelist = {'mc_DYJets_merged':6025.2 ,
              'mc_TTJets_mad_merged':831.8,
              'mc_ZZ_merged':31.8,
              'mc_WW_merged':65.9,
              'mc_WZ_merged':118.7}

rdfilelist = [#'DoubleMuon',
               'data_MuonEG_merged']
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
init_cut = "(step == 2 && isTight)"

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

    j=1
    for i in mcfilelist:
        j=j+1
        rootfilename = i+".root"
        print rootfilename
        samplename = i.strip().split("_")[0]
        tt = ROOT.TFile(rootfilename)
        tree = tt.ttll.Get("top")
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
        hs.Add(histo)
        leg.AddEntry(histo,samplename,"f")

    for i in rdfilelist:
        rootfilename = i+".root"
        print rootfilename
        samplename = i.strip().split("_")[0]
        tt = ROOT.TFile(rootfilename)
        tree = tt.ttll.Get("top")
        histo = copy.deepcopy(hist_maker(samplename, title, bin_set, x_name, y_name, tree, plotvar, tcut))
        h_rd.Add(histo)
        
    canvas = ROOT.TCanvas(plotvar,plotvar)
    hs.Draw()
    h_rd.Draw("psame")
    leg.Draw("same")
    canvas.SaveAs(title+".root")
    canvas.SaveAs(title+".eps")

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

    j=1
    for i in mcfilelist:
        j=j+1
        rootfilename = i+".root"
        print rootfilename
        samplename = i.strip().split("_")[0]
        tt = ROOT.TFile(rootfilename)
        tree = tt.ttll.Get("top")
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
        hs.Add(histo)
        leg.AddEntry(histo,samplename,"f")

    for i in rdfilelist:
        rootfilename = i+".root"
        print rootfilename
        samplename = i.strip().split("_")[0]
        tt = ROOT.TFile(rootfilename)
        tree = tt.ttll.Get("top")
        histo = copy.deepcopy(hist_maker(samplename, title, bin_set, x_name, y_name, tree, plotvar, tcut))
        h_rd.Add(histo)
        
    canvas = ROOT.TCanvas(title,title)
    hs.Draw()
    h_rd.Draw("psame")
    leg.Draw("same")
    canvas.SaveAs(title+".root")
    canvas.SaveAs(title+".eps")

    del leg

#dilepmass = ROOT.gDirectory.Get("dilepmass")
