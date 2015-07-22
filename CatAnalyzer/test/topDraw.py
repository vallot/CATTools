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

plot = 1
datalumi = 55.926
#datalumi = 49.502
plotvar = "ll_mass"
tcut = "step >= 1 && channel == 0"
bin_set = [38, 20, 400]
x_name = "a.u."
y_name = "M [GeV]"
if plot == 2:
    plotvar = "njet"
    tcut = "step >= 2 && channel == 0"
    bin_set = [10, 0, 10]
if plot == 3:
    plotvar = "nbjet"
    tcut = "step >= 4 && channel == 0"
    bin_set = [10, 0, 10]
    

mcfilelist = {'DYJetsToLL_M-50':6025.2 ,
              'TTJets_TuneCUETP8M1_13TeV-madgraphMLM':831.8,
              'ZZ_TuneCUETP8M1_13TeV':31.8,
              'WW_TuneCUETP8M1_13TeV':65.9,
              'WZ_TuneCUETP8M1_13TeV':118.7,
              'WJetsToLNu_TuneCUETP8M1_13TeV':61526.7,
              'ST_tW_antitop_5f_inclusiveDecays_13TeV':35.6,
              'ST_tW_top_5f_inclusiveDecays_13TeV':35.6,
              }

rdfilelist = [#'DoubleMuon',
               #'SingleMuon'
               #'DoubleEG'
               'SingleElectron'
               #'MuonEG'
               ]

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
    histo = copy.deepcopy(hist_maker(samplename, plotvar, bin_set, x_name, y_name, tree, plotvar, tcut))
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
    histo = copy.deepcopy(hist_maker(samplename, plotvar, bin_set, x_name, y_name, tree, plotvar, tcut))
    h_rd.Add(histo)
    
canvas = ROOT.TCanvas(plotvar,plotvar)
hs.Draw()
h_rd.Draw("psame")
leg.Draw("same")
canvas.SaveAs(plotvar+".root")
canvas.SaveAs(plotvar+".eps")

#dilepmass = ROOT.gDirectory.Get("dilepmass")

