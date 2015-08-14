import ROOT,copy
import array
import sys

def hist_maker(name, title, bin_set, x_name, y_name, tr, br, cut):
    if title == "ll_m":
        massbin = array.array('d', [20, 30, 45, 75, 105, 125, 165, 260, 400])
        hist = ROOT.TH1F(name, title, 8, massbin)
    else: hist = ROOT.TH1F(name, title, bin_set[0], bin_set[1], bin_set[2])
    hist.SetStats(0)
    tr.Project(name, br, cut)
    return hist

plotvar = sys.argv[1]

datalumi = 40.2

mcfilelist = ['TTJets_TuneCUETP8M1_13TeV-madgraphMLM',
              'WJetsToLNu_TuneCUETP8M1_13TeV',
              'ST_tW_antitop_5f_inclusiveDecays_13TeV',
              'ST_tW_top_5f_inclusiveDecays_13TeV',
              'ZZ_TuneCUETP8M1_13TeV',
              'WW_TuneCUETP8M1_13TeV',
              'WZ_TuneCUETP8M1_13TeV',
              'DYJetsToLL_M-50']
mcval_l = [831.8, 61526.7, 35.6, 35.6, 31.8, 65.9, 118.7, 6025.2]
color_l = [2, 0, 3, 5, 7, 7, 7, 4]

rdfilelist = ['MuonEG']

plotvar_l = ["ll_m", "njet", "MET", "nbjet"]
bin_set_l = [[38, 20, 400], [10, 0, 10], [30,0,300], [6, 0, 6]]
x_name_l = ["ElMu M(ll) [GeV]   ", "ElMu Jet Multiplicity   ", "ElMu Missing Et[GeV/c]   ", "ElMu b Jet Multiplicity   "]

bin_set = bin_set_l[plotvar_l.index(plotvar)]
x_name = x_name_l[plotvar_l.index(plotvar)]
y_name = "Number of Events"

tcut = "step >=1 && channel == 1"

hs = ROOT.THStack(plotvar, plotvar)
if plotvar == "ll_m":
    massbin = array.array('d', [20, 30, 45, 75, 105, 125, 165, 260, 400])
    h_rd = ROOT.TH1F(plotvar, plotvar, 8, massbin)
else: h_rd = ROOT.TH1F(plotvar, plotvar, bin_set[0], bin_set[1], bin_set[2])
h_rd.SetMarkerStyle(20)

leg = ROOT.TLegend(0.57,0.5,0.88,0.88)
leg.SetTextSize(0.035)
leg.SetTextFont(42)
leg.SetLineColor(0)
leg.SetFillColor(0)

mchist_l = []
samplenames = []
for i, f in enumerate(mcfilelist):
    rootfilename = f+".root"
    print rootfilename
    samplename = f.strip().split("_")[0]
    samplenames.append(samplename)
    tt = ROOT.TFile(rootfilename)
    tree = tt.ttll.Get("tree")
    """
	# untill better way to get nentries
    tempdraw = plotvar +" >> temp" +samplename
    tree.Draw(tempdraw)
    temphist = ROOT.gDirectory.Get("temp" +samplename)
    # untill better way to get nentries
    scale = mcval_l[i]*datalumi / temphist.GetEntries()
    """
    scale = mcval_l[i]*datalumi / tree.GetEntries()
    if samplename == "TTJets":
        histo1 = copy.deepcopy(hist_maker(samplename, plotvar, bin_set, x_name, y_name, tree, plotvar, tcut+" && parton_channel == 3"))
        histo2 = copy.deepcopy(hist_maker(samplename+" others", plotvar, bin_set, x_name, y_name, tree, plotvar, tcut+" && parton_channel != 3"))
        histo1.SetFillColor(color_l[i])
        histo2.SetFillColor(6)
        histo1.Scale(scale)
        histo2.Scale(scale)
        num = histo1.GetEntries()+histo2.GetEntries()
        mchist_l.append(histo1) 
        mchist_l.append(histo2) 
    else:
        histo = copy.deepcopy(hist_maker(samplename, plotvar, bin_set, x_name, y_name, tree, plotvar, tcut))
        histo.SetFillColor(color_l[i])
        histo.Scale(scale)
        mchist_l.append(histo)

for i in rdfilelist:
    rootfilename = i+".root"

    print rootfilename
    samplename = i.strip().split("_")[0]
    tt = ROOT.TFile(rootfilename)
    tree = tt.ttll.Get("tree")
    histo = copy.deepcopy(hist_maker(samplename, plotvar, bin_set, x_name, y_name, tree, plotvar, tcut))
    h_rd.Add(histo)

canvas = ROOT.TCanvas(plotvar,plotvar, 600, 600)
#canvas.SetLogy(1)
#hs.SetMaximum(300)

for i in mchist_l:
    hs.Add(i)

leg.AddEntry(h_rd,"Data","lp")
leg.AddEntry(mchist_l[8], samplenames[7], "f")
leg.AddEntry(mchist_l[5], "WW/WZ/ZZ", "f")
leg.AddEntry(mchist_l[4], "Single top W", "f")
leg.AddEntry(mchist_l[3], "Single topbar W", "f")
leg.AddEntry(mchist_l[2], samplenames[1], "f")
leg.AddEntry(mchist_l[1], samplenames[0]+" others", "f")
leg.AddEntry(mchist_l[0], samplenames[0]+str(num), "f")
      
hs.Draw()
hs.SetTitle("")
hs.GetXaxis().SetTitle(x_name)
hs.GetXaxis().SetTitleSize(0.044)
hs.GetXaxis().SetLabelSize(0.03)
hs.GetYaxis().SetTitle(y_name)
hs.GetYaxis().CenterTitle()
hs.GetYaxis().SetTitleSize(0.044)
hs.GetYaxis().SetLabelSize(0.03)
h_rd.Draw("epsame")

tex1 = ROOT.TLatex(0.56,0.92,"%.1f pb^{-1},   #sqrt{s} = 13 TeV"%(datalumi))
tex1.SetNDC()
tex1.SetTextFont(42)
tex1.SetTextSize(0.042)
tex1.Draw()

tex2 = ROOT.TLatex(0.1,0.92,"CMS Preliminary")
tex2.SetNDC()
tex2.SetTextFont(42)
tex2.SetTextSize(0.042)
tex2.Draw("same")

leg.Draw("same")
#canvas.SaveAs(plotvar+".root")
canvas.SaveAs(plotvar+".png")

