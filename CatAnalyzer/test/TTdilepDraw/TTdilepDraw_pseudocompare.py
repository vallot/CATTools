import ROOT, copy

def hist_maker(name, title, tr, br, cut):
    hist = ROOT.TH2F(name, title, 4, -1, 3, 4, -1, 3)
    hist.SetStats(0)
    tr.Project(name, br, cut)
    return hist


mcfile = 'TTJets_TuneCUETP8M1_13TeV-powheg-pythia8'
rootfilename = mcfile+".root"
samplename = mcfile.strip().split("_")[0]
tt = ROOT.TFile(rootfilename)
tree = tt.ttll.Get("tree")

cnv = ROOT.TCanvas()
plotvar = "phasein_pseudolep:phasein_pseudojet"
cut = "pseudoTop_channel != 10000 && channel == 1 && parton_channel == 3 && ((parton_mode1 == 2 && parton_mode2 == 3) || (parton_mode1 == 3 && parton_mode2 == 2))"
histo = copy.deepcopy(hist_maker("histo", plotvar, tree, plotvar, cut))
histo.SetTitle("No pseudoTop cut")
#histo.SetTitle("pseudoTop em channel")
#histo.SetTitle("pseudoTop not em channel")

for i in range(4):
	histo.GetXaxis().SetBinLabel(i+1, str(i-1))
	histo.GetYaxis().SetBinLabel(i+1, str(i-1))
histo.GetXaxis().SetLabelSize(0.06)
histo.GetXaxis().SetTitleSize(0.045)
histo.GetXaxis().SetTitleOffset(1)
histo.GetXaxis().SetTitle("# of phase-in pseudoTop leptons")
histo.GetYaxis().SetLabelSize(0.06)
histo.GetYaxis().SetTitleSize(0.045)
histo.GetYaxis().SetTitleOffset(0.7)
histo.GetYaxis().SetTitle("# of phase-in pseudoTop jets")

histo.Scale(1/float(histo.GetEntries())*100)
ROOT.gStyle.SetPaintTextFormat("3.4f%%")
histo.SetMarkerSize(1.5)
histo.Draw("text")
cnv.SetGrid()
cnv.Print("pseudo.png")
