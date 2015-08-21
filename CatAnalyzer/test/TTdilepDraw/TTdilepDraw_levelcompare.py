import ROOT,copy
import array
import sys

def hist_maker(name, title, tr, br, cut):
    hist = ROOT.TH2F(name, title, 7, 0, 7, 7, 0, 7)
    hist.SetStats(0)
    tr.Project(name, br, cut)
    return hist

levels = ["gen", "parton", "pseudoTop"]

tt = ROOT.TFile("ttpowheg.root")
#tt = ROOT.TFile("TTJets_TuneCUETP8M1_13TeV-madgraphMLM.root")
tree = tt.ttll.Get("tree")
print tree.GetEntries()
axisname = ["Hadron", "Muon", "Elec", "Tau"]

for i, lv in enumerate(levels):
	canvas = ROOT.TCanvas("compare"+lv, "compare"+lv, 600, 600)

	plotvar1 = "{0}_mode1:{0}_mode2".format(lv)
	plotvar2 = "{0}_mode2:{0}_mode1".format(lv)
	cut = "lepinPhase == 1 && jetinPhase == 1"
	histo = copy.deepcopy(hist_maker("histo", plotvar1, tree, plotvar1, cut))
	histo_inverted = copy.deepcopy(hist_maker("histo_inverted", plotvar2, tree, plotvar2, cut))

	histo.Add(histo_inverted)
	
	########### rebinning(there's gotta be better way...)
	for j in range(7):
		histo.SetBinContent(2, j+1, histo.GetBinContent(1, j+1))
		histo.SetBinContent(j+1, 2, histo.GetBinContent(j+1, 1))
	histo.SetBinContent(2, 2, histo.GetBinContent(1, 1))
	for j in range(4):
		histo.SetBinContent(5, j+1, histo.GetBinContent(5, j+1)+histo.GetBinContent(6, j+1)+histo.GetBinContent(7, j+1))
		histo.SetBinContent(j+1, 5, histo.GetBinContent(j+1, 5)+histo.GetBinContent(j+1, 6)+histo.GetBinContent(j+1, 7))

	lastcell = 0.
	for j in range(3):
		for z in range(3):
			lastcell += histo.GetBinContent(5+j, 5+z)
	histo.SetBinContent(5, 5, lastcell)

	for j in range(2, 6):
		histo.SetBinContent(j, j, histo.GetBinContent(j, j)/2.)
	########### rebinning(there's gotta be better way...)

	for j in range(4):
		histo.GetXaxis().SetBinLabel(j+2, axisname[j])
		histo.GetYaxis().SetBinLabel(j+2, axisname[j])
	histo.GetXaxis().SetRangeUser(1, 5)
	histo.GetYaxis().SetRangeUser(1, 5)
	
	histo.SetStats(0)
	#histo.Scale(1/float(histo.GetEntries())*100)
	#ROOT.gStyle.SetPaintTextFormat("3.4f%%")
	histo.Draw("text")

	canvas.SetGrid()
	canvas.SaveAs("compare_"+lv+".png")
