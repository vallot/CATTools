import ROOT,copy
import array
import sys

channels = ["", " && channel == 0", " && channel == 1"," && channel == 2"]
filename = ["", "_mm", "_em","_ee"]

tt = ROOT.TFile("TTJets_TuneCUETP8M1_13TeV-madgraphMLM.root")
tree = tt.ttll.Get("tree")
axisname = ["Hadron", "Muon", "Elec", "Tau"]

for i, ch in enumerate(channels):
	canvas = ROOT.TCanvas("parton_ch"+ch, "parton_ch"+ch, 600, 600)
	tree.Draw("parton_mode1:parton_mode2 >> histo1", "step >= 1"+ch, "text")
	histo1 = ROOT.gDirectory.Get("histo1")
	histo = histo1.Clone()

	for j in range(7):
		histo.SetBinContent(2, j+1, histo1.GetBinContent(1, j+1))
		histo.SetBinContent(j+1, 2, histo1.GetBinContent(j+1, 1))
	histo.SetBinContent(2, 2, histo1.GetBinContent(1, 1))
	for j in range(4):
		histo.SetBinContent(5, j+1, histo.GetBinContent(5, j+1)+histo.GetBinContent(6, j+1)+histo.GetBinContent(7, j+1))
		histo.SetBinContent(j+1, 5, histo.GetBinContent(j+1, 5)+histo.GetBinContent(j+1, 6)+histo.GetBinContent(j+1, 7))

	lastcell = 0.
	for j in range(3):
		for z in range(3):
			lastcell += histo1.GetBinContent(5+j, 5+z)
	histo.SetBinContent(5, 5, lastcell)

	histo.GetXaxis().SetRangeUser(1, 5)
	histo.GetYaxis().SetRangeUser(1, 5)
	histo.GetXaxis().SetLabelSize(0.05)
	histo.GetYaxis().SetLabelSize(0.05)
	
	for j in range(4):
		histo.GetXaxis().SetBinLabel(j+2, axisname[j])
		histo.GetYaxis().SetBinLabel(j+2, axisname[j])

	histo.SetTitle("parton mode"+filename[i])
	histo.SetStats(0)
	histo.Draw("text")

	canvas.SetGrid()
	canvas.SaveAs("partonmode_analysis"+filename[i]+".png")

