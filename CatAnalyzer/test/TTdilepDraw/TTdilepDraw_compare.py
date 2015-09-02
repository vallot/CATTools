import ROOT,copy
import array
import sys

channels = ["gen", "parton", "pseudoTop"]
filename = ["", "_mm", "_em","_ee"]

tt = ROOT.TFile("TT_TuneCUETP8M1_13TeV-powheg.root")
tree = tt.ttll.Get("tree")
print tree.GetEntries()
axisname = ["Hadron", "Muon", "Elec", "Tau"]

for i, ch in enumerate(channels):
	canvas = ROOT.TCanvas("compare"+ch, "compare"+ch, 600, 600)
	str1 = "{0}_mode1:{0}_mode2 >> {0}_histo1".format(ch)
	str2 = "{0}_mode2:{0}_mode1 >> {0}_histo2".format(ch)
	tree.Draw(str1, "", "text")
	tree.Draw(str2, "", "text")

	histo1 = ROOT.gDirectory.Get(ch+"_histo1")
	histo2 = ROOT.gDirectory.Get(ch+"_histo2")
	histo = histo1.DrawClone()

	histo.Add(histo2)
	
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

	for j in range(4):
		histo.GetXaxis().SetBinLabel(j+2, axisname[j])
		histo.GetYaxis().SetBinLabel(j+2, axisname[j])
	histo.GetXaxis().SetRangeUser(1, 5)
	histo.GetYaxis().SetRangeUser(1, 5)
	
	histo.SetStats(0)
	histo.Draw("text")

	canvas.SetGrid()
	canvas.SaveAs("compare_"+ch+".png")

