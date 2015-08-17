import ROOT,copy
import array
import sys

def histMaker(tree, name, tcut, plotvar, bin_set):
	hist = ROOT.TH1F(name, plotvar, bin_set[0], bin_set[1], bin_set[2])
	hist.SetStats(0)
	tree.Project(name, plotvar, tcut)
	return hist

datalumi = 40.2
crosssection = 831.8

#input values here
plotvar = "ll_m"
bin_set = [38, 20, 400]
x_name = "Mll"
tcuts = []
tcut0 = "step >=1 && parton_channel != 3"
tcuts.append("step >=1 && parton_channel == 2 && (parton_mode1 > 3 || parton_mode2 > 3)")
tcuts.append("step >=1 && parton_channel == 2 && (parton_mode1 == 0 || parton_mode2 == 0)")
tcuts.append("step >=1 && parton_channel == 1 && ((parton_mode1 == 0 && parton_mode2 != 0)|| (parton_mode2 == 0 && parton_mode1 != 0))")
tcuts.append("step >=1 && parton_channel == 1 && parton_mode1 > 3 && parton_mode2 > 3")
tcuts.append("step >=1 && parton_channel == 1 && parton_mode1 == 0 && parton_mode2 == 0")
tcuts.append("step >=1 && parton_channel == 3 && inPhase == 0")
name = ["lep+tau", "lep+hadron", "tau+hadron", "full tau", "full hadron", "out of phase"]
channels = ["", " && channel == 0", " && channel == 1"," && channel == 2"]
filename = ["", "mm_", "em_","ee_"]
#input values here

y_name = "Number of Events"
tt = ROOT.TFile("TTJets_TuneCUETP8M1_13TeV-madgraphMLM.root")
tree = tt.ttll.Get("tree")
scale = crosssection*datalumi/tree.GetEntries()
canvas = ROOT.TCanvas(plotvar, plotvar, 600, 600)

for j, ch in enumerate(channels):
	hs = ROOT.THStack(plotvar, filename[j]+plotvar)
	leg = ROOT.TLegend(0.55, 0.55, 0.88, 0.88)

	histl = copy.deepcopy(histMaker(tree, "ttothers"+filename[j], tcut0+ch, plotvar, bin_set))
	histl.SetLineColor(1)
	histl.SetLineWidth(3)
	histl.SetLineStyle(2)
	histl.Scale(scale)
	histlname = "ttothers"
	leg.AddEntry(histl, histlname, "lp")

	for i, tcut in enumerate(tcuts):
		histo = copy.deepcopy(histMaker(tree, name[i]+filename[j], tcut+ch, plotvar, bin_set))
		num= int(histo.GetEntries())
		histo.SetLineColor(i+2)
		histo.SetFillColor(i+2)
		histo.Scale(scale)
		hs.Add(histo)
		leg.AddEntry(histo, name[i]+"  "+str(num), "f")

#	hs.SetMaximum(0.4)
#	if j == 0:
#		hs.SetMaximum(1)
	hs.Draw()
	hs.GetXaxis().SetTitle(x_name)
	hs.GetXaxis().SetTitleSize(0.044)
	hs.GetXaxis().SetLabelSize(0.03)
	hs.GetYaxis().SetTitle(y_name)
	hs.GetYaxis().CenterTitle()
	hs.GetYaxis().SetTitleSize(0.04)
	hs.GetYaxis().SetLabelSize(0.03)
	histl.Draw("same")
	leg.Draw()
	canvas.SaveAs("ttothers_"+filename[j]+plotvar+".png")

