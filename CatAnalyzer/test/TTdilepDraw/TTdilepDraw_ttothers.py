import ROOT,copy
import array
import sys

def histMaker(tree, name, tcut, plotvar, bin_set):
	hist = ROOT.TH1F(name, plotvar, bin_set[0], bin_set[1], bin_set[2])
	hist.SetStats(0)
	tree.Project(name, plotvar, tcut)
	return hist

########## user input start ############
plotvar = "ll_m"
x_name = "Mll"
bin_set = [20, 20, 400]
tcuts = []
tcut_total = "step == 5"
tcut_signal = "step == 5 && parton_channel == 3"
tcut_others = "step == 5 && parton_channel != 3"
tcuts.append("step == 5 && ((parton_channel == 2 && (parton_mode1 > 3 || parton_mode2 > 3)) || (parton_channel == 1 && parton_mode1 > 3 && parton_mode2 > 3))")
tcuts.append("step == 5 && ((parton_channel == 2 && (parton_mode1 == 0 || parton_mode2 == 0)) || (parton_channel == 1 && ((parton_mode1 == 0 && parton_mode2 != 0)|| (parton_mode2 == 0 && parton_mode1 != 0))))")
tcuts.append("step == 5 && parton_channel == 1 && parton_mode1 == 0 && parton_mode2 == 0")
tcuts.append("step == 5 && parton_channel == 3 && lepinPhase == 0")
name = ["dilep_tau", "1lep + 1jet", "full hadron", "out of phase"]
channels = ["", " && channel == 0", " && channel == 1"," && channel == 2"]
filename = ["", "mm_", "em_","ee_"]
########### user input end #############
 
y_name = "Number of Events"
tt = ROOT.TFile("TT_TuneCUETP8M1_13TeV-powheg.root")
tree = tt.ttll.Get("tree")

for j, ch in enumerate(channels):
	lst = []
	canvas = ROOT.TCanvas(plotvar+ch, plotvar+ch, 600, 600)
	canvas.SetLogy()
	leg = ROOT.TLegend(0.58, 0.6, 0.88, 0.88)

	total = copy.deepcopy(histMaker(tree, "total", tcut_total, plotvar, bin_set))
	tmp = ROOT.TH1F("tmp","tmp",380,20,400)
	for q in range(80):
		tmp.Fill(20)
	tmp.SetTitle("TTbar others Fraction_s5")
	tmp.GetXaxis().SetTitle("Mll")
	tmp.GetYaxis().SetTitle("%")
	tmp.SetLineColor(0)
	tmp.SetStats(0)
	tmp.Draw()

	h_signal = copy.deepcopy(histMaker(tree, "ttsignal"+filename[j], tcut_signal+ch, plotvar, bin_set))
	num = str(round(h_signal.GetEntries()/float(total.GetEntries())*100, 2))+"%"
	h_signal.SetLineColor(1)
	h_signal.SetLineWidth(4)
	h_signal.Divide(total)
	h_signal.Scale(100)
	h_signal.Draw("same")
#	h_signal.GetYaxis().SetRangeUser(0, 2)
	leg.AddEntry(h_signal,"ttsignal"+"  "+num, "l")

	h_others = copy.deepcopy(histMaker(tree, "ttothers"+filename[j], tcut_others+ch, plotvar, bin_set))
	num = str(round(h_others.GetEntries()/float(total.GetEntries())*100, 2))+"%"
	h_others.SetLineColor(1)
	h_others.SetLineWidth(2)
	h_others.SetLineStyle(2)
	h_others.Divide(total)
	h_others.Scale(100)
	h_others.Draw("same")
	leg.AddEntry(h_others,"ttothers"+" "+num, "l")

	for i, tcut in enumerate(tcuts):
		histo = copy.deepcopy(histMaker(tree, name[i]+filename[j], tcut+ch, plotvar, bin_set))
		num = str(round(histo.GetEntries()/float(total.GetEntries())*100, 2))+"%"
		histo.SetLineColor(i+2)
		histo.SetLineWidth(2)
		histo.Divide(total)
		histo.Scale(100)
		lst.append(histo)
		leg.AddEntry(histo, name[i]+"  "+str(num), "l")

	for z in lst:
		z.Draw("same")
	leg.Draw()

	canvas.SaveAs("ttothersRatio_s5_"+filename[j]+plotvar+".png")

