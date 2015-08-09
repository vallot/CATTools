import ROOT,copy
import array
import sys

def histMaker(tree, name, tcut, plotvar, bin_set):
	hist = ROOT.TH1F(name, plotvar, bin_set[0], bin_set[1], bin_set[2])
	hist.SetStats(0)
	tree.Project(name, plotvar, tcut)
	return hist

datalumi = 7.3
crosssection = 831.8

#input values here
plotvar = "ll_m"
bin_set = [38, 20, 400]
x_name = "Mll"
tautolep = []
ttcut = "step >=1 && parton_channel == 3"
ttotherscut = "step >=1 && parton_channel != 3"
tautolep.append("step >=1 && channel == 0 && (parton_mode1 == 2 && parton_mode2 == 2) && inPhase != 0")#->tt
tautolep.append("step >=1 && channel == 0 && (parton_mode1 == 5 && parton_mode2 == 5) && inPhase != 0")#->tt
tautolep.append("step >=1 && channel == 0 && ((parton_mode1 == 5 && parton_mode2 == 2) || (parton_mode1 == 2 && parton_mode2 == 5)) && inPhase != 0")#->tt
tautolep.append("step >=1 && channel == 2 && (parton_mode1 == 3 && parton_mode2 == 3) && inPhase != 0")#->tt
tautolep.append("step >=1 && channel == 2 && (parton_mode1 == 6 && parton_mode2 == 6) && inPhase != 0")#->tt
tautolep.append("step >=1 && channel == 2 && ((parton_mode1 == 6 && parton_mode2 == 3) || (parton_mode1 == 3 && parton_mode2 == 6)) && inPhase != 0")#->tt
tautolep.append("step >=1 && channel == 1 && ((parton_mode1 == 3 && parton_mode2 == 2) || (parton_mode1 == 2 && parton_mode2 == 3)) && inPhase != 0")#->tt
tautolep.append("step >=1 && channel == 1 && ((parton_mode1 == 5 && parton_mode2 == 6) || (parton_mode1 == 6 && parton_mode2 == 5)) && inPhase != 0")#->tt
tautolep.append("step >=1 && channel == 1 && ((parton_mode1 == 2 && parton_mode2 == 6) || (parton_mode1 == 6 && parton_mode2 == 2)) && inPhase != 0")#->tt
tautolep.append("step >=1 && channel == 1 && ((parton_mode1 == 3 && parton_mode2 == 5) || (parton_mode1 == 5 && parton_mode2 == 3)) && inPhase != 0")#->tt
#input values here

y_name = "Number of Events"
tt = ROOT.TFile("TTJets_TuneCUETP8M1_13TeV-madgraphMLM.root")
tree = tt.ttll.Get("tree")
scale = crosssection*datalumi/tree.GetEntries()
canvas = ROOT.TCanvas(plotvar, plotvar, 600, 600)
leg = ROOT.TLegend(0.45, 0.65, 0.88, 0.88)

hs_before = ROOT.THStack(plotvar, plotvar+"_before")
hs_corrected = ROOT.THStack(plotvar, plotvar+"_corrected")

histo = copy.deepcopy(histMaker(tree, "tt_before", ttcut, plotvar, bin_set))
histo.SetLineColor(38)
histo.SetFillColor(38)
histo.SetFillStyle(3004)
histo.Scale(scale)
ttEntry = histo.GetEntries()
hs_before.Add(histo)

histo = copy.deepcopy(histMaker(tree, "ttothers_before", ttotherscut, plotvar, bin_set))
histo.SetLineColor(30)
histo.SetFillColor(30)
histo.SetFillStyle(3004)
histo.Scale(scale)
ttothersEntry = histo.GetEntries()
hs_before.Add(histo)
totalEntry = hs_before.GetStack().Last().GetEntries()
percent1 = str(round((ttEntry/float(totalEntry))*100))+"%"
percent2 = str(round((ttothersEntry/float(totalEntry))*100))+"%"
leg.AddEntry(histo, "tt_before"+"  "+percent1, "f")
leg.AddEntry(histo, "ttothers_before"+"  "+percent2, "f")

for i, cut in enumerate(tautolep):
	histo = copy.deepcopy(histMaker(tree, str(i), cut, plotvar, bin_set))
	histo.Scale(scale)
	hs_corrected.Add(histo)

last = hs_corrected.GetStack().Last()
last.SetLineColor(2)
last.SetLineStyle(2)
last.SetLineWidth(3)
last.SetTitle("ttothers_corrected")
percent3 = str(round((last.GetEntries()/float(totalEntry))*100))+"%"
leg.AddEntry(last, last.GetTitle()+"  "+percent3, "lp")

hs_before.SetTitle("ttothers_corrected_ll_m")
#hs_before.SetName("ttothers_corrected_ll_m")
hs_before.Draw()
last.Draw("same")
leg.Draw()
canvas.SaveAs("tt_correction_"+plotvar+".png")

