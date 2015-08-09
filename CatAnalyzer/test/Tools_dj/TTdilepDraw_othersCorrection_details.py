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
bin_set = [10, 20, 400]
x_name = "Mll"
chmatchcut = []
tautolep1cut = []
tautolep2cut = []
totalcut = "step >=1"
ttcut = "step >=1 && parton_channel == 3"
chmatchcut.append("step >=1 && channel == 0 && (parton_mode1 == 2 && parton_mode2 == 2)")
chmatchcut.append("step >=1 && channel == 1 && ((parton_mode1 == 3 && parton_mode2 == 2) || (parton_mode1 == 2 && parton_mode2 == 3))")
chmatchcut.append("step >=1 && channel == 2 && (parton_mode1 == 3 && parton_mode2 == 3)")
tautolep1cut.append("step >=1 && channel == 0 && ((parton_mode1 == 5 && parton_mode2 == 2) || (parton_mode1 == 2 && parton_mode2 == 5))")
tautolep1cut.append("step >=1 && channel == 1 && ((parton_mode1 == 2 && parton_mode2 == 6) || (parton_mode1 == 6 && parton_mode2 == 2))")
tautolep1cut.append("step >=1 && channel == 1 && ((parton_mode1 == 3 && parton_mode2 == 5) || (parton_mode1 == 5 && parton_mode2 == 3))")
tautolep1cut.append("step >=1 && channel == 2 && ((parton_mode1 == 6 && parton_mode2 == 3) || (parton_mode1 == 3 && parton_mode2 == 6))")
tautolep2cut.append("step >=1 && channel == 0 && (parton_mode1 == 5 && parton_mode2 == 5)")
tautolep2cut.append("step >=1 && channel == 1 && ((parton_mode1 == 5 && parton_mode2 == 6) || (parton_mode1 == 6 && parton_mode2 == 5))")
tautolep2cut.append("step >=1 && channel == 2 && (parton_mode1 == 6 && parton_mode2 == 6)")
#input values here

y_name = "a.u. (%)"
#y_name = "Number of Events"
tt = ROOT.TFile("TTJets_TuneCUETP8M1_13TeV-madgraphMLM.root")
tree = tt.ttll.Get("tree")
scale = crosssection*datalumi/tree.GetEntries()
canvas = ROOT.TCanvas(plotvar, plotvar, 600, 600)
leg = ROOT.TLegend(0.5, 0.15, 0.88, 0.4)

histo1 = copy.deepcopy(histMaker(tree, "total", totalcut, plotvar, bin_set))
totalEntry = histo1.GetEntries()
totalIntegral = histo1.Integral()

histo2 = copy.deepcopy(histMaker(tree, "tt_before", ttcut, plotvar, bin_set))
histo2.SetLineColor(17)
histo2.SetFillColor(17)
histo2.SetFillStyle(3002)
ttEntry = histo2.GetEntries()
ttIntegral = histo2.Integral()
histo2.Divide(histo1)
histo2.Scale(100)
percent1 = str(round((ttEntry/float(totalEntry))*100, 2))+"%"
leg.AddEntry(histo2, "ttbar_before"+"  "+percent1, "f")

hs_corrected = ROOT.THStack(plotvar, plotvar)
hs_chmatch = ROOT.THStack(plotvar, plotvar+"chmatch")
hs_tautolep1 = ROOT.THStack(plotvar, plotvar+"tautolep1")
hs_tautolep2 = ROOT.THStack(plotvar, plotvar+"tautolep2")
hs_corrected2 = ROOT.THStack(plotvar, plotvar+"_corrected")

for i, cut in enumerate(chmatchcut):
	histo = copy.deepcopy(histMaker(tree, str(i), cut, plotvar, bin_set))
	hs_chmatch.Add(histo)

tautolep1cut += chmatchcut
for i, cut in enumerate(tautolep1cut):
	histo = copy.deepcopy(histMaker(tree, str(i), cut, plotvar, bin_set))
	hs_tautolep1.Add(histo)

tautolep2cut += tautolep1cut
for i, cut in enumerate(tautolep2cut):
	histo = copy.deepcopy(histMaker(tree, str(i), cut, plotvar, bin_set))
	hs_tautolep2.Add(histo)

tautolep = tautolep2cut
for i, cut in enumerate(tautolep):
	histo = copy.deepcopy(histMaker(tree, str(i), cut+" && inPhase != 0", plotvar, bin_set))
	hs_corrected2.Add(histo)

chmatch = hs_chmatch.GetStack().Last()
tautolep1 = hs_tautolep1.GetStack().Last()
tautolep2 = hs_tautolep2.GetStack().Last()
last = hs_corrected2.GetStack().Last()

chmatch.SetLineColor(15)
chmatch.SetLineStyle(2)
chmatch.SetLineWidth(2)
percent2 = str(round((chmatch.GetEntries()/float(totalEntry))*100, 2))+"%"
chmatch.Divide(histo1)
chmatch.Scale(100)
leg.AddEntry(chmatch, "channel matched  "+percent2, "l")
tautolep1.SetLineColor(13)
tautolep1.SetLineStyle(3)
tautolep1.SetLineWidth(2)
percent3 = str(round((tautolep1.GetEntries()/float(totalEntry))*100, 2))+"%"
tautolep1.Divide(histo1)
tautolep1.Scale(100)
leg.AddEntry(tautolep1, "1 tau to lep  "+percent3, "l")
tautolep2.SetLineColor(12)
#tautolep2.SetLineStyle(4)
tautolep2.SetLineWidth(2)
percent4 = str(round((tautolep2.GetEntries()/float(totalEntry))*100, 2))+"%"
tautolep2.Divide(histo1)
tautolep2.Scale(100)
leg.AddEntry(tautolep2, "2 tau to lep  "+percent4, "l")
last.SetLineColor(1)
last.SetLineWidth(4)
last.Scale(100)
percent5 = str(round((last.GetEntries()/float(totalEntry))*100, 2))+"%"
last.Divide(histo1)
leg.AddEntry(last, "in phase   "+percent5, "l")

histo2.SetTitle("ttbar_corrected_ll_m")
yaxis = histo2.GetYaxis()
xaxis = histo2.GetXaxis()
yaxis.SetRangeUser(70., 100.)
yaxis.SetTitle(y_name)
xaxis.SetTitle(x_name)
xaxis.SetTitleSize(0.044)
xaxis.SetLabelSize(0.03)
yaxis.CenterTitle()
yaxis.SetTitleSize(0.04)
yaxis.SetLabelSize(0.03)
histo2.Draw()
tautolep2.Draw("same")
tautolep1.Draw("same")
chmatch.Draw("same")
last.Draw("same")
leg.Draw()

canvas.SaveAs("tt_correction_"+plotvar+".png")

