import ROOT,copy
import json
import os 
from array import array

def hist_maker(name, title, bins, x_name, y_name, tr, br, cut):
    hist = ROOT.TH1F(name, title, len(bins)-1, array('f', bins))
    hist.SetStats(0)
    tr.Project(name, br, cut)
    if "1" in title:
        hist2 = ROOT.TH1F(name+"2", title, len(bins)-1, array('f', bins))
        br = plotvar.replace("1", "2")
        tr.Project(name+"2", br, cut)
        hist.Add(hist2)
    return hist

def stack_design(hs, x_name, y_name):
	hs.SetTitle("")
	hs.GetXaxis().SetTitle(x_name)
	hs.GetXaxis().SetTitleSize(0.044)
	hs.GetXaxis().SetLabelSize(0.03)
	hs.GetYaxis().SetTitle(y_name)
	hs.GetYaxis().CenterTitle()
	hs.GetYaxis().SetTitleSize(0.048)
	hs.GetYaxis().SetLabelSize(0.03)
	hs.SetMaximum(hs.GetMaximum()+100)
	return hs
	
############ external infos ############
datalumi = 226.1

mcfilelist = ['TT_powheg',
              'WJets',
              'SingleTbar_tW',
              'SingleTop_tW',
              'ZZ',
              'WW',
              'WZ',
              'DYJets']
rdfilelist = ['MuonEG']

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset.json" % os.environ['CMSSW_BASE']))
color_l = [2, 6, 0, 3, 5, 7, 7, 7, 4]

massbin = [20, 30, 45, 75, 105, 125, 165, 260, 400]
njetbin = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
metbin = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300]
nbjetbin = [0, 1, 2, 3, 4, 5, 6]
topptbin = [0, 70, 140, 240, 400]
toprapibin = [0.0, 0.25, 0.65, 1.15, 1.7, 2.5]
ptbin = [20, 40, 70, 120, 180, 400]
etabin = [-2.5, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2.5]
jetptbin = [30, 50, 80, 130, 210, 500]
############ external infos ############

############ user input infos ############
step = 6#untill add b_step=6 to tree
############ user input infos ############

#control plot info
plotvar_l = ["ll_m", "njet", "MET", "nbjet"]
bins = [massbin, njetbin, metbin, nbjetbin] 
x_name = ["MuEl M(ll) [GeV/c^{2}]", "MuEl Jet Multiplicity", "MuEl Missing Et[GeV/c]", "MuEl b Jet Multiplicity"]
tcut = "step >= %d && channel == 1 && tri == 1 && filtered"%step

#analysis plot info
if step == 6:
	plotvar_l = ["top1_pt", "top1_rapi", "ttbar_pt", "ttbar_rapi", "lep1_pt", "lep1_eta", "jet1_pt", "jet1_eta"]
	bins = [topptbin, toprapibin, topptbin, toprapibin, ptbin, etabin, jetptbin, etabin]
	x_name = ["MuEl top p_{T} [GeV/c]", "MuEl top Rapidity", "MuEl TTbar p_{T} [GeV/c]", "MuEl TTbar Rapidity", "MuEl S6 lepton p_{T} [GeV/c]", "MuEl lepton #eta", "MuEl S6 Jet p_{T} [GeV/c]", "MuEl Jet #eta"]
	tcut = "step >= %d && channel == 1 && tri == 1 && filtered && top1_pt > 0"%5

y_name = "Number of Events"
signal_tcut = "&& (parton_channel == 2 && ((parton_mode1 == 1 && parton_mode2 == 2) || (parton_mode1 == 2 && parton_mode2 == 1)))"

#define texts
tex1 = ROOT.TLatex(0.52,0.91,"%.1f pb^{-1}, #sqrt{s} = 13 TeV 25ns"%(datalumi))
tex1.SetNDC()
tex1.SetTextFont(42)
tex1.SetTextSize(0.036)
tex1.SetTextColorAlpha(1, 0.)

tex2 = ROOT.TLatex(0.09,0.91,"CMS Preliminary")
tex2.SetNDC()
tex2.SetTextFont(61)
tex2.SetTextSize(0.044)

#making control plot start
for j, plotvar in enumerate(plotvar_l):
	mchist_l = []
	print plotvar

	#declare stack(for mc)
	hs = ROOT.THStack(plotvar, plotvar)

	#controlplot_MC
	for i, file in enumerate(mcfilelist):
		rootfilename = file+".root"
		tt = ROOT.TFile(rootfilename)
		tree = tt.ttll.Get("tree")
		for data in datasets:
			if data["name"] == file:
				scale = datalumi/data["lumi"]
		if file == "TT_powheg":
			histo1 = copy.deepcopy(hist_maker(file, plotvar, bins[j], x_name[j], y_name, tree, plotvar, tcut+signal_tcut))
			histo2 = copy.deepcopy(hist_maker(file+" others", plotvar, bins[j], x_name[j], y_name, tree, plotvar, tcut))
			histo2.Add(histo1, -1)
			histo1.Scale(scale)
			histo2.Scale(scale)
			mchist_l.append(histo1) 
			mchist_l.append(histo2)
		else:
			histo = copy.deepcopy(hist_maker(file, plotvar, bins[j], x_name[j], y_name, tree, plotvar, tcut))
			histo.Scale(scale)
			mchist_l.append(histo)

	#controlplot_Data
	for file in rdfilelist:
		rootfilename = file+".root"
		tt = ROOT.TFile(rootfilename)
		tree = tt.ttll.Get("tree")
		h_rd = copy.deepcopy(hist_maker(file, plotvar, bins[j], x_name[j], y_name, tree, plotvar, tcut))
		h_rd.SetMarkerStyle(20)

	#drawing canvas
	canvas = ROOT.TCanvas(plotvar, plotvar, 600, 600)

	for i, h in enumerate(mchist_l):
		h.SetFillColor(color_l[i])
		h.SetLineColor(1)
		hs.Add(h)

	hs.Draw()
	stack_design(hs, x_name[j], y_name)

	h_rd.Draw("epsame")
	tex1.Draw()
	#tex2.Draw("same")

	#drawing legend
	leg = ROOT.TLegend(0.5,0.5,0.82,0.88)
	leg.SetTextSize(0.035)
	leg.SetTextFont(42)
	leg.SetLineColor(0)
	leg.SetFillColor(0)
	leg.AddEntry(h_rd,"Data","lp")
	leg.AddEntry(mchist_l[8], mcfilelist[7], "f")
	leg.AddEntry(mchist_l[5], "WW/WZ/ZZ", "f")
	leg.AddEntry(mchist_l[4], mcfilelist[3], "f")
	leg.AddEntry(mchist_l[3], mcfilelist[2], "f")
	leg.AddEntry(mchist_l[2], mcfilelist[1], "f")
	leg.AddEntry(mchist_l[1], mcfilelist[0]+" others", "f")
	leg.AddEntry(mchist_l[0], mcfilelist[0], "f")
	leg.Draw("same")

	if "1" in plotvar:
		plotvar = plotvar.replace("1", "")
	canvas.SaveAs("s%d_%s.png"%(step, plotvar))


#print number in each steps
print
rootfilename = rdfilelist[0]+".root"
tree = tt.ttll.Get("tree")
for i in range(0,6):
	tree.Draw("ll_m >> hist", "step >=%d && channel == 1 && tri ==1 && filtered"%(i))
	hist = ROOT.gDirectory.Get("hist")
	print "s%d   %d"%(i, hist.GetEntries())
tree.Draw("ll_m >> hist", "step >=5 && channel == 1 && tri ==1 && filtered && top1_pt >= 0")
hist = ROOT.gDirectory.Get("hist")
print "s%d   %d"%(i, hist.GetEntries())


