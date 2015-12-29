#!/usr/bin/env python
import ROOT, copy, getopt, sys
from CATTools.CatAnalyzer.histoHelper import *


try:
    opts, args = getopt.getopt(sys.argv[1:],"hdc:w:b:p:x:y:a:s:",["binning","plotvar","channel"])
except getopt.GetoptError:          
    print 'Usage : ./topDraw.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./topDraw.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog>'
        sys.exit()
    elif opt in ("-a", "--channel"):
        channel = int(arg)
    elif opt in ("-b", "--binning"):
        binning = eval(arg)
    elif opt in ("-p", "--plotvar"):
        plotvar = arg

def hist_maker(name, title, tr, br, cut, bin):
	if len(bin) == 3:
		hist = ROOT.TH1F(name, title, bin[0], bin[1], bin[2])
	else:
		hist = ROOT.TH1F(name, title, len(bin)-1, array.array('f',bin))
	tr.Project(name, br, cut)
	return copy.deepcopy(hist)

def hist_maker2D(name, title, tr, br, cut, bin):
	if len(bin) == 3:
		hist = ROOT.TH2F(name, title, bin[0], bin[1], bin[2], bin[0], bin[1], bin[2])
	else:
		hist = ROOT.TH2F(name, title, len(bin)-1, array.array('f',bin), len(bin)-1, array.array('f',bin))
	tr.Project(name, br, cut)
	return copy.deepcopy(hist)



signalfile = 'TT_powheg'
bkgfilelist = ['WJets', 'SingleTbar_tW', 'SingleTop_tW', 'ZZ', 'WW', 'WZ', 'DYJets', 'DYJets_10to50']
rdfilelist = ['MuonEG_Run2015','DoubleEG_Run2015','DoubleMuon_Run2015']

rootfileDir = "/cms/scratch/tt8888tt/cattools_v746/src/CATTools/CatAnalyzer/test/v7-4-6/"
channel_name = ['Combined', 'MuEl', 'ElEl', 'MuMu']

#weight = '1.'
weight = '(weight/lepweight)'
stepchcut = "step==6&&channel==%s"%channel

phase_cut = '(gentop1_pt!=-9&&gentop2_pt!=-9)'
if channel == 1: ttother_cut = "(parton_channel==2 && ((parton_mode1==1 && parton_mode2==2) || (parton_mode1==2 && parton_mode2==1)))"
elif channel == 2: ttother_cut = "(parton_channel==2 && (parton_mode1==2 && parton_mode2==2))"
elif channel == 3: ttother_cut = "(parton_channel==2 && (parton_mode1==1 && parton_mode2==1))"

#naming
vars = plotvar.split(',')
name = vars[0].split('_')[0][0]+vars[0].split('_')[1]
if 'tt' in plotvar: name = 't'+name
if 'jet' in plotvar: name = 'b'+name
name = name.replace('rapi','y')

vars = plotvar.split(',')

#Signal
rootfilename = rootfileDir+signalfile+".root"
tt = ROOT.TFile(rootfilename)
tree = tt.cattree.Get("nom")

out_rt = ROOT.TFile("%s_%s_mc_%s.root"%(signalfile.lower(),channel_name[channel],name), "RECREATE")
cnv = ROOT.TCanvas()

h1_Gen = ROOT.TH1F('g_'+name, 'Gen', len(binning)-1, array.array('f',binning))
h1_Parton = ROOT.TH1F('p_'+name, 'Parton', len(binning)-1, array.array('f',binning))
h1_Reco = ROOT.TH1F('c_'+name, 'Reco', len(binning)-1, array.array('f',binning))
h1_Reco_noPhaseCut = ROOT.TH1F('s6_'+name, 'Reco_noPhaseCut', len(binning)-1, array.array('f',binning))
h2_GenReco = ROOT.TH2F('h2_'+name, 'GenReco', len(binning)-1, array.array('f',binning), len(binning)-1, array.array('f',binning))

for var in vars:
	h1_Gen.Add(hist_maker('hist', plotvar, tree, 'gen'+var, "(channel==%s&&%s&&%s)*%s"%(channel,ttother_cut,phase_cut,weight), binning))
	h1_Parton.Add(hist_maker('hist', plotvar, tree, 'parton'+var, "(%s&&%s&&%s)*%s"%(stepchcut,ttother_cut,phase_cut,weight), binning))
	h1_Reco.Add(hist_maker('hist', plotvar, tree, var, "(%s&&%s&&%s)*%s"%(stepchcut,ttother_cut,phase_cut,weight), binning))
	h1_Reco_noPhaseCut.Add(hist_maker('hist', plotvar, tree, var, "(%s)*%s"%(stepchcut,weight), binning))
	h2_GenReco.Add(hist_maker2D('hist', plotvar, tree, '%s:gen%s'%(var,var), "(%s&&%s&&%s&&pseudoTop_channel==%s)*%s"%(stepchcut,ttother_cut,phase_cut,channel,weight), binning))

out_rt.Write()
out_rt.Close()

#Background
for mcname in bkgfilelist:
	print mcname
	rootfilename = rootfileDir+mcname+".root"
	tt = ROOT.TFile(rootfilename)
	tree = tt.cattree.Get("nom")

	out_rt = ROOT.TFile("%s_%s_mc_%s.root"%(mcname.lower(),channel_name[channel],name), "RECREATE")
	cnv = ROOT.TCanvas()

	h1_Reco_noPhaseCut = ROOT.TH1F('s6_'+name, 'Reco_noPhaseCut', len(binning)-1, array.array('f',binning))

	for var in vars:
		h1_Reco_noPhaseCut.Add(hist_maker('hist', plotvar, tree, var, "(%s)*%s"%(stepchcut,weight), binning))

	out_rt.Write()
	out_rt.Close()

#Data
rootfilename = rootfileDir+rdfilelist[channel]+".root"
tt = ROOT.TFile(rootfilename)
tree = tt.cattree.Get("nom")

out_rt = ROOT.TFile("%s_%s_data_%s.root"%(rdfilelist[channel].lower(),channel_name[channel],name), "RECREATE")
cnv = ROOT.TCanvas()

h1_Reco_noPhaseCut = ROOT.TH1F('s6_'+name, 'Reco_noPhaseCut', len(binning)-1, array.array('f',binning))
for var in vars:
	h1_Reco_noPhaseCut.Add(hist_maker('hist', plotvar, tree, var, "(%s)*%s"%(stepchcut,weight), binning))

out_rt.Write()
out_rt.Close()
