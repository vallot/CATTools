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

signalfile = 'TT_powheg'
bkgfilelist = ['WJets', 'SingleTbar_tW', 'SingleTop_tW', 'ZZ', 'WW', 'WZ', 'DYJets', 'DYJets_10to50']
rdfilelist = ['MuonEG_Run2015','DoubleEG_Run2015','DoubleMuon_Run2015']
rootfileDir = "file:/xrootd/store/user/tt8888tt/v7-6-2/"
channel_name = ['MuEl', 'ElEl', 'MuMu']

#set weight
bkweight = '(genweight*puweight)'
weight = '(genweight*puweight*lepweight*csvweights[0])'

#cuts
stepch_tcut = 'step==6&&channel==%s'%channel
filter_tcut = 'tri==1&&filtered==1'
phase_tcut  = 'gentop1_pt!=-9&&gentop2_pt!=-9&&pseudoTop_channel==%s'%channel
if   channel == 1: ttother_tcut = "(parton_channel==2 && ((parton_mode1==1 && parton_mode2==2) || (parton_mode1==2 && parton_mode2==1)))"
elif channel == 2: ttother_tcut = "(parton_channel==2 && (parton_mode1==2 && parton_mode2==2))"
elif channel == 3: ttother_tcut = "(parton_channel==2 && (parton_mode1==1 && parton_mode2==1))"

#naming
vars = plotvar.split(',')
name = vars[0].split('_')[0][0]+vars[0].split('_')[1]
if 'tt' in plotvar: name = 't'+name
if 'jet' in plotvar: name = 'b'+name
name = name.replace('rapi','y')
title = plotvar

#adding Background to s6
bkg_s6 = []
for mcname in bkgfilelist:
	rootfilename = rootfileDir+mcname+".root"
	tt = ROOT.TFile(rootfilename)
	tree = tt.cattree.Get("nom")
	if len(binning) == 3:
		tmp_bk = ROOT.TH1F('tmp_s6_'+name, 'Reco_noPhaseCut', binning[0], binning[1], binning[2])
	else:
		tmp_bk = ROOT.TH1F('tmp_s6_'+name, 'Reco_noPhaseCut', len(binning)-1, array.array('f',binning))
	for var in vars:
		tmp_bk.Add(getTH1(title, binning, tree, var, "(%s&&%s)*%s"%(stepch_tcut,filter_tcut,weight)))
	bkg_s6.append(copy.deepcopy(tmp_bk))

#Signal
rootfilename = rootfileDir+signalfile+".root"
tt = ROOT.TFile(rootfilename)
tree = tt.cattree.Get("nom")

mc_rt = ROOT.TFile("%s_%s_mc_%s.root"%(signalfile.lower(),channel_name[channel-1],name), "RECREATE")
cnv = ROOT.TCanvas()

if len(binning) == 3: 
	h1_Gen = ROOT.TH1F('c_'+name, 'Gen', binning[0], binning[1], binning[2])
	h1_Parton = ROOT.TH1F('p_'+name, 'Parton', binning[0], binning[1], binning[2])
	h1_Reco = ROOT.TH1F('c6_'+name, 'Reco', binning[0], binning[1], binning[2])
	h1_Reco_noPhaseCut = ROOT.TH1F('s6_'+name, 'Reco_noPhaseCut', binning[0], binning[1], binning[2])
	h2_GenReco = ROOT.TH2F('h2_'+name, 'GenReco', binning[0], binning[1], binning[2], binning[0], binning[1], binning[2])
	h1_Signal_rate = ROOT.TH1F('h_'+name+'_rate', 'Signal_rate', binning[0], binning[1], binning[2])
else:
	h1_Gen = ROOT.TH1F('c_'+name, 'Gen', len(binning)-1, array.array('f',binning))
	h1_Parton = ROOT.TH1F('p_'+name, 'Parton', len(binning)-1, array.array('f',binning))
	h1_Reco = ROOT.TH1F('c6_'+name, 'Reco', len(binning)-1, array.array('f',binning))
	h1_Reco_noPhaseCut = ROOT.TH1F('s6_'+name, 'Reco_noPhaseCut', len(binning)-1, array.array('f',binning))
	h2_GenReco = ROOT.TH2F('h2_'+name, 'GenReco', len(binning)-1, array.array('f',binning), len(binning)-1, array.array('f',binning))
	h1_Signal_rate = ROOT.TH1F('h_'+name+'_rate', 'Signal_rate', len(binning)-1, array.array('f',binning))

for var in vars:
	h1_Gen.Add(getTH1(title, binning, tree, 'gen'+var, "(%s&&%s)*%s"%(ttother_tcut,phase_tcut,bkweight)))
	h1_Parton.Add(getTH1(title, binning, tree, 'parton'+var, "(%s)*%s"%(ttother_tcut,bkweight)))
	h1_Reco.Add(getTH1(title, binning, tree, var, "(%s&&%s&&%s&&%s)*%s"%(stepch_tcut,filter_tcut,ttother_tcut,phase_tcut,weight)))
	h1_Reco_noPhaseCut.Add(getTH1(title, binning, tree, var, "(%s&&%s)*%s"%(stepch_tcut,filter_tcut,weight)))
	h2_GenReco.Add(getTH2(title, binning, tree, 'gen%s:%s'%(var,var), "(%s&&%s&&%s&&%s)*%s"%(stepch_tcut,filter_tcut,ttother_tcut,phase_tcut,weight)))

for h in bkg_s6:
	h1_Reco_noPhaseCut.Add(h)

h1_Signal_rate.Divide(h1_Reco, h1_Reco_noPhaseCut)

mc_rt.Write()
mc_rt.Close()

#Data
rootfilename = rootfileDir+rdfilelist[channel-1]+".root"
tt = ROOT.TFile(rootfilename)
tree = tt.cattree.Get("nom")

data_rt = ROOT.TFile("data_%s_%s.root"%(channel_name[channel-1],name), "RECREATE")
cnv = ROOT.TCanvas()

if len(binning) == 3:
	h1_Reco_noPhaseCut = ROOT.TH1F('s6_'+name, 'Reco_noPhaseCut', binning[0], binning[1], binning[2])
else:
	h1_Reco_noPhaseCut = ROOT.TH1F('s6_'+name, 'Reco_noPhaseCut', len(binning)-1, array.array('f',binning))
for var in vars:
	h1_Reco_noPhaseCut.Add(getTH1(title, binning, tree, var, "(%s&&%s)*%s"%(stepch_tcut,filter_tcut,weight)))

data_rt.Write()
data_rt.Close()

