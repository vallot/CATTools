#!/usr/bin/env python
from ROOT import *
import json, os, sys, math, getopt 
from CATTools.CatAnalyzer.histoHelper import *
gROOT.SetBatch(True)

file0 = TFile.Open("jet_charge_tree_all.root")
tree = file0.Get("cattree/bjet_color")

c1 = makeCanvas("Lepton pdgId vs Jet charge",False)

h1 = TH2I("lep_jetq","Lepton pdgId VS Jet charge",30,-15,15,10,-5,5)
tree.Project("lep_jetq","jet_charge:lep_pdgId")

h1.SetXaxis().SetTitle("Lepton pdgId")
h1.SetYaxis().SetTitle("Jet Charge")

h1.Draw("colz")
c1.SaveAs("lep_jetq.png")


