#!/usr/bin/env python
from ROOT import *
import json, os, sys, math, getopt 
from CATTools.CatAnalyzer.histoHelper import *
gROOT.SetBatch(True)

file0 = TFile.Open("jet_color_tree.root")
tree = file0.Get("cattree/bjet_color")
c1 = makeCanvas("c1",False)

tree.Draw("lep_pt[0]>>temp_h1")
h1 = file0.Get("temp_h1")
print h1.GetXaxis().GetNbins()
c1.SaveAs("h1.eps")
