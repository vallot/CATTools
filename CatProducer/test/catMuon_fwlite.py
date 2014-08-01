#! /usr/bin/env python

import ROOT
import sys
from DataFormats.FWLite import Events, Handle
from math import *

events = Events (['../prod/CAT.root'])

handleMu1  = Handle ("std::vector<cat::Muon>")
handleMu2  = Handle ("std::vector<cat::Muon>")

# for now, label is just a tuple of strings that is initialized just
# like and edm::InputTag
labelMu1 = ("catMuons","","PAT")
labelMu2 = ("catMuonsWeighted","","PAT")

ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.gROOT.SetStyle('Plain') # white background
histdefault = ROOT.TH1F ("histdefault", "relIso03", 30, 0, 0.3)
histweighted = ROOT.TH1F ("histweighted", "relIso03", 30, 0, 0.3)

# loop over events
count= 0
for event in events:
    count+=1 
    if count % 1000 == 0 :
	print count
    event.getByLabel (labelMu1, handleMu1)
    event.getByLabel (labelMu2, handleMu2)
    # get the product
    muons1 = handleMu1.product()
    muons2 = handleMu2.product()
    
    for m1,m2 in zip(muons1,muons2)  :
        if( m1.isTightMuon() ) :
          histdefault.Fill(m1.relIso()); 
        if( m2.isTightMuon() ) :
          histweighted.Fill(m2.relIso()); 
    
from ROOT import TCanvas, TLegend
c = TCanvas('c','Example',200,10,700,500)

histweighted.Draw()
histweighted.SetStats(0)
histweighted.GetXaxis().SetTitle("Relative Isolation")
histweighted.GetYaxis().SetTitle("Number of Muons")
histdefault.Draw("same")
histweighted.SetLineColor(2)

l = TLegend(0.70,0.70,0.80,0.80)
l.AddEntry(histdefault,"dbeta","L")
l.AddEntry(histweighted,"weighted","L")
l.SetTextSize(0.04)
l.SetFillColor(0)
l.SetLineColor(0)
l.Draw()

c.Update()
c.Draw()
c.Print("relIso.pdf")
