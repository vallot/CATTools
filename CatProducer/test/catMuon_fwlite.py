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

#ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.gROOT.SetStyle('Plain') # white background
histdefault = ROOT.TH1F ("histdefault", "neutral", 100, 0, 1)
histweighted = ROOT.TH1F ("histweighted", "neutral", 100, 0, 1)

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
        histdefault.Fill(m1.neutralHadronIso()); 
        histweighted.Fill(m2.neutralHadronIso()); 
        if m1.neutralHadronIso() > m2.neutralHadronIso():
	  print "Muon neutral iso : default vs weighed = ", m1.neutralHadronIso(),m2.neutralHadronIso()
        if m1.photonIso() > m2.photonIso():
	  print "Muon photon iso : default vs weighed = ", m1.photonIso(),m2.photonIso()
    
from ROOT import TCanvas
c = TCanvas('c','Example',200,10,700,500)

histdefault.Draw()
#histweighted.Draw("same")
c.Update()
#c.Draw()
