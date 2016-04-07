#!/usr/bin/env python

inputfile="../cattools/cattree_TT_powheg.root"



from ROOT import * #, CATTools.CatAnalyzer.CMS_lumi, json, os,copy

file = TFile.Open(inputfile)
tree = file.Get("cattree/nom")
totEvent = tree.GetEntries()
out = TFile("d0.root","RECREATE")

ntuple = TNtuple("nt","nt","step:pt:eta:phi:mass:dR:relPt:vProb:LXY:L3D:DCA")

for idx, t in enumerate(tree) :
  print idx,totEvent
  for d0_idx, d0 in enumerate(t.d0) :
    ntuple.Fill(t.step, d0.Pt(), d0.Eta(), d0.Phi(), d0.M(), t.d0_dRTrue[d0_idx], t.d0_relPtTrue[d0_idx], t.d0_vProb[d0_idx], t.d0_LXY[d0_idx], t.d0_L3D[d0_idx], t.d0_dca[d0_idx])
ntuple.Write()
out.Close()
file.Close()
