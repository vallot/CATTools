#!/usr/bin/env python

inputfile="../cattools/cattree_TT_powheg.root"



from ROOT import * #, CATTools.CatAnalyzer.CMS_lumi, json, os,copy
import numpy as np



file = TFile.Open(inputfile)
tree = file.Get("cattree/nom")
totEvent = tree.GetEntries()
out = TFile("d0.root","RECREATE")

ntuple  = TNtuple("nt","nt",  "pt:eta:mass:vProb:LXY:L3D:DCA:trackQuality/I:lep1DR:lep2DR:lepMinDRMass/I:lepLowMass/I:correctM/I")
ntuple2 = TNtuple("nt2","nt2","step/I:pt:eta:phi:mass:dR:relPt:vProb:LXY:L3D:DCA:isFromB:isFromTop:trackQuality")

for idx, t in enumerate(tree) :
  print idx,totEvent
  if ( t.step<5 ) : continue
  for d0_idx, d0 in enumerate(t.d0) :
    lepMinDRMass1 = -9
    lepMinDRMass2 = -9
    lepLowMass1 = -9
    lepLowMass2 = -9
    correctM = -9
    wrongM = -9

    lep1DR = t.lep1.DeltaR(d0)
    lep2DR = t.lep1.DeltaR(d0)

    lep1_M = (t.lep1+d0).M()
    lep2_M = (t.lep2+d0).M()
    if ( lep1DR < lep2DR) :
      lepMinDRMass = 1 
    else :
      lepMinDRMass = 2

    if ( lep1_M< lep2_M) :
      lepLowMass = 1
    else :
      lepLowMass = 2

    if ( t.d0_isFromTop[d0_idx] ==6 ) :
      if ( t.lep1_pid <0 ) :
        correctM = 1
      else :
        correctM = 2
    elif (t.d0_isFromTop[d0_idx] == -6 ) :
      if ( t.lep1_pid>0 ) :
        correctM = 1
      else :
        wrongM = 2 
       
    ntuple.Fill( d0.Pt(), d0.Eta(), d0.M(), t.d0_vProb[d0_idx], t.d0_LXY[d0_idx], t.d0_L3D[d0_idx], t.d0_dca[d0_idx],t.d0_trackQuality[d0_idx], lep1DR, lep2DR, lepMinDRMass, lepLowMass, correctM) 
    ntuple2.Fill(t.step, d0.Pt(), d0.Eta(), d0.Phi(), d0.M(), t.d0_dRTrue[d0_idx], t.d0_relPtTrue[d0_idx], t.d0_vProb[d0_idx], t.d0_LXY[d0_idx], t.d0_L3D[d0_idx], t.d0_dca[d0_idx],t.d0_isFromB[d0_idx], t.d0_isFromTop[d0_idx], t.d0_trackQuality[d0_idx])
ntuple.Write()
out.Close()
file.Close()
