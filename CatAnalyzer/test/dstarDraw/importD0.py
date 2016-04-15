#!/usr/bin/env python

inputfile="../cattools/cattree_TT_powheg.root"



from ROOT import * #, CATTools.CatAnalyzer.CMS_lumi, json, os,copy
import numpy as np



file = TFile.Open(inputfile)
tree = file.Get("cattree/nom")
totEvent = tree.GetEntries()
out = TFile("d0.root","RECREATE")

ntuple0  = TNtuple("nt0","nt0",  "step:pt:eta:phi:mass:dR:relPt:vProb:LXY:L3D:DCA:isFromB:isFromTop:trackQuality")
ntuple1  = TNtuple("nt1","nt1",  "pt:eta:mass:vProb:LXY:L3D:DCA:trackQuality:lep1DR:lep2DR:lepMinDRMass:lepLowMass:correctM")
ntuple2  = TNtuple("nt2","nt2",  "pt:eta:mass:vProb:LXY:L3D:DCA:trackQuality:lep1DR:lep2DR:lepMaxDRMass:lepHighMass:WrongM")

for idx, t in enumerate(tree) :
  print idx,totEvent
  if ( t.step<5 ) : continue
  for d0_idx, d0 in enumerate(t.d0) :
    lepMinDRMass = -9
    lepMaxDRMass = -9
    lepLowMass  = -9
    lepHighMass = -9
    correctM = -9
    wrongM = -9

    lep1DR = t.lep1.DeltaR(d0)
    lep2DR = t.lep1.DeltaR(d0)

    lep1_M = (t.lep1+d0).M()
    lep2_M = (t.lep2+d0).M()
    if ( lep1DR < lep2DR) :
      lepMinDRMass = lep1_M 
      lepMaxDRMass = lep2_M 
    else :
      lepMinDRMass = lep2_M
      lepMaxDRMass = lep1_M 

    if ( lep1_M< lep2_M) :
      lepLowMass  = lep1_M
      lepHighMass = lep2_M
    else :
      lepLowMass  = lep2_M
      lepHighMass = lep1_M

    if ( t.d0_isFromTop[d0_idx] ==6 ) :
      if ( t.lep1_pid <0 ) :
        correctM = lep1_M
        wrongM = lep2_M
      else :
        correctM = lep2_M
        wrongM = lep1_M
    elif (t.d0_isFromTop[d0_idx] == -6 ) :
      if ( t.lep1_pid>0 ) :
        correctM = lep1_M
        wrongM = lep2_M
      else :
        correctM = lep2_M 
        wrongM = lep1_M
       
    ntuple0.Fill(t.step, d0.Pt(), d0.Eta(), d0.Phi(), d0.M(), t.d0_dRTrue[d0_idx], t.d0_relPtTrue[d0_idx], t.d0_vProb[d0_idx], t.d0_LXY[d0_idx], t.d0_L3D[d0_idx], t.d0_dca[d0_idx],t.d0_isFromB[d0_idx], t.d0_isFromTop[d0_idx], t.d0_trackQuality[d0_idx])
    ntuple1.Fill( d0.Pt(), d0.Eta(), d0.M(), t.d0_vProb[d0_idx], t.d0_LXY[d0_idx], t.d0_L3D[d0_idx], t.d0_dca[d0_idx],t.d0_trackQuality[d0_idx], lep1DR, lep2DR, lepMinDRMass, lepLowMass, correctM) 
    ntuple2.Fill( d0.Pt(), d0.Eta(), d0.M(), t.d0_vProb[d0_idx], t.d0_LXY[d0_idx], t.d0_L3D[d0_idx], t.d0_dca[d0_idx],t.d0_trackQuality[d0_idx], lep1DR, lep2DR, lepMaxDRMass, lepHighMass, wrongM) 
ntuple0.Write()
ntuple1.Write()
ntuple2.Write()
out.Close()
file.Close()
