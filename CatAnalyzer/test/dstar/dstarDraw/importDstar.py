#!/usr/bin/env python

inputfile="../cattools/cattree_TT_powheg.root"



from ROOT import * #, CATTools.CatAnalyzer.CMS_lumi, json, os,copy
import numpy as np



file = TFile.Open(inputfile)
tree = file.Get("cattree/nom")
totEvent = tree.GetEntries()
out = TFile("dstar.root","RECREATE")

ntuple0  = TNtuple("nt0","nt0",  "pt:eta:phi:mass:dR:relPt:vProb:LXY:L3D:DCA:isFromB:isFromTop:trackQuality:diffMass")
ntuple1  = TNtuple("nt1","nt1",  "pt:eta:mass:vProb:LXY:L3D:DCA:trackQuality:lep1DR:lep2DR:lepMinDRMass:lepLowMass:correctM")
ntuple2  = TNtuple("nt2","nt2",  "pt:eta:mass:vProb:LXY:L3D:DCA:trackQuality:lep1DR:lep2DR:lepMaxDRMass:lepHighMass:WrongM")
ntuple3  = TNtuple("nt3","nt3",  "lep1_pid:lep2_pid:dstar_q:dR:relPt:LXY:L3D:isFromTop")

for idx, t in enumerate(tree) :
  print idx,totEvent
  if ( t.step<5 ) : continue
  #if idx >1000 : break
  for dstar_idx, dstar in enumerate(t.dstar) :
    lepMinDRMass = -9
    lepMaxDRMass = -9
    lepLowMass  = -9
    lepHighMass = -9
    correctM = -9
    wrongM = -9

    lep1DR = t.lep1.DeltaR(dstar)
    lep2DR = t.lep1.DeltaR(dstar)

    lep1_M = (t.lep1+dstar).M()
    lep2_M = (t.lep2+dstar).M()
    print idx,"dstar pt,eta,phi,M",dstar.Pt(),dstar.Eta(),dstar.Phi,dstar.M()
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

    if ( t.dstar_isFromTop[dstar_idx] ==6 ) :
      if ( t.lep1_pid <0 ) :
        correctM = lep1_M
        wrongM = lep2_M
      else :
        correctM = lep2_M
        wrongM = lep1_M
    elif (t.dstar_isFromTop[dstar_idx] == -6 ) :
      if ( t.lep1_pid>0 ) :
        correctM = lep1_M
        wrongM = lep2_M
      else :
        correctM = lep2_M 
        wrongM = lep1_M

    ntuple0.Fill( dstar.Pt(), dstar.Eta(), dstar.Phi(), dstar.M(), t.dstar_dRTrue[dstar_idx], t.dstar_relPtTrue[dstar_idx], t.dstar_vProb[dstar_idx], t.dstar_LXY[dstar_idx], t.dstar_L3D[dstar_idx], t.dstar_dca[dstar_idx],t.dstar_isFromB[dstar_idx], t.dstar_isFromTop[dstar_idx], t.dstar_trackQuality[dstar_idx], t.dstar_diffMass[dstar_idx] )
    ntuple1.Fill( dstar.Pt(), dstar.Eta(), dstar.M(), t.dstar_vProb[dstar_idx], t.dstar_LXY[dstar_idx], t.dstar_L3D[dstar_idx], t.dstar_dca[dstar_idx],t.dstar_trackQuality[dstar_idx], lep1DR, lep2DR, lepMinDRMass, lepLowMass, correctM) 
    ntuple2.Fill( dstar.Pt(), dstar.Eta(), dstar.M(), t.dstar_vProb[dstar_idx], t.dstar_LXY[dstar_idx], t.dstar_L3D[dstar_idx], t.dstar_dca[dstar_idx],t.dstar_trackQuality[dstar_idx], lep1DR, lep2DR, lepMaxDRMass, lepHighMass, wrongM)
    ntuple3.Fill(t.lep1_pid, t.lep2_pid, t.dstar_dau3_q[dstar_idx], t.dstar_dRTrue[dstar_idx], t.dstar_relPtTrue[dstar_idx], t.dstar_LXY[dstar_idx], t.dstar_L3D[dstar_idx], t.dstar_isFromTop[dstar_idx])
ntuple0.Write()
ntuple1.Write()
ntuple2.Write()
ntuple3.Write()
out.Close()
file.Close()
