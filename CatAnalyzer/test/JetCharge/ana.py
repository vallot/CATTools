#!/usr/bin/env python
from ROOT import *
import json, os, sys, math, getopt 
from CATTools.CatAnalyzer.histoHelper import *
gROOT.SetBatch(True)

gStyle.SetPalette(1)
def anaTree( tree, label="true" ) :
  c0 = makeCanvas("Jet pdgId",False)
  h0 = getTH1("Jet pdgId ; pdgCode of Jet origin",[50,-25,25],tree, "jet_pdgId","1")
  h0.Draw("hist")
  c0.SaveAs("jet_pdgId%s.png"%label)
  

  c1 = makeCanvas("Lepton pdgId vs Jet charge",False)

  h1 = getTH2("Lepton pdgId VS Jet charge",[[30,-15,15],[10,-5,5]],tree, "jet_charge:lep_pdgId","1")

  h1.GetXaxis().SetTitle("Lepton pdgId")
  h1.GetYaxis().SetTitle("Jet Charge")
  h1.Draw("colz")
  c1.SaveAs("lep_jetq_%s.png"%label)

  default_cut = "jet_pdgId[0]!=0 && jet_pdgId[1]!=0"
  h2 = []

  separation = []

  for idx, pdgid in enumerate([5,-5]) :
    c2 = makeCanvas("Jet Charge from pdgId:%d"%(pdgid),False)
    h2.append(getTH1("Jet charge from pdgId:%d"%(pdgid),[10,-5,5],tree,"jet_charge[0]","%s&&jet_pdgId[0]==%d"%(default_cut, pdgid)))
    h2[idx].Sumw2()
    fr = h2[idx].Fit('gaus',"S")
    separation.append(fr.GetParams()[1])
    h2[idx].Draw() 
    c2.SaveAs("jetq_pdgid_%d_%s.png"%(pdgid,label))

  c3 = makeCanvas("Jet Charge from b quark")
  h2[0].SetLineColor(kRed)
  h2[0].SetFillColor(kRed)
  h2[0].SetFillStyle(3004)
  h2[1].SetLineColor(kBlue)
  h2[1].SetFillColor(kBlue)
  h2[1].SetFillStyle(3005)
  h2[0].Draw("hist")
  h2[1].Draw("samehist")
  c3.BuildLegend( 0.4,0.2,0.8,0.4)
  c3.SaveAs("jetq_from_b_%s.png"%label)

  h3 = []
  h3.append( getTH1( "delta_jetq_b",[10,-5,5],tree,"(jet_charge[0]-jet_charge[1])","%s&&jet_pdgId[0]==5"%(default_cut)))
  h3.append( getTH1( "delta_jetq_bbar",[10,-5,5],tree,"(jet_charge[0]-jet_charge[1])","%s&&jet_pdgId[0]==-5"%(default_cut)))
  h3[0].Sumw2()
  h3[1].Sumw2()
  fr = h3[0].Fit('gaus',"S")
  separation.append(fr.GetParams()[1])
  fr = h3[1].Fit('gaus',"S")
  separation.append(fr.GetParams()[1])

  for idx, pdgid in enumerate( [5,-5]) :
    c35 = makeCanvas("Diff jet Charge Jet1_pdgId: %d"%(pdgid),False)
    h3[idx].Draw()
    c35.SaveAs("diffJetCharge_pdgId_%d_%s.png"%(pdgid,label))


  h3[0].SetLineColor(kRed)
  h3[0].SetFillColor(kRed)
  h3[0].SetFillStyle(3004)
  h3[1].SetLineColor(kBlue)
  h3[1].SetFillColor(kBlue)
  h3[1].SetFillStyle(3005)

  c4 = makeCanvas("#delta (jet1_q-jet2_q) at b quark")
  h3[0].Draw("hist")
  h3[1].Draw("samehist")

    

  c4.BuildLegend( 0.4,0.2,0.8,0.4)
  c4.SaveAs("delta_jetq_%s.png"%label)

  h4 = []
  #h1 = getTH2("Lepton pdgId VS Jet charge",[[30,-15,15],[10,-5,5]],tree, "jet_charge:lep_pdgId","1")
  #h4 = TH2I("q_vs_true","charge vs true",2,0,2,2,0,2)
  h4.append(getTH2("q_vs_true_only_bquark",[[2,0,2],[2,0,2]],tree,"(jet_charge[0]<jet_charge[1]):(jet_pdgId[0]==5)","abs(jet_pdgId[0])==5&&abs(jet_pdgId[1])==5 && %s"%(default_cut)))
  h4.append(getTH2("q_vs_true_with_no_btag",[[2,0,2],[2,0,2]],tree,"(jet_charge[0]<jet_charge[1]):(jet_pdgId[0]==5)","%s"%default_cut))
  h4.append(getTH2("q_vs_true_with_jet1_btag",[[2,0,2],[2,0,2]],tree,"(jet_charge[0]<jet_charge[1]):(jet_pdgId[0]==5)","%s && jet_btag[0] "%default_cut))
  h4.append(getTH2("q_vs_true_with_jet1,2_btag",[[2,0,2],[2,0,2]],tree,"(jet_charge[0]<jet_charge[1]):(jet_pdgId[0]==5)","%s && jet_btag[0]&&jet_btag[1]"%default_cut))

  for idx in range(len(h4)) :
    h4[idx].GetXaxis().SetBinLabel(1,"true bbar")
    h4[idx].GetXaxis().SetBinLabel(2,"true b")
    h4[idx].GetYaxis().SetBinLabel(2,"Neg Charge(b)")
    h4[idx].GetYaxis().SetBinLabel(1,"Pos Charge(bbar)")

  gStyle.SetPaintTextFormat("2.4g")
  for idx, contents in enumerate(["only_bquark","no_btag","jet1_btag","jet1,2_btag"]) :
    c5 = TCanvas("q_true_%s"%(contents),"charge vs true",1000,600)
    h4[idx].Scale(1/ h4[idx].GetEntries()*100)
    h4[idx].Draw("colztext")
    c5.SaveAs("q_vs_true_%s_%s.png"%(contents, label))


if __name__ == "__main__" : 

  file0 = TFile.Open("data/jetcharge_all.root")
  tree = file0.Get("cattree/true_bjet_charge")
  recoTree = file0.Get("cattree/reco_bjet_charge")

  anaTree( tree )
  anaTree( recoTree, "reco")

