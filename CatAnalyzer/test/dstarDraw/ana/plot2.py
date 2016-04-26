#!/usr/bin/env python
import os
from ROOT import *
from CATTools.CatAnalyzer.histoHelper import *
infile = TFile.Open("../dstar.root")

tree = infile.Get("nt1")
#tree2 = infile.Get("nt2")

#tree.Print()

#parameter = [ "DCA","LXY","L3D", "pt","abs(eta)", ]


parameter_good = ['lepMinDRMass', 'lepLowMass','correctM']
parameter_bad =  ['lepMaxDRMass', 'lepHighMass','WrongM']
color = [kBlue, kRed, kBlack]
xmax = [0.2,1.0,3.0, 200.0, 5.0]

#signal_cut = "dR>0.0 && dR<0.1 && abs(relPt)>0.0&& abs(relPt)<0.1"
#bkg_cut = "(dR>0.1 || abs(relPt)>0.1)"  

optimal_cut = "LXY>0.1 && L3D>0.2 && correctM>0"


c1 = TCanvas("l_svx","l_svx")
cut = optimal_cut
a = "Lepton+SecVtx with deltaR",
tree.Project("l_secvtx_dR","lepMinDRMass:correctM",cut,"colz")

c1.BuildLegend()
c1.SaveAs("lepMin.png")

c1 = TCanvas("l_svx","l_svx")
tree.Draw("lepLowMass:correctM",cut,'colz')
c1.BuildLegend()
c1.SaveAs("lepLow.png")

