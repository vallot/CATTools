#!/usr/bin/env python
import math, array, ROOT, copy
import PhysicsTools.PythonAnalysis.rootplot.core as rootplotcore
import os
from ROOT import *
from CATTools.CatAnalyzer.histoHelper import *

def comparePlot(title, par, sig_hist, bkg_hist, x_title, y_title) :
  c1 = makeCanvas(title,False)
  sig_hist.SetLineColor(kBlue)
  sig_hist.SetFillStyle(3004)
  sig_hist.SetFillColor(kBlue)
  sig_hist.Scale( 1./sig_hist.GetEntries())
  sig_hist = setDefTH1Style(sig_hist, x_title, y_title)
  sig_hist.SetTitle("True "+par)

  bkg_hist.SetLineColor(kRed)
  bkg_hist.SetFillStyle(3005)
  bkg_hist.SetFillColor(kRed)
  bkg_hist.Scale( 1./bkg_hist.GetEntries())
  bkg_hist = setDefTH1Style(bkg_hist, x_title, y_title)
  bkg_hist.SetTitle("Fake "+par)

  if ( sig_hist.GetMaximum() > bkg_hist.GetMaximum()) : 
    bkg_hist.SetMaximum( sig_hist.GetMaximum() )
  bkg_hist.SetMaximum( bkg_hist.GetMaximum()*1.2)

  bkg_hist.Draw("hist")
  sig_hist.Draw("samehist")
  c1.BuildLegend()
  c1.SaveAs(title+".png")
  

infile = TFile.Open("../d0.root")
tree = infile.Get("nt0")

parameter = [ "DCA","LXY","L3D", "pt","abs(eta)","@pt.size()" ]
title = ["Distance of Closest approach","Lengh between PrimaryVertex and Secondary Vertex on XY Plane","Lengh between PrimaryVertex and Secondary Vertex on 3D space[cm]","Transvers Momentum","PeudoRapidy","Secondary Vertex Multiplicity"]
xaxis_title = ["DCA[cm]","LXY[cm]","L3D[cm]","pT[GeV/c^{2}]","|\eta|","nSV"]
yaxis_title = ["Normalized Entries"]*6
xmax = [0.2,0.5,1, 200.0, 5.0,10.0]
binning = [ [],[40,0.0,0.5], [40, 0.0, 1.0], [] ,[] ,[]]

signal_cut = "dR>0.0 && dR<0.1 && abs(relPt)>0.0&& abs(relPt)<0.1 && abs(isFromTop)==6"
bkg_cut = "!(%s)"%(signal_cut)


for idx, par in enumerate(parameter) :
  common_cut = "%s<%f&&%s>0.0&&"%(par,xmax[idx],par)
  sigcut = common_cut+signal_cut
  bkgcut = common_cut+bkg_cut
  sig_hist = getTH1( title[idx], binning[idx],tree, par, sigcut)
  bkg_hist = getTH1( title[idx], binning[idx],tree, par, bkgcut)
  comparePlot( title[idx], par,  sig_hist, bkg_hist, xaxis_title[idx], yaxis_title[idx]) 

