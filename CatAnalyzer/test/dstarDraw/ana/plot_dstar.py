#!/usr/bin/env python
import math, array, ROOT, copy
import PhysicsTools.PythonAnalysis.rootplot.core as rootplotcore
import os
from ROOT import *
from CATTools.CatAnalyzer.histoHelper import *
import CATTools.CatAnalyzer.CMS_lumi
datalumi=2.17
CMS_lumi.lumi_sqrtS = "%.1f fb^{-1}, #sqrt{s} = 13 TeV"%(datalumi)
datalumi = datalumi*1000 # due to fb
CMS_lumi.extraText   = "Private work"


def comparePlot(title, par, sig_hist, bkg_hist, x_title, y_title) :
  c1 = makeCanvas(title,False)
  sig_hist.SetLineColor(kBlue)
  sig_hist.SetFillStyle(3004)
  sig_hist.SetFillColor(kBlue)
  sig_hist = setDefTH1Style(sig_hist, x_title, y_title)
  sig_hist.SetTitle("True "+par)

  bkg_hist.SetLineColor(kRed)
  bkg_hist.SetFillStyle(3005)
  bkg_hist.SetFillColor(kRed)
  bkg_hist = setDefTH1Style(bkg_hist, x_title, y_title)
  bkg_hist.SetTitle("Fake "+par)


  if ( sig_hist.GetEntries() !=0 ) : sig_hist.Scale( 1./sig_hist.GetEntries())
  if ( bkg_hist.GetEntries() != 0 ) : bkg_hist.Scale( 1./bkg_hist.GetEntries())


  if ( sig_hist.GetMaximum() > bkg_hist.GetMaximum()) : 
    bkg_hist.SetMaximum( sig_hist.GetMaximum() )
  bkg_hist.SetMaximum( bkg_hist.GetMaximum()*1.2)

  bkg_hist.Draw("hist")
  sig_hist.Draw("samehist")
  c1.BuildLegend()
  c1.cd()
  iPos = 11
  if( iPos==0 ):
      cmsLumi.relPosX = 0.12
  CMS_lumi.CMS_lumi(c1, 0, iPos)

  c1.Modified()
  c1.Update()
  c1.SaveAs(title+".png")


  

infile = TFile.Open("../dstar.root")
tree = infile.Get("nt0")

parameter = [ "DCA","LXY","L3D", "pt","abs(eta)","@pt.size()","diffMass","mass" ]
title = ["Distance of Closest approach","Lengh between PrimaryVertex and Secondary Vertex on XY Plane","Lengh between PrimaryVertex and Secondary Vertex on 3D space[cm]","Transvers Momentum","PeudoRapidy","Secondary Vertex Multiplicity", "MassDifference", "Mass"]
xaxis_title = ["DCA[cm]","LXY[cm]","L3D[cm]","pT[GeV/c^{2}]","|\eta|","nSV","M_{D*} - M_{D0}","Mass [GeV/c^2]" ]
yaxis_title = ["Entries"]*8
xmax = [0.2,0.5,1, 200.0, 5.0,10.0, 0.25, 2.2]
binning = [ [],[40,0.0,0.5], [40, 0.0, 1.0], [] ,[] ,[],[], []]

signal_cut = "dR>0.0 && dR<0.1 && abs(relPt)>0.0&& abs(relPt)<0.1 && abs(isFromTop)==6"
bkg_cut = "!(%s)"%(signal_cut)
#applied_cut = "LXY>0.1 && L3D>0.2&&"
applied_cut = ""

for idx, par in enumerate(parameter) :
  common_cut = "%s<%f&&%s>0.0&&"%(par,xmax[idx],par)
  common_cut = common_cut + applied_cut
  sigcut = common_cut+signal_cut
  bkgcut = common_cut+bkg_cut
  sig_hist = getTH1( title[idx], binning[idx],tree, par, sigcut)
  bkg_hist = getTH1( title[idx], binning[idx],tree, par, bkgcut)
  comparePlot( title[idx], par,  sig_hist, bkg_hist, xaxis_title[idx], yaxis_title[idx]) 

