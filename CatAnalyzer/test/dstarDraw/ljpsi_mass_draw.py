#!/usr/bin/env python
import glob
from ROOT import *
from array import array
filelist= glob.glob("invMass*.root")

 
mean_gauss=[]
error_gauss = []
for file in filelist :
  title = file.split(".root")[0]
  file0 = TFile(file)
  hist1 = file0.Get("ljpsi1_mass") 
  hist2 = file0.Get("ljpsi2_mass")
 
  clone_hist = hist1.Clone()
  clone_hist2 = hist2.Clone()

  #clone_hist.Add( clone_hist2)
  clone_hist.Scale( 252.89 * 19.78 * 1000 / nevent[title])
  #clone_hist.Rebin(5)
  # Gauss
  #clone_hist.Fit("gaus")
  tf1 = TF1("f1","gaus",30,80)
  for x in range(5) :
    clone_hist.Fit(tf1)
  c1 = TCanvas("c1","c1",600,600)
  clone_hist.Draw()
  fitresult = TVirtualFitter.GetFitter()
  sig = fitresult.GetParameter(1)
  error = gitreuslt.GetErrors(1)
  mean_gauss.append(sig)
  mean_error.append(error)
  c1.SaveAs(title+"_gauss.png")

mass_various = [ 166.5, 169.5, 171.5, 173.5, 175.5, 178.5]
x_array = array("f", mass_various)


c1 = TCanvas("c1","c1",600,600)
y_array = array("f",mean_gauss)
h1 = TGraph( len(x_array), x_array, y_array)
h1.Draw("AL*")
c1.SaveAs("top_mass_gauss.png")
  
