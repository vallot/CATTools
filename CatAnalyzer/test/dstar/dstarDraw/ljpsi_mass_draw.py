#!/usr/bin/env python
import glob
from ROOT import *
import array as ar
import sys
from CATTools.CatAnalyzer.histoHelper import *
import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os, getopt, sys
import ROOT.TRandom3 as ran3

random = ran3.TRandom3()
random.SetSeed(1)

"""
def poissonFluctuation(hist, NumPE=10000, resetError=False ) :
  hist_for_PE_mean = TH1F("dist_mean_pe","Distribution of PE's mean value",120,0,120)
  hist_for_PE_err = TH1F("dist_err_pe","Distribution of PE's mean err value",100,0,50)
  sample_hist  = None
  for n in range( NumPE ) :
    clone2 = hist.Clone()
    for x in range(clone2.GetNbinsX()) :
      fluctuation = random.Poisson( clone2.GetBinContent(x+1))
      clone2.SetBinContent(x+1, fluctuation   )
      if ( resetError) : clone2.SetBinError(x+1, TMath.Sqrt(clone2.GetBinContent(x+1)))
    clone2.Fit( tf1 )
    if ( n==0 ) : 
      sample_hist = clone2
      c1 = makeCanvas(clone2.GetTitle(),False)
      c1 = setMargins(c1, False)
      setDefAxis( clone2.GetXaxis(), "M_{l+SV} Mean [GeV/c^2]",1)
      setDefAxis( clone2.GetYaxis(), "Entries",1)
      clone2.Draw()
      c1.Update()
      st = clone2.FindObject("stats")
      st.SetY1NDC(0.8)
      st.SetY2NDC(0.9)
      CMS_lumi.CMS_lumi( c1,0,11)
      c1.SaveAs(clone2.GetName()+".png")
    fitresult = TVirtualFitter.GetFitter()
    hist_for_PE_mean.Fill( fitresult.GetParameter(1))
    hist_for_PE_err.Fill( fitresult.GetParError(1))

    mean_list.append( fitresult.GetParameter(1))
    mean_err_list.append( fitresult.GetParError(1))

    sigma_list.append( fitresult.GetParameter(2))
    sigma_err_list.append( fitresult.GetParError(2))
    
  sig = sum(sig_list)/float(len(sig_list))
  error = sum(err_list)/float(len(err_list))
  return sample_hist, mean_hist, mean_err_hist, mean, mean_err, sigma, sigma_err 
"""    




if ( len(sys.argv) != 2 ) : 
  print "Please input inv file."
  sys.exit(-1)
type=sys.argv[1].split('invMass_')[-1]
infile = sys.argv[1]

massList = [  
"invMass_1665",
"invMass_1695",
"invMass_1715",
"invMass_1735",
"invMass_1755",
"invMass_1785",
"invMass_nominal"
]

data =["Run2015"]


mean_gauss=[]
error_gauss = []
x_error_gauss=[]

datalumi = 2.17
CMS_lumi.lumi_sqrtS = "%.1f fb^{-1}, #sqrt{s} = 13 TeV"%(datalumi)
datalumi = datalumi*1000 # due to fb
CMS_lumi.extraText   = "Private work"

CMS_lumi.relPosX = 0.06
CMS_lumi.cmsTextSize = 0.6

file0 = TFile(infile)




for mass in massList :
  hist1= file0.Get(mass)
  #hist1.SetStats(False)
  clone_hist = hist1.Clone()
  tf1 = TF1("f1","gaus",50,80)
  if ( mass == "invMass_nominal") :
    sig_list =[]
    err_list =[]
    """ 
    for n in range( 10000 ) :
      clone2 = clone_hist.Clone()
      for x in range(clone2.GetNbinsX()) :   
        fluctuation = random.Poisson( clone2.GetBinContent(x+1))
        clone2.SetBinContent(x+1, fluctuation   )
        clone2.SetBinError(x+1, TMath.Sqrt(clone2.GetBinContent(x+1)))
      clone2.Fit( tf1 )
      if ( n == 0 ) :
        c1 = makeCanvas("Pseudo-Data from ttbar MC",False)
        c1 = setMargins(c1, False)
        setDefAxis( clone2.GetXaxis(), "M_{l+SV} Mean [GeV/c^2]",1)
        setDefAxis( clone2.GetYaxis(), "Entries",1)
        clone2.Draw()
        c1.Update()
        st = clone2.FindObject("stats")
        st.SetY1NDC(0.8)
        st.SetY2NDC(0.9)
        CMS_lumi.CMS_lumi( c1,0,11)
        c1.SaveAs("PseudoData.png")
        
      fitresult = TVirtualFitter.GetFitter()
      hist_for_PE_mean.Fill( fitresult.GetParameter(1))
      hist_for_PE_err.Fill( fitresult.GetParError(1))
      sig_list.append( fitresult.GetParameter(1))
      err_list.append( fitresult.GetParError(1))
    sig = sum(sig_list)/float(len(sig_list))
    error = sum(err_list)/float(len(err_list))
    print "Sig : %f / %f , Err : %f / %f"%(sig, sum(sig_list)/float(len(sig_list)),error, sum(err_list)/float(len(err_list)))

    c1 = makeCanvas(hist_for_PE_mean.GetTitle(),False)
    gStyle.SetOptStat(1111)
    c1 = setMargins(c1, False)
    setDefAxis( hist_for_PE_mean.GetXaxis(), "M_{l+SV} Mean on Pseudo Data",1)
    setDefAxis( hist_for_PE_mean.GetYaxis(), "Entries",1)
    hist_for_PE_mean.SetStats(True)
    hist_for_PE_mean.Draw()
    c1.Update()
    st = hist_for_PE_mean.FindObject("stats")
    st.SetY1NDC(0.8)
    st.SetY2NDC(0.9)
    CMS_lumi.CMS_lumi( c1,0,11)
    c1.SaveAs(mass+"_PE_mean.png")

    c1 = makeCanvas(hist_for_PE_err.GetTitle(),False)
    gStyle.SetOptStat(1111)
    c1 = setMargins(c1, False)
    setDefAxis( hist_for_PE_err.GetXaxis(), "M_{l+svx} Error on Pseudo Data",1)
    setDefAxis( hist_for_PE_err.GetYaxis(), "Entries",1)
    gStyle.SetOptStat(11)
    hist_for_PE_err.SetStats(True)
    hist_for_PE_err.Draw()
    c1.Update()
    st = hist_for_PE_err.FindObject("stats")
    st.SetY1NDC(0.8)
    st.SetY2NDC(0.9)
    CMS_lumi.CMS_lumi( c1,0,11)
    c1.SaveAs(mass+"_PE_err.png")
    gStyle.SetOptStat(0)
    """
         
  clone_hist.Fit(tf1)
  c1 = makeCanvas(mass,False)
  c1 = setMargins(c1, False)
  c1.SetTopMargin(30)
  setDefAxis( clone_hist.GetXaxis(), "M_{l+SV} [GeV/c^2]",1)
  setDefAxis( clone_hist.GetYaxis(), "Entries",1)
  clone_hist.SetTitle("Invariant Mass of Lepton + Secondary Vertex")
  gStyle.SetTitleFontSize(0.2);
  clone_hist.Draw()
  leg = TLegend(0.4,0.15, 0.7,0.25)
  if ( mass != "invMass_nominal") : leg.AddEntry( clone_hist, "M_{top}=%.1fGeV/c^2"%(float(mass.split("invMass_")[-1])/10.)  )
  else : leg.AddEntry(clone_hist, "M_{top} nominal")
  leg.SetTextSize(0.03)
  leg.Draw()
  c1.Update()
  st = clone_hist.FindObject("stats")
  st.SetY1NDC(0.8)
  st.SetY2NDC(0.9)
  CMS_lumi.CMS_lumi( c1, 0, 11 )
  c1.Modified()
  c1.Update()

  if ( mass != "invMass_nominal") :
    fitresult = TVirtualFitter.GetFitter()
    sig = fitresult.GetParameter(1)
    error = fitresult.GetParError(1)

  mean_gauss.append(sig)
  error_gauss.append(error)

  x_error_gauss.append(0)
  c1.SaveAs(mass+"_gauss.png")

mass_various = [ 166.5, 169.5, 171.5, 173.5, 175.5, 178.5]
x_array = ar.array("f", mass_various)
x_error_array = ar.array("f",x_error_gauss[:-1])
error_array = ar.array("f",error_gauss[:-1])





## Below lines for Real Data analysis.
# Blinding!!

rdhist= file0.Get("Run2015")

if ( rdhist.Integral() / rdhist.GetNbinsX() *2  < 10 ) :
  #rdhist.Rebin(2) 
  pass

tf_rd = TF1("f2","gaus",50,80)
rdhist.Fit(tf_rd)
c1 = makeCanvas("Run2015",False)
c1 = setMargins(c1, False)
setDefAxis( rdhist.GetXaxis(), "M_{l+SV} [GeV/c^2]",1)
setDefAxis( rdhist.GetYaxis(), "Entries",1)
rdhist.SetTitle("Invariant Mass of Lepton + Seconary Vertex on Run2015")

#rdhist.SetStats(False)
gStyle.SetTitleFontSize(0.2)

rdhist.Draw()
c1.Update()
st = rdhist.FindObject("stats")
st.SetY1NDC(0.8)
st.SetY2NDC(0.9)
CMS_lumi.CMS_lumi( c1,0,11)
c1.Modified()
c1.Update()

c1.SaveAs("run2015_"+type+".png")

fitresult_rd = TVirtualFitter.GetFitter()
mean = fitresult_rd.GetParameter(1)
error = fitresult_rd.GetParError(1)

c1 = makeCanvas("M_{t} vs M_{l+secvtx}",False)
y_array = ar.array("f",mean_gauss)
h1 = TGraphErrors( len(x_array), x_array, y_array, x_error_array, error_array)
h1.SetName("M_{t} MC samples")
h1.SetMarkerStyle( 23)
h1.SetMarkerColor ( kRed )
h1.SetMarkerSize( 1.5 )
#h1.Draw("ALP")
h1.Fit("pol1")
fitresult2 = TVirtualFitter.GetFitter()
y_0  = fitresult2.GetParameter(0)
slope    = fitresult2.GetParameter(1)


top_mass = (mean-y_0)/ slope
top_mass_up =(mean+error-y_0)/ slope
top_mass_error = top_mass_up - top_mass
print "Top mass ",top_mass,"+-", top_mass_error

arr_lsvtxMass = ar.array("f",[mean])
arr_topMass = ar.array("f",[top_mass])
y_error_array = ar.array("f",[error])
x_error_array = ar.array("f",[top_mass_error])


h2 = TGraphErrors( 1, arr_topMass, arr_lsvtxMass, x_error_array, y_error_array)
h2.SetName("Run2015")
h2.SetMarkerColor( kBlue )
h2.SetMarkerStyle( 22)
h2.SetMarkerSize(2) 

x_array = ar.array("f", mass_various[:-1])
x_error_array = ar.array("f",x_error_gauss[:-1])
error_array = ar.array("f",error_gauss[:-1])


nominal_mass = (mean_gauss[-1]-y_0) / slope
nominal_mass_up     = (mean_gauss[-1]-y_0+ error_gauss[-1]) / slope
nominal_mass_error  = nominal_mass_up - nominal_mass


h3 =TGraphErrors(1, ar.array("f",[nominal_mass]), ar.array("f",[mean_gauss[-1]]) , ar.array("f",[x_error_array[-1]]), ar.array("f",[nominal_mass_error]))
h3.SetName("M_{t} MC nominal")
h3.SetMarkerColor( kBlue)
h3.SetMarkerStyle(25)
h3.SetMarkerSize(2)



#h2.Draw("LP")
mg = TMultiGraph();
mg.SetTitle("M_{t} vs Invariant mass of (Lepton + Secondary Vertex); M_{t} [GeV/c^{2}] ; Invariant M_{l+D}[GeV/c^{2}]")
mg.Add(h1)
mg.Add(h2) 
#mg.Add(h3)
mg.Draw("AP")
mg.GetXaxis().SetRangeUser(165,180)


mg.GetYaxis().SetRangeUser( min(mean_gauss[:-2])-5 ,max(mean_gauss[:-2])+5)

setDefAxis( mg.GetXaxis(), "M_{top} [GeV/c^2]",1)
setDefAxis( mg.GetYaxis(), "Invariant M_{l+D}[GeV/c^{2}]",1)


leg = TLegend(0.2,0.7, 0.5,0.8)
leg.AddEntry( h1, "M_{top} MC samples")
leg.AddEntry( h2, "Run2015")
#leg.AddEntry( h3, "MC nominal sample")
leg.Draw()

pavet = TPaveText()
pavet.AddText("Measured M_{top} : %3.2f #pm %.2f GeV/c^2 (Stat.)"%(top_mass,top_mass_error))
#pavet.AddText("Nominal M_{top} : %3.2f #pm %.2f GeV/c^2 (Stat.)"%(nominal_mass,nominal_mass_error))
pavet.SetTextSize(0.03)
pavet.SetX1NDC(0.4)
pavet.SetX2NDC(0.9)
pavet.SetY1NDC(0.15)
pavet.SetY2NDC(0.25)

pavet.Draw()
c1.Update()

st = h1.FindObject("stats")
st.SetY1NDC(0.3)
st.SetY2NDC(0.4)

CMS_lumi.lumiTextSize = 0.7
CMS_lumi.cmsTextSize = 0.9
CMS_lumi.CMS_lumi( c1,0,11)
c1.Modified()
c1.Update()
c1.SaveAs("top_mass_gauss_"+type+".png")
