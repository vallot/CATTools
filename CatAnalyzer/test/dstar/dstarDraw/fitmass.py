#!/usr/bin/env python
import types, ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os, getopt, sys, copy
from CATTools.CatAnalyzer.histoHelper import *
import DYestimation
ROOT.gROOT.SetBatch(True)

ROOT.gSystem.Load("libRooFit")


nFitType = 0


def myGetNameFit(): 
  if nFitType == 0 : # Landau
    return "landau"
  if nFitType == 1 : # Gaussian
    return "gaus"


def myGetFitFunc(name, title, xvar, centralvar, sigmavar):
  if nFitType == 0 : # Landau
    return ROOT.RooLandau(name, title, xvar, centralvar, sigmavar)
  if nFitType == 1 : # Gaussian
    return ROOT.RooGaussian(name, title, xvar, centralvar, sigmavar)


def myFitInvMass(hist, strName, x, binning, range, dicStyle):
  tf1 = ROOT.TF1("f1_fit", myGetNameFit(), binning[ 1 ], binning[ 2 ])
  
  fMinRange = binning[ 1 ]
  fMaxRange = binning[ 2 ]
  
  if len(range) == 0 : 
    hist.Fit(tf1)
  else : 
    fMinRange = range[ 0 ]
    fMaxRange = range[ 1 ]
    
    hist.Fit(tf1, "", "", fMinRange, fMaxRange)
  
  #x = ROOT.RooRealVar("invmass",hist.GetXaxis().GetTitle(), binning[ 1 ], binning[ 2 ])
  xfitarg = ROOT.RooArgList(x, "invmass")
  dh = ROOT.RooDataHist("dh", "data histogram", xfitarg, hist)
  
  dSigma = tf1.GetParameter(2)
  dSigmaError = tf1.GetParError(2) * 10.0

  CB_central = ROOT.RooRealVar("mean" + dicStyle[ "suffix" ], "mean" + dicStyle[ "suffix" ], 
    tf1.GetParameter(1), fMinRange, fMaxRange)
  CB_sigma   = ROOT.RooRealVar("sigma" + dicStyle[ "suffix" ], "sigma" + dicStyle[ "suffix" ], 
    dSigma, dSigma - dSigmaError, dSigma + dSigmaError)
  
  #sig_pdf   = ROOT.RooLandau("sig_fit", "signal p.d.f", x, CB_central, CB_sigma)
  sig_pdf   = myGetFitFunc("sig_fit", "signal p.d.f", x, CB_central, CB_sigma)
  model = sig_pdf
  
  fitResult = model.fitTo(dh, ROOT.RooFit.Extended(False), ROOT.RooFit.Save())
  
  top_mass_frame = x.frame()
  strNameHisto = "roofithisto_" + strName
  dh.plotOn(top_mass_frame, ROOT.RooFit.Name(strNameHisto), 
    ROOT.RooFit.MarkerColor(dicStyle[ "color" ]), ROOT.RooFit.MarkerStyle(dicStyle[ "marker" ]))
  
  #model.paramOn(top_mass_frame, ROOT.RooFit.Format("NELU", ROOT.RooFit.AutoPrecision(2)), 
  #  ROOT.RooFit.Layout(0.1, 0.4, 0.9))

  strNameCurve = "fittingcurve_" + strName
  model.plotOn(top_mass_frame, ROOT.RooFit.Name(strNameCurve), 
    ROOT.RooFit.LineColor(dicStyle[ "color" ]), ROOT.RooFit.LineStyle(dicStyle[ "line" ]))
  model.plotOn(top_mass_frame, ROOT.RooFit.Components(ROOT.RooArgSet(sig_pdf)),
    ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(dicStyle[ "line" ]))
  
  return {"frame":top_mass_frame, "histo":strNameHisto, "graph":strNameCurve, 
    "peak_val":CB_central.getVal(), "peak_err":CB_central.getError(), 
    "sigma_val":CB_sigma.getVal(), "sigma_err":CB_sigma.getError(), "chi2":top_mass_frame.chiSquare()}


################################################################
## Getting inputs
################################################################
strFitType = "Landau"
rangeFit = []

try:
    opts, args = getopt.getopt(sys.argv[2:], "hdnot:r:",["type","range"])
except getopt.GetoptError:          
    print 'Usage : ./fitmass.py <json> -t <type> -r <range>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./fitmass.py <json> -t <type> -r <range>'
        sys.exit()
    elif opt in ("-t", "--type"):
        strFitType = arg
    elif opt in ("-r", "--range"):
        rangeFit = eval(arg)

if strFitType == "Landau" :
  nFitType = 0
elif strFitType == "Gaussian" :
  nFitType = 1
  print "Gaussian go"
else :
  print "Error: wrong fitting type"
  sys.exit(2)

################################################################
## Reading the configuration file
################################################################
try:
  if len(sys.argv) is 1 : 
    print "Usage : %s [rootfilename]"
    sys.exit(2)

  fileDicHists = open(sys.argv[ 1 ])
  dicHists = json.load(fileDicHists)
  fileDicHists.close()
except IOError:
  print "Usage : %s [rootfilename]  # rootfilename must be correct"
  sys.exit(2)

listPartPath = sys.argv[ 1 ].split("/")
strPathPos = ""

for i in range(len(listPartPath) - 1) :
  strPathPos = strPathPos + listPartPath[ i ] + "/"

################################################################
## Initializing the result root file
################################################################
rootHists = ROOT.TFile.Open(strPathPos + dicHists["rootfilename"])
rootFits = ROOT.TFile.Open("fres/fittings_" + dicHists["rootfilename"],"RECREATE")

binning = dicHists["binning"]
x_name = dicHists["x_name"]
y_name = dicHists["y_name"]

canvasMain = makeCanvas("canvasMain")
strKeyData = ""

dicHistoStyle = {
  "data":    {"label":"Data",   "marker":ROOT.kFullDotLarge, "color":ROOT.kBlack, "line":ROOT.kDashed}, 
  "1695":    {"label":"169.5",  "marker":ROOT.kMultiply,     "color":8,           "line":ROOT.kSolid}, 
  "nominal": {"label":"172.5", "marker":ROOT.kCircle,       "color":ROOT.kBlue,  "line":ROOT.kSolid}, 
  "1755":    {"label":"175.5",  "marker":ROOT.kPlus,         "color":ROOT.kPink,  "line":ROOT.kSolid}
}

################################################################
##  Getting peaks of TT samples
################################################################
for strKey in dicHists.keys():
  if type(dicHists[ strKey ]) != types.DictType or "type" not in dicHists[ strKey ] : continue
  if dicHists[ strKey ][ "type" ] not in ["TT_onlytt", "TT_withbkg", "data"] : continue
  
  if dicHists[ strKey ][ "type" ] == "data" : strKeyData = strKey
  
  # If str is missing, since strKey is in unicode type, it yields an error.
  histCurr = ROOT.TH1D(rootHists.Get(str(strKey)))
  
  strStyleCurr = ""
  
  for strStyle in dicHistoStyle.keys() : 
    if strStyle in strKey : 
      strStyleCurr = strStyle
      break
  
  strSuffix = ""
  if dicHists[ strKey ][ "type" ] in ["TT_onlytt", "TT_withbkg"] : 
    strSuffix = "_%0.1f"%(dicHists[ strKey ][ "mass" ])
  
  dicHistoStyle[ strStyleCurr ][ "suffix" ] = strSuffix
  dicHists[ strKey ][ "label" ] = dicHistoStyle[ strStyleCurr ][ "label" ]
  
  # -- Fitting by RooFit (no bkg)
  print ""
  print "#####################################################"
  print "### Fitting (%s) by RooFit is beginning"%(strKey)
  print "#####################################################"
  
  x = ROOT.RooRealVar("invmass_" + strKey, histCurr.GetXaxis().GetTitle(), binning[ 1 ], binning[ 2 ])
  
  nFitTypeOrg = nFitType
  
  if "Gen" not in dicHists or dicHists[ "Gen" ] == 0: 
    print "##### Using the given distribution #####"
  else:
    print "##### Using Gaussian distribution for Gen #####"
    nFitType = 1

  dicHists[ strKey ][ "fitres" ] = myFitInvMass(histCurr, strKey, x, binning, 
    rangeFit, dicHistoStyle[ strStyleCurr ])

  nFitType = nFitTypeOrg

  # -- Dumping
  print "------ Peak (%s) by RooFit : %f, Err : %f, chi : %f"%(strKey, 
    dicHists[ strKey ][ "fitres" ][ "peak_val" ], 
    dicHists[ strKey ][ "fitres" ][ "peak_err" ], 
    dicHists[ strKey ][ "fitres" ][ "chi2" ])
  print "------      (sigma (%s) : %f, %f)"%(strKey, dicHists[ strKey ][ "fitres" ][ "sigma_val" ], 
    dicHists[ strKey ][ "fitres" ][ "sigma_err" ])
  
  # -- Setting the name and shapes of histogram
  frameCurr = dicHists[ strKey ][ "fitres" ][ "frame" ]
  #dicHists[ strKey ][ "fitres" ][ "x" ] = x
  #frameCurr = x.frame()
  
  frameCurr.SetName(strKey)
  frameCurr.SetTitle(histCurr.GetTitle())
  
  setDefAxis(frameCurr.GetXaxis(), x_name, 1.0)
  setDefAxis(frameCurr.GetYaxis(), y_name, 1.2)
  frameCurr.SetMinimum(0)
  
  # -- Saving the histogram and the fitting curve
  rootFits.cd()
  frameCurr.Write()

################################################################
## Drawing all multiple plots
################################################################
for strType in [ "TT_onlytt", "TT_withbkg" ] : 
  # -- Making a legend
  leg = ROOT.TLegend(0.76,0.73,0.76+0.20,0.91)
  leg.SetBorderSize(0)
  leg.SetNColumns(2)
  leg.SetTextSize(0.040)
  leg.SetTextFont(42)
  leg.SetLineColor(0)
  leg.SetFillColor(0)
  leg.SetFillStyle(0)
  
  # -- Finding the maximum
  dMax = 0.0
  
  for strKey in dicHists.keys():
    if type(dicHists[ strKey ]) != types.DictType or "type" not in dicHists[ strKey ] : continue
    if dicHists[ strKey ][ "type" ] != strType : continue
    
    dMaxCurr = dicHists[ strKey ][ "fitres" ][ "frame" ].GetMaximum()
    if dMax < dMaxCurr : dMax = dMaxCurr
  
  # -- Drawing the data plot first
  dicHists[ strKeyData ][ "fitres" ][ "frame" ].SetMaximum(dMax)
  
  canvasMain.cd()
  dicHists[ strKeyData ][ "fitres" ][ "frame" ].Draw()
  
  # -- Adding label into the legend
  strNameCurve = dicHists[ strKeyData ][ "fitres" ][ "graph" ]
  curveFit = dicHists[ strKeyData ][ "fitres" ][ "frame" ].findObject(strNameCurve)
  leg.AddEntry(curveFit, " ", "l")
  
  strNameHisto = dicHists[ strKeyData ][ "fitres" ][ "histo" ]
  histRooFit = dicHists[ strKeyData ][ "fitres" ][ "frame" ].findObject(strNameHisto)
  leg.AddEntry(histRooFit, dicHists[ strKeyData ][ "label" ], "p")
  
  for strKey in dicHists.keys():
    if type(dicHists[ strKey ]) != types.DictType or "type" not in dicHists[ strKey ] : continue
    if dicHists[ strKey ][ "type" ] != strType : continue
    
    frameX = dicHists[ strKey ][ "fitres" ][ "frame" ]
    
    frameX.SetMarkerStyle(ROOT.kStar)
    
    canvasMain.cd()
    frameX.Draw("same")
    
    # -- Adding label into the legend
    strNameCurve = dicHists[ strKey ][ "fitres" ][ "graph" ]
    curveFit = dicHists[ strKey ][ "fitres" ][ "frame" ].findObject(strNameCurve)
    leg.AddEntry(curveFit, " ", "l")
    
    strNameHisto = dicHists[ strKey ][ "fitres" ][ "histo" ]
    histRooFit = dicHists[ strKey ][ "fitres" ][ "frame" ].findObject(strNameHisto)
    leg.AddEntry(histRooFit, dicHists[ strKey ][ "label" ], "p")

    leg.Draw("same")
  
  strFilenameFits = "fres/fitting_plots_%s_%s.png"%(dicHists["rootfilename"], strType)
  canvasMain.SaveAs(strFilenameFits)
  print strFilenameFits

################################################################
##  Prepare to draw the linear plot
################################################################
dBinMinPeak = 164.0
dBinMaxPeak = 184.0
dSizeBin = 0.1

dDatMinPeak =  1048576.0
dDatMaxPeak = -1048576.0

for strKey in dicHists.keys() :
  if type(dicHists[ strKey ]) != types.DictType or "type" not in dicHists[ strKey ] : continue
  if dicHists[ strKey ][ "type" ] not in ["TT_onlytt", "TT_withbkg"] : continue
  
  dDatVal = dicHists[ strKey ][ "fitres" ][ "peak_val" ]
  dDatErr = dicHists[ strKey ][ "fitres" ][ "peak_err" ]
  
  if dDatMinPeak > dDatVal - dDatErr : dDatMinPeak = dDatVal - dDatErr
  if dDatMaxPeak < dDatVal + dDatErr : dDatMaxPeak = dDatVal + dDatErr

dDatMean = ( dDatMaxPeak + dDatMinPeak ) / 2
dDatMinPeak = dDatMinPeak - ( dDatMean - dDatMinPeak ) * 1.0
dDatMaxPeak = dDatMaxPeak + ( dDatMaxPeak - dDatMean ) * 1.0
print "Range in linear plot : ", dDatMinPeak, dDatMean, dDatMaxPeak

dTan = 0.0
dYP  = 0.0

################################################################
##  Plotting the calibration curve by RooFit
################################################################
histPeak_nb = ROOT.TH1D("histPeak_onlyTT", "Peaks (no bkg)", 
  int(( dBinMaxPeak - dBinMinPeak ) / dSizeBin), dBinMinPeak, dBinMaxPeak)
histPeak_bk = ROOT.TH1D("histPeak_withbkg", "Peaks", 
  int(( dBinMaxPeak - dBinMinPeak ) / dSizeBin), dBinMinPeak, dBinMaxPeak)

listDicPeak = [
  {"name":"no_bkg", "hist":histPeak_nb, "type":"TT_onlytt", 
    "label":"Without background", "marker":ROOT.kMultiply, "color":ROOT.kBlue}, 
  {"name":"bkg", "hist":histPeak_bk, "type":"TT_withbkg", 
    "label":"With background", "marker":ROOT.kCircle, "color":ROOT.kRed}
]

for dicPeak in listDicPeak : 
  print "##### Fitting by RooFit (%s)"%(dicPeak[ "name" ])
  
  histPeak = dicPeak[ "hist" ]
  
  histPeak.SetLineColor(1)

  for strKey in dicHists.keys() :
    if type(dicHists[ strKey ]) != types.DictType or "type" not in dicHists[ strKey ] : continue
    if dicHists[ strKey ][ "type" ] != dicPeak[ "type" ] : continue
    
    dMass = dicHists[ strKey ][ "mass" ]
    dDatVal = dicHists[ strKey ][ "fitres" ][ "peak_val" ]
    dDatErr = dicHists[ strKey ][ "fitres" ][ "peak_err" ]
    
    print dMass, "(" + strKey + ")", dDatVal, dDatErr
    
    nIdxBin = int(( dMass - dBinMinPeak ) / dSizeBin)
    histPeak.SetBinContent(nIdxBin, dDatVal)
    histPeak.SetBinError(nIdxBin, dDatErr)

  histPeak.SetMinimum(dDatMinPeak)
  histPeak.SetMaximum(dDatMaxPeak)

  tf1 = ROOT.TF1("f1_peaks_RooFit_" + dicPeak[ "name" ], "pol1", dBinMinPeak, dBinMaxPeak)
  histPeak.Fit(tf1)
  
  if dicPeak[ "name" ] == "bkg" : 
    dTan = tf1.GetParameter("p1")
    dYP  = tf1.GetParameter("p0")

  setDefAxis(histPeak.GetXaxis(), "M_{top} [GeV/c^2]", 1)
  setDefAxis(histPeak.GetYaxis(), x_name, 0.9)

  rootFits.cd()
  histPeak.Draw()
  histPeak.Write()

if dTan != 0 : 
    print "Top mass = %lf +/- %lf"%(( dicHists[ strKeyData ][ "fitres" ][ "peak_val" ] - dYP ) / dTan, 
      dicHists[ strKeyData ][ "fitres" ][ "peak_err" ] / dTan)
else : 
    print "Top mass = (Unknown because of zeto-tan)"

################################################################
##  Drawing all calibration curves together
################################################################
canvasMain.cd()

# -- Making a legend
leg = ROOT.TLegend(0.64,0.79,0.64+0.30,0.93)
leg.SetBorderSize(0)
#leg.SetNColumns(2)
leg.SetTextSize(0.040)
leg.SetTextFont(42)
leg.SetLineColor(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)

strSign = ""

for dicPeak in listDicPeak : 
  histPeak = dicPeak[ "hist" ]
  
  histPeak.SetMarkerStyle(dicPeak[ "marker" ])
  histPeak.SetMarkerColor(dicPeak[ "color" ])
  try : 
    histPeak.GetFunction("f1_peaks_RooFit_" + dicPeak[ "name" ]).SetLineColor(dicPeak[ "color" ])
  except : 
    print dicPeak, histPeak.GetFunction("f1_peaks_RooFit_" + dicPeak[ "name" ])
  
  ROOT.gStyle.SetOptFit(0)
  histPeak.Draw(strSign + "E1")
  strSign = "same"

  leg.AddEntry(histPeak, dicPeak[ "label" ], "lp")

polyData = ROOT.TH1D("histPeak_data", "Data", 
  int(( dBinMaxPeak - dBinMinPeak ) / dSizeBin), dBinMinPeak, dBinMaxPeak)

for i in range(int(( dBinMaxPeak - dBinMinPeak ) / dSizeBin) + 1) : 
  polyData.SetBinContent(i, dicHists[ strKeyData ][ "fitres" ][ "peak_val" ])
  polyData.SetBinError(i, dicHists[ strKeyData ][ "fitres" ][ "peak_err" ])

polyData.SetFillColorAlpha(ROOT.kBlack, 0.3)
polyData.SetMarkerStyle(ROOT.kDot)

polyData.Draw("sameE3")

leg.AddEntry(polyData, "Data", "f")

leg.Draw("same")

canvasMain.SaveAs("fres/calibration_cuve_%s.png"%(dicHists["rootfilename"]))

################################################################
##  Everything is over; closing the file
################################################################

rootFits.Write()
rootFits.Close()
rootHists.Close()

print ''
print ''
