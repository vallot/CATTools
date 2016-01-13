#!/usr/bin/env python

import sys, os
sys.argv.append("-b")
from ROOT import *

srcbase = "pass3/hists"
outbase = "pass3/result"

os.makedirs(outbase)

meas = RooStats.HistFactory.Measurement("xsec", "xsec")

meas.SetOutputFilePrefix("%s/meas" % outbase)
meas.SetPOI("xsecOverSM")
meas.AddConstantParam("Lumi")
meas.AddConstantParam("alpha_systTTLL")

meas.SetLumi(2.11*1000)
meas.SetLumiRelErr(0.10)
meas.SetExportOnly(False)

fRD = "pass3/hists/central/Data.root"

infoMCs = {
    "TTLL":"t_bar_t__Jets_LL",
    "TTOthers":"t_bar_t__Jets_Others",
    "SingleTop":"Single_top",
    "Dibosons":"Dibosons",
    "Tribosons":"Tribosons",
    "WJets":"W_Jets",
    "DYJets":"Z__gamma_rightarrow_ll",
}

chan = {}
samples = {}
for mode in ("ee", "mm", "em"):
    #hFullPath = "eventsTTLL/%s/step5/event_st" % mode
    hFullPath = "eventsTTLL/%s/step4/jet1_btag" % mode
    #hFullPath = "eventsTTLL/%s/step4/z_m_noveto" % mode

    chan[mode] = RooStats.HistFactory.Channel(mode)
    hName = os.path.basename(hFullPath)
    hPath = os.path.dirname(hFullPath) + "/"
    chan[mode].SetData(hName, fRD, hPath)
    chan[mode].SetStatErrorConfig(0.10, "Poisson")

    for name in infoMCs:
        fPrefix = infoMCs[name]
        samples[mode+name] = RooStats.HistFactory.Sample(name, hName, "pass3/hists/central/%s.root" % fPrefix, hPath)
        if name == "TTLL":
            samples[mode+name].AddNormFactor("xsecOverSM", 1, 0.85, 1.05)
        else:
            #samples[mode+name].AddOverallSys("syst_%s" % name, 0.95, 1.05) ## xsec error
            #samples[mode+name].ActivateStatError()
            pass

        ## Systematic uncertainty variations
        for syst in os.listdir(srcbase):
            fNameDn = "%s/%s/dn/%s.root" % (srcbase, syst, fPrefix)
            fNameUp = "%s/%s/up/%s.root" % (srcbase, syst, fPrefix)
            if not os.path.exists(fNameDn) or not os.path.exists(fNameUp): continue

            samples[mode+name].AddHistoSys("syst_"+syst, hName, fNameDn, hPath, hName, fNameUp, hPath)

        chan[mode].AddSample(samples[mode+name])

    meas.AddChannel(chan[mode])

meas.CollectHistograms()
meas.PrintTree()

RooStats.HistFactory.MakeModelAndMeasurementFast(meas)
