#!/usr/bin/env python

from ROOT import *
import json
import sys, os

jsdir = "%s/src/CATTools/CatAnalyzer/data" % os.environ["CMSSW_BASE"]
info = {}
for x in json.loads(open(jsdir+"/dataset.json").read()):
    info[x["name"]] = x

srcdir = "pass1_makeHist"
for srcpath, dirs, files in os.walk(srcdir):
    rootFiles = [x for x in files if x.endswith('.root')]
    if len(rootFiles) == 0: continue

    ## Double check with dataset JSON
    name = srcpath.split('/')[1]
    if name not in info:
        print "Cannot find dataset", name, "from the dataset JSON"
        continue
    dsName = info[name]["DataSetName"]
    if dsName.endswith("AOD"): continue # skip real data

    xsec = info[name]["xsec"]

    ## Prepare output directory
    outpath = srcpath.replace(srcdir, 'pass2_mergeHist')
    if not os.path.exists(outpath): os.makedirs(outpath)
    for file in files:
        print "@@ Producing plots for the sample", srcpath
        fsrc = TFile(os.path.join(srcpath, file))
        fout = TFile(os.path.join(outpath, file), "RECREATE")

        ## Get the sum of weights
        hWeight = fsrc.Get("ttll/overall/weight")
        sumW = hWeight.Integral()*hWeight.GetMean()

        ## Visit all cut flows and do the normalization
        src_moddir = fsrc.Get("ttll")
        out_moddir = fout.mkdir("ttll")
        for chName in [x.GetName() for x in src_moddir.GetListOfKeys()]:
            if chName == "overall": continue
            src_chdir = src_moddir.Get(chName)
            out_chdir = out_moddir.mkdir(chName)
            for stepName in [x.GetName() for x in src_chdir.GetListOfKeys()]:
                src_stepdir = src_chdir.Get(stepName)
                ## Cut flow table
                if src_stepdir.IsA().InheritsFrom("TH1"):
                    out_chdir.cd()
                    hout = src_stepdir.Clone()
                    if sumW !=0: hout.Scale(xsec/sumW)
                    hout.Write()
                ## Control plots under cut step directory
                else:
                    out_chdir.mkdir(stepName).cd()
                
                    for histName in [x.GetName() for x in src_stepdir.GetListOfKeys()]:
                        hout = src_stepdir.Get(histName).Clone()
                        if sumW != 0: hout.Scale(xsec/sumW)
                        hout.Write()

