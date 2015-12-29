#!/usr/bin/env python

from ROOT import *
import json
import sys, os
import string
from multiprocessing import Pool

srcbase, outbase = "pass1", "pass2"

def merge(outfileName, mergeInfo):
    print "@@ Producing plots for the sample", outfileName

    type = mergeInfo['type']
    samples = mergeInfo['samples']
    fout = TFile(outfileName, "RECREATE")
    out_moddir = fout.mkdir("ttll")

    for i, x in enumerate(samples):
        fsrc = TFile(x['file'])

        scale = 1.0
        if type == 'MC':
            ## Get the sum of weights
            if fsrc.GetDirectory("agen") != None:
                hWeight = fsrc.Get("agen/hWeight")
                ## Get weight distribution before gen filter - only if exists
            else:
                hWeight = fsrc.Get("ttll/overall/weight")
            if hWeight != None:
                sumW = hWeight.Integral()*hWeight.GetMean()
                if sumW != 0.0: scale = x['xsec']/sumW
            else:
                print "Cannot find weight histogram!!! Skipping cross section normalization."

        ## Visit all cut flows and do the normalization
        src_moddir = fsrc.Get("ttll")
        for chName in [x.GetName() for x in src_moddir.GetListOfKeys()]:
            if chName == "overall": continue

            if i == 0: out_moddir.mkdir(chName)
            src_chdir = src_moddir.Get(chName)
            out_chdir = out_moddir.Get(chName)

            for stepName in [x.GetName() for x in src_chdir.GetListOfKeys()]:
                src_stepdir = src_chdir.Get(stepName)
                ## Cut flow table
                if src_stepdir.IsA().InheritsFrom("TH1"):
                    hsrc = src_stepdir
                    hsrc.Scale(scale)

                    out_chdir.cd()
                    if i == 0:
                        hout = hsrc.Clone()
                        hout.Reset()
                    hout = out_chdir.Get(stepName)
                    hout.Add(hsrc)
                    hout.Write("", TObject.kOverwrite)
                ## Control plots under cut step directory
                else:
                    if i == 0: out_chdir.mkdir(stepName)
                    out_stepdir = out_chdir.Get(stepName)
                
                    for histName in [x.GetName() for x in src_stepdir.GetListOfKeys()]:
                        hsrc = src_stepdir.Get(histName)
                        hsrc.Scale(scale)

                        out_stepdir.cd()
                        if i == 0:
                            hout = hsrc.Clone()
                            hout.Reset()
                        hout = out_stepdir.Get(histName)
                        hout.Add(hsrc)
                        hout.Write("", TObject.kOverwrite)

if __name__ == '__main__':
    p = Pool(30)

    mergeList = json.loads(open(outbase+"/samples.json").read())

    ## Start mergin
    for outfileName in mergeList:
        ## Prepare output directory
        outdir = os.path.dirname(outfileName)
        mergeInfo = mergeList[outfileName]

        if not os.path.exists(outdir): os.makedirs(outdir)

        p.apply_async(merge, [outfileName, mergeInfo])

    p.close()
    p.join()
