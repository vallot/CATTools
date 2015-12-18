#!/usr/bin/env python

from ROOT import *
import json
import sys, os
import string

srcbase, outbase = "pass1", "pass2"
jsdir = "%s/src/CATTools/CatAnalyzer/data" % os.environ["CMSSW_BASE"]

## Build (filesystem safe) title -> dataset info list mapping
info = {}
for x in json.loads(open(jsdir+"/dataset.json").read()):
    title = x["title"]
    name = x["name"]
    xsec = x["xsec"]

    safeTitle = list(title)
    validChars = list(string.ascii_letters) + list(string.digits) + ['-_.']
    for i, c in enumerate(safeTitle):
        if c not in validChars: safeTitle[i] = '_'
    safeTitle = ''.join(safeTitle)

    info[name] = {'title':title, 'safetitle':safeTitle, 'xsec':xsec, 'DataSetName':x["DataSetName"]}

## Get directory structure to reorganize input
mergeList = {}
for path, dirs, files in os.walk(srcbase):
    rootFiles = [x for x in files if x.endswith('.root')]
    if len(rootFiles) == 0: continue

    ## Split sample name and uncertainty variations from the path name
    pp = path[len(srcbase)+1:].split('/')
    name = pp[0]

    ## Special care for ttbar signal samples due to gen level splitting
    nameInInfo = name
    if name.startswith("TT"):
        if name.endswith("LL"): nameInInfo = name[:-2]
        elif name.endswith("LLOthers"): nameInInfo = name[:-8]

    xsec = info[nameInInfo]['xsec']
    safeTitle = info[nameInInfo]['safetitle']
    datasetName = info[nameInInfo]['DataSetName']
    if datasetName.endswith("AOD"): type = "data"
    else: type = "MC"

    ## Put back suffix
    if name.startswith("TT"):
        if name.endswith("LL"): safeTitle += "_LL"
        elif name.endswith("LLOthers"): safeTitle += "_LLOthers"

    for f in rootFiles:
        outpath = [outbase]
        if len(pp) > 1: outpath.extend(pp[1:])
        outpath.extend([f[:-5], safeTitle+'.root'])
        outfile = '/'.join(outpath)

        if outfile not in mergeList:
            mergeList[outfile] = {'type':type, 'samples':[],}

        mergeList[outfile]['samples'].append({
            'xsec':xsec,
            'file':os.path.join(path, f),
        })

## Start mergin
for outfileName in mergeList:
    ## Prepare output directory
    outdir = os.path.dirname(outfileName)
    mergeInfo = mergeList[outfileName]
    type = mergeInfo['type']
    samples = mergeInfo['samples']

    if not os.path.exists(outdir): os.makedirs(outdir)

    print "@@ Producing plots for the sample", outfileName
    fout = TFile(outfileName, "RECREATE")
    out_moddir = fout.mkdir("ttll")
    for i, x in enumerate(samples):
        fsrc = TFile(x['file'])

        scale = 1.0
        if type == 'MC':
            ## Get the sum of weights
            hWeight = fsrc.Get("ttll/overall/weight")
            sumW = hWeight.Integral()*hWeight.GetMean()
            if sumW != 0.0: scale = x['xsec']/sumW

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
                    hout.Write()
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
                        hout.Write()

