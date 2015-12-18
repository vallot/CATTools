#!/usr/bin/env python

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

    safeTitle = list(title)
    validChars = list(string.ascii_letters) + list(string.digits) + ['-_.']
    for i, c in enumerate(safeTitle):
        if c not in validChars: safeTitle[i] = '_'
    safeTitle = ''.join(safeTitle)

    info[name] = {
        'title':title, 'safetitle':safeTitle, 
        'xsec':x['xsec'], 'DataSetName':x["DataSetName"], "colour":x["colour"]}

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
    color = info[nameInInfo]['colour']
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
            mergeList[outfile] = {'type':type, 'color':color, 'samples':[],}

        mergeList[outfile]['samples'].append({
            'xsec':xsec,
            'file':os.path.join(path, f),
        })

if not os.path.exists(outbase): os.mkdir(outbase)

js = json.dumps(mergeList, indent=4, sort_keys=True)
js = js.replace("},", "},\n")
f = open("%s/samples.json" % outbase, "w")
f.write(js)
f.close()

print "%s/samples.json is created." % outbase
print "You can edit json file manually if needed"
