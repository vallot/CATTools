#!/usr/bin/env python

import json
import sys, os
import string

blacklist = [
    "ZZTo", "ZZto", "WWTo", "WZTo", "GluGluTo", "WpWp", "WW_dps",
    "HToMuMu",
]
srcbase, outbase = "pass1", "pass2"

jsdir = "%s/src/CATTools/CatAnalyzer/data/dataset" % os.environ["CMSSW_BASE"]

## Build (filesystem safe) title -> dataset info list mapping
info = {}
blacklisted = []
for x in json.loads(open(jsdir+"/dataset.json").read()):
    title = x["title"]
    name = x["name"]

    if any([b in name for b in blacklist]):
        blacklisted.append(name)
        continue

    safeTitle = list(title)
    validChars = list(string.ascii_letters) + list(string.digits) + ['-_.']
    for i, c in enumerate(safeTitle):
        if c not in validChars: safeTitle[i] = '_'
    safeTitle = ''.join(safeTitle)

    info[name] = {
        'title':title, 'safetitle':safeTitle, 
        'xsec':x['xsec'], 'DataSetName':x["DataSetName"], "colour":x["colour"]}

## Get directory structure to reorganize input
srcList = {}
mergeList = {}
for path, dirs, files in os.walk(srcbase):
    rootFiles = [x for x in files if x.endswith('.root') and not x.startswith('out_')]
    if len(rootFiles) == 0: continue

    pp = path.split('/')[1:]
    if len(pp) == 0: continue
    if pp[0] not in srcList: srcList[pp[0]] = {}
    if len(pp) == 1: srcList[pp[0]]['central'] = rootFiles
    elif len(pp) > 1: srcList[pp[0]][pp[1]] = rootFiles

    ## Split sample name and uncertainty variations from the path name
    pp = path[len(srcbase)+1:].split('/')
    name = pp[0]

    ## Special care for ttbar signal samples due to gen level splitting
    nameInInfo = name
    if nameInInfo in blacklisted: continue
    if name.startswith("TT"):
        if name.endswith("_LL"): nameInInfo = name[:-3]
        elif name.endswith("_Others"): nameInInfo = name[:-7]

    xsec = info[nameInInfo]['xsec']
    color = info[nameInInfo]['colour']
    safeTitle = info[nameInInfo]['safetitle']
    datasetName = info[nameInInfo]['DataSetName']
    if datasetName.endswith("AOD"): type = "data"
    else: type = "MC"

    ## Put back suffix
    if name.startswith("TT"):
        if name.endswith("_LL"): safeTitle += "_LL"
        elif name.endswith("_Others"): safeTitle += "_Others"

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
f = open("%s/samples.json" % outbase, "w")
f.write(js)
f.close()

js = json.dumps(srcList, indent=4, sort_keys=True)
f = open("%s/samples.json" % srcbase, "w")
f.write(js)
f.close()

print "NOTE: you can give blacklist by editing this script"
print "Following samples are blacklisted"
print " ","\n  ".join(blacklisted)
print "%s/samples.json is created." % outbase
print "You can edit json file manually if needed"
print "WARNING!!! some samples are known to be overapping. Please remove them manually!!!"
