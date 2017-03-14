#!/usr/bin/env python

## Define list of samples and their recipes
samples = {
    'data':[
        "SingleElectron_Run2016B", "SingleElectron_Run2016C", "SingleElectron_Run2016D", "SingleElectron_Run2016E", "SingleElectron_Run2016F", "SingleElectron_Run2016G", "SingleElectron_Run2016H_v2", "SingleElectron_Run2016H_v3",
        "SingleMuon_Run2016B", "SingleMuon_Run2016C", "SingleMuon_Run2016D", "SingleMuon_Run2016E", "SingleMuon_Run2016F", "SingleMuon_Run2016G", "SingleMuon_Run2016H_v2", "SingleMuon_Run2016H_v3",
    ],
    'bkg':[
        "DYJets", "DYJets_10to50",
        "SingleTop_s", "SingleTop_t", "SingleTbar_t", "SingleTop_tW", "SingleTbar_tW", 
        "WJets",
        "WW", "WZ", "ZZ", 
        "WWW", "WWZ", "WZZ", "ZZZ", 
    ],
    'sig':["TTLJ_powheg",],
    'sig.ttOthers':["TT_powheg",],
    'sig.ttNoFilter':["ttW", "ttZ",],
}

import sys, os
outDir = "pass1"
if os.path.exists(outDir):
    print "Directory already exists. remove or rename it."
    sys.exit()
os.mkdir(outDir)
os.system("ln -s ../analyze_sig_cfg.py %s/analyze_sig_cfg.py" % outDir)
os.system("ln -s ../analyze_bkg_cfg.py %s/analyze_bkg_cfg.py" % outDir)
os.system("ln -s ../analyze_data_cfg.py %s/analyze_data_cfg.py" % outDir)

## Load JSON file and categorize datasets
import json
dataDir = "%s/src/CATTools/CatAnalyzer/data/dataset" % os.environ["CMSSW_BASE"]
js = json.loads(open("%s/dataset.json" % dataDir).read())

## Write script to run create-batch
cmds = {}
for type, sampleNames in samples.iteritems():
    modifiers = []
    if '.' in type:
        modifiers = type.split('.')[1:]
        type = type.split('.')[0]

    suffix = ""
    baseargs = []
    if 'ttOthers' in modifiers:
        baseargs += ["filterParton.invert=True"]
        suffix = ".Others"
    elif 'ttNoFilter' in modifiers:
        baseargs += ['filterParton.nLepton=-1']

    for name in sampleNames:
        jobName = "%s%s" % (name, suffix)

        submitCmd  = "create-batch --cfg analyze_%s_cfg.py --maxFiles 25 " % type
        submitCmd += " --jobName %s --fileList %s/dataset_%s.txt " % (jobName, dataDir, name)

        args = baseargs
        if len(args) > 0: submitCmd += (" --args '%s'" % ' '.join(args))

        cmds[jobName] = submitCmd

out = open("%s/submit.sh" % outDir, "w")
for cmd in sorted([x for x in cmds]): print>>out, cmds[cmd]
out.close()

os.system("chmod +x %s/submit.sh" % outDir)

print "Prepared to submit cmsRun jobs."
print "Submit jobs using helper script under pass1 directory yourself."
print "> cd pass1"
print "> ./submit.sh"
print ""
print "or, use the xargs magic to submit them all"
print "> cat submit.sh | sed -e 's;create-batch;;g' -e 's; *$;;g' | xargs -l -P20 create-batch"
print ""

