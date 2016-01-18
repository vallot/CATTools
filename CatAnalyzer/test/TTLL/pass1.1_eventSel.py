#!/usr/bin/env python

import sys, os
outDir = "pass1"
if os.path.exists(outDir):
    print "Directory already exists. remove or rename it."
    sys.exit()
os.mkdir(outDir)
os.system("ln -s ../analyze_sig_cfg.py %s/analyze_sig_cfg.py" % outDir)
os.system("ln -s ../analyze_bkg_cfg.py %s/analyze_bkg_cfg.py" % outDir)
os.system("ln -s ../analyze_data_cfg.py %s/analyze_data_cfg.py" % outDir)
os.system("ln -s ../customise_saveEvent_cfg.py %s/customise_saveEvent_cfg.py" % outDir)

## Load JSON file and categorize datasets
import json
dataDir = "%s/src/CATTools/CatAnalyzer/data" % os.environ["CMSSW_BASE"]
js = json.loads(open("%s/dataset.json" % dataDir).read())

def isBlacklisted(name):
    for x in [
        "GluGluToZZ", "HToMuMu", "WpWp", "WW_dps",
        "WWTo2L2Nu_powheg", "WZTo", "ZZTo", "ZZto",
        "SingleElectron_Run2015", "SingleMuon_Run2015", ]:
        if x in name: return True
    return False

sigList = []
bkgList = []
datList = []

for j in js:
    datasetName = j['DataSetName']
    name = j['name']
    title = j['title']
    if datasetName.endswith("/MINIAOD"):
        if datasetName.startswith("/Jet"): continue
        datList.append(j)
    elif title.startswith('t#bar{t}'): sigList.append(j)
    else:
        if datasetName.startswith("/QCD"): continue
        bkgList.append(j)

## Write script to run create-batch

## Start from the central samples
out = open("%s/submit_dat_central.sh" % outDir, "w")
for d in datList:
    name = d['name']
    if isBlacklisted(name): continue
    submitCmd  = "create-batch --cfg analyze_data_cfg.py --maxFiles 100"
    submitCmd += " --jobName %s --fileList %s/dataset/dataset_%s.txt" % (name, dataDir, name)

    print>>out, (submitCmd + " --jobName %s/central --customise customise_saveEvent_cfg.py" % name)
out.close()

out = open("%s/submit_bkg_central.sh" % outDir, "w")
for d in bkgList:
    name = d['name']
    if isBlacklisted(name): continue
    submitCmd  = "create-batch --cfg analyze_bkg_cfg.py --maxFiles 100"
    submitCmd += " --jobName %s --fileList %s/dataset/dataset_%s.txt" % (name, dataDir, name)

    print>>out, (submitCmd + " --jobName %s/central --customise customise_saveEvent_cfg.py" % name)

for d in sigList: ## TTbar others are treated as background
    name = d['name']
    if isBlacklisted(name): continue
    submitCmd  = "create-batch --cfg analyze_sig_cfg.py --maxFiles 100"
    submitCmd += " --fileList %s/dataset/dataset_%s.txt" % (dataDir, name)

    arg0 = 'filterPartonTTLL.invert=True'
    print>>out, (submitCmd + " --jobName %s_Others/central --customise customise_saveEvent_cfg.py --args '%s'" % (name, arg0))
out.close()

out = open("%s/submit_sig_central.sh" % outDir, "w")
for d in sigList:
    name = d['name']
    if isBlacklisted(name): continue
    submitCmd  = "create-batch --cfg analyze_sig_cfg.py --maxFiles 100"
    submitCmd += " --fileList %s/dataset/dataset_%s.txt" % (dataDir, name)

    arg0 = 'filterPartonTTLL.invert=False'
    print>>out, (submitCmd + " --jobName %s_LL/central --customise customise_saveEvent_cfg.py --args '%s'" % (name, arg0))
out.close()




## Then continue to the systematic variations
out_sig = open("%s/submit_sig_unc.sh" % outDir, "w")
out_bkg = open("%s/submit_bkg_unc.sh" % outDir, "w")
out_dat = open("%s/submit_dat_unc.sh" % outDir, "w")

## Loop over common systematics
systAny = {
    'mu_pt/up':'eventsTTLL.muon.scaleDirection=1',
    'mu_pt/dn':'eventsTTLL.muon.scaleDirection=-1',
    'el_pt/up':'eventsTTLL.electron.scaleDirection=1',
    'el_pt/dn':'eventsTTLL.electron.scaleDirection=-1',
    'jet_cor/up':'eventsTTLL.jet.scaleDirection=1',
    'jet_cor/dn':'eventsTTLL.jet.scaleDirection=-1',
}
for systName in systAny:
    syst = systAny[systName]

    for d in datList:
        name = d['name']
        if isBlacklisted(name): continue
        submitCmd  = "create-batch --cfg analyze_data_cfg.py --maxFiles 100"
        submitCmd += " --fileList %s/dataset/dataset_%s.txt" % (dataDir, name)
        print>>out_dat, (submitCmd + (" --jobName %s/%s --args '%s'" % (name, systName, syst)))

    for d in bkgList:
        name = d['name']
        if isBlacklisted(name): continue
        submitCmd  = "create-batch --cfg analyze_bkg_cfg.py --maxFiles 100"
        submitCmd += " --fileList %s/dataset/dataset_%s.txt" % (dataDir, name)
        print>>out_bkg, (submitCmd + (" --jobName %s/%s --args '%s'" % (name, systName, syst)))

    for d in sigList:
        name = d['name']
        if '_scaleup' in name or '_scaledown' in name: continue ## Skip this variations for scale up/down samples
        if isBlacklisted(name): continue
        submitCmd  = "create-batch --cfg analyze_sig_cfg.py --maxFiles 100"
        submitCmd += " --fileList %s/dataset/dataset_%s.txt" % (dataDir, name)

        print>>out_sig, (submitCmd + (" --jobName %s_LL/%s --args '%s'" % (name, systName, syst)))
        print>>out_bkg, (submitCmd + (" --jobName %s_Others/%s --args '%s filterPartonTTLL.invert=True'" % (name, systName, syst)))

## Then loop over MC specific systematics
systMC = {
    'jet_res/up':'eventsTTLL.jet.resolDirection=1',
    'jet_res/dn':'eventsTTLL.jet.resolDirection=-1',
    'pileup/up':'eventsTTLL.vertex.pileupWeight="pileupWeight:up"',
    'pileup/dn':'eventsTTLL.vertex.pileupWeight="pileupWeight:dn"',
    'mu_eff/up':'eventsTTLL.muon.efficiencySFDirection=1',
    'mu_eff/dn':'eventsTTLL.muon.efficiencySFDirection=-1',
    'el_eff/up':'eventsTTLL.electron.efficiencySFDirection=1',
    'el_eff/dn':'eventsTTLL.electron.efficiencySFDirection=-1',
}
for systName in systMC:
    syst = systMC[systName]

    for d in bkgList:
        name = d['name']
        if isBlacklisted(name): continue
        submitCmd  = "create-batch --cfg analyze_bkg_cfg.py --maxFiles 100"
        submitCmd += " --fileList %s/dataset/dataset_%s.txt" % (dataDir, name)
        print>>out_bkg, (submitCmd + (" --jobName %s/%s --args '%s'" % (name, systName, syst)))

    for d in sigList:
        name = d['name']
        if '_scaleup' in name or '_scaledown' in name: continue ## Skip this variations for scale up/down samples
        if isBlacklisted(name): continue
        submitCmd  = "create-batch --cfg analyze_sig_cfg.py --maxFiles 100"
        submitCmd += " --fileList %s/dataset/dataset_%s.txt" % (dataDir, name)

        print>>out_sig, (submitCmd + (" --jobName %s_LL/%s --args '%s'" % (name, systName, syst)))
        print>>out_bkg, (submitCmd + (" --jobName %s_Others/%s --args '%s filterPartonTTLL.invert=True'" % (name, systName, syst)))

## Scale up/down systematic uncertainty from LHE weight
## This uncertainty have to be combined with envelope
## Let us assume index1-10 are for the scale variations (muF & muR)
for i in range(1,9): # total 8 scale variations, 3 muF x 3 muR and one for central weight
    systName = "gen_scale/%d" % i
    syst = 'eventsTTLL.genWeight.src="genWeight:pdfWeights" eventsTTLL.genWeight.index=%d' % i
    for d in bkgList:
        name = d['name']
        if isBlacklisted(name): continue
        submitCmd  = "create-batch --cfg analyze_bkg_cfg.py --maxFiles 100"
        submitCmd += " --fileList %s/dataset/dataset_%s.txt" % (dataDir, name)
        print>>out_bkg, (submitCmd + (" --jobName %s/%s --args '%s'" % (name, systName, syst)))

    ## Modification for ttbar samples
    syst += ' agen.weight="genWeight:pdfWeights" agen.weightIndex=%d' % i
    for d in sigList:
        name = d['name']
        if isBlacklisted(name): continue
        submitCmd  = "create-batch --cfg analyze_sig_cfg.py --maxFiles 100"
        submitCmd += " --fileList %s/dataset/dataset_%s.txt" % (dataDir, name)

        print>>out_sig, (submitCmd + (" --jobName %s_LL/%s --args '%s'" % (name, systName, syst)))
        print>>out_bkg, (submitCmd + (" --jobName %s_Others/%s --args '%s filterPartonTTLL.invert=True'" % (name, systName, syst)))

## PDF weights
## 110 variations for aMC@NLO, 248 for POWHEG
## NOTE: there is weight vector, but we don't do it for LO generator here.
## -- for note madgraph weightSize = 445 + 8 variations
for d in bkgList:
    name = d['name']
    if isBlacklisted(name): continue
    if '_aMC' in name: weightSize = 110
    elif '_powheg' in name: weightSize = 248
    else: continue

    submitCmd  = "create-batch --cfg analyze_sig_cfg.py --maxFiles 100"
    submitCmd += " --fileList %s/dataset/dataset_%s.txt" % (dataDir, name)

    for i in range(9, weightSize+1):
        arg = 'eventsTTLL.genWeight.src="genWeight:pdfWeights" eventsTTLL.genWeight.index=%d' % i
        print>>out_bkg, (submitCmd + (" --jobName %s/gen_PDF/%d --args '%s'" % (name, i, arg)))

for d in sigList:
    name = d['name']
    if isBlacklisted(name): continue
    if '_scale' in name: continue
    if '_aMC' in name: weightSize = 110
    elif '_powheg' in name: weightSize = 248
    else: continue

    submitCmd  = "create-batch --cfg analyze_sig_cfg.py --maxFiles 100"
    submitCmd += " --fileList %s/dataset/dataset_%s.txt" % (dataDir, name)

    for i in range(9, weightSize+1):
        arg  = 'eventsTTLL.genWeight.src="genWeight:pdfWeights" eventsTTLL.genWeight.index=%d' % i
        arg += ' agen.weight="genWeight:pdfWeights" agen.weightIndex=%d' % i
        print>>out_sig, (submitCmd + (" --jobName %s_LL/gen_PDF/%d --args '%s'" % (name, i, arg)))
        print>>out_bkg, (submitCmd + (" --jobName %s_Others/gen_PDF/%d --args '%s filterPartonTTLL.invert=True'" % (name, i, arg)))

out_sig.close()
out_bkg.close()
out_dat.close()

os.system("chmod +x %s/submit_sig_central.sh" % outDir)
os.system("chmod +x %s/submit_bkg_central.sh" % outDir)
os.system("chmod +x %s/submit_dat_central.sh" % outDir)
os.system("chmod +x %s/submit_sig_unc.sh" % outDir)
os.system("chmod +x %s/submit_bkg_unc.sh" % outDir)
os.system("chmod +x %s/submit_dat_unc.sh" % outDir)

print "Prepared to submit cmsRun jobs."
print "Submit jobs using helper script under pass1 directory yourself."
print "> cd pass1"
print "> ./submit_dat_central.sh"
print "> ./submit_sig_central.sh"
print "> ./submit_bkg_central.sh"
print ">"
print "> ./submit_dat_unc.sh"
print "> ./submit_sig_unc.sh"
print "> ./submit_bkg_unc.sh"
print ""
print "or, use the xargs magic to submit them all"
print "> cat submit*.sh | sed -e 's;create-batch;;g' | xargs -L1 -P20 create-batch"
print ""

