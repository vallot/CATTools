#!/usr/bin/env python

## Define recipes
recipes = {
    'nominal':{'nominal':[]},
    'antiIso':{'antiIso':['mu.muon.applyAntiIso=True', 'el.electron.applyAntiIso=True']},
    'syst':{
        'mu_pt/up':['mu.muon.scaleDirection=1', 'el.muon.scaleDirection=1'],
        'mu_pt/dn':['mu.muon.scaleDirection=-1', 'el.muon.scaleDirection=-1'],
        'el_pt/up':['el.electron.scaleDirection=1', 'el.electron.scaleDirection=1'],
        'el_pt/dn':['el.electron.scaleDirection=-1', 'el.electron.scaleDirection=-1'],
        'jet_cor/up':['mu.jet.scaleDirection=1', 'el.jet.scaleDirection=1'],
        'jet_cor/dn':['mu.jet.scaleDirection=-1', 'el.jet.scaleDirection=-1'],
    },
    'systMC':{
        'jet_res/up':['mu.jet.resolDirection=1', 'el.jet.resolDirection=1'],
        'jet_res/dn':['mu.jet.resolDirection=-1', 'mu.jet.resolDirection=-1'],
    },
}

## Define list of samples and their recipes
samples = {
    'data nominal syst antiIso':[
        ## Real data
        ## cmsRun analyze_dat_cfg.py
        ## do common systematic variations
        "SingleElectron_Run2016B", "SingleElectron_Run2016C", "SingleElectron_Run2016D", "SingleElectron_Run2016E", "SingleElectron_Run2016F", "SingleElectron_Run2016G", "SingleElectron_Run2016H_v2", "SingleElectron_Run2016H_v3",
        "SingleMuon_Run2016B", "SingleMuon_Run2016C", "SingleMuon_Run2016D", "SingleMuon_Run2016E", "SingleMuon_Run2016F", "SingleMuon_Run2016G", "SingleMuon_Run2016H_v2", "SingleMuon_Run2016H_v3",
    ],
    'bkg nominal syst systMC antiIso':[
        ## MC samples without ttbar genFilter
        ## cmsRun analyze_bkg_cfg.py
        ## do common systematic variations
        ## do MC common systematic variations
        ## No ttbar specific gen level analyzer, no ttbar genFilters
        "DYJets", "DYJets_10to50",
        "SingleTop_s", "SingleTop_t", "SingleTbar_t", "SingleTop_tW", "SingleTbar_tW",
        "WJets", "WJets_MG",
        "WW", "WZ", "ZZ",
        "WWW", "WWZ", "WZZ", "ZZZ",
    ],
    'sig nominal syst systMC antiIso':[
        ## MC Signal samples, select ttbar-dilepton at Gen.Level
        ## cmsRun analyze_sig_cfg.py
        ## do common systematic variations
        ## do MC common systematic variations
        ## do ttbar specific gen level analyzer, do ttbar-dilepton genFilters
        ## do PDF variations
        "TTLJ_powheg",
        #"ttbb",
    ],
    'sig.noll nominal syst systMC antiIso':[
        ## ttbar-others
        ## cmsRun analyze_sig_cfg.py
        ## do common systematic variations
        ## do MC common systematic variations
        ## do ttbar specific gen level analyzer, do ttbar-dilepton genFilters
        "TT_powheg",
        #"ttW", "ttZ",
    ],
    'sig nominal':[
        ## Variation on generator choice
        ## produce nominals only
        ## cmsRun analyze_sig_cfg.py
        ## do ttbar specific gen level analyzer, do ttbar-dilepton genFilters
        "TTJets_MG", "TTJets_aMC",
        "TT_powheg_herwig",

        "TT_powheg_FSRdown","TT_powheg_FSRup",
        "TT_powheg_ISRdown","TT_powheg_ISRup",
        "TT_powheg_UEdown","TT_powheg_UEup",
        "TT_powheg_mtop1665","TT_powheg_mtop1695","TT_powheg_mtop1715",
        "TT_powheg_mtop1735","TT_powheg_mtop1755","TT_powheg_mtop1785",
    ],
    'sig.nofilter nominal syst systMC':[
        "TT_AntitopLeptonicDecay_TH_1L3B_Eta_Hct", "TT_TopLeptonicDecay_TH_1L3B_Eta_Hct",
        "TT_AntitopLeptonicDecay_TH_1L3B_Eta_Hut", "TT_TopLeptonicDecay_TH_1L3B_Eta_Hut",
    ],
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
for key in samples:
    type, recipesToRun = key.split()[0], key.split()[1:]
    modifiers = []
    if '.' in type:
        modifiers = type.split('.')[1:]
        type = type.split('.')[0]

    suffix = ""
    baseargs = []
    if 'noll' in modifiers:
        baseargs += ["filterParton.invert=True"]
        suffix = ".Others"
    if 'nofilter' in modifiers:
        baseargs += ["filterParton.nLepton=-1"]

    for name in samples[key]:
        for recipeToRun in recipesToRun:
            for recipe in recipes[recipeToRun]:
                jobName = "%s%s/%s" % (name, suffix, recipe)

                submitCmd  = "create-batch --cfg analyze_%s_cfg.py --maxFiles 25 " % type
                submitCmd += " --jobName %s --fileList %s/dataset_%s.txt " % (jobName, dataDir, name)

                args = baseargs+recipes[recipeToRun][recipe]
                if len(args) > 0: submitCmd += (" --args '%s'" % ' '.join(args))

                cmds[jobName] = submitCmd

out = open("%s/submit_nominal.sh" % outDir, "w")
for cmd in sorted([x for x in cmds if 'nominal' in x and 'antiIso' not in x]): print>>out, cmds[cmd]
out.close()
out = open("%s/submit_unc.sh" % outDir, "w")
for cmd in sorted([x for x in cmds if 'nominal' not in x and 'antiIso' not in x]): print>>out, cmds[cmd]
out.close()
out = open("%s/submit_antiIso.sh" % outDir, "w")
for cmd in sorted([x for x in cmds if 'antiIso' in x]): print>>out, cmds[cmd]
out.close()

os.system("chmod +x %s/submit_nominal.sh" % outDir)
os.system("chmod +x %s/submit_unc.sh" % outDir)
os.system("chmod +x %s/submit_antiIso.sh" % outDir)

print "Prepared to submit cmsRun jobs."
print "Submit jobs using helper script under pass1 directory yourself."
print "> cd pass1"
print "> ./submit.sh"
print ""
print "or, use the xargs magic to submit them all"
print "> cat submit*.sh | sed -e 's;create-batch;;g' -e 's; *$;;g' | xargs -l -P20 create-batch"
print ""

