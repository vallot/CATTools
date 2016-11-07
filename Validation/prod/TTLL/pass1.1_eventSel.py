#!/usr/bin/env python

## Define recipes
recipes = {
    'nominal':{'nominal':[]},
    'syst':{
        'mu_pt/up':['eventsTTLL.muon.scaleDirection=1'],
        'mu_pt/dn':['eventsTTLL.muon.scaleDirection=-1'],
        'el_pt/up':['eventsTTLL.electron.scaleDirection=1'],
        'el_pt/dn':['eventsTTLL.electron.scaleDirection=-1'],
        'jet_cor/up':['eventsTTLL.jet.scaleDirection=1'],
        'jet_cor/dn':['eventsTTLL.jet.scaleDirection=-1'],
    },
    'systMC':{
        'jet_res/up':['eventsTTLL.jet.resolDirection=1'],
        'jet_res/dn':['eventsTTLL.jet.resolDirection=-1'],
        'pileup/up':['eventsTTLL.vertex.pileupWeight="pileupWeight:up"'],
        'pileup/dn':['eventsTTLL.vertex.pileupWeight="pileupWeight:dn"'],
        'mu_eff/up':['eventsTTLL.muon.efficiencySFDirection=1'],
        'mu_eff/dn':['eventsTTLL.muon.efficiencySFDirection=-1'],
        'el_eff/up':['eventsTTLL.electron.efficiencySFDirection=1'],
        'el_eff/dn':['eventsTTLL.electron.efficiencySFDirection=-1'],
    },
    'scaleup':{},
    'scaledn':{},
    'pdf':{},
}

## Scale up/down systematic uncertainty from LHE weight
## This uncertainty have to be combined with envelope
## Let us assume index1-10 are for the scale variations (muF & muR)
## total 8 scale variations, 3 muF x 3 muR and one for nominal weight
## Skip unphysical scale variation combinations, (muF=2, muR=0.5) and (muF=0.5, muR=2) should be skipped
## and combine with PS level scale variations samples
for i in range(3):
    for ss in ("up", "dn"):
        if ss == 'up':
            s  = ['eventsTTLL.genWeight.src="flatGenWeights:scaleup"',
                  'agen.weight="flatGenWeights:scaleup"']
        else: s = ['eventsTTLL.genWeight.src="flatGenWeights:scaledown"',
                   'agen.weight="flatGenWeights:scaledown"']

        s.extend(['eventsTTLL.genWeight.index=%d' % i, 'agen.weightIndex=%d' % i])
        recipes['scale'+ss]['gen_scale/%s_%d' % (ss, i)] = s[:]

## PDF weights
## Weight vector size differs to include different PDF considerations
## -> 110 variations for aMC@NLO, 248 for POWHEG
## Basically, (1+8 scale variations) + (1+100 NNPDF variations) + (other PDF variations) + (1+N hdamp variations)
## NOTE: there is weight vector, but we don't do it for LO generator here.
for i in range(100):
    r = ['eventsTTLL.genWeight.src="flatGenWeights:pdf"', 'eventsTTLL.genWeight.index=%d' % i,
         'agen.weight="flatGenWeights:pdf"', 'agen.weightIndex=%d' % i]
    recipes['pdf']['gen_pdf/%d' % i] = r[:] 

## Define list of samples and their recipes
samples = {
    'data nominal syst':[
        ## Real data
        ## cmsRun analyze_dat_cfg.py
        ## do common systematic variations
        "DoubleEG_Run2016B", "DoubleEG_Run2016C", 
        "DoubleMuon_Run2016B", "DoubleMuon_Run2016C", 
        "MuonEG_Run2016B", "MuonEG_Run2016C", 
    ],
    'bkg nominal syst systMC':[
        ## MC samples without ttbar genFilter
        ## cmsRun analyze_bkg_cfg.py
        ## do common systematic variations
        ## do MC common systematic variations
        ## No ttbar specific gen level analyzer, no ttbar genFilters
        "DYJets", "DYJets_10to50", "DYJets_MG", "DYJets_MG_10to50", 
        "SingleTop_s", "SingleTop_t", "SingleTbar_t", "SingleTop_tW", "SingleTbar_tW", 
        "WJets_MG", 
        "WW", "WZ", "ZZ", 
        "WWW", "WWZ", "WZZ", "ZZZ", 
    ],
    'sig nominal syst systMC scaleup scaledn pdf':[
        ## MC Signal samples, select ttbar-dilepton at Gen.Level
        ## cmsRun analyze_sig_cfg.py
        ## do common systematic variations
        ## do MC common systematic variations
        ## do ttbar specific gen level analyzer, do ttbar-dilepton genFilters
        ## do PDF variations
        "TTLL_powheg", 
        "ttbb",
    ],
    'sig nominal syst systMC':[
        "ttW", "ttZ",
    ],
    'sig.noll nominal syst systMC':[
        ## ttbar-others
        ## cmsRun analyze_sig_cfg.py
        ## do common systematic variations
        ## do MC common systematic variations
        ## do ttbar specific gen level analyzer, do ttbar-dilepton genFilters
        "TT_powheg", 
        "ttW", "ttZ",
    ],
    'sig nominal scaleup':[
        ## MC signal samples, for the systematic unc. variations
        ## Therefore, produce nominals only
        ## cmsRun analyze_sig_cfg.py
        ## do ttbar specific gen level analyzer, do ttbar-dilepton genFilters
        "TTLL_powheg_scaleup", "TT_powheg_scaleup",
    ],
    'sig nominal scaledn':[
        "TTLL_powheg_scaledown", "TT_powheg_scaledown",
    ],
    'sig nominal':[
        ## Variation on generator choice
        ## produce nominals only
        ## cmsRun analyze_sig_cfg.py
        ## do ttbar specific gen level analyzer, do ttbar-dilepton genFilters
        "TTJets_MG", "TTJets_aMC", 
        "TT_powheg_herwig", 

        "TT_powheg_mtop1695", "TT_powheg_mtop1755",

        "TT_powheg_noCR", "TT_powheg_mpiOFF", 

        "TTLL_powheg_alphaS", "TT_powheg_alphaS", 
        "TT_powheg_herwig_mpiOFF", 
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
        baseargs += ["filterPartonTTLL.invert=True"]
        suffix = "_Others"

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
for cmd in sorted([x for x in cmds if 'nominal' in x]): print>>out, cmds[cmd]
out.close()
out = open("%s/submit_unc.sh" % outDir, "w")
for cmd in sorted([x for x in cmds if 'nominal' not in x]): print>>out, cmds[cmd]
out.close()

os.system("chmod +x %s/submit_nominal.sh" % outDir)
os.system("chmod +x %s/submit_unc.sh" % outDir)

print "Prepared to submit cmsRun jobs."
print "Submit jobs using helper script under pass1 directory yourself."
print "> cd pass1"
print "> ./submit.sh"
print ""
print "or, use the xargs magic to submit them all"
print "> cat submit.sh | sed -e 's;create-batch;;g' | xargs -L1 -P20 create-batch"
print ""

