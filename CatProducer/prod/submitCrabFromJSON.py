#!/usr/bin/env python
import sys, os

if len(sys.argv) < 2:
    print "Usage: %s vX-Y-Z"
    print "or, you can put your own 'catGetDatasetInfo'-compatible json file"
    sys.exit(1)

reqName = ""
if sys.argv[1].endswith('.json'):
    if len(sys.argv) < 3:
        print "!!! Please give me the version name of the production, vX-Y-Z"
        sys.exit(1)
    js = sys.argv[1]
    reqName = sys.argv[2]
else:
    reqName = sys.argv[1]
    os.system("catGetDatasetInfo %s" % reqName)
    js = "%s/src/CATTools/CatAnalyzer/data/dataset/dataset.json"%os.environ['CMSSW_BASE']
if not os.path.exists(js):
    print "!!! Cannot open json file,", js
    sys.exit(2)

import CATTools.CatProducer.catDefinitions_cfi as cat
import json
datasets = json.load(open(js))
queuesRD, queuesMC = [], []
for d in datasets:
    name = d['DataSetName']
    if len(name) == 0: continue
    pdName, sdName = name.split('/')[1], name.split('/')[2]
    if os.path.exists('crab_%s_%s' % (reqName, pdName)): continue
    elif os.path.exists('crab_%s_%s_%s' % (reqName, pdName, sdName)): continue
    if len(d['path']) != 0: continue

    runOnMC = (d['type'] != 'Data')

    ## Options for this sample
    opts = {"useMiniAOD":"True"}
    if runOnMC:
        opts['globalTag'] = cat.globalTag_mc
        opts['runOnMC'] = "True"
        if name.startswith('/TT') or name.startswith("/tt"): opts['runGenTop'] = "True"
    else:
        opts['globalTag'] = cat.globalTag_rd
        opts['runOnMC'] = "False"

    ## Override options if it is set from the dataset JSON file
    if 'options' in d:
        for opt in d['options'].split(','):
            key, val = opt.split('=')
            opts[key] = val

    ## List-fy opts to use in crab API
    opts = ["%s=%s" % (key, opts[key]) for key in opts]

    if d['type'] == 'Data': queuesRD.append((name, opts))
    else: queuesMC.append((name, opts))

## Now job configuration is almost ready. Use CrabAPI to configure jobs
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.transferLogs    = False
config.General.transferOutputs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'PAT2CAT_cfg.py'

config.section_("Data")
config.Data.publication  = False
config.Data.allowNonValidInputDataset = True

config.section_("Site")
config.Site.storageSite = 'T3_KR_KISTI'
config.Data.outLFNDirBase = '/store/group/CAT/' 

## Submit real data jobs
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 40
config.Data.lumiMask = os.environ['CMSSW_BASE']+'/src/CATTools/CatProducer/data/LumiMask/'+cat.lumiJSON
for name, opts in queuesRD:
    label = name.split('/')[1]+'_'+name.split('/')[2]
    config.General.requestName = '%s_%s' % (reqName, label)
    config.Data.inputDataset = name
    config.Data.outputDatasetTag = '%s_%s' % (reqName, name.split('/')[2])
    print config
    
## Submit MC jobs
config.Data.splitting='FileBased'
config.Data.unitsPerJob=1 
config.Data.lumiMask = ''
for name, opts in queuesMC:
    label = name.split('/')[1]
    config.General.requestName = '%s_%s' % (reqName, label)
    config.Data.inputDataset = name
    config.Data.outputDatasetTag = '%s_%s' % (reqName, name.split('/')[2])
    print config

