#!/usr/bin/env python
import sys, os

if len(sys.argv) < 2 or '-h' in sys.argv or '--help' in sys.argv:
    print "Usage: %s vX-Y-Z" % sys.argv[0]
    print "or, you can put your own 'catGetDatasetInfo'-compatible json file"
    print "       %s dataset.json vX-Y-Z" % sys.argv[0]
    sys.exit(0)

doSubmit = True 
if '--dryrun' in sys.argv:
    doSubmit = False
    sys.argv.remove('--dryrun')

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
    name = d['name']
    dataset = d['DataSetName']
    if len(dataset) == 0: continue
    pdName, sdName = dataset.split('/')[1], dataset.split('/')[2]
    if os.path.exists('crab_%s_%s' % (reqName, name)): continue
    if len(d['path']) != 0: continue

    runOnMC = (d['type'] != 'Data')

    ## Options for this sample
    opts = {"useMiniAOD":"True"}
    if runOnMC:
        opts['globalTag'] = cat.globalTag_mc
        opts['runOnMC'] = "True"
    else:
        opts['globalTag'] = cat.globalTag_rd
        opts['runOnMC'] = "False"

    ## Override options if it is set from the dataset JSON file
    if 'options' in d:
        for opt in d['options'].split(','):
            key, val = opt.split('=')
            opts[key] = val

    ## List-fy opts to use in crab API
    opts = [str("%s=%s" % (key, opts[key])) for key in opts]

    if d['type'] == 'Data': queuesRD.append((name, dataset, opts))
    else: queuesMC.append((name, dataset, opts))

## Now job configuration is almost ready. Use CrabAPI to configure jobs
from WMCore.Configuration import Configuration
def initConfig(runOnMC, reqName):
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
    config.Data.outLFNDirBase = '/store/group/CAT/%s' % reqName

    config.Data.splitting='FileBased'
    config.Data.unitsPerJob=1 

    return config

## Start job submission
#from CRABAPI.RawCommand import crabCommand

## Submit real data jobs
import time
from tempfile import mkstemp
for name, dataset, opts in queuesRD:
    print "@@@ Creating", name
    config = initConfig(False, reqName)
    config.General.requestName = '%s_%s' % (reqName, name)
    config.Data.inputDataset = dataset
    config.Data.outputDatasetTag = dataset.split('/')[2]
    config.JobType.pyCfgParams = opts
    fd, fName = mkstemp(suffix='.py')
    f = os.fdopen(fd, 'w')
    print>>f, config
    f.close()
    print "@@@ Submitting", name
    #crabCommand('submit', config=fName, dryrun=(not doSubmit))
    if doSubmit: os.system('crab submit %s' % fName)
    else: os.system('crab submit --dryrun %s' % fName)
    time.sleep(1)
    os.remove(fName)
    if os.path.exists(fName+'c'): os.remove(fName+'c') ## remove .pyc file

## Submit MC jobs
for name, dataset, opts in queuesMC:
    print "@@@ Creating", name
    config = initConfig(True, reqName)
    config.General.requestName = '%s_%s' % (reqName, name)
    config.Data.inputDataset = dataset
    config.Data.outputDatasetTag = dataset.split('/')[2]
    config.JobType.pyCfgParams = opts
    fd, fName = mkstemp(suffix='.py')
    f = os.fdopen(fd, 'w')
    print>>f, config
    f.close()
    print "@@@ Submitting", name
    #crabCommand('submit', config=fName, dryrun=(not doSubmit))
    if doSubmit: os.system('crab submit %s' % fName)
    else: os.system('crab submit --dryrun %s' % fName)
    time.sleep(1)
    os.remove(fName)
    if os.path.exists(fName+'c'): os.remove(fName+'c') ## remove .pyc file
