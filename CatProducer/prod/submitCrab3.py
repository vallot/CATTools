#!/usr/bin/env python
import os,json,sys,shutil,time,getopt
import CATTools.CatProducer.catDefinitions_cfi as cat

def submitjob(requestName, dataset, globalTag, lumiMask, submit):
    print "creating job"
    print dataset

    isMiniAOD="False"
    isMC = True
    datatype = dataset.strip().split("/")[-1]
    if datatype == "MINIAOD" or datatype == "MINIAODSIM" :
        isMiniAOD="True"

    if datatype == "AOD" or datatype == "MINIAOD" :
        isMC = False        
        if globalTag == None:
            globalTag = cat.globalTag_rd

    if globalTag == None:
        globalTag = cat.globalTag_mc
    if lumiMask == None:
        lumiMask = '../data/LumiMask/%s.txt'%cat.lumiJSONSilver
        
    dataSplitting   = " Data.splitting='FileBased' "
    dataUnitsPerJob = " Data.unitsPerJob=1 "
    dataLumiMask    = ''    
    ## Special option for TTbar signal samples
    doGenTop = False
    if (dataset.startswith('/TT') or dataset.startswith('/tt')):
        doGenTop = True
    ### dirty way for now since crab3 doesnt allow lists to be passed by cmd line
    pyCfgParams     = "config.JobType.pyCfgParams = ['runOnMC=True','useMiniAOD=%s','globalTag=%s','runGenTop=%s']"%(isMiniAOD,globalTag,doGenTop)
    ### MC or Data?
    if isMC == False:
        dataSplitting   = " Data.splitting='LumiBased' "
        dataUnitsPerJob = " Data.unitsPerJob=40 "
        dataLumiMask    = " Data.lumiMask='%s'"%(lumiMask)
        pyCfgParams     = "config.JobType.pyCfgParams = ['runOnMC=False','useMiniAOD=%s','globalTag=%s']"%(isMiniAOD,globalTag)

    ## pyCfgParams cannot be set from cmd line yet
    shutil.copy2('crabConfig.py', 'crab.py')
    print pyCfgParams
    with open("crab.py", "a") as myfile:
        myfile.write(pyCfgParams)

    ### adjusting label for requestName
    dataset = dataset.strip()
    if isMC :
        label = dataset.split("/")[1]
    else :
        label = dataset.split("/")[1]+"_"+dataset.split("/")[2]

    if requestName:
        dataRequestName = '%s_%s'%(requestName,label)
        outputDatasetTag = '%s_%s'%(requestName,dataset.split("/")[2])

    if submit:
        sendjob = "crab submit -c crab.py Data.outputDatasetTag='%s' General.requestName='%s' Data.inputDataset='%s'"%(outputDatasetTag,dataRequestName,dataset) + dataSplitting + dataUnitsPerJob + dataLumiMask
    else :
        sendjob = "crab submit --dryrun -c crab.py Data.outputDatasetTag='%s' General.requestName='%s' Data.inputDataset='%s'"%(outputDatasetTag,dataRequestName,dataset) + dataSplitting + dataUnitsPerJob + dataLumiMask

    print sendjob
    print "submiting job"
    os.system(sendjob)
    #lines = open("crab.py")
    #print lines.read()
    os.remove("crab.py")
    time.sleep(5)
    
submitBlock = None
requestName = ""
datasets = []
inputFile =None
submit = False
lumiMask =""
globalTag =""

try:
    opts, args = getopt.getopt(sys.argv[1:],"hsi:n:l:g:b:",["requestName","inputFile","lumiMask","globalTag","submitBlock"])
except getopt.GetoptError:          
    print 'Usage : ./submitCrab3.py -n <requestName> -i <inputFile> -l <lumiMask> -g <globalTag>'
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./submitCrab3.py -n <requestName> -i <inputFile> -l <lumiMask> -g <globalTag>'
        sys.exit()
    elif opt in ("-n", "--requestName"):
        requestName = arg
    elif opt in ("-s"):
        submit = True
    elif opt in ("-i", "--inputFile"):
        inputFile = arg
        if os.path.isfile(inputFile):
            lines = open(inputFile)
            datasets = lines.readlines()
        else:
            datasets.append(inputFile)
    elif opt in ("-l", "--lumiMask"):
        lumiMask = arg
    elif opt in ("-b", "--submitBlock"):
        submitBlock = arg
    elif opt in ("-g", "--globalTag"):
        globalTag = arg

if requestName == "" :
    print "requestName(-n) is mandantory"
    sys.exit(-1)

if inputFile is None:
    catGetDatasetInfo = 'catGetDatasetInfo %s'%(requestName)
    os.system(catGetDatasetInfo)
    datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/dataset.json"%os.environ['CMSSW_BASE']))
    for d in datasets:
        dataset = d['DataSetName']
        if len( d['path']) == 0:
            #print d['path'], len( d['path'])
            submitjob(requestName, dataset, None,None, submit)
        
        #if submitBlock == '1' and 'QCD' in dataset:
        #    continue
        #if submitBlock == '2' and 'QCD' not in dataset:
        #    continue

else:
    for dataset in datasets:
        if len(dataset) < 10:
            continue
        if dataset.startswith("#"):
            continue
        submitjob(requestName, dataset, globalTag, lumiMask, submit)
        

if not submit:
    print "Dry run, not submitting job and only printing crab3 command"
    print "Add -s to submit job"
    print 'Usage : ./submitCrab3.py -n <requestName> -i <inputFile> -s'
