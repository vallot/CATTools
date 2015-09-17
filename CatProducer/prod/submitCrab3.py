#!/usr/bin/env python
import os,sys,getopt,shutil

requestName = ""
datasets = []
inputFile =""
submit = False
lumiMask =""
globalTag =""
crabcommand ='crab submit -c crab.py'

try:
    opts, args = getopt.getopt(sys.argv[1:],"hsi:n:l:g:",["requestName","inputFile","lumiMask","globalTag"])
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
    elif opt in ("-g", "--globalTag"):
        globalTag = arg

if len(globalTag) == 0:
    print "need to define globalTag, -g <globalTag>"
    sys.exit()

print datasets
for dataset in datasets:
    if len(dataset) < 10:
        continue
    if dataset.startswith("#"):
        continue
    
    isMiniAOD="False"
    datatype = dataset.strip().split("/")[-1]
    if datatype == "MINIAOD" or datatype == "MINIAODSIM" :
        isMiniAOD="True"

    isMC = True
    dataSplitting   = " Data.splitting='FileBased' "
    dataUnitsPerJob = " Data.unitsPerJob=1 "
    dataLumiMask    = ""
    ### dirty way for now since crab3 doesnt allow lists to be passed by cmd line
    pyCfgParams     = "config.JobType.pyCfgParams = ['runOnMC=True','useMiniAOD=%s','globalTag=%s']"%(isMiniAOD,globalTag)
    ### MC or Data?
    if datatype == "AOD" or datatype == "MINIAOD" :
        if len(lumiMask) == 0:
            print "NOTE! no lumiMask was selected!"
            print "NOTE! no lumiMask was selected!"
            print "NOTE! no lumiMask was selected!"
            print "NOTE! no lumiMask was selected!"
            print "NOTE! no lumiMask was selected!"
            print "NOTE! no lumiMask was selected!"
            print "NOTE! no lumiMask was selected!"
            print "NOTE! no lumiMask was selected!"
            print "NOTE! no lumiMask was selected!"
            print "NOTE! no lumiMask was selected!"
        isMC = False
        #dataSplitting   = " Data.splitting='LumiBased' "
        #dataUnitsPerJob = " Data.unitsPerJob=10 "
        dataLumiMask    = " Data.lumiMask='%s'"%(lumiMask)
        pyCfgParams     = "config.JobType.pyCfgParams = ['runOnMC=False','useMiniAOD=%s','globalTag=%s']"%(isMiniAOD,globalTag)
    ## Special option for TTbar signal samples
    if isMC and (dataset.startswith('/TT') or dataset.startswith('/tt')):
        pyCfgParams = pyCfgParams[:-2] + (",'runGenTop=True']")

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
        publishDataName = '%s_%s'%(requestName,dataset.split("/")[2])
    
    sendjob = crabcommand + " Data.publishDataName='%s' General.requestName='%s' Data.inputDataset='%s'"%(publishDataName,dataRequestName,dataset) + dataSplitting + dataUnitsPerJob + dataLumiMask
    print sendjob
    if submit:
        print "submiting job"
        os.system(sendjob)

    #lines = open("crab.py")
    #print lines.read()
    os.remove("crab.py")
    
if not submit:
    print "Dry run, not submitting job and only printing crab3 command"
    print "Add -s to submit job"
    print 'Usage : ./submitCrab3.py -n <requestName> -i <inputFile> -s'
