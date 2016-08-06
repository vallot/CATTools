import os, time

InputDB        = "SamplesTest.info"
FileHeader     = "Tree_LepJets_KinReco_v8-0-0_Spring16-80X_5913pb-1"
OutputLocation = "/xrootd/store/user/brochero/v8-0-0/"

DelayTime = 120. # Time in seconds

def NumberOfCondorJobs ():
    print "condor_q brochero > .tempCondor.info" 
    os.system("condor_q brochero > .tempCondor.info") 
    with open(".tempCondor.info", "rb") as fcondor:
        fcondor.seek(-2, 2)             # Jump to the second last byte.
        while fcondor.read(1) != b"\n": # Until EOL is found...
            fcondor.seek(-2, 1)         # ...jump back the read byte plus one more.
        tempjr = []
        tempjr = fcondor.readline().split() 
        condorJobsRun = int(tempjr[0])
    print "rm -rf .tempCondor.info" 
    os.system("rm -rf .tempCondor.info") 
    return condorJobsRun;

def NumberOfFiles (str):
    fSamLoc = open(str,'r')
    nfSamLoc = 0
    for slline in fSamLoc:
        if ".root" in slline:
            nfSamLoc += 1
    return nfSamLoc;


print "Reading database for Lepton+Jets in " + InputDB

fr = open(InputDB,'r')

NSamples = 0
for line in fr:
    if ".txt" in line:
        NSamples += 1

print "Number of sample to be processed -> " + str(NSamples)
fr.seek(0)

SamNam = []
SamLoc = []
nfilesperjob = 0
nJobs = 0
isample = 0
nsrunning = 0 
groupJobs = False

for line in fr:
    if ".txt" in line:
        tempsn = []
        tempsn = line.rstrip().split()
        SamNam.append(str(tempsn[0]))
        SamLoc.append(str(tempsn[1]))
        del tempsn
        nfSamLoc = NumberOfFiles (SamLoc[nsrunning])
        if nfSamLoc < 1000: maxf = 1
        else:               maxf =  int(round(nfSamLoc/1000.))    

        print  str(nfSamLoc) + " root files. Max number of files per job " + str(maxf)  
        print "./create-batch --jobName " + SamNam[nsrunning] + " --fileList " + SamLoc[nsrunning] + " --maxFiles " + str(maxf) + " --cfg ttbbLepJetsAnalyzer_cfg.py --queue batch6 --transferDest /xrootd/store/user/brochero/"
        os.system("./create-batch --jobName " + SamNam[nsrunning] + " --fileList " + SamLoc[nsrunning] + " --maxFiles " + str(maxf) + " --cfg ttbbLepJetsAnalyzer_cfg.py --queue batch6 --transferDest /xrootd/store/user/brochero/")

        nJobs += nfSamLoc/maxf
        print "Number of jobs " + str(nJobs)
        if nJobs > 1000:
            groupJobs = True

        jobsRunning = True
        while jobsRunning:
            condorJobsRun = NumberOfCondorJobs ()    
            if groupJobs:
                print "Limit of jobs reached..."
                print "Number of Condor Jobs = " + str(condorJobsRun)
                if condorJobsRun is 0:
                    groupJobs = False
                    nJobs = 0
                    print "Condor jobs are done."
                else:
                    print "Condor is still running " + str(condorJobsRun) + " job(s). Waiting " + str(round(DelayTime/60.)) + " min to check again..."
                    time.sleep(DelayTime) # Time in seconds 
            elif nJobs > 0:
                if isample < (NSamples-1):
                    jobsRunning = False
                    nsrunning += 1
                    print "More Jobs can be sent..."
                else:
                    print "Waiting " + str(round(DelayTime/60.)) + " min for the last " + str(condorJobsRun) + " job(s)..."
                    nJobs = condorJobsRun
                    time.sleep(DelayTime)
            elif nJobs == 0:
                print "All jobs done for " + str(len(SamNam)) + " samples."
                print str(len(SamNam))
                for index in range(len(SamNam)):
                    print "Merging " + SamNam[index] + " sample."
                    print "hadd " + OutputLocation + FileHeader + "_" + SamNam[index] + ".root " + SamNam[index] + "/*.root "
                    os.system("hadd " + OutputLocation + FileHeader + "_" + SamNam[index] + ".root " + SamNam[index] + "/*.root ")
                    print "Removing " + SamNam[index] + " directory..."
                    print "rm -rf " + SamNam[index]
                    os.system("rm -rf " + SamNam[index])
                nsrunning = 0
                jobsRunning = False
                del SamNam[:]                
                del SamLoc[:]                
    isample += 1
fr.close()
print "All samples done. Output files are in " + OutputLocation
