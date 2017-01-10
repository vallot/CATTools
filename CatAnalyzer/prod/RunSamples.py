import os, time, socket

InputDB        = "SamTra.info"
FileHeader     = "Tree_LepJets_KFCSVOrder01NoSkim_v8-0-1_Spring16-80X_15920pb-1"
OutputLocation = "/xrootd/store/user/brochero/v8-0-1/"

DelayTime = 120. # Time in seconds
maxNjobs = 4000  # Maximum number of jobs running simultaneously
def NumberOfCondorJobs (str):
    condorNF = ".tempCondor_" + socket.gethostname() + "_" + str + "_" + time.strftime('%Hh%Mm%Ss') + ".info"
    print "condor_q brochero > " + condorNF
    os.system("condor_q brochero > " + condorNF)
    with open(condorNF, "rb") as fcondor:
        fcondor.seek(-2, 2)             # Jump to the second last byte.
        while fcondor.read(1) != b"\n": # Until EOL is found...
            fcondor.seek(-2, 1)         # ...jump back the read byte plus one more.
        tempjr = []
        tempjr = fcondor.readline().split()
        condorJobsRun = int(tempjr[0])
    print "rm -rf " + condorNF
    os.system("rm -rf " + condorNF)
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
SamArg = []
nfilesperjob = 0
nJobs = 0
isample = 0
nsrunning = 0 
groupJobs = False

for line in fr:
    if ".txt" in line:
        tempsn = []
        tempsn = line.rstrip().split()
        if len(tempsn) > 1:
            SamNam.append(str(tempsn[0]))
            SamLoc.append(str(tempsn[1]))
            tSamArg = ""
            if len(tempsn) > 2:
                for narg in range(2, len(tempsn)):
                    tSamArg += " " + str(tempsn[narg])
            SamArg.append(str(tSamArg))
        else:
            print "There is no enough parameter!!!"
            quit()
        del tempsn
        nfSamLoc = NumberOfFiles (SamLoc[nsrunning])
        if nfSamLoc < 1000: maxf = 1
        else:               maxf =  int(round(nfSamLoc/1000.))    

        print  str(nfSamLoc) + " root files. Max number of files per job " + str(maxf)  
        CreateJob = str("./create-batch --jobName " + SamNam[nsrunning] + " --fileList " + SamLoc[nsrunning] + " --maxFiles " + str(maxf) + " --cfg ttbbLepJetsAnalyzer_cfg.py --queue batch6 --transferDest /xrootd/store/user/brochero/")
        if SamArg[nsrunning] is not "":
            CreateJob += " --args ' " + SamArg[nsrunning] + " ' "
        print CreateJob
        os.system(CreateJob)

        nJobs += nfSamLoc/maxf
        print "Number of jobs " + str(nJobs)
        if nJobs > maxNjobs:
            groupJobs = True

        jobsRunning = True
        while jobsRunning:
            condorJobsRun = NumberOfCondorJobs (SamNam[nsrunning])    
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
                del SamArg[:]                
        isample += 1
fr.close()
print "All samples done. Output files are in " + OutputLocation
