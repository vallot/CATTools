#!/usr/bin/env python
import os, time, socket, sys

UserName       = os.environ["USER"]
BaseDir        = os.environ["CMSSW_BASE"]+"/src"
InputDB        = str(sys.argv[1])
FileHeader     = "Tree_LepJets_Cat_v8-0-6_Spring16-80X_36814pb-1"
OutputLocation = "/xrootd/store/user/brochero/v8-0-6/"

DelayTime = 120. # Time in seconds
maxNjobs = 2000  # Maximum number of jobs running simultaneously

def NumberOfCondorJobs (str):
    condorNF = ".tempCondor_" + socket.gethostname() + "_" + str + "_" + time.strftime('%Hh%Mm%Ss') + ".info"
    print "condor_q %s > %s" % (UserName, condorNF)
    os.system("condor_q %s > %s" % (UserName, condorNF))
    if not os.path.isfile(condorNF):
        condorJobsRun = 1000 # If it cannot read it returs 1000 jobs!
    else:
        with open(condorNF, "rb") as fcondor:
            fcondor.seek(-2, 2)             # Jump to the second last byte.
            while fcondor.read(1) != b"\n": # Until EOL is found...
                fcondor.seek(-2, 1)         # ...jump back the read byte plus one more.
            tempjr = []
            tempjr = fcondor.readline().split()
            if tempjr[1] == "jobs;":
                condorJobsRun = int(tempjr[0])
            else: condorJobsRun = 1000 # If it cannot read it returs 1000 jobs!
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
            SamLoc.append(BaseDir+"/CATTools/CatAnalyzer/data/dataset_v8-0-6/"+str(tempsn[1]))
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
        if nfSamLoc < 500: maxf = 1
        else:              maxf =  int(round(nfSamLoc/500.))    

        print  str(nfSamLoc) + " root files. Max number of files per job " + str(maxf)  
        CreateJob = str("create-batch --jobName " + SamNam[nsrunning] + " --fileList " + SamLoc[nsrunning] + " --maxFiles " + str(maxf) + " --cfg ttbbLepJetsAnalyzer_cfg.py --queue batch6")
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
                    print "hadd -f " + OutputLocation + FileHeader + "_" + SamNam[index] + ".root " + SamNam[index] + "/*.root "
                    os.system("hadd -f " + OutputLocation + FileHeader + "_" + SamNam[index] + ".root " + SamNam[index] + "/*.root ")
                    BackupDir = str(OutputLocation + "BackupFiles")
                    if not os.path.isdir(BackupDir):
                        print "Creating " + BackupDir
                        os.system("mkdir " + BackupDir)
                    print "Moving " + SamNam[index] + " directory to " + BackupDir
                    os.system("mv " + SamNam[index] + "  " + BackupDir)
                    print "Files kept in the BackupFiles directory!! Delete them your self! "
                nsrunning = 0
                jobsRunning = False
                del SamNam[:]                
                del SamLoc[:]                
                del SamArg[:]                
        isample += 1
fr.close()
print "All samples done. Output files are in " + OutputLocation
