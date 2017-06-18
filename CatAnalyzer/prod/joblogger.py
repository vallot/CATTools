#!/usr/bin/env python
import os
import re
import sys
import time
import commands

from os.path import join, getsize

def sorted_ls(path):
    mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
    return list(sorted(os.listdir(path), key=mtime))

def doWork(finalpath, action):

  sample = os.path.basename(finalpath)

  if action == "report" or action == "submit":
    joblist = os.listdir(finalpath+"/Log")
    for i in joblist:
      num = i.split("_")
      file = finalpath+"/Res/Tree_"+sample+"_"+num[1]+".root"
      filesize = 0
      if os.path.isfile( file ):
        filesize = os.path.getsize( file )
      #another test
      error = ""
      jobDir = sorted_ls( finalpath+"/Log/Job_"+num[1] )
      logDir = []
      for jobDirList in jobDir:
        if jobDirList.startswith("LSF"):
          logDir.append(jobDirList)  
      if len(logDir) > 0:
        log = finalpath+"/Log/Job_"+num[1]+"/"+logDir[-1]+"/STDOUT"
        for line in open(log):
          if "Segmentation fault" in line:
            error = "Segmentation fault"
          if "Probable system resources exhausted" in line:
            error = "Probable system resources exhausted"
          if "exhausted the virtual memory available to the process" in line:
            error = "exhausted the virtual memory available to the process"
      else:
        error = "No log file available"
      ## print out
      if filesize < 10000 or error != "":
        if action == "report":
          if filesize == 0:
            print "Bad job number : " + num[1] + " : file was not created"  
          else: 
            if error == "":
              print "Bad job number : " + num[1] + " : file size is too small. it must be corrupted." 
            else:
              print "Bad job number : [" +logDir[-1]+"]" + num[1] + " : " + error
        if action == "submit":
          print "Submitting job : " + num[1] 
          currdir = commands.getoutput('pwd')
          os.chdir("%s/Log/Job_%s" % (finalpath,num[1]) )
          os.system("bsub -q 8nm < ./batchScript.sh")
          os.chdir(currdir)

  if action == "merge":
    pathdir = os.path.dirname(finalpath)
    print "Merging files in " + finalpath
    destinationfile = pathdir+"/Tree_"+sample+".root"
    if os.path.isfile( destinationfile ):
      print destinationfile + " exists."
    else:
      #filelist = os.listdir(finalpath+"/Res")
      filelist = os.listdir(finalpath)
      filelist.sort()
      nlist =  len(filelist)

      if nlist > 500 :
        #make list of 500 files
        n = 1 
        list = [] 
        tmp = ""
        for f in filelist:
          #tmp += finalpath+"/Res/"+f+" "
          tmp += finalpath+"/"+f+" "
          k = n%500
          if k == 0 or n == nlist:
            list.append(tmp)
            tmp = ""
          n += 1 
        ## merge 500 files first to tmp file.
        j = 0 
        nmerge = len(list)
        print "Will creat "+str(nmerge)+" tmpfile.root files"
        for l in list:
          tmpfile = pathdir+"/tmp_"+sample+"_"+str(j)+".root"
          os.system("hadd -f "+tmpfile+" "+l)
          j += 1

        os.system("hadd -f "+destinationfile+" "+pathdir+"/tmp_"+sample+"*.root")
        os.system("rm -rf "+pathdir+"/tmp_"+sample+"*.root")
      else:
        #os.system("hadd -f "+destinationfile+" "+finalpath+"/Res/Tree_*.root")
        os.system("hadd -f "+destinationfile+" "+finalpath+"/Tree_*.root")
    #x = raw_input("Remove directory %s (y/n)?" % (finalpath))
    #if x == "y":
    #  os.system("rm -rf "+finalpath+"/Log")
    #  os.system("rm -rf "+finalpath+"/Res")
    #  os.system("rm -rf "+finalpath)

def jobReport( path , action):
   list = os.listdir(path)
   for sample in list:
     finalpath = path + "/" + sample
     if os.path.isdir(finalpath):
       print "Checking " + finalpath
       #logdir = os.path.exists(finalpath+"/Log")
       logdir = os.path.exists(finalpath)
       if logdir == True:
         doWork(finalpath, action)
       else:
         print "No sample Log directory!"
  

if __name__ == '__main__':

    from optparse import OptionParser

    parser = OptionParser()
    parser.usage = """
    %prog [options] <sample_dir>
    
    Prints the list of bad jobs.
    Using the options, you can get a log of what happened during each bad job,
    and you can resubmit these jobs.
    """

    parser.add_option("-r", "--report", dest="report",
                      action = 'store_true',
                      default=False,
                      help='Print report for bad jobs.')

    parser.add_option("-s", "--submit", dest="submit",
                      action = 'store_true',
                      default=False,
                      help='Resubmission command')

    parser.add_option("-m", "--merge", dest="merge",
                      action = 'store_true',
                      default=False,
                      help='hadd -f merged.root *.root')


    (options,args) = parser.parse_args()
 
    path = args[0]

    if options.report:
        jobReport(path, "report")
    if options.submit:
        jobReport(path, "submit")
    if options.merge:
        jobReport(path, "merge")

