#!/usr/bin/env python

import ROOT,os,sys,copy,json
ROOT.gROOT.SetBatch(True)

execmd = "h2muDraw.py"
cmssw = os.environ['CMSSW_BASE']

outDir = "plot"
if not os.path.isdir(outDir):os.mkdir(outDir)
try:
    info_json = json.load(open("%s/info.json" % outDir))
except:
    print"""
    !! WARNNING

    please try to run "set_config.py" file which makes json file at first.
    and then you can run this file.
    """
try:json_used = sys.argv[1]
except:json_used = 'Golden'
if json_used == 'Silver':
    print"""
    Used Silver JSON
    For ggF
    """
else:print"""
    Used Golden JSON
    For VBF
    if you want to use Silver JSON, argument is 'Silver'.
    """
sh_s = """#!/bin/bash
workD=%s
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $workD
eval `scramv1 runtime -sh`
cd -
mkdir -p %s
host=`hostname`
./%s
"""

jds_s = """
executable = %s
universe = vanilla
getenv = True

num= %s
dirSaved= %s
f_name= %s

error = $(dirSaved)/error_$(num)
log = $(dirSaved)/log_$(num)

should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = h2muDraw.py
transfer_output_files = $(f_name).png,
transfer_output_remaps = "$(f_name).png = $(dirSaved)/$(f_name).png"

queue
""" 
for i,info in enumerate(info_json):
    plotvar,f_name,cut,binning,x_name,y_name = info['plotvar'],info['f_name'],info['cut'],info['binning'],info['x_name'],info['y_name']
    if plotvar == 'll_m':
        args = "-c \'%s\' -b %s -p \'%s\' -x \'%s\' -y \'%s\' -f \'%s\' -j \'%s\' -d"%(cut,binning,plotvar,x_name,y_name,f_name,json_used)
    else:
        args = "-c \'%s\' -b %s -p \'%s\' -x \'%s\' -y \'%s\' -f \'%s\' -j \'%s\'"%(cut,binning,plotvar,x_name,y_name,f_name,json_used)
    
    #if not i<1:continue
    name_sh, name_jds = "%s/tmp_%s.sh"%(outDir,i),"%s/tmp_%s.jds"%(outDir,i)
    f_sh = sh_s%(cmssw, outDir,execmd+" "+args)
    tmp_f_sh = open(name_sh,"w")
    tmp_f_sh.write(f_sh)
    tmp_f_sh.close()
    f_jds = jds_s%(name_sh,i,"%s"%(outDir),f_name)
    if json_used == 'Silver':f_jds = jds_s%(name_sh,i,"%s"%(outDir),'Silver_'+f_name)
    tmp_f_jds = open(name_jds,"w")
    tmp_f_jds.write(f_jds)
    tmp_f_jds.close()
    #os.system("chmod 755 %s"%name_jds)
    #os.system("chmod 755 %s"%name_sh)
    #os.system("condor_submit %s"%name_jds)
    #os.system("rm -f tmp.jds")
    #os.system("rm -f tmp.sh")
    print "%d/%d"%(i,len(info_json))
    print ("./h2muDraw.py %s"%args)
    os.system("./h2muDraw.py %s&"%args)
