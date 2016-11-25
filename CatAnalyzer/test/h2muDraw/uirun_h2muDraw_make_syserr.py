#!/usr/bin/env python

import ROOT,os,sys,copy,json
ROOT.gROOT.SetBatch(True)

execmd = "h2muDraw_syserr.py"
import set_config as conf
outDir = conf.getoutDir()

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
for i,info in enumerate(info_json):
    if i<25:continue
    plotvar,f_name,cut,binning,x_name,y_name = info['plotvar'],info['f_name'],info['cut'],info['binning'],info['x_name'],info['y_name']
    if plotvar == 'll_m':
        args = "-c \'%s\' -b %s -p \"%s\" -x \'%s\' -y \'%s\' -f \'%s\' -j \'%s\' -d"%(cut,binning,plotvar,x_name,y_name,f_name,json_used)
    else:
        args = "-c \'%s\' -b %s -p \"%s\" -x \'%s\' -y \'%s\' -f \'%s\' -j \'%s\'"%(cut,binning,plotvar,x_name,y_name,f_name,json_used)
    print "./%s %s &"%(execmd,args)
    os.system("./%s %s &"%(execmd,args))
    print "%d/%d"%(i+1,len(info_json))
       
#print ("./h2muDraw.py -c \"%s\" -b %s -p %s -x %s -y %s -f %s/%s&"%(cut,binning,plotvar,x_name,y_name,outDir,f_name))
