#!/usr/bin/python

import os,time

multicrabdir = []
for fn in os.listdir('.'):
    if os.path.isdir(fn):
        if fn.startswith('multicrab_'):
            multicrabdir.append(fn)
            print "found multicrab job dir:", fn

for n in range(50):
    for d in multicrabdir:
        command_submit = "multicrab -submit 500 -c "+d
        print "loop", n, " ,", command_submit
        os.system(command_submit)
        command_status = "multicrab -status -c "+d
        print "loop", n, " ,", command_status        
        os.system(command_status)
        command_get = "multicrab -get -c "+d
        print "loop", n, " ,", command_get
        os.system(command_get)
        command_resubmit = "multicrab -resubmit bad -c "+d
        print "loop", n, " ,", command_resubmit
        os.system(command_resubmit)
        time.sleep(2000)
