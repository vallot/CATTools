#!/usr/bin/env python
import os, glob

lists = os.listdir(".")

dirs = []
for file in lists :
  if ( os.path.isdir(file) ) :
    dirs.append(file)


for dir in dirs :
  cmd = "hadd jcTree_%s.root %s/jet_charge_tree_*.root"%(dir,dir)
  print cmd
  os.system(cmd)
