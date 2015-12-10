#!/usr/bin/env python 
import os, sys
from ROOT import *
from CATTools.CatAnalyzer.histoHelper import *
ROOT.gROOT.SetBatch(True)
from merge2Canvas import mergecanvas

plotvar = ["ll_m","ll_pt","met","nvertex","lep1_pt","lep2_pt","lep1_eta","lep2_eta"]
muid = ["tight","medium"]
jetcat = ["0jet","1jet","2jet"]

for i in plotvar:
    for j in muid:
        for k in jetcat:
            #Silver except
            mergecanvas("%s_%s_%s"%(i,j,k),"","Silver")
            #only Silver accepted
            mergecanvas("%s_%s_%s"%(i,j,k),"Silver")

mergeDir = "re_canvas"
if not os.path.isdir(mergeDir):os.mkdir(mergeDir)
os.system("mv *.png %s"%mergeDir)
