#!/usr/bin/env python 

import os, sys
from ROOT import *
from CATTools.CatAnalyzer.histoHelper import *
ROOT.gROOT.SetBatch(True)

def mergecanvas(arg,accept="",excpt="something"):
    outDir = "%s/src/CATTools/CatAnalyzer/test/h2muDraw/plot" %( os.environ['CMSSW_BASE'] )

    if not arg:
        print "Usage : python merge2Canvas.py <plotvar_muid_njet>"
        print """
        Example : python merge2Canvas.py ll_m_medium_0jet 
                                         _tight _1jet
                                                _2jet
        """ 
        sys.exit(0)

    imagelist = []
    #setMargins(canvas,False)
    for i in os.listdir(outDir):
        if (arg in i) and ('.png' in i):
            if not accept in i:continue
            if excpt in i:continue
            print i
            img = TImage.Open("%s/%s"%(outDir,i))
            #img.Vectorize()
            #name = i.split(".png")[0].split("_")[-1]
            imagelist.append(copy.deepcopy(img))
    imagelist.sort()
    print imagelist
    f_name = ""
    if len(accept)>1:
        f_name = "%s_"%accept
    #sys.exit(2)
    if not '2jet' in arg:
        cat = ['tight', 'loose']
        order = ['BB','BO','BE','OO','OE','EE']
        ratio = 0.7
        W_ref = int(1800*ratio)
        H_ref = int(1300*ratio)
    else:
        cat = ['VBF_tight','ggF_tight','loose']
        order = False
        ratio = 0.7
        W_ref = int(1600*ratio)
        H_ref = int(800*ratio)

    for i in cat:
        canvas = TCanvas("canv","canv",W_ref,H_ref)
        if order:
            canvas.Divide(3,2)
            for j in imagelist:
                for order_i,k in enumerate(order):
                    if (i in j.GetName().split("_")[-2:-1]) and (k in j.GetName()):
                        print j.GetName(), j.GetName().split("_"), i
                        canvas.cd(order_i+1)
                        j.Draw("X")
            canvas.Modified()
            canvas.Update()
            canvas.SaveAs("%s%s_%s.png"%(f_name,arg,i))
        else:
            canvas.Divide(3,1)
            for j_i,j in enumerate(imagelist):
                 canvas.cd(j_i+1)
                 j.Draw("X")
            canvas.Modified()
            canvas.Update()
            canvas.SaveAs("%s%s.png"%(f_name,arg))

#for test
#mergecanvas("ll_m_medium_0jet","","Silver")
#mergecanvas("ll_m_medium_0jet","Silver")

