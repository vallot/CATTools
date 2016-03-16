#!/usr/bin/env python

import ROOT,os,sys,copy
ROOT.gROOT.SetBatch(True)

ds = []
outDir = "%s/src/CATTools/CatAnalyzer/test/h2muDraw/plot" %( os.environ['CMSSW_BASE'] )
if not os.path.isdir(outDir):os.mkdir(outDir)

info = [['nvertex','nvertex','(step>=5&&filtered==1)','[30,0,30]','no. vertex','events'],
        ['ll_m','ll_m','(step>=5&&filtered==1)','[200,0,200]','mass [GeV]','events'],
        ['ll_pt','ll_pt','(step>=5&&filtered==1)','[100,0,100]','diMuon p_{T} [GeV]','events'],
        ['lep1_pt','lep1_pt','(step>=5&&filtered==1)','[100,0,100]','leading muon p_{T} [GeV]','events'],
        ['lep2_pt','lep2_pt','(step>=5&&filtered==1)','[100,0,100]','sub-leading muon p_{T} [GeV]','events'],
        ['lep1_eta','lep1_eta','(step>=5&&filtered==1)','[100,-3,3]','#eta','events'],
        ['lep2_eta','lep2_eta','(step>=5&&filtered==1)','[100,-3,3]','#eta','events'],
        ['met','met','(step>=5&&filtered==1)','[100,0,100]','met [GeV]','events'],
       ]

muid = ["tight","medium"]
muid_cut = ["(isTight==1)","(isMedium==1)"]

# ====== jetcat init. ======
jet0_tight = "(cat_jet == 1)"
jet0_loose = "(cat_jet == 2)"
jet1_tight = "(cat_jet == 3)"
jet1_loose = "(cat_jet == 4)"
jet2_vbf = "(cat_jet == 5 || cat_jet == 7)"
jet2_ggf = "(cat_jet == 6)"
jet2_loose = "(cat_jet == 8)"

cat_jet = ["0jet_tight","0jet_loose","1jet_tight","1jet_loose","2jet_VBF_tight","2jet_ggF_tight","2jet_loose"]
cat_jet_cut = [jet0_tight, jet0_loose, jet1_tight, jet1_loose, jet2_vbf, jet2_ggf, jet2_loose]

BB = "(cat_eta == 1)"
BO = "(cat_eta == 2)"
BE = "(cat_eta == 3)"
OO = "(cat_eta == 4)"
OE = "(cat_eta == 5)"
EE = "(cat_eta == 6)"

cat_eta = ["BB","BO","BE","OO","OE","EE"]
cat_eta_cut = [BB,BO,BE,OO,OE,EE]
# ====== jetcat endl. ======

tmp_info = copy.deepcopy(info)
for info_loop in range(len(info)):
    for i_i,i in enumerate(muid_cut):
        tmp_info[info_loop][2] += "*%s"%i
        tmp_info[info_loop][1] += "_%s"%muid[i_i]
        info.append(tmp_info[info_loop])
        tmp_info = copy.deepcopy(info)
        for j_i, j in enumerate(cat_jet_cut):
            tmp_info[info_loop][2] += "*%s*%s"%(i,j)
            tmp_info[info_loop][1] += "_%s_%s"%(muid[i_i],cat_jet[j_i])
            info.append(tmp_info[info_loop])
            tmp_info = copy.deepcopy(info)
            for k_i, k in enumerate(cat_eta_cut):
                if j_i<4:
                    tmp_info[info_loop][2] += "*%s*%s*%s"%(i,j,k)
                    tmp_info[info_loop][1] += "_%s_%s_%s"%(muid[i_i],cat_jet[j_i],cat_eta[k_i])
                    info.append(tmp_info[info_loop])
                    tmp_info = copy.deepcopy(info)
print '== start =='
info.sort()
import json
f_json = open("%s/info.json"%outDir,"w")
for i in range(len(info)):
    ds.append({'plotvar':info[i][0],
               'f_name':info[i][1],
               'cut':info[i][2],
               'binning':info[i][3],
               'x_name':info[i][4],
               'y_name':info[i][5],
    })
    sorted(ds[i])
print ds
print>>f_json, json.dumps(ds)
f_json.close()
print '== done =='
print 'After run this script, run \"run_h2muDraw.py\"'
