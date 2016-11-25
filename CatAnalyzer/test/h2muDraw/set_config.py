#!/usr/bin/env python

import ROOT,os,sys,copy
ROOT.gROOT.SetBatch(True)

ds = []
outDir = "plot"
if not os.path.isdir(outDir):os.mkdir(outDir)

info = [
        #['nvertex','nvertex','(ll_m>50&&step>=5&&filtered==1)','[30,0,30]','no. vertex','events'],
        ['dilep.M()','ll_m','(dilep.M()>50&&step>=5&&filtered==1)','[300,0,300]','mass [GeV]','events'],
        #['dilep.Pt()','ll_pt','(ll_m>50&&step>=5&&filtered==1)','[100,0,100]','diMuon p_{T} [GeV]','events'],
        #['lep1.Pt()','lep1_pt','(ll_m>50&&step>=5&&filtered==1)','[100,0,100]','leading muon p_{T} [GeV]','events'],
        #['lep2.Pt()','lep2_pt','(ll_m>50&&step>=5&&filtered==1)','[100,0,100]','sub-leading muon p_{T} [GeV]','events'],
        #['lep1.Eta()','lep1_eta','(ll_m>50&&step>=5&&filtered==1)','[100,-3,3]','#eta','events'],
        #['lep2.Eta()','lep2_eta','(ll_m>50&&step>=5&&filtered==1)','[100,-3,3]','#eta','events'],
        #['met','met','(ll_m>50&&step>=5&&filtered==1)','[100,0,100]','met [GeV]','events'],
       ]

muid = ["tight","medium"]
muid_cut = ["(isTight==1)","(isMedium==1)"]

# ====== jetcat init. ======
jet2_vbf = "(cat == 1)"
jet2_ggf = "(cat == 2)"
jet2_loose = "(cat == 3)"
jet01_tight = "(cat == 4)"
jet01_loose = "(cat == 5)"

cat = ["2jet_VBF_tight","2jet_ggF_tight","2jet_loose","01jet_tight","01jet_loose"]
cat_cut = [jet2_vbf, jet2_ggf, jet2_loose, jet01_tight, jet01_loose]

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
        for j_i, j in enumerate(cat_cut):
            tmp_info[info_loop][2] += "*%s*%s"%(i,j)
            tmp_info[info_loop][1] += "_%s_%s"%(muid[i_i],cat[j_i])
            info.append(tmp_info[info_loop])
            tmp_info = copy.deepcopy(info)
            for k_i, k in enumerate(cat_eta_cut):
                if j_i>2:
                    tmp_info[info_loop][2] += "*%s*%s*%s"%(i,j,k)
                    tmp_info[info_loop][1] += "_%s_%s_%s"%(muid[i_i],cat[j_i],cat_eta[k_i])
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
import pprint
pprint.pprint(ds)
print>>f_json, json.dumps(ds,indent=4,sort_keys=True)
f_json.close()
print '== done =='
print 'After run this script, run \"run_h2muDraw.py\"'
