#!/usr/bin/env python
import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os, getopt, sys
import array
from CATTools.CatAnalyzer.histoHelper import *
from CATTools.CatAnalyzer.histoHelper_h2mu import *
from ROOT import TLorentzVector
#import DYestimation
ROOT.gROOT.SetBatch(True)

#import pdb
#pdb.set_trace()

class drawSysErr:

    def __init__(self,_rfname,_versus):
        if (_rfname==None or _versus==None):
            self.rfname=None
            self.versus=None
            sys.exit(0)
        else :
            self.rfname=_rfname
            self.versus=_versus
            self.color=[ROOT.kBlack,ROOT.kRed,ROOT.kBlue]
            self.binning=[300,0,300]
            self.result=[]
            self.len=0
            self.entries=[]
            self.ofname=""
            self.bin=0
            self.ratiorange=1

    def __del__(self):
        self.getresult()
        with open(self.txt,"r") as f:
            print f.read()
        print "=== complete plotting samples ==="

    def comparison(self,vs,i):
        count=len(self.ltcut)
        if vs=="pu":
            if i<3:
                return 1
            else:
                return 0
        if vs=="mu":
            if i==0 or i==count or i==count*2:
                return 1
            else:
                return 0
        if vs=="mueff":
            if i>2 and i<5:
                return 1
            else:
                return 0
    def choose(self):
        question="""
    Choose the number among the things below.\n
"""

        f=1
        for i in range(len(self.rfname)):
            for j in range(len(self.versus)):
                question+="    %d) sample : %s, comparison : %s\n"%(f,self.rfname[i],self.versus[j])
                f+=1
        question+="    %d) exit\n"%(f)
        number = input(question)
        if number==f:
            return 0
        else:
            self.number = number
            return 1
    def setcmsLumi(self,_cmsLumi):
        self.cmsLumi = _cmsLumi

    def setnumber(self,_number):
        self.number = _number

    def setlistcolor(self,_lcolor):
        self.color = _lcolor

    def setbinning(self,_binning):
        self.binning = _binning
    def setratiorange(self,_ratiorange):
        self.ratiorange = _ratiorange

    def getresult(self):
        self.txt = "sys_err_%d.txt"%(self.bin)
        f_txt = open(self.txt,"w")
        syserr={}
        print>>f_txt, (" %d GeV "%self.bin).center(50,"*")
        row = '{:>%d}; {:>14}; {:>14}; '%(self.len)
        words = ["sample","N_original","MC_stat"]
        for i in range(len(self.versus)):
            row += '{:>14}; {:>14}; '
            words.append("%s_up"%self.versus[i])
            words.append("%s_dn"%self.versus[i])

        print>>f_txt, row.format(*words)
        #initialize dictionary
        re_word = ''
        for i,word in enumerate(words):
            syserr[word]=0
            if i<1: re_word += '{%s:>%d}; '%(word,self.len)
            else: re_word += '{%s:>14}; '%word
        #print re_word
        for i,_rfname in enumerate(self.rfname):
            count = 0
            for j,_entries in enumerate(self.entries):
                if _entries['rfname'] != _rfname:continue
                syserr["sample"],syserr["N_original"],syserr["MC_stat"] = _entries['rfname'],"%.6f"%_entries['N0'],"+- %.4f%%"%_entries['N0_stat']
                for k,_versus in enumerate(self.versus):
                    if _entries['vs'] != _versus:continue
                    syserr["%s_up"%_versus],syserr["%s_dn"%_versus] = "%.6f"%_entries['N1'],"%.6f"%_entries['N2']
                count += 1
            if (count == count and syserr["sample"]!=0):
                print>>f_txt, re_word.format(**syserr)
            for key in syserr:
                syserr[key] = 0    
        f_txt.close()
                    


    def exe(self,_ltname,_ltcut,_bin=125):
        if not self.number:self.choose()
        rfname = self.rfname[(self.number-1)/len(_versus)]
        versus = self.versus[(self.number+1)%len(_versus)]
        if self.len < len(rfname):
            self.len = len(rfname)
        print rfname, versus
        self.ofname = "sys_err_%s_%s"%(rfname,versus)
        ofname = self.ofname

        rfile = ROOT.TFile("histo_%s.root"%rfname)
        self.ltname = _ltname
        self.ltcut = _ltcut
        lcolor = self.color
        binning = self.binning 
        self.bin = _bin
        tr_min=1
        tr_max=len(self.ltname)*len(self.ltcut)
        ratiorange=self.ratiorange

        lh=[]
        lhratio=[]


        for i in range(len(self.ltname)):
            for j in range(len(self.ltcut)):
                h=rfile.Get(self.ltname[i]+"_"+self.ltcut[j])
                h.SetFillColor(0)
                h.SetLineColor(lcolor[j]+i)
                h.SetMarkerColor(lcolor[j]+i)
                h.SetMarkerSize(0.5)
                if len(binning) == 1:
                    #hnew = h.Clone()
                    hnew = h.Rebin(binning[0],"hnew")
                else:
                    hnew = h.Rebin(len(binning)-1,"hnew",array.array('d',binning))
                rename = self.ltname[i]+"_"+self.ltcut[j]
                hnew.SetName(rename.split("*")[-1])
                lh.append(copy.deepcopy(hnew))
             
                hratio=hnew.Clone(self.ltname[i]+"_"+self.ltcut[j])
                hratio.Reset()
                h.SetLineColor(lcolor[j]+i)
                lhratio.append(copy.deepcopy(hratio))

        rfile.Close()

        doRatio = True
        #c=ROOT.TCanvas("c","c",800,600)
        canv = makeCanvas("c", doRatio)
        pads=[canv]
        pads = rootplotcore.divide_canvas(canv, 0.3)
        pads[0].cd()
        h1=lh[0].Clone("h1")
        h1.GetXaxis().SetLabelSize(0)
        if len(binning)>1:
            h1.GetYaxis().SetTitle("entries")
        else:
            h1.GetYaxis().SetTitle("entries/%d GeV"%binning[0])

        #f_txt = open("%s.txt"%(ofname),"w")
        leg=ROOT.TLegend(0.6,0.7,0.9,0.9)
        h1.Draw("hist e")
        leg.SetHeader(rfname)
        rename0=self.ltname[0]+"_"+self.ltcut[0]
        #leg.AddEntry(lh[0],rename0.split("*")[0],"l")    
        leg.AddEntry(lh[0],"normal_weight","l")    
        entries = []
        entries.append(lh[0].GetBinContent(lh[0].FindBin(self.bin)))
        if lh[0].GetBinContent(lh[0].FindBin(self.bin))==0:
            errpercent=0.0
        else:
            errpercent=( 100 * lh[0].GetBinError(lh[0].FindBin(self.bin)) )/ lh[0].GetBinContent(lh[0].FindBin(self.bin))
        #errpercent=lh[0].GetBinError(lh[0].FindBin(self.bin))
        entries.append(errpercent)
        leg.AddEntry(0,"%.4f around %s GeV"%(entries[0],self.bin),"")
        for i in range(tr_min,tr_max):
            if (self.comparison(versus,i)):
                lh[i].Draw("hist same")        
                leg.AddEntry(lh[i],lh[i].GetName(),"l")
                entry = lh[i].GetBinContent(lh[i].FindBin(self.bin))
                leg.AddEntry(0,"%.4f around %s GeV"%(entry,self.bin),"")
                entries.append(entry)
            else:continue
        leg.SetTextSize(0.02)
        legheader=leg.GetListOfPrimitives().First()
        legheader.SetTextAlign(22)
        legheader.SetTextSize(0.04)
        leg.Draw("same")
        pads[0].SetLogy()

        tmpdic={
                'rfname':rfname,
                'vs':versus,
                'N0':entries[0],
                'N0_stat':entries[1],
                'N1':entries[2],
                'N2':entries[3],
               }
        if len(entries)>4:
            for i in range(4,len(entries)):
                tmpdic['N%d'%(i-1)] = entries[i]

        self.entries.append(tmpdic)     
  
        pads[1].cd()
        pads[1].SetGridy()
        f=0
        for i in range(tr_min,tr_max):
            #print lh[i], lh[0], lhratio[i]
            lhratio[i].Divide(lh[i],lh[0],1.,1.,"B")
            #print>>f_txt, "===== %s ====="%(lh[i].GetName())
            #for j in range(len(bins)):
            #    nbin=lhratio[i].FindBin(bins[j])
            #    print>>f_txt, "== mass value : %d ==\n"%bins[j]
            #    print>>f_txt, " sys_err : %f\n"%(lhratio[i].GetBinContent(nbin)-1)
            if (self.comparison(versus,i)):
                if (f==0):
                    lhratio[i].Draw("ep")
                else:
                    lhratio[i].Draw("esame")
            else:continue
            lhratio[i].SetLabelSize(0.03,"Y")
            lhratio[i].GetYaxis().SetTitle("%s/nom"%versus)
            lhratio[i].GetXaxis().SetTitle("M_{#mu#mu} [GeV/c^{2}]")
            lhratio[i].SetMaximum(1.+ratiorange)
            lhratio[i].SetMinimum(1.-ratiorange)
            f+=1
        #f_txt.close()
        for p in pads:
            p.RedrawAxis()
            p.SetGridx()
            p.Modified()
            p.Update()
        iPos = 0
        if( iPos==0 ):
            self.cmsLumi.relPosX = 0.1
        self.cmsLumi.CMS_lumi(pads[0], 0, iPos)

        canv.cd()
        canv.Modified()
        canv.Update()
        canv.SaveAs("%s_%s.png"%(ofname,self.bin))
        canv.Clear()

    
if __name__ == "__main__":
    
    datalumi = 15920

    CMS_lumi.lumi_sqrtS = "%.2f fb^{-1}, #sqrt{s} = 13 TeV 25ns "%(float(datalumi)/1000)
   
    _rfname = [
               "background_ll_m",
               #"higgsx30",
               "higgs_ll_m"
              ]
    _versus = [
               "pu",
               "mu",
               "mueff"
              ]
    #_ltname = ["/nom","/mu_u","/mu_d","/jes_u","/jes_d","/jer_u","/jer_d"]
    _ltname = ["/nom","/mu_u","/mu_d"]

    # the number of both ltcut and ltcolor has to same each other
    _ltcut = ["weight*(mueffweight)","(genweight)*(mueffweight)*(puweight_up)","(genweight)*(mueffweight)*(puweight_dn)","(genweight)*(puweight)*(mueffweight_up)","(genweight)*(puweight)*(mueffweight_dn)"]
    #_ltcut = ["weight","(genweight)*(puweight_up)","(genweight)*(puweight_dn)","mueffweight","mueffweight_up","mueffweight_dn"]
    _lcolor = [ROOT.kBlack,ROOT.kRed,ROOT.kBlue,ROOT.kOrange,ROOT.kViolet,ROOT.kTeal]

    #_binning = [0,30,40,50,60,70,80,90,100,110,120,130,140,160,180,230,300]
    _binning = [10]
    _bins = [115,125,135,150]
    _ratiorange = 0.3
    
    """
    outDir = "plot"
    import json,pprint
    _totrfname = []
    with open('%s/info.json'%outDir) as f_json:
        data = json.load(f_json)
        for i in data:
            fname = i['f_name']
            if not "ll_m" in fname:continue
            if "medium" in fname:continue
            jet01=fname.split("01jet_")
            if len(jet01)>1:
                if len(jet01[1].split("_"))>1:
                    print fname
                    continue
            for j in _rfname:
                _totrfname.append(j+"_"+i['f_name'])
    print _totrfname
    """
    """
    A = drawSysErr(_rfname,_versus)
    A.setlistcolor(_lcolor)
    A.setbinning(_binning)
    while A.choose():
        A.exe(_ltname,_ltcut)
    """
    for _bin in _bins:
        A = drawSysErr(_rfname,_versus)
        A.setcmsLumi(CMS_lumi)
        A.setlistcolor(_lcolor)
        A.setbinning(_binning)
        A.setratiorange(_ratiorange)
        for j in range(1,len(_rfname)*len(_versus)+1):
            A.setnumber(j)
            A.exe(_ltname,_ltcut,_bin) 
        del A
    """
    for _bin in _bins:
        A = drawSysErr(_totrfname,_versus)
        A.setlistcolor(_lcolor)
        A.setbinning(_binning)
        for j in range(1,len(_totrfname)*len(_versus)+1):
            A.setnumber(j)
            A.exe(_ltname,_ltcut,_bin) 
        del A
    """
