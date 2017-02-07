import math, array, ROOT, copy, CMS_lumi, tdrstyle
import PhysicsTools.PythonAnalysis.rootplot.core as rootplotcore
tdrstyle.setTDRStyle()
from histoHelper import *

def fBreitWigner(x,par):
    '''
    par[0] = constant
    par[1] = mean
    par[2] = gamma
    '''
    pi = ROOT.TMath.Pi()
    return par[2]/((x[0]-par[1])*(x[0]-par[1]) + par[2]*par[2]/4) * (2/pi) * par[0]

def fParameterization(x,par):
    '''
    par[0] = constant
    par[1] = mean
    par[2] = gamma
    par[3] = lambda
    par[4] = beta
    '''
    return par[0] * ROOT.TMath.Exp(-par[3]*x[0]) * ( ( par[4]*par[2] / ( (x[0]-par[1])*(x[0]-par[1]) + par[2]*par[2]/4 ) ) + ( (1-par[4])/ (x[0]*x[0]) ) )


def drawBWFit(name, data, x_min, x_max, doLog=False, draw=False):
    c = ROOT.TCanvas(name,name,800,600)
    hd= data.Clone("hdata")
    hd.GetXaxis().SetRangeUser(x_min-1,x_max+1)
    hd.SetTitle("%s;M(#mu#mu)[GeV];events"%(name))

    smean = data.GetMean()
    sgamma = data.GetRMS()
    max = data.GetMaximum()
    print x_min,x_max, smean, sgamma

    fbw = ROOT.TF1("fbw",fBreitWigner,x_min,x_max,3)
    #fbw = ROOT.TF1("fbw","[1]/((x-[0])*(x-[0]) + [1]*[1]/4) * (2*TMath::Pi())",x_min,x_max)
    #fbw.SetParameters(1000,smean,sgamma)
    fbw.SetParLimits(1,smean-10,smean+10)
    fbw.SetParLimits(2,1,x_max-x_min)
    hd.Fit("fbw","R")
    hd.Draw()
    hd.SetMaximum(max*1.1)
    #fbw.Draw("same")
    c.SaveAs(name)
    gmean = fbw.GetParameter(1)
    ggamma = fbw.GetParameter(2)
    gmeanerr = fbw.GetParError(1)
    ggammaerr = fbw.GetParError(2)
    print gmean,gmeanerr,ggamma,ggammaerr
    return gmean,gmeanerr,ggamma,ggammaerr

def parameterization(name,cmsLumi,data, mclist, x_min, x_max, mean, meanerr, gamma, gammaerr,doLog=True,doRatio=True):
    canv=makeCanvas("c",doRatio)
    pads=[canv]
    pads = rootplotcore.divide_canvas(canv, 0.3)
    pads[0].cd()
    pads[0].SetGridy()

    fp = ROOT.TF1("fp",fParameterization,x_min,x_max,5)
    fp.SetNpx(x_max-x_min)

    leg = ROOT.TLegend(0.71,0.60,0.90,0.80)
    leg.SetBorderSize(0)
    #leg.SetNColumns(2)
    leg.SetTextSize(0.029)
    leg.SetTextFont(42)
    leg.SetLineColor(0)
    leg.SetFillColor(0)
    leg.AddEntry(data,"Data","lp")
    hs1 = ROOT.THStack("hs_%s_mc"%(name), "hs_%s_mc"%(name))
    hs2 = ROOT.THStack("hs_%s_mc"%(name), "hs_%s_mc"%(name))
    leghist = []
    for i, mc in enumerate(mclist):
        hnew = mc.Clone("hnew"+mc.GetName())
        print mc.GetTitle()
        hnew.Sumw2(False)
        if ("GG_H" in mc.GetTitle()) or ("VBF_H" in mc.GetTitle()):
            hs2.Add(hnew)
        else:
            hs1.Add(hnew)
        inversed = mclist[len(mclist)-1-i]
        if not any(inversed.GetTitle() == s for s in leghist):
            #leg.AddEntry(inversed, inversed.GetTitle(), "f")
            leghist.append(inversed.GetTitle())
    hsum = hs1.GetStack().Last()
    nmc = hsum.Integral(x_min,x_max)
    print nmc

    # fit error
    fp_up=ROOT.TF1("fp_up",fParameterization,x_min,x_max,5)
    fp_dn=ROOT.TF1("fp_dn",fParameterization,x_min,x_max,5)
    fp_up.SetObjectStat(0)
    fp_dn.SetObjectStat(0)
    fp_up.SetNpx(x_max-x_min)
    fp_dn.SetNpx(x_max-x_min)

    fp.SetLineColor(ROOT.kAzure)
    fp.SetFillColor(ROOT.kAzure+10)
    #fp.SetParLimits(0,nmc,nmc)
    fp.SetParLimits(1,mean-meanerr,mean+meanerr)
    fp.SetParLimits(2,gamma-gammaerr,gamma+gammaerr)
    fp.SetLineWidth(2)
    fp.SetParNames("N_{bg}","m_{Z}","#gamma","#Lambda","#beta")
    from time import localtime, strftime
    #leg.SetHeader(strftime("%Y-%m-%d %H:%M", localtime()))
    #legheader=leg.GetListOfPrimitives().First()
    #legheader.SetTextAlign(22)
    #legheader.SetTextSize(0.04)
    hs2sum = hs2.GetStack().Last()
    hs2sum.SetLineColor(ROOT.kOrange)
    hs2sum.SetLineWidth(3)
    hs2sum.SetFillColor(0)
    leg.AddEntry(fp,"fit","lf")
    leg.AddEntry(hs2sum,"Higgs x 30","l")

    #data.Draw("AXIS")
    data.Fit("fp","R0+")
    data.Draw("ex0")

    par = ""
    parerr = ""
    for i in range(5):
        par = fp.GetParameter(i)
        parerr =fp.GetParError(i)
        print "fp_up.par[%d] : %f"%(i,par+parerr)
        print "fp_dn.par[%d] : %f"%(i,par-parerr)
        fp_up.FixParameter(i,par+parerr)
        fp_dn.FixParameter(i,par-parerr)
        
    fp_up.SetLineColor(ROOT.kOrange)
    fp_dn.SetLineColor(ROOT.kOrange)
    fp_up.SetFillColor(ROOT.kAzure+10)
    fp_up.SetLineWidth(2)
    fp_dn.SetLineWidth(2)

    data.Fit("fp_up","R0+")
    data.Fit("fp_dn","R0+")

    #fp.Draw("same")
    #fp_up.Draw("same")
    #fp_dn.Draw("same")
    hs2sum.Draw("lsame")
    #hs1.Draw("same")
    leg.Draw("same")
    pads[0].Update()

    x, y_up, y_dn, y = shade(pads[0],fp_up,fp_dn,fp)
    data.Draw("e same x0")
    hs2sum.Draw("lsame")
    #hs1.Draw("same")
    leg.Draw("same")
    

    print "fp : %f"%fp.Integral(120,130)
    print "fp_up : %f"%fp_up.Integral(120,130)
    print "fp_dn : %f"%fp_dn.Integral(120,130)

    print "fp : %f"%fp.Integral(120,130)
    print "fp_up percent : %f%%"%((fp_up.Integral(120,130)-fp.Integral(120,130))/fp.Integral(120,130)*100)
    print "fp_dn percent : %f%%"%((fp_dn.Integral(120,130)-fp.Integral(120,130))/fp.Integral(120,130)*100)

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)
    if doLog:
        #data.SetMinimum(10**-2)
        pads[0].SetLogy()
    data.GetXaxis().SetLabelSize(0)
    data.GetYaxis().SetLabelSize(0.03)
    data.GetYaxis().SetTitle("Events/ 1 GeV/c^{2}")
    data.SetAxisRange(110,160,"X")
    
    #pads[0].Modified()
    #setNameIntoCanvas(pads[0],name)

    pads[1].cd()
    pads[1].SetGridy()
    ratiorange=0.3
    lhratio=[]
    yratio_up=[]
    yratio_dn=[]
    for i in range(len(x)):
        #print i,x[i],y_dn[i],y[i],y_up[i]
        yratio_up.append(y_up[i]/y[i])
        yratio_dn.append(y_dn[i]/y[i])
      
    #print len(x)
    from array import array
    gr=ROOT.TGraph(len(x),array('d',x),array('d',y))
    gr_ratio_up=ROOT.TGraph(len(x),array('d',x),array('d',yratio_up))
    gr_ratio_dn=ROOT.TGraph(len(x),array('d',x),array('d',yratio_dn))
    hfit=copy.deepcopy(fp.GetHistogram())
    #if hfit.GetSumw2N() == 0:
    #    hfit.Sumw2()
    hfit_up=copy.deepcopy(fp_up.GetHistogram())
    hfit_dn=copy.deepcopy(fp_dn.GetHistogram())

    lh = [hfit,copy.deepcopy(data),hfit_up,hfit_dn]
    lh[1].GetXaxis().SetLabelSize(0.03)
    if lh[0].GetSumw2N() == 0:
        print "lh0 has no sumw2!!!"
        lh[0].Sumw2(1)
    #lh[0].GetXaxis().SetLimits(110,160)
    #for i in range(len(lh)):
        #lh[i].SetAxisRange(110,160,"X")
        #lh[1].GetXaxis().SetLimits(lh[1].FindBin(110),lh[1].FindBin(300))
        #lh[i]=lh[i].Clone()
        #lh[i].SetBins(50, 110, 160)
    htmp = ROOT.TH1F("tmp","tmp",50,110,160)
    k=0
    for i in range(110,300):
        k+=1
        if i>120 and i<130:continue
        htmp.SetBinContent(k,lh[1].GetBinContent(lh[1].FindBin(i)))
        print k,lh[1].GetBinContent(lh[1].FindBin(i))
    x,y,exl,eyl,exh,eyh=[],[],[],[],[],[]
    for i in range(110,300):
        y1,y2,y3,y4,y5 = lh[0].GetBinContent(lh[0].FindBin(i)), lh[1].GetBinContent(lh[1].FindBin(i)), lh[2].GetBinContent(lh[2].FindBin(i)), lh[3].GetBinContent(lh[3].FindBin(i)), htmp.GetBinContent(htmp.FindBin(i))
        print """ 
    lh[0].FindBin(i) = %d, lh[1].FindBin(i) = %d, lh[2].FindBin(i) = %d, lh[3].FindBin(i) = %d, htmp.FindBin(i) = %d
        """%(lh[0].FindBin(i), lh[1].FindBin(i), lh[2].FindBin(i), lh[3].FindBin(i),htmp.FindBin(i))
        print """ 
    lh[0].GetBinContent(lh[0].FindBin(i)) = %d, lh[1].GetBinContent(lh[1].FindBin(i)) = %d, lh[2].GetBinContent(lh[2].FindBin(i)) = %d, lh[3].GetBinContent(lh[3].FindBin(i)) = %d,  htmp.GetBinContent(htmp.FindBin(i)) = %d
        """%(y1,y2,y3,y4,y5)
        x.append(ROOT.Double(i))
        y.append(ROOT.Double(1))
        exl.append(ROOT.Double(0))
        exh.append(ROOT.Double(0))
        eyl.append(ROOT.Double(abs((y1-y4)/y1)))
        eyh.append(ROOT.Double(abs((y1-y3)/y1)))        
        #eyl.append(ROOT.Double(y4/y1))
        #eyh.append(ROOT.Double(y3/y1))        
    print x,y,exl,exh,eyl,eyh
    print len(x),len(y),len(exl),len(exh),len(eyl),len(eyh)
    band = ROOT.TGraphAsymmErrors(len(x),array('d',x),array('d',y),array('d',exl),array('d',exh),array('d',eyl),array('d',eyh))
    xline = ROOT.TGraph(len(x),array('d',x),array('d',y))
    band.SetFillColor(fp.GetFillColor())
    band.SetLineColor(fp.GetLineColor())
    band.SetLineWidth(fp.GetLineWidth()) 
    xline.SetLineColor(fp.GetLineColor())
    xline.SetLineWidth(fp.GetLineWidth()) 
    #band.SetFillColor(fp_up.GetFillColor())
    lh[1] = htmp.Clone()    
    lh[1].Sumw2(1)
    for i in range(len(lh)):
        lh[i].SetAxisRange(110,160,"X")
    #lh = [gr.GetHistogram(),data,gr_ratio_up.GetHistogram(),gr_ratio_dn.GetHistogram()]
    lcolor = [2,1,ROOT.kBlue-9,ROOT.kBlue-9]
    #lfcolor = [0,0,ROOT.kBlue-9,ROOT.kBlue-9]
    lfcolor = [0,0,ROOT.kAzure+10,0]
    ldrawopt = ["","pesame","hist same","hist same"]
    lhratio=[]
    dataratio=False
    for i in range(len(lh)):
        hratio=lh[1].Clone("hratio_%d"%i)
        hratio.Reset()
        if i == 1:
            dataratio=True
            hratio.Sumw2(0)
            pratio=ROOT.TH1F("tempratio","tempratio",300,0,300)
        else:
            dataratio=False
        htmp1=lh[1].Clone("htmp1_%d"%i)
        htmp2=lh[1].Clone("htmp2_%d"%i)
        htmp1.Reset()
        htmp2.Reset()
        #hratio.Sumw2(1)
        #hratio.Divide(lh[0],lh[i],1.,1.,"B")
        for j in range(110,161):
            valtmp1=lh[0].GetBinContent(lh[0].FindBin(j))
            valtmp2=lh[i].GetBinContent(lh[i].FindBin(j))
            htmp1.Fill(j,valtmp1)
            htmp2.Fill(j,valtmp2)
            errtmp1=lh[0].GetBinError(lh[0].FindBin(j))
            errtmp2=lh[i].GetBinError(lh[0].FindBin(j))
            reerrtmp1=errtmp1/valtmp1
            if valtmp2==0:
                reerrtmp2=0
            else:
                reerrtmp2=errtmp2/valtmp2
            print "mass : %d , htmp1 : %d , htmp2 : %d , errtmp1 : %.4f, errtmp2 : %.4f, reerrtmp1 : %.4f, reerrtmp2 : %.4f"%(j,lh[0].GetBinContent(lh[0].FindBin(j)),lh[1].GetBinContent(lh[1].FindBin(j)),errtmp1,errtmp2,reerrtmp1,reerrtmp2)
            print "valtmp2/valtmp1 = %.4f , errtmp2/valtmp1 = %.4f"%(valtmp2/valtmp1,errtmp2/valtmp1)
            if dataratio:
                hratio.SetBinContent(j,float(valtmp2/valtmp1))
                hratio.SetBinError(j,float(errtmp2/valtmp1))
                pratio.SetBinContent(j,float(valtmp2/valtmp1))
                pratio.SetBinError(j,float(errtmp2/valtmp1))
            else:
                dataratio=False

        if dataratio:dataratio=True
        else:
            hratio.Divide(htmp1,htmp2,1.,1.,"B")
            
        hratio.SetLineColor(lcolor[i])
        hratio.SetFillColor(lfcolor[i])
        lhratio.append(copy.deepcopy(hratio))
    pads[1].cd()
    hratio.Reset()
    hratio.Draw("")
    ROOT.gDirectory.Add(band)
    band.Draw("4")
    xline.Draw("l")
    #lhratio[1].Draw("%s"%ldrawopt[1])
    pratio.Draw("epsame")
    for i in range(len(lh)+2):
        canvas=ROOT.TCanvas()
        canvas.cd()
        if i<len(lh):lh[i].Draw()
        elif i==len(lh):
            band.Draw("a4")
            xline.Draw("l")
            #band.Draw("L")
        elif i==len(lh)+1:
            lhratio[1].Draw()
        canvas.SaveAs("hratio_%d.png"%i)
        pads[1].cd()
        
    print """
    data.GetNbinsX() = %d
    hfit.GetNbinsX() = %d
    hfit_up.GetNbinsX() = %d
    hfit_dn.GetNbinsX() = %d
    hratio.GetNbinsX() = %d
    """%(lh[1].GetNbinsX(),lh[0].GetNbinsX(),lh[2].GetNbinsX(),lh[3].GetNbinsX(),hratio.GetNbinsX()) 

    hratio.SetLabelSize(0.03,"Y")
    hratio.SetLabelSize(0.04,"X")
    hratio.GetYaxis().SetTitle("#frac{data}{fit}")
    hratio.GetYaxis().SetNdivisions(9)
    hratio.GetXaxis().SetTitle("M_{#mu#mu} [GeV/c^{2}]")
    #hratio[0].SetAxisRange(110,160,"X")
    hratio.SetMaximum(1.+ratiorange)
    hratio.SetMinimum(1.-ratiorange)

    for p in pads:
        p.RedrawAxis()
        p.SetGridx()
        p.Modified()
        p.Update()

    iPos = 0
    if( iPos==0 ):
        cmsLumi.relPosX = 0.1
    cmsLumi.CMS_lumi(pads[0], 0, iPos)

    canv.cd()
    canv.Modified()
    canv.Update()

    canv.SaveAs(name)

def setNameIntoCanvas(pad,name):
    if "." in name:
        fname=name.split(".")[0]
    else:
        fname=name
    print fname

    H = pad.GetWh()
    W = pad.GetWw() 
    l = pad.GetLeftMargin()
    t = pad.GetTopMargin()
    r = pad.GetRightMargin()
    b = pad.GetBottomMargin()
    #pad.DrawFrame(0,0,1,1)

    print "="*50
    print H, W
    print l, t, r, b
    print "="*50

    #need to fix the code for putting text in canvas.
    latex=ROOT.TLatex()
    latex.SetTextSize(0.045)
    latex.DrawLatexNDC(r,1-t,fname)
    #latex.DrawLatex(W/2,H/2,fname)
    #pad.Modified()
    pad.Update()
    return pad

def setLastHist(mclist): 
    hs = ROOT.THStack("hs_mc", "hs_mc")
    hratio = mclist[0].Clone("hratio")
    hratio.Reset()
    leghist = []
    for i, mc in enumerate(mclist):
        #hnew = mc.Clone("hnew"+mc.GetName())
        #hnew.Sumw2(False)
        #hs.Add(hnew)
        hratio.Add(mc)
        inversed = mclist[len(mclist)-1-i]
        if not any(inversed.GetTitle() == s for s in leghist):
            #leg.AddEntry(inversed, inversed.GetTitle(), "f")
            leghist.append(inversed.GetTitle())
    return hratio

def shade(canv, f1, f2, f3):
    ## Rene's code
    # ref : https://root.cern.ch/phpBB3/viewtopic.php?t=5520
    # shade the area between f1 and f2 and draw f3 on top
    # create a TGraph to store the function values
    # shaded area is the fill/color/style of f1

    mg = ROOT.TMultiGraph()
    gr = ROOT.TGraph()
    gr0 = ROOT.TGraph()
    gr.SetFillColor(f1.GetFillColor())
    gr0.SetFillColor(f3.GetFillColor())
    gr0.SetLineColor(f3.GetLineColor())
    gr0.SetLineWidth(f3.GetLineWidth()) 

    #get picture range
    xmin = canv.GetUxmin()
    xmax = canv.GetUxmax()
    ymin = canv.GetUymin()
    ymax = canv.GetUymax()
    print xmin, xmax, ymin, ymax
    gr.SetTitle("banded graph")
    lx = []
    lyf1 = []
    lyf2 = []
    lyf3 = []
   
    #process first function
    print f3.ClassName()
    if "TF" in f3.ClassName():
        npx = f3.GetNpx()
    else:
        npx = f3.GetN()
    npoints=0
    dx = (xmax-xmin)/npx
    x = xmin+0.5*dx
    while (x <= xmax):
        y = f1.Eval(x)
        y0 = f3.Eval(x)
        if (y < ymin): y = ymin
        if (y > ymax): y = ymax
        gr.SetPoint(npoints,x,y)
        gr0.SetPoint(npoints,x,y0)
        lx.append(x)
        lyf1.append(f1.Eval(x))
        lyf2.append(f2.Eval(x))
        lyf3.append(f3.Eval(x))
        npoints+=1
        x += dx

    #process second function
    x = xmax-0.5*dx
    while (x >= xmin):
        y = f2.Eval(x)
        if (y < ymin): y = ymin
        if (y > ymax): y = ymax
        gr.SetPoint(npoints,x,y)
        npoints+=1
        x -= dx
    
    xp = ROOT.Double() 
    yp = ROOT.Double()
    for i in range(npoints):
        gr.GetPoint(i,xp,yp)
        #print i,xp,yp
    #ROOT.gDirectory.Add(gr)
    #ROOT.gDirectory.Add(gr0)
    ROOT.gDirectory.Add(mg)
    #gr.Draw("cf") #draw graph with fill area option
    #gr0.Draw("l")
    #f3.Draw("lsame") #superimpose function
    mg.Add(gr,"cf")
    mg.Add(gr0,"cl")
    mg.Draw("")
    canv.Update()
    return lx, lyf1, lyf2, lyf3
    
def DrawErrorBand(graph):
    isErrorBand = graph.GetErrorYhigh(0) != -1 and graph.GetErrorYlow(0) != -1
    npoints     = graph.GetN()
 
    if not isErrorBand:
        graph.Draw("l same")
        return
 
    # Declare individual TGraph objects used in drawing error band
    central, min, max = ROOT.TGraph(), ROOT.TGraph(), ROOT.TGraph()
    shapes = []
    for i in range((npoints-1)*4):
        shapes.append(ROOT.TGraph())
 
    # Set ownership of TGraph objects
    ROOT.SetOwnership(central, False)
    ROOT.SetOwnership(    min, False)
    ROOT.SetOwnership(    max, False)
    for shape in shapes:
        ROOT.SetOwnership(shape, False)
 
    # Get data points from TGraphAsymmErrors
    x, y, ymin, ymax = [], [], [], []
    for i in range(npoints):
        tmpX, tmpY = ROOT.Double(0), ROOT.Double(0)
        graph.GetPoint(i, tmpX, tmpY)
        x.append(tmpX)
        y.append(tmpY)
        ymin.append(tmpY - graph.GetErrorYlow(i))
        ymax.append(tmpY + graph.GetErrorYhigh(i))
 
    # Fill central, min and max graphs
    for i in range(npoints):
        central.SetPoint(i, x[i], y[i])
    min.SetPoint(i, x[i], ymin[i])
    max.SetPoint(i, x[i], ymax[i])
 
    # Fill shapes which will be shaded to create the error band
    for i in range(npoints-1):
        for version in range(4):
            shapes[i+(npoints-1)*version].SetPoint((version+0)%4, x[i],   ymax[i])
            shapes[i+(npoints-1)*version].SetPoint((version+1)%4, x[i+1], ymax[i+1])
            shapes[i+(npoints-1)*version].SetPoint((version+2)%4, x[i+1], ymin[i+1])
            shapes[i+(npoints-1)*version].SetPoint((version+3)%4, x[i],   ymin[i])
 
    # Set attributes to those of input graph
    central.SetLineColor(graph.GetLineColor())
    central.SetLineStyle(graph.GetLineStyle())
    central.SetLineWidth(graph.GetLineWidth())
    min.SetLineColor(graph.GetLineColor())
    min.SetLineStyle(graph.GetLineStyle())
    max.SetLineColor(graph.GetLineColor())
    max.SetLineStyle(graph.GetLineStyle())
    for shape in shapes:
        shape.SetFillColor(graph.GetFillColor())
        shape.SetFillStyle(graph.GetFillStyle())
 
    # Draw
    for shape in shapes:
        shape.Draw("f same")
    min.Draw("l same")
    max.Draw("l same")
    central.Draw("l same")
    ROOT.gPad.RedrawAxis()


