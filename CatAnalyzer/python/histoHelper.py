import math, array, ROOT, copy, CMS_lumi, tdrstyle
import PhysicsTools.PythonAnalysis.rootplot.core as rootplotcore
tdrstyle.setTDRStyle()

def getTH1(title, binning, tree, plotvar, cut, scale = 0.):
    if len(binning) == 3:
        hist = ROOT.TH1D("name", title, binning[0], binning[1], binning[2])
    else:
        hist = ROOT.TH1D("name", title, len(binning)-1, array.array('f', binning))
    tree.Project("name", plotvar, cut)
    if hist.GetSumw2N() == 0:
        hist.Sumw2()
    if scale != 0:
        hist.Scale(scale)
    return copy.deepcopy(hist)

def getTH2(title, binning, tree, plotvar, cut, scale = 0.):
    if len(binning) == 3:
        hist = ROOT.TH2D("name", title, binning[0], binning[1], binning[2], binning[0], binning[1], binning[2])
    elif len(binning) == 2 and len(binning[0]) ==3 :
        hist = ROOT.TH2D("name", title, binning[0][0], binning[0][1], binning[0][2], binning[1][0], binning[1][1], binning[1][2])
    else:
        hist = ROOT.TH2D("name", title, len(binning)-1, array.array('f', binning), len(binning)-1, array.array('f', binning))
    tree.Project("name", plotvar, cut)
    if hist.GetSumw2N() == 0:
        hist.Sumw2()
    if scale != 0:
        hist.Scale(scale)
    return copy.deepcopy(hist)

def makeTH1(filename, treename, title, binning, plotvar, cut, scale = 0.):
    tfile = ROOT.TFile(filename)
    tree  = tfile.Get(treename)
    
    if len(binning) == 3:
        hist = ROOT.TH1D("temp", title, binning[0], binning[1], binning[2])
    else:
        hist = ROOT.TH1D("temp", title, len(binning)-1, array.array('f', binning))      
        
    for var in plotvar.split(','):
        hist.Add(getTH1(title, binning, tree, var, cut, scale))
        
    return copy.deepcopy(hist)

def getEntries(filename, treename):
    tfile = ROOT.TFile(filename)
    tree  = tfile.Get(treename)
    return tree.GetEntriesFast()

def getWeightedEntries(filename, treename, plotvar, weight):
    weighthist = makeTH1(filename, treename, '', [1, 0, 1], plotvar, weight)    
    return weighthist.Integral(-1,2)

def divide_canvas(canvas, ratio_fraction):
    margins = [ROOT.gStyle.GetPadTopMargin(), ROOT.gStyle.GetPadBottomMargin()]
    useable_height = 1 - (margins[0] + margins[1])
    canvas.Clear()
    pad = ROOT.TPad('mainPad', 'mainPad', 0., 0., 1., 1.)
    pad.SetFillStyle(4000)
    pad.Draw()
    pad.SetBottomMargin(margins[1] + ratio_fraction * useable_height)
    pad_ratio = ROOT.TPad('ratioPad', 'ratioPad', 0., 0., 1., 1.);
    pad_ratio.SetFillStyle(4000)
    pad_ratio.Draw()
    pad_ratio.SetTopMargin(margins[0] + (1 - ratio_fraction) * useable_height)
    return pad, pad_ratio

def makeCanvas(name, doRatio=False):
    H_ref = 600;
    if doRatio:
        H_ref = 800;
    W_ref = 800;
    canvas = ROOT.TCanvas(name,name,W_ref,H_ref)    
    return canvas

def setMargins(canvas, doRatio=False):
    H_ref = 600;
    if doRatio:
        H_ref = 800;
    W_ref = 800;
    W = W_ref
    H  = H_ref
    T = 0.08*H_ref
    B = 0.12*H_ref 
    L = 0.12*W_ref
    R = 0.04*W_ref
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin( L/W )
    canvas.SetRightMargin( R/W )
    #canvas.SetTopMargin( T/H )
    #canvas.SetBottomMargin( B/H )
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    return canvas

def setDefAxis(axis, title, offset):
    axis.SetTitle(title)
    axis.SetTitleOffset(offset)
    axis.SetTitleColor(1)
    axis.SetTitleFont(42)
    axis.SetTitleSize(0.043)
    axis.SetLabelColor(1)
    axis.SetLabelFont(42)
    axis.SetLabelOffset(0.007)
    axis.SetLabelSize(0.03)
    axis.SetAxisColor(1)
    axis.SetTickLength(0.03)
    axis.SetNdivisions(510)
    #axis.SetStripDecimals(True)
    #axis.SetPadTickX(1)
    #axis.SetPadTickY(1)

def setDefTH1Style(th1, x_name, y_name):
    setDefAxis(th1.GetYaxis(),y_name, 1.2)
    setDefAxis(th1.GetXaxis(),x_name, 1)
    ROOT.gStyle.SetStripDecimals(True)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.cd()
    return th1
    
def drawTH1(name, cmsLumi, mclist, data, x_name, y_name, doLog=False, doRatio=True, ratioRange=0.45, siglist=None, legx=0.68, legfontsize=0.030):
    leg = ROOT.TLegend(legx,0.68,legx+0.2,0.91)
    leg.SetBorderSize(0)
    #leg.SetNColumns(2)
    leg.SetTextSize(legfontsize)
    leg.SetTextFont(42)
    leg.SetLineColor(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.AddEntry(data,"Data","lp")
    
    leghist = []
    
    if siglist is not None:
        #leg.AddEntry(sig, sig.GetTitle(), "l")
        #leghist.append(sig.GetTitle())
        for i, sig in enumerate(siglist):
            leg.AddEntry(sig, sig.GetTitle(), "l")
            leghist.append(sig.GetTitle())

    hs = ROOT.THStack("hs_%s_mc"%(name), "hs_%s_mc"%(name))
    hratio = mclist[0].Clone("hratio")
    hratio.Reset()
    
    for i, mc in enumerate(mclist):
        hnew = mc.Clone("hnew"+mc.GetName())
        hnew.Sumw2(False)
        hs.Add(hnew)
        hratio.Add(mc)
        inversed = mclist[len(mclist)-1-i]
        if not any(inversed.GetTitle() == s for s in leghist):
            leg.AddEntry(inversed, inversed.GetTitle(), "f")
            leghist.append(inversed.GetTitle())
                        
    hratio.Divide(data,hratio,1.,1.,"B")

    tdrstyle.setTDRStyle()

    setDefTH1Style(data, x_name, y_name)
    data.SetMaximum(data.GetMaximum()*1.8)
    if doLog:
        #data.SetMaximum(10**7)
        #data.SetMinimum(10**-3)
        data.SetMaximum(data.GetMaximum()*100)
    else:
        data.GetYaxis().SetTitleSize(0.04)
        data.GetYaxis().SetLabelSize(0.024)
        data.GetYaxis().SetTitleOffset(1.35)
        
    ratio_fraction = 0
    if doRatio:
        ratio_fraction = 0.3        
        data.GetXaxis().SetLabelSize(0)
        data.GetXaxis().SetTitleSize(0)
        setDefTH1Style(hratio, x_name, "Data/MC")
        hratio.GetYaxis().CenterTitle()
        hratio.GetYaxis().SetNdivisions(5)
            
    canv = makeCanvas(name, doRatio)
    pads=[canv]
    pads = rootplotcore.divide_canvas(canv, ratio_fraction)
    pads[0].cd()
    
    setMargins(pads[0],doRatio)
    if doLog:
        pads[0].SetLogy()

    data.Draw()
    hs.Draw("same")
    
    if siglist is not None:
        for i, sig in enumerate(siglist):
            sig.Draw("samehist")
    
    data.Draw("esamex0")
    leg.Draw("same")

    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.DrawLatex(0.25, 0.85, name.split('_')[0])
    #tex.DrawLatex(canv.GetLeftMargin()*1.4, 1-canv.GetTopMargin()*2.8, name.split('_')[0])
    
    pads[0].Update()

    if doRatio:
        pads[1].cd()
        pads[1].SetGridy()
        setMargins(pads[1],doRatio)
        hratio.SetLineColor(1)
        hratio.Draw("e")
        hratio.SetMaximum(1.+ratioRange)
        hratio.SetMinimum(1.-ratioRange)

    for p in pads:
        p.RedrawAxis()
        p.Modified()
        p.Update()

    canv.cd()
    iPos = 0
    if( iPos==0 ):
        cmsLumi.relPosX = 0.1
    cmsLumi.CMS_lumi(pads[0], 4, iPos)

    canv.Modified()
    canv.Update()
    canv.SaveAs(name+".png")
    canv.SaveAs(name+".pdf")

def drellYanEstimation(mc_ee_in, mc_ee_out, mc_mm_in, mc_mm_out,
                       rd_ee_in, rd_mm_in, rd_em_in, kMM, kEE):
    #kMM = math.sqrt(rd_mm_in/rd_ee_in)/2.
    #kEE = math.sqrt(rd_ee_in/rd_mm_in)/2.

    rMC_mm = mc_mm_out/mc_mm_in
    rMC_ee = mc_ee_out/mc_ee_in
    print "rMC_mm  ", rMC_mm
    print "rMC_ee  ", rMC_ee
    
    nOutEst_mm = rMC_mm*(rd_mm_in - rd_em_in*kMM)
    nOutEst_ee = rMC_ee*(rd_ee_in - rd_em_in*kEE)
    return nOutEst_ee/mc_ee_out,nOutEst_mm/mc_mm_out

def findDataSet(name, datasets):
    for data in datasets:
        if data["name"] == name:
            return data
    return None

def adderrs(err1, err2, sign=1.):
    return math.sqrt(err1**2+sign*err2**2)

def table(mchistList, errList, signal_hist, signal_err):
    nums = {}
    errs = {}
    total = total_err = 0

    titles = list(set([mc.GetTitle() for mc in mchistList]))
    for t in titles:
        nums[t] = 0
        errs[t] = 0

    for i, mc in enumerate(mchistList):
        nbins = mc.GetSize()-2
        nums[mc.GetTitle()] += mc.Integral(0,nbins+1)
        errs[mc.GetTitle()] = adderrs(errs[mc.GetTitle()], errList[i])

        total += mc.Integral(0,nbins+1)
        total_err = adderrs(total_err, errList[i])
    
    nums['total'] = total
    errs['total'] = total_err

    bkg = total - signal_hist.Integral(0,signal_hist.GetSize()-1)
    bkg_err = adderrs(total_err, signal_err, -1)
    nums['bkg'] = bkg
    errs['bkg'] = bkg_err

    return nums, errs

def set_palette(name="", ncontours=999):
    """Set a color palette from a given RGB list
    stops, red, green and blue should all be lists of the same length
    see set_decent_colors for an example"""

    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.34, 0.00]
    # elif name == "whatever":
        # (define more palettes)
    else:
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]

    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)

    npoints = len(s)
    TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    gStyle.SetNumberContours(ncontours)

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

def parameterization(name, data, mclist, x_min, x_max, mean, meanerr, gamma, gammaerr,doLog=True,doRatio=True):
    canv=makeCanvas("c",doRatio)
    pads=[canv]
    pads = rootplotcore.divide_canvas(canv, 0.3)
    pads[0].cd()
    pads[0].SetGridy()

    fp = ROOT.TF1("fp",fParameterization,x_min,x_max,5)

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
    hratio = mclist[0].Clone("hratio")
    hratio.Reset()
    leghist = []
    for i, mc in enumerate(mclist):
        hnew = mc.Clone("hnew"+mc.GetName())
        print mc.GetTitle()
        hnew.Sumw2(False)
        if ("GG_H" in mc.GetTitle()) or ("VBF_H" in mc.GetTitle()):
            hs2.Add(hnew)
        else:
            hs1.Add(hnew)
        hratio.Add(mc)
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

    fp.SetLineColor(ROOT.kAzure)
    fp.SetFillColor(ROOT.kAzure+10)
    #fp.SetParLimits(0,nmc,nmc)
    fp.SetParLimits(1,mean-meanerr,mean+meanerr)
    fp.SetParLimits(2,gamma-gammaerr,gamma+gammaerr)
    fp.SetLineWidth(2)
    fp.SetParNames("N_{bg}","m_{Z}","#gamma","#Lambda","#beta")
    from time import localtime, strftime
    leg.SetHeader(strftime("%Y-%m-%d %H:%M", localtime()))
    #legheader=leg.GetListOfPrimitives().First()
    #legheader.SetTextAlign(22)
    #legheader.SetTextSize(0.04)
    hs2sum = hs2.GetStack().Last()
    hs2sum.SetLineColor(ROOT.kOrange)
    hs2sum.SetLineWidth(3)
    hs2sum.SetFillColor(0)
    leg.AddEntry(fp,"fit","lf")
    leg.AddEntry(hs2sum,"Higgs x 30","l")

    data.Fit("fp","R0")
    #data.Draw("AXIS")
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
    setNameIntoCanvas(pads[0],name)

    pads[1].cd()
    pads[1].SetGridy()
    ratiorange=3
    lhratio=[]
    yratio_up=[]
    yratio_dn=[]
    for i in range(len(x)):
        print i,x[i],y_dn[i],y[i],y_up[i]
        yratio_up.append(y_up[i]/y[i])
        yratio_dn.append(y_dn[i]/y[i])
      
    print len(x)
    from array import array
    gr=ROOT.TGraph(len(x),array('d',x),array('d',y))
    gr_ratio_up=ROOT.TGraph(len(x),array('d',x),array('d',yratio_up))
    gr_ratio_dn=ROOT.TGraph(len(x),array('d',x),array('d',yratio_dn))
        
    #lh = [ROOT.gROOT.GetFunction("fp").GetHistogram(),data,ROOT.gROOT.GetFunction("fp_up").GetHistogram(),ROOT.gROOT.GetFunction("fp_dn").GetHistogram()]
    lh = [gr.GetHistogram(),data,gr_ratio_up.GetHistogram(),gr_ratio_dn.GetHistogram()]
    lcolor = [2,1,ROOT.kBlue-9,ROOT.kBlue-9]
    ldrawopt = ["","epsame","f","f"]
    for i in range(len(lh)):
        hratio=data.Clone("hratio_%d"%i)
        hratio.Reset()
        if i<2:hratio.Divide(lh[0],lh[i],1.,1.,"B")
        else:hratio=lh[i]
        hratio.SetLineColor(lcolor[i])
        hratio.Draw(ldrawopt[i])

    hratio.SetLabelSize(0.03,"Y")
    #hratio.GetYaxis().SetTitle("#frac{data}{fit}")
    hratio.GetXaxis().SetTitle("M_{#mu#mu} [GeV/c^{2}]")
    hratio.SetAxisRange(110,160,"X")
    hratio.SetMaximum(0.+ratiorange)
    hratio.SetMinimum(0.-ratiorange)

    for p in pads:
        p.RedrawAxis()
        p.SetGridx()
        p.Modified()
        p.Update()

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

    gr = ROOT.TGraph()
    gr.SetFillColor(f1.GetFillColor())
    print gr.GetFillColor()
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
    npx = f3.GetNpx()
    npoints=0
    dx = (xmax-xmin)/npx
    x = xmin+0.5*dx
    while (x <= xmax):
        y = f1.Eval(x)
        if (y < ymin): y = ymin
        if (y > ymax): y = ymax
        gr.SetPoint(npoints,x,y)
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
        print i,xp,yp
    ROOT.gDirectory.Add(gr)
    gr.Draw("cf") #draw graph with fill area option
    f3.Draw("lsame") #superimpose function
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


