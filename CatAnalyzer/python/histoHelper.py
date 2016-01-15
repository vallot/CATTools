import math, array, ROOT, copy, CMS_lumi, tdrstyle
import PhysicsTools.PythonAnalysis.rootplot.core as rootplotcore
tdrstyle.setTDRStyle()

def getTH1(title, binning, tree, plotvar, cut, scale = 0.):
    if len(binning) == 3:
        hist = ROOT.TH1D("name", title, binning[0], binning[1], binning[2])
    else:
        hist = ROOT.TH1D("name", title, len(binning)-1, array.array('f', binning))
    #tree.Project("name", 'min(%s,%s-0.001)'%(plotvar,binning[2]), cut)
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

def makeCanvas(name, doRatio):
    H_ref = 600;
    if doRatio:
        H_ref = 800;
    W_ref = 800;
    canvas = ROOT.TCanvas(name,name,W_ref,H_ref)    
    return canvas

def setMargins(canvas, doRatio):
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
    th1.GetYaxis().CenterTitle()
    ROOT.gStyle.SetStripDecimals(True)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.cd()
    return th1
    
def drawTH1(name, cmsLumi, mclist, data, x_name, y_name, doLog=False, doRatio=True, ratioRange=0.45):
    #leg = ROOT.TLegend(0.58,0.78,0.8,0.9)
    leg = ROOT.TLegend(0.71,0.68,0.88,0.91)
    leg.SetBorderSize(0)
    #leg.SetNColumns(2)
    leg.SetTextSize(0.029)
    leg.SetTextFont(42)
    leg.SetLineColor(0)
    leg.SetFillColor(0)
    leg.AddEntry(data,"Data","lp")

    hs = ROOT.THStack("hs_%s_mc"%(name), "hs_%s_mc"%(name))
    hratio = mclist[0].Clone("hratio")
    hratio.Reset()
    leghist = []
    for i, mc in enumerate(mclist):
        hnew = mc.Clone("hnew"+mc.GetName())
        hnew.Sumw2(False)
        hs.Add(hnew)
        hratio.Add(mc)
        inversed = mclist[len(mclist)-1-i]
        if not any(inversed.GetTitle() == s for s in leghist):
            leg.AddEntry(inversed, inversed.GetTitle(), "f")
            leghist.append(inversed.GetTitle())
                        
    #hratio.Divide(data,mclist[0],1.,1.,"B")
    hratio.Divide(data,hratio,1.,1.,"B")

    tdrstyle.setTDRStyle()

    setDefTH1Style(data, x_name, y_name)
    data.SetMaximum(data.GetMaximum()*1.8)
    if doLog:
        #data.SetMaximum(10**7)
        data.SetMaximum(data.GetMaximum()*100)
        
    ratio_fraction = 0
    if doRatio:
        ratio_fraction = 0.3        
        data.GetXaxis().SetLabelSize(0)
        data.GetXaxis().SetTitleSize(0)
        data.GetYaxis().CenterTitle()
        setDefTH1Style(hratio, x_name, "Data/MC")
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
    #hs.Draw("samenostack")
    data.Draw("esamex0")
    leg.Draw("same")
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
    iPos = 11
    if( iPos==0 ):
        cmsLumi.relPosX = 0.12    
    cmsLumi.CMS_lumi(canv, 0, iPos)
    
    canv.Modified()
    canv.Update()
    canv.SaveAs(name)

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
