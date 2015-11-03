import ROOT, copy, CMS_lumi, tdrstyle

def getTH1(name, binning, tr, br, cut):
    hist = ROOT.TH1F(name, name, binning[0], binning[1], binning[2])
    tr.Project(name, br, cut)
    hist.Sumw2()
    return copy.deepcopy(hist)

def makeCanvas(name, doRatio):
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
        
    canvas = ROOT.TCanvas(name,name,50,50,W,H)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin( L/W )
    canvas.SetRightMargin( R/W )
    canvas.SetTopMargin( T/H )
    canvas.SetBottomMargin( B/H )
    canvas.SetTickx(0)
    canvas.SetTicky(0)
    if doRatio:
        canvas.Divide(1,2)
    return copy.deepcopy(canvas)

def drawTH1(name, cmsLumi, mclist, data, x_name, y_name, doLog=True, doRatio=False):
    iPos = 11
    if( iPos==0 ):
        cmsLumi.relPosX = 0.12

    leg = ROOT.TLegend(0.7,0.7,0.82,0.88)
    leg.SetTextSize(0.035)
    leg.SetTextFont(42)
    leg.SetLineColor(0)
    leg.SetFillColor(0)
    leg.AddEntry(data,"Data","lp")

    hs = ROOT.THStack("hs_%s_mc"%(name), "hs_%s_mc"%(name))
    
    for mc in mclist:
        hnew = mc.Clone("hnew"+mc.GetName())
        hnew.Sumw2(False)
        hs.Add(hnew)
        leg.AddEntry(mc, mc.GetName(), "f")
    
    canvas = makeCanvas(name, doRatio)
    tdrstyle.setTDRStyle()    
    cmsLumi.CMS_lumi(canvas, 0, iPos)
    if doLog:
        canvas.SetLogy()
    
    data.Draw()
    data.GetXaxis().SetTitle(x_name)
    data.GetYaxis().SetTitle(y_name)
    data.GetYaxis().SetTitleOffset(1)
    
    hs.Draw("same")
    data.Draw("esamex0")
    leg.Draw("same")
    canvas.Update()
    canvas.RedrawAxis()
    
    canvas.SaveAs("%s.png"%(name))
