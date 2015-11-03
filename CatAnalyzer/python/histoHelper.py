import ROOT, copy, CMS_lumi, tdrstyle
tdrstyle.setTDRStyle()    

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

def drawTH1(name, cmsLumi, mclist, data, x_name, y_name, doLog=True):
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
    
    canvas = makeCanvas(name, False)
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
    
def drawTH1withRatio(name, cmsLumi, mclist, data, x_name, y_name, doLog=True, doRatio=True):
    iPos = 11
    if( iPos==0 ):
        cmsLumi.relPosX = 0.12
    canv = makeCanvas(name, doRatio)

    pads=[canv]
    
    if doRatio:
        ratio = 0.4
        canv.Clear()
        
        pads = [ROOT.TPad("%s_main"%canv.GetName(),"%s_main"%canv.GetName(),0.,1.- ratio,1.,1.),
                ROOT.TPad("%s_ratio"%canv.GetName(),"%s_ratio"%canv.GetName(),0.,0.,1.,1.- ratio)]
        for p in pads:
            p.Draw()
            
        pads[0].SetBottomMargin(ratio*0.03)
        pads[1].SetBottomMargin(pads[1].GetBottomMargin()+pads[1].GetTopMargin()+ratio*0.2)
        pads[1].SetTopMargin(ratio*0.03)
        pads[1].SetTickx()
        print pads[1].GetName()

    #cmsLumi.CMS_lumi(pads[0], 0, iPos)
        
    leg = ROOT.TLegend(0.7,0.7,0.82,0.88)
    leg.SetTextSize(0.035)
    leg.SetTextFont(42)
    leg.SetLineColor(0)
    leg.SetFillColor(0)
    leg.AddEntry(data,"Data","lp")
    data.GetXaxis().SetTitle(x_name)
    data.GetYaxis().SetTitle(y_name)
    data.GetYaxis().SetTitleOffset(1)

    hs = ROOT.THStack("hs_%s_mc"%(name), "hs_%s_mc"%(name))
    hratio = mclist[0].Clone("hratio")
    hratio.Reset()
    
    for mc in mclist:
        hnew = mc.Clone("hnew"+mc.GetName())
        hnew.Sumw2(False)
        hs.Add(hnew)
        hratio.Add(mc)
        leg.AddEntry(mc, mc.GetName(), "f")
    
    #if doLog:
    #    canvas.SetLogy()
    pads[0].cd()
    #data.Draw("same")
    #hs.Draw("same")
    #data.Draw("esamex0")
    #leg.Draw("same")
    
    if doRatio:
        pads[1].cd()
        hratio.Divide(hratio,data,1.,1.,"B")
        hratio.Draw("e")
        
    for p in pads:
        print p
        p.RedrawAxis()
        p.Modified()
        p.Update()
    canv.SaveAs("%s.png"%(name))

    print "finished"
