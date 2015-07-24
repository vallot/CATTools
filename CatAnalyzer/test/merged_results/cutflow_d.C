{
 
  /////////////////////////////////////////
    
    // Stacked plot (data + MC)
    // by Youn Jung Roh (youn.jung.roh@cern.ch
    // 2015 July 21st
    
  /////////////////////////////////////////

  #include <TMath.h>
  #include <TH2F.h>
  #include <TH1F.h>
  #include <iomanip>
    
  gStyle->SetCanvasBorderMode(0); gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadColor(kWhite);    gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);  gStyle->SetFrameFillColor(0);
    
  gStyle->SetFrameFillStyle(0);   gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);   gStyle->SetFrameLineWidth(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11100011);
    
  		
  int ctt=2; int ctto =6; int cwj=0, cdy=4, ctw=5, ctbw=3, cdb=7;
    
  for(int levt=1; levt<2; levt++){

    if(levt==0) int leg_with_evt = 0;
    if(levt==1) int leg_with_evt = 1;
    
    for(int ch1=2; ch1<3; ch1++){
        
      char ch[100];
      if(ch1==0)sprintf(ch, "MuMu");
      if(ch1==1)sprintf(ch, "ElEl");
      if(ch1==2)sprintf(ch, "MuEl");
      if(ch1==3)sprintf(ch, "Combined");
 
      //for(int pl=0; pl<8; pl++){
        
      int pl=1; if(1){

      int sbin=0;
      double max=1.5;
        
      if(pl==0){
          
        int bin=8;
        Double_t xAxis_mll[9] = {20, 30, 50, 76, 106, 130, 170, 260, 400};
        char title1[100]; sprintf(title1, "s1_llm");
        char title2[100]; sprintf(title2, "c1_llm");
        char xAxis[100]; sprintf(xAxis, Form("%s (S1) M(ll) [GeV/c^{2}] ",ch));
        TH1F *psdata     = new TH1F("data","",bin, xAxis_mll); psdata->Sumw2();
        TH1F *mc     = new TH1F("mc","",bin, xAxis_mll); mc->Sumw2();
        TH1F *ratio     = new TH1F("ratio","",bin, xAxis_mll); ratio->Sumw2();

      }
        
      if(pl==1){
        int bin=8;
        Double_t xAxis_mll[9] = {20, 30, 50, 76, 106, 130, 170, 260, 400};
        char title1[100]; sprintf(title1, "s2_llm");
        char title2[100]; sprintf(title2, "c2_llm");
        char xAxis[100]; sprintf(xAxis, Form("%s (S2) M(ll) [GeV/c^{2}] ",ch));
        TH1F *psdata     = new TH1F("data","",bin, xAxis_mll); psdata->Sumw2();
        TH1F *mc     = new TH1F("mc","",bin, xAxis_mll); mc->Sumw2();
        TH1F *ratio     = new TH1F("ratio","",bin, xAxis_mll); ratio->Sumw2();
      }

      if(pl==2){
        int bin=10; max=400;
        char title1[100]; sprintf(title1, "s2_njet");
        char title2[100]; sprintf(title2, "c2_njet");
        char xAxis[100]; sprintf(xAxis, Form("%s (S2) Jet Multiplicity",ch));
        TH1F *psdata     = new TH1F("data","",bin, 0,10); psdata->Sumw2();
        TH1F *mc     = new TH1F("mc","",bin, 0,10); mc->Sumw2();
        TH1F *ratio     = new TH1F("ratio","",bin, 0,10); ratio->Sumw2();
      }

      if(pl==3){
        int bin=10; sbin=2; max=400;
        char title1[100]; sprintf(title1, "s3_njet");
        char title2[100]; sprintf(title2, "c3_njet");
        char xAxis[100]; sprintf(xAxis, Form("%s (S3) Jet Multiplicity",ch));
        TH1F *psdata     = new TH1F("data","",bin, 0,10); psdata->Sumw2();
        TH1F *mc     = new TH1F("mc","",bin, 0,10); mc->Sumw2();
        TH1F *ratio     = new TH1F("ratio","",bin, 0,10); ratio->Sumw2();
      }

      if(pl==4){
        int bin=30; sbin=0;
        char title1[100]; sprintf(title1, "s3_met");
        char title2[100]; sprintf(title2, "c3_met");
        char xAxis[100]; sprintf(xAxis, Form("%s (S3)  Missing E_{T} [GeV]",ch));
        TH1F *psdata     = new TH1F("data","",bin, 0,300); psdata->Sumw2();
        TH1F *mc     = new TH1F("mc","",bin, 0,300); mc->Sumw2();
        TH1F *ratio     = new TH1F("ratio","",bin, 0,300); ratio->Sumw2();
      }
        
      if(pl==5){
        int bin=30; sbin=0;
        char title1[100]; sprintf(title1, "s4_met");
        char title2[100]; sprintf(title2, "c4_met");
        char xAxis[100]; sprintf(xAxis, Form("%s (S4)  Missing E_{T} [GeV]",ch));
        TH1F *psdata     = new TH1F("data","",bin, 0,300); psdata->Sumw2();
        TH1F *mc     = new TH1F("mc","",bin, 0,300); mc->Sumw2();
        TH1F *ratio     = new TH1F("ratio","",bin, 0,300); ratio->Sumw2();
      }

      if(pl==6){
        int bin=6; sbin=0; max=220;
        char title1[100]; sprintf(title1, "s6_nbjet");
        char title2[100]; sprintf(title2, "c6_nbjet");
        char xAxis[100]; sprintf(xAxis, Form("%s (S4) b Jet Multiplicity",ch));
        TH1F *psdata     = new TH1F("data","",bin, 0,6); psdata->Sumw2();
        TH1F *mc     = new TH1F("mc","",bin, 0,6); mc->Sumw2();
        TH1F *ratio     = new TH1F("ratio","",bin, 0,6); ratio->Sumw2();
      }

      if(pl==7){
        int bin=6; sbin=1; max=200;
        char title1[100]; sprintf(title1, "s5_nbjet");
        char title2[100]; sprintf(title2, "c5_nbjet");
        char xAxis[100]; sprintf(xAxis, Form("%s (S5) b Jet Multiplicity",ch));
        TH1F *psdata     = new TH1F("data","",bin, 0,6); psdata->Sumw2();
        TH1F *mc     = new TH1F("mc","",bin, 0,6); mc->Sumw2();
        TH1F *ratio     = new TH1F("ratio","",bin, 0,6); ratio->Sumw2();
      }
        
          
          
          
          
          
      char ofile[100]; sprintf(ofile,"p_%s_%s_%d.png",title1,ch,leg_with_evt);
    
        
      TFile *s0   = new TFile(Form("em_%s_mc.root",ch)); //data
      TFile *s1   = new TFile(Form("ttbar_powheg_%s_mc.root",ch)); //mc
      TFile *s2   = new TFile(Form("tw_%s_mc.root",ch));
      TFile *s3   = new TFile(Form("tbw_%s_mc.root",ch));
      TFile *s4   = new TFile(Form("dy_%s_mc.root",ch));
      TFile *s5   = new TFile(Form("wjet_%s_mc.root",ch));
      TFile *s6   = new TFile(Form("ttbar_other_%s_mc.root",ch));
	  TFile *s7   = new TFile(Form("ww_%s_mc.root",ch));
	  TFile *s8   = new TFile(Form("wz_%s_mc.root",ch));
	  TFile *s9   = new TFile(Form("zz_%s_mc.root",ch));
          
      TH1F *n1 = (TH1F*) s1->Get("h_nevent"); //ca1->Sumw2();
      TH1F *n2 = (TH1F*) s2->Get("h_nevent"); //ca2->Sumw2();
      TH1F *n3 = (TH1F*) s3->Get("h_nevent"); //ca3->Sumw2();
      TH1F *n4 = (TH1F*) s4->Get("h_nevent"); //ca4->Sumw2();
      TH1F *n5 = (TH1F*) s5->Get("h_nevent"); //ca5->Sumw2();
      TH1F *n6 = (TH1F*) s6->Get("h_nevent"); //ca6->Sumw2();
	  TH1F *n7 = (TH1F*) s7->Get("h_nevent");
	  TH1F *n8 = (TH1F*) s8->Get("h_nevent");
	  TH1F *n9 = (TH1F*) s9->Get("h_nevent");

      int ZJets = n4->GetEntries();//2829164;
      int SingleToptW = n2->GetEntries();//986100;
      int SingleTopBartW = n3->GetEntries();//971800;
      int WJets = n5->GetEntries();//10017930;
      int TTbarFullLepMGDecays =n1->GetEntries();//25446993;
      int TTbarSemiLeptMGDecays =n6->GetEntries();//25446993;
      int TTbarHadronicMGDecays =n6->GetEntries();//25446993;
	  int WW=n7->GetEntries();
	  int WZ=n8->GetEntries();
	  int ZZ=n9->GetEntries();
	
      float lumi1 = (float(TTbarFullLepMGDecays)/831.76); //26 //23.6196
      float lumi2 = (float(SingleToptW)/35.6);
      float lumi3 = (float(SingleTopBartW)/35.6);
      float lumi4 = (float(ZJets)/6025.2);
      float lumi5 = (float(WJets)/61526.7);
      float lumi6 = (float(TTbarHadronicMGDecays)/831.76);
	  float lumi7 = (float(WW)/118.7);
	  float lumi8 = (float(WZ)/65.9);
	  float lumi9 = (float(ZZ)/31.8);

      float intlumi=41;//49.502;
    
      TH1F *pca1 = (TH1F*) s0->Get(title1);
      TH1F *ca0 = (TH1F*) s1->Get(title2); //ca1->Sumw2();
      TH1F *ca1 = (TH1F*) s1->Get(title1); //ca1->Sumw2();
      TH1F *ca2 = (TH1F*) s2->Get(title1); //ca2->Sumw2();
      TH1F *ca3 = (TH1F*) s3->Get(title1); //ca3->Sumw2();
      TH1F *ca4 = (TH1F*) s4->Get(title1); //ca4->Sumw2();
      TH1F *ca5 = (TH1F*) s5->Get(title1); //ca5->Sumw2();
      TH1F *ca6 = (TH1F*) s6->Get(title1); //ca6->Sumw2();
      TH1F *ca7 = (TH1F*) s7->Get(title1);
      TH1F *ca8 = (TH1F*) s8->Get(title1);
      TH1F *ca9 = (TH1F*) s9->Get(title1);

      pca1->SetBinContent(bin,pca1->GetBinContent(bin)+pca1->GetBinContent(bin+1));
      ca0->SetBinContent(bin,ca0->GetBinContent(bin)+ca0->GetBinContent(bin+1));
      ca1->SetBinContent(bin,ca1->GetBinContent(bin)+ca1->GetBinContent(bin+1));
      ca2->SetBinContent(bin,ca2->GetBinContent(bin)+ca2->GetBinContent(bin+1));
      ca3->SetBinContent(bin,ca3->GetBinContent(bin)+ca3->GetBinContent(bin+1));
      ca4->SetBinContent(bin,ca4->GetBinContent(bin)+ca4->GetBinContent(bin+1));
      ca5->SetBinContent(bin,ca5->GetBinContent(bin)+ca5->GetBinContent(bin+1));
      ca6->SetBinContent(bin,ca6->GetBinContent(bin)+ca6->GetBinContent(bin+1));
	  ca7->SetBinContent(bin,ca7->GetBinContent(bin)+ca7->GetBinContent(bin+1));
	  ca8->SetBinContent(bin,ca8->GetBinContent(bin)+ca8->GetBinContent(bin+1));
	  ca9->SetBinContent(bin,ca9->GetBinContent(bin)+ca9->GetBinContent(bin+1));
        
      ca0->Sumw2();
      ca1->Sumw2();
      ca2->Sumw2();
      ca3->Sumw2();
      ca4->Sumw2();
      ca5->Sumw2();
      ca6->Sumw2();
	  ca7->Sumw2();
	  ca8->Sumw2();
	  ca9->Sumw2();
        
      ca0->Scale(intlumi/lumi1);
      ca1->Scale(intlumi/lumi1);
      ca2->Scale(intlumi/lumi2);
      ca3->Scale(intlumi/lumi3);
      ca4->Scale(intlumi/lumi4);
      ca5->Scale(intlumi/lumi5);
      ca6->Scale(intlumi/lumi6);
	  ca7->Scale(intlumi/lumi7);
	  ca8->Scale(intlumi/lumi8);
	  ca9->Scale(intlumi/lumi9);

      TH1F *mcall = (TH1F*) ca1->Clone("mcall");
      mcall->Add(ca2);
      mcall->Add(ca3);
      mcall->Add(ca4);
      mcall->Add(ca5);
      mcall->Add(ca6);
	  mcall->Add(ca7);
	  mcall->Add(ca8);
	  mcall->Add(ca9);

     // pca1->Scale(mcall->Integral()/pca1->Integral());

      TCanvas* c3  = new TCanvas("c3", "c3", 600,800);
      c3->Divide(1,2,0,0,0);
        
      Int_t chatch = 1756;
      TColor *color = new TColor(chatch, 0.3, 0.3, 0.3, "", 0.25); // alpha = 0.5
      c3->cd(1); //gPad->SetLogy();
      p11_1 = (TPad*)c3->GetPad(1);
      p11_1->SetPad(0.04,0.3,0.98,0.96);
      //p11_1->SetLogy();
      p11_1->SetRightMargin(0.05);
      p11_1->SetTopMargin(0.05);
 
      THStack *hs2 = new THStack("hs2","");
      //gPad->SetLogx();

      Int_t ci = TColor::GetColor("#ffcc33");
      //ca1->SetTitle(title1);
      //ca1->SetMaximum(ca1->GetMaximum()*10);
      //ca0->SetMinimum(10);
      //ca0->SetMinimum(10);
      //ca0->SetMarkerStyle(20);
      //ca0->SetMarkerSize(1.0);
      //ca1->Draw();
        
      ca1->SetFillColor(ctto);
      ca2->SetFillColor(ctw);
      ca3->SetFillColor(ctbw);
      ca4->SetFillColor(cdy);
      ca5->SetFillColor(cwj);
      ca6->SetFillColor(ctto);
      ca7->SetFillColor(cdb);
      ca8->SetFillColor(cdb);
      ca9->SetFillColor(cdb);
    
      ca1->SetLineColor(ctto);
      ca2->SetLineColor(ctw);
      ca3->SetLineColor(ctbw);
      ca4->SetLineColor(cdy);
      ca5->SetLineColor(cwj);
      ca6->SetLineColor(ctto);
      ca7->SetLineColor(cdb);
	  ca8->SetLineColor(cdb);
	  ca9->SetLineColor(cdb);

      ca1->SetMaximum(TMath::Max(ca4->GetMaximum(),pca1->GetMaximum())*3);
    
      hs2->Add(ca1);//ttbar
      hs2->Add(ca6);//others
      hs2->Add(ca5);//wj
      hs2->Add(ca3);//tbw
      hs2->Add(ca2);//tw
      hs2->Add(ca7);
      hs2->Add(ca8);
      hs2->Add(ca9);
      hs2->Add(ca4);//dy
        
      hs2->Draw("hist");
        
      ca0->SetLineColor(ctt);
      ca0->SetFillColor(ctt);
    
      ca0->Draw("Same hist");
        
          
      for(int b=sbin; b<bin; b++){
        
            
        ratio->Divide(pca1,mcall,1,1,"B");
            
            
      }
    
      pca1->SetMarkerStyle(20);
      pca1->SetMarkerSize(1.0);
      pca1->Draw("Sames E");
    
      TLatex *   tex = new TLatex(0.9,0.955,"#sqrt{s} = 13 TeV");
  
      tex->SetNDC();
      tex->SetTextAlign(31);
      tex->SetTextFont(42);
      
      tex->SetTextSize(0.045);
      tex->SetLineWidth(2);
      //tex->Draw();

      tex = new TLatex(0.9,0.96,Form("%.1f pb^{-1},   #sqrt{s} = 13 TeV",intlumi));
      tex->SetNDC();
      tex->SetTextAlign(31);
      tex->SetTextFont(42);
      tex->SetTextSize(0.05);
      tex->SetLineWidth(2);
      tex->Draw();

      tex = new TLatex(0.42,0.96,"CMS Preliminary");
      tex->SetNDC();
      tex->SetTextAlign(31);
      tex->SetTextFont(42);
      tex->SetTextSize(0.05);
      tex->SetLineWidth(2);
      tex->Draw();
      
        
      char data[100]; char data2[100];
      int evt = ca0->GetEntries()*(intlumi/lumi1);
      int da = pca1->GetEntries();
      int evt_all = ca1->GetEntries()*(intlumi/lumi1);
      int ntw = ca2->GetEntries()*(intlumi/lumi2);
      int ntbw = ca3->GetEntries()*(intlumi/lumi3);
      int ndy= ca4->GetEntries()*(intlumi/lumi4);
      int ntto = (evt_all-evt)+(ca6->GetEntries()*(intlumi/lumi6));
      int nwj = ca5->GetEntries()*(intlumi/lumi5);
      int nww = ca7->GetEntries()*(intlumi/lumi7);
      int nwz = ca8->GetEntries()*(intlumi/lumi8);
      int nzz = ca9->GetEntries()*(intlumi/lumi9);
        
      sprintf(data,"Data = %d",da);
        
      int evt2 = evt+ntto+ntw+ntbw+ndy+nwj+nww+nwz+nzz;
      sprintf(data2," MC Event = %d",evt2);
      //sprintf(data2,"MuMu Step 1",evt2);

      tex = new TLatex(0.81,0.8,data);
      tex->SetNDC();
      tex->SetTextAlign(31);
      tex->SetTextSize(0.03);
      tex->SetLineWidth(2);
      tex->Draw();
  
      tex = new TLatex(0.81,0.76,data2);
      tex->SetNDC();
      tex->SetTextAlign(31);
      tex->SetTextSize(0.03);
      tex->SetLineWidth(2);
      tex->Draw();


      tex = new TLatex(0.04,0.3341969,"Number of Events");
      tex->SetNDC();
      tex->SetTextFont(42);
      tex->SetTextSize(0.06);
      tex->SetLineWidth(2);
      tex->SetTextAngle(90);
      tex->Draw();

      TLegend *leg = new TLegend(0.6,0.4,0.93,0.93,NULL,"brNDC");
      leg->SetTextSize(0.045);
      
      leg->SetTextFont(42);
      leg->SetLineColor(0);
      leg->SetLineStyle(1);
      leg->SetLineWidth(1);
      leg->SetFillColor(0);
      leg->SetFillStyle(1001);
 
      if(!leg_with_evt){
      TLegendEntry *entry=leg->AddEntry("NULL",Form("Data"),"l p");
        
      entry->SetLineStyle(1);
      entry->SetLineWidth(1);
      entry->SetMarkerColor(1);
      entry->SetMarkerStyle(20);
      entry->SetMarkerSize(1);
    
      entry=leg->AddEntry("NULL",Form("ZJet"),"f");
      ci = TColor::GetColor("#ffcc33");
      entry->SetFillColor(cdy);////////////////////////sort of yellow
      entry->SetFillStyle(1001);
      ci = TColor::GetColor("#663300");
      entry->SetLineColor(ci);
      entry->SetLineStyle(1);
      entry->SetLineWidth(1);
      entry->SetMarkerColor(1);
      entry->SetMarkerStyle(21);
      entry->SetMarkerSize(1);

      entry=leg->AddEntry("NULL",Form("WW/WZ/ZZ"),"f");
      entry->SetFillColor(cdb);/////////////////////// light blue
      entry->SetFillStyle(1001);
      ci = TColor::GetColor("#663300");
      entry->SetLineColor(ci);
      entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1);


  entry=leg->AddEntry("NULL",Form("Single top W"),"f");
  entry->SetFillColor(ctw);/////////////////////// light blue
  entry->SetFillStyle(1001);
  ci = TColor::GetColor("#663300");
  entry->SetLineColor(ci);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1);

  entry=leg->AddEntry("NULL",Form("Single topbar W "),"f");
  entry->SetFillColor(ctbw);/////////////////////// light blue
  entry->SetFillStyle(1001);
  ci = TColor::GetColor("#663300");
  entry->SetLineColor(ci);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1);

    entry=leg->AddEntry("NULL",Form("WJets "),"f");
    entry->SetFillColor(cwj);////////////////////////// red
    entry->SetFillStyle(1001);
    ci = TColor::GetColor("#663300");
    entry->SetLineColor(ci);
    entry->SetLineStyle(1);
    entry->SetLineWidth(1);
    entry->SetMarkerColor(1);
    entry->SetMarkerStyle(20);
    entry->SetMarkerSize(1);
    
  entry=leg->AddEntry("NULL",Form("ttbar others"),"f");
  entry->SetFillColor(ctto);////////////////////////// red
  entry->SetFillStyle(1001);
  ci = TColor::GetColor("#663300");
  entry->SetLineColor(ci);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1);  

  entry=leg->AddEntry("NULL",Form("ttbar "),"f");
  entry->SetFillColor(ctt);////////////////////////// red
  entry->SetFillStyle(1001);
  ci = TColor::GetColor("#663300");
  entry->SetLineColor(ci);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1);
        }
 

        if(leg_with_evt){
            TLegendEntry *entry=leg->AddEntry("NULL",Form("Data"),"l p");
            
            //entry->SetFillStyle(1001);
            
            //entry->SetLineColor(ci);
            entry->SetLineStyle(1);
            entry->SetLineWidth(1);
            entry->SetMarkerColor(1);
            entry->SetMarkerStyle(20);
            entry->SetMarkerSize(1);
	    entry->SetTextSize(0.025);

            
            
        entry=leg->AddEntry("NULL",Form("ZJet  %d", ndy),"f");
        ci = TColor::GetColor("#ffcc33");
        entry->SetFillColor(cdy);////////////////////////sort of yellow
        entry->SetFillStyle(1001);
        ci = TColor::GetColor("#663300");
        entry->SetLineColor(ci);
        entry->SetLineStyle(1);
        entry->SetLineWidth(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(21);
        entry->SetMarkerSize(1);
	entry->SetTextSize(0.025);

	entry=leg->AddEntry("NULL",Form("WW/WZ/ZZ  %d",nww+nwz+nzz),"f");
	entry->SetFillColor(cdb);/////////////////////// light blue
	entry->SetFillStyle(1001);
	ci = TColor::GetColor("#663300");
	entry->SetLineColor(ci);
	entry->SetLineStyle(1);
	entry->SetLineWidth(1);
	entry->SetMarkerColor(1);
	entry->SetMarkerStyle(20);
	entry->SetMarkerSize(1);
	entry->SetTextSize(0.025);


        entry=leg->AddEntry("NULL",Form("Single top W  %d", ntw),"f");
        entry->SetFillColor(ctw);/////////////////////// light blue
        entry->SetFillStyle(1001);
        ci = TColor::GetColor("#663300");
        entry->SetLineColor(ci);
        entry->SetLineStyle(1);
        entry->SetLineWidth(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(20);
        entry->SetMarkerSize(1);
	entry->SetTextSize(0.025);
        
        entry=leg->AddEntry("NULL",Form("Single topbar W  %d", ntbw),"f");
        entry->SetFillColor(ctbw);/////////////////////// light blue
        entry->SetFillStyle(1001);
        ci = TColor::GetColor("#663300");
        entry->SetLineColor(ci);
        entry->SetLineStyle(1);
        entry->SetLineWidth(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(20);
        entry->SetMarkerSize(1);
	entry->SetTextSize(0.025);
        
        entry=leg->AddEntry("NULL",Form("WJets  %d", nwj),"f");
        entry->SetFillColor(cwj);////////////////////////// red
        entry->SetFillStyle(1001);
        ci = TColor::GetColor("#663300");
        entry->SetLineColor(ci);
        entry->SetLineStyle(1);
        entry->SetLineWidth(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(20);
        entry->SetMarkerSize(1);
	entry->SetTextSize(0.025);
        
        entry=leg->AddEntry("NULL",Form("ttbar others haddsemi + out %d + %d", ntto, (evt_all-evt)),"f");
        entry->SetFillColor(ctto);////////////////////////// red
        entry->SetFillStyle(1001);
        ci = TColor::GetColor("#663300");
        entry->SetLineColor(ci);
        entry->SetLineStyle(1);
        entry->SetLineWidth(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(20);
        entry->SetMarkerSize(1);
	entry->SetTextSize(0.025);
        
        entry=leg->AddEntry("NULL",Form("ttbar  %d", evt),"f");
        entry->SetFillColor(ctt);////////////////////////// red
        entry->SetFillStyle(1001);
        ci = TColor::GetColor("#663300");
        entry->SetLineColor(ci);
        entry->SetLineStyle(1);
        entry->SetLineWidth(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(20);
        entry->SetMarkerSize(1);
	entry->SetTextSize(0.025);
        
        entry=leg->AddEntry("NULL",data,"");
	entry->SetTextSize(0.025);
        entry=leg->AddEntry("NULL",data2,"");
	entry->SetTextSize(0.025);
    }
        

  leg->Draw();
        
        c3->cd(2);
        p11_2 = (TPad*)c3->GetPad(2);
        p11_2->SetPad(0.04,0.02,0.98,0.3);
        p11_2->SetBottomMargin(0.35);
        p11_2->SetRightMargin(0.05);
        gPad->SetGridy();
      
        ratio->GetYaxis()->SetTitleSize(0.10);
        ratio->GetYaxis()->SetLabelSize(0.08);
        ratio->GetXaxis()->SetLabelSize(0.1);
        ratio->SetMarkerStyle(20);
        //ratio->SetMarkerSize(1);
        ratio->SetMinimum(0.7);
        ratio->SetMaximum(1.3);
        ratio->Draw("same");
        
        
        
        TGraphErrors *thegraph = new TGraphErrors(ratio);
        thegraph->SetFillStyle(3004);
        
        thegraph->SetFillColor(1);
        //thegraph->SetLineColor(1);
        //thegraph->SetFillStyle(1001);
        
        //thegraph->SetFillColor(chatch);
        //thegraph->SetLineColor(chatch);
        
        thegraph->Draw("e2SAME");
        
        
        TLatex *  tex = new TLatex(0.04116466,0.4,"Data/MC");
        tex->SetNDC();
        tex->SetTextSize(0.12);
        tex->SetTextFont(42);
        tex->SetTextAngle(90);
        tex->SetLineWidth(2);
        tex->Draw();

	if(ch1==3){
	  if(pl == 0 ){tex = new TLatex(0.3,0.0891969,xAxis);}
	  else if(pl == 1){tex = new TLatex(0.3,0.0891969,xAxis);}
	  else if(pl == 2){tex = new TLatex(0.3,0.0891969,xAxis);}
	  else if(pl == 3){tex = new TLatex(0.3,0.0891969,xAxis);}
	  else if(pl == 4){tex = new TLatex(0.2,0.0891969,xAxis);}
	  else if(pl == 5){tex = new TLatex(0.2,0.0891969,xAxis);}
	  else if(pl == 6){tex = new TLatex(0.2,0.0891969,xAxis);}
	  else if(pl == 7){tex = new TLatex(0.2,0.0891969,xAxis);}
	}
	else{
	  if(pl == 0 ){tex = new TLatex(0.4,0.0891969,xAxis);}
	  else if(pl == 1){tex = new TLatex(0.4,0.0891969,xAxis);}
	  else if(pl == 2){tex = new TLatex(0.4,0.0891969,xAxis);}
	  else if(pl == 3){tex = new TLatex(0.4,0.0891969,xAxis);}
	  else if(pl == 4){tex = new TLatex(0.3,0.0891969,xAxis);}
	  else if(pl == 5){tex = new TLatex(0.3,0.0891969,xAxis);}
	  else if(pl == 6){tex = new TLatex(0.35,0.0891969,xAxis);}
	  else if(pl == 7){tex = new TLatex(0.35,0.0891969,xAxis);}
	}


//        tex = new TLatex(0.5,0.0891969,xAxis);
        tex->SetNDC();
        tex->SetTextSize(0.15);
        tex->SetTextFont(42);

        tex->SetLineWidth(2);
        tex->Draw();

        
        c3->Update(); c3->SaveAs(ofile);
    
    
  
    
}
}




//gROOT->ProcessLine(".q");

}
}













