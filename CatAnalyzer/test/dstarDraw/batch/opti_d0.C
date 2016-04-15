int opti_d0(float relPt_val=0.05, float dR_val =0.15 )
{

  gROOT->SetBatch(true);
  //TString inputfile("../cattools/cattree_TTJets_aMC.root");
  //TString inputfile("../cattools/cattree_WW.root");
  TString inputfile("d0.root");



  RooRealVar x("mass","d0 mass",1.6,2.2);
  RooRealVar step("step","step",0,8);
  RooRealVar dR("dR","dR between Gen and reco",0,0.5);
  RooRealVar relPt("relPt","relPt",-1,1);


  auto gaus_mean = RooRealVar("mean","mean",1.864, 1.6,2.2);
  auto gaus_sigma = RooRealVar("sigma","sigma",0,1);
  auto gaus_pdf = RooGaussian("sig","signal p.d.f",x,gaus_mean, gaus_sigma);


  auto expo_tau = RooRealVar("exp_const","exp_const",0,-100,100);
  auto exp_pdf = RooExponential("bkg","bkg p.d.f",x,expo_tau);


  auto file = TFile::Open(inputfile);
  TTree* tree = (TTree*)file->Get("nt");
  TFile* out = new TFile("mcTruth.root","RECREATE");
  TH2F* h1 = new TH2F("mcTruth_sig","Significance of mc Truth;abs(relPt); dR ",20,0,1,20,0,0.5); 
  TH2F* h2 = new TH2F("ratio_sig_bkg","sig/bkg  ;abs(relPt); dR ",20,0,1,20,0,0.5); 
  TH2F* h3 = new TH2F("mcTruth_sig2","sig/bkg*signifi  ;abs(relPt); dR ",20,0,1,20,0,0.5); 
  TH2F* h4 = new TH2F("sig_bkg_event","sig/bkg*nsig;abs(relPt); dR ",20,0,1,20,0,0.5); 

  RooDataSet data("data","dataset with d0 mass",RooArgSet(x,step,dR,relPt),RooFit::Import(*tree));

  system("touch relPt_dR.png");
  //for( float rel_pt = 0.05 ; rel_pt <1.0 ; rel_pt = rel_pt+(1/20.)) { 
  //for( float dR = 0.025 ; dR <0.5 ; dR = dR+(0.5/20.)) {
  {
      auto cut = TString::Format("step>=5 && dR<%f && abs(relPt)<%f && dR>=0.0 && abs(relPt) >=0.0",dR_val, relPt_val);
      std::cout<<cut<<std::endl;
      RooDataSet* d1 = dynamic_cast<RooDataSet*>(data.reduce(cut.Data()));
        //RooRealVar fsig("fsig","Signal fraction",0.5,0,1);

        int totEvent=tree->GetEntries();

      auto nsig = RooRealVar("nsig","nsig",0,totEvent);
      auto nbkg = RooRealVar("nbkg","nbkg",0,totEvent);

      auto model = RooAddPdf("model1","model1",RooArgList(gaus_pdf, exp_pdf),RooArgList(nsig,nbkg));

      
      auto* fitResult = model.fitTo(*d1, RooFit::Save());
      if ( fitResult != nullptr) {
        //fitResult->Print();
        std::cout<< "rel_pt : "<<relPt_val << "dR : "<<dR_val << " "<<nsig.getVal() << " "<<nsig.getError()<< "    "<<nbkg.getVal() << " "<<nbkg.getError()<<std::endl;
        if ( abs(gaus_mean.getVal() - 1.8640)>0.05 ) { h1->Fill( relPt_val, dR_val, 0) ; h1->Write() ; out->Close(); return -1;}
        if ( abs(gaus_sigma.getVal())>0.1 ) { h1->Fill( relPt_val, dR_val, 0) ; h1->Write() ; out->Close();return -1;}
        float significance = nsig.getVal() / TMath::Sqrt( nsig.getVal()+ nbkg.getVal());
        h1->Fill( relPt_val-0.025, dR_val-0.0125, significance);
        h2->Fill( relPt_val-0.025, dR_val-0.0125, nsig.getVal()/nbkg.getVal());
        h3->Fill( relPt_val-0.025, dR_val-0.0125, nsig.getVal()/nbkg.getVal()*significance);
        h4->Fill( relPt_val-0.025, dR_val-0.0125, (nsig.getVal()+nbkg.getVal())*nsig.getVal()/nbkg.getVal());
        delete fitResult;
      }

      RooPlot* d0_mass_frame = x.frame();
      d1->plotOn(d0_mass_frame);

      model.plotOn(d0_mass_frame);
      model.plotOn(d0_mass_frame, RooFit::Components(exp_pdf),RooFit::LineStyle(kDashed));
      model.plotOn(d0_mass_frame, RooFit::Components(gaus_pdf),RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));

      TString canvas_name = TString::Format("DeltaR_%0.4f_ReleativePt_%0.4f",dR_val, relPt_val);
      TCanvas* c1 = new TCanvas(canvas_name.Data(),canvas_name.Data(),600,600);
      d0_mass_frame->Draw();
      d0_mass_frame->SetTitle(canvas_name.Data());
      c1->SaveAs( "relPt_dR.png");
      delete c1;
      delete d1;
    }
    //}
 // }
  h1->Write();
  h2->Write();
  h3->Write();
  h4->Write();
  out->Close();
  return 0;
}
