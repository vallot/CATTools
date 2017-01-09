#include <TH1F.h>
#include <TFile.h>
#include <TString.h>

#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooHistPdf.h>
#include <RooDataHist.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooSimultaneous.h>
#include <RooCategory.h>
#include <RooDataSet.h>

#include <string>
#include <iostream>
#include <map>

using namespace std;

//#define WSCALE 1

void fit()
{
  const char* srcdir = "hists";
  //const char* varName = "bjets_n";
  //const char* varName = "event_mT";
  //const char* varName = "event_mlj";
  const char* varName = "event_m3";

  const double lumi = 36.4*1e3; // # in invpb
  const double xs0 = 831.76;

  const std::map<std::string, TFile*> files = {
    {"rdata", new TFile(Form("%s/data.root" , srcdir))},
    {"wjets", new TFile(Form("%s/WJet.root" , srcdir))},
    {"mcbkg", new TFile(Form("%s/bkg.root"  , srcdir))},
    {"ttbar", new TFile(Form("%s/ttbar.root", srcdir))},
  };

  std::map<std::string, TH1*> hs_el, hs_mu;
  std::map<std::string, double> n0_el, n0_mu;
  for ( auto key=files.begin(); key!=files.end(); ++key ) {
    const auto name = key->first;
    const auto file = key->second;

    TH1* h_el = (TH1*)file->Get(Form("el/%s", varName));
    TH1* h_mu = (TH1*)file->Get(Form("mu/%s", varName));
    hs_el[name] = h_el;
    hs_mu[name] = h_mu;
    n0_el[name] = h_el->Integral();
    n0_mu[name] = h_mu->Integral();
  }

  // Normalization constants, glbal parameters and PoI
  RooRealVar var_lumi("lumi", "luminosity", lumi, "fb^{-1}");
  RooRealVar var_xs0("xs0", "cross section", xs0, "pb");// ## Default xs at the generation
  RooRealVar var_xs("xs", "cross section", xs0, 0.5*xs0, 1.5*xs0, "pb");

  auto hTemp = hs_el["rdata"];
  auto tempAxis = hTemp->GetXaxis();
  RooBinning bin_x(tempAxis->GetXmin(), tempAxis->GetXmax());
  for ( int i=1; i<=hTemp->GetNbinsX()+1; ++i ) {
    bin_x.addBoundary(tempAxis->GetBinLowEdge(i));
  }
  RooRealVar var_x(varName, varName, tempAxis->GetXmin(), tempAxis->GetXmax());
  var_x.setBinning(bin_x);

  // Build RooFit histograms
  //RooDataHist dh_el_rdata("dh_el_rdata", "dh_el_rdata", RooArgList(var), hs_el["rdata"]);
  RooDataHist dh_el_wjets("dh_el_wjets", "dh_el_wjets", RooArgList(var_x), hs_el["wjets"]);
  RooDataHist dh_el_mcbkg("dh_el_mcbkg", "dh_el_mcbkg", RooArgList(var_x), hs_el["mcbkg"]);
  RooDataHist dh_el_ttbar("dh_el_ttbar", "dh_el_ttbar", RooArgList(var_x), hs_el["ttbar"]);

  //RooDataHist dh_mu_rdata("dh_mu_rdata", "dh_mu_rdata", RooArgList(var), hs_mu["rdata"]);
  RooDataHist dh_mu_wjets("dh_mu_wjets", "dh_mu_wjets", RooArgList(var_x), hs_mu["wjets"]);
  RooDataHist dh_mu_mcbkg("dh_mu_mcbkg", "dh_mu_mcbkg", RooArgList(var_x), hs_mu["mcbkg"]);
  RooDataHist dh_mu_ttbar("dh_mu_ttbar", "dh_mu_ttbar", RooArgList(var_x), hs_mu["ttbar"]);

  RooRealVar var_weight("var_weight", "var_weight", 0, std::numeric_limits<double>::min(),
                                                       std::numeric_limits<double>::max());
  RooCategory cat("cat", "cat");
  cat.defineType("mu");
  cat.defineType("el");

  RooDataSet dh_comb("dh_comb", "dh_comb", RooArgSet(var_x, cat, var_weight), "var_weight");
  
  for ( int i=1, n=hs_el["rdata"]->GetNbinsX(); i<=n; ++i ) {
    const double x = hs_el["rdata"]->GetXaxis()->GetBinCenter(i);
    const double n_el = hs_el["rdata"]->GetBinContent(i);
    const double n_mu = hs_mu["rdata"]->GetBinContent(i);
    var_x.setVal(x);
    var_weight.setVal(n_el);
    cat.setLabel("el");
    dh_comb.add(RooArgSet(var_x, cat, var_weight), n_el);
    var_x.setVal(x);
    var_weight.setVal(n_mu);
    cat.setLabel("mu");
    dh_comb.add(RooArgSet(var_x, cat, var_weight), n_mu);
  }

  // Build RooHistPDFs
  RooHistPdf rp_el_wjets("rp_el_wjets", "el_wjets", RooArgSet(var_x), dh_el_wjets);
  RooHistPdf rp_el_mcbkg("rp_el_mcbkg", "el_mcbkg", RooArgSet(var_x), dh_el_mcbkg);
  RooHistPdf rp_el_ttbar("rp_el_ttbar", "el_ttbar", RooArgSet(var_x), dh_el_ttbar);

  RooHistPdf rp_mu_wjets("rp_mu_wjets", "mu_wjets", RooArgSet(var_x), dh_mu_wjets);
  RooHistPdf rp_mu_mcbkg("rp_mu_mcbkg", "mu_mcbkg", RooArgSet(var_x), dh_mu_mcbkg);
  RooHistPdf rp_mu_ttbar("rp_mu_ttbar", "mu_ttbar", RooArgSet(var_x), dh_mu_ttbar);

  // Normalization scales
  RooConstVar var_norm_el_wjets("norm_el_wjets", "el_wjets", n0_el["wjets"]);
  RooConstVar var_norm_el_mcbkg("norm_el_mcbkg", "el_mcbkg", n0_el["mcbkg"]);
  RooConstVar var_norm_el_ttbar("norm_el_ttbar", "el_ttbar", n0_el["ttbar"]);

  RooConstVar var_norm_mu_wjets("norm_mu_wjets", "mu_wjets", n0_mu["wjets"]);
  RooConstVar var_norm_mu_mcbkg("norm_mu_mcbkg", "mu_mcbkg", n0_mu["mcbkg"]);
  RooConstVar var_norm_mu_ttbar("norm_mu_ttbar", "mu_ttbar", n0_mu["ttbar"]);

  // Number of ttbar, express in terms of xs
  // Histograms are already normalized by (xs at generation), at lumi=1/pb
  RooFormulaVar var_n_mu_ttbar("n_mu_ttbar", "xs/xs0*lumi*norm_mu_ttbar",
                               RooArgList(var_xs, var_xs0, var_lumi, var_norm_mu_ttbar));
  RooFormulaVar var_n_el_ttbar("n_el_ttbar", "xs/xs0*lumi*norm_el_ttbar",
                               RooArgList(var_xs, var_xs0, var_lumi, var_norm_el_ttbar));

  // Number of... others
#ifndef WSCALE
  RooRealVar var_scale_wjets("scale_wjets", "wjets scaling", 1);
#else
  RooRealVar var_scale_wjets("scale_wjets", "wjets scaling", 1, 0, 2);
#endif
  RooFormulaVar var_n_mu_wjets("n_mu_wjets", "lumi*norm_mu_wjets*scale_wjets",
                               RooArgList(var_lumi, var_norm_mu_wjets, var_scale_wjets));
  RooFormulaVar var_n_mu_mcbkg("n_mu_mcbkg", "lumi*norm_mu_mcbkg", RooArgList(var_lumi, var_norm_mu_mcbkg));

  RooFormulaVar var_n_el_wjets("n_el_wjets", "lumi*norm_el_wjets*scale_wjets",
                               RooArgList(var_lumi, var_norm_el_wjets, var_scale_wjets));
  RooFormulaVar var_n_el_mcbkg("n_el_mcbkg", "lumi*norm_el_mcbkg", RooArgList(var_lumi, var_norm_el_mcbkg));


  // Define sumPDF
  RooAddPdf rf_mu_sumPDF("mu_sumPDF", "sumPDF", RooArgList(rp_mu_mcbkg, rp_mu_wjets, rp_mu_ttbar),
                         RooArgList(var_n_mu_mcbkg, var_n_mu_wjets, var_n_mu_ttbar));
  RooAddPdf rf_el_sumPDF("el_sumPDF", "sumPDF", RooArgList(rp_el_mcbkg, rp_el_wjets, rp_el_ttbar),
                         RooArgList(var_n_el_mcbkg, var_n_el_wjets, var_n_el_ttbar));

  // Simultaneous fit
  RooSimultaneous rf_simPDF("simPDF", "simultaenous PDF", cat);
  rf_simPDF.addPdf(rf_mu_sumPDF, "mu");
  rf_simPDF.addPdf(rf_el_sumPDF, "el");
  rf_simPDF.fitTo(dh_comb);

  // Fitting is done. Continue to draw them
  auto frame_el = var_x.frame();
  dh_comb.plotOn(frame_el, RooFit::Cut("cat==cat::el"), RooFit::Binning(bin_x));
  rf_simPDF.plotOn(frame_el, RooFit::Slice(cat, "el"), RooFit::ProjWData(dh_comb),
                             RooFit::DrawOption("F"), RooFit::FillColor(kBlue));
  rf_simPDF.plotOn(frame_el, RooFit::Slice(cat, "el"), RooFit::ProjWData(dh_comb), RooFit::Components("rp_el_ttbar,rp_el_wjets"),
                             RooFit::DrawOption("F"), RooFit::FillColor(kGreen));
  rf_simPDF.plotOn(frame_el, RooFit::Slice(cat, "el"), RooFit::ProjWData(dh_comb), RooFit::Components("rp_el_ttbar"),
                             RooFit::DrawOption("F"), RooFit::FillColor(kRed));
  //rf_simPDF.plotOn(frame_el, RooFit::Slice(cat, "el"), RooFit::ProjWData(dh_comb), RooFit::Components("rp_el_wjets"));
  //rf_simPDF.plotOn(frame_el, RooFit::Slice(cat, "el"), RooFit::ProjWData(dh_comb), RooFit::Components("rp_el_mcbkg"));
  dh_comb.plotOn(frame_el, RooFit::Cut("cat==cat::el"), RooFit::Binning(bin_x), RooFit::DataError(RooAbsData::Poisson));

  auto frame_mu = var_x.frame();
  dh_comb.plotOn(frame_mu, RooFit::Cut("cat==cat::mu"), RooFit::Binning(bin_x));
  rf_simPDF.plotOn(frame_mu, RooFit::Slice(cat, "mu"), RooFit::ProjWData(dh_comb),
                             RooFit::DrawOption("F"), RooFit::FillColor(kBlue));
  rf_simPDF.plotOn(frame_mu, RooFit::Slice(cat, "mu"), RooFit::ProjWData(dh_comb), RooFit::Components("rp_mu_ttbar,rp_mu_wjets"),
                             RooFit::DrawOption("F"), RooFit::FillColor(kGreen));
  rf_simPDF.plotOn(frame_mu, RooFit::Slice(cat, "mu"), RooFit::ProjWData(dh_comb), RooFit::Components("rp_mu_ttbar"),
                             RooFit::DrawOption("F"), RooFit::FillColor(kRed));
  //rf_simPDF.plotOn(frame_mu, RooFit::Slice(cat, "mu"), RooFit::ProjWData(dh_comb), RooFit::Components("rp_mu_wjets"));
  //rf_simPDF.plotOn(frame_mu, RooFit::Slice(cat, "mu"), RooFit::ProjWData(dh_comb), RooFit::Components("rp_mu_mcbkg"));
  dh_comb.plotOn(frame_mu, RooFit::Cut("cat==cat::mu"), RooFit::Binning(bin_x), RooFit::DataError(RooAbsData::Poisson));

  TCanvas* c1 = new TCanvas("c1", "c1", 1000, 500);
  c1->Divide(2,1);
  c1->cd(1);
  frame_el->Draw();
  c1->cd(2);
  frame_mu->Draw();

  auto nll = rf_simPDF.createNLL(dh_comb, RooFit::Extended(true));
#ifndef WSCALE
  auto nllFrame = var_xs.frame();
  nll->plotOn(nllFrame, RooFit::ShiftToZero());
  auto nllCurve = nllFrame->getCurve();
  const double maxNLL = std::max(nllCurve->Eval(var_xs.getMin()), nllCurve->Eval(var_xs.getMax()));
  nllFrame->SetMaximum(maxNLL*1.1);
  nllFrame->SetMinimum(-maxNLL*0.1);

  TCanvas* c2 = new TCanvas("c2", "c2", 500, 500);
  nllFrame->Draw();
#else
  TCanvas* c3 = new TCanvas("c3", "c3", 500, 500);
  RooMinuit m(*nll);
  m.migrad();
  m.hesse();
  m.minos();
  auto nllFrame2 = m.contour(var_xs, var_scale_wjets, 1, 2, 3);
  double x = 0, ex = 0;
  x = var_xs.getVal();
  ex = var_xs.getError();
  nllFrame2->SetAxisRange(x-4*ex, x+4*ex, "X");
  x = var_scale_wjets.getVal();
  ex = var_scale_wjets.getError();
  nllFrame2->SetAxisRange(x-4*ex, x+4*ex, "Y");
  nllFrame2->Draw();
#endif

  cout << "\n\n*********************************\n";
  cout << "Cross section = " << var_xs.getVal() << " + " << var_xs.getErrorHi() << " - " << var_xs.getErrorLo() << endl;
  cout << "scale W+jets = " << var_scale_wjets.getVal() << " + " << var_scale_wjets.getErrorHi() << " - " << var_scale_wjets.getErrorLo() << endl;
  cout << "*********************************\n\n";

}

