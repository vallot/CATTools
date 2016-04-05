//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 22 11:57:23 2016 by ROOT version 6.02/13
// from TTree nom/nom
// found on file: cattree.root
//////////////////////////////////////////////////////////

#ifndef dstarAna_h
#define dstarAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
using namespace std;
class dstarAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nvertex;
   Int_t           step;
   Int_t           channel;
   Int_t           njet;
   Int_t           nbjet;
   Bool_t          step1;
   Bool_t          step2;
   Bool_t          step3;
   Bool_t          step4;
   Bool_t          step5;
   Bool_t          step6;
   Float_t         tri;
   Bool_t          filtered;
   Float_t         met;
   Float_t         weight;
   vector<float>   *pdfWeights;
   vector<float>   *scaleWeights;
   Float_t         topPtWeight;
   Float_t         puweight;
   Float_t         puweight_up;
   Float_t         puweight_dn;
   Float_t         genweight;
   Float_t         mueffweight;
   Float_t         mueffweight_up;
   Float_t         mueffweight_dn;
   Float_t         eleffweight;
   Float_t         eleffweight_up;
   Float_t         eleffweight_dn;
   Float_t         btagweight;
   Float_t         btagweight_up;
   Float_t         btagweight_dn;
   vector<float>   *csvweights;
   Float_t         lep1_pt;
   Float_t         lep1_eta;
   Float_t         lep1_phi;
   Float_t         lep2_pt;
   Float_t         lep2_eta;
   Float_t         lep2_phi;
   Float_t         ll_pt;
   Float_t         ll_eta;
   Float_t         ll_phi;
   Float_t         ll_m;
   Float_t         partontop1_pt;
   Float_t         partontop1_eta;
   Float_t         partontop1_rapi;
   Float_t         partontop1_m;
   Float_t         partontop2_pt;
   Float_t         partontop2_eta;
   Float_t         partontop2_rapi;
   Float_t         partontop2_m;
   Float_t         partonttbar_pt;
   Float_t         partonttbar_eta;
   Float_t         partonttbar_dphi;
   Float_t         partonttbar_rapi;
   Float_t         partonttbar_m;
   Int_t           parton_channel;
   Int_t           parton_mode1;
   Float_t         partonlep1_pt;
   Float_t         partonlep1_eta;
   Float_t         partonlep2_pt;
   Float_t         partonlep2_eta;
   Float_t         partonjet1_pt;
   Float_t         partonjet1_eta;
   Float_t         partonjet2_pt;
   Float_t         partonjet2_eta;
   Int_t           parton_mode2;
   Bool_t          partonInPhase;
   Bool_t          partonInPhaseLep;
   Bool_t          partonInPhaseJet;
   Float_t         gentop1_pt;
   Float_t         gentop1_eta;
   Float_t         gentop1_rapi;
   Float_t         gentop1_m;
   Float_t         gentop2_pt;
   Float_t         gentop2_eta;
   Float_t         gentop2_rapi;
   Float_t         gentop2_m;
   Float_t         genttbar_pt;
   Float_t         genttbar_eta;
   Float_t         genttbar_dphi;
   Float_t         genttbar_rapi;
   Float_t         genttbar_m;
   Int_t           pseudoTop_channel;
   Float_t         genlep1_pt;
   Float_t         genlep1_eta;
   Float_t         genlep2_pt;
   Float_t         genlep2_eta;
   Float_t         genjet1_pt;
   Float_t         genjet1_eta;
   Float_t         genjet2_pt;
   Float_t         genjet2_eta;
   Bool_t          pseudoInPhase;
   Float_t         jet1_pt;
   Float_t         jet2_pt;
   Float_t         jet1_eta;
   Float_t         jet2_eta;
   Float_t         jet1_CSVInclV2;
   Float_t         jet2_CSVInclV2;
   Float_t         top1_pt;
   Float_t         top1_eta;
   Float_t         top1_phi;
   Float_t         top1_rapi;
   Float_t         top1_m;
   Float_t         top2_pt;
   Float_t         top2_eta;
   Float_t         top2_phi;
   Float_t         top2_rapi;
   Float_t         top2_m;
   Float_t         ttbar_pt;
   Float_t         ttbar_eta;
   Float_t         ttbar_dphi;
   Float_t         ttbar_rapi;
   Float_t         ttbar_m;
   Int_t           is3lep;
   vector<float>   *d0_pt;
   vector<float>   *d0_eta;
   vector<float>   *d0_phi;
   vector<float>   *d0_m;
   vector<bool>    *d0_true;
   vector<bool>    *d0_fit;
   vector<float>   *d0_L3D;
   vector<float>   *d0_LXY;
   vector<float>   *d0_dRTrue;
   vector<float>   *d0_relPtTrue;
   vector<float>   *d0_dau1_pt;
   vector<float>   *d0_dau1_eta;
   vector<float>   *d0_dau1_phi;
   vector<float>   *d0_dau1_m;
   vector<float>   *d0_dau1_q;
   vector<float>   *d0_dau2_pt;
   vector<float>   *d0_dau2_eta;
   vector<float>   *d0_dau2_phi;
   vector<float>   *d0_dau2_m;
   vector<float>   *d0_dau2_q;
   vector<float>   *dstar_pt;
   vector<float>   *dstar_eta;
   vector<float>   *dstar_phi;
   vector<float>   *dstar_m;
   vector<bool>    *dstar_true;
   vector<bool>    *dstar_fit;
   vector<float>   *dstar_L3D;
   vector<float>   *dstar_LXY;
   vector<float>   *dstar_dRTrue;
   vector<float>   *dstar_relPtTrue;
   vector<float>   *dstar_dau1_pt;
   vector<float>   *dstar_dau1_eta;
   vector<float>   *dstar_dau1_phi;
   vector<float>   *dstar_dau1_m;
   vector<float>   *dstar_dau1_q;
   vector<float>   *dstar_dau2_pt;
   vector<float>   *dstar_dau2_eta;
   vector<float>   *dstar_dau2_phi;
   vector<float>   *dstar_dau2_m;
   vector<float>   *dstar_dau2_q;
   vector<float>   *dstar_dau3_pt;
   vector<float>   *dstar_dau3_eta;
   vector<float>   *dstar_dau3_phi;
   vector<float>   *dstar_dau3_m;
   vector<float>   *dstar_dau3_q;

   // List of branches
   TBranch        *b_nvertex;   //!
   TBranch        *b_step;   //!
   TBranch        *b_channel;   //!
   TBranch        *b_njet;   //!
   TBranch        *b_nbjet;   //!
   TBranch        *b_step1;   //!
   TBranch        *b_step2;   //!
   TBranch        *b_step3;   //!
   TBranch        *b_step4;   //!
   TBranch        *b_step5;   //!
   TBranch        *b_step6;   //!
   TBranch        *b_tri;   //!
   TBranch        *b_filtered;   //!
   TBranch        *b_met;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_pdfWeights;   //!
   TBranch        *b_scaleWeights;   //!
   TBranch        *b_topPtWeight;   //!
   TBranch        *b_puweight;   //!
   TBranch        *b_puweight_up;   //!
   TBranch        *b_puweight_dn;   //!
   TBranch        *b_genweight;   //!
   TBranch        *b_mueffweight;   //!
   TBranch        *b_mueffweight_up;   //!
   TBranch        *b_mueffweight_dn;   //!
   TBranch        *b_eleffweight;   //!
   TBranch        *b_eleffweight_up;   //!
   TBranch        *b_eleffweight_dn;   //!
   TBranch        *b_btagweight;   //!
   TBranch        *b_btagweight_up;   //!
   TBranch        *b_btagweight_dn;   //!
   TBranch        *b_csvweights;   //!
   TBranch        *b_lep1_pt;   //!
   TBranch        *b_lep1_eta;   //!
   TBranch        *b_lep1_phi;   //!
   TBranch        *b_lep2_pt;   //!
   TBranch        *b_lep2_eta;   //!
   TBranch        *b_lep2_phi;   //!
   TBranch        *b_ll_pt;   //!
   TBranch        *b_ll_eta;   //!
   TBranch        *b_ll_phi;   //!
   TBranch        *b_ll_m;   //!
   TBranch        *b_partontop1_pt;   //!
   TBranch        *b_partontop1_eta;   //!
   TBranch        *b_partontop1_rapi;   //!
   TBranch        *b_partontop1_m;   //!
   TBranch        *b_partontop2_pt;   //!
   TBranch        *b_partontop2_eta;   //!
   TBranch        *b_partontop2_rapi;   //!
   TBranch        *b_partontop2_m;   //!
   TBranch        *b_partonttbar_pt;   //!
   TBranch        *b_partonttbar_eta;   //!
   TBranch        *b_partonttbar_dphi;   //!
   TBranch        *b_partonttbar_rapi;   //!
   TBranch        *b_partonttbar_m;   //!
   TBranch        *b_parton_channel;   //!
   TBranch        *b_parton_mode1;   //!
   TBranch        *b_partonlep1_pt;   //!
   TBranch        *b_partonlep1_eta;   //!
   TBranch        *b_partonlep2_pt;   //!
   TBranch        *b_partonlep2_eta;   //!
   TBranch        *b_partonjet1_pt;   //!
   TBranch        *b_partonjet1_eta;   //!
   TBranch        *b_partonjet2_pt;   //!
   TBranch        *b_partonjet2_eta;   //!
   TBranch        *b_parton_mode2;   //!
   TBranch        *b_partonInPhase;   //!
   TBranch        *b_partonInPhaseLep;   //!
   TBranch        *b_partonInPhaseJet;   //!
   TBranch        *b_gentop1_pt;   //!
   TBranch        *b_gentop1_eta;   //!
   TBranch        *b_gentop1_rapi;   //!
   TBranch        *b_gentop1_m;   //!
   TBranch        *b_gentop2_pt;   //!
   TBranch        *b_gentop2_eta;   //!
   TBranch        *b_gentop2_rapi;   //!
   TBranch        *b_gentop2_m;   //!
   TBranch        *b_genttbar_pt;   //!
   TBranch        *b_genttbar_eta;   //!
   TBranch        *b_genttbar_dphi;   //!
   TBranch        *b_genttbar_rapi;   //!
   TBranch        *b_genttbar_m;   //!
   TBranch        *b_pseudoTop_channel;   //!
   TBranch        *b_genlep1_pt;   //!
   TBranch        *b_genlep1_eta;   //!
   TBranch        *b_genlep2_pt;   //!
   TBranch        *b_genlep2_eta;   //!
   TBranch        *b_genjet1_pt;   //!
   TBranch        *b_genjet1_eta;   //!
   TBranch        *b_genjet2_pt;   //!
   TBranch        *b_genjet2_eta;   //!
   TBranch        *b_pseudoInPhase;   //!
   TBranch        *b_jet1_pt;   //!
   TBranch        *b_jet2_pt;   //!
   TBranch        *b_jet1_eta;   //!
   TBranch        *b_jet2_eta;   //!
   TBranch        *b_jet1_CSVInclV2;   //!
   TBranch        *b_jet2_CSVInclV2;   //!
   TBranch        *b_top1_pt;   //!
   TBranch        *b_top1_eta;   //!
   TBranch        *b_top1_phi;   //!
   TBranch        *b_top1_rapi;   //!
   TBranch        *b_top1_m;   //!
   TBranch        *b_top2_pt;   //!
   TBranch        *b_top2_eta;   //!
   TBranch        *b_top2_phi;   //!
   TBranch        *b_top2_rapi;   //!
   TBranch        *b_top2_m;   //!
   TBranch        *b_ttbar_pt;   //!
   TBranch        *b_ttbar_eta;   //!
   TBranch        *b_ttbar_dphi;   //!
   TBranch        *b_ttbar_rapi;   //!
   TBranch        *b_ttbar_m;   //!
   TBranch        *b_is3lep;   //!
   TBranch        *b_d0_pt;   //!
   TBranch        *b_d0_eta;   //!
   TBranch        *b_d0_phi;   //!
   TBranch        *b_d0_m;   //!
   TBranch        *b_d0_true;   //!
   TBranch        *b_d0_fit;   //!
   TBranch        *b_d0_L3D;   //!
   TBranch        *b_d0_LXY;   //!
   TBranch        *b_d0_dRTrue;   //!
   TBranch        *b_d0_relPtTrue;   //!
   TBranch        *b_d0_dau1_pt;   //!
   TBranch        *b_d0_dau1_eta;   //!
   TBranch        *b_d0_dau1_phi;   //!
   TBranch        *b_d0_dau1_m;   //!
   TBranch        *b_d0_dau1_q;   //!
   TBranch        *b_d0_dau2_pt;   //!
   TBranch        *b_d0_dau2_eta;   //!
   TBranch        *b_d0_dau2_phi;   //!
   TBranch        *b_d0_dau2_m;   //!
   TBranch        *b_d0_dau2_q;   //!
   TBranch        *b_dstar_pt;   //!
   TBranch        *b_dstar_eta;   //!
   TBranch        *b_dstar_phi;   //!
   TBranch        *b_dstar_m;   //!
   TBranch        *b_dstar_true;   //!
   TBranch        *b_dstar_fit;   //!
   TBranch        *b_dstar_L3D;   //!
   TBranch        *b_dstar_LXY;   //!
   TBranch        *b_dstar_dRTrue;   //!
   TBranch        *b_dstar_relPtTrue;   //!
   TBranch        *b_dstar_dau1_pt;   //!
   TBranch        *b_dstar_dau1_eta;   //!
   TBranch        *b_dstar_dau1_phi;   //!
   TBranch        *b_dstar_dau1_m;   //!
   TBranch        *b_dstar_dau1_q;   //!
   TBranch        *b_dstar_dau2_pt;   //!
   TBranch        *b_dstar_dau2_eta;   //!
   TBranch        *b_dstar_dau2_phi;   //!
   TBranch        *b_dstar_dau2_m;   //!
   TBranch        *b_dstar_dau2_q;   //!
   TBranch        *b_dstar_dau3_pt;   //!
   TBranch        *b_dstar_dau3_eta;   //!
   TBranch        *b_dstar_dau3_phi;   //!
   TBranch        *b_dstar_dau3_m;   //!
   TBranch        *b_dstar_dau3_q;   //!

   dstarAna(TTree *tree=0);
   virtual ~dstarAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef dstarAna_cxx
dstarAna::dstarAna(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("cattree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("cattree.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("cattree.root:/cattree");
      dir->GetObject("nom",tree);

   }
   Init(tree);
}

dstarAna::~dstarAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t dstarAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t dstarAna::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void dstarAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   pdfWeights = 0;
   scaleWeights = 0;
   csvweights = 0;
   d0_pt = 0;
   d0_eta = 0;
   d0_phi = 0;
   d0_m = 0;
   d0_true = 0;
   d0_fit = 0;
   d0_L3D = 0;
   d0_LXY = 0;
   d0_dRTrue = 0;
   d0_relPtTrue = 0;
   d0_dau1_pt = 0;
   d0_dau1_eta = 0;
   d0_dau1_phi = 0;
   d0_dau1_m = 0;
   d0_dau1_q = 0;
   d0_dau2_pt = 0;
   d0_dau2_eta = 0;
   d0_dau2_phi = 0;
   d0_dau2_m = 0;
   d0_dau2_q = 0;
   dstar_pt = 0;
   dstar_eta = 0;
   dstar_phi = 0;
   dstar_m = 0;
   dstar_true = 0;
   dstar_fit = 0;
   dstar_L3D = 0;
   dstar_LXY = 0;
   dstar_dRTrue = 0;
   dstar_relPtTrue = 0;
   dstar_dau1_pt = 0;
   dstar_dau1_eta = 0;
   dstar_dau1_phi = 0;
   dstar_dau1_m = 0;
   dstar_dau1_q = 0;
   dstar_dau2_pt = 0;
   dstar_dau2_eta = 0;
   dstar_dau2_phi = 0;
   dstar_dau2_m = 0;
   dstar_dau2_q = 0;
   dstar_dau3_pt = 0;
   dstar_dau3_eta = 0;
   dstar_dau3_phi = 0;
   dstar_dau3_m = 0;
   dstar_dau3_q = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nvertex", &nvertex, &b_nvertex);
   fChain->SetBranchAddress("step", &step, &b_step);
   fChain->SetBranchAddress("channel", &channel, &b_channel);
   fChain->SetBranchAddress("njet", &njet, &b_njet);
   fChain->SetBranchAddress("nbjet", &nbjet, &b_nbjet);
   fChain->SetBranchAddress("step1", &step1, &b_step1);
   fChain->SetBranchAddress("step2", &step2, &b_step2);
   fChain->SetBranchAddress("step3", &step3, &b_step3);
   fChain->SetBranchAddress("step4", &step4, &b_step4);
   fChain->SetBranchAddress("step5", &step5, &b_step5);
   fChain->SetBranchAddress("step6", &step6, &b_step6);
   fChain->SetBranchAddress("tri", &tri, &b_tri);
   fChain->SetBranchAddress("filtered", &filtered, &b_filtered);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("pdfWeights", &pdfWeights, &b_pdfWeights);
   fChain->SetBranchAddress("scaleWeights", &scaleWeights, &b_scaleWeights);
   fChain->SetBranchAddress("topPtWeight", &topPtWeight, &b_topPtWeight);
   fChain->SetBranchAddress("puweight", &puweight, &b_puweight);
   fChain->SetBranchAddress("puweight_up", &puweight_up, &b_puweight_up);
   fChain->SetBranchAddress("puweight_dn", &puweight_dn, &b_puweight_dn);
   fChain->SetBranchAddress("genweight", &genweight, &b_genweight);
   fChain->SetBranchAddress("mueffweight", &mueffweight, &b_mueffweight);
   fChain->SetBranchAddress("mueffweight_up", &mueffweight_up, &b_mueffweight_up);
   fChain->SetBranchAddress("mueffweight_dn", &mueffweight_dn, &b_mueffweight_dn);
   fChain->SetBranchAddress("eleffweight", &eleffweight, &b_eleffweight);
   fChain->SetBranchAddress("eleffweight_up", &eleffweight_up, &b_eleffweight_up);
   fChain->SetBranchAddress("eleffweight_dn", &eleffweight_dn, &b_eleffweight_dn);
   fChain->SetBranchAddress("btagweight", &btagweight, &b_btagweight);
   fChain->SetBranchAddress("btagweight_up", &btagweight_up, &b_btagweight_up);
   fChain->SetBranchAddress("btagweight_dn", &btagweight_dn, &b_btagweight_dn);
   fChain->SetBranchAddress("csvweights", &csvweights, &b_csvweights);
   fChain->SetBranchAddress("lep1_pt", &lep1_pt, &b_lep1_pt);
   fChain->SetBranchAddress("lep1_eta", &lep1_eta, &b_lep1_eta);
   fChain->SetBranchAddress("lep1_phi", &lep1_phi, &b_lep1_phi);
   fChain->SetBranchAddress("lep2_pt", &lep2_pt, &b_lep2_pt);
   fChain->SetBranchAddress("lep2_eta", &lep2_eta, &b_lep2_eta);
   fChain->SetBranchAddress("lep2_phi", &lep2_phi, &b_lep2_phi);
   fChain->SetBranchAddress("ll_pt", &ll_pt, &b_ll_pt);
   fChain->SetBranchAddress("ll_eta", &ll_eta, &b_ll_eta);
   fChain->SetBranchAddress("ll_phi", &ll_phi, &b_ll_phi);
   fChain->SetBranchAddress("ll_m", &ll_m, &b_ll_m);
   fChain->SetBranchAddress("partontop1_pt", &partontop1_pt, &b_partontop1_pt);
   fChain->SetBranchAddress("partontop1_eta", &partontop1_eta, &b_partontop1_eta);
   fChain->SetBranchAddress("partontop1_rapi", &partontop1_rapi, &b_partontop1_rapi);
   fChain->SetBranchAddress("partontop1_m", &partontop1_m, &b_partontop1_m);
   fChain->SetBranchAddress("partontop2_pt", &partontop2_pt, &b_partontop2_pt);
   fChain->SetBranchAddress("partontop2_eta", &partontop2_eta, &b_partontop2_eta);
   fChain->SetBranchAddress("partontop2_rapi", &partontop2_rapi, &b_partontop2_rapi);
   fChain->SetBranchAddress("partontop2_m", &partontop2_m, &b_partontop2_m);
   fChain->SetBranchAddress("partonttbar_pt", &partonttbar_pt, &b_partonttbar_pt);
   fChain->SetBranchAddress("partonttbar_eta", &partonttbar_eta, &b_partonttbar_eta);
   fChain->SetBranchAddress("partonttbar_dphi", &partonttbar_dphi, &b_partonttbar_dphi);
   fChain->SetBranchAddress("partonttbar_rapi", &partonttbar_rapi, &b_partonttbar_rapi);
   fChain->SetBranchAddress("partonttbar_m", &partonttbar_m, &b_partonttbar_m);
   fChain->SetBranchAddress("parton_channel", &parton_channel, &b_parton_channel);
   fChain->SetBranchAddress("parton_mode1", &parton_mode1, &b_parton_mode1);
   fChain->SetBranchAddress("partonlep1_pt", &partonlep1_pt, &b_partonlep1_pt);
   fChain->SetBranchAddress("partonlep1_eta", &partonlep1_eta, &b_partonlep1_eta);
   fChain->SetBranchAddress("partonlep2_pt", &partonlep2_pt, &b_partonlep2_pt);
   fChain->SetBranchAddress("partonlep2_eta", &partonlep2_eta, &b_partonlep2_eta);
   fChain->SetBranchAddress("partonjet1_pt", &partonjet1_pt, &b_partonjet1_pt);
   fChain->SetBranchAddress("partonjet1_eta", &partonjet1_eta, &b_partonjet1_eta);
   fChain->SetBranchAddress("partonjet2_pt", &partonjet2_pt, &b_partonjet2_pt);
   fChain->SetBranchAddress("partonjet2_eta", &partonjet2_eta, &b_partonjet2_eta);
   fChain->SetBranchAddress("parton_mode2", &parton_mode2, &b_parton_mode2);
   fChain->SetBranchAddress("partonInPhase", &partonInPhase, &b_partonInPhase);
   fChain->SetBranchAddress("partonInPhaseLep", &partonInPhaseLep, &b_partonInPhaseLep);
   fChain->SetBranchAddress("partonInPhaseJet", &partonInPhaseJet, &b_partonInPhaseJet);
   fChain->SetBranchAddress("gentop1_pt", &gentop1_pt, &b_gentop1_pt);
   fChain->SetBranchAddress("gentop1_eta", &gentop1_eta, &b_gentop1_eta);
   fChain->SetBranchAddress("gentop1_rapi", &gentop1_rapi, &b_gentop1_rapi);
   fChain->SetBranchAddress("gentop1_m", &gentop1_m, &b_gentop1_m);
   fChain->SetBranchAddress("gentop2_pt", &gentop2_pt, &b_gentop2_pt);
   fChain->SetBranchAddress("gentop2_eta", &gentop2_eta, &b_gentop2_eta);
   fChain->SetBranchAddress("gentop2_rapi", &gentop2_rapi, &b_gentop2_rapi);
   fChain->SetBranchAddress("gentop2_m", &gentop2_m, &b_gentop2_m);
   fChain->SetBranchAddress("genttbar_pt", &genttbar_pt, &b_genttbar_pt);
   fChain->SetBranchAddress("genttbar_eta", &genttbar_eta, &b_genttbar_eta);
   fChain->SetBranchAddress("genttbar_dphi", &genttbar_dphi, &b_genttbar_dphi);
   fChain->SetBranchAddress("genttbar_rapi", &genttbar_rapi, &b_genttbar_rapi);
   fChain->SetBranchAddress("genttbar_m", &genttbar_m, &b_genttbar_m);
   fChain->SetBranchAddress("pseudoTop_channel", &pseudoTop_channel, &b_pseudoTop_channel);
   fChain->SetBranchAddress("genlep1_pt", &genlep1_pt, &b_genlep1_pt);
   fChain->SetBranchAddress("genlep1_eta", &genlep1_eta, &b_genlep1_eta);
   fChain->SetBranchAddress("genlep2_pt", &genlep2_pt, &b_genlep2_pt);
   fChain->SetBranchAddress("genlep2_eta", &genlep2_eta, &b_genlep2_eta);
   fChain->SetBranchAddress("genjet1_pt", &genjet1_pt, &b_genjet1_pt);
   fChain->SetBranchAddress("genjet1_eta", &genjet1_eta, &b_genjet1_eta);
   fChain->SetBranchAddress("genjet2_pt", &genjet2_pt, &b_genjet2_pt);
   fChain->SetBranchAddress("genjet2_eta", &genjet2_eta, &b_genjet2_eta);
   fChain->SetBranchAddress("pseudoInPhase", &pseudoInPhase, &b_pseudoInPhase);
   fChain->SetBranchAddress("jet1_pt", &jet1_pt, &b_jet1_pt);
   fChain->SetBranchAddress("jet2_pt", &jet2_pt, &b_jet2_pt);
   fChain->SetBranchAddress("jet1_eta", &jet1_eta, &b_jet1_eta);
   fChain->SetBranchAddress("jet2_eta", &jet2_eta, &b_jet2_eta);
   fChain->SetBranchAddress("jet1_CSVInclV2", &jet1_CSVInclV2, &b_jet1_CSVInclV2);
   fChain->SetBranchAddress("jet2_CSVInclV2", &jet2_CSVInclV2, &b_jet2_CSVInclV2);
   fChain->SetBranchAddress("top1_pt", &top1_pt, &b_top1_pt);
   fChain->SetBranchAddress("top1_eta", &top1_eta, &b_top1_eta);
   fChain->SetBranchAddress("top1_phi", &top1_phi, &b_top1_phi);
   fChain->SetBranchAddress("top1_rapi", &top1_rapi, &b_top1_rapi);
   fChain->SetBranchAddress("top1_m", &top1_m, &b_top1_m);
   fChain->SetBranchAddress("top2_pt", &top2_pt, &b_top2_pt);
   fChain->SetBranchAddress("top2_eta", &top2_eta, &b_top2_eta);
   fChain->SetBranchAddress("top2_phi", &top2_phi, &b_top2_phi);
   fChain->SetBranchAddress("top2_rapi", &top2_rapi, &b_top2_rapi);
   fChain->SetBranchAddress("top2_m", &top2_m, &b_top2_m);
   fChain->SetBranchAddress("ttbar_pt", &ttbar_pt, &b_ttbar_pt);
   fChain->SetBranchAddress("ttbar_eta", &ttbar_eta, &b_ttbar_eta);
   fChain->SetBranchAddress("ttbar_dphi", &ttbar_dphi, &b_ttbar_dphi);
   fChain->SetBranchAddress("ttbar_rapi", &ttbar_rapi, &b_ttbar_rapi);
   fChain->SetBranchAddress("ttbar_m", &ttbar_m, &b_ttbar_m);
   fChain->SetBranchAddress("is3lep", &is3lep, &b_is3lep);
   fChain->SetBranchAddress("d0_pt", &d0_pt, &b_d0_pt);
   fChain->SetBranchAddress("d0_eta", &d0_eta, &b_d0_eta);
   fChain->SetBranchAddress("d0_phi", &d0_phi, &b_d0_phi);
   fChain->SetBranchAddress("d0_m", &d0_m, &b_d0_m);
   fChain->SetBranchAddress("d0_true", &d0_true, &b_d0_true);
   fChain->SetBranchAddress("d0_fit", &d0_fit, &b_d0_fit);
   fChain->SetBranchAddress("d0_L3D", &d0_L3D, &b_d0_L3D);
   fChain->SetBranchAddress("d0_LXY", &d0_LXY, &b_d0_LXY);
   fChain->SetBranchAddress("d0_dRTrue", &d0_dRTrue, &b_d0_dRTrue);
   fChain->SetBranchAddress("d0_relPtTrue", &d0_relPtTrue, &b_d0_relPtTrue);
   fChain->SetBranchAddress("d0_dau1_pt", &d0_dau1_pt, &b_d0_dau1_pt);
   fChain->SetBranchAddress("d0_dau1_eta", &d0_dau1_eta, &b_d0_dau1_eta);
   fChain->SetBranchAddress("d0_dau1_phi", &d0_dau1_phi, &b_d0_dau1_phi);
   fChain->SetBranchAddress("d0_dau1_m", &d0_dau1_m, &b_d0_dau1_m);
   fChain->SetBranchAddress("d0_dau1_q", &d0_dau1_q, &b_d0_dau1_q);
   fChain->SetBranchAddress("d0_dau2_pt", &d0_dau2_pt, &b_d0_dau2_pt);
   fChain->SetBranchAddress("d0_dau2_eta", &d0_dau2_eta, &b_d0_dau2_eta);
   fChain->SetBranchAddress("d0_dau2_phi", &d0_dau2_phi, &b_d0_dau2_phi);
   fChain->SetBranchAddress("d0_dau2_m", &d0_dau2_m, &b_d0_dau2_m);
   fChain->SetBranchAddress("d0_dau2_q", &d0_dau2_q, &b_d0_dau2_q);
   fChain->SetBranchAddress("dstar_pt", &dstar_pt, &b_dstar_pt);
   fChain->SetBranchAddress("dstar_eta", &dstar_eta, &b_dstar_eta);
   fChain->SetBranchAddress("dstar_phi", &dstar_phi, &b_dstar_phi);
   fChain->SetBranchAddress("dstar_m", &dstar_m, &b_dstar_m);
   fChain->SetBranchAddress("dstar_true", &dstar_true, &b_dstar_true);
   fChain->SetBranchAddress("dstar_fit", &dstar_fit, &b_dstar_fit);
   fChain->SetBranchAddress("dstar_L3D", &dstar_L3D, &b_dstar_L3D);
   fChain->SetBranchAddress("dstar_LXY", &dstar_LXY, &b_dstar_LXY);
   fChain->SetBranchAddress("dstar_dRTrue", &dstar_dRTrue, &b_dstar_dRTrue);
   fChain->SetBranchAddress("dstar_relPtTrue", &dstar_relPtTrue, &b_dstar_relPtTrue);
   fChain->SetBranchAddress("dstar_dau1_pt", &dstar_dau1_pt, &b_dstar_dau1_pt);
   fChain->SetBranchAddress("dstar_dau1_eta", &dstar_dau1_eta, &b_dstar_dau1_eta);
   fChain->SetBranchAddress("dstar_dau1_phi", &dstar_dau1_phi, &b_dstar_dau1_phi);
   fChain->SetBranchAddress("dstar_dau1_m", &dstar_dau1_m, &b_dstar_dau1_m);
   fChain->SetBranchAddress("dstar_dau1_q", &dstar_dau1_q, &b_dstar_dau1_q);
   fChain->SetBranchAddress("dstar_dau2_pt", &dstar_dau2_pt, &b_dstar_dau2_pt);
   fChain->SetBranchAddress("dstar_dau2_eta", &dstar_dau2_eta, &b_dstar_dau2_eta);
   fChain->SetBranchAddress("dstar_dau2_phi", &dstar_dau2_phi, &b_dstar_dau2_phi);
   fChain->SetBranchAddress("dstar_dau2_m", &dstar_dau2_m, &b_dstar_dau2_m);
   fChain->SetBranchAddress("dstar_dau2_q", &dstar_dau2_q, &b_dstar_dau2_q);
   fChain->SetBranchAddress("dstar_dau3_pt", &dstar_dau3_pt, &b_dstar_dau3_pt);
   fChain->SetBranchAddress("dstar_dau3_eta", &dstar_dau3_eta, &b_dstar_dau3_eta);
   fChain->SetBranchAddress("dstar_dau3_phi", &dstar_dau3_phi, &b_dstar_dau3_phi);
   fChain->SetBranchAddress("dstar_dau3_m", &dstar_dau3_m, &b_dstar_dau3_m);
   fChain->SetBranchAddress("dstar_dau3_q", &dstar_dau3_q, &b_dstar_dau3_q);
   Notify();
}

Bool_t dstarAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void dstarAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t dstarAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef dstarAna_cxx
