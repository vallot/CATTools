#include "CATTools/CatAnalyzer/interface/TTEventSelector.h"
#include <cctype>
#include <algorithm>

using namespace std;
using namespace cat;
using namespace TopEventCommonGlobal;


TTEventSelector::TTEventSelector(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& iC ) {
  TTEventSelector(iConfig, iC);
}
TTEventSelector::TTEventSelector(const edm::ParameterSet& iConfig, edm::ConsumesCollector& iC ) {
  typedef std::vector<double> vdouble;
  recoFiltersToken_ = iC.consumes<int>(iConfig.getParameter<edm::InputTag>("recoFilters"));
  nGoodVertexToken_ = iC.consumes<int>(iConfig.getParameter<edm::InputTag>("nGoodVertex"));
  lumiSelectionToken_ = iC.consumes<int>(iConfig.getParameter<edm::InputTag>("lumiSelection"));

  csvWeight.initCSVWeight(false, "csvv2");
  bTagWeightL.init(3, "csvv2", BTagEntry::OP_LOOSE , 1);
  bTagWeightM.init(3, "csvv2", BTagEntry::OP_MEDIUM, 1);
  bTagWeightT.init(3, "csvv2", BTagEntry::OP_TIGHT , 1);


  jetToken_  = iC.consumes<cat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
  metToken_  = iC.consumes<cat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"));
  vtxToken_  = iC.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));

  const auto muonSet = iConfig.getParameter<edm::ParameterSet>("muon");
  muonToken_ = iC.consumes<cat::MuonCollection>(muonSet.getParameter<edm::InputTag>("src"));
  const auto muonSFSet = muonSet.getParameter<edm::ParameterSet>("effSF");
  muonSF_.set(muonSFSet.getParameter<vdouble>("pt_bins"),
              muonSFSet.getParameter<vdouble>("abseta_bins"),
              muonSFSet.getParameter<vdouble>("values"),
              muonSFSet.getParameter<vdouble>("errors"));
  muonPtCut_  = iConfig.getParameter<double>("MuonPtCut");
  muonEtaCut_ = iConfig.getParameter<double>("MuonEtaCut");
  muonIsoCut_ = iConfig.getParameter<double>("MuonIsoCut");
  muonIDCut_  = iConfig.getParameter<string>("MuonIDCut");
  std::transform(muonIDCut_.begin(), muonIDCut_.end(), muonIDCut_.begin(),::tolower);

  const auto elecSet = iConfig.getParameter<edm::ParameterSet>("electron");
  elecToken_ = iC.consumes<cat::ElectronCollection>(elecSet.getParameter<edm::InputTag>("src"));
  const auto elecSFSet = elecSet.getParameter<edm::ParameterSet>("effSF");
  elecSF_.set(elecSFSet.getParameter<vdouble>("pt_bins"),
              elecSFSet.getParameter<vdouble>("abseta_bins"),
              elecSFSet.getParameter<vdouble>("values"),
              elecSFSet.getParameter<vdouble>("errors"));
  electronPtCut_  = iConfig.getParameter<double>("ElectronPtCut" );
  electronEtaCut_ = iConfig.getParameter<double>("ElectronEtaCut");
  electronIsoCut_ = iConfig.getParameter<double>("ElectronIsoCut");
  electronIDCut_  = iConfig.getParameter<string>("ElectronIDCut");

  jetPtCut_  = iConfig.getParameter<double>("JetPtCut");
  jetEtaCut_ = iConfig.getParameter<double>("JetEtaCut");

  bJetCSV_   = iConfig.getParameter<string>("bJetCSV");
  std::transform(bJetCSV_.begin(), bJetCSV_.end(), bJetCSV_.begin(),::tolower);

  /*
  auto solverPSet = iConfig.getParameter<edm::ParameterSet>("solver");
  auto algoName = solverPSet.getParameter<std::string>("algo");
  std::transform(algoName.begin(), algoName.end(), algoName.begin(), ::toupper);
  if      ( algoName == "CMSKIN" ) solver_.reset(new CMSKinSolver(solverPSet));
  else if ( algoName == "DESYMASSLOOP" ) solver_.reset(new DESYMassLoopSolver(solverPSet));
  else if ( algoName == "DESYSMEARED" ) solver_.reset(new DESYSmearedSolver(solverPSet));
  else if ( algoName == "MT2"    ) solver_.reset(new MT2Solver(solverPSet));
  else if ( algoName == "MAOS"   ) solver_.reset(new MAOSSolver(solverPSet));
  else if ( algoName == "DEFAULT" ) solver_.reset(new TTDileptonSolver(solverPSet));
  else {
    cerr << "The solver name \"" << solverPSet.getParameter<std::string>("algo") << "\" is not known please check spellings.\n";
    cerr << "Fall back to the default dummy solver\n";
    solver_.reset(new TTDileptonSolver(solverPSet)); // A dummy solver
  }
  // PseudoTop
  auto solverPSetPT = iConfig.getParameter<edm::ParameterSet>("solverPseudoTop");
  auto algoNamePT = solverPSetPT.getParameter<std::string>("algo");
  std::transform(algoNamePT.begin(), algoNamePT.end(), algoNamePT.begin(), ::toupper);
  if      ( algoNamePT == "CMSKIN" ) solverPT_.reset(new CMSKinSolver(solverPSetPT));
  else if ( algoNamePT == "DESYMASSLOOP" ) solverPT_.reset(new DESYMassLoopSolver(solverPSetPT));
  else if ( algoNamePT == "DESYSMEARED" ) solverPT_.reset(new DESYSmearedSolver(solverPSetPT));
  else if ( algoNamePT == "MT2"    ) solverPT_.reset(new MT2Solver(solverPSetPT));
  else if ( algoNamePT == "MAOS"   ) solverPT_.reset(new MAOSSolver(solverPSetPT));
  else if ( algoNamePT == "DEFAULT" ) solverPT_.reset(new TTDileptonSolver(solverPSetPT));
  else {
    cerr << "The solver name \"" << solverPSetPT.getParameter<std::string>("algoPT") << "\" is not known please check spellings.\n";
    cerr << "Fall back to the default dummy solver\n";
    solverPT_.reset(new TTDileptonSolver(solverPSetPT)); // A dummy solver
  }
  */
}

 
const reco::Candidate* TTEventSelector::getLast(const reco::Candidate* p) const
{
  for ( size_t i=0, n=p->numberOfDaughters(); i<n; ++i )
  {
    const reco::Candidate* dau = p->daughter(i);
    if ( p->pdgId() == dau->pdgId() ) return getLast(dau);
  }
  return p;
}

cat::JetCollection TTEventSelector::selectJets(const cat::JetCollection& jets, const LeptonPtrs& recolep, sys_e sys)
{
  TopEventInfo& evInfo_ = TopEventInfo::getInstance();
  // Initialize SF_btag
  float Jet_SF_CSV[19];
  for (unsigned int iu=0; iu<19; iu++) Jet_SF_CSV[iu] = 1.0;

  cat::JetCollection seljets;
  for (auto& j : jets) {
    cat::Jet jet(j);
    if (sys == sys_jes_u) jet.setP4(j.p4() * j.shiftedEnUp());
    if (sys == sys_jes_d) jet.setP4(j.p4() * j.shiftedEnDown());
    if (sys == sys_jer_u) jet.setP4(j.p4() * j.smearedResUp());
    if (sys == sys_jer_d) jet.setP4(j.p4() * j.smearedResDown());

    if (jet.pt() < jetPtCut_) continue;
    if (std::abs(jet.eta()) > jetEtaCut_ )  continue;
    if (!jet.LooseId()) continue;

    bool hasOverLap = false;
    for (auto lep : recolep){
      if (deltaR(jet.p4(),lep->p4()) < 0.4) hasOverLap = true;
    }
    if (hasOverLap) continue;
    // printf("jet with pt %4.1f\n", jet.pt());
    //if (sys == sys_btag_u) evInfo_.btagweight *= jet.scaleFactorCSVv2(cat::Jet::BTAGCSV_LOOSE, 1);
    //else if (sys == sys_btag_d) evInfo_.btagweight *= jet.scaleFactorCSVv2(cat::Jet::BTAGCSV_LOOSE, -1);
    //else evInfo_.btagweight *= jet.scaleFactorCSVv2(cat::Jet::BTAGCSV_LOOSE, 0);
    for (unsigned int iu=0; iu<19; iu++) Jet_SF_CSV[iu] *= csvWeight.getSF(jet, iu);
    seljets.push_back(jet);
  }
  for (unsigned int iu=0; iu<19; iu++) evInfo_.csvweights.push_back(Jet_SF_CSV[iu]);

  evInfo_.btagweight = Jet_SF_CSV[0];
  // if      ( sys == sys_btag_u ) evInfo_.btagweight = bTagWeightL.eventWeight(seljets, 1);
  // else if ( sys == sys_btag_d ) evInfo_.btagweight = bTagWeightL.eventWeight(seljets, 2);
  // else                          evInfo_.btagweight = bTagWeightL.eventWeight(seljets, 0);

  return seljets;
}

cat::JetCollection TTEventSelector::selectBJets(const JetCollection& jets) const
{
  cat::JetCollection selBjets;
  for (auto& jet : jets) {
    if ( bJetCSV_.compare(std::string("tight"))==0 && jet.bDiscriminator(BTAG_CSVv2) < WP_BTAG_CSVv2T) continue;
    else if ( bJetCSV_.compare(std::string("medium"))==0 && jet.bDiscriminator(BTAG_CSVv2) < WP_BTAG_CSVv2M) continue;
    else if ( bJetCSV_.compare(std::string("loose"))==0 && jet.bDiscriminator(BTAG_CSVv2) < WP_BTAG_CSVv2L) continue;
    //if (jet.bDiscriminator(BTAG_CSVv2) < WP_BTAG_CSVv2M) continue;//forsync
    //printf("b jet with pt %4.1f\n", jet.pt());
    selBjets.push_back(jet);
  }
  return selBjets;
}



 
