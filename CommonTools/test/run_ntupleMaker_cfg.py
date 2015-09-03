import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.load("Configuration.StandardSequences.Services_cff")
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource",
      fileNames = cms.untracked.vstring(
#'/store/user/youn/cat710_phy14_ttbar_2025_aod/catTuple_972.root'
#in Kisti v736
"root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-3-6_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v4/150820_215807/0000/catTuple_1.root",
#'file:catTuple_1.root'
#'file:/cms/home/youn/work/cattool/tag711/cat/src/CATTools/CatProducer/prod/catTuple.root'
#          '/store/user/jlee/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/141219_091640/0000/catTuple_1.root',
#          '/store/user/jlee/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/141219_091640/0000/catTuple_2.root',
#          '/store/user/jlee/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/141219_091640/0000/catTuple_3.root',
#          '/store/user/jlee/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/141219_091640/0000/catTuple_4.root',
#          '/store/user/jlee/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/141219_091640/0000/catTuple_5.root',
      )
)

process.nEventsTotal = cms.EDProducer("EventCountProducer")
process.simpleTriggerMaker = cms.EDProducer("SimpleTriggerMaker",
    InputTriggerLabel  = cms.InputTag("catTrigger"),
    hltPathNames = cms.vstring(
     "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
     "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
     "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
     "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
     "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*"
    ),
)

process.ntuple = cms.EDAnalyzer("GenericNtupleMaker",
    failureMode = cms.untracked.string("keep"), # choose one among keep/skip/error
    eventCounters = cms.vstring("nEventsTotal"), #"nEventsTotal", "nEventsClean", "nEventsPAT"),
    int = cms.PSet(
        HLTMu17TrkIsoVVLMu8TrkIsoVVLDZ            =   cms.PSet(src = cms.InputTag("simpleTriggerMaker", "HLTMu17TrkIsoVVLMu8TrkIsoVVLDZ",                         )),
        HLTMu17TrkIsoVVLTkMu8TrkIsoVVLDZ          =   cms.PSet(src = cms.InputTag("simpleTriggerMaker", "HLTMu17TrkIsoVVLTkMu8TrkIsoVVLDZ",                       )),
        HLTEle17Ele12CaloIdLTrackIdLIsoVLDZ        =   cms.PSet(src = cms.InputTag("simpleTriggerMaker", "HLTEle17Ele12CaloIdLTrackIdLIsoVLDZ",                   )),
        HLTMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVL  =   cms.PSet(src = cms.InputTag("simpleTriggerMaker", "HLTMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVL",             )),
        HLTMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVL   =   cms.PSet(src = cms.InputTag("simpleTriggerMaker", "HLTMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVL"               )),
   
        #HLTDoubleEle33CaloIdLGsfTrkIdVL             =   cms.PSet(src = cms.InputTag("catTrigger", "HLTDoubleEle33CaloIdLGsfTrkIdVL"             )),
        #HLTEle12CaloIdLTrackIdLIsoVL                =   cms.PSet(src = cms.InputTag("catTrigger", "HLTEle12CaloIdLTrackIdLIsoVL"                )),
        #HLTEle16Ele12Ele8CaloIdLTrackIdL            =   cms.PSet(src = cms.InputTag("catTrigger", "HLTEle16Ele12Ele8CaloIdLTrackIdL"            )),
        #HLTEle17CaloIdLTrackIdLIsoVL                =   cms.PSet(src = cms.InputTag("catTrigger", "HLTEle17CaloIdLTrackIdLIsoVL"                )),
        #HLTEle17Ele12CaloIdLTrackIdLIsoVLDZ         =   cms.PSet(src = cms.InputTag("catTrigger", "HLTEle17Ele12CaloIdLTrackIdLIsoVLDZ"         )),
        #HLTEle23Ele12CaloIdLTrackIdLIsoVL           =   cms.PSet(src = cms.InputTag("catTrigger", "HLTEle23Ele12CaloIdLTrackIdLIsoVL"           )),
        #HLTEle23Ele12CaloIdLTrackIdLIsoVLDZ         =   cms.PSet(src = cms.InputTag("catTrigger", "HLTEle23Ele12CaloIdLTrackIdLIsoVLDZ"         )),
        #HLTEle27eta2p1WPLooseGsfTriCentralPFJet30   =   cms.PSet(src = cms.InputTag("catTrigger", "HLTEle27eta2p1WPLooseGsfTriCentralPFJet30"   )),
        #HLTMu17Mu8DZ                                =   cms.PSet(src = cms.InputTag("catTrigger", "HLTMu17Mu8DZ"                                )),
        #HLTMu17TkMu8DZ                              =   cms.PSet(src = cms.InputTag("catTrigger", "HLTMu17TkMu8DZ"                              )),
        #HLTMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVL   =   cms.PSet(src = cms.InputTag("catTrigger", "HLTMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVL"   )),
        #HLTMu17TrkIsoVVLMu8TrkIsoVVL                =   cms.PSet(src = cms.InputTag("catTrigger", "HLTMu17TrkIsoVVLMu8TrkIsoVVL"                )),
        #HLTMu17TrkIsoVVLMu8TrkIsoVVLDZ              =   cms.PSet(src = cms.InputTag("catTrigger", "HLTMu17TrkIsoVVLMu8TrkIsoVVLDZ"              )),
        #HLTMu17TrkIsoVVLTkMu8TrkIsoVVL              =   cms.PSet(src = cms.InputTag("catTrigger", "HLTMu17TrkIsoVVLTkMu8TrkIsoVVL"              )),
        #HLTMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVL    =   cms.PSet(src = cms.InputTag("catTrigger", "HLTMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVL"    )),

        nGoodPV           =   cms.PSet(src = cms.InputTag("catVertex"   , "nGoodPV"          )),
        nPV               =   cms.PSet(src = cms.InputTag("catVertex"   , "nPV"              )),

        pdfWeightId1 =   cms.PSet(src = cms.InputTag("pdfWeight", "id1" )),
        pdfWeightId2 =   cms.PSet(src = cms.InputTag("pdfWeight", "id2" )),

        nTrueInteraction  =   cms.PSet(src = cms.InputTag("pileupWeight", "nTrueInteraction" )),
    ),
    float = cms.PSet(
        puWeight   = cms.PSet(src = cms.InputTag("pileupWeight")),
        puWeightUp = cms.PSet(src = cms.InputTag("pileupWeight", "up")),
        puWeightDn = cms.PSet(src = cms.InputTag("pileupWeight", "dn")),

        pdfWeightQ  =   cms.PSet(src = cms.InputTag("pdfWeight", "Q" )),
        pdfWeightX1 =   cms.PSet(src = cms.InputTag("pdfWeight", "x1" )),
        pdfWeightX2 =   cms.PSet(src = cms.InputTag("pdfWeight", "x2" )),

    ),
    floats = cms.PSet(
        pdfWeight = cms.PSet(src = cms.InputTag("pdfWeight")),
    ),
    cands = cms.PSet(
        muons = cms.PSet(
            src = cms.InputTag("catMuons"),
            #index = cms.untracked.int32(0),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                #relIso = cms.string("relIso"),
                relIso03 = cms.string("relIso(0.3)"),
                relIso04 = cms.string("relIso(0.4)"),
                isTracker = cms.string("isTrackerMuon"),
                isGlobal = cms.string("isGlobalMuon"),
                isLoose = cms.string("isLooseMuon"),
                isTight = cms.string("isTightMuon"),
                dxy = cms.string("dxy"),
                dz = cms.string("dz"),
                q = cms.string("charge"),
                isPF = cms.string("isPFMuon"),
                normalizedChi2 = cms.string("normalizedChi2"),
                numberOfValidMuonHits = cms.string("numberOfValidMuonHits"),
                numberOfValidPixelHits = cms.string("numberOfValidPixelHits"),
                trackerLayersWithMeasurement = cms.string("trackerLayersWithMeasurement"),
                numberOfMatchedStations = cms.string("numberOfMatchedStations"),
                #matched = cms.string("mcMatched"),
            ),
            selections = cms.untracked.PSet(),
        ),
        electrons = cms.PSet(
            src = cms.InputTag("catElectrons"),
            #index = cms.untracked.int32(0),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                #relIso = cms.string("relIso"),
                #idLoose = cms.string("electronID('eidLoose')"),
                #idTight = cms.string("electronID('eidTight')"),
                #idMedium2 = cms.string("electronID('cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium')"),
                #idVeto2 = cms.string("electronID('cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto')"),
                #idVeto = cms.string("electronID('cutBasedElectronID-CSA14-PU20bx25-V0-standalone-veto')"),
                #idMedium2= cms.string("electronID('cutBasedElectronID_CSA14_PU20bx25_V0_standalone_medium')"),
                #idMedium= cms.string("electronID('cutBasedElectronID-CSA14-PU20bx25-V0-standalone-medium')"),
                #mva = cms.string("electronID('mvaTrigV0')"),
                idHEEP51 = cms.string("electronID('heepElectronID-HEEPV51')"),
                idLoose = cms.string("electronID('cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose')"), 
                idMedium = cms.string("electronID('cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium')"), 
                idTight = cms.string("electronID('cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight')"), 
                idTeto = cms.string("electronID('cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto')"),
                relIso03 = cms.string("relIso(0.3)"),
                relIso04 = cms.string("relIso(0.4)"),
                chIso03 = cms.string("chargedHadronIso(0.3)"),
                nhIso03 = cms.string("neutralHadronIso(0.3)"),
                phIso03 = cms.string("photonIso(0.3)"),
                puChIso03 = cms.string("puChargedHadronIso(0.3)"),
                chIso04 = cms.string("chargedHadronIso(0.4)"),
                nhIso04 = cms.string("neutralHadronIso(0.4)"),
                phIso04 = cms.string("photonIso(0.4)"),
                puChIso04 = cms.string("puChargedHadronIso(0.4)"), 
                #rhoIso03 = cms.string("rho"),
                scEta = cms.string("scEta"),
                #dxy = cms.string("dxy"),
                #dz = cms.string("dz"),
                q = cms.string("charge"),
                #isGsfCtfScPixChargeConsistent = cms.string("isGsfCtfScPixChargeConsistent"),
            ),
            selections = cms.untracked.PSet(
                isPassBaseId = cms.string("passConversionVeto && isPF && gsfTrack.hitPattern.numberOfLostHits('MISSING_INNER_HITS') <= 0"),
            ),
        ),
        jets = cms.PSet(
            src = cms.InputTag("catJets"),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                vtxMass = cms.string("vtxMass"),
                CSVInclV2 = cms.string("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')"),
                #CSVInclV2 = cms.string("bDiscriminator('combinedInclusiveSecondaryVertexV2BJetTags')"),
                partonFlavour = cms.string("partonFlavour"),
                hadronFlavour = cms.string("hadronFlavour"),
            ),
            selections = cms.untracked.PSet(
                isLoose = cms.string("LooseId"),
                isPFId = cms.string("pileupJetId"),
            ),
        ),
        jetsPuppi = cms.PSet(
            src = cms.InputTag("catJetsPuppi"),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                vtxMass = cms.string("vtxMass"),
                CSVInclV2 = cms.string("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')"),
                #CSVInclV2 = cms.string("bDiscriminator('combinedInclusiveSecondaryVertexV2BJetTags')"),
                partonFlavour = cms.string("partonFlavour"),
                hadronFlavour = cms.string("hadronFlavour"),
            ),
            selections = cms.untracked.PSet(
                isLoose = cms.string("LooseId"),
                isPFId = cms.string("pileupJetId"),
            ),
        ),
        met = cms.PSet(
            src = cms.InputTag("catMETs"),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                phi = cms.string("phi"),
            ),
            selections = cms.untracked.PSet(),
        ),
        metNoHF = cms.PSet(
            src = cms.InputTag("catMETsNoHF"),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                phi = cms.string("phi"),
            ),
            selections = cms.untracked.PSet(),
        ),
        metPfMva = cms.PSet(
            src = cms.InputTag("catMETsPfMva"),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                phi = cms.string("phi"),
            ),
            selections = cms.untracked.PSet(),
        ),
        metPuppi = cms.PSet(
            src = cms.InputTag("catMETsPuppi"),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                phi = cms.string("phi"),
            ),
            selections = cms.untracked.PSet(),
        ),
        #gentop = cms.PSet(
        #    src = cms.InputTag("catGenTops"),
        #    #index = cms.untracked.int32(0),
        #    exprs = cms.untracked.PSet(
        #        lepton1_pt      = cms.string("lepton1().Pt()"),
        #        lepton1_eta      = cms.string("lepton1().Eta()"),
        #        lepton1_phi      = cms.string("lepton1().Phi()"),
        #        lepton2_pt      = cms.string("lepton2().Pt()"),
        #        lepton2_eta      = cms.string("lepton2().Eta()"),
        #        lepton2_phi      = cms.string("lepton2().Phi()"),
        #        allHadronic      = cms.string("allHadronic"),
        #        semiLeptonic     = cms.string("semiLeptonic"),
        #        diLeptonicMuoMuo = cms.string("diLeptonicMuoMuo"),
        #        diLeptonicMuoEle = cms.string("diLeptonicMuoEle"),
        #        diLeptonicEleEle = cms.string("diLeptonicEleEle"),
        #        diLeptonicTauMuo = cms.string("diLeptonicTauMuo"),
        #        diLeptonicTauEle = cms.string("diLeptonicTauEle"),
        #        diLeptonicTauTau = cms.string("diLeptonicTauTau"),
        #        NbJets1           = cms.string("NbJets(1)"),
        #        NbJets201         = cms.string("NbJets20(1)"),
        #        NbJets251         = cms.string("NbJets25(1)"),
        #        NbJets301         = cms.string("NbJets30(1)"),
        #        NbJets401         = cms.string("NbJets40(1)"),
        #        NaddbJets1        = cms.string("NaddbJets(1)"),
        #        NaddbJets201      = cms.string("NaddbJets20(1)"),
        #        NaddbJets401      = cms.string("NaddbJets40(1)"),
        #        NcJets1           = cms.string("NcJets(1)"),
        #        NcJets101         = cms.string("NcJets10(1)"),
        #        NcJets151         = cms.string("NcJets15(1)"),
        #        NcJets201         = cms.string("NcJets20(1)"),
        #        NcJets251         = cms.string("NcJets25(1)"),
        #        NcJets301         = cms.string("NcJets30(1)"),
        #        NcJets401         = cms.string("NcJets40(1)"),
        #        NbJets           = cms.string("NbJets(0)"),
        #        NbJets20         = cms.string("NbJets20(0)"),
        #        NbJets25         = cms.string("NbJets25(0)"),
        #        NbJets30         = cms.string("NbJets30(0)"),
        #        NbJets40         = cms.string("NbJets40(0)"),
        #        NaddbJets        = cms.string("NaddbJets(0)"),
        #        NaddbJets20      = cms.string("NaddbJets20(0)"),
        #        NaddbJets40      = cms.string("NaddbJets40(0)"),
        #        NcJets           = cms.string("NcJets(0)"),
        #        NcJets10         = cms.string("NcJets10(0)"),
        #        NcJets15         = cms.string("NcJets15(0)"),
        #        NcJets20         = cms.string("NcJets20(0)"),
        #        NcJets25         = cms.string("NcJets25(0)"),
        #        NcJets30         = cms.string("NcJets30(0)"),
        #        NcJets40         = cms.string("NcJets40(0)"),
        #        NJets            = cms.string("NJets"),
        #        NJets10          = cms.string("NJets10"),
        #        NJets20          = cms.string("NJets20"),
        #        NJets25          = cms.string("NJets25"),
        #        NJets30          = cms.string("NJets30"),
        #        NJets40          = cms.string("NJets40"),
        #        NaddJets20       = cms.string("NaddJets20"),
        #        NaddJets40       = cms.string("NaddJets40"),
        #        NbQuarksTop      = cms.string("NbQuarksTop"),
        #        NbQuarksNoTop    = cms.string("NbQuarksNoTop"),
        #        NbQuarks         = cms.string("NbQuarks"),
        #        NbQuarks20       = cms.string("NbQuarks20"),
        #        NbQuarks40       = cms.string("NbQuarks40"), 
        #        NaddbQuarks20    = cms.string("NaddbQuarks20"),
        #        NaddbQuarks40    = cms.string("NaddbQuarks40"),
        #        NcQuarks         = cms.string("NcQuarks"),
        #    ),
        #    selections = cms.untracked.PSet(),
        #),
        slimmedGenJets = cms.PSet(
            src = cms.InputTag("slimmedGenJets",""),
            #index = cms.untracked.int32(0),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                #pdgId = cms.string("pdgId"),
                #q = cms.string("charge"),
                #status = cms.string("status"),
            ),
            selections = cms.untracked.PSet(),
        ),
        partonTop = cms.PSet(
            src = cms.InputTag("partonTop"),
            #index = cms.untracked.int32(0),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                pdgId = cms.string("pdgId"),
                q = cms.string("charge"),
                #status = cms.string("status"),
            ),
            selections = cms.untracked.PSet(),
        ),
        pseudoTopJet = cms.PSet(
            src = cms.InputTag("pseudoTop","jets"),
            #index = cms.untracked.int32(0),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                pdgId = cms.string("pdgId"),
                q = cms.string("charge"),
                #status = cms.string("status"),
            ),
            selections = cms.untracked.PSet(),
        ),
        pseudoTopLepton = cms.PSet(
            src = cms.InputTag("pseudoTop","leptons"),
            #index = cms.untracked.int32(0),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                pdgId = cms.string("pdgId"),
                q = cms.string("charge"),
                #status = cms.string("status"),
            ),
            selections = cms.untracked.PSet(),
        ),
        pseudoTopNu = cms.PSet(
            src = cms.InputTag("pseudoTop","neutrinos"),
            #index = cms.untracked.int32(0),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                #m   = cms.string("mass"),
                pdgId = cms.string("pdgId"),
                #q = cms.string("charge"),
                #status = cms.string("status"),
            ),
            selections = cms.untracked.PSet(),
        ),
        pseudoTop = cms.PSet(
            src = cms.InputTag("pseudoTop"),
            #index = cms.untracked.int32(0),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                pdgId = cms.string("pdgId"),
                q = cms.string("charge"),
                #status = cms.string("status"),
            ),
            selections = cms.untracked.PSet(),
        ),
    ),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)

process.load("CATTools.CatProducer.pseudoTop_cff")
process.p = cms.Path(
    process.nEventsTotal*
    process.partonTop*
    process.simpleTriggerMaker*
    process.ntuple
)

