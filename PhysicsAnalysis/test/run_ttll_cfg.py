import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.load("Configuration.StandardSequences.Services_cff")
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
for i in xrange(1,101):
    process.source.fileNames.append('/store/user/jlee/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/141219_091640/0000/catTuple_%d.root' % i)
#    '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU40bx25_POSTLS170_V7-v2/00000/0A30732D-FA26-E411-9A59-E0CB4E29C4FD.root',
#'/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/16964A21-7344-E311-BBE6-00A0D1EE8ECC.root',
#]

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("out.root"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_*_*_Ana",
    )
)

#process.outPath = cms.EndPath(process.out)

process.ttbar = cms.EDProducer("TTbarDileptonProducer",
#    solver = cms.string("Default"),
#    solver = cms.string("CMSKIN"),
#    solver = cms.string("NUWGT"),
    solver = cms.string("MT2"),
    muons = cms.InputTag("catMuons"),
    electrons = cms.InputTag("catElectrons"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag("catMETs"),
    genParticles = cms.InputTag("prunedGenParticles"),
)

process.ntuple = cms.EDAnalyzer("GenericNtupleMaker",
    failureMode = cms.untracked.string("error"), # choose one among keep/skip/error
    eventCounters = cms.vstring(), #"nEventsTotal", "nEventsClean", "nEventsPAT"),
    int = cms.PSet(
        channel = cms.PSet(src = cms.InputTag("ttbar:channel")),
    ),
    double = cms.PSet(
    ),
    doubles = cms.PSet(
        ttbar_mLL = cms.PSet(src = cms.InputTag("ttbar:mLL")),
        ttbar_dphi = cms.PSet(src = cms.InputTag("ttbar:dphi")),
        ttbar_mAddJJ = cms.PSet(src = cms.InputTag("ttbar:mAddJJ")),
        ttbar_mLB = cms.PSet(src = cms.InputTag("ttbar:mLB")),
        ttbar_aux = cms.PSet(src = cms.InputTag("ttbar:aux")),
    ),
    cands = cms.PSet(
        ttbar = cms.PSet(
            src = cms.InputTag("ttbar"),
            index = cms.untracked.int32(0),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),

                top1pt  = cms.string("daughterPtr(0).pt"),
                top1eta = cms.string("daughterPtr(0).eta"),
                top1phi = cms.string("daughterPtr(0).phi"),
                top1m   = cms.string("daughterPtr(0).mass"),

                top2pt  = cms.string("daughterPtr(1).pt"),
                top2eta = cms.string("daughterPtr(1).eta"),
                top2phi = cms.string("daughterPtr(1).phi"),
                top2m   = cms.string("daughterPtr(1).mass"),

                w1m   = cms.string("daughterPtr(0).daughter(0).mass"),
                w2m   = cms.string("daughterPtr(1).daughter(0).mass"),
                b1Pt = cms.string("daughterPtr(0).daughter(1).pt"),
                b2Pt = cms.string("daughterPtr(1).daughter(1).pt"),
                b1Tag = cms.string("daughterPtr(0).daughter(1).bDiscriminator('combinedInclusiveSecondaryVertexV2BJetTags')"),
                b2Tag = cms.string("daughterPtr(1).daughter(1).bDiscriminator('combinedInclusiveSecondaryVertexV2BJetTags')"),

                lep1Pt = cms.string("daughter(0).daughter(0).daughter(0).pt"),
                lep2Pt = cms.string("daughter(1).daughter(0).daughter(0).pt"),
                nu1Pt = cms.string("daughter(0).daughter(0).daughter(1).pt"),
                nu2Pt = cms.string("daughter(1).daughter(0).daughter(1).pt"),
            ),
            selections = cms.untracked.PSet(
            )
        ),
    ),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

process.p = cms.Path(
    process.ttbar
  * process.ntuple
)

