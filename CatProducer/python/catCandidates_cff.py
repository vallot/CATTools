import FWCore.ParameterSet.Config as cms
from CATTools.CatProducer.catEventContent_cff import *

def addCatCommonObjects(process):
    process.load("CATTools.CatProducer.producers.triggerProducer_cfi")
    process.load("CATTools.CatProducer.producers.vertexProducer_cfi")

    process.load("CATTools.CatProducer.producers.muonProducer_cfi")
    process.load("CATTools.CatProducer.producers.electronProducer_cfi")
    process.load("CATTools.CatProducer.producers.jetProducer_cfi")
    process.load("CATTools.CatProducer.producers.metProducer_cfi")

    process.load("CATTools.CatProducer.producers.photonProducer_cfi")
    process.load("CATTools.CatProducer.producers.tauProducer_cfi")
    process.load("CATTools.CatProducer.producers.fatjetProducer_cfi")

    process.RandomNumberGeneratorService.catJets = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(1),
    )
    process.RandomNumberGeneratorService.catFatJets = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(1),
    )
    process.RandomNumberGeneratorService.catJetsPuppi = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(1),
    )

    process.catOut.outputCommands.extend(catEventContent)
    process.catObjectTask.add(
        process.catTrigger, process.catVertex,
        process.catMuons, process.catElectrons, process.catJets, process.catMETs,
        process.catPhotons, process.catTaus, process.catFatJets,
    )
    return process

def addCatCommonMCObjects(process):
    process.load("CATTools.CatProducer.pileupWeight_cff")
    process.load("CATTools.CatProducer.producers.genWeight_cff")

    process.catOut.outputCommands.extend(catEventContentMC)
    process.catObjectTask.add(
        process.pileupWeight, process.genWeight
    )

    return process

def addCatSecVertexObjects(process):
    process.load("CATTools.CatProducer.producers.secondaryVertexProducer_cfi")

    process.catOut.outputCommands.extend(catEventContentSecVertexs)
    process.catObjectTask.add(
        process.catSecVertexs
    )
    return process

def addCatDstarObjects(process):
    process.load("CATTools.CatProducer.producers.dstarProducer_cfi")

    process.catOut.outputCommands.extend(['keep *_catDstars_*_*',])
    process.catObjectTask.add(
        process.catDstars
    )
    return process

def addCatGenTopObjects(process):
    process.load("CATTools.CatProducer.producers.genTopProducer_cfi") #please do not remove it.
    process.load("CATTools.CatProducer.producers.genJetHadronMatch_cfi")
    process.load("TopQuarkAnalysis.TopTools.GenTtbarCategorizer_cfi")

    # for GenTtbarCategories
    from CATTools.CatProducer.Tools.tools import getHFTool
    genHFTool(process, True)
    process.catOut.outputCommands.extend(catEventContentTOPMC)
    process.catObjectTask.add(
        process.catGenTops
    )
    return process

def addCatParticleTopObjects(process):
    process.load("CATTools.CatProducer.mcTruthTop.particleTop_cff")
    process.catOut.outputCommands.extend(catEventContentTOPParticleLevel)
    return process
