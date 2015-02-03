import FWCore.ParameterSet.Config as cms

from CATTools.CatProducer.muonProducer_cfi import *
from CATTools.CatProducer.electronProducer_cfi import *
from CATTools.CatProducer.photonProducer_cfi import *
from CATTools.CatProducer.jetProducer_cfi import *
from CATTools.CatProducer.metProducer_cfi import *
from CATTools.CatProducer.genJetProducer_cfi import *
from CATTools.CatProducer.genTopProducer_cfi import *
from CATTools.CatProducer.mcParticleProducer_cfi import *
from CATTools.CatProducer.secondaryVertexProducer_cfi import *

makeCatCandidates =  cms.Sequence( 
    catMuons*
    catElectrons*
    catJets*
    catMETs*
    catGenJets
    #catMCParticles
    #catPhotons*
    #catTaus*
    # MC information below
    # dont need for now
    #catGenTops*
    #catSecVertexs
) 
