import FWCore.ParameterSet.Config as cms

from CATTools.CatProducer.CATMuonProducer_cfi import *
from CATTools.CatProducer.CATElectronProducer_cfi import *
from CATTools.CatProducer.CATPhotonProducer_cfi import *
from CATTools.CatProducer.CATJetProducer_cfi import *
from CATTools.CatProducer.CATMetProducer_cfi import *
from CATTools.CatProducer.CATGenJetProducer_cfi import *
from CATTools.CatProducer.CATGenTopProducer_cfi import *
from CATTools.CatProducer.CATMcParticleProducer_cfi import *

makeCatCandidates =  cms.Sequence( 
    catMuons*
    catElectrons*
    # turning off for now since pf photons dont work
    #catPhotons*
    catJets*
    catMETs*
    # MC information below
    catGenJets*
    # dont need for now
    #catGenTops*
    catMCParticles
) 
