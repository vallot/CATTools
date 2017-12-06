import FWCore.ParameterSet.Config as cms

#from TopQuarkAnalysis.TopEventProducers.producers.pseudoTop_cfi import *
from GeneratorInterface.RivetInterface.particleLevel_cfi import *
from GeneratorInterface.RivetInterface.genParticles2HepMC_cff import *
from GeneratorInterface.RivetInterface.mergedGenParticles_cfi import *
genParticles2HepMC.genParticles = "mergedGenParticles"
