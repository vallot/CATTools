import FWCore.ParameterSet.Config as cms

def ntupler_load(process, name):
    process.load("CATTools.CatAnalyzer.analyzers.ntuple_cfi")
    for typeName in ['int', 'ints', 'float', 'floats']:
        for key, val in getattr(process.ntuple, typeName).parameters_().iteritems():
            setattr(getattr(process.ntuple, typeName), key, cms.InputTag(val.value().replace('eventSel', name)))
    for key, val in process.ntuple.cands.parameters_().iteritems():
        getattr(process.ntuple.cands, key).src = cms.InputTag(val.src.value().replace('eventSel', name))
    return process

def ntupler_addVarsTTLL(process, name):
    #process.load("CATTools.CatAnalyzer.ttll.ttbarDileptonKinSolutionProducer_cfi")
    from CATTools.CatAnalyzer.ttll.ttbarDileptonKinSolutionProducer_cfi import ttbarDileptonKin
    from CATTools.CatAnalyzer.ttll.ttbarDileptonKinSolutionAlgos_cff import ttbarDileptonKinAlgoPSetDESYSmeared as solverAlgo
    process.ttllKin = ttbarDileptonKin.clone()
    del ttbarDileptonKin
    process.ttllKin.solver = solverAlgo

    process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
        ttllKin = cms.PSet(
            initialSeed = cms.untracked.uint32(123456),
            engineName = cms.untracked.string('TRandom3')
        )
    )

    process.ttllLeptons = cms.EDProducer("CandViewMerger",
        src = cms.VInputTag(cms.InputTag(name, "electrons"), cms.InputTag(name, "muons")),
    )
    process.ttllKin.jets = name+":jets"
    process.ttllKin.leptons = "ttllLeptons"
    process.ttllKin.met = name+":met"
    process.ttllKin.metphi = name+":metphi"

    #process.ntuple.floats.weights_genOthers  = cms.InputTag("flatGenWeights", "others")
    process.ntuple.floats.kin_mLL = cms.InputTag("ttllKin", "mLL")
    process.ntuple.floats.kin_mLB = cms.InputTag("ttllKin", "mLB")
    #process.ntuple.floats.kin_mAddJJ = cms.InputTag("ttllKin", "mAddJJ")
    #process.ntuple.floats.kin_dphi = cms.InputTag("ttllKin", "dphi")
    process.ntuple.floats.kin_quality = cms.InputTag("ttllKin", "quality")

    process.ntuple.cands.kin = cms.PSet(
        src = cms.InputTag("ttllKin"),
        exprsF = cms.untracked.PSet(
            pt  = cms.string("pt"),
            eta = cms.string("eta"),
            phi = cms.string("phi"),
            m   = cms.string("mass"),
        ),
    )
    return process

def ntupler_addVarsTTGen(process):
    process.ntuple.ints.partonTop_modes = cms.InputTag("partonTop", "modes")

    process.ntuple.cands.partonTop = cms.PSet(
        src = cms.InputTag("partonTop"),
        exprsF = cms.untracked.PSet(
            pt  = cms.string("pt"),
            eta = cms.string("eta"),
            phi = cms.string("phi"),
            m   = cms.string("mass"),
        ),
        exprsI = cms.untracked.PSet(
            q = cms.string("charge"),
            pdgId = cms.string("pdgId"),
        ),
    )
    process.ntuple.cands.pseudoTop = cms.PSet(
        src = cms.InputTag("pseudoTop"),
        exprsF = cms.untracked.PSet(
            pt  = cms.string("pt"),
            eta = cms.string("eta"),
            phi = cms.string("phi"),
            m   = cms.string("mass"),
        ),
        exprsI = cms.untracked.PSet(
            q = cms.string("charge"),
            pdgId = cms.string("pdgId"),
        ),
    )

    return process

def ntupler_addVarsGen(process, name):
    process.ntuple.float.weight = cms.InputTag(name, "weight")
    process.ntuple.float.weight_gen = cms.InputTag("flatGenWeights", "")
    process.ntuple.float.weight_pu   = cms.InputTag("pileupWeight")
    process.ntuple.float.weight_puUp = cms.InputTag("pileupWeight", "up")
    process.ntuple.float.weight_puDn = cms.InputTag("pileupWeight", "dn")

    from CATTools.CatAnalyzer.csvWeights_cfi import csvWeights
    setattr(process, 'csvWeights'+name, csvWeights.clone(src = cms.InputTag(name, "jets")))
    del csvWeights
    process.ntuple.float.weight_CSV = cms.InputTag("csvWeights"+name)
    process.ntuple.floats.weights_CSV = cms.InputTag("csvWeights"+name, "syst")

    process.ntuple.floats.weights_pdf = cms.InputTag("flatGenWeights", "pdf")
    process.ntuple.floats.weights_scaleUp = cms.InputTag("flatGenWeights", "scaleup")
    process.ntuple.floats.weights_scaleDn = cms.InputTag("flatGenWeights", "scaledown")

    return process

def ntupler_addVarsGenTop(process):
    process.ntuple.cands.genTop = cms.PSet(
        src = cms.InputTag("catGenTops"),
        index = cms.untracked.int32(0),
        exprsF = cms.untracked.PSet(
            top1_pt  = cms.string("topquark1().pt()"),
            top1_eta = cms.string("topquark1().eta()"),
            top1_phi = cms.string("topquark1().phi()"),
            top1_m   = cms.string("topquark1().mass()"),
            top2_pt  = cms.string("topquark2().pt()"),
            top2_eta = cms.string("topquark2().eta()"),
            top2_phi = cms.string("topquark2().phi()"),
            top2_m   = cms.string("topquark2().mass()"),

            lep1_pt  = cms.string("lepton1().pt()"),
            lep1_eta = cms.string("lepton1().eta()"),
            lep1_phi = cms.string("lepton1().phi()"),
            lep1_m   = cms.string("lepton1().mass()"),
            lep2_pt  = cms.string("lepton2().pt()"),
            lep2_eta = cms.string("lepton2().eta()"),
            lep2_phi = cms.string("lepton2().phi()"),
            lep2_m   = cms.string("lepton2().mass()"),

            bjet1_pt  = cms.string("bJetsFromTop1().pt()"),
            bjet1_eta = cms.string("bJetsFromTop1().eta()"),
            bjet1_phi = cms.string("bJetsFromTop1().phi()"),
            bjet1_m   = cms.string("bJetsFromTop1().mass()"),
            bjet2_pt  = cms.string("bJetsFromTop2().pt()"),
            bjet2_eta = cms.string("bJetsFromTop2().eta()"),
            bjet2_phi = cms.string("bJetsFromTop2().phi()"),
            bjet2_m   = cms.string("bJetsFromTop2().mass()"),

            addjet1_pt  = cms.string("addJets1().pt()"),
            addjet1_eta = cms.string("addJets1().eta()"),
            addjet1_phi = cms.string("addJets1().phi()"),
            addjet1_m   = cms.string("addJets1().mass()"),
            addjet2_pt  = cms.string("addJets2().pt()"),
            addjet2_eta = cms.string("addJets2().eta()"),
            addjet2_phi = cms.string("addJets2().phi()"),
            addjet2_m   = cms.string("addJets2().mass()"),

            addjets_dR = cms.string("dRaddJets()"),
        ),
        exprsI = cms.untracked.PSet(
            bjets20_n = cms.string("NbJets20()"),
            bjets40_n = cms.string("NbJets40()"),
            addbjets20_n = cms.string("NaddbJets20()"),
            addbjets40_n = cms.string("NaddbJets40()"),
            addcjets20_n = cms.string("NaddcJets20()"),
            addcjets40_n = cms.string("NaddcJets40()"),
            wjets_n = cms.string("NWJets()"),
        ),
        exprsB = cms.untracked.PSet(
            semiLeptonic = cms.string("semiLeptonic()"),
            diLeptonic = cms.string("diLeptonic()"),
        ),
    )
    return process
