import FWCore.ParameterSet.Config as cms

def ntupler_load(process, selectorName, ntuplerName = "ntuple"):
    from CATTools.CatAnalyzer.analyzers.ntuple_cfi import ntuple
    for typeName in ['int', 'ints', 'float', 'floats']:
        for key, val in getattr(ntuple, typeName).parameters_().iteritems():
            setattr(getattr(ntuple, typeName), key, cms.InputTag(val.value().replace('eventSel', selectorName)))
    for key, val in ntuple.cands.parameters_().iteritems():
        getattr(ntuple.cands, key).src = cms.InputTag(val.src.value().replace('eventSel', selectorName))

    setattr(process, ntuplerName, ntuple.clone())
    return process

def ntupler_addVarsTTLL(process, selectorName, ntuplerName = "ntuple"):
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
        src = cms.VInputTag(cms.InputTag(selectorName, "electrons"), cms.InputTag(selectorName, "muons")),
    )
    process.ttllKin.jets = selectorName+":jets"
    process.ttllKin.leptons = "ttllLeptons"
    process.ttllKin.met = selectorName+":met"
    process.ttllKin.metphi = selectorName+":metphi"

    #getattr(process, ntuplerName).floats.weights_genOthers  = cms.InputTag("flatGenWeights", "others")
    getattr(process, ntuplerName).floats.kin_mLL = cms.InputTag("ttllKin", "mLL")
    getattr(process, ntuplerName).floats.kin_mLB = cms.InputTag("ttllKin", "mLB")
    #getattr(process, ntuplerName).floats.kin_mAddJJ = cms.InputTag("ttllKin", "mAddJJ")
    #getattr(process, ntuplerName).floats.kin_dphi = cms.InputTag("ttllKin", "dphi")
    getattr(process, ntuplerName).floats.kin_quality = cms.InputTag("ttllKin", "quality")

    getattr(process, ntuplerName).cands.kin = cms.PSet(
        src = cms.InputTag("ttllKin"),
        exprsF = cms.untracked.PSet(
            pt  = cms.string("pt"),
            eta = cms.string("eta"),
            phi = cms.string("phi"),
            m   = cms.string("mass"),
        ),
    )
    return process

def ntupler_addVarsTTGen(process, ntuplerName = "ntuple"):
    getattr(process, ntuplerName).ints.partonTop_modes = cms.InputTag("partonTop", "modes")

    getattr(process, ntuplerName).cands.partonTop = cms.PSet(
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
    getattr(process, ntuplerName).cands.pseudoTop = cms.PSet(
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

def ntupler_addVarsGen(process, selectorName, ntuplerName = "ntuple"):
    getattr(process, ntuplerName).float.weight = cms.InputTag(selectorName, "weight")
    getattr(process, ntuplerName).float.weight_gen = cms.InputTag("flatGenWeights", "")
    getattr(process, ntuplerName).float.weight_pu   = cms.InputTag("pileupWeight")
    getattr(process, ntuplerName).float.weight_puUp = cms.InputTag("pileupWeight", "up")
    getattr(process, ntuplerName).float.weight_puDn = cms.InputTag("pileupWeight", "dn")

    from CATTools.CatAnalyzer.csvWeights_cfi import csvWeights
    setattr(process, 'csvWeights'+selectorName, csvWeights.clone(src = cms.InputTag(selectorName, "jets")))
    del csvWeights
    getattr(process, ntuplerName).float.weight_CSV = cms.InputTag("csvWeights"+selectorName)
    getattr(process, ntuplerName).floats.weights_CSV = cms.InputTag("csvWeights"+selectorName, "syst")

    getattr(process, ntuplerName).floats.weights_pdf = cms.InputTag("flatGenWeights", "pdf")
    getattr(process, ntuplerName).floats.weights_scaleUp = cms.InputTag("flatGenWeights", "scaleup")
    getattr(process, ntuplerName).floats.weights_scaleDn = cms.InputTag("flatGenWeights", "scaledown")

    return process

def ntupler_addVarsGenTop(process, ntuplerName = "ntuple"):
    getattr(process, ntuplerName).cands.genTop = cms.PSet(
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
