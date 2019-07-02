import FWCore.ParameterSet.Config as cms

def computeAverageSF(set1, lumi1, set2, lumi2):
    sfSet = set1.clone()

    w1 = lumi1/(lumi1+lumi2)
    w2 = lumi2/(lumi1+lumi2)
    for i in range(len(sfSet.values)):
        wv1, wv2 = w1*set1.values[i], w2*set2.values[i]
        we1, we2 = w1*set1.errors[i], w2*set2.errors[i]
        sfSet.values[i] = wv1+wv2
        sfSet.errors[i] = (we1*we1 + we2*we2)**0.5

    return sfSet

def combineSF(set1, set2, additionalUnc1=0, additionalUnc2=0):
    from bisect import bisect_right

    ## build pt bins (x axis)
    binsX1, binsX2 = set1.pt_bins[:], set2.pt_bins[:]
    nbinsX1, nbinsX2 = len(binsX1), len(binsX2)
    binsX = sorted(list(set().union(binsX1, binsX2)))

    ## build (abs) eta bins (y axis)
    isAbsEta1 = hasattr(set1, 'abseta_bins')
    isAbsEta2 = hasattr(set2, 'abseta_bins')
    binNameY, binNameY1, binNameY2 = 'eta_bins', 'eta_bins', 'eta_bins'
    if isAbsEta1: binNameY1 = 'abseta_bins'
    if isAbsEta2: binNameY2 = 'abseta_bins'
    if isAbsEta1 == isAbsEta2 == True: binNameY = 'abseta_bins'

    binsY1 = getattr(set1, binNameY1)[:]
    binsY2 = getattr(set2, binNameY2)[:]
    nbinsY1, nbinsY2 = len(binsY1), len(binsY2)
    if isAbsEta1 == isAbsEta2:
        binsY = sorted(list(set().union(binsY1, binsY2)))
    else:
        print "Expanding bins"
        binsY = sorted(list(set().union(binsY1, [-y for y in binsY1],
                                        binsY2, [-y for y in binsY2])))

    values, errors = [], []
    for y in binsY[:-1]:
        y1, y2 = y, y
        if isAbsEta1 and not isAbsEta2: y1 = abs(y)
        if isAbsEta2 and not isAbsEta1: y2 = abs(y)

        binY1 = bisect_right(binsY1, y1)
        binY2 = bisect_right(binsY2, y2)
        for x in binsX[:-1]:
            binX1 = bisect_right(binsX1, x)
            binX2 = bisect_right(binsX2, x)

            value1, value2, error1, error2 = 1, 1, 0, 0

            if 0 < binX1 < nbinsX1 and 0 < binY1 < nbinsY1:
                bin1 = (binX1-1)+(binY1-1)*(nbinsX1-1)
                value1 = set1.values[bin1]
                error1 = set1.errors[bin1]
            if 0 < binX2 < nbinsX2 and 0 < binY2 < nbinsY2:
                bin2 = (binX2-1)+(binY2-1)*(nbinsX2-1)
                value2 = set2.values[bin2]
                error2 = set2.errors[bin2]

            addErr1 = additionalUnc1*value1
            addErr2 = additionalUnc2*value2

            values.append(value1*value2)
            errors.append((error1**2+error2**2+addErr1**2+addErr2**2)**0.5)

    sfSet = cms.PSet(
        pt_bins = cms.vdouble(binsX),
        values = cms.vdouble(values),
        errors = cms.vdouble(errors),
    )
    setattr(sfSet, binNameY, cms.vdouble(binsY))

    return sfSet

dummySF = cms.PSet(
    pt_bins = cms.vdouble(0,1e9),
    abseta_bins = cms.vdouble(0,1e9),
    eta_bins = cms.vdouble(-1e9,1e9),
    values = cms.vdouble(1.0),
    errors = cms.vdouble(0.0),
)

##Muon SF reference for 2018: https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2018
#The numbers are retrieved from json file in gitlab
#https://gitlab.cern.ch/cms-muonPOG/MuonReferenceEfficiencies/blob/master/EfficienciesStudies/2018/jsonfiles/RunABCD_SF_ISO.json
#https://gitlab.cern.ch/cms-muonPOG/MuonReferenceEfficiencies/blob/master/EfficienciesStudies/2018/jsonfiles/RunABCD_SF_ID.json
#Uncertainties are square-summed

#NUM_TightID_DEN_TrackerMuons, rounded keeping 7 significant figures
muonSFTightId102X  = cms.PSet(
    pt_bins = cms.vdouble(20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 120.0,),
    abseta_bins = cms.vdouble(0.0, 0.9, 1.2, 2.1, 2.4,),
    values = cms.vdouble(
      0.9913892, 0.9914348, 0.9913756, 0.9914412, 0.9916862, 0.9905799, 
      0.9859874, 0.9839878, 0.9838238, 0.9842166, 0.9830003, 0.9816889,
      0.9914728, 0.9911847, 0.9907980, 0.9899933, 0.9907797, 0.9886367,
      0.9762337, 0.9735941, 0.9739650, 0.9738877, 0.9737076, 0.9677501,
    ),
    errors = cms.vdouble(
      0.007422989, 0.001385032, 0.00195549, 0.0007994815, 0.000647879, 0.002601976,
      0.006915646, 0.002156646, 0.002721865, 0.0008127334, 0.001451969, 0.002571885,
      0.006162521, 0.002709858, 0.001138584, 0.0009250287, 0.001043604, 0.0008427016,
      0.0102527, 0.001318936, 0.002176843, 0.001642799, 0.003151512, 0.005727141,
    ),
)

#"NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta" for 2018
muonSFTightIso102X = cms.PSet(
    pt_bins = cms.vdouble(20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 120.0,),
    abseta_bins = cms.vdouble(0.0, 0.9, 1.2, 2.1, 2.4,),
    values = cms.vdouble(
      0.9878177, 0.9888813, 0.9935667, 0.9966412, 0.9972884, 0.9986075,
      0.9862650, 0.9910231, 0.9937724, 0.9967368, 0.9974498, 0.9994680,
      1.0287069, 1.0199128, 1.0102499, 1.0051965, 1.0029043, 1.0020165,
      1.0646216, 1.0432520, 1.0249592, 1.0124054, 1.0071798, 1.0058323,
    ),
    errors = cms.vdouble(
      0.003172181, 0.001493586, 0.0004602945, 0.0002128453, 0.0003927085, 0.0005351785,
      0.005238489, 0.002669965, 0.0008072914, 0.000345866, 0.0007226002, 0.001017408,
      0.002466541, 0.001392536, 0.0006276717, 0.0002431922, 0.0004399094, 0.0006245023,
      0.004033066, 0.002196287, 0.0008832271, 0.0005292326, 0.001031393, 0.001686611,
    ),
)


## ELectron SF for 2018: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#102X_series_Dataset_2018_Autumn
electronSFReco102X = cms.PSet(
    pt_bins = cms.vdouble(10, 20, 45, 75, 100, 500),
    eta_bins = cms.vdouble(-2.5, -2, -1.566, -1.444, -1, -0.5, 0, 0.5, 1, 1.444, 1.566, 2, 2.5),
    values = cms.vdouble(
        1.019979, 0.9865564, 0.9876797, 1.002056, 0.9836735,
        0.9833507, 0.9866393, 0.988764, 1.018367, 0.9979675,
        1.256778, 0.986095, 0.9676724, 1.043812, 1.036677,
        1.030702, 0.987526, 0.9896907, 1.00818, 1.007179,
        0.9712154, 0.9917696, 0.9897751, 1.008155, 1.006104,
        0.9749728, 0.9896801, 0.9897541, 1.003055, 1.003043,
        0.9749728, 0.9968944, 0.9969136, 1.021561, 1.003043,
        0.9712154, 0.998968, 0.9979466, 1.019467, 1.006104,
        1.030702, 0.9875, 0.9896694, 1.010278, 1.007179,
        1.256778, 0.9731935, 0.9803922, 1.035332, 1.036677,
        0.9833507, 0.9876797, 0.9908164, 1.010183, 0.9979675,
        1.019979, 0.983556, 0.9877301, 1.007165, 0.9836735,
        ),
    errors = cms.vdouble(
        0.01505554, 0.001457952, 0.001448989, 0.01135466, 0.01194357,
        0.01778149, 0.00145271, 0.001443812, 0.01064251, 0.008682931,
        0.06373602, 0.002598568, 0.003261537, 0.02954206, 0.02847978,
        0.02117676, 0.001471606, 0.001459457, 0.006876683, 0.006880209,
        0.007151603, 0.001457201, 0.001448989, 0.005267125, 0.004549477,
        0.02101546, 0.001461719, 0.001451965, 0.004572736, 0.004535635,
        0.02101546, 0.001461719, 0.001451965, 0.004572736, 0.004535635,
        0.007151603, 0.001457201, 0.001448989, 0.00503263, 0.004549477,
        0.02117676, 0.001471606, 0.001459457, 0.006876683, 0.006880209,
        0.06373602, 0.002598568, 0.003468648, 0.02857233, 0.02847978,
        0.01778149, 0.00145271, 0.001443812, 0.009670574, 0.008682931,
        0.01505554, 0.001457952, 0.001448989, 0.01124703, 0.01194357,
        ),
)

electronSFCutBasedTight102X = cms.PSet(
    pt_bins = cms.vdouble(10, 20, 35, 50, 100, 200, 500,),
    eta_bins = cms.vdouble(-2.5, -2, -1.566, -1.444, -0.8, 0, 0.8, 1.444, 1.566, 2, 2.5,),
    values = cms.vdouble(
        1.203593, 1.162617, 1.091988, 1.063686, 1.027744, 0.9655938,
        1.022075, 1.006202, 1.010485, 1.0, 1.004592, 1.015625,
        1, 1, 1, 1, 1, 1,
        0.9944904, 0.9450172, 0.9519231, 0.9555273, 0.9642032, 0.9499444, 
        0.9713604, 0.9456522, 0.9522581, 0.9529554, 0.9564732, 1.026196,
        1.021687, 0.9596899, 0.9587097, 0.9588875, 0.9851258, 0.9636363,
        1.0, 0.9527972, 0.9583333, 0.9639176, 0.9964243, 0.9876266,
        1, 1, 1, 1, 1, 1,
        0.9889868, 1.011041, 0.9986911, 1.004932, 0.9942792, 0.9671262,
        1.268293, 1.141791, 1.084948, 1.037383, 1.008353, 1.059595,
        ),
    errors = cms.vdouble(
        0.04611582, 0.03274636, 0.01799809, 0.0266549, 0.03122403, 0.08148414,
        0.0476775, 0.03030152, 0.003983473, 0.007563576, 0.01912701, 0.05232746,
        1, 1, 1, 1, 1, 1,
        0.04182999, 0.02418594, 0.005900555, 0.006679659, 0.02871655, 0.04401765,
        0.02410042, 0.01665703, 0.00328969, 0.01147343, 0.01611903, 0.03158081,
        0.02410042, 0.01665703, 0.00328969, 0.01147343, 0.01611903, 0.02839926,
        0.04193757, 0.02418594, 0.005900555, 0.006679659, 0.02882864, 0.04396343,
        1, 1, 1, 1, 1, 1,
        0.0476775, 0.03030152, 0.003983473, 0.007730507, 0.01912701, 0.0488459,
        0.04611582, 0.03274636, 0.01799809, 0.0266549, 0.03098541, 0.08194443,
        ),
)

electronSFHLTZvtx102X = cms.PSet( #Dummy for 2018
    pt_bins = cms.vdouble(10, 1000,),
    eta_bins = cms.vdouble(-2.5, 2.5,),
    values = cms.vdouble(
        1.0,
        ),
    errors = cms.vdouble(
        0.0, 
        ),
)

#### Trigger SF
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs2017
#/afs/cern.ch/user/a/abdatta/public/SingleEG_JetHT_Trigger_Scale_Factors_ttHbb_Legacy2018_v1.root
trigSF_El32_El28HT150_ttH_legacy18_v1 = cms.PSet(
    pt_bins = cms.vdouble(30, 37, 45, 55, 100, 500),
    eta_bins = cms.vdouble(-2.5, -1.566, -1.4442, -0.8, 0, 0.8, 1.4442, 1.566, 2.5),
    values = cms.vdouble(
      0.7777987, 0.9242307, 0.9434291, 0.9495266, 0.9347277,
      1.0, 1.0, 1.0, 1.0, 1.0,
      0.900519, 0.9559557, 0.9584732, 0.9589031, 0.9675712,
      0.9105919, 0.9590567, 0.9653778, 0.9771409, 0.9836695,
      0.8792987, 0.9545439, 0.9532964, 0.9550337, 0.9572403,
      0.8400212, 0.9305611, 0.9508112, 0.9613527, 0.9652394,
      1.0, 1.0, 1.0, 1.0, 1.0,
      0.805282, 0.9581472, 0.9550291, 0.9627885, 0.944918,
    ),
    errors = cms.vdouble(
      0.03689664, 0.0253604, 0.02268788, 0.01329217, 0.01606752,
      1.0, 1.0, 1.0, 1.0, 1.0,
      0.02231461, 0.01762826, 0.01494909, 0.008288923, 0.0102015,
      0.01663414, 0.01414303, 0.01098092, 0.007149299, 0.008524632,
      0.01763328, 0.014016, 0.01151442, 0.008910196, 0.008249124,
      0.02272108, 0.01964398, 0.01577521, 0.009659816, 0.009887309,
      1.0, 1.0, 1.0, 1.0, 1.0,
      0.03402347, 0.02508284, 0.0201781, 0.01309658, 0.01788023,
    ),
)

#2018 before (< 316361) and after (>= 316361) HLT update. 
trigSF_IsoMu24A = cms.PSet(
    pt_bins = cms.vdouble(26, 30, 40, 50, 60, 120, 200, 300, 1200,),
    abseta_bins = cms.vdouble(0.0, 0.9, 1.2, 2.1, 2.4,),
    values = cms.vdouble(
      0.9698456, 0.9755923, 0.9765972, 0.9760637, 0.9751763, 0.9631842, 0.9618107, 0.9548403,
      0.9473602, 0.9566572, 0.9598887, 0.9599884, 0.9579104, 0.9551055, 0.9241194, 0.9619218,
      0.9685892, 0.9729108, 0.9738206, 0.9756057, 0.9760235, 0.9758267, 0.9664002, 0.9373623,
      0.9627936, 0.9956144, 1.010219, 1.009003, 1.02275, 0.989849, 1.0521, 1.025299,
    ),
    errors = cms.vdouble(
      0.001012076, 0.0003187262, 0.00023558, 0.0005508166, 0.0007756378, 0.007316778, 0.009773176, 0.01645118,
      0.002660559, 0.0006708286, 0.0004398361, 0.000930777, 0.001480966, 0.005692229, 0.01516908, 0.03259916,
      0.001870185, 0.0005944767, 0.0003954684, 0.0008655241, 0.001579233, 0.005378611, 0.01354485, 0.04598159,
      0.004382902, 0.001529621, 0.001165707, 0.002526876, 0.004402612, 0.02084879, 0.06024973, 0.139445,
    ),
)

trigSF_IsoMu24BCD = cms.PSet(
    pt_bins = cms.vdouble(26, 30, 40, 50, 60, 120, 200, 300, 1200,),
    abseta_bins = cms.vdouble(0.0, 0.9, 1.2, 2.1, 2.4,),
    values = cms.vdouble(
      0.9755329, 0.980781, 0.9810476, 0.9803042, 0.9786978, 0.9762971, 0.9823214, 0.9524137,
      0.9664956, 0.9749003, 0.9767722, 0.9752547, 0.9721834, 0.9712423, 0.9438655, 0.9384831,
      1.017329, 1.016138, 1.012924, 1.009952, 1.007131, 1.004202, 0.9930608, 1.000356,
      0.98497, 0.9992023, 1.006886, 1.00484, 1.015947, 1.000126, 0.994727, 0.9649939,
    ),
    errors = cms.vdouble(
      0.000640226, 0.0001956416, 0.0001451622, 0.0003306695, 0.0005331751, 0.006705951, 0.005111584, 0.01058117,
      0.001622344, 0.0003989414, 0.000250946, 0.0005373003, 0.0008404671, 0.003513026, 0.009799997, 0.019118,
      0.001207691, 0.0003863534, 0.0002508327, 0.0005667949, 0.0009652937, 0.003628991, 0.008000816, 0.01735078,
      0.002910666, 0.001012022, 0.0007806635, 0.001763289, 0.003272583, 0.02504356, 0.03431472, 0.09741847,
    ),
)

## Combined Scale factors
muonSFTight102X = combineSF(muonSFTightId102X, muonSFTightIso102X)
electronSFCutTight102X = combineSF(electronSFReco102X, electronSFCutBasedTight102X)
trigSF_IsoMu24 = computeAverageSF(trigSF_IsoMu24A, 8.951, trigSF_IsoMu24BCD, 50.811)
