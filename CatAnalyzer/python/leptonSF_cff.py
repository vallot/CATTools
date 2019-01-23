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

##Muon SF reference for 2017: https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2017
muonSFTightIdOnly94X  = cms.PSet(
    pt_bins = cms.vdouble(20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 120.0,),
    abseta_bins = cms.vdouble(0.0, 0.9, 1.2, 2.1, 2.4,),
    values = cms.vdouble(
      0.99107776277569515, 0.98741046826208401, 0.9907753279135898, 0.98924835889520468, 0.98555451603347632, 0.98980573770933888,
      0.99273892755152437, 0.98506393976251205, 0.98653594641822473, 0.984913093101493, 0.98390563847600077, 0.98406040314346799,
      0.99242527198773844, 0.98908844612849334, 0.99464690698838409, 0.99265288251551831, 0.99063642229435289, 0.99204643221439792,
      0.9758095839531763, 0.97451535941798839, 0.97874105001587464, 0.978189122919501, 0.96735684160978941, 0.97663118567312024,
    ),
    errors = cms.vdouble(
      0.00486221, 0.00234205, 0.00172048, 0.000553118, 0.00278313, 0.00337485,
      0.0584377, 0.021916, 0.00234107, 0.020152, 0.00168673, 0.0124799,
      0.00779607, 0.0149208, 0.0126635, 0.00995829, 0.00112734, 0.00281051,
      0.0045305, 0.00298196, 0.00328128, 0.0019898, 0.00425963, 0.00957463,
    ),
)
#"NUM_TightRelIso_DEN_TightIDandIPCut" for 2017
muonSFTightIsoOnly94X = cms.PSet(
    pt_bins = cms.vdouble(20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 120.0,),
    abseta_bins = cms.vdouble(0.0, 0.9, 1.2, 2.1, 2.4,),
    values = cms.vdouble(
      0.99278134880281932, 0.99733223261925152, 0.99757673354437448, 0.99821234048998209, 0.99844669557462695, 0.9995299697442489,
      0.99573295192897304, 0.98896867902501906, 0.99473797534413366, 0.99701599571535326, 0.99934624944433836, 1.0000177267113954,
      0.9914010546323393, 0.99406676632373769, 0.99552377293660199, 0.99746791885529507, 0.99859688733839247, 0.99964911760449982,
      0.98882971036295109, 0.9923197624578961, 0.99499586727665668, 0.99792605855678485, 0.99936392435049604, 0.99834681165529127,
    ),
    errors = cms.vdouble(
      0.00377464, 0.00190303, 0.000742077, 0.000281225, 0.000572324, 0.000640803,
      0.00650588, 0.00334989, 0.00103829, 0.000419564, 0.0018559, 0.00131999,
      0.00273885, 0.00148174, 0.000557557, 0.000231935, 0.000573589, 0.000714106,
      0.00372993, 0.00201685, 0.00079892, 0.000413348, 0.000897969, 0.0014165,
    ),
)


## ELectron SF for 2017: https://twiki.cern.ch/twiki/bin/view/CMS/Egamma2017DataRecommendations#Efficiency_Scale_Factors
electronSFMVAWP80RecoOnly94Xv2 = cms.PSet(
    pt_bins = cms.vdouble(20, 45, 75, 100, 500),
    eta_bins = cms.vdouble(-2.5, -2, -1.566, -1.444, -1, -0.5, 0, 0.5, 1, 1.444, 1.566, 2, 2.5),
    values = cms.vdouble(
        0.9886246, 0.9846154, 1.00102, 1.007157,
        0.9907503, 0.9908069, 1.006091, 0.9918946,
        0.9815242, 0.9590958, 1.046739, 0.9837487,
        0.987526, 0.9886715, 1.005128, 1.001032,
        0.9897119, 0.9907975, 1.002039, 1.001019,
        0.9855521, 0.9887295, 1.006141, 0.9868554,
        0.983454, 0.9866255, 1.006141, 0.9868554,
        0.9865702, 0.9886948, 1.002039, 1.001019,
        0.9843587, 0.982438, 1.005128, 1.001032,
        0.9848131, 0.9727075, 1.046739, 0.9837487,
        0.988718, 0.9908164, 1.006091, 0.9918946,
        0.9917611, 0.9857143, 1.00102, 1.007157,
        ),
    errors = cms.vdouble(
        0.001459457, 0.001446766, 0.0065338, 0.009646197,
        0.001451965, 0.001443812, 0.005076142, 0.006487461,
        0.002859715, 0.004136677, 0.02200919, 0.03426086,
        0.001472372, 0.001458704, 0.005128205, 0.007441798,
        0.001457952, 0.001449732, 0.003675384, 0.00509684,
        0.001460964, 0.001451965, 0.003690431, 0.004289829,
        0.001460964, 0.001451965, 0.003690431, 0.004289829,
        0.001457952, 0.001449732, 0.003675384, 0.00509684,
        0.001472372, 0.001458704, 0.005128205, 0.007441798,
        0.003127273, 0.004369787, 0.02200919, 0.03426086,
        0.001451965, 0.001443812, 0.005076142, 0.006487461,
        0.001459457, 0.001446766, 0.0065338, 0.009646197,
        ),
)

electronSFMVAWP80IsoIdOnly94Xv2 = cms.PSet(
    pt_bins = cms.vdouble(10, 20, 35, 50, 100, 200, 500,),
    eta_bins = cms.vdouble(-2.5, -2, -1.566, -1.444, -0.8, 0, 0.8, 1.444, 1.566, 2, 2.5,),
    values = cms.vdouble(
        0.89441, 0.8888889, 0.8977833, 0.9032648, 0.9578164, 0.9098143,
        0.9239131, 0.9116915, 0.926097, 0.9430168, 0.9656319, 0.9230769,
        1, 1, 1, 1, 1, 1,
        0.9851412, 0.953271, 0.9446495, 0.9463647, 0.9579832, 0.9485294,
        0.976923, 0.977027, 0.9625935, 0.9600484, 0.9668639, 0.962877,
        0.969697, 0.9746667, 0.9594096, 0.958134, 0.9706573, 0.9110867,
        0.9818731, 0.9567568, 0.9479554, 0.9474313, 0.9733334, 0.9311981,
        1, 1, 1, 1, 1, 1,
        0.9497965, 0.9216418, 0.9295612, 0.9455556, 0.9647577, 0.9093887,
        0.9029734, 0.8897849, 0.893696, 0.88929, 0.9271845, 0.9591584,
        ),
    errors = cms.vdouble(
        0.04670853, 0.02164992, 0.01013675, 0.01063642, 0.03075603, 0.121146,
        0.03299417, 0.0261194, 0.01042464, 0.008503809, 0.08210882, 0.0842841,
        1, 1, 1, 1, 1, 1,
        0.09301971, 0.03423142, 0.0124226, 0.008687482, 0.01716556, 0.1273453,
        0.0205824, 0.02219844, 0.00672618, 0.009398616, 0.02604336, 0.1205491,
        0.0205824, 0.02219844, 0.00672618, 0.009398616, 0.02596826, 0.1203642,
        0.09301971, 0.03423142, 0.0124226, 0.008687482, 0.01716556, 0.1274177,
        1, 1, 1, 1, 1, 1,
        0.03299417, 0.0261194, 0.01042464, 0.008503809, 0.08210882, 0.08493492,
        0.04664833, 0.02164992, 0.01013675, 0.01063642, 0.03087054, 0.1217787,
        ),
)

electronSFCutBasedTightIDOnly94Xv2 = cms.PSet(
    pt_bins = cms.vdouble(10, 20, 35, 50, 100, 200, 500,),
    eta_bins = cms.vdouble(-2.5, -2, -1.566, -1.444, -0.8, 0, 0.8, 1.444, 1.566, 2, 2.5,),
    values = cms.vdouble(
        0.9725686, 0.9276094, 0.9196051, 0.9136598, 0.9892473, 0.9257991,
        0.9592668, 0.9583333, 0.9547803, 0.9613993, 0.9740991, 0.9724062,
        1, 1, 1, 1, 1, 1,
        1.056995, 0.9880952, 0.9607843, 0.9591837, 0.9800704, 0.9520174,
        1.008791, 0.990712, 0.9658344, 0.9563636, 0.976082, 0.9576547,
        0.9869565, 0.9847328, 0.9673629, 0.9613993, 0.976298, 0.9599133,
        1.031088, 0.9897436, 0.9661972, 0.9568528, 0.9708624, 0.9575626,
        1, 1, 1, 1, 1, 1,
        0.9674134, 0.958209, 0.9498069, 0.9521531, 0.9920635, 0.9482577,
        0.9830918, 0.9448161, 0.9230769, 0.9247449, 0.9569266, 1.08918,
        ),
    errors = cms.vdouble(
        0.08555808, 0.02676687, 0.01428862, 0.01875479, 0.02069147, 0.1037975,
        0.05818755, 0.03030507, 0.01570057, 0.01069339, 0.05700198, 0.1181784,
        1, 1, 1, 1, 1, 1,
        0.06833439, 0.0364394, 0.01733, 0.009541985, 0.02559112, 0.05700315,
        0.04977479, 0.02316144, 0.01568978, 0.01550877, 0.04159333, 0.03288074,
        0.04949505, 0.02316144, 0.01568978, 0.01550877, 0.04157304, 0.03293771,
        0.06833439, 0.0364394, 0.01733, 0.009541985, 0.02563613, 0.05709869,
        1, 1, 1, 1, 1, 1,
        0.05818755, 0.03030507, 0.01570057, 0.01069339, 0.05702298, 0.1176582,
        0.08549427, 0.02676687, 0.01428862, 0.01875479, 0.02058619, 0.1061907,
        ),
)

electronSFHLTZvtx94X = cms.PSet(
    pt_bins = cms.vdouble(10, 500,),
    eta_bins = cms.vdouble(-2.5, 2.5,),
    values = cms.vdouble(
        0.991,
        ),
    errors = cms.vdouble(
        0.001, 
        ),
)

#### Trigger SF
#### 2017 First
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs2017
trigSF_El35_El28HT150_ttH_v5 = cms.PSet(
    pt_bins = cms.vdouble(30, 37, 45, 55, 100, 500),
    eta_bins = cms.vdouble(-2.5, -1.566, -1.4442, -0.8, 0, 0.8, 1.4442, 1.566, 2.5),
    values = cms.vdouble(
      0.7139729, 0.8447125, 0.9022676, 0.9130857, 0.9120451,
      1.0, 1.0, 1.0, 1.0, 1.0,
      0.7392274, 0.9299934, 0.9216896, 0.9450684, 0.9684056,
      0.8691187, 0.9304942, 0.9576964, 0.9394679, 0.9587997,
      0.8331928, 0.9234828, 0.9195148, 0.9293138, 0.9373967,
      0.7717557, 0.9024622, 0.9291959, 0.9294295, 0.9403721,
      1.0, 1.0, 1.0, 1.0, 1.0,
      0.6207152, 0.8319214, 0.8981477, 0.8827884, 0.8901988,
    ),
    errors = cms.vdouble(
      0.04828459, 0.04693317, 0.03059252, 0.01935367, 0.02818359,
      1.0, 1.0, 1.0, 1.0, 1.0,
      0.04742585, 0.03061207, 0.02393724, 0.01672756, 0.01694731,
      0.02420729, 0.02000441, 0.03078572, 0.008935409, 0.01075451,
      0.02210147, 0.0179609, 0.01646352, 0.009328098, 0.01109023,
      0.05940178, 0.06457422, 0.04676345, 0.02045333, 0.01399639,
      1.0, 1.0, 1.0, 1.0, 1.0,
      0.08215168, 0.03667785, 0.03386165, 0.01843878, 0.02542128,
    ),
)

trigSF_IsoMu27 = cms.PSet(
    pt_bins = cms.vdouble(29, 32, 40, 50, 60, 120, 200, 1200,),
    abseta_bins = cms.vdouble(0.0, 0.9, 1.2, 2.1, 2.4,),
    values = cms.vdouble(
      0.9704915, 0.9728889, 0.9721969, 0.9719104, 0.9705052, 0.969812, 0.9626321,
      0.9515225, 0.9544191, 0.9525473, 0.9518248, 0.9449776, 0.9343395, 0.9266601,
      0.9900779, 0.9965652, 0.9939828, 0.9901482, 0.9873659, 0.9921802, 0.9828429,
      0.8622288, 0.9108149, 0.9317348, 0.9442655, 0.9484941, 0.9726911, 0.8435198,
    ),
    errors = cms.vdouble(
      0.0007590597, 0.0002695662, 0.0001989391, 0.0004486466, 0.000878429, 0.003061058, 0.005976491,
      0.001908455, 0.0005574671, 0.000354451, 0.0007887377, 0.001520861, 0.00551032, 0.01390627,
      0.001663918, 0.000554091, 0.0003458035, 0.0007651448, 0.001404214, 0.005650038, 0.01995085,
      0.003523487, 0.001257815, 0.00091109, 0.002090694, 0.004220159, 0.02028692, 0.04057829,
    ),
)

## Combined Scale factors
#2017 SF - ongoing
muonSFTight94X = combineSF(muonSFTightIdOnly94X, muonSFTightIsoOnly94X)
electronSFMVAWP8094Xv2 = combineSF( combineSF( electronSFMVAWP80IsoIdOnly94Xv2, electronSFMVAWP80RecoOnly94Xv2 ), electronSFHLTZvtx94X)
