import FWCore.ParameterSet.Config as cms

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

## Muon SF reference https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults
## SF for Run2016BCDE, before HIP issue
muonSFTrackingOnly = cms.PSet(
    eta_bins = cms.vdouble(-2.400000, -2.100000, -1.600000,-1.200000,-0.900000,-0.600000,-0.300000,-0.200000,0.200000,0.300000,0.600000,0.900000,1.200000,1.600000,2.100000,2.400000,),
    values = cms.vdouble(
        1.000085,
        0.997879,
        0.998387,
        0.998666,
        0.999174,
        0.999416,
        0.999280,
        0.999405,
        0.999678,
        0.999544,
        0.999311,
        0.998382,
        0.997051,
        0.998565,
        0.998243,
    ),
    errors = cms.vdouble(
        0.000871,
        0.000246,
        0.000206,
        0.000227,
        0.000111,
        0.000106,
        0.000182,
        0.000081,
        0.000190,
        0.000107,
        0.000118,
        0.000227,
        0.000218,
        0.000226,
        0.000849,
    ),
)

## SF for Run2016G-H, after HIP issue
muonSFTrackingGHOnly = cms.PSet(
    eta_bins = cms.vdouble(-2.400000, -2.100000, -1.600000,-1.200000,-0.900000,-0.600000,-0.300000,-0.200000,0.200000,0.300000,0.600000,0.900000,1.200000,1.600000,2.100000,2.400000,),
    values = cms.vdouble(
        0.982741,
        0.992149,
        0.994382,
        0.995412,
        0.995896,
        0.995767,
        0.994363,
        0.994700,
        0.995835,
        0.997718,
        0.997370,
        0.996884,
        0.993423,
        0.991586,
        0.976802,
    ),
    errors = cms.vdouble(
        0.001011,
        0.000252,
        0.000227,
        0.000266,
        0.000123,
        0.000119,
        0.000233,
        0.000099,
        0.000224,
        0.000109,
        0.000120,
        0.000255,
        0.000236,
        0.000241,
        0.001033,
    ),
)


## SF for Run2016BCDE, before HIP issue
muonSFTightIdOnly = cms.PSet(
    pt_bins = cms.vdouble(20.000000, 25.000000, 30.000000, 40.000000, 50.000000, 60.000000, 120.000000,),
    abseta_bins = cms.vdouble(0.000000, 0.900000, 1.200000, 2.100000, 2.400000,),
    values = cms.vdouble(
        0.980779, 0.979747, 0.981756, 0.982723, 0.979586, 0.992805,
        0.961656, 0.959323, 0.965162, 0.967988, 0.969637, 0.967575,
        0.982584, 0.982259, 0.984453, 0.987816, 0.985261, 0.988935,
        0.970229, 0.969708, 0.967787, 0.970770, 0.967764, 0.963107,
    ),
    errors = cms.vdouble(
        0.002778, 0.001092, 0.000380, 0.000290, 0.000900, 0.002036,
        0.003445, 0.002070, 0.002116, 0.000499, 0.001376, 0.002804,
        0.001669, 0.000842, 0.000324, 0.000274, 0.000907, 0.002019,
        0.003056, 0.001790, 0.000959, 0.003717, 0.002790, 0.005782,
    ),
)

## SF for Run2016G-H, after HIP issue
muonSFTightGHIdOnly = cms.PSet(
    pt_bins = cms.vdouble(20.000000, 25.000000, 30.000000, 40.000000, 50.000000, 60.000000, 120.000000,),
    abseta_bins = cms.vdouble(0.000000, 0.900000, 1.200000, 2.100000, 2.400000,),
    values = cms.vdouble(
        0.993173, 0.986990, 0.987596, 0.989777, 0.984749, 0.991370,
        0.985596, 0.984686, 0.983914, 0.983265, 0.980582, 0.983879,
        0.990863, 0.990917, 0.992066, 0.993847, 0.985655, 0.988584,
        0.981501, 0.979109, 0.971526, 0.974776, 0.967651, 0.963199,
    ),
    errors = cms.vdouble(
        0.002642, 0.046316, 0.000369, 0.000300, 0.000785, 0.001626,
        0.015529, 0.001699, 0.000600, 0.000467, 0.001423, 0.002968,
        0.001871, 0.000941, 0.020594, 0.000274, 0.001019, 0.002380,
        0.003092, 0.053287, 0.002042, 0.006083, 0.002817, 0.005718,
    ),
)
## SF for Run2016BCDE, before HIP issue
muonSFTightIsoOnly = cms.PSet(
    pt_bins = cms.vdouble(20.000000, 25.000000, 30.000000, 40.000000, 50.000000, 60.000000, 120.000000,),
    abseta_bins = cms.vdouble(0.000000, 0.900000, 1.200000, 2.100000, 2.400000,),
    values = cms.vdouble(
        0.987723, 0.993733, 0.994052, 0.995378, 0.996878, 0.998548,
        0.993162, 1.001004, 0.999295, 0.997179, 0.999354, 0.999297,
        0.989717, 0.994223, 0.997134, 0.997531, 0.997972, 0.999017,
        0.975341, 0.985961, 0.993083, 0.996770, 0.997456, 1.001517,
    ),
    errors = cms.vdouble(
        0.002236, 0.001190, 0.000385, 0.000215, 0.000363, 0.000463,
        0.003570, 0.002053, 0.000699, 0.000005, 0.000599, 0.000771,
        0.001684, 0.000995, 0.000381, 0.000159, 0.000338, 0.000456,
        0.002945, 0.001762, 0.000687, 0.000437, 0.000766, 0.001126,
    ),
)

## SF for Run2016G-H, after HIP issue
muonSFTightGHIsoOnly = cms.PSet(
    pt_bins = cms.vdouble(20.000000, 25.000000, 30.000000, 40.000000, 50.000000, 60.000000, 120.000000,),
    abseta_bins = cms.vdouble(0.000000, 0.900000, 1.200000, 2.100000, 2.400000,),
    values = cms.vdouble(
        0.981144, 0.992791, 0.993452, 0.995200, 0.996716, 0.999064,
        0.997662, 0.999726, 0.999598, 0.998324, 0.998887, 0.998908,
        0.993687, 0.997960, 0.998980, 0.998603, 0.998800, 0.999453,
        0.994233, 0.999098, 1.000063, 1.000221, 1.000003, 1.001529,
    ),
    errors = cms.vdouble(
        0.002283, 0.001220, 0.000396, 0.000130, 0.000374, 0.000472,
        0.003563, 0.002081, 0.000705, 0.000300, 0.000611, 0.000794,
        0.001676, 0.001003, 0.000380, 0.000008, 0.000340, 0.000461,
        0.002858, 0.001746, 0.000692, 0.000260, 0.000787, 0.001123,
    ),
)


## Muon SF reference https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffsRun2
## Retrieve data from the cmsdoc web page:
##   https://cmsdoc.cern.ch/cms/Physics/muon/ReferenceEfficiencies/Run2015/25ns/MuonID_Z_RunCD_Reco76X_Feb15.json
##   https://cmsdoc.cern.ch/cms/Physics/muon/ReferenceEfficiencies/Run2015/25ns/MuonIso_Z_RunCD_Reco76X_Feb15.json
muonSFTight76X = cms.PSet(
    # Values of "MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1 + MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1"
    pt_bins = cms.vdouble(20, 25, 30, 40, 50, 60, 120),
    abseta_bins = cms.vdouble(0, 0.9, 1.2, 2.1, 2.4),
    values = cms.vdouble(
        0.975283, 0.980129, 0.986676, 0.985789, 0.980877, 0.985445,
        0.984484, 0.970783, 0.980673, 0.980057, 0.977943, 0.984907,
        0.986784, 0.990845, 0.993367, 0.990324, 0.988209, 0.991915,
        0.971296, 0.971883, 0.976009, 0.977, 0.973655, 0.963135,
    ),
    errors = cms.vdouble(
        0.0157468, 0.0136165, 0.0120878, 0.0118247, 0.0125879, 0.0135921,
        0.0181504, 0.0154891, 0.0127943, 0.0121561, 0.0134922, 0.0152142,
        0.0145121, 0.0132486, 0.0120305, 0.0116593, 0.0124603, 0.0137997,
        0.0176168, 0.0153647, 0.0131657, 0.0126734, 0.0148934, 0.0197265,
    ),

)

## Electron SF from https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2#Efficiencies_and_scale_factors
## root file is retrieved from 
##  Reco: http://fcouderc.web.cern.ch/fcouderc/EGamma/scaleFactors/Moriond17/approval/RECO/passingRECO/egammaEffi.txt_EGM2D.root 
##   ID: http://fcouderc.web.cern.ch/fcouderc/EGamma/scaleFactors/Moriond17/approval/EleID/passingMedium80X/egammaEffi.txt_EGM2D.root 
electronSFRecoOnly = cms.PSet(
    pt_bins = cms.vdouble(25.000000, 500.000000,),
    eta_bins = cms.vdouble(-2.500000, -2.450000, -2.400000, -2.300000, -2.200000, -2.000000, -1.800000, -1.630000, -1.566000, -1.444000, -1.200000, -1.000000, -0.600000, -0.400000, -0.200000, 0.000000, 0.200000, 0.400000, 0.600000, 1.000000, 1.200000, 1.444000, 1.566000, 1.630000, 1.800000, 2.000000, 2.200000, 2.300000, 2.400000, 2.450000, 2.500000,),
    values = cms.vdouble(
        1.317604,
        1.113780,
        1.024625,
        1.013641,
        1.007277,
        0.994819,
        0.994786,
        0.991632,
        0.963129,
        0.989701,
        0.985656,
        0.981595,
        0.984678,
        0.981614,
        0.980433,
        0.984552,
        0.988764,
        0.987743,
        0.987743,
        0.987743,
        0.987680,
        0.967598,
        0.989627,
        0.992761,
        0.991761,
        0.997940,
        1.001037,
        0.989507,
        0.970519,
        0.907,
        ),
    errors = cms.vdouble(
        0.018239,
        0.011067,
        0.008158,
        0.007133,
        0.004203,
        0.006493,
        0.005166,
        0.005512,
        0.026030,
        0.003599,
        0.005064,
        0.003312,
        0.006129,
        0.006358,
        0.005302,
        0.005302,
        0.006358,
        0.006129,
        0.003312,
        0.005064,
        0.003599,
        0.026030,
        0.005512,
        0.005166,
        0.006493,
        0.004203,
        0.007133,
        0.008158,
        0.011067,
        0.018239,
        ),
)
electronSFCutBasedIDMediumWPIdOnly = cms.PSet(
    pt_bins = cms.vdouble(10.000000,20.000000,35.000000,50.000000,90.000000,150.000000,500.000000),
    eta_bins = cms.vdouble(-2.500000,-2.000000,-1.566000,-1.444000,-0.800000,0.000000,0.800000,1.444000,1.566000,2.000000,2.500000),
    values = cms.vdouble(
        0.820574, 0.914201, 0.955919, 0.974940, 1.040745, 1.040745,
        0.822650, 0.943662, 0.979369, 0.991908, 1.018059, 1.018059,
        1.027119, 1.000000, 0.985653, 0.987838, 1.086614, 1.086614,
        0.959108, 0.958556, 0.969767, 0.972129, 0.993363, 0.993363,
        0.922027, 0.946164, 0.956725, 0.958427, 0.985588, 0.985588,
        0.939571, 0.971467, 0.979882, 0.979569, 1.018038, 1.018038,
        0.966480, 0.972973, 0.976608, 0.982022, 1.018059, 1.018059,
        0.993399, 0.966270, 0.970000, 0.986395, 1.015605, 1.015605,
        0.853659, 0.930496, 0.973430, 0.985006, 0.997750, 0.997750,
        0.815385, 0.891273, 0.942928, 0.967098, 1.015222, 1.015222,
        ),
    errors = cms.vdouble(
        0.023166, 0.008260, 0.007552, 0.005656, 0.029675, 0.029866,
        0.018241, 0.012760, 0.005869, 0.005066, 0.011225, 0.011722,
        0.068255, 0.179759, 0.012664, 0.014144, 0.035429, 0.035590,
        0.050900, 0.011640, 0.004441, 0.012475, 0.008380, 0.009034,
        0.026763, 0.012595, 0.004033, 0.011094, 0.010382, 0.010917,
        0.026763, 0.012595, 0.004033, 0.011094, 0.010382, 0.010917,
        0.050900, 0.011640, 0.004441, 0.012475, 0.008380, 0.009034,
        0.068255, 0.179759, 0.012664, 0.014221, 0.035289, 0.035450,
        0.018241, 0.012760, 0.005869, 0.005066, 0.011225, 0.011722,
        0.023166, 0.008260, 0.007552, 0.005491, 0.029675, 0.029866,
        ),
)
## Electorn SF from http://fcouderc.web.cern.ch/fcouderc/EGamma/scaleFactors/moriond2016_76X/eleCutBasedID/CutBasedID_MediumWP_76X_18Feb.txt_egammaPlots.pdf
electronSFCutBasedIDMediumWP76X = cms.PSet(
    # Values of "CutBasedID_MediumWP"
    pt_bins = cms.vdouble(10, 20, 30, 40, 50, 200),
    abseta_bins = cms.vdouble(0, 0.8, 1.444, 1.566, 2, 2.5),
    values = cms.vdouble(
        1.00726, 0.971182, 0.984906, 0.985899, 0.988598,
        1.0902, 0.983359, 0.987179, 0.98692, 0.986159,
        1.08642, 0.963054, 0.949123, 0.981612, 0.997257,
        0.984444, 0.936809, 0.975066, 0.992806, 1.00579,
        1.03557, 0.986446, 0.963351, 1.00611, 1.00949,
    ),
    errors = cms.vdouble(
        0.0656124, 0.0279776, 0.0291216, 0.0175478, 0.00603364,
        0.0391175, 0.0256741, 0.0171527, 0.0277844, 0.00540994,
        0.0800687, 0.261884 , 0.0511487, 0.036227 , 0.113466  ,
        0.0447765, 0.0186922, 0.0274652, 0.0404665, 0.0125727 ,
        0.0303603, 0.0263447, 0.0445603, 0.0103732, 0.0137829 ,
    ),
)

## Retrieve data from the cmsdoc web page:
##   https://cmsdoc.cern.ch/cms/Physics/muon/ReferenceEfficiencies/Run2015/25ns/MuonID_Z_RunCD_Reco74X_Dec1.json
muonSFTight74X = cms.PSet(
    # Values of "NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1 + NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1"
    pt_bins = cms.vdouble(20, 25, 30, 40, 50, 60, 120),
    abseta_bins = cms.vdouble(0, 0.9, 1.2, 2.1, 2.4),
    values = cms.vdouble(
        0.979479, 0.984375, 0.987068, 0.986086, 0.983674, 0.985021,
        0.978054, 0.980927, 0.984458, 0.981568, 0.975983, 0.984876,
        0.99541, 0.991193, 0.994704, 0.992697, 0.988925, 0.995536,
        0.983608, 0.974278, 0.976037, 0.980598, 0.969666, 0.969941,
    ),
    errors = cms.vdouble(
        0.0191201, 0.0168928, 0.0151444, 0.0148111, 0.015527, 0.0162596,
        0.021688, 0.0189092, 0.0159175, 0.0150692, 0.0164092, 0.0178394,
        0.0179975, 0.0166027, 0.0151404, 0.0147407, 0.0154785, 0.0167413,
        0.0219747, 0.0194149, 0.0165014, 0.0159522, 0.0180924, 0.0243999,
    ),
)

electronSFWP9074X = cms.PSet(
    ## SF from https://indico.cern.ch/event/369259/contribution/3/attachments/1204731/1755463/update_SFs_triggering_MVA_ID_systematics.pdf
    eta_bins = cms.vdouble(-2.5, -1.5, -1.0, 0, 1.0, 1.5, 2.5),
    pt_bins = cms.vdouble(15, 25, 35, 45, 55, 1e9),
    values = cms.vdouble(
        0.96, 0.95, 0.98, 0.99, 0.99, 0.97,
        0.98, 0.97, 0.97, 0.99, 0.99, 0.98,
        0.98, 0.99, 0.99, 0.99, 0.99, 0.98,
        0.98, 0.99, 0.99, 0.99, 0.99, 0.99,
        0.99, 0.99, 1.00, 1.00, 0.99, 0.99,
    ),
    errors = cms.vdouble(
        0.02, 0.06, 0.01, 0.01, 0.02, 0.01,
        0.01, 0.01, 0.02, 0.01, 0.01, 0.01,
        0.00, 0.00, 0.00, 0.00, 0.00, 0.01,
        0.01, 0.00, 0.01, 0.00, 0.01, 0.01,
        0.01, 0.02, 0.01, 0.00, 0.02, 0.02,
    ),
)

electronSFCutBasedIDMediumWP74X = cms.PSet(
    ## SF from https://indico.cern.ch/event/491528/contribution/2/attachments/1231399/1805319/CutBasedId_EGM_19Feb.pdf
    ## Actual numbers are taken from https://arun-afs.web.cern.ch/arun-afs/Final_Fits_Data_74X/CutBasedID_MediumWP_fromTemplates_withSyst_Final.txt
    # Values of "CutBasedID_MediumWP"
    pt_bins = cms.vdouble(10, 20, 30, 40, 50, 200),
    eta_bins = cms.vdouble(-2.5, -2, -1.566, -1.444, -0.8, 0, 0.8, 1.444, 1.566, 2, 2.5),
    values = cms.vdouble(
        1.25131, 1.0536, 1.01235, 1.00875, 1.006,
        1.09254, 0.969453, 0.971279, 0.983313, 0.985109,
        1.29412, 1.04043, 1.02482, 0.990155, 0.948856,
        1.32673, 1.05939, 0.990897, 0.974057, 0.958097,
        1.23371, 1.02262, 0.978616, 0.969838, 0.967489,
        1.22889, 1.01799, 0.986164, 0.974448, 0.970819,
        1.30827, 1.07742, 1.00654, 0.976359, 0.964773,
        1.35498, 1.18801, 0.990991, 0.956153, 0.974394,
        1.06283, 0.980392, 0.961892, 0.978469, 0.994253,
        1.20769, 1.02649, 1.0096, 1.00125, 1.00481,
    ),
    errors = cms.vdouble(
        0.712863, 0.0517635, 0.0242599, 0.0222102, 0.0337348,
        0.662784, 0.0757762, 0.575897, 0.0392023, 0.0671983,
        3.75487, 1.17119, 0.124767, 0.0150145, 0.0555887,
        1.2277, 0.20447, 0.0280944, 0.00855398, 0.0355446,
        0.813571, 0.108281, 0.00853765, 0.0217419, 0.0302874,
        0.761126, 0.580145, 0.0313306, 0.0196618, 0.0253676,
        0.874693, 0.237932, 0.650881, 0.0103168, 0.0286228,
        4.59432, 1.84669, 0.1764, 0.0506619, 0.161535,
        0.797724, 0.0632952, 0.0377029, 0.0313089, 0.0300417,
        0.617985, 0.0693681, 0.0293348, 0.024152, 0.0894748,
    ),
)

## Combined Scale factors
## id syst 1%+0.5%(quadrature), iso syst 1%
muonSFTight = combineSF(muonSFTightIdOnly, muonSFTightIsoOnly, (0.01**2+0.005**2)**0.5, 0.01)
muonSFTightGH = combineSF(muonSFTightGHIdOnly, muonSFTightGHIsoOnly, (0.01**2+0.005**2)**0.5, 0.01)
electronSFCutBasedIDMediumWP = combineSF(electronSFRecoOnly, electronSFCutBasedIDMediumWPIdOnly)

