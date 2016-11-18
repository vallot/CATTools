#!/bin/bash

# For D0
echo "######## For D0 ########"

./run6_d0.sh
./run6_d0.sh "noTopPtW"
./run6_d0.sh "TTnominal=TT_powheg_herwig"
./run6_d0.sh "TTnominal=TT_powheg_mpiOFF"

./wait.sh "noNotification"
echo "Nominal is done"
echo "noTopPtW is done"
echo "herwig is done"
echo "mpiOFF is done"

./run6_d0.sh "TTnominal=TT_powheg_noCR"
./run6_d0.sh "TTnominal=TT_powheg_scaledown"
./run6_d0.sh "TTnominal=TT_powheg_scaleup"

./wait.sh "noNotification"
echo "noCR is done"
echo "scaledown is done"
echo "scaleup is done"

./run6_d0.sh "JES_Up"
./run6_d0.sh "JES_Down"
./run6_d0.sh "JER_Up"

./wait.sh "noNotification"
echo "JES up is done"
echo "JES down is done"
echo "JER up is done"

./run6_d0.sh "JER_Down"
./run6_d0.sh "El_Up"
./run6_d0.sh "El_Down"

./wait.sh "noNotification"
echo "JER down is done"
echo "Electron up is done"
echo "Electron down is done"

./run6_d0.sh "Mu_Up"
./run6_d0.sh "Mu_Down"

./massPlot.py -s 0 -b [30,0,150] -p 'd0_lepSV_correctM' -x 'M_{l+D0} [GeV/c^2]' >> invmass/res_003.txt &

./wait.sh
echo "Muon up is done"
echo "Muon down is done"
echo "Gen level is done"


# For D*
echo "######## For D* ########"
./run6_dstar.sh; ./wait.sh "noNotification"; echo "Nominal is done"

./run6_dstar.sh "noTopPtW"; ./wait.sh "noNotification"; echo "noTopPtW is done"
./run6_dstar.sh "TTnominal=TT_powheg_herwig"; ./wait.sh "noNotification"; echo "herwig is done"
./run6_dstar.sh "TTnominal=TT_powheg_mpiOFF"; ./wait.sh "noNotification"; echo "mpiOFF is done"
./run6_dstar.sh "TTnominal=TT_powheg_noCR"; ./wait.sh "noNotification"; echo "noCR is done"
./run6_dstar.sh "TTnominal=TT_powheg_scaledown"; ./wait.sh "noNotification"; echo "scaledown is done"
./run6_dstar.sh "TTnominal=TT_powheg_scaleup"; ./wait.sh "noNotification"; echo "scaleup is done"

./run6_dstar.sh "JES_Up"; ./wait.sh "noNotification"; echo "JES up is done"
./run6_dstar.sh "JES_Down"; ./wait.sh "noNotification"; echo "JES down is done"
./run6_dstar.sh "JER_Up"; ./wait.sh "noNotification"; echo "JER up is done"
./run6_dstar.sh "JER_Down"; ./wait.sh "noNotification"; echo "JER down is done"
./run6_dstar.sh "El_Up"; ./wait.sh "noNotification"; echo "Electron up is done"
./run6_dstar.sh "El_Down"; ./wait.sh "noNotification"; echo "Electron down is done"
./run6_dstar.sh "Mu_Up"; ./wait.sh "noNotification"; echo "Muon up is done"
./run6_dstar.sh "Mu_Down"; ./wait.sh "noNotification"; echo "Muon down is done"

./massPlot.py -s 0 -b [30,0,150] -p 'dstar_lepSV_correctM' -x 'M_{l+D^{*}} [GeV/c^2]' >> invmass/res_130.txt &
./wait.sh
echo "Gen level is done"


