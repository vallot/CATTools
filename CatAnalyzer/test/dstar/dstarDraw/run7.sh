#!/bin/bash

date

./run6.sh; ./wait.sh "noNotification"; echo "Nominal is done"

./run6.sh "noTopPtW"; ./wait.sh "noNotification"; echo "noTopPtW is done"
./run6.sh "TTnominal=TT_powheg_herwig"; ./wait.sh "noNotification"; echo "herwig is done"
./run6.sh "TTnominal=TT_powheg_mpiOFF"; ./wait.sh "noNotification"; echo "mpiOFF is done"
./run6.sh "TTnominal=TT_powheg_noCR"; ./wait.sh "noNotification"; echo "noCR is done"
./run6.sh "TTnominal=TT_powheg_scaledown"; ./wait.sh "noNotification"; echo "scaledown is done"
./run6.sh "TTnominal=TT_powheg_scaleup"; ./wait.sh "noNotification"; echo "scaleup is done"

./run6.sh "JES_Up"; ./wait.sh "noNotification"; echo "JES up is done"
./run6.sh "JES_Down"; ./wait.sh "noNotification"; echo "JES down is done"
./run6.sh "JER_Up"; ./wait.sh "noNotification"; echo "JER up is done"
./run6.sh "JER_Down"; ./wait.sh "noNotification"; echo "JER down is done"
./run6.sh "El_Up"; ./wait.sh "noNotification"; echo "Electron up is done"
./run6.sh "El_Down"; ./wait.sh "noNotification"; echo "Electron down is done"
./run6.sh "Mu_Up"; ./wait.sh "noNotification"; echo "Muon up is done"
./run6.sh "Mu_Down"; ./wait.sh "noNotification"; echo "Muon down is done"

./massPlot.py -s 0 -b [30,0,150] -p 'd0_lepSV_correctM' -x 'M_{l+D0} [GeV/c^2]' >> invmass/res_003.txt &
./massPlot.py -s 0 -b [30,0,150] -p 'dstar_lepSV_correctM' -x 'M_{l+D^{*}} [GeV/c^2]' >> invmass/res_130.txt &
./wait.sh
echo "Gen level is done"

date


