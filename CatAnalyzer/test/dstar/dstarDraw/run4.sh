#!/bin/bash
eval `scramv1 runtime -sh`

cut='(abs(d0_dRTrue)<0.1&&abs(d0_relPtTrue)<0.1)'
acut='(abs(d0_dRTrue)>0.1||abs(d0_relPtTrue)>0.1)'

cut2='(abs(dstar_dRTrue)<0.1&&abs(dstar_relPtTrue)<0.1)'
acut2='(abs(dstar_dRTrue)>0.1||abs(dstar_relPtTrue)>0.1)'
./compare.py -c "${cut}"  -a "${acut}" -s 5 -b [60,0,0.1] -p 'd0_dca' -x 'DCA of D0 Cands' -f 'compare' & 
./compare.py -c "${cut2}"  -a "${acut2}" -s 5 -b [60,0,0.1] -p 'dstar_dca' -x 'DCA of D* Cands' -f 'compare'> /dev/null &
./compare.py -c "${cut}"  -a "${acut}" -s 5 -b [60,0,0.3] -p 'd0_vProb' -x 'Vertex Probability of D0 Cands' -f 'compare'> /dev/null &
./compare.py -c "${cut2}"  -a "${acut2}" -s 5 -b [60,0,0.3] -p 'dstar_vProb' -x 'Vertex Probability of D* Cands' -f 'compare'> /dev/null &
./compare.py -c "${cut}"  -a "${acut}" -s 5 -b [60,0,0.5] -p 'd0_LXY' -x 'Lxy of D0 Cands' -f 'compare'> /dev/null &
./compare.py -c "${cut2}"  -a "${acut2}" -s 5 -b [60,0,0.5] -p 'dstar_LXY' -x 'Lxy of D* Cands' -f 'compare'> /dev/null &
./compare.py -c "${cut}"  -a "${acut}" -s 5 -b [60,0,1] -p 'd0_L3D' -x 'L_{3D} of D0 Cands' -f 'compare'> /dev/null &
./compare.py -c "${cut2}"  -a "${acut2}" -s 5 -b [60,0,1] -p 'dstar_L3D' -x 'L_{3D} of D* Cands' -f 'compare'> /dev/null &


