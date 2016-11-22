#!/bin/bash

treename='nom'

if [ $# != "0" ]; then
	treename=$1
fi

eval `scramv1 runtime -sh`

cut='(abs(d0_dRTrue)<0.1&&abs(d0_relPtTrue)<0.1)'
cut2='(abs(d0_dRTrue)<0.1&&abs(d0_relPtTrue)<0.1)&&abs(dstar_dRTrue)<0.1&&abs(dstar_dRTrue)<0.1&&abs(dstar_isFromTop)==6'

relcut='(d0_LXY>0.1&&d0_L3D>0.2)'
relcut2='(d0_LXY>0.1&&d0_L3D>0.2&&dstar_LXY>0.1&&dstar_L3D>0.2)'

./dstarDraw.py -t $treename -c "${relcut}" -s 5 -b [20,1.6,2.2] -p d0.M\(\) -x 'Mass of D0 Cands [GeV/c^2]' -f 'cut_applied' > /dev/null &
./dstarDraw.py -t $treename -c "${relcut2}" -s 5 -b [20,1.9,2.1] -p dstar.M\(\) -x 'Mass of Dstar Cands [GeV/c^2]' -f 'cut_applied' > /dev/null &
./dstarDraw.py -t $treename -c "${relcut2}" -s 5 -b [20,0.135,0.17] -p dstar_diffMass -x 'M_{K #pi #pi}-M_{K #pi} [GeV/c^2]' -f 'cut_applied'> /dev/null &


