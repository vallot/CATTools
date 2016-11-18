#!/bin/bash

treename='nom'

if [ $# != "0" ]; then
	treename=$1
fi

eval `scramv1 runtime -sh`
for step in 1 5 
do
./dstarDraw.py -t $treename -s $step -b [20,1.6,2.2] -p d0.M\(\) -x 'Mass of D0 Cands [GeV/c^{2}]' > /dev/null &
./dstarDraw.py -t $treename -s $step -b [20,1.9,2.1] -p dstar.M\(\) -x 'Mass of Dstar Cands [GeV/c^{2}]' > /dev/null &
./dstarDraw.py -t $treename -s $step -b [20,0.135,0.17] -p 'dstar_diffMass' -x 'M_{K #pi #pi}-M_{K #pi} [GeV/c^{2}]' > /dev/null &
./dstarDraw.py -t $treename -s $step -b [60,0,0.1] -p 'd0_dca' -x 'DCA of D0 Cands [cm]' > /dev/null &
./dstarDraw.py -t $treename -s $step -b [60,0,0.1] -p 'dstar_dca' -x 'DCA of D* Cands [cm]' > /dev/null &
./dstarDraw.py -t $treename -s $step -b [60,0,0.3] -p 'd0_vProb' -x 'Vertex Probability of D0 Cands' > /dev/null &
./dstarDraw.py -t $treename -s $step -b [60,0,0.3] -p 'dstar_vProb' -x 'Vertex Probability of D* Cands' > /dev/null &
./dstarDraw.py -t $treename -s $step -b [60,0,0.5] -p 'd0_LXY' -x 'Lxy of D0 Cands [cm]' > /dev/null &
./dstarDraw.py -t $treename -s $step -b [60,0,0.5] -p 'dstar_LXY' -x 'Lxy of D* Cands [cm]' > /dev/null &
./dstarDraw.py -t $treename -s $step -b [60,0,1] -p 'd0_L3D' -x 'L_{3D} of D0 Cands [cm]' > /dev/null &
./dstarDraw.py -t $treename -s $step -b [60,0,1] -p 'dstar_L3D' -x 'L_{3D} of D* Cands [cm]' > /dev/null &
done 

