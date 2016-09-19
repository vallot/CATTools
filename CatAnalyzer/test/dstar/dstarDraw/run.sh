#!/bin/bash
eval `scramv1 runtime -sh`
for step in {1,5} 
do
./dstarDraw.py -s $step -b [60,20,320] -p dilep.M\(\) -x 'M(ll) [GeV/c^{2}]' -d -o > tmp &
./dstarDraw.py -s $step -b [10,0,10] -p njet -x 'Jet Multiplicity' -d -o > tmp &
./dstarDraw.py -s $step -b [20,0,200] -p met -x 'Missing Et [GeV]' -d -o > tmp &
./dstarDraw.py -s $step -b [6,0,6] -p nbjet -x 'b Jet Multiplicity' -d -o > tmp &
done

bash ./run2.sh

bash ./run3.sh
bash ./run5.sh
#bash ./run6.sh
