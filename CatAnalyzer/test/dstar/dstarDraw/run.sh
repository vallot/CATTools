#!/bin/bash

treename='nom'

if [ $# != "0" ]; then
	treename=$1
fi

eval `scramv1 runtime -sh`
for step in {1,2,3,4,5} 
do
./dstarDraw.py -t $treename -s $step -b [60,20,320] -p dilep.M\(\) -x 'M(ll) [GeV/c^{2}]' -d -o > /dev/null &
./dstarDraw.py -t $treename -s $step -b [10,0,10] -p njet -x 'Jet Multiplicity' -d -o > /dev/null &
./dstarDraw.py -t $treename -s $step -b [20,0,200] -p met -x 'Missing Et [GeV]' -d -o > /dev/null &
./dstarDraw.py -t $treename -s $step -b [6,0,6] -p nbjet -x 'b Jet Multiplicity' -d -o > /dev/null &
done

#bash ./run2.sh $treename

#bash ./run3.sh $treename
#bash ./run5.sh $treename
#bash ./run6.sh
