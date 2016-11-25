#!/bin/bash

eval `scramv1 runtime -sh`

rm invmass/res_[0-9]*.txt

type='null'

if [ $# != "0" ]; then
	type=$1
fi

#cut='d0_L3D>0.2 && d0_LXY>0.1&& dstar_L3D>0.2 && dstar_LXY>0.1 && abs(dstar_diffMass-0.145)<0.01'
#cut='abs(dstar_diffMass-0.145)<0.02 && dstar_LXY>0.09'
cut_dstar='d0_L3D>0.2&&d0_LXY>0.1&&dstar_L3D>0.2&&dstar_LXY>0.1&&abs(dstar_diffMass-0.145)<0.01'
cut_dstar_nomassConstain='d0_L3D>0.2&&d0_LXY>0.1&&dstar_L3D>0.2&&dstar_LXY>0.1'
cut_dstar_noCut='abs(dstar_diffMass-0.145)<0.01'

./massPlot.py -c "${cut_dstar}"  -s 5 -b [30,0,150] -p 'dstar_lepSV_lowM' -x 'M_{l+D^{*}} [GeV/c^2]' -t $type >> invmass/res_101.txt  & 
./massPlot.py -c "${cut_dstar}"  -s 5 -b [30,0,150] -p 'dstar_lepSV_dRM' -x 'M_{l+D^{*}} [GeV/c^2]' -t $type >> invmass/res_102.txt &
./massPlot.py -c "${cut_dstar}"  -s 5 -b [30,0,150] -p 'dstar_opCharge_M' -x 'M_{l+D^{*}} [GeV/c^2]' -t $type >> invmass/res_103.txt &

./massPlot.py -c "${cut_dstar_nomassConstain}"  -s 5 -b [10,0,150] -p 'dstar_lepSV_lowM' -x 'M_{l+D^{*}} [GeV/c^2]' -f 'noMasscut' -t $type >> invmass/res_111.txt & 
./massPlot.py -c "${cut_dstar_nomassConstain}"  -s 5 -b [10,0,150] -p 'dstar_lepSV_dRM' -x 'M_{l+D^{*}} [GeV/c^2]' -f 'noMasscut' -t $type >> invmass/res_112.txt &
./massPlot.py -c "${cut_dstar_nomassConstain}"  -s 5 -b [10,0,150] -p 'dstar_opCharge_M' -x 'M_{l+D^{*}} [GeV/c^2]' -f 'noMasscut' -t $type >> invmass/res_113.txt &

./massPlot.py -c "${cut_dstar_noCut}"  -s 5 -b [10,0,150] -p 'dstar_lepSV_lowM' -x 'M_{l+D^{*}} [GeV/c^2]' -f 'onlyMassCut' -t $type >> invmass/res_121.txt & 
./massPlot.py -c "${cut_dstar_noCut}"  -s 5 -b [10,0,150] -p 'dstar_lepSV_dRM' -x 'M_{l+D^{*}} [GeV/c^2]' -f 'onlyMassCut' -t $type >> invmass/res_122.txt &
./massPlot.py -c "${cut_dstar_noCut}"  -s 5 -b [10,0,150] -p 'dstar_opCharge_M' -x 'M_{l+D^{*}} [GeV/c^2]' -f 'onlyMassCut' -t $type >> invmass/res_123.txt &

#./massPlot.py -s 0 -b [30,0,150] -p 'dstar_lepSV_correctM' -x 'M_{l+D^{*}} [GeV/c^2]' -t $type >> invmass/res_130.txt &


