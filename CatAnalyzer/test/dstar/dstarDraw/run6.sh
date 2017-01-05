#!/bin/bash
eval `scramv1 runtime -sh`

#cut='d0_L3D>0.2 && d0_LXY>0.1&& dstar_L3D>0.2 && dstar_LXY>0.1 && abs(dstar_diffMass-0.145)<0.01'
#cut='abs(dstar_diffMass-0.145)<0.02 && dstar_LXY>0.09'
cut_d0='d0_L3D>0.2&&d0_LXY>0.1&&abs(d0.M()-1.8648)<0.050'
cut_dstar='d0_L3D>0.2&&d0_LXY>0.1&&dstar_L3D>0.2&&dstar_LXY>0.1&&abs(dstar_diffMass-0.145)<0.01'
cut_dstar_nomassConstain='d0_L3D>0.2&&d0_LXY>0.1&&dstar_L3D>0.2&&dstar_LXY>0.1'
cut_dstar_noCut='abs(dstar_diffMass-0.145)<0.01'

./massPlot.py -c "${cut_d0}"  -s 5 -b [30,0,120] -p 'd0_lepSV_lowM1' -x 'M_{l+D0}' & 
./massPlot.py -c "${cut_d0}"  -s 5 -b [30,0,120] -p 'd0_lepSV_dRM1' -x 'M_{l+D0}'  &

./massPlot.py -c "${cut_dstar}"  -s 5 -b [30,0,120] -p 'dstar_lepSV_lowM1' -x 'M_{l+D^{*}}'  & 
./massPlot.py -c "${cut_dstar}"  -s 5 -b [30,0,120] -p 'dstar_lepSV_dRM1' -x 'M_{l+D^{*}}'  &
./massPlot.py -c "${cut_dstar}"  -s 5 -b [30,0,120] -p 'dstar_lepSV_correctM' -x 'M_{l+D^{*}}'  &
./massPlot.py -c "${cut_dstar}"  -s 5 -b [30,0,120] -p 'dstar_opCharge_M' -x 'M_{l+D^{*}}'  &

./massPlot.py -c "${cut_dstar_nomassConstain}"  -s 5 -b [10,0,120] -p 'dstar_lepSV_lowM1' -x 'M_{l+D^{*}}' -f 'noMasscut' & 
./massPlot.py -c "${cut_dstar_nomassConstain}"  -s 5 -b [10,0,120] -p 'dstar_lepSV_dRM1' -x 'M_{l+D^{*}}' -f 'noMasscut' &
./massPlot.py -c "${cut_dstar_nomassConstain}"  -s 5 -b [10,0,120] -p 'dstar_lepSV_correctM' -x 'M_{l+D^{*}}' -f 'noMasscut' &
./massPlot.py -c "${cut_dstar_nomassConstain}"  -s 5 -b [10,0,120] -p 'dstar_opCharge_M' -x 'M_{l+D^{*}}' -f 'noMasscut' &

./massPlot.py -c "${cut_dstar_noCut}"  -s 5 -b [10,0,120] -p 'dstar_lepSV_lowM1' -x 'M_{l+D^{*}}' -f 'onlyMassCut' & 
./massPlot.py -c "${cut_dstar_noCut}"  -s 5 -b [10,0,120] -p 'dstar_lepSV_dRM1' -x 'M_{l+D^{*}}' -f 'onlyMassCut' &
./massPlot.py -c "${cut_dstar_noCut}"  -s 5 -b [10,0,120] -p 'dstar_lepSV_correctM' -x 'M_{l+D^{*}}' -f 'onlyMassCut' &
./massPlot.py -c "${cut_dstar_noCut}"  -s 5 -b [10,0,120] -p 'dstar_opCharge_M' -x 'M_{l+D^{*}}' -f 'onlyMassCut' &

