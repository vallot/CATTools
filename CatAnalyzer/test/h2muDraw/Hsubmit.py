#./DYdraw.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x-name> -y <y_name> -f <f_name> -d <dolog> 
"""
presel_cut  = 'dilep.M()>20&&jet1.Pt()>40&&jet2.Pt()>30&&step==5'
VBFT_cut    = '%s&&dijet.M()>650&&abs(dijet.Eta())>3.5'%(presel_cut)
GGFT_cut    = '%s&&!(%s)&&dijet.M()>250&&dilep.Pt()>50'%(presel_cut, VBFT_cut)
BFL_cut     = '%s&&!(%s)&&!(%s)'%(presel_cut, VBFT_cut, GGFT_cut)
JetT_cut    = '!(%s)&&dilep.M()>20&&step==5&&dilep.Pt()>=25'%(presel_cut)
JetL_cut    = '!(%s)&&dilep.Pt()<25&&dilep.M()>20&&step==5'%(presel_cut)
"""

std_cut     = 'dilep.M()>2&&step==5'

weight      = 'weight*(mueffweight)'

json_used   = 'Golden'

x_name_l    = ["Invariant Mass [GeV]", "Transverse Momentum [GeV]"]

plotvar_l   = ["dilep.M\(\)", "dilep.Pt\(\)"]

binset_l    = ["[300,0,300]", "[250,0,300]", "[200,0,300]"]


lst         = []
for i in [1,2,3,4,5]:
    for b in (plotvar_l):
        for c in (x_name_l):
            for a in (binset_l):

                if (b,c) == ("dilep.M\(\)", "Invariant Mass [GeV]"): 
                    cmd ="./h2muDraw.py -c \'%s&&cat==%s\' -b %s -p %s -x '%s' -w '%s' -j '%s'"%(std_cut, i, a, b, c, weight, json_used)
                    if i == 1:
                        cmd = cmd +' -f \'VBFT_M\''
                    if i == 2:
                        cmd = cmd +' -f \'GGFT_M\'' 
                    if i == 3:
                        cmd = cmd +' -f \'VBFL_M\'' 
                    if i == 4:
                        cmd = cmd +' -f \'JetT_M\''  
                    if i == 5: 
                        cmd = cmd +' -f \'JetL_M\'' 

                elif b == 'dilep.Pt\(\)' and c == "Transverse Momentum [GeV]":
                    cmd ="./h2muDraw.py -c \'%s&&cat==%s\' -b %s -p %s -x '%s' -w '%s' -j '%s'"%(std_cut, i, a, b, c, weight, json_used)
                    if i == 1:
                        cmd = cmd +' -f \'VBFT_Pt\'' 
                    if i == 2:
                        cmd = cmd +' -f \'GGFT_Pt\'' 
                    if i == 3:
                        cmd = cmd +' -f \'VBFL_Pt\'' 
                    if i == 4:
                        cmd = cmd +' -f \'JetT_Pt\'' 
                    if i == 5: 
                        cmd = cmd +' -f \'JetL_Pt\''  
                else:
                    break
                cmd = cmd +' > /dev/null &'
                lst.append(cmd)   

for l in lst:
    print l


