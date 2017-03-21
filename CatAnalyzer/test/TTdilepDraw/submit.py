#./topDraw.py.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog>

plotvar_l = ["dilep.M\(\)", "njet", "met", "nbjet",
			 "lep1.Pt\(\),lep2.Pt\(\)", "lep1.Eta\(\),lep2.Eta\(\)",
			 "pseudojet1.Pt\(\),pseudojet2.Pt\(\)", "pseudojet1.Eta\(\),pseudojet2.Eta\(\)",
			 "pseudotop1.Pt\(\),pseudotop2.Pt\(\)", "pseudotop1.Rapidity\(\),pseudotop2.Rapidity\(\)",
			 "pseudottbar.Pt\(\)", "pseudottbar.Rapidity\(\)", "pseudottbar.M\(\)", "abs\(pseudottbar_dphi\)"]
x_name_l  = ["M(ll) [GeV]", "Jet Multiplicity", "Missing Et [GeV]", "b Jet Multiplicity",
             "p_{T}^{lep} [GeV]", "lepton #eta",
             "p_{T}^{jet} [GeV]", "jet #eta",
             "p_{T}^{t} [GeV]", "y^{t}",
             "p_{T}^{t#bar{t}} [GeV]", "y^{t#bar{t}}", "M^{t#bar{t}} [GeV]", "#Delta#phi^{t#bar{t}} [Rad]"]
binset_l  = ["[60,20,320]", "[10,0,10]", "[20,0,200]", "[6,0,6]",
			 "[9,20,290]", "[10,-2.5,2.5]",
			 "[9,30,300]", "[10,-2.5,2.5]",
			 "[20,0,500]", "[10,-2.5,2.5]",
			 "[10,0,400]", "[10,-2.5,2.5]", "[10,300,1200]", "[10,0.0,3.141592]"]

lst=[]
for channel in [0,1,2,3]:
    for step in [1,2,3,4,5,6]:
        for b, p, x in zip(binset_l, plotvar_l, x_name_l):
            cmd = "./topDraw.py -a %d -s %d -b %s -p %s -x '%s' -o"%(channel, step, b, p, x)
            if step <= 5: cmd = cmd+' -d'
            if step >= 5:
                cmd = cmd+' -w \'genweight*puweight*mueffweight*eleffweight*tri*btagweight\''
                cmd = cmd+' -c \'tri!=0&&filtered==1&&is3lep==2&&pseudojet1.Pt()>30&&pseudojet2.Pt()>30&&lep1.Pt()>20&&lep2.Pt()>20\''
            cmd = cmd+' > /dev/null &'
            lst.append(cmd)

for l in lst:
    print l


#print """
"""
./topDraw.py -a 1 -s 1 -b [60,20,320] -p dilep.M\(\) -x 'M(ll) [GeV]' -o -d > /dev/null &
./topDraw.py -a 1 -s 2 -b [60,20,320] -p dilep.M\(\) -x 'M(ll) [GeV]' -o -d > /dev/null &
./topDraw.py -a 1 -s 2 -b [10,0,10] -p njet -x 'Jet Multiplicity' -o -d > /dev/null &
./topDraw.py -a 1 -s 3 -b [10,0,10] -p njet -x 'Jet Multiplicity' -o -d > /dev/null &
./topDraw.py -a 1 -s 3 -b [20,0,200] -p met -x 'Missing Et [GeV]' -o -d > /dev/null &
./topDraw.py -a 1 -s 4 -b [20,0,200] -p met -x 'Missing Et [GeV]' -o -d > /dev/null &
./topDraw.py -a 1 -s 4 -b [6,0,6] -p nbjet -x 'b Jet Multiplicity' -o -d > /dev/null &
./topDraw.py -a 1 -s 5 -b [6,0,6] -p nbjet -x 'b Jet Multiplicity' -o -d -w 'genweight*puweight*mueffweight*eleffweight*tri*btagweight' -c 'tri!=0&&filtered==1&&is3lep==2&&pseudojet1.Pt()>30&&pseudojet2.Pt()>30&&lep1.Pt()>20&&lep2.Pt()>20' > /dev/null &
./topDraw.py -a 1 -s 6 -b [9,20,290] -p lep1.Pt\(\),lep2.Pt\(\) -x 'p_{T}^{lep} [GeV]' -o -w 'genweight*puweight*mueffweight*eleffweight*tri*btagweight' -c 'tri!=0&&filtered==1&&is3lep==2&&pseudojet1.Pt()>30&&pseudojet2.Pt()>30&&lep1.Pt()>20&&lep2.Pt()>20' > /dev/null &
./topDraw.py -a 1 -s 6 -b [10,-2.5,2.5] -p lep1.Eta\(\),lep2.Eta\(\) -x 'lepton #eta' -o -w 'genweight*puweight*mueffweight*eleffweight*tri*btagweight' -c 'tri!=0&&filtered==1&&is3lep==2&&pseudojet1.Pt()>30&&pseudojet2.Pt()>30&&lep1.Pt()>20&&lep2.Pt()>20' > /dev/null &
./topDraw.py -a 1 -s 6 -b [9,30,300] -p pseudojet1.Pt\(\),pseudojet2.Pt\(\) -x 'p_{T}^{jet} [GeV]' -o -w 'genweight*puweight*mueffweight*eleffweight*tri*btagweight' -c 'tri!=0&&filtered==1&&is3lep==2&&pseudojet1.Pt()>30&&pseudojet2.Pt()>30&&lep1.Pt()>20&&lep2.Pt()>20' > /dev/null &
./topDraw.py -a 1 -s 6 -b [10,-2.5,2.5] -p pseudojet1.Eta\(\),pseudojet2.Eta\(\) -x 'jet #eta' -o -w 'genweight*puweight*mueffweight*eleffweight*tri*btagweight' -c 'tri!=0&&filtered==1&&is3lep==2&&pseudojet1.Pt()>30&&pseudojet2.Pt()>30&&lep1.Pt()>20&&lep2.Pt()>20' > /dev/null &
./topDraw.py -a 1 -s 6 -b [20,0,500] -p pseudotop1.Pt\(\),pseudotop2.Pt\(\) -x 'p_{T}^{t} [GeV]' -o -w 'genweight*puweight*mueffweight*eleffweight*tri*btagweight' -c 'tri!=0&&filtered==1&&is3lep==2&&pseudojet1.Pt()>30&&pseudojet2.Pt()>30&&lep1.Pt()>20&&lep2.Pt()>20' > /dev/null &
./topDraw.py -a 1 -s 6 -b [10,-2.5,2.5] -p pseudotop1.Rapidity\(\),pseudotop2.Rapidity\(\) -x 'y^{t}' -o -w 'genweight*puweight*mueffweight*eleffweight*tri*btagweight' -c 'tri!=0&&filtered==1&&is3lep==2&&pseudojet1.Pt()>30&&pseudojet2.Pt()>30&&lep1.Pt()>20&&lep2.Pt()>20' > /dev/null &
./topDraw.py -a 1 -s 6 -b [10,0,400] -p pseudottbar.Pt\(\) -x 'p_{T}^{t#bar{t}} [GeV]' -o -w 'genweight*puweight*mueffweight*eleffweight*tri*btagweight' -c 'tri!=0&&filtered==1&&is3lep==2&&pseudojet1.Pt()>30&&pseudojet2.Pt()>30&&lep1.Pt()>20&&lep2.Pt()>20' > /dev/null &
./topDraw.py -a 1 -s 6 -b [10,-2.5,2.5] -p pseudottbar.Rapidity\(\) -x 'y^{t#bar{t}}' -o -w 'genweight*puweight*mueffweight*eleffweight*tri*btagweight' -c 'tri!=0&&filtered==1&&is3lep==2&&pseudojet1.Pt()>30&&pseudojet2.Pt()>30&&lep1.Pt()>20&&lep2.Pt()>20' > /dev/null &
./topDraw.py -a 1 -s 6 -b [10,300,1200] -p pseudottbar.M\(\) -x 'M^{t#bar{t}} [GeV]' -o -w 'genweight*puweight*mueffweight*eleffweight*tri*btagweight' -c 'tri!=0&&filtered==1&&is3lep==2&&pseudojet1.Pt()>30&&pseudojet2.Pt()>30&&lep1.Pt()>20&&lep2.Pt()>20' > /dev/null &
./topDraw.py -a 1 -s 6 -b [10,0.0,3.141592] -p abs\(pseudottbar_dphi\) -x '#Delta#phi^{t#bar{t}} [Rad]' -o -w 'genweight*puweight*mueffweight*eleffweight*tri*btagweight' -c 'tri!=0&&filtered==1&&is3lep==2&&pseudojet1.Pt()>30&&pseudojet2.Pt()>30&&lep1.Pt()>20&&lep2.Pt()>20' > /dev/null &

"""
