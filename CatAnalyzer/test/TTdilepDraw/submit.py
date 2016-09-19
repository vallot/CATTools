#./topDraw.py.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog>

plotvar_l = ["dilep.M\(\)", "njet", "met", "nbjet",
			 "lep1.Pt\(\),lep2.Pt\(\)", "lep1.Eta\(\),lep2.Eta\(\)",
			 "pseudojet1.Pt\(\),pseudojet2.Pt\(\)", "pseudojet1.Eta\(\),pseudojet2.Eta\(\)",
			 "pseudotop1.Pt\(\),pseudotop2.Pt\(\)", "pseudotop1.Rapidity\(\),pseudotop2.Rapidity\(\)",
			 "pseudottbar.Pt\(\)", "pseudottbar.Rapidity\(\)", "pseudottbar.M\(\)", "abs\(pseudottbar_dphi\)"]
x_name_l  = ["M(ll) [GeV]", "Jet Multiplicity", "Missing Et [GeV]", "b Jet Multiplicity",
             "p_{T}^{l} [GeV]", "lepton #eta",
             "p_{T}^{j} [GeV]", "jet #eta",
             "p_{T}^{t} [GeV]", "y^{t}",
             "p_{T}^{t#bar{t}} [GeV]", "y^{t#bar{t}}", "M^{t#bar{t}} [GeV]", "#Delta#phi^{t#bar{t}}"]
binset_l  = ["[60,20,320]", "[10,0,10]", "[20,0,200]", "[6,0,6]",
			 "[9,20,290]", "[10,-2.5,2.5]",
			 "[9,30,300]", "[10,-2.5,2.5]",
			 "[20,0,500]", "[10,-2.5,2.5]",
			 "[10,0,400]", "[10,-2.5,2.5]", "[10,300,1200]", "[10,0.0,3.141592]"]
"""
binset_l  = ["[60,20,320]", "[10,0,10]", "[20,0,200]", "[6,0,6]",
			 "[20,30,40,60,80,120,180,400]", "[-2.5,-1.5,-1,-0.5,0,0.5,1,1.5,2.5]",
			 "[30,50,80,130,210,500]", "[-2.5,-1.5,-1,-0.5,0,0.5,1,1.5,2.5]",
			 "[0,60,120,200,500]", "[-2.5,-1.5,-0.5,0.5,1.5,2.5]",
			 "[0,20,50,100,200,400]", "[-2.5,-1.5,-1,-0.5,0,0.5,1,1.5,2.5]", "[200,400,520,700,1200]", "[0.0,1.57,2.61,3.016,3.141592]"]
"""

lst=[]
for channel in [0,1,2,3]:
    for step in [1,2,3,4,5,6]:
        for b, p, x in zip(binset_l, plotvar_l, x_name_l):
            cmd = "./topDraw.py -a %d -s %d -b %s -p %s -x '%s' -o"%(channel, step, b, p, x)
            if step <= 5: cmd = cmd+' -d'
            if step >= 5:
                cmd = cmd+' -w \'genweight*puweight*mueffweight*eleffweight*tri*btagweight\''
                cmd = cmd+' -c \'tri!=0&&filtered==1&&is3lep==2&&pseudojet1.Pt()>30&&pseudojet2.Pt()>30&&lep1.Pt()>20&&lep2.Pt()>20\''
            cmd = cmd+' > tmp &'
            lst.append(cmd)

for l in lst:
    print l
