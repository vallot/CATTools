#./topDraw.py.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog>

plotvar_l = ["ll_m", "njet", "met", "nbjet",
			 "lep1_pt,lep2_pt", "lep1_eta,lep2_eta",
			 "jet1_pt,jet2_pt", "jet1_eta,jet2_eta",
			 "top1_pt,top2_pt", "top1_rapi,top2_rapi",
			 "ttbar_pt", "ttbar_rapi", "ttbar_m", 'ttbar_dphi']
x_name_l  = ["M(ll) [GeV/c^{2}]", "Jet Multiplicity", "Missing Et [GeV]", "b Jet Multiplicity",
			 "lepton p_{T} [GeV/c]", "lepton #eta",
			 "Jet p_{T} [GeV/c]", "Jet #eta",
			 "Top p_{T} [GeV/c]", "Top Rapidity",
			 "TTbar p_{T} [GeV/c]", "TTbar Rapidity", "TTbar Mass [GeV/c^{2}]", "TTbar #Delta#phi"]
binset_l  = ["[60,20,320]", "[10,0,10]", "[20,0,200]", "[6,0,6]",
			 "[10,0,300]", "[10,-2.5,2.5]",
			 "[10,0,300]", "[10,-2.5,2.5]",
			 "[20,0,500]", "[10,-2.5,2.5]",
			 "[10,0,400]", "[10,-2.5,2.5]", "[10,200,1200]", "[10,0.0,3.141592]"]
			 #"[20,30,40,60,80,120,180,400]", "[-2.5,-1.5,-1,-0.5,0,0.5,1,1.5,2.5]",
			 #"[30,50,80,130,210,500]", "[-2.5,-1.5,-1,-0.5,0,0.5,1,1.5,2.5]",
			 #"[0,60,120,200,500]", "[-2.5,-1.5,-0.5,0.5,1.5,2.5]",
			 #"[0,20,50,100,200,400]", "[-2.5,-1.5,-1,-0.5,0,0.5,1,1.5,2.5]", "[200,400,520,700,1200]", "[0.0,1.57,2.61,3.016,3.141592]"]

#os.system(command)
for channel in range(1,4):
	for step in range(1,6):
		for p, plotvar in enumerate(plotvar_l[:4]):
			if step < 5:
				print "./topDraw.py -a %d -s %d -b %s -p %s -x '%s' -d -o > tmp &"%(channel, step, binset_l[p], plotvar, x_name_l[p])
			else:
				print "./topDraw.py -a %d -s %d -b %s -p %s -x '%s' -d -o -w 'genweight*puweight*mueffweight*eleffweight*tri*btagweight' > tmp &"%(channel, step, binset_l[p], plotvar, x_name_l[p])
	print

for i in range(1,4):
	for j in range(6,7):
		for p, plotvar in enumerate(plotvar_l[4:]):
			if step < 5:
				print "./topDraw.py -a %d -s %d -b %s -p %s -x '%s' > tmp &"%(i, j, binset_l[4+p], plotvar, x_name_l[4+p])
			else:
				print "./topDraw.py -a %d -s %d -b %s -p %s -x '%s' -w 'genweight*puweight*mueffweight*eleffweight*tri*btagweight' > tmp &"%(i, j, binset_l[4+p], plotvar, x_name_l[4+p])
	print

for j in range(1,6):
	for p, plotvar in enumerate(plotvar_l[:4]):
		print "./topDraw_combined.py -s %d -b %s -p %s -x '%s' -d -o > tmp &"%(j, binset_l[p], plotvar, x_name_l[p])
print

for j in range(6,7):
	for p, plotvar in enumerate(plotvar_l[4:]):
		print "./topDraw_combined.py -s %d -b %s -p %s -x '%s' > tmp &"%(j, binset_l[4+p], plotvar, x_name_l[4+p])
print

for i in range(1,4):
	for p, plotvar in enumerate(plotvar_l[4:]):
		print "./topDraw_top.py -a %d -p %s -b %s > tmp &"%(i, plotvar, binset_l[4+p])
print

for p, plotvar in enumerate(plotvar_l[4:]):
	print "./topDraw_top_combined.py -p %s -b %s > tmp &"%(plotvar, binset_l[4+p])
print

for i in range(1,4):
	for p, plotvar in enumerate(plotvar_l[4:]):
		print "./topDraw_sys_bkg.py -a %d -p %s -b %s > tmp &"%(i, plotvar, binset_l[4+p])
print

for p, plotvar in enumerate(plotvar_l[4:]):
	print "./topDraw_sys_combined.py -p %s -b %s > tmp &"%(plotvar, binset_l[4+p])
print

