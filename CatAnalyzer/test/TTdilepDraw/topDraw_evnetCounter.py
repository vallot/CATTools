#!/usr/bin/env python
from CATTools.CatAnalyzer.histoHelper_TtbarDilep import *

# step infor: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TTbarXSecSynchronization
def printbyCut(cut):
	print "%8d"%(getWeightedEntries(rootfile, tname, 'tri', cut)),
	

rootfile = '/cms/scratch/tt8888tt/cattools_v746/src/CATTools/CatAnalyzer/test/cattree.root'
#rootfile = '/cms/scratch/tt8888tt/cattools_v746/src/CATTools/CatAnalyzer/test/v7-4-6/TT_powheg.root'

tname = "cattree/nom"
fstlow = ['    ', 'elmu', 'elel', 'mumu']

print fstlow[0],
for s in ['0a','0b','0c']+range(1, 7):
	print "%8s"%("step"+str(s)),
print

for channel in range(1, 4):
	print fstlow[channel],

	step0_cut = 'step>=%i'%(-1)
	printbyCut(step0_cut)
	printbyCut(step0_cut+'&&tri==1')
	printbyCut(step0_cut+'&&tri==1&&filtered==1')

	for step in range(1, 7):
		if channel == 1 : partonCh_cut = '&&((parton_mode1 == 1 && parton_mode2 == 2) || (parton_mode1 == 2 && parton_mode2 == 1))'
		elif channel == 2 : partonCh_cut = '&&(parton_mode1 == 2 && parton_mode2 == 2)'
		elif channel == 3 : partonCh_cut = '&&(parton_mode1 == 1 && parton_mode2 == 1)'

		stepch_cut = 'step>=%i&&channel==%i'%(step,channel)
		#stepch_cut = stepch_cut + partonCh_cut	
		printbyCut(stepch_cut+'&&tri==1&&filtered==1')
	print
print

