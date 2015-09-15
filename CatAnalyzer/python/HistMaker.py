#!/usr/bin/env python

from ROOT import *
from array import *

class HistMaker(object):
    def __init__(self, inFileNames, outFileName,
                 modName = "ntuple", treeName = "event", eventCount = "hNEvent"):
        self.weightVar = "1"
        self.precut = "1"
        self.h1 = {}
        self.h2 = {}
        self.cutSteps = []
        self.ntupleStep = None
        self.isMC = False
        self.isOK = False

        try:
            self.nEventTotal = float(eventCount)
            eventCounterName = None
        except:
            eventCounterName = eventCount
            self.nEventTotal = 0

        ## Load input tree and event number for scalers
        self.chain = TChain("%s/%s" % (modName, treeName))
        if hasattr(inFileNames, 'split'): inFileNames = inFileNames.split(',')
        for fName in inFileNames:
            f = TFile.Open(fName)
            if None == f: continue
            self.chain.Add(fName)
            self.isOK = True

            if eventCounterName is not None:
                hNEvent = f.Get("%s/%s" % (modName, eventCounterName))
                if hNEvent != None: self.nEventTotal += hNEvent.GetBinContent(1)
 
        self.scale = max(1, self.nEventTotal)
        self.outFile = TFile(outFileName, "RECREATE")
        if self.outFile == None or not self.outFile.IsOpen(): self.isOK = False

    def setWeightVar(self, weightVar):
        self.weightVar = weightVar
        if weightVar != "1": self.isMC = True

    def setPrecut(self, precut):
        self.precut = precut

    def addCutStep(self, name, cut, histNames):
        if cut == '': cut = "1"
        if type(histNames) == type(''): histNames = [x.strip() for x in histNames.split(',')]
        self.cutSteps.append( (name, cut, histNames) )

    def addH1(self, *args):
        name, varexp, title = args[:3]
        binOptions = args[3:]
        if len(binOptions) == 1:
            bins = array('d', binOptions[0])
        elif len(binOptions) == 3:
            nbins, min, max = binOptions
            bins = array('d', [min+(max-min)*i/nbins for i in range(nbins+1)])
        else:
            print "!!! Wrong input to addH1()"
            return
        self.h1[name] = (varexp, title, bins)

    def addH2(self, *args):
        name, varexp, title = args[:3]
        binOptions = args[3:]
        if len(binOptions) == 2:
            binsX = array('d', binOptions[0])
            binsY = array('d', binOptions[1])
        elif len(binOptions) == 6:
            nbinsX, minX, maxX = binOptions[:3]
            nbinsY, minY, maxY = binOptions[3:]
            binsX = array('d', [minX+(maxX-minX)*i/nbinsX for i in range(nbinsX+1)])
            binsY = array('d', [minY+(maxY-minY)*i/nbinsY for i in range(nbinsY+1)])
        else:
            print "!!! Wrong input to addH2()"
            return
        self.h2[name] = (varexp, title, binsX, binsY)

    def storeNtuple(self, step):
        self.ntupleStep = step
       
    def process(self):
        if not self.isOK:
            print "HistMaker is not initialized, ", self.outFile
            return

        self.outFile.cd()
        nCutStep = len(self.cutSteps)
        hNEvent = TH1F("hNEvent", "hNEvent", nCutStep+2, 1, nCutStep+3)
        hWeight = TH1F("hWeight", "hWeight", nCutStep+2, 1, nCutStep+3)

        hNEvent.GetXaxis().SetBinLabel(1, "Total")
        hWeight.GetXaxis().SetBinLabel(1, "Total")
        hNEvent.SetBinContent(1, self.nEventTotal)
        hWeight.SetBinContent(1, self.nEventTotal)

        stackedCut = self.precut[:]
        nPassed = self.chain.Draw(">>eventList", stackedCut)
        hNEvent.GetXaxis().SetBinLabel(2, "Precut")
        hNEvent.SetBinContent(2, nPassed)
        hWeight.GetXaxis().SetBinLabel(2, "Precut")
        if self.isMC: self.chain.Draw("2>>+hWeight", "(%s)*(%s)" % (self.weightVar, stackedCut), "goff")
        else: self.chain.Draw("2>>+hWeight", "%s" % stackedCut, "goff")

        for i, cutStepInfo in enumerate(self.cutSteps):
            name, cut, histNames = cutStepInfo

            self.outFile.cd()
            stackedCut = "(%s) && (%s)" % (stackedCut, cut)
            nPassed = self.chain.Draw(">>eventList", stackedCut)
            hNEvent.GetXaxis().SetBinLabel(i+3, name)
            hNEvent.SetBinContent(i+3, nPassed)
            hWeight.GetXaxis().SetBinLabel(i+3, name)
            if self.isMC: self.chain.Draw("%d>>+hWeight" % (i+3), "(%s)*(%s)" % (self.weightVar, stackedCut), "goff")
            else: self.chain.Draw("%d>>+hWeight" % (i+3), "%s" % stackedCut, "goff")

            print "Cut step %d/%d (%s) : %.1f %.1f" % (i+1, len(self.cutSteps), name, hNEvent.GetBinContent(i+3), hWeight.GetBinContent(i+3))

            cutStepDir = self.outFile.mkdir(name)
            cutStepDir.cd()

            for histName in histNames:
                if histName in self.h1:
                    varexp, title, bins = self.h1[histName]
                    h = TH1F(histName, title, len(bins)-1, bins)
                    h.Sumw2()
                    if self.isMC: h.SetOption("hist")
                    xmax = h.GetXaxis().GetXmax()-1e-9*h.GetXaxis().GetBinWidth(len(bins))
                    if self.isMC: self.chain.Draw("min(%s,%f)>>%s" % (varexp, xmax, histName), "(%s)*(%s)" % (self.weightVar, stackedCut), "goff")
                    else: self.chain.Draw("min(%s,%f)>>%s" % (varexp, xmax, histName), "%s" % stackedCut, "goff")
                elif histName in self.h2:
                    varexp, title, binsX, binsY = self.h2[histName]
                    h = TH2F(histName, title, len(binsX)-1, binsX, len(binsY)-1, binsY)
                    h.Sumw2()
                    if self.isMC: self.chain.Draw("%s>>%s" % (varexp, histName), "(%s)*(%s)" % (self.weightVar, stackedCut), "goff")
                    else: self.chain.Draw("%s>>%s" % (varexp, histName), "%s" % stackedCut, "goff")
                else:
                    print "Histogram", histName, "in cut step", name, "not defined."
                    continue

                h.Scale(1./self.scale)
                h.Write()
                h = None

            if self.ntupleStep != None and name == self.ntupleStep:
                print "Store ntuple"
                ntuple = self.chain.CopyTree(stackedCut)
                ntuple.Write()
                print "Stored ntuple"

        self.outFile.cd()
        hNEvent.Write()
        hWeight.Write()
        self.outFile.Close()

