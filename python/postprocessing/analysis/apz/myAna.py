#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducer import *

import argparse
parser = argparse.ArgumentParser("")
parser.add_argument('--local', dest='local', action="store_true", default=False, help="Run local")
parser.add_argument('--jobNum', dest='jobNum', type=int, default=1, help="")
opt = parser.parse_args()


class MyAnalysis(Module):
    def __init__(self):
	self.writeHistFile=True

    def beginJob(self,histFile=None,histDirName=None):
	Module.beginJob(self,histFile,histDirName)

	#self.h_vpt=ROOT.TH1F('sumpt',   'sumpt',   100, 0, 1000)
        #self.addObject(self.h_vpt)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("sel1",  "O");
        self.out.branch("sel2",  "O");
        self.out.branch("sel3",  "O");
        self.out.branch("nMu",  "I");
        self.out.branch("mll",  "F");
        self.out.branch("mllg",  "F");
        self.out.branch("dR_mu1g",  "F");
        self.out.branch("dR_mu2g",  "F");

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        muons = Collection(event, "Muon")
        photons = Collection(event, "Photon")
        #electrons = Collection(event, "Electron")
        #jets = Collection(event, "Jet")

        nMuons = len(muons)
        sel1, sel2, sel3 = 0,0,0
        myMuons  = [x for x in muons if x.pt > 6  and x.softId == 1]
        myPhotons = [x for x in photons if x.pt > 30  and x.cutBased17Bitmap >= 2 and abs(x.eta) < 1.6 ]
        mll, mllg, dR_mu1g, dR_mu2g = -11, -11, -11, -11

        if len(myMuons)>=2 and len(myPhotons)>=1:
            lep1, lep2, gamma = ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector()
            lep1 = myMuons[0].p4()
            lep2 = myMuons[1].p4()
            gamma = myPhotons[0].p4()
            mll = (lep1 + lep2).M()
            mllg = (lep1 + lep2 + gamma).M()
            if abs(myPhotons[0].eta) < 1.4442 and myMuons[0].pt > 23 and myMuons[0].pfRelIso04_all < 0.4:
                sel1 = 1
            if myMuons[0].pt > 24 and myMuons[0].pfRelIso04_all < 0.4:
                sel2 = 1
            if myPhotons[0].pt > 31 and abs(myPhotons[0].eta) < 1.5 and myMuons[0].pt > 24 and myMuons[0].pfRelIso04_all < 0.4:
                sel3 = 1

        self.out.fillBranch("nMu",nMuons)
        self.out.fillBranch("mll",mll)
        self.out.fillBranch("mllg",mllg)
        self.out.fillBranch("dR_mu1g",dR_mu1g)
        self.out.fillBranch("dR_mu2g",dR_mu2g)
        self.out.fillBranch("sel1",sel1)
        self.out.fillBranch("sel2",sel2)
        self.out.fillBranch("sel3",sel3)

        return True


#preselection="Sum$(Muon_pt > 5.5 && Muon_softId) >=2"
#preselection="Sum$(Muon_pt > 22 && Muon_softId && Muon_pfRelIso04_all < 0.4) >=1"
#preselection="(Sum$(Photon_pt > 28 && abs(Photon_eta) < 1.6 && Photon_cutBasedBitmap >= 2) >= 1)"
preselection='''Sum$(Muon_pt > 5.5 && Muon_softId) >= 2 &&
 Sum$(Muon_pt > 22  && Muon_softId && Muon_pfRelIso04_all < 0.4) >=1 &&
 Sum$(Photon_pt > 28 && abs(Photon_eta) < 1.6) >= 1'''

files=["/eos/cms/store/user/andrey/f.root"]
files=["/eos/cms/store/user/andrey/f_Nano14Dec2018_data.root"]
#files=["root://cms-xrd-global.cern.ch://store/mc/RunIISummer16NanoAODv4/DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/50000/F78ECFA1-42BF-9A49-8462-A8B868EFF690.root"]
#this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles,runsAndLumis


if opt.local:
    p=PostProcessor("./outDir", files, cut=preselection.replace('\n',' '), branchsel='keep_and_drop_input.txt', outputbranchsel='keep_and_drop_output.txt',
                    modules=[MyAnalysis()], histFileName="./outDir/histOut.root",histDirName="plots")
else:
    p=PostProcessor(".", inputFiles(), cut=preselection.replace('\n',' '), branchsel='keep_and_drop_input.txt', outputbranchsel='keep_and_drop_output.txt',
                    modules=[muonScaleRes2016(),MyAnalysis()], provenance=True,fwkJobReport=True,jsonInput=runsAndLumis())

p.run()

print "ls -lR:"
os.system("ls -lR")
print "DONE"
