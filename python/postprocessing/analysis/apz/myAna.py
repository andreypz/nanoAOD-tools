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
parser.add_argument('--jobNum', dest='jobNum', type=int, default=0, help="Job number (defines the input file)")
parser.add_argument('--evt', dest='evt', type=int, default=None, help="Maximum events to process (in local mode only)")
parser.add_argument('--era', dest='era', type=int, default=None, choices = [2016, 2017, 2018], help="Era")
parser.add_argument('--muCorr', dest='muCorr', action="store_true", default=False, help="Use muon pt correction")

opt = parser.parse_args()


class MyAnalysis(Module):
    def __init__(self, era, muCorr):
	self.writeHistFile=True
        self.era = era
        self.muCorr = muCorr

    def beginJob(self,histFile=None,histDirName=None):
	Module.beginJob(self,histFile,histDirName)

	#self.h_vpt=ROOT.TH1F('sumpt',   'sumpt',   100, 0, 1000)
        #self.addObject(self.h_vpt)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("a_sel1",  "O");
        self.out.branch("a_sel2",  "O");
        self.out.branch("a_sel3",  "O");

        self.out.branch("a_sel1_t",  "O");
        self.out.branch("a_sel2_t",  "O");
        self.out.branch("a_sel3_t",  "O");
        self.out.branch("a_sel4_t",  "O");

        self.out.branch("a_sel1_j",  "O");
        self.out.branch("a_sel2_j",  "O");
        self.out.branch("a_sel3_j",  "O");

        self.out.branch("a_trg_1Mu",  "O");
        self.out.branch("a_trg_2Mu",  "O");
        self.out.branch("a_trg_MuPho","O");
        self.out.branch("a_mll",  "F");

        self.out.branch("a_nMu",   "I");
        self.out.branch("a_nPho",  "I");
        self.out.branch("a_nPho_t","I");
        self.out.branch("a_nJet",  "I");

        self.out.branch("a_mllg",  "F");
        self.out.branch("a_Pt_mu1",  "F");
        self.out.branch("a_dR_mu1g",  "F");
        self.out.branch("a_dR_mu2g",  "F");

        self.out.branch("a_mllj",  "F");
        self.out.branch("a_dR_mu1j",  "F");
        self.out.branch("a_dR_mu2j",  "F");

        self.out.branch("a_mllg_t",  "F");
        self.out.branch("a_dR_mu1g_t",  "F");
        self.out.branch("a_dR_mu2g_t",  "F");

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        muons = Collection(event, "Muon")
        photons = Collection(event, "Photon")
        #electrons = Collection(event, "Electron")
        jets = Collection(event, "Jet")

        if self.era == 2016:
            trg_MuPho = event.HLT_Mu17_Photon30_CaloIdL_L1ISO

            trg_2Mu = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL | event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ | \
                      event.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL | event.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ

            trg_1Mu = event.HLT_IsoMu24 | event.HLT_IsoTkMu24

        if self.era == 2017:
            try:
                trg_MuPho = event.HLT_Mu17_Photon30_IsoCaloId
                # it only exists starting from 2017C-era
            except:
                trg_MuPho = 0
            trg_2Mu = event.HLT_DoubleMu20_7_Mass0to30_Photon23
            trg_1Mu = event.HLT_IsoMu27

        if self.era == 2018:
            trg_MuPho = event.HLT_Mu17_Photon30_IsoCaloId
            trg_2Mu = event.HLT_DoubleMu20_7_Mass0to30_Photon23
            trg_1Mu = event.HLT_IsoMu24

        sel1, sel2, sel3 = 0,0,0
        sel1_t, sel2_t, sel3_t, sel4_t = 0,0,0,0
        sel1_j, sel2_j, sel3_j = 0,0,0
        myMuons  = [x for x in muons if x.pt > 6  and x.softId == 1]
        nMuons = len(myMuons)

        if nMuons<2:
            return False

        lep1 = myMuons[0].p4()
        lep2 = myMuons[1].p4()
        mll = (lep1 + lep2).M()

        mu1Pt = lep1.Pt()

        if self.era==2016:
            #myPhotons = [x for x in photons if x.pt > 30  and x.cutBased >= 2 and abs(x.eta) < 1.6]
            myPhotons = [x for x in photons if x.pt > 30  and x.cutBased >= 2 and abs(x.eta) < 1.6 and lep1.DeltaR(x.p4()) > 0.3 and lep2.DeltaR(x.p4()) > 0.3]
        else:
            myPhotons = [x for x in photons if x.pt > 30  and x.cutBasedBitmap >= 2 and abs(x.eta) < 1.6 and lep1.DeltaR(x.p4()) > 0.3 and lep2.DeltaR(x.p4()) > 0.3]

        myPhotons_tight_MVA = [x for x in photons if x.pt > 30  and x.mvaID_WP90==True and abs(x.eta) < 1.6 and lep1.DeltaR(x.p4()) > 0.3 and lep2.DeltaR(x.p4()) > 0.3]
        myJets = [x for x in jets if x.pt > 30 and x.jetId >= 2 and abs(x.eta) < 1.6 and lep1.DeltaR(x.p4()) > 0.3 and lep2.DeltaR(x.p4()) > 0.3]

        nPho   = len(myPhotons)
        nPho_t = len(myPhotons_tight_MVA)
        nJet   = len(myJets)
                
        mllg, dR_mu1g, dR_mu2g = -11, -11, -11
        mllg_t, dR_mu1g_t, dR_mu2g_t = -11, -11, -11
        mllj, dR_mu1j, dR_mu2j = -11, -11, -11

        # print(nPho, nPho_t, nJet)
        if nPho==0 and nPho_t==0 and nJet==0:
            # print("No photons or jets. So we pass to next event")
            return False

        if nPho >= 1:
            #if lep1.Pt() < lep2.Pt() :
            #    print ("Warning: worng pt-ordering of muons!!")
            
            gamma = myPhotons[0].p4()
            mllg = (lep1 + lep2 + gamma).M()
            dR_mu1g = lep1.DeltaR(gamma)
            dR_mu2g = lep2.DeltaR(gamma)
            if abs(gamma.Eta()) < 1.4442 and myMuons[0].pt > 23 and myMuons[0].pfRelIso04_all < 0.4:
                sel1 = 1
            if myMuons[0].pt > 24 and myMuons[0].pfRelIso04_all < 0.4:
                sel2 = 1
            if gamma.Pt() > 31 and abs(gamma.Eta()) < 1.5 and myMuons[0].pt > 24 and myMuons[0].pfRelIso04_all < 0.4:
                sel3 = 1

        if nPho_t >= 1:
            gamma_t = myPhotons_tight_MVA[0].p4()
            mllg_t = (lep1 + lep2 + gamma_t).M()
            dR_mu1g_t = lep1.DeltaR(gamma_t)
            dR_mu2g_t = lep2.DeltaR(gamma_t)
            if abs(gamma_t.Eta()) < 1.4442 and myMuons[0].pt > 23 and myMuons[0].pfRelIso04_all < 0.4:
                sel1_t = 1
            if myMuons[0].pt > 24 and myMuons[0].pfRelIso04_all < 0.4:
                sel2_t = 1
            if gamma_t.Pt() > 31 and abs(gamma_t.Eta()) < 1.5 and myMuons[0].pt > 24 and myMuons[0].pfRelIso04_all < 0.4:
                sel3_t = 1
            if gamma_t.Pt() > 31 and abs(gamma_t.Eta()) < 1.5 and myMuons[0].pt > 26 and myMuons[0].pfRelIso04_all < 0.1:
                sel4_t = 1

        if nJet >= 1:
            jet = myJets[0].p4()
            mllj = (lep1 + lep2 + jet).M()
            dR_mu1j = lep1.DeltaR(jet)
            dR_mu2j = lep2.DeltaR(jet)
            if abs(jet.Eta()) < 1.4442 and myMuons[0].pt > 23 and myMuons[0].pfRelIso04_all < 0.4:
                sel1_j = 1
            if myMuons[0].pt > 24 and myMuons[0].pfRelIso04_all < 0.4:
                sel2_j = 1
            if jet.Pt() > 31 and abs(jet.Eta()) < 1.5 and myMuons[0].pt > 24 and myMuons[0].pfRelIso04_all < 0.4:
                sel3_j = 1

        self.out.fillBranch("a_nMu",   nMuons)
        self.out.fillBranch("a_nPho",  nPho)
        self.out.fillBranch("a_nPho_t",nPho_t)
        self.out.fillBranch("a_nJet",  nJet)

        self.out.fillBranch("a_Pt_mu1",mu1Pt)

        self.out.fillBranch("a_mll",mll)
        self.out.fillBranch("a_mllg",mllg)
        self.out.fillBranch("a_dR_mu1g",dR_mu1g)
        self.out.fillBranch("a_dR_mu2g",dR_mu2g)
        self.out.fillBranch("a_sel1",sel1)
        self.out.fillBranch("a_sel2",sel2)
        self.out.fillBranch("a_sel3",sel3)

        self.out.fillBranch("a_mllg_t",mllg_t)
        self.out.fillBranch("a_dR_mu1g_t",dR_mu1g_t)
        self.out.fillBranch("a_dR_mu2g_t",dR_mu2g_t)
        self.out.fillBranch("a_sel1_t",sel1_t)
        self.out.fillBranch("a_sel2_t",sel2_t)
        self.out.fillBranch("a_sel3_t",sel3_t)
        self.out.fillBranch("a_sel4_t",sel4_t)


        self.out.fillBranch("a_mllj",mllj)
        self.out.fillBranch("a_dR_mu1j",dR_mu1j)
        self.out.fillBranch("a_dR_mu2j",dR_mu2j)
        self.out.fillBranch("a_sel1_j",sel1_j)
        self.out.fillBranch("a_sel2_j",sel2_j)
        self.out.fillBranch("a_sel3_j",sel3_j)

        self.out.fillBranch("a_trg_1Mu",trg_1Mu)
        self.out.fillBranch("a_trg_2Mu",trg_2Mu)
        self.out.fillBranch("a_trg_MuPho",trg_MuPho)


        return True


#preselection="Sum$(Muon_pt > 5.5 && Muon_softId) >=2"
#preselection="Sum$(Muon_pt > 22 && Muon_softId && Muon_pfRelIso04_all < 0.4) >=1"
#preselection="(Sum$(Photon_pt > 28 && abs(Photon_eta) < 1.6 && Photon_cutBasedBitmap >= 2) >= 1)"
preselection='''Sum$(Muon_pt > 5.5 && Muon_softId) >= 2 &&
 Sum$(Muon_pt > 22  && Muon_softId && Muon_pfRelIso04_all < 0.4) >=1'''

#preselection='''Sum$(Muon_pt > 5.5 && Muon_softId) >= 2 &&
# Sum$(Muon_pt > 22  && Muon_softId && Muon_pfRelIso04_all < 0.4) >=1 &&
# Sum$(Photon_pt > 28 && abs(Photon_eta) < 1.6) >= 1'''

files=["/eos/cms/store/user/andrey/f.root"]

files = [ "/user/andreypz/MuEGamma_2016G_81A6.root",
          "/user/andreypz/MuEGamma_2017F_C377.root",
          "/user/andreypz/MuEGamma_2018C_0CEF.root",
          "/user/andreypz/MuEGamma_2017B_5C2E.root",
          "/user/andreypz/ZH_HCC_ZLL_NanoV6_2016_FB6D.root",
          "/user/andreypz/ZH_Htautau_NanoV6_2017_3867.root"
          "/eos/cms/store/user/andrey/f_Nano14Dec2018_data.root",
          #"root://xrootd-cms.infn.it//store/data/Run2016G/MuonEG/NANOAOD/Nano25Oct2019-v1/40000/FB976070-45BA-5C44-B5C7-1355201581A6.root",
      ]

#this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles,runsAndLumis



myModules = []
if opt.muCorr:
    if opt.era==2016:
        myModules.append(muonScaleRes2016())
    if opt.era==2017:
        myModules.append(muonScaleRes2017())
    if opt.era==2018:
        myModules.append(muonScaleRes2018())

myModules.append(MyAnalysis(opt.era, opt.muCorr))

print("List of mudules to run:", myModules)

if opt.local:
    p=PostProcessor("./outDir", [files[opt.jobNum]], cut=preselection.replace('\n',' '), branchsel='keep_and_drop_input.txt', outputbranchsel='keep_and_drop_output.txt',
                    modules=myModules, histFileName="./outDir/histOut_"+str(opt.jobNum)+".root",histDirName="plots",maxEntries=opt.evt)
else:
    p=PostProcessor(".", inputFiles(), cut=preselection.replace('\n',' '), branchsel='keep_and_drop_input.txt', outputbranchsel='keep_and_drop_output.txt',
                    modules=myModules, provenance=True,fwkJobReport=True,jsonInput=runsAndLumis())

p.run()

#print "ls -lR:"
#os.system("ls -lR")
print "DONE"
