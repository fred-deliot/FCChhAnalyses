from __future__ import division
from heppy.framework.analyzer import Analyzer
from heppy.statistics.tree import Tree
from heppy.analyzers.ntuple import *
from heppy.particles.tlv.resonance import Resonance2 as Resonance
from heppy.particles.tlv.particle import Particle
from heppy.FCChhAnalyses.analyzers.TRFbtag import *

import math
import ROOT
from ROOT import *
import collections
#from array import array
import array
import os

class TreeProducer(Analyzer):

    def beginLoop(self, setup):
        super(TreeProducer, self).beginLoop(setup)
        self.rootfile = TFile('/'.join([self.dirName,
                                        'tree.root']),
                              'recreate')
        self.tree = Tree( 'events', '')
        
        self.tree.var('weight', float)
        self.tree.var('weight_0tagex', float)
        self.tree.var('weight_1tagex', float)
        self.tree.var('weight_2tagex', float)
        self.tree.var('weight_3tagex', float)
        self.tree.var('weight_4tagex', float)
        self.tree.var('missingET', float)
        self.tree.var('numberOfElectrons', int)
        self.tree.var('numberOfMuons', int)
	self.tree.var('numberOfLeptons',int)

	self.tree.var('numberOfJets', int)
	self.tree.var('numberOfBJets', int)

	self.tree.var('Ht',float)

        self.tree.var('Jet1_dR_lep', float)
        self.tree.var('Jet2_dR_lep', float)
        self.tree.var('Jet3_dR_lep', float)
        self.tree.var('Jet4_dR_lep', float)
        self.tree.var('Jet5_dR_lep', float)

	self.tree.var('SSee',int)
	self.tree.var('SSem',int)
	self.tree.var('SSmm',int)
	self.tree.var('SSeee',int)
	self.tree.var('SSemm',int)
	self.tree.var('SSmmm',int)

        bookParticle(self.tree, 'Electron1')
        bookParticle(self.tree, 'Electron2')
        bookParticle(self.tree, 'Electron3')

        bookParticle(self.tree, 'Muon1')
        bookParticle(self.tree, 'Muon2')
        bookParticle(self.tree, 'Muon3')

        bookParticle(self.tree, 'Jet1')
        bookParticle(self.tree, 'Jet2')
        bookParticle(self.tree, 'Jet3')
        bookParticle(self.tree, 'Jet4')
        bookParticle(self.tree, 'Jet5')


    def corrMET(self, jet1, pdg1, jet2, pdg2, jet3, pdg3, jet4, pdg4, met):
        dphi1 = abs(jet1.p4().DeltaPhi(met.p4()))
        dphi2 = abs(jet2.p4().DeltaPhi(met.p4()))
        dphi3 = abs(jet3.p4().DeltaPhi(met.p4()))
        dphi4 = abs(jet4.p4().DeltaPhi(met.p4()))

        metp4 = ROOT.TLorentzVector()
        px = met.p4().Px()
        py = met.p4().Py()

        dphi_min = min(dphi1, dphi2, dphi3, dphi4)
        if (dphi_min == dphi1):
            pz = jet1.p4().Pz()/2.
            e = math.sqrt(px**2 + py**2 + pz**2)
            metp4.SetPxPyPzE(px, py, pz, e) 
            jetcorr1   = Particle(pdg1, 0, jet1.p4() + metp4, 1)
            jetcorr2   = Particle(pdg2, 0, jet2.p4(), 1)
            jetcorr3   = Particle(pdg3, 0, jet3.p4(), 1)
            jetcorr4   = Particle(pdg4, 0, jet4.p4(), 1)
        elif (dphi_min == dphi2):
            pz = jet2.p4().Pz()/2.
            e = math.sqrt(px**2 + py**2 + pz**2)
            metp4.SetPxPyPzE(px, py, pz, e) 
            jetcorr1  = Particle(pdg1, 0, jet1.p4(), 1)
            jetcorr2  = Particle(pdg2, 0, jet2.p4() + metp4, 1)
            jetcorr3  = Particle(pdg3, 0, jet3.p4(), 1)
            jetcorr4  = Particle(pdg4, 0, jet4.p4(), 1)
        elif (dphi_min == dphi3):
            pz = jet3.p4().Pz()/2.
            e = math.sqrt(px**2 + py**2 + pz**2)
            metp4.SetPxPyPzE(px, py, pz, e)
            jetcorr1  = Particle(pdg1, 0, jet1.p4(), 1)
            jetcorr2  = Particle(pdg2, 0, jet2.p4(), 1)
            jetcorr3  = Particle(pdg3, 0, jet3.p4() + metp4, 1)
            jetcorr4  = Particle(pdg4, 0, jet4.p4(), 1)
        else :
            pz = jet4.p4().Pz()/2.
            e = math.sqrt(px**2 + py**2 + pz**2)
            metp4.SetPxPyPzE(px, py, pz, e)
            jetcorr1  = Particle(pdg1, 0, jet1.p4(), 1)
            jetcorr2  = Particle(pdg2, 0, jet2.p4(), 1)
            jetcorr3  = Particle(pdg3, 0, jet3.p4(), 1)
            jetcorr4  = Particle(pdg4, 0, jet4.p4() + metp4, 1)
        return jetcorr1,jetcorr2,jetcorr3,jetcorr4

    def fillMass(self, jet1, jet2, jet3, jet4):
        mj1j2j3j4 = ROOT.TLorentzVector()
        j1 = ROOT.TLorentzVector(); j2 = ROOT.TLorentzVector(); j3 = ROOT.TLorentzVector(); j4 = ROOT.TLorentzVector()
        j1.SetPtEtaPhiE(jet1.pt(), jet1.eta(), jet1.phi(), jet1.e())
        j2.SetPtEtaPhiE(jet2.pt(), jet2.eta(), jet2.phi(), jet2.e())
        j3.SetPtEtaPhiE(jet3.pt(), jet3.eta(), jet3.phi(), jet3.e())
        j4.SetPtEtaPhiE(jet4.pt(), jet4.eta(), jet4.phi(), jet4.e())
        mj1j2j3j4 = j1+j2+j3+j4
        return mj1j2j3j4.M()

     
    def process(self, event):
        self.tree.reset()
#        jets_trk08    = getattr(event, self.cfg_ana.jets_trk08_20)
        jets_pf04     = getattr(event, self.cfg_ana.jets_pf04)
        jets_pf04_pdg = event.jets_pf04_pdg #        jets_pf08     = getattr(event, self.cfg_ana.jets_pf08_30)

        electrons = getattr(event, self.cfg_ana.electrons)
        muons = getattr(event, self.cfg_ana.muons)

	Ht = 0

# pt cut on lepton collections in analysis.py)
# flage if this is a same-sign lepton or a multilepton events

	nlep = len(electrons) + len(muons)
	if (nlep<2):
	   return 

	SSee = 0
	SSem = 0
	SSmm = 0
	SSeee = 0
	SSeem = 0
	SSemm = 0
	SSmmm = 0
	if (len(electrons)>2):
	   SSeee = 1
	if (len(electrons)==2 and len(muons)>0):
	   SSeem = 1
	if (len(electrons)==2 and len(muons)==0):
	   if (electrons[0].pdgid()*electrons[1].pdgid()>0):
	     SSee =1
	if (len(electrons)==1 and len(muons)==1):
	   if (electrons[0].pdgid()*muons[0].pdgid()>0):
	     SSem = 1
	if (len(muons)==2 and len(electrons)==0):
	   if (muons[0].pdgid()*muons[1].pdgid()>0):
	      SSmm = 1
        if (len(electrons)==1 and len(muons)>1):
          SSemm = 1
	if (len(muons)>2 and SSeee==0 and SSeem==0 and SSemm==0):
	  SSmmm = 1
	
	for el in electrons:
	  Ht += el.pt()

	for mu in muons:
	  Ht += mu.pt()
	

# number of jets and b-jets (pt cut on jet collection in analysis.py)

        njets = len(jets_pf04)
        nbjets = 0
        for jet in jets_pf04:
	   Ht += jet.pt()
           if (jet.tags['bf'] > 0):
              nbjets +=1

# minimum number of jets and b-jets
	if (njets<5 or nbjets<2):
	  return

	self.tree.fill('numberOfLeptons' ,nlep)

	self.tree.fill('SSee' ,SSee)
        self.tree.fill('SSem' ,SSem)
        self.tree.fill('SSmm' ,SSmm)
        self.tree.fill('SSeee' ,SSeee)
        self.tree.fill('SSemm' ,SSemm)
        self.tree.fill('SSmmm' ,SSmm)

	fillParticle(self.tree, 'Jet1',jets_pf04[0])
        fillParticle(self.tree, 'Jet2',jets_pf04[1])
        fillParticle(self.tree, 'Jet3',jets_pf04[2])
        fillParticle(self.tree, 'Jet4',jets_pf04[3])
        fillParticle(self.tree, 'Jet5',jets_pf04[4])


        if ( len(electrons) >=1 ): fillParticle(self.tree, 'Electron1', electrons[0])
        if ( len(electrons) >=2 ): fillParticle(self.tree, 'Electron2', electrons[1])
        if ( len(electrons) >=3 ): fillParticle(self.tree, 'Electron3', electrons[2])

        if ( len(muons) >=1 ): fillParticle(self.tree, 'Muon1', muons[0])
        if ( len(muons) >=2 ): fillParticle(self.tree, 'Muon2', muons[1])
        if ( len(muons) >=3 ): fillParticle(self.tree, 'Muon3', muons[2])


	self.tree.fill('Ht',Ht)

        self.tree.fill('numberOfJets',njets)
        self.tree.fill('numberOfBJets',nbjets)

        self.tree.fill('weight' , event.weight )
        self.tree.fill('missingET', event.met.pt())
        self.tree.fill('numberOfElectrons', len(electrons))
        self.tree.fill('numberOfMuons', len(muons))



#        if ( len(jets_trk08)>=4 and  len(jets_pf08)>=4):
	if ( len(jets_pf04)>=5):
            #print len(jets_trk08), len(jets_pf08), len(electrons), len(muons)

            # compute the closest lepton/jet dR
            Jet_dR=[]
            for i_jet in range(5):
              the_Jet_dR = 999
              j = ROOT.TLorentzVector()
              j.SetPtEtaPhiE(jets_pf04[i_jet].pt(), jets_pf04[i_jet].eta(), jets_pf04[i_jet].phi(), jets_pf04[i_jet].e())
              for i_e in range(len(electrons)):
                e = ROOT.TLorentzVector()
                e.SetPtEtaPhiE(electrons[i_e].pt(), electrons[i_e].eta(), electrons[i_e].phi(), electrons[i_e].e())
                if j.DeltaR(e)<the_Jet_dR : the_Jet_dR=j.DeltaR(e)
              for i_m in range(len(muons)):
                m = ROOT.TLorentzVector()
                m.SetPtEtaPhiE(muons[i_m].pt(), muons[i_m].eta(), muons[i_m].phi(), muons[i_m].e())
                if j.DeltaR(m)<the_Jet_dR : the_Jet_dR=j.DeltaR(m)
              Jet_dR.append(the_Jet_dR)
            self.tree.fill('Jet1_dR_lep' , Jet_dR[0] )
            self.tree.fill('Jet2_dR_lep' , Jet_dR[1] )
            self.tree.fill('Jet3_dR_lep' , Jet_dR[2] )
            self.tree.fill('Jet4_dR_lep' , Jet_dR[3] )
            self.tree.fill('Jet5_dR_lep' , Jet_dR[4] )



            # TRF / truth b-tagging -> need at least 2 jets_pf04
            use_DELPHES=False
            weight_0tagex=0.
            weight_1tagex=0.
            weight_2tagex=0.
            weight_3tagex=0.
            weight_4tagex=0.
            jet=[]
            ipdg=0
            for i in range(len(jets_pf04)):
              if use_DELPHES==True:
                ipdg = jets_pf04[i].tags['flav']
                if ipdg!=4 and ipdg!=5 : ipdg=0
              else:
                ipdg = jets_pf04_pdg[i].flavour
              jet.append([jets_pf04[i],ipdg])
            weight_0tagex=getNbTagEx(0,jet,4)
            weight_1tagex=getNbTagEx(1,jet,4)
            weight_2tagex=getNbTagEx(2,jet,4)
            weight_3tagex=getNbTagEx(3,jet,4)
            weight_4tagex=getNbTagEx(4,jet,4)
            self.tree.fill('weight_0tagex', weight_0tagex)
            self.tree.fill('weight_1tagex', weight_1tagex)
            self.tree.fill('weight_2tagex', weight_2tagex)
            self.tree.fill('weight_3tagex', weight_3tagex)
            self.tree.fill('weight_4tagex', weight_4tagex)


        self.tree.tree.Fill()

    def write(self, setup):
        self.rootfile.Write()
        self.rootfile.Close()

