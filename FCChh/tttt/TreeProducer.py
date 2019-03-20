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

        bookParticle(self.tree, 'Electron1')
        bookParticle(self.tree, 'Electron2')
        bookParticle(self.tree, 'Electron3')

        bookParticle(self.tree, 'Muon1')
        bookParticle(self.tree, 'Muon2')
        bookParticle(self.tree, 'Muon3')

	
	### variables for SSML channel 

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
	self.tree.var('SSeem',int)
	self.tree.var('SSemm',int)
	self.tree.var('SSmmm',int)
	self.tree.var('SS4l',int)

        bookParticle(self.tree, 'Jet1')
        bookParticle(self.tree, 'Jet2')
        bookParticle(self.tree, 'Jet3')
        bookParticle(self.tree, 'Jet4')
        bookParticle(self.tree, 'Jet5')


	### variables for lepton+jets channel

	 #trk08 no SD
        self.tree.var('Jet1_trk08_tau1', float)       
        self.tree.var('Jet1_trk08_tau2', float)
        self.tree.var('Jet1_trk08_tau3', float)
        self.tree.var('Jet2_trk08_tau1', float)
        self.tree.var('Jet2_trk08_tau2', float)
        self.tree.var('Jet2_trk08_tau3', float)
        self.tree.var('Jet3_trk08_tau1', float)
        self.tree.var('Jet3_trk08_tau2', float)
        self.tree.var('Jet3_trk08_tau3', float)
        self.tree.var('Jet4_trk08_tau1', float)
        self.tree.var('Jet4_trk08_tau2', float)
        self.tree.var('Jet4_trk08_tau3', float)
        self.tree.var('Jet1_trk08_tau32', float)
        self.tree.var('Jet1_trk08_tau31', float)
        self.tree.var('Jet1_trk08_tau21', float)
        self.tree.var('Jet2_trk08_tau32', float)
        self.tree.var('Jet2_trk08_tau31', float)
        self.tree.var('Jet2_trk08_tau21', float)
        self.tree.var('Jet3_trk08_tau32', float)
        self.tree.var('Jet3_trk08_tau31', float)
        self.tree.var('Jet3_trk08_tau21', float)
        self.tree.var('Jet4_trk08_tau32', float)
        self.tree.var('Jet4_trk08_tau31', float)
        self.tree.var('Jet4_trk08_tau21', float)

        bookParticle(self.tree, 'Jet1_trk08_Corr_MetCorr')
        bookParticle(self.tree, 'Jet2_trk08_Corr_MetCorr')
        bookParticle(self.tree, 'Jet3_trk08_Corr_MetCorr')
        bookParticle(self.tree, 'Jet4_trk08_Corr_MetCorr')

        #trk08 SD
        bookParticle(self.tree, 'Jet1_trk08_SD_Corr')
        bookParticle(self.tree, 'Jet2_trk08_SD_Corr')
        bookParticle(self.tree, 'Jet3_trk08_SD_Corr')
        bookParticle(self.tree, 'Jet4_trk08_SD_Corr')

        bookParticle(self.tree, 'Jet1_trk08_SD_MetCorr')
        bookParticle(self.tree, 'Jet2_trk08_SD_MetCorr')
        bookParticle(self.tree, 'Jet3_trk08_SD_MetCorr')
        bookParticle(self.tree, 'Jet4_trk08_SD_MetCorr')
        
        bookParticle(self.tree, 'Jet1_trk08_SD_Corr_MetCorr')
        bookParticle(self.tree, 'Jet2_trk08_SD_Corr_MetCorr')
        bookParticle(self.tree, 'Jet3_trk08_SD_Corr_MetCorr')
        bookParticle(self.tree, 'Jet4_trk08_SD_Corr_MetCorr')

        self.tree.var('Jet1_trk08_SD_Corr_m', float)
        self.tree.var('Jet2_trk08_SD_Corr_m', float)
        self.tree.var('Jet3_trk08_SD_Corr_m', float)
        self.tree.var('Jet4_trk08_SD_Corr_m', float)

        self.tree.var('rapiditySeparation_trk08_j1j2', float)
        self.tree.var('rapiditySeparation_trk08_j1j3', float)
        self.tree.var('rapiditySeparation_trk08_j1j4', float)
        self.tree.var('transverseMomentumAsymmetry_trk08_j1j2', float)
        self.tree.var('transverseMomentumAsymmetry_trk08_j1j3', float)
        self.tree.var('transverseMomentumAsymmetry_trk08_j1j4', float)

        self.tree.var('Mj1j2j3j4_trk08', float)
        self.tree.var('Mj1j2j3j4_trk08_Corr', float)
        self.tree.var('Mj1j2j3j4_trk08_MetCorr', float)
        self.tree.var('Mj1j2j3j4_trk08_Corr_MetCorr', float)

        self.tree.var('Mj1j2j3j4_pf08', float)
        self.tree.var('Mj1j2j3j4_pf08_MetCorr', float)

        self.tree.var('Jet1_trk08_dR_lep', float)
        self.tree.var('Jet2_trk08_dR_lep', float)
        self.tree.var('Jet3_trk08_dR_lep', float)
        self.tree.var('Jet4_trk08_dR_lep', float)



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
        jets_trk08    = getattr(event, self.cfg_ana.jets_trk08_20)
        jets_pf04     = getattr(event, self.cfg_ana.jets_pf04)
        jets_pf04_pdg = event.jets_pf04_pdg
	jets_pf08     = getattr(event, self.cfg_ana.jets_pf08_30)

        electrons = getattr(event, self.cfg_ana.electrons)
        muons = getattr(event, self.cfg_ana.muons)

	Ht = 0
        SSee = 0
        SSem = 0
        SSmm = 0
        SSeee = 0
        SSeem = 0
        SSemm = 0
        SSmmm = 0
	SS4l = 0

	# pt cut on lepton collections in analysis.py

	#####################################################################################
	##### selection for the SSML channel 
        #####################################################################################


	# flag if this is a same-sign lepton or a multilepton events

	nlep = len(electrons) + len(muons)
	if (nlep>=2):

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

#	  sum = SSeee+SSeem+SSemm+SSmmm+SSee+SSem+SSmm
#	  if (sum>1):
#	     print  'SSeee:',SSeee,' SSeem:',SSeem,' SSemm:',SSemm,' SSmmm:',SSmmm, ' SSee:',SSee,' SSem:', SSem, ' SSmm:',SSmm

	  if (nlep>=4):
	     SS4l = 1
	
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
	  if (njets>=5 and nbjets>=2):

	     self.tree.fill('numberOfLeptons' ,nlep)

	     self.tree.fill('SSee' ,SSee)
             self.tree.fill('SSem' ,SSem)
             self.tree.fill('SSmm' ,SSmm)
             self.tree.fill('SSeee' ,SSeee)
	     self.tree.fill('SSeem' ,SSeem)
             self.tree.fill('SSemm' ,SSemm)
             self.tree.fill('SSmmm' ,SSmmm)
	     self.tree.fill('SS4l' ,SS4l)

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


             pdg1 = 0
             pdg2 = 0
             pdg3 = 0
             pdg4 = 0
	     pdg5 = 0
	     if (jets_pf04[0].tags['bf'] > 0):
	        pdg1 = 5
             if (jets_pf04[1].tags['bf'] > 0):
                pdg2 = 5
             if (jets_pf04[2].tags['bf'] > 0):
                pdg3 = 5
             if (jets_pf04[3].tags['bf'] > 0):
                pdg4 = 5
             if (jets_pf04[4].tags['bf'] > 0):
                pdg5 = 5

	     myjet1 = Particle(pdg1, 0, jets_pf04[0].p4(), 1)
             myjet2 = Particle(pdg2, 0, jets_pf04[1].p4(), 1)
             myjet3 = Particle(pdg3, 0, jets_pf04[2].p4(), 1)
             myjet4 = Particle(pdg4, 0, jets_pf04[3].p4(), 1)
             myjet5 = Particle(pdg5, 0, jets_pf04[4].p4(), 1)

             fillParticle(self.tree, 'Jet1',myjet1)
             fillParticle(self.tree, 'Jet2',myjet2)
             fillParticle(self.tree, 'Jet3',myjet3)
             fillParticle(self.tree, 'Jet4',myjet4)
             fillParticle(self.tree, 'Jet5',myjet5)

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


        #####################################################################################
        ##### selection for the lepton+jets channel
        #####################################################################################


        if ( len(jets_trk08)>=4 and  len(jets_pf08)>=4):
          if jets_trk08[0].pt()>200.:
            #print len(jets_trk08), len(jets_pf08), len(electrons), len(muons)

            j1 = ROOT.TLorentzVector(); j2 = ROOT.TLorentzVector(); j3 = ROOT.TLorentzVector(); j4 = ROOT.TLorentzVector()
            j1.SetPtEtaPhiE(jets_trk08[0].pt(), jets_trk08[0].eta(), jets_trk08[0].phi(), jets_trk08[0].e())
            j2.SetPtEtaPhiE(jets_trk08[1].pt(), jets_trk08[1].eta(), jets_trk08[1].phi(), jets_trk08[1].e())
            j3.SetPtEtaPhiE(jets_trk08[2].pt(), jets_trk08[2].eta(), jets_trk08[2].phi(), jets_trk08[2].e())
            j4.SetPtEtaPhiE(jets_trk08[3].pt(), jets_trk08[3].eta(), jets_trk08[3].phi(), jets_trk08[3].e())

            # compute the closest lepton/jet dR
            Jet_dR=[]
            for i_jet in range(4):
              the_Jet_dR = 999
              j = ROOT.TLorentzVector()
              j.SetPtEtaPhiE(jets_trk08[i_jet].pt(), jets_trk08[i_jet].eta(), jets_trk08[i_jet].phi(), jets_trk08[i_jet].e())
              for i_e in range(len(electrons)):
                e = ROOT.TLorentzVector()
                e.SetPtEtaPhiE(electrons[i_e].pt(), electrons[i_e].eta(), electrons[i_e].phi(), electrons[i_e].e())
                if j.DeltaR(e)<the_Jet_dR : the_Jet_dR=j.DeltaR(e)
              for i_m in range(len(muons)):
                m = ROOT.TLorentzVector()
                m.SetPtEtaPhiE(muons[i_m].pt(), muons[i_m].eta(), muons[i_m].phi(), muons[i_m].e())
                if j.DeltaR(m)<the_Jet_dR : the_Jet_dR=j.DeltaR(m)
              Jet_dR.append(the_Jet_dR)
            self.tree.fill('Jet1_trk08_dR_lep' , Jet_dR[0] )
            self.tree.fill('Jet2_trk08_dR_lep' , Jet_dR[1] )
            self.tree.fill('Jet3_trk08_dR_lep' , Jet_dR[2] )
            self.tree.fill('Jet4_trk08_dR_lep' , Jet_dR[3] )

            self.tree.fill('weight' , event.weight )
            self.tree.fill('missingET', event.met.pt())
            self.tree.fill('numberOfElectrons', len(electrons))
            self.tree.fill('numberOfMuons', len(muons))

            self.tree.fill('rapiditySeparation_trk08_j1j2', abs(jets_trk08[0].eta() - jets_trk08[1].eta()))
            self.tree.fill('rapiditySeparation_trk08_j1j3', abs(jets_trk08[0].eta() - jets_trk08[2].eta()))
            self.tree.fill('rapiditySeparation_trk08_j1j4', abs(jets_trk08[0].eta() - jets_trk08[3].eta()))
            self.tree.fill('transverseMomentumAsymmetry_trk08_j1j2', (jets_trk08[0].pt() - jets_trk08[1].pt())/(jets_trk08[0].pt() + jets_trk08[1].pt()))
            self.tree.fill('transverseMomentumAsymmetry_trk08_j1j3', (jets_trk08[0].pt() - jets_trk08[2].pt())/(jets_trk08[0].pt() + jets_trk08[2].pt()))
            self.tree.fill('transverseMomentumAsymmetry_trk08_j1j4', (jets_trk08[0].pt() - jets_trk08[3].pt())/(jets_trk08[0].pt() + jets_trk08[3].pt()))

            self.tree.fill('Jet1_trk08_tau1' , jets_trk08[0].tau1 )
            self.tree.fill('Jet1_trk08_tau2' , jets_trk08[0].tau2 )
            self.tree.fill('Jet1_trk08_tau3' , jets_trk08[0].tau3 )
            self.tree.fill('Jet2_trk08_tau1' , jets_trk08[1].tau1 )
            self.tree.fill('Jet2_trk08_tau2' , jets_trk08[1].tau2 )
            self.tree.fill('Jet2_trk08_tau3' , jets_trk08[1].tau3 )
            self.tree.fill('Jet3_trk08_tau1' , jets_trk08[2].tau1 )
            self.tree.fill('Jet3_trk08_tau2' , jets_trk08[2].tau2 )
            self.tree.fill('Jet3_trk08_tau3' , jets_trk08[2].tau3 )
            self.tree.fill('Jet4_trk08_tau1' , jets_trk08[3].tau1 )
            self.tree.fill('Jet4_trk08_tau2' , jets_trk08[3].tau2 )
            self.tree.fill('Jet4_trk08_tau3' , jets_trk08[3].tau3 )

            Jet1_trk08_tau31 = -999.0
            Jet1_trk08_tau21 = -999.0
            Jet1_trk08_tau32 = -999.0
            Jet2_trk08_tau31 = -999.0
            Jet2_trk08_tau21 = -999.0
            Jet2_trk08_tau32 = -999.0
            Jet3_trk08_tau31 = -999.0
            Jet3_trk08_tau21 = -999.0
            Jet3_trk08_tau32 = -999.0
            Jet4_trk08_tau31 = -999.0
            Jet4_trk08_tau21 = -999.0
            Jet4_trk08_tau32 = -999.0

            if (jets_trk08[0].tau1 != 0.0):
                Jet1_trk08_tau31 = jets_trk08[0].tau3/jets_trk08[0].tau1
                Jet1_trk08_tau21 = jets_trk08[0].tau2/jets_trk08[0].tau1 
            if (jets_trk08[0].tau2 != 0.0):
                Jet1_trk08_tau32 = jets_trk08[0].tau3/jets_trk08[0].tau2

            if (jets_trk08[1].tau1 != 0.0):
                Jet2_trk08_tau31 = jets_trk08[1].tau3/jets_trk08[1].tau1
                Jet2_trk08_tau21 = jets_trk08[1].tau2/jets_trk08[1].tau1
            if (jets_trk08[1].tau2 != 0.0):
                Jet2_trk08_tau32 = jets_trk08[1].tau3/jets_trk08[1].tau2

            if (jets_trk08[2].tau1 != 0.0):
                Jet3_trk08_tau31 = jets_trk08[2].tau3/jets_trk08[2].tau1
                Jet3_trk08_tau21 = jets_trk08[2].tau2/jets_trk08[2].tau1
            if (jets_trk08[2].tau2 != 0.0):
                Jet3_trk08_tau32 = jets_trk08[2].tau3/jets_trk08[2].tau2

            if (jets_trk08[3].tau1 != 0.0):
                Jet4_trk08_tau31 = jets_trk08[3].tau3/jets_trk08[3].tau1
                Jet4_trk08_tau21 = jets_trk08[3].tau2/jets_trk08[3].tau1
            if (jets_trk08[3].tau2 != 0.0):
                Jet4_trk08_tau32 = jets_trk08[3].tau3/jets_trk08[3].tau2

            self.tree.fill('Jet1_trk08_tau31', Jet1_trk08_tau31)
            self.tree.fill('Jet1_trk08_tau21', Jet1_trk08_tau21)
            self.tree.fill('Jet1_trk08_tau32', Jet1_trk08_tau32)
            self.tree.fill('Jet2_trk08_tau31', Jet2_trk08_tau31)
            self.tree.fill('Jet2_trk08_tau21', Jet2_trk08_tau21)
            self.tree.fill('Jet2_trk08_tau32', Jet2_trk08_tau32)
            self.tree.fill('Jet3_trk08_tau31', Jet3_trk08_tau31)
            self.tree.fill('Jet3_trk08_tau21', Jet3_trk08_tau21)
            self.tree.fill('Jet3_trk08_tau32', Jet3_trk08_tau32)
            self.tree.fill('Jet4_trk08_tau31', Jet4_trk08_tau31)
            self.tree.fill('Jet4_trk08_tau21', Jet4_trk08_tau21)
            self.tree.fill('Jet4_trk08_tau32', Jet4_trk08_tau32)

	    # here is btag, need matching in DR
            Jet1_trk08_dR_pf04 = 999
            Jet2_trk08_dR_pf04 = 999
            Jet3_trk08_dR_pf04 = 999
            Jet4_trk08_dR_pf04 = 999
	    for j in jets_pf04:
                pf04= ROOT.TLorentzVector()
                pf04.SetPtEtaPhiE(j.pt(), j.eta(), j.phi(), j.e())
                if j.tags['bf'] > 0:
                    if pf04.DeltaR(j1)<Jet1_trk08_dR_pf04:
                        Jet1_trk08_dR_pf04=pf04.DeltaR(j1)
                    if pf04.DeltaR(j2)<Jet2_trk08_dR_pf04:
                        Jet2_trk08_dR_pf04=pf04.DeltaR(j2)
                    if pf04.DeltaR(j3)<Jet3_trk08_dR_pf04:
                        Jet3_trk08_dR_pf04=pf04.DeltaR(j3)
                    if pf04.DeltaR(j4)<Jet4_trk08_dR_pf04:
                        Jet4_trk08_dR_pf04=pf04.DeltaR(j4)

            pdg1 = 0
            pdg2 = 0
            pdg3 = 0
            pdg4 = 0
            if Jet1_trk08_dR_pf04 < 0.3:
                pdg1 = 5
            if Jet2_trk08_dR_pf04 < 0.3:
                pdg2 = 5
            if Jet3_trk08_dR_pf04 < 0.3:
                pdg3 = 5
            if Jet4_trk08_dR_pf04 < 0.3:
                pdg4 = 5


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

            ######################
            # trkjet08 mass info #
            ######################

            #MATCHING PF08 and trk08 for CORRECTION
            Jet1_trk08_dR_pf08 = 999
            Jet2_trk08_dR_pf08 = 999
            Jet3_trk08_dR_pf08 = 999
            Jet4_trk08_dR_pf08 = 999
            Jet1_pf08 = None
            Jet2_pf08 = None
            Jet3_pf08 = None
            Jet4_pf08 = None
            for j in jets_pf08:
                pf08= ROOT.TLorentzVector()
                pf08.SetPtEtaPhiE(j.pt(), j.eta(), j.phi(), j.e())
                if pf08.DeltaR(j1)<Jet1_trk08_dR_pf08:
                    Jet1_trk08_dR_pf08=pf08.DeltaR(j1)
                    Jet1_pf08=j
                if pf08.DeltaR(j2)<Jet2_trk08_dR_pf08:
                    Jet2_trk08_dR_pf08=pf08.DeltaR(j2)
                    Jet2_pf08=j
                if pf08.DeltaR(j3)<Jet3_trk08_dR_pf08:
                    Jet3_trk08_dR_pf08=pf08.DeltaR(j3)
                    Jet3_pf08=j
                if pf08.DeltaR(j4)<Jet4_trk08_dR_pf08:
                    Jet4_trk08_dR_pf08=pf08.DeltaR(j4)
                    Jet4_pf08=j

            #print 'jet1 dr ',Jet1_trk08_dR_pf08,'  pf08   ',Jet1_pf08,'  trk08  ',jets_trk08[0]
            #print 'jet2 dr ',Jet2_trk08_dR_pf08,'  pf08   ',Jet2_pf08,'  trk08  ',jets_trk08[1]

	    corr1 = Jet1_pf08.p4().Pt()/j1.Pt()
	    corr2 = Jet2_pf08.p4().Pt()/j2.Pt()
            corr3 = Jet3_pf08.p4().Pt()/j3.Pt()
            corr4 = Jet4_pf08.p4().Pt()/j4.Pt()

            #print 'corr 1  ',corr1,'   corr2  ',corr2
            #NORMAL TRK08 SD corrected jet
	    p4sd1 = ROOT.TLorentzVector(); p4sd2 = ROOT.TLorentzVector(); p4sd3 = ROOT.TLorentzVector(); p4sd4 = ROOT.TLorentzVector()
	    p4sd1.SetPtEtaPhiM(jets_trk08[0].subjetsSoftDrop[0].p4().Pt()*corr1, 
	    			jets_trk08[0].eta(), 
				jets_trk08[0].phi(), 
				jets_trk08[0].subjetsSoftDrop[0].p4().M()*corr1)
	    
	    p4sd2.SetPtEtaPhiM(jets_trk08[1].subjetsSoftDrop[0].p4().Pt()*corr2, 
	    			jets_trk08[1].eta(), 
				jets_trk08[1].phi(), 
				jets_trk08[1].subjetsSoftDrop[0].p4().M()*corr2)

            p4sd3.SetPtEtaPhiM(jets_trk08[2].subjetsSoftDrop[0].p4().Pt()*corr3,
                                jets_trk08[2].eta(),
                                jets_trk08[2].phi(),
                                jets_trk08[2].subjetsSoftDrop[0].p4().M()*corr3)

            p4sd4.SetPtEtaPhiM(jets_trk08[3].subjetsSoftDrop[0].p4().Pt()*corr4,
                                jets_trk08[3].eta(),
                                jets_trk08[3].phi(),
                                jets_trk08[3].subjetsSoftDrop[0].p4().M()*corr4)
	    
            sdjet1_corr = Particle(pdg1, 0, p4sd1, 1)
            sdjet2_corr = Particle(pdg2, 0, p4sd2, 1)
            sdjet3_corr = Particle(pdg3, 0, p4sd3, 1)
            sdjet4_corr = Particle(pdg4, 0, p4sd4, 1)
            fillParticle(self.tree, 'Jet1_trk08_SD_Corr', sdjet1_corr)
            fillParticle(self.tree, 'Jet2_trk08_SD_Corr', sdjet2_corr)
            fillParticle(self.tree, 'Jet3_trk08_SD_Corr', sdjet3_corr)
            fillParticle(self.tree, 'Jet4_trk08_SD_Corr', sdjet4_corr)

            #NORMAL TRK08 SD jet
	    #sdjet1 = Particle(pdg1, 0, jets_trk08[0].subjetsSoftDrop[0].p4(), 1)
            #sdjet2 = Particle(pdg2, 0, jets_trk08[1].subjetsSoftDrop[0].p4(), 1)
            #fillParticle(self.tree, 'Jet1_trk08_SD', sdjet1)
            #fillParticle(self.tree, 'Jet2_trk08_SD', sdjet2)

            #CORRECTED TRK08 jet
	    p4jet1_corr = ROOT.TLorentzVector(); p4jet2_corr = ROOT.TLorentzVector(); p4jet3_corr = ROOT.TLorentzVector(); p4jet4_corr = ROOT.TLorentzVector()
            p4jet1_corr.SetPtEtaPhiM(jets_trk08[0].pt()*corr1, jets_trk08[0].eta(), jets_trk08[0].phi(), jets_trk08[0].m()*corr1)
	    p4jet2_corr.SetPtEtaPhiM(jets_trk08[1].pt()*corr2, jets_trk08[1].eta(), jets_trk08[1].phi(), jets_trk08[1].m()*corr2)
            p4jet3_corr.SetPtEtaPhiM(jets_trk08[2].pt()*corr3, jets_trk08[2].eta(), jets_trk08[2].phi(), jets_trk08[2].m()*corr3)
            p4jet4_corr.SetPtEtaPhiM(jets_trk08[3].pt()*corr4, jets_trk08[3].eta(), jets_trk08[3].phi(), jets_trk08[3].m()*corr4)

            jet1_corr = Particle(pdg1, 0, p4jet1_corr, 1)
            jet2_corr = Particle(pdg2, 0, p4jet2_corr, 1)
            jet3_corr = Particle(pdg3, 0, p4jet3_corr, 1)
            jet4_corr = Particle(pdg4, 0, p4jet4_corr, 1)
            #fillParticle(self.tree, 'Jet1_trk08_Corr', jet1_corr)
            #fillParticle(self.tree, 'Jet2_trk08_Corr', jet2_corr)

 
            # associate MET to one jet or another based on softdrop
            sdjetmet1, sdjetmet2, sdjetmet3, sdjetmet4 = self.corrMET(jets_trk08[0].subjetsSoftDrop[0], pdg1, jets_trk08[1].subjetsSoftDrop[0], pdg2, jets_trk08[2].subjetsSoftDrop[0], pdg3, jets_trk08[3].subjetsSoftDrop[0], pdg4, event.met)
            fillParticle(self.tree, 'Jet1_trk08_SD_MetCorr', sdjetmet1)
            fillParticle(self.tree, 'Jet2_trk08_SD_MetCorr', sdjetmet2)
            fillParticle(self.tree, 'Jet3_trk08_SD_MetCorr', sdjetmet3)
            fillParticle(self.tree, 'Jet4_trk08_SD_MetCorr', sdjetmet4)


            sdjetmet1, sdjetmet2, sdjetmet3, sdjetmet4 = self.corrMET(sdjet1_corr, pdg1, sdjet2_corr, pdg2, sdjet3_corr, pdg3, sdjet4_corr, pdg4, event.met)
            fillParticle(self.tree, 'Jet1_trk08_SD_Corr_MetCorr', sdjetmet1)
            fillParticle(self.tree, 'Jet2_trk08_SD_Corr_MetCorr', sdjetmet2)
            fillParticle(self.tree, 'Jet3_trk08_SD_Corr_MetCorr', sdjetmet3)
            fillParticle(self.tree, 'Jet4_trk08_SD_Corr_MetCorr', sdjetmet4)



            # form masses
            self.tree.fill( 'Mj1j2j3j4_trk08',self.fillMass(jets_trk08[0],jets_trk08[1],jets_trk08[2],jets_trk08[3]))
            self.tree.fill( 'Mj1j2j3j4_trk08_Corr',self.fillMass(jet1_corr,jet2_corr,jet3_corr,jet4_corr))
            jetmet1, jetmet2, jetmet3, jetmet4 = self.corrMET(jets_trk08[0], pdg1, jets_trk08[1], pdg2, jets_trk08[2], pdg3, jets_trk08[3], pdg4, event.met)
            self.tree.fill( 'Mj1j2j3j4_trk08_MetCorr',self.fillMass(jetmet1,jetmet2,jetmet3,jetmet4))
            #fillParticle(self.tree, 'Jet1_trk08_MetCorr', jetmet1)

            jetmet1, jetmet2, jetmet3, jetmet4 = self.corrMET(jet1_corr, pdg1, jet2_corr, pdg2, jet3_corr, pdg3, jet4_corr, pdg4, event.met)
            self.tree.fill( 'Mj1j2j3j4_trk08_Corr_MetCorr',self.fillMass(jetmet1,jetmet2,jetmet3,jetmet4))
            fillParticle(self.tree, 'Jet1_trk08_Corr_MetCorr', jetmet1)
            fillParticle(self.tree, 'Jet2_trk08_Corr_MetCorr', jetmet2)
            fillParticle(self.tree, 'Jet3_trk08_Corr_MetCorr', jetmet3)
            fillParticle(self.tree, 'Jet4_trk08_Corr_MetCorr', jetmet4)

            self.tree.fill( 'Mj1j2j3j4_pf08', self.fillMass(jets_pf08[0],jets_pf08[1],jets_pf08[2],jets_pf08[3]))
            jetmet1, jetmet2, jetmet3, jetmet4 = self.corrMET(jets_pf08[0], pdg1, jets_pf08[1], pdg2, jets_pf08[2], pdg3, jets_pf08[3], pdg4, event.met)
            self.tree.fill( 'Mj1j2j3j4_pf08_MetCorr', self.fillMass(jetmet1,jetmet2,jetmet3,jetmet4))


            if ( len(electrons) >=1 ): fillParticle(self.tree, 'Electron1', electrons[0])
            if ( len(electrons) >=2 ): fillParticle(self.tree, 'Electron2', electrons[1])
            if ( len(electrons) >=3 ): fillParticle(self.tree, 'Electron3', electrons[2])

            if ( len(muons) >=1 ): fillParticle(self.tree, 'Muon1', muons[0])
            if ( len(muons) >=2 ): fillParticle(self.tree, 'Muon2', muons[1])
            if ( len(muons) >=3 ): fillParticle(self.tree, 'Muon3', muons[2])


        self.tree.tree.Fill()

    def write(self, setup):
        self.rootfile.Write()
        self.rootfile.Close()

