#!/usr/bin/env python3

import numpy as np
import argparse
import os
from array import array
import pyjetty.alihfjets.dev.hfjet.process.base.process_io_data_hf as hfdio
from pyjetty.mputils.mputils import perror, pinfo, pwarning
from pyjetty.mputils import treewriter, jet_analysis

import yaml
import fastjet as fj
import fjext
import fjtools
import fjcontrib

import ROOT
ROOT.gROOT.SetBatch(True)

class HFAnalysisInvMass(hfdio.HFAnalysis):
	def __init__(self, **kwargs):
		self.fout = None
		super(HFAnalysisInvMass, self).__init__(**kwargs)

		self.fout = ROOT.TFile(self.name+'.root', 'recreate')
		self.fout.cd()
		
		self.hNevents = ROOT.TH1F("hNevents", "hNevents", 2, array('f', [-0.5, 0.5, 1.5]) )
		self.hNevents.Sumw2()

		title = [ '#it{p}_{T}^{ch jet}', '#it{p}_{T}^{D^{0}}', 'Invmass', '#it{#lambda}_{#it{#alpha}=1}', 
			'#it{#lambda}_{#it{#alpha}=1.5}', '#it{#lambda}_{#it{#alpha}=2}', '#it{#lambda}_{#it{#alpha}=3}']

		self.nbins=array('i', (60, 60, 370, 160, 160, 160, 160))
		self.binlow=array('d', (0, 0, 1.7, 0, 0, 0, 0))
		self.binhigh=array('d', (60, 60, 2.07, 0.8, 0.8, 0.8, 0.8))
		self.fsparseJet=ROOT.THnSparseD("hsparsejet","hsparsejet;Jet_pt;D_pt;D_m;a10,a15,a20,a30",7,self.nbins,self.binlow,self.binhigh)
		self.fsparseJetvalue=array('d',(0,0,0,0,0,0,0))
		self.fsparseJet.Sumw2()
		
		self.fsparseJetSD=ROOT.THnSparseD("hsparsejet_groomed","hsparsejet_groomed;Jet_pt;D_pt;D_m;a10,a15,a20,a30",7,self.nbins,self.binlow,self.binhigh)
		self.fsparseJetSDvalue=array('d',(0,0,0,0,0,0,0))
		self.fsparseJetSD.Sumw2()

		for i in range(0,7):
			self.fsparseJet.GetAxis(i).SetTitle(title[i])
			self.fsparseJetSD.GetAxis(i).SetTitle(title[i])

		print(self.config_file)
	
		with open(self.config_file, 'r') as stream:
			self.config = yaml.load(stream,Loader=yaml.FullLoader)

	def process_events(self, df):
		print("processing event info")
		self.hNevents.Fill(1,df["ev_id"].nunique())
	
	def analysis(self, df):
		#print(df)

		m_array = np.full((self.df_tracks['ParticlePt'].values.size), 0.1396)	

		djmm = fjtools.DJetMatchMaker()
		djmm.set_ch_pt_eta_phi_m(self.df_tracks['ParticlePt'].values, self.df_tracks['ParticleEta'].values, self.df_tracks['ParticlePhi'].values, m_array)
		djmm.set_Ds_pt_eta_phi_m(df['pt_cand'].values, df['eta_cand'].values, df['phi_cand'].values, df['inv_mass'].values)
		djmm.set_daughters0_pt_eta_phi(df['pt_prong0'].values, df['eta_prong0'].values, df['phi_prong0'].values)
		djmm.set_daughters1_pt_eta_phi(df['pt_prong1'].values, df['eta_prong1'].values, df['phi_prong1'].values)


		#run for each D candidate to build jet
		for id0, d0 in enumerate(djmm.Ds):

			#daughter tracks matching
			_parts_and_ds = djmm.match(0.005, id0)
			#replacing daughter tracks with matched D0 candidate
			#including D0 	
			_parts_and_ds.push_back(d0)
	
			#jet reconstruction with D0 and charged particle
			jetR=0.4
			ja = jet_analysis.JetAnalysis(jet_R = jetR,jet_RecombinationScheme=fj.E_scheme, particle_eta_max=0.9, explicit_ghosts=False)
			ja.analyze_event(_parts_and_ds)
			
			if len(ja.jets) < 1:
				continue

			jets = ja.jets_as_psj_vector()

			#filtering D0 jets
			djets = djmm.filter_D0_jets(jets)


			if len(djets) > 0:
				j = djets[0]
				dcand = djmm.get_Dcand_in_jet(j)
				#number of constitutents > 1
				#print("number of constituents ", len(j.constituents()))
				#if len(j.constituents())<=2:
				#	continue	
				self.fsparseJetvalue[0] = j.perp()
				self.fsparseJetvalue[1] = dcand[0].perp()
				self.fsparseJetvalue[2] = dcand[0].m()
				self.fsparseJetvalue[3] = fjext.lambda_beta_kappa(j,  1.0, 1.0 ,0.4)
				self.fsparseJetvalue[4] = fjext.lambda_beta_kappa(j,  1.5, 1.0 ,0.4)
				self.fsparseJetvalue[5] = fjext.lambda_beta_kappa(j,  2.0, 1.0 ,0.4)
				self.fsparseJetvalue[6] = fjext.lambda_beta_kappa(j,  3.0, 1.0 ,0.4)
									
				self.fsparseJet.Fill(self.fsparseJetvalue)
	

				#jets with the SD groomed axis################################ 	
				gshop = fjcontrib.GroomerShop(j, jetR, fj.cambridge_algorithm)
				sd_zcut=0.2
				sd_beta=0
				jet_groomed_lund = gshop.soft_drop(sd_beta,sd_zcut,jetR)
				jet_groomed = jet_groomed_lund.pair()
				
				#ang with the SD groomed axis################################  
				a10 = fjext.lambda_beta_kappa(jet_groomed,  1.0, 1.0 ,0.4)
				a15 = fjext.lambda_beta_kappa(jet_groomed,  1.5, 1.0 ,0.4)
				a20 = fjext.lambda_beta_kappa(jet_groomed,  2.0, 1.0 ,0.4)
				a30 = fjext.lambda_beta_kappa(jet_groomed,  3.0, 1.0 ,0.4)
				
				#if SD condition fails put in the negative bin
				if jet_groomed_lund.Delta() < 0:
					a10 = a15 = a20 = a30 = -0.005
				
				#else fill the THnSparse
				self.fsparseJetSDvalue[0] = jet_groomed.perp()
				self.fsparseJetSDvalue[1] = dcand[0].perp()
				self.fsparseJetSDvalue[2] = dcand[0].m()
				self.fsparseJetSDvalue[3] = a10
				self.fsparseJetSDvalue[4] = a15
				self.fsparseJetSDvalue[5] = a20
				self.fsparseJetSDvalue[6] = a30

				self.fsparseJetSD.Fill(self.fsparseJetSDvalue)
	
			if len(djets) > 1:
				perror("more than one jet per D candidate?")

		return True
	

	def finalize(self):
		self.hNevents.Write()
		self.fsparseJet.Write()	
		self.fsparseJetSD.Write()
		self.fout.Close()
		pinfo(self.fout.GetName(), 'written.')


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='D0 analysis on alice data', prog=os.path.basename(__file__))
	parser.add_argument('-c', '--configFile', help='Path of config file for analysis', type=str,
			metavar='configFile',default='/home/preeti/analysis/pyjetty/pyjetty/alihfjets/dev/hfjet/config/hf_ang/configcuts_ptbin.yaml', required=True)
	parser.add_argument('-f', '--inputFile', action='store',
                      type=str, metavar='inputFile',
                      default='AnalysisResults.root',
                      help='Path of ROOT file containing TTrees')
	parser.add_argument('-o', '--output', help="output name / file name in the end", type=str, default='test_hfana')
	args = parser.parse_args()
	
	print('Configuring...')
	if not os.path.exists(args.configFile):
		print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
		sys.exit(0)

	print('inputFile: \'{0}\''.format(args.inputFile))
	print('----------------------------------------------------------------')
	if not os.path.exists(args.inputFile):
		print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
		sys.exit(0)

	hfaio = hfdio.HFAnalysisIO()
	
	hfa = HFAnalysisInvMass(config_file=args.configFile, name = args.output)
	print('applying event selection cuts')
	hfa.event_selection.add_selection_range_abs('z_vtx_reco', 10)
	hfa.event_selection.add_selection_equal('is_ev_rej', 0)	
	print('applying D0 selection cuts')
	#topomatic cut suggested by D2H
	hfa.d0_selection.add_selection_range_abs('max_norm_d0d0exp',2)

	hfa.d0_selection.add_selection_range('pt_cand', 2, 1e3)
	hfa.d0_selection.add_selection_range_abs('eta_cand', 0.8)
	hfa.d0_selection.add_selection_nsig('nsigTPC_Pi_0', 'nsigTOF_Pi_0', 'nsigTPC_K_1', 'nsigTOF_K_1', 'nsigTPC_Pi_1', 'nsigTOF_Pi_1', 'nsigTPC_K_0', 'nsigTOF_K_0', 3, -900)
		

	hfaio.add_analysis(hfa)
	
	hfaio.execute_analyses_on_inputfile(args.inputFile)
	
	hfa.finalize()
