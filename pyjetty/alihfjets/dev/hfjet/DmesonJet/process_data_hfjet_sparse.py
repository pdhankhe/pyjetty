#!/usr/bin/env python3

import numpy as np
import argparse
import os
from array import array
import pyjetty.alihfjets.dev.hfjet.process_io_data_hf as hfdio
from pyjetty.mputils import perror, pinfo, pwarning, treewriter, jet_analysis

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

		#self.nbins=array('i',(60,60,60,350,60,60))
		#self.binlow=array('d',(0,0,0,1.72,0,0))
		#self.binhigh=array('d',(60,60,60,2.07,0.6,0.6))

		#self.fsparseJet=ROOT.THnSparseD("hsparsejet","hsparsejet;Jet_pt;JetWTA_pt;D_pt;D_m;dR;dR_WTA",6,self.nbins,self.binlow,self.binhigh)
		#self.fsparseJetvalue=array('d',(0,0,0,0,0,0))

		self.nbins=array('i',(60,60,350,160,160,160,160))
		self.binlow=array('d',(0,0,1.72,0,0,0,0))
		self.binhigh=array('d',(60,60,2.07,0.8,0.8,0.8,0.8))
		self.fsparseJet=ROOT.THnSparseD("hsparsejet","hsparsejet;Jet_pt;D_pt;D_m;a10,a15,a20,a30",7,self.nbins,self.binlow,self.binhigh)
		self.fsparseJetvalue=array('d',(0,0,0,0,0,0,0))
		self.fsparseJet.Sumw2()

	
		print(self.config_file)
	
		with open(self.config_file, 'r') as stream:
			self.config = yaml.load(stream,Loader=yaml.FullLoader)	

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
			ja = jet_analysis.JetAnalysis(jet_R = jetR,jet_RecombinationScheme=fj.E_scheme, particle_eta_max=0.9,explicit_ghosts=False)
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
			#if len(j.constituents())<=1:
				#	continue
				
				#if (j.pt())>50:
				#	print("jet pt greater than 50 GeV/c")
				#	continue

				#jets with the winner take all axis################################		
				jet_def_wta = fj.JetDefinition(fj.cambridge_algorithm, 2*jetR)
				jet_def_wta.set_recombination_scheme(fj.WTA_pt_scheme)
				#print('WTA jet definition is:', jet_def_wta)
				reclusterer_wta =  fjcontrib.Recluster(jet_def_wta)
				jet_wta = reclusterer_wta.result(j)
				################################

				self.fsparseJetvalue[0]=j.perp()
				self.fsparseJetvalue[1]=dcand[0].perp()
				self.fsparseJetvalue[2]=dcand[0].m()
				self.fsparseJetvalue[3]=fjext.lambda_beta_kappa(j,  1.0, 1.0 ,0.4)
				self.fsparseJetvalue[4]=fjext.lambda_beta_kappa(j,  1.5, 1.0 ,0.4)
				self.fsparseJetvalue[5]=fjext.lambda_beta_kappa(j,  2.0, 1.0 ,0.4)
				self.fsparseJetvalue[6]=fjext.lambda_beta_kappa(j,  3.0, 1.0 ,0.4)
				#self.fsparseJetvalue[3]=jet_wta.perp()
				#self.fsparseJetvalue[4]=j.delta_R(dcand[0])
				#self.fsparseJetvalue[5]=jet_wta.delta_R(dcand[0])
									
				self.fsparseJet.Fill(self.fsparseJetvalue)
	
			if len(djets) > 1:
				perror("more than one jet per D candidate?")

		return True
	

	def finalize(self):

		self.fsparseJet.Write()	
		self.hNevents.Write()
		self.fout.Close()
		pinfo(self.fout.GetName(), 'written.')


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='D0 analysis on alice data', prog=os.path.basename(__file__))
	parser.add_argument('-c', '--configFile', help='Path of config file for analysis', type=str,metavar='configFile',default='config/configcuts.yaml', required=True)
	parser.add_argument('-f', '--flist', help='file list to process', type=str, default=None, required=True)
	parser.add_argument('-n', '--nfiles', help='max n files to process', type=int, default=0, required=False)
	parser.add_argument('-o', '--output', help="output name / file name in the end", type=str, default='test_hfana')
	args = parser.parse_args()
	
	if not os.path.exists(args.configFile):
		print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
		sys.exit(0)


	hfaio = hfdio.HFAnalysisIO()
	
	hfa = HFAnalysisInvMass(config_file=args.configFile, name = args.output)

	hfa.event_selection.add_selection_range_abs('z_vtx_reco', 10)
	hfa.event_selection.add_selection_equal('is_ev_rej', 0)	
	
	#topomatic cut suggested by D2H
	hfa.d0_selection.add_selection_range_abs('max_norm_d0d0exp',2)

	hfa.d0_selection.add_selection_range('pt_cand', 2, 1e3)
	hfa.d0_selection.add_selection_range_abs('eta_cand', 0.8)
	hfa.d0_selection.add_selection_nsig('nsigTPC_Pi_0', 'nsigTOF_Pi_0', 'nsigTPC_K_1', 'nsigTOF_K_1', 'nsigTPC_Pi_1', 'nsigTOF_Pi_1', 'nsigTPC_K_0', 'nsigTOF_K_0', 3, -900)
		

	hfaio.add_analysis(hfa)
	
	hfaio.execute_analyses_on_file_list(args.flist, args.nfiles)
	
	hfa.finalize()
