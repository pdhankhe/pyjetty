#!/usr/bin/env python3

import numpy as np
import argparse
import os
from array import array
import pandas as pd

from pyjetty.alice_analysis.process.base import process_io, process_utils, jet_info, process_base
import pyjetty.alihfjets.dev.hfjet.process.user.hf_EEC.process_io_mc_hf as hfdio
from pyjetty.mputils.mputils import perror, pinfo, pwarning
from pyjetty.mputils import treewriter, jet_analysis

import fastjet as fj
import fjext
import fjtools
import fjcontrib
import math as ma
import yaml
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

		self.name='truth_jet_prompt'
		self.histo_truth_prompt = ROOT.TH2F(self.name, self.name,55,0,55,40,0,40)
		self.histo_truth_prompt.Sumw2()
		
		self.name='truth_jet_fd'
		self.histo_truth_fd = ROOT.TH2F(self.name, self.name,55,0,55,40,0,40)
		self.histo_truth_fd.Sumw2()

		self.name='particle_jet_prompt'
		self.histo_particle_prompt = ROOT.TH2F(self.name, self.name,55,0,55,40,0,40)
		self.histo_particle_prompt.Sumw2()

		self.name='particle_jet_fd'
		self.histo_particle_fd = ROOT.TH2F(self.name, self.name,55,0,55,40,0,40)
		self.histo_particle_fd.Sumw2()

		obs_bins_mass = np.concatenate((np.linspace(0, 0.9, 10), np.linspace(1, 9.8, 45), np.linspace(10, 14.5, 10),
		np.linspace(15, 19, 5), np.linspace(20, 60, 9)))

		pt_bins_jet = np.linspace(0, 60, 60)
		pt_bins_dmeson = np.linspace(0, 60, 60)
		dmeson_mass_bin = np.linspace(1.7, 2.07, 370)

		dim = 7
		nbins  = [len(pt_bins_jet)-1, len(pt_bins_dmeson)-1, len(pt_bins_jet)-1, len(pt_bins_dmeson)-1, len(dmeson_mass_bin)-1, len(obs_bins_mass)-1, len(obs_bins_mass)-1]
		min_li = [pt_bins_jet[0],     pt_bins_dmeson[0], pt_bins_jet[0],     pt_bins_dmeson[0],    dmeson_mass_bin[0],     obs_bins_mass[0], obs_bins_mass[0]    ]
		max_li = [pt_bins_jet[-1],    pt_bins_dmeson[-1], pt_bins_jet[-1],    pt_bins_dmeson[-1],    dmeson_mass_bin[-1],    obs_bins_mass[-1], obs_bins_mass[-1]   ]
		nbins = (nbins)
		xmin = (min_li)
		xmax = (max_li)
		nbins_array = array('i', nbins)
		xmin_array = array('d', xmin)
		xmax_array = array('d', xmax)

		print(self.config_file)
	
		with open(self.config_file, 'r') as stream:
			self.config = yaml.load(stream,Loader=yaml.FullLoader)

	def process_events(self, df):
		print(df)
		self.hNevents.Fill(1,df["ev_id"].nunique())
	

	def analysis(self, df, isMC):
		
		isJet=False

		m_array = np.full((self.df_tracks['ParticlePt'].values.size), 0.1396)
		djmm = fjtools.DJetMatchMaker()
		djmm.set_ch_pt_eta_phi_m(self.df_tracks['ParticlePt'].values, self.df_tracks['ParticleEta'].values, self.df_tracks['ParticlePhi'].values, m_array)
		
		if isMC:
			m_cand_gen_array = np.full((df['pt_cand'].values.size), 1.864)
			djmm.set_Ds_pt_eta_phi_m(df['pt_cand'].values, df['eta_cand'].values, df['phi_cand'].values,m_cand_gen_array)

		else:
			djmm.set_Ds_pt_eta_phi_m(df['pt_cand'].values, df['eta_cand'].values, df['phi_cand'].values, df['inv_mass'].values)
			djmm.set_daughters0_pt_eta_phi(df['pt_prong0'].values, df['eta_prong0'].values, df['phi_prong0'].values)
			djmm.set_daughters1_pt_eta_phi(df['pt_prong1'].values, df['eta_prong1'].values, df['phi_prong1'].values)
		
		
		array_index_df= df.index.values
		self.array_col_df=self.cand_identifier+self.jet_identifier
		ana_df=pd.DataFrame(columns=self.array_col_df,index=array_index_df)
		#ana_df=pd.DataFrame(columns=['dR'],index=array_index_df)			
	
		#run for each D candidate to build jet
		for id0, d0 in enumerate(djmm.Ds):
			#daughter tracks matching
			
			if isMC:
				_parts_and_ds=djmm.ch
				for i in range(0,len(_parts_and_ds)):
					for j in range(i+1,len(_parts_and_ds)):
						daughtersSum=_parts_and_ds[i]+_parts_and_ds[j]
						diff=daughtersSum-d0
						if(ma.sqrt((diff.px()*diff.px()))<0.001 and ma.sqrt((diff.py()*diff.py()))<0.001 and ma.sqrt((diff.pz()*diff.pz()))<0.001):
							_parts_and_ds[i]=_parts_and_ds[i]* 1.e-6
							_parts_and_ds[j]=_parts_and_ds[j]* 1.e-6
			else:
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
				isJet=True
				j = djets[0]
				dcand = djmm.get_Dcand_in_jet(j)
				
				# Extract information for EEC   
				constituents = fj.sorted_by_pt(j.constituents())
				c_select = fj.vectorPJ()
				trk_thrd = 1 # track pt threshold  
				# apply pT threshold on jet constituents
				for c in constituents:
					if c.pt() < trk_thrd:
						break
					#print("constituent used for pair =", c)
					c_select.append(c)
				dphi_cut = -9999
				deta_cut = -9999

				################################
				ana_df.at[array_index_df[id0],self.cand_identifier] = df[self.cand_identifier].values[id0]
				ana_df.at[array_index_df[id0],self.jet_identifier] = [j.pt(),j.eta(),j.phi()]	
				ana_df.at[array_index_df[id0],'ungroomed_mass'] = 1
				ana_df.at[array_index_df[id0],'groomed_mass'] = 1
							
			if len(djets) > 1:
				perror("more than one jet per D candidate?")
				continue

		return ana_df

	def fill_generated_info(self,_gen_df):
		print("filling generated information")
		_gen_prompt_df = _gen_df[_gen_df['ismcprompt']==1]
		_gen_fd_df = _gen_df[_gen_df['ismcfd']==1]

		for index_gen, row in _gen_prompt_df.iterrows():
			self.histo_truth_prompt.Fill(row['jet_pt'],row['pt_cand'])
		for index_gen, row in _gen_fd_df.iterrows():
			self.histo_truth_fd.Fill(row['jet_pt'],row['pt_cand'])

	def fill_reco_info(self,_rec_df,_gen_df):
		print("calculate eff")

	def finalize(self):
		self.fout.Close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='D0 analysis on alice data', prog=os.path.basename(__file__))
	parser.add_argument('-c', '--configFile', help='Path of config file for analysis', type=str,metavar='configFile',default='config/configcuts.yaml', required=True)
	parser.add_argument('-f', '--flist', help='file list to process', type=str, default=None, required=True)
	parser.add_argument('-o', '--output', help="output name / file name in the end", type=str, default='test_hfana')
	args = parser.parse_args()

	if not os.path.exists(args.configFile):
		print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
		sys.exit(0)	

	hfaio = hfdio.HFAnalysisIO()
	
	hfa = HFAnalysisInvMass(config_file=args.configFile, name = args.output)

	hfa.event_selection.add_selection_range_abs('z_vtx_reco', 10)
	hfa.event_selection.add_selection_equal('is_ev_rej', 0)	
	
	#topomatic cut suggested by D2H.
	hfa.d0_selection.add_selection_range_abs('max_norm_d0d0exp',2)
	
	hfa.d0_selection.add_selection_range('pt_cand', 2, 1e3)
	hfa.d0_selection.add_selection_range_abs('eta_cand', 0.8)
	hfa.d0_selection.add_selection_nsig('nsigTPC_Pi_0', 'nsigTOF_Pi_0', 'nsigTPC_K_1', 'nsigTOF_K_1', 'nsigTPC_Pi_1', 'nsigTOF_Pi_1', 'nsigTPC_K_0', 'nsigTOF_K_0', 3, -900)

	hfa.d0_gen_selection.add_selection_range('dau_in_acc',0, 1.1)
	hfa.d0_gen_selection.add_selection_range_abs('eta_cand',0.8)
	hfa.d0_gen_selection.add_selection_range('pt_cand',2,1e3) 
	

	hfaio.add_analysis(hfa)
	
	hfaio.execute_analyses_on_file_list(args.flist)
	
	hfa.finalize()
