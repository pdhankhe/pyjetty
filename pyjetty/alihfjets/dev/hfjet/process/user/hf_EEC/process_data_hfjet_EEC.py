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
import ecorrel

import ROOT
ROOT.gROOT.SetBatch(True)

def linbins(xmin, xmax, nbins):
	lspace = np.linspace(xmin, xmax, nbins+1)
	arr = array('f', lspace)
	return arr

def logbins(xmin, xmax, nbins):
	lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
	arr = array('f', lspace)
	return arr


class HFAnalysisInvMass(hfdio.HFAnalysis):
	def __init__(self, **kwargs):
		self.fout = None
		super(HFAnalysisInvMass, self).__init__(**kwargs)

		self.fout = ROOT.TFile(self.name+'.root', 'recreate')
		self.fout.cd()
		
		self.hNevents = ROOT.TH1F("hNevents", "hNevents", 2, array('f', [-0.5, 0.5, 1.5]) )
		self.hNevents.Sumw2()

		title = [ '#it{p}_{T}^{ch jet}', '#it{p}_{T}^{D^{0}}', 'Invmass', '#it{R}_{L}']
		title_jetinfo = [ '#it{p}_{T}^{ch jet}', '#it{p}_{T}^{D^{0}}', 'Invmass', '#it{N}_{const}']
		
		RL_bins = np.linspace(0, 1, 500) 

		pt_bins_jet = np.linspace(0, 60, 60)
		pt_bins_dmeson = np.linspace(0, 60, 60)
		dmeson_mass_bin = np.linspace(1.7, 2.07, 370)

		dim = 4
		nbins  = [len(pt_bins_jet)-1, len(pt_bins_dmeson)-1, len(dmeson_mass_bin)-1, len(RL_bins)-1]
		min_li = [pt_bins_jet[0],     pt_bins_dmeson[0],     dmeson_mass_bin[0],     RL_bins[0]    ]
		max_li = [pt_bins_jet[-1],    pt_bins_dmeson[-1],    dmeson_mass_bin[-1],    RL_bins[-1]   ]

		nbins = (nbins)
		xmin = (min_li)
		xmax = (max_li)

		nbins_array = array('i', nbins)
		xmin_array = array('d', xmin)
		xmax_array = array('d', xmax)

		self.fsparseJet = ROOT.THnSparseD("hsparsejet","hsparsejet; Jet_pt; D_pt; D_m; RL", dim,  nbins_array, xmin_array, xmax_array)
		self.fsparseJetvalue = array( 'd', ( 0, 0, 0, 0))
		self.fsparseJet.Sumw2()
		
		self.fsparseJetInfo = ROOT.THnSparseD("hsparsejet_info","hsparsejet_info; Jet_pt; D_pt; D_m; Nconst.", dim, nbins_array, xmin_array, xmax_array)
		self.fsparseJetInfovalue = array( 'd', ( 0, 0, 0, 0))
		self.fsparseJetInfo.Sumw2()

		for i in range(0,dim):
			self.fsparseJet.GetAxis(i).SetTitle(title[i])
			self.fsparseJetInfo.GetAxis(i).SetTitle(title_jetinfo[i])

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

				self.fsparseJetInfovalue[0] = j.perp()
				self.fsparseJetInfovalue[1] = dcand[0].perp()
				self.fsparseJetInfovalue[2] = dcand[0].m()
				self.fsparseJetInfovalue[3] = len(c_select)

				self.fsparseJetInfo.Fill(self.fsparseJetInfovalue)
		
				new_corr = ecorrel.CorrelatorBuilder(c_select, dcand, j.perp(), 2, 1, dphi_cut, deta_cut)
				
				for index in range(new_corr.correlator(2).rs().size()):
					
					self.fsparseJetvalue[0] = j.perp()
					self.fsparseJetvalue[1] = dcand[0].perp()
					self.fsparseJetvalue[2] = dcand[0].m()
					self.fsparseJetvalue[3] = new_corr.correlator(2).rs()[index]
					#fill thnsparse with per pair and weight 				
					self.fsparseJet.Fill(self.fsparseJetvalue, new_corr.correlator(2).weights()[index])
	
			if len(djets) > 1:
				perror("more than one jet per D candidate?")

		return True
	

	def finalize(self):
		self.hNevents.Write()
		self.fsparseJet.Write()	
		self.fsparseJetInfo.Write()
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
