#!/usr/bin/env python3

import numpy as np
import argparse
import os
from array import array
import pyjetty.alihfjets.dev.hfjet.process_io_data as hfdio
from pyjetty.mputils import perror, pinfo, pwarning, treewriter, jet_analysis

import yaml
import fastjet as fj
import fjext
import fjtools
import fjcontrib
from pyjetty.alice_analysis.process.user.substructure import process_data_base



import ROOT
ROOT.gROOT.SetBatch(True)

class HFAnalysisInvMass(hfdio.HFAnalysis):
	def __init__(self, **kwargs):
		self.fout = None
		super(HFAnalysisInvMass, self).__init__(**kwargs)

		self.fout = ROOT.TFile(self.name+'.root', 'recreate')
		self.fout.cd()
		
		self.nbins=array('i',(60,60,350,160))
		self.binlow=array('d',(0,0,1.72,0))
		self.binhigh=array('d',(60,60,2.07,0.8))

		self.fsparseJet=ROOT.THnSparseD("hsparsejet","hsparsejet;Jet_pt;D_pt;D_m;dR",4,self.nbins,self.binlow,self.binhigh)
		self.fsparseJetvalue=array('d',(0,0,0,0))
		self.fsparseJet.Sumw2()

		self.nbinsWTA=array('i',(60,60,60,350,160,160))
		self.binlowWTA=array('d',(0,0,0,1.72,0,0))
		self.binhighWTA=array('d',(60,60,60,2.07,0.8,0.8))

		self.fsparseJetWTA=ROOT.THnSparseD("hsparsejetWTA","hsparsejetWTA;Jet_pt;JetWTA_pt;D_pt;D_m;dR_WTA;dR_STD_WTA",6,self.nbinsWTA,self.binlowWTA,self.binhighWTA)
		self.fsparseJetWTAvalue=array('d',(0,0,0,0,0,0))
		self.fsparseJetWTA.Sumw2()

		self.nbinsSD=array('i',(60,60,60,350,160,160,160))
		self.binlowSD=array('d',(0,0,0,1.72,-0.1,-0.1,-0.1))
		self.binhighSD=array('d',(60,60,60,2.07,0.7,0.7,0.7))

		self.fsparseJetSD=ROOT.THnSparseD("hsparsejetSD","hsparsejetSD;Jet_pt;JetSD_pt;D_pt;D_m;dR_SD,dR_STD_SD,dR_WTA_SD",7,self.nbinsSD,self.binlowSD,self.binhighSD)
		self.fsparseJetSDvalue=array('d',(0,0,0,0,0,0,0))
		self.fsparseJetSD.Sumw2()
			
		print(self.config_file)
	
		with open(self.config_file, 'r') as stream:
			self.config = yaml.load(stream,Loader=yaml.FullLoader)	
	
	def sd_jet(self,j,jetR,reclustering_algorithm,sd_zcut,sd_beta):
		gshop = fjcontrib.GroomerShop(j, jetR, fj.cambridge_algorithm)
		jet_groomed_lund = gshop.soft_drop(sd_zcut,sd_beta,jetR)
		jet_groomed = jet_groomed_lund.pair()
		deltaR = j.delta_R(jet_groomed)
		if jet_groomed_lund.Delta() < 0:
			deltaR=-1
		return deltaR

	def analysis(self, df):
		#print(df)

		m_array = np.full((self.df_tracks['ParticlePt'].values.size), 0.1396)	

		djmm = fjtools.DJetMatchMaker()
		djmm.set_ch_pt_eta_phi_m(self.df_tracks['ParticlePt'].values, self.df_tracks['ParticleEta'].values, self.df_tracks['ParticlePhi'].values, m_array)
		djmm.set_Ds_pt_eta_phi_m(df['pt_cand'].values, df['eta_cand'].values, df['phi_cand'].values, df['inv_mass'].values)
		djmm.set_daughters0_pt_eta_phi(df['pt_prong0'].values, df['eta_prong0'].values, df['phi_prong0'].values)
		djmm.set_daughters1_pt_eta_phi(df['pt_prong1'].values, df['eta_prong1'].values, df['phi_prong1'].values)

		print(df)
		print(self.df_tracks)
		jetR=0.4
		jets_an = jet_analysis.JetAnalysis(jet_R = jetR,jet_RecombinationScheme=fj.E_scheme, particle_eta_max=0.9,jet_pt_min=5.0)
		_parts_and_ds=djmm.ch
		jets_an.analyze_event(_parts_and_ds)
		if len(jets_an.jets) < 1:
			return False
		jets_selected = jets_an.jets_as_psj_vector()
		print(jets_selected)
		for ja in jets_selected:
			#ja = jets.jets_as_psj_vector()
			print(ja)
			#if len(ja.jets) < 1:
			#	return False


			print("number of constitutent"+str(len(ja.constituents())))
			#if len(ja.constituents())<=1:
			#	return False
	
			constituents = ja.constituents()	
			sorted_by_pt_constituents = fj.sorted_by_pt(constituents)
			lead_part = sorted_by_pt_constituents[0]
			self.fsparseJetvalue[0]=ja.perp()
			self.fsparseJetvalue[1]=lead_part.perp()
			self.fsparseJetvalue[2]=lead_part.m()
			self.fsparseJetvalue[3]=ja.delta_R(lead_part)
			self.fsparseJet.Fill(self.fsparseJetvalue)
	
			#jets with the winner take all axis################################ 
			jet_def_wta = fj.JetDefinition(fj.cambridge_algorithm, 2*jetR)
			jet_def_wta.set_recombination_scheme(fj.WTA_pt_scheme)
			reclusterer_wta =  fjcontrib.Recluster(jet_def_wta)
			jet_wta = reclusterer_wta.result(ja)
			################################

			self.fsparseJetWTAvalue[0]=ja.perp()
			self.fsparseJetWTAvalue[1]=jet_wta.perp()
			self.fsparseJetWTAvalue[2]=lead_part.perp()
			self.fsparseJetWTAvalue[3]=lead_part.m()
			self.fsparseJetWTAvalue[4]=jet_wta.delta_R(lead_part)
			self.fsparseJetWTAvalue[5]=ja.delta_R(jet_wta)

			self.fsparseJetWTA.Fill(self.fsparseJetWTAvalue)
	
			#jets with the SD axis################################ 		
			gshop = fjcontrib.GroomerShop(ja, jetR, fj.cambridge_algorithm)
			sd_zcut=0.1
			sd_beta=0	
			#gshop.soft_drop(beta, zcut, jetR)
			jet_groomed_lund = gshop.soft_drop(sd_beta,sd_zcut,jetR)
			jet_groomed = jet_groomed_lund.pair()
				
			deltaR=jet_groomed.delta_R(lead_part)
			deltaR_STD_SD=ja.delta_R(jet_groomed)
			deltaR_WTA_SD=jet_wta.delta_R(jet_groomed)
		
			print("number of constitutent"+str(len(ja.constituents())))	
			#print("number of constitutent after"+str(len(jet_groomed.constituents())))
			print("after grooming="+str(jet_groomed_lund.Delta()))	
			if jet_groomed_lund.Delta() < 0:
				deltaR=-0.005
				deltaR_STD_SD=-0.005
				deltaR_WTA_SD=-0.005
				print("deltaR_STD_SD="+str(deltaR_STD_SD))
				print("deltaR_WTA_SD="+str(deltaR_WTA_SD))
			################################
			
			self.fsparseJetSDvalue[0]=ja.perp()
			self.fsparseJetSDvalue[1]=jet_groomed.perp()
			self.fsparseJetSDvalue[2]=lead_part.perp()
			self.fsparseJetSDvalue[3]=lead_part.m()
			self.fsparseJetSDvalue[4]=deltaR
			self.fsparseJetSDvalue[5]=deltaR_STD_SD
			self.fsparseJetSDvalue[6]=deltaR_WTA_SD

			self.fsparseJetSD.Fill(self.fsparseJetSDvalue)


		return True
	

	def finalize(self):
		self.fsparseJet.Write()
		self.fsparseJetWTA.Write()
		self.fsparseJetSD.Write()
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
