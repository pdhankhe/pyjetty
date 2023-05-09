#!/usr/bin/env python3

import numpy as np
import argparse
import os
from array import array
import pandas as pd
import pyjetty.alihfjets.dev.hfjet.MC_process_io_mc_hf as hfdio
from pyjetty.alice_analysis.process.base import process_io, process_utils, jet_info, process_base
from pyjetty.mputils import perror, pinfo, pwarning, treewriter, jet_analysis


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
		self.twjc = treewriter.RTreeWriter(tree_name='d0jc', fout=self.fout)
		self.twjc_gen = treewriter.RTreeWriter(tree_name='d0jc_gen', fout=self.fout)

		self.name='truth_jetpt_D_pt_cand'
		self.histo_truth = ROOT.TH3F(self.name, self.name,55,0,55,40,0,40,2,0,2)
		self.histo_truth.Sumw2()
		self.fhistovalue=array('d',(0,0,0))

		self.name2='THnsparese_reco_level_jetpt'
		self.nbins=array('i',(55,55,40,35,100,220,50,10000,2,2,60,60))
		self.binlow=array('d',(0,0,0,1.7,0,-1.1,0,-0.02,0,0,0,0))
		self.binhigh=array('d',(55,55,40,2.05,1,1.1,50,0.02,2,2,0.6,0.6))
		self.fsparseJet=ROOT.THnSparseD("hsparsejet","hsparsejet;Jet_gen_pt;Jet_reco_pt;D_pt;D_m;D_cosp;D_cos_t_star;norm_d_lxy;imp_par_pro;cand_type_gen;cand_type_reco;dR_gen;dR_reco;",12,self.nbins,self.binlow,self.binhigh)
		self.fsparseJet.Sumw2()
		self.fsparseJetvalue=array('d',(0,0,0,0,0,0,0,0,0,0,0,0))

		print(self.config_file)
	
		with open(self.config_file, 'r') as stream:
			self.config = yaml.load(stream,Loader=yaml.FullLoader)
	
	def analysis(self, df):
		#print("number of candidate per event"+ str(len(df)))
		#print(df)
		m_array = np.full((self.df_tracks['ParticlePt'].values.size), 0.1396)	

		djmm = fjtools.DJetMatchMaker()
		djmm.set_ch_pt_eta_phi_m(self.df_tracks['ParticlePt'].values, self.df_tracks['ParticleEta'].values, self.df_tracks['ParticlePhi'].values, m_array)
		djmm.set_Ds_pt_eta_phi_m(df['pt_cand'].values, df['eta_cand'].values, df['phi_cand'].values, df['inv_mass'].values)
		djmm.set_daughters0_pt_eta_phi(df['pt_prong0'].values, df['eta_prong0'].values, df['phi_prong0'].values)
		djmm.set_daughters1_pt_eta_phi(df['pt_prong1'].values, df['eta_prong1'].values, df['phi_prong1'].values)
		
		array_rec_index_df= df.index.values
		array_rec_col_df=[]
	
		#print(*djmm.Ds)
		pd_reco=pd.DataFrame(columns=array_rec_col_df,index=array_rec_index_df)
		#run for each D candidate to build jet
		for id0, d0 in enumerate(djmm.Ds):
			#daughter tracks matching
			_parts_and_ds = djmm.match(0.005, id0)
			#replacing daughter tracks with matched D0 candidate
			#including D0 	
			_parts_and_ds.push_back(d0)
			#print("daughter 1")
			#print(*djmm.daughters0)
			daughter1=djmm.daughters0[id0]	

			#print("daughter 2")
			#print(*djmm.daughters1)
			daughter2=djmm.daughters1[id0]
			dr_daughters=daughter1.delta_R(daughter2)
			print("angular distance between daughters: "+str(dr_daughters))
			print("number of particle for jet reconstrcution"+str(len(_parts_and_ds)))	
			#jet reconstruction with D0 and charged particle
			jetR=0.4
			ja = jet_analysis.JetAnalysis(jet_R = jetR,jet_RecombinationScheme=fj.pt_scheme, particle_eta_max=0.9,explicit_ghosts=False)
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

				#constituents=j.constituents()
				#print(*constituents)	
				#print(constituents[0].pt())
				#highest_pt_constituent=fj.sorted_by_pt(constituents)
				#print("number of constitutents"+str(len(j.constituents())))
				#print(*j.constituents())
				#print(*highest_pt_constituent)
				#jets with the winner take all axis################################		
				jet_def_wta = fj.JetDefinition(fj.cambridge_algorithm, 2*jetR)
				jet_def_wta.set_recombination_scheme(fj.WTA_pt_scheme)
				#print('WTA jet definition is:', jet_def_wta)
				reclusterer_wta =  fjcontrib.Recluster(jet_def_wta)
				jet_wta = reclusterer_wta.result(j)
				################################
				D_cand_type = df['cand_type'].values
				
				
				pd_reco.at[array_rec_index_df[id0],'run_number']=df['run_number'].values[id0]
				pd_reco.at[array_rec_index_df[id0],'ev_id']=df['ev_id'].values[id0]
				pd_reco.at[array_rec_index_df[id0],'ev_id_ext']=df['ev_id_ext'].values[id0]
				pd_reco.at[array_rec_index_df[id0],'cand_type']=df['cand_type'].values[id0]
				pd_reco.at[array_rec_index_df[id0],'pt_cand']=df['pt_cand'].values[id0]
				pd_reco.at[array_rec_index_df[id0],'eta_cand']=df['eta_cand'].values[id0]
				pd_reco.at[array_rec_index_df[id0],'phi_cand']=df['phi_cand'].values[id0]
				pd_reco.at[array_rec_index_df[id0],'inv_mass']=df['inv_mass'].values[id0]
				pd_reco.at[array_rec_index_df[id0],'jet_pt']=j.pt()
				pd_reco.at[array_rec_index_df[id0],'jet_eta']=j.eta()
				pd_reco.at[array_rec_index_df[id0],'jet_phi']=j.phi()
				pd_reco.at[array_rec_index_df[id0],'dR']=j.delta_R(dcand[0])

			if len(djets) > 1:
				perror("more than one jet per D candidate?")
				continue

		#print(pd_reco)
		return pd_reco

	def analysis_gen(self, df):
		m_gen_array = np.full((self.df_gen_tracks['ParticlePt'].values.size), 0.1396)
		m_cand_gen_array = np.full((df['pt_cand'].values.size), 1.864)
		djmm_gen = fjtools.DJetMatchMaker()
		djmm_gen.set_ch_pt_eta_phi_m(self.df_gen_tracks['ParticlePt'].values, self.df_gen_tracks['ParticleEta'].values, self.df_gen_tracks['ParticlePhi'].values, m_gen_array)
		djmm_gen.set_Ds_pt_eta_phi_m(df['pt_cand'].values, df['eta_cand'].values, df['phi_cand'].values,m_cand_gen_array)


		array_gen_index_df= df.index.values
		array_gen_col_df=[]
		pd_gen=pd.DataFrame(columns=array_gen_col_df,index=array_gen_index_df)
		for id0, d0 in enumerate(djmm_gen.Ds):
			#_parts_and_ds = djmm_gen.ch
			_parts_and_ds = djmm_gen.ch
			#print("number of particles for jet reco: "+str(len(_parts_and_ds)))
			print(*_parts_and_ds)
			for i in range(0,len(_parts_and_ds)):
				for j in range(i+1,len(_parts_and_ds)):
						daughtersSum=_parts_and_ds[i]+_parts_and_ds[j]
						diff=daughtersSum-d0				
						if(ma.sqrt((diff.px()*diff.px()))<0.001 and ma.sqrt((diff.py()*diff.py()))<0.001 and ma.sqrt((diff.pz()*diff.pz()))<0.001):
							print("match")
							print("daughter 1: "+str(_parts_and_ds[i]))
							print("daugter 2: "+str(_parts_and_ds[j]))
							_parts_and_ds[i]=_parts_and_ds[i]* 1.e-6
							_parts_and_ds[j]=_parts_and_ds[j]* 1.e-6

			
			print("number of particles for jet reco after replacing daughters: "+str(len(_parts_and_ds)))
			_parts_and_ds.push_back(d0)
			print(*_parts_and_ds)
			ja_gen = jet_analysis.JetAnalysis(jet_R = 0.4,jet_RecombinationScheme=fj.pt_scheme, particle_eta_max=0.9,explicit_ghosts=False)

			ja_gen.analyze_event(_parts_and_ds)
			if len(ja_gen.jets) < 1:
				continue

			jets = ja_gen.jets_as_psj_vector()
			djets = djmm_gen.filter_D0_jets(jets)

			if len(djets) > 0:
				j_gen = djets[0]
				dcand = djmm_gen.get_Dcand_in_jet(j_gen)
				#print("D candidate selected")
				#print(dcand[0])	
				#D_cand_type = df['cand_type'].values
				#self.twjc_gen.fill_branches(jet = j_gen,dR = j_gen.delta_R(dcand[0]),D = dcand[0],Dmeson_cand_type=float(D_cand_type[id0]))
				#self.twjc_gen.fill_tree()
				#constituents=j_gen.constituents()
				#highest_pt_constituent=fj.sorted_by_pt(constituents)
				#print(*highest_pt_constituent)
				pd_gen.at[array_gen_index_df[id0],'run_number']=df['run_number'].values[id0]
				pd_gen.at[array_gen_index_df[id0],'ev_id']=df['ev_id'].values[id0]
				pd_gen.at[array_gen_index_df[id0],'ev_id_ext']=df['ev_id_ext'].values[id0]
				pd_gen.at[array_gen_index_df[id0],'cand_type']=df['cand_type'].values[id0]
				pd_gen.at[array_gen_index_df[id0],'pt_cand']=df['pt_cand'].values[id0]
				pd_gen.at[array_gen_index_df[id0],'eta_cand']=df['eta_cand'].values[id0]
				pd_gen.at[array_gen_index_df[id0],'phi_cand']=df['phi_cand'].values[id0]
				pd_gen.at[array_gen_index_df[id0],'jet_pt']=j_gen.pt()
				pd_gen.at[array_gen_index_df[id0],'jet_eta']=j_gen.eta()
				pd_gen.at[array_gen_index_df[id0],'jet_phi']=j_gen.phi()
				pd_gen.at[array_gen_index_df[id0],'dR']=j_gen.delta_R(dcand[0])
				pd_gen.at[array_gen_index_df[id0],'dau_in_acc']=df['dau_in_acc'].values[id0]

			if len(djets) > 1:
				perror("more than one jet per D candidate?")
				continue
		print(pd_gen)
		return pd_gen

	def ResponseMatrix(self,df_fjparticles_reco,df_fjparticles_gen):

		print("response matrix generated here")
		array_rec_index_df=df_fjparticles_reco.index.values
		array_gen_index_df=df_fjparticles_gen.index.values
	
		for i, index_gen in enumerate(array_gen_index_df):
			#print(index_gen)
			Jetpt_gen=df_fjparticles_gen['jet_pt'].values[i]
			Dpt_gen=df_fjparticles_gen['pt_cand'].values[i]			
			CandType_gen=df_fjparticles_gen['cand_type'].values[i]
			pr = int(CandType_gen) & 0b1000
			fdn = int(CandType_gen) & 0b10000
			
			if pr:
				cand_type=1.5
			
			if fdn:
				cand_type=0.5
			
			if pr or fdn:
				self.histo_truth.Fill(Jetpt_gen,Dpt_gen,cand_type)
			
			for j, index_reco in enumerate(array_rec_index_df):
				if index_gen==index_reco:
					#print("index matched run jet matching")
					#print("gen index"+str(index_gen)+": reco index"+str(index_reco))
	
					phi_gen=df_fjparticles_gen['jet_phi'].values[i]
					eta_gen=df_fjparticles_gen['jet_eta'].values[i]
					phi_reco=df_fjparticles_reco['jet_phi'].values[j]
					eta_reco=df_fjparticles_reco['jet_eta'].values[j]
					phi_cand_gen=df_fjparticles_gen['phi_cand'].values[i]
					eta_cand_gen=df_fjparticles_gen['eta_cand'].values[i]
					phi_cand_reco=df_fjparticles_reco['phi_cand'].values[j]
					eta_cand_reco=df_fjparticles_reco['eta_cand'].values[j]					

					self.fsparseJetvalue[0]=Jetpt_gen
					self.fsparseJetvalue[1]=df_fjparticles_reco['jet_pt'].values[j]
					self.fsparseJetvalue[2]=df_fjparticles_reco['pt_cand'].values[j]
					self.fsparseJetvalue[3]=df_fjparticles_reco['inv_mass'].values[j]
					self.fsparseJetvalue[8]=cand_type
					Cand_type_reco=df_fjparticles_reco['cand_type'].values[j]
					self.fsparseJetvalue[10]=df_fjparticles_gen['dR'].values[i]
					self.fsparseJetvalue[11]=df_fjparticles_reco['dR'].values[j]
					delta_R=np.sqrt( (eta_gen - eta_reco)**2 + (phi_gen - phi_reco)**2 )
					delta_cand_R=np.sqrt((eta_cand_gen - eta_cand_reco)**2 + (phi_cand_gen - phi_cand_reco)**2 )

					pr_reco = int(Cand_type_reco) & 0b1000
					fdn_reco= int(Cand_type_reco) & 0b10000				
					
					if pr_reco:
						self.fsparseJetvalue[9]=1.5
					if fdn_reco:
						self.fsparseJetvalue[9]=0.5

					if pr_reco or fdn_reco:
						if delta_R<(0.6*0.4):
							if delta_cand_R<(0.1):
								print("matched")
								print("gen index"+str(index_gen)+": reco index"+str(index_reco))
								#print(cand_type)
								print("gen candidate type: "+str(cand_type)+"reco candidate type: "+str(self.fsparseJetvalue[9]))
								self.fsparseJet.Fill(self.fsparseJetvalue)
							#print("gen index:"+str(index_gen)+" recon index:"+str(index_gen)+"cand type gen"+str(Jetpt_gen)+"can type reco"+str(self.fsparseJetvalue[1])+"cand type:"+str(self.fsparseJetvalue[9])+"gen can eta phi:"+ str(eta_cand_gen)+","+str(phi_cand_gen)+"reco can eta phi:"+str(eta_cand_reco)+","+str(phi_cand_reco)+"Dcan pt gen: " +str(Dpt_gen)+"Dcan pt reco: "+str(self.fsparseJetvalue[2]))	
								#print("generated jet info")
							#print("cand_type  ; pt_cand;  eta_cand ;phi_cand   ; jet_pt  ; jet_eta ; jet_phi ; jet constitutent N;  highest pt inside jet ")
								#print(df_fjparticles_gen.values[i])
								#print("reconstructed jet info")
							#print("cand_type  ; pt_cand;  eta_cand ;phi_cand   ; jet_pt  ; jet_eta ; jet_phi ; jet constitutent N;  highest pt inside jet ")
								#print(df_fjparticles_reco.values[j])
			
	def finalize(self):
	
		self.histo_truth.Write()
		self.fsparseJet.Write()
		self.fout.Write()
		self.fout.Close()

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
	
	#topomatic cut suggested by D2H.
	hfa.d0_selection.add_selection_range_abs('max_norm_d0d0exp',2)

	hfa.d0_selection.add_selection_range('pt_cand', 2, 1e3)
	hfa.d0_selection.add_selection_range_abs('eta_cand', 0.8)
	hfa.d0_selection.add_selection_nsig('nsigTPC_Pi_0', 'nsigTOF_Pi_0', 'nsigTPC_K_1', 'nsigTOF_K_1', 'nsigTPC_Pi_1', 'nsigTOF_Pi_1', 'nsigTPC_K_0', 'nsigTOF_K_0', 3, -900)


	hfa.d0_gen_selection.add_selection_range('dau_in_acc',0, 1.1)
	hfa.d0_gen_selection.add_selection_range_abs('eta_cand',0.8)
	hfa.d0_gen_selection.add_selection_range('pt_cand',2,1e3) 
	

	hfaio.add_analysis(hfa)
	
	hfaio.execute_analyses_on_file_list(args.flist, args.nfiles)
	
	hfa.finalize()
