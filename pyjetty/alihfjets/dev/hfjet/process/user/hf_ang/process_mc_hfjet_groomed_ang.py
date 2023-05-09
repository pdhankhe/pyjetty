#!/usr/bin/env python3

import numpy as np
import argparse
import os
from array import array
import pandas as pd

from pyjetty.alice_analysis.process.base import process_io, process_utils, jet_info, process_base
import pyjetty.alihfjets.dev.hfjet.process.base.process_io_mc_hf as hfdio
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

	
		self.name='THnSparse_reflection_signal'
		self.nbins=array('i',(55,55,40,40,370))
		self.binlow=array('d',(0,0,0,0,1.7))
		self.binhigh=array('d',(55,55,40,40,2.07))
		title = [ '#it{p}_{T,truth}^{ch jet}', '#it{p}_{T,det}^{ch jet}', '#it{p}_{T,truth}^{D}', '#it{p}_{T,det}^{D}', 'Invmass']
		self.fsparse_reflection=ROOT.THnSparseD(self.name,"Jet_gen_pt;Jet_reco_pt;D_pt_gen;D_pt_rec;D_m",5,self.nbins,self.binlow,self.binhigh)
		self.fsparse_reflection.Sumw2()
		for i in range(0,5):
                        self.fsparse_reflection.GetAxis(i).SetTitle(title[i])
		self.fsparse_reflection_value=array('d',(0,0,0,0,0))

		title_RM = [ '#it{p}_{T,truth}^{ch jet}', '#it{p}_{T,det}^{ch jet}', '#it{p}_{T,truth}^{D}', '#it{p}_{T,det}^{D}', 'Invmass',
				 '#it{#lambda}_{#it{#alpha},truth}', '#it{#lambda}_{#it{#alpha},det}']
		self.nbins_RM=array('i',(55,55,40,40,370,160,160))
		self.binlow_RM=array('d',(0,0,0,0,1.7,0,0))
		self.binhigh_RM=array('d',(55,55,40,40,2.07,0.8,0.8))
        
		self.name_prompt='THnSparse_prompt_signal_alpha_1'
		self.fsparsejet_prompt_alpha1=ROOT.THnSparseD(self.name_prompt,"Jet_gen_pt;Jet_reco_pt;D_pt_gen;D_pt_rec;invmass;dLB1K1_gen;dLB1K1_reco",7,self.nbins_RM,self.binlow_RM,self.binhigh_RM)
		self.fsparsejet_prompt_alpha1.Sumw2()
		self.fsparse_prompt_alpha1_value=array('d',(0,0,0,0,0,0,0))
        
		self.name_prompt='THnSparse_prompt_signal_alpha_15'
		self.fsparsejet_prompt_alpha15=ROOT.THnSparseD(self.name_prompt,"Jet_gen_pt;Jet_reco_pt;D_pt_gen;D_pt_rec;invmass;dLB15K1_gen;dLB15K1_reco",7,self.nbins_RM,self.binlow_RM,self.binhigh_RM)
		self.fsparsejet_prompt_alpha15.Sumw2()
		self.fsparse_prompt_alpha15_value=array('d',(0,0,0,0,0,0,0))

		self.name_prompt='THnSparse_prompt_signal_alpha_2'
		self.fsparsejet_prompt_alpha2=ROOT.THnSparseD(self.name_prompt,"Jet_gen_pt;Jet_reco_pt;D_pt_gen;D_pt_rec;invmass;dLB2K1_gen;dLB2K1_reco",7,self.nbins_RM,self.binlow_RM,self.binhigh_RM)
		self.fsparsejet_prompt_alpha2.Sumw2()
		self.fsparse_prompt_alpha2_value=array('d',(0,0,0,0,0,0,0))
        
        
		self.name_prompt='THnSparse_prompt_signal_alpha_3'
		self.fsparsejet_prompt_alpha3=ROOT.THnSparseD(self.name_prompt,"Jet_gen_pt;Jet_reco_pt;D_pt_gen;D_pt_rec;invmass;dLB3K1_gen;dLB3K1_reco",7,self.nbins_RM,self.binlow_RM,self.binhigh_RM)
		self.fsparsejet_prompt_alpha3.Sumw2()
		self.fsparse_prompt_alpha3_value=array('d',(0,0,0,0,0,0,0))
        
        
		self.name_nonprompt='THnSparse_nonprompt_signal_alpha_1'
		self.fsparsejet_nonprompt_alpha1=ROOT.THnSparseD(self.name_nonprompt,"Jet_gen_pt;Jet_reco_pt;D_pt_gen;D_pt_rec;invmass;dLB1K1_gen;dLB1K1_reco",7,self.nbins_RM,self.binlow_RM,self.binhigh_RM)
		self.fsparsejet_nonprompt_alpha1.Sumw2()
		self.fsparse_nonprompt_alpha1_value=array('d',(0,0,0,0,0,0,0))


		self.name_nonprompt='THnSparse_nonprompt_signal_alpha_15'
		self.fsparsejet_nonprompt_alpha15=ROOT.THnSparseD(self.name_nonprompt,"Jet_gen_pt;Jet_reco_pt;D_pt_gen;D_pt_rec;invmass;dLB15K1_gen;dLB15K1_reco",7,self.nbins_RM,self.binlow_RM,self.binhigh_RM)
		self.fsparsejet_nonprompt_alpha15.Sumw2()
		self.fsparse_nonprompt_alpha15_value=array('d',(0,0,0,0,0,0,0))

		self.name_nonprompt='THnSparse_nonprompt_signal_alpha_2'
		self.fsparsejet_nonprompt_alpha2=ROOT.THnSparseD(self.name_nonprompt,"Jet_gen_pt;Jet_reco_pt;D_pt_gen;D_pt_rec;invmass;dLB2K1_gen;dLB2K1_reco",7,self.nbins_RM,self.binlow_RM,self.binhigh_RM)
		self.fsparsejet_nonprompt_alpha2.Sumw2()
		self.fsparse_nonprompt_alpha2_value=array('d',(0,0,0,0,0,0,0))
		
		self.name_nonprompt='THnSparse_nonprompt_signal_alpha_3'
		self.fsparsejet_nonprompt_alpha3=ROOT.THnSparseD(self.name_nonprompt,"Jet_gen_pt;Jet_reco_pt;D_pt_gen;D_pt_rec;invmass;dLB3K1_gen;dLB3K1_reco",7,self.nbins_RM,self.binlow_RM,self.binhigh_RM)
		self.fsparsejet_nonprompt_alpha3.Sumw2()
		self.fsparse_nonprompt_alpha3_value=array('d',(0,0,0,0,0,0,0))


		for i in range(0,7):
			self.fsparsejet_prompt_alpha1.GetAxis(i).SetTitle(title_RM[i])
			self.fsparsejet_prompt_alpha15.GetAxis(i).SetTitle(title_RM[i])
			self.fsparsejet_prompt_alpha2.GetAxis(i).SetTitle(title_RM[i])
			self.fsparsejet_prompt_alpha3.GetAxis(i).SetTitle(title_RM[i])
			self.fsparsejet_nonprompt_alpha1.GetAxis(i).SetTitle(title_RM[i])
			self.fsparsejet_nonprompt_alpha15.GetAxis(i).SetTitle(title_RM[i])
			self.fsparsejet_nonprompt_alpha2.GetAxis(i).SetTitle(title_RM[i])
			self.fsparsejet_nonprompt_alpha3.GetAxis(i).SetTitle(title_RM[i])


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
			ja = jet_analysis.JetAnalysis(jet_R = jetR,jet_RecombinationScheme=fj.E_scheme, particle_eta_max=0.9,explicit_ghosts=True)	
			ja.analyze_event(_parts_and_ds)

			if len(ja.jets) < 1:
				continue
			jets = ja.jets_as_psj_vector()
			#filtering D0 jets
			djets = djmm.filter_D0_jets(jets)

			if len(djets) > 0:
				j = djets[0]
				dcand = djmm.get_Dcand_in_jet(j)
				################################
				

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
				
				ana_df.at[array_index_df[id0],self.cand_identifier] = df[self.cand_identifier].values[id0]
				ana_df.at[array_index_df[id0],self.jet_identifier] = [jet_groomed.pt(),jet_groomed.eta(),jet_groomed.phi()]	
				ana_df.at[array_index_df[id0],'dLB1K1'] = a10
				ana_df.at[array_index_df[id0],'dLB15K1'] = a15
				ana_df.at[array_index_df[id0],'dLB2K1'] = a20
				ana_df.at[array_index_df[id0],'dLB3K1'] = a30
			
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

		_par_prompt_df = _rec_df[_rec_df['ismcprompt']==1]
		_par_fd_df = _rec_df[_rec_df['ismcfd']==1]
		
		for index_par, row in _par_prompt_df.iterrows():
			self.histo_particle_prompt.Fill(row['jet_pt'],row['pt_cand'])

		for index_par, row in _par_fd_df.iterrows():
			self.histo_particle_fd.Fill(row['jet_pt'],row['pt_cand'])

		self.df_matching=pd.merge(_gen_df,_rec_df,left_index=True, right_index=True)
		
		print("Matching started")
		print(self.df_matching)
		for index_matching, row in self.df_matching.iterrows():
			phi_gen=row['jet_phi_x']
			eta_gen=row['jet_eta_x']
			phi_reco=row['jet_phi_y']
			eta_reco=row['jet_eta_y']
			delta_jet_R=np.sqrt( (eta_gen - eta_reco)**2 + (phi_gen - phi_reco)**2 )
			phi_cand_gen=row['phi_cand_x']
			eta_cand_gen=row['eta_cand_x']
			phi_cand_reco=row['phi_cand_y']
			eta_cand_reco=row['eta_cand_y']
			delta_cand_R=np.sqrt((eta_cand_gen - eta_cand_reco)**2 + (phi_cand_gen - phi_cand_reco)**2 )

			cand_type_reflection = row['ismcrefl_y']
			cand_type_prompt = row['ismcprompt_y']
			cand_type_nonprompt = row['ismcfd_y']
			
			self.fsparse_reflection_value[0]=row['jet_pt_x']
			self.fsparse_reflection_value[1]=row['jet_pt_y']
			self.fsparse_reflection_value[2]=row['pt_cand_x']
			self.fsparse_reflection_value[3]=row['pt_cand_y']
			self.fsparse_reflection_value[4]=row['inv_mass']

			self.fsparse_prompt_alpha1_value[0]=row['jet_pt_x']
			self.fsparse_prompt_alpha1_value[1]=row['jet_pt_y']
			self.fsparse_prompt_alpha1_value[2]=row['pt_cand_x']
			self.fsparse_prompt_alpha1_value[3]=row['pt_cand_y']
			self.fsparse_prompt_alpha1_value[4]=row['inv_mass']
			self.fsparse_prompt_alpha1_value[5]=row['dLB1K1_x']
			self.fsparse_prompt_alpha1_value[6]=row['dLB1K1_y']

			self.fsparse_prompt_alpha15_value[0]=row['jet_pt_x']
			self.fsparse_prompt_alpha15_value[1]=row['jet_pt_y']
			self.fsparse_prompt_alpha15_value[2]=row['pt_cand_x']
			self.fsparse_prompt_alpha15_value[3]=row['pt_cand_y']
			self.fsparse_prompt_alpha15_value[4]=row['inv_mass']
			self.fsparse_prompt_alpha15_value[5]=row['dLB15K1_x']
			self.fsparse_prompt_alpha15_value[6]=row['dLB15K1_y']
            
            
			self.fsparse_prompt_alpha2_value[0]=row['jet_pt_x']
			self.fsparse_prompt_alpha2_value[1]=row['jet_pt_y']
			self.fsparse_prompt_alpha2_value[2]=row['pt_cand_x']
			self.fsparse_prompt_alpha2_value[3]=row['pt_cand_y']
			self.fsparse_prompt_alpha2_value[4]=row['inv_mass']
			self.fsparse_prompt_alpha2_value[5]=row['dLB2K1_x']
			self.fsparse_prompt_alpha2_value[6]=row['dLB2K1_y']
            
            
			self.fsparse_prompt_alpha3_value[0]=row['jet_pt_x']
			self.fsparse_prompt_alpha3_value[1]=row['jet_pt_y']
			self.fsparse_prompt_alpha3_value[2]=row['pt_cand_x']
			self.fsparse_prompt_alpha3_value[3]=row['pt_cand_y']
			self.fsparse_prompt_alpha3_value[4]=row['inv_mass']
			self.fsparse_prompt_alpha3_value[5]=row['dLB3K1_x']
			self.fsparse_prompt_alpha3_value[6]=row['dLB3K1_y']
            
            
			self.fsparse_nonprompt_alpha1_value[0]=row['jet_pt_x']
			self.fsparse_nonprompt_alpha1_value[1]=row['jet_pt_y']
			self.fsparse_nonprompt_alpha1_value[2]=row['pt_cand_x']
			self.fsparse_nonprompt_alpha1_value[3]=row['pt_cand_y']
			self.fsparse_nonprompt_alpha1_value[4]=row['inv_mass']
			self.fsparse_nonprompt_alpha1_value[5]=row['dLB1K1_x']
			self.fsparse_nonprompt_alpha1_value[6]=row['dLB1K1_y']
   
   
			self.fsparse_nonprompt_alpha15_value[0]=row['jet_pt_x']
			self.fsparse_nonprompt_alpha15_value[1]=row['jet_pt_y']
			self.fsparse_nonprompt_alpha15_value[2]=row['pt_cand_x']
			self.fsparse_nonprompt_alpha15_value[3]=row['pt_cand_y']
			self.fsparse_nonprompt_alpha15_value[4]=row['inv_mass']
			self.fsparse_nonprompt_alpha15_value[5]=row['dLB15K1_x']
			self.fsparse_nonprompt_alpha15_value[6]=row['dLB15K1_y']
            
            
            
			self.fsparse_nonprompt_alpha2_value[0]=row['jet_pt_x']
			self.fsparse_nonprompt_alpha2_value[1]=row['jet_pt_y']
			self.fsparse_nonprompt_alpha2_value[2]=row['pt_cand_x']
			self.fsparse_nonprompt_alpha2_value[3]=row['pt_cand_y']
			self.fsparse_nonprompt_alpha2_value[4]=row['inv_mass']
			self.fsparse_nonprompt_alpha2_value[5]=row['dLB2K1_x']
			self.fsparse_nonprompt_alpha2_value[6]=row['dLB2K1_y']
            
            
			self.fsparse_nonprompt_alpha3_value[0]=row['jet_pt_x']
			self.fsparse_nonprompt_alpha3_value[1]=row['jet_pt_y']
			self.fsparse_nonprompt_alpha3_value[2]=row['pt_cand_x']
			self.fsparse_nonprompt_alpha3_value[3]=row['pt_cand_y']
			self.fsparse_nonprompt_alpha3_value[4]=row['inv_mass']
			self.fsparse_nonprompt_alpha3_value[5]=row['dLB3K1_x']
			self.fsparse_nonprompt_alpha3_value[6]=row['dLB3K1_y']
   
   
			#reflection sample
			#non prompt 
			#prompt
			
			if delta_jet_R<(0.6*0.4):
				if delta_cand_R<(0.1):
					#print("matched")
					if cand_type_reflection==1:
						#print("reflection candidate is "+str(row['ismcrefl_y']))
						self.fsparse_reflection.Fill(self.fsparse_reflection_value)

					if cand_type_prompt==1:
						#print("prompt candidate is "+str(row['ismcprompt_y']))
						self.fsparsejet_prompt_alpha1.Fill(self.fsparse_prompt_alpha1_value)
						self.fsparsejet_prompt_alpha15.Fill(self.fsparse_prompt_alpha15_value)
						self.fsparsejet_prompt_alpha2.Fill(self.fsparse_prompt_alpha2_value)
						self.fsparsejet_prompt_alpha3.Fill(self.fsparse_prompt_alpha3_value)
					if cand_type_nonprompt==1:
						#print("nonprompt candidate is "+str(row['ismcfd_y']))
						self.fsparsejet_nonprompt_alpha1.Fill(self.fsparse_nonprompt_alpha1_value)
						self.fsparsejet_nonprompt_alpha15.Fill(self.fsparse_nonprompt_alpha15_value)
						self.fsparsejet_nonprompt_alpha2.Fill(self.fsparse_nonprompt_alpha2_value)
						self.fsparsejet_nonprompt_alpha3.Fill(self.fsparse_nonprompt_alpha3_value)


	def finalize(self):
		self.hNevents.Write()
		self.histo_truth_prompt.Write()
		self.histo_truth_fd.Write()
		self.histo_particle_prompt.Write()
		self.histo_particle_fd.Write()
		self.fsparse_reflection.Write()
		self.fsparsejet_prompt_alpha1.Write()
		self.fsparsejet_prompt_alpha15.Write()
		self.fsparsejet_prompt_alpha2.Write()
		self.fsparsejet_prompt_alpha3.Write()
		self.fsparsejet_nonprompt_alpha1.Write()
		self.fsparsejet_nonprompt_alpha15.Write()
		self.fsparsejet_nonprompt_alpha2.Write()
		self.fsparsejet_nonprompt_alpha3.Write()
		self.fout.Write()
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
