#!/usr/bin/env python3

import numpy as np
import argparse
import os
from array import array
import pandas as pd

from pyjetty.alice_analysis.process.base import process_io, process_utils, jet_info, process_base
#import pyjetty.alihfjets.dev.hfjet.DmesonJet.analysis.EEC.testcode.process_io_mc_hf as hfdio
import pyjetty.alihfjets.dev.hfjet.process.base.process_io_mc_hf as hfdio
from pyjetty.mputils.mputils import perror, pinfo, pwarning
from pyjetty.mputils import treewriter, jet_analysis

import fastjet as fj
import fjext
import fjtools
import tqdm
import fjcontrib
import math as ma
import yaml
import ROOT
import time
import ecorrel
ROOT.gROOT.SetBatch(True)

class HFAnalysisInvMass(hfdio.HFAnalysis, process_base.ProcessBase):

    def __init__(self, **kwargs):
        super(HFAnalysisInvMass, self).__init__(**kwargs)
        #self.fout = ROOT.TFile(self.name+'.root', 'recreate')
        #self.fout.cd()
        self.initialize_config()
        print(self.config_file)


    def HFAnalysis(self):
            print("dataframe will be created jet analysis")
            start_time = time.time()
            print('--- {} seconds ---'.format(time.time() - start_time))
            
            # Initialize configuration and histograms
            self.initialize_config()
            self.initialize_user_output_objects()
            
            # Find jets and fill histograms
            print('Find jets...')
            self.analyzeEvents()
            
            
            # Plot histograms
            print('Save histograms...')
            print('--- {} seconds ---'.format(time.time() - start_time))
            process_base.ProcessBase.save_output_objects(self)
      
    def initialize_config(self):
        print("initializing the configuration")
        with open(self.config_file, 'r') as stream:
            self.config = yaml.load(stream,Loader=yaml.FullLoader)
    
        
        # Set configuration for analysis
        self.jetR_list = self.config["jetR"]
        # Mass assumption for jet reconstruction
        self.track_mass = self.config["track_mass"]
        self.track_random_mass = self.config["track_random_mass"]   # randomly assign K or p mass to tracks
        # Detector kinematics
        self.eta_max = 0.9
    
    #---------------------------------------------------------------
    # Initialize histograms
    #---------------------------------------------------------------
    def initialize_user_output_objects(self):
        title = [ '#it{p}_{T}^{ch jet}', '#it{p}_{T}^{D^{0}}', 'Invmass', 'Candidate Info', '#it{R}_{L}']
        title_jetinfo = [ '#it{p}_{T}^{ch jet}', '#it{p}_{T}^{D^{0}}', 'Invmass', '#it{N}_{const}']
        RL_bins  = np.logspace(np.log10(1E-4), np.log10(1), 51)
        print("print RL bins")
        #print(RL_bins)
        pt_bins_jet = np.linspace(0, 60, 61)
        pt_bins_dmeson = np.linspace(0, 60, 61)
        dmeson_mass_bin = np.linspace(1.7, 2.07, 371)
        cand_type_bins = np.linspace(0,3,4)

        dim = 5 #set this to 5??
        nbins  = [len(pt_bins_jet)-1, len(pt_bins_dmeson)-1, len(dmeson_mass_bin)-1, 3, 50] #should nbins be 3, not 2, for cand_type_bins??
        min_li = [pt_bins_jet[0],     pt_bins_dmeson[0],     dmeson_mass_bin[0],     cand_type_bins[0],     RL_bins[0]     ]
        max_li = [pt_bins_jet[-1],    pt_bins_dmeson[-1],    dmeson_mass_bin[-1],    cand_type_bins[-1],    RL_bins[-1]    ]

        nbins = (nbins)
        xmin = (min_li)
        xmax = (max_li)
        
        nbins_array = array('i', nbins)
        xmin_array = array('d', xmin)
        xmax_array = array('d', xmax)

        self.fsparseJet_reco = ROOT.THnSparseD("hsparsejet_reco","hsparsejet_reco; Jet_pt; D_pt; D_m; cand_type; RL", dim,  nbins_array, xmin_array, xmax_array)
        self.fsparseJet_gen = ROOT.THnSparseD("hsparsejet_gen","hsparsejet_gen; Jet_pt; D_pt; D_m; cand_type; RL", dim,  nbins_array, xmin_array, xmax_array)
        self.fsparseJetCandInfo_reco = ROOT.THnSparseD("hsparsejetcand_reco","hsparsejetcand_reco; Jet_pt; D_pt; D_m; cand_type", dim-1,  nbins_array[:-1], xmin_array[:-1], xmax_array[:-1])
        self.fsparseJetCandInfo_gen = ROOT.THnSparseD("hsparsejetcand_gen","hsparsejetcand_gen; Jet_pt; D_pt; D_m; cand_type", dim-1,  nbins_array[:-1], xmin_array[:-1], xmax_array[:-1])
        
        self.fsparseJetvalue = array( 'd', ( 0, 0, 0, 0, 0)) # added one more zero here
        self.fsparseJetCandInfo = array( 'd', ( 0, 0, 0, 0))
        
        self.fsparseJet_reco.Sumw2()
        self.fsparseJet_gen.Sumw2()
        self.fsparseJetCandInfo_reco.Sumw2()
        self.fsparseJetCandInfo_gen.Sumw2()
        
        for i in range(0,dim):
            self.fsparseJet_reco.GetAxis(i).SetTitle(title[i])
            self.fsparseJet_gen.GetAxis(i).SetTitle(title[i])
            if i < dim-1:
                self.fsparseJetCandInfo_reco.GetAxis(i).SetTitle(title[i])
                self.fsparseJetCandInfo_gen.GetAxis(i).SetTitle(title[i])
            
            if i == 0 or i == 1:
                self.fsparseJet_reco.SetBinEdges(i, pt_bins_jet)
                self.fsparseJet_gen.SetBinEdges(i, pt_bins_jet)
                self.fsparseJetCandInfo_reco.SetBinEdges(i, pt_bins_jet)
                self.fsparseJetCandInfo_gen.SetBinEdges(i, pt_bins_jet)
            if i == 2:
                self.fsparseJet_reco.SetBinEdges(i, dmeson_mass_bin)
                self.fsparseJet_gen.SetBinEdges(i, dmeson_mass_bin)
                self.fsparseJetCandInfo_reco.SetBinEdges(i, dmeson_mass_bin)
                self.fsparseJetCandInfo_gen.SetBinEdges(i, dmeson_mass_bin)
            if i == 3:
                self.fsparseJet_reco.SetBinEdges(i, cand_type_bins)
                self.fsparseJet_gen.SetBinEdges(i, cand_type_bins)
                self.fsparseJetCandInfo_reco.SetBinEdges(i, cand_type_bins)
                self.fsparseJetCandInfo_gen.SetBinEdges(i, cand_type_bins)
            if i == 4:
                self.fsparseJet_reco.SetBinEdges(i, RL_bins)
                self.fsparseJet_gen.SetBinEdges(i, RL_bins)

                
    def analyzeEvents(self):
        self.hfio.execute_analyses_on_file_list(self.input_file)
    
        
    def D_jet_reconstruction_event_by_event(self, df, IsGen):
    
        m_array = np.full((self.df_tracks['ParticlePt'].values.size), 0.1396)
        djmm = fjtools.DJetMatchMaker()
        djmm.set_ch_pt_eta_phi_m(self.df_tracks['ParticlePt'].values, self.df_tracks['ParticleEta'].values, self.df_tracks['ParticlePhi'].values, m_array)
        
        if IsGen:
            m_cand_gen_array = np.full((df['pt_cand'].values.size), 1.864)
            djmm.set_Ds_pt_eta_phi_m(df['pt_cand'].values, df['eta_cand'].values, df['phi_cand'].values,m_cand_gen_array)
        else:
            djmm.set_Ds_pt_eta_phi_m(df['pt_cand'].values, df['eta_cand'].values, df['phi_cand'].values, df['inv_mass'].values)
            djmm.set_daughters0_pt_eta_phi(df['pt_prong0'].values, df['eta_prong0'].values, df['phi_prong0'].values)
            djmm.set_daughters1_pt_eta_phi(df['pt_prong1'].values, df['eta_prong1'].values, df['phi_prong1'].values)
        
        for id0, d0 in enumerate(djmm.Ds):
            #replacing daughter tracks with matched D0 candidate
            if IsGen:
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
                
            #including D0
            _parts_and_ds.push_back(d0)
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
                print(j)
                print(dcand)
                
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
                cand_type_prompt = df['ismcprompt'].values[id0]
                cand_type_non_prompt = df['ismcfd'].values[id0]
                cand_type = -1
                
                if cand_type_prompt == 1: #prompt
                    cand_type = 1
                elif cand_type_non_prompt == 1: #non prompt
                    cand_type = 2
                else: #reflection
                    cand_type = 0
                
                self.fsparseJetCandInfo[0] = j.perp()
                self.fsparseJetCandInfo[1] = dcand[0].perp()
                self.fsparseJetCandInfo[2] = dcand[0].m()
                self.fsparseJetCandInfo[3] = cand_type
                
                if IsGen:
                    self.fsparseJetCandInfo_gen.Fill(self.fsparseJetCandInfo)
                else:
                    self.fsparseJetCandInfo_reco.Fill(self.fsparseJetCandInfo)
                
                new_corr = ecorrel.CorrelatorBuilder(c_select, dcand, j.perp(), 2, 1, dphi_cut, deta_cut)
#                print("df[ismcprompt]", df['ismcprompt'])
#                print("df[ismcprompt][0]", df['ismcprompt'].values[0])
                
                
#                if cand_type_prompt==1: # filling information only for prompt candidate for now
                for index in range(new_corr.correlator(2).rs().size()): # loop over pairs
                    self.fsparseJetvalue[0] = j.perp()
                    self.fsparseJetvalue[1] = dcand[0].perp()
                    self.fsparseJetvalue[2] = dcand[0].m()
                    self.fsparseJetvalue[3] = cand_type
                    self.fsparseJetvalue[4] = new_corr.correlator(2).rs()[index]
                    #fill thnsparse with per pair and weight
                    if IsGen:
                        self.fsparseJet_gen.Fill(self.fsparseJetvalue, new_corr.correlator(2).weights()[index])
                    else:
                        self.fsparseJet_reco.Fill(self.fsparseJetvalue, new_corr.correlator(2).weights()[index])
                    
            if len(djets) > 1:
                perror("more than one jet per D candidate?")
                continue
      
   

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='D0 analysis on alice data', prog=os.path.basename(__file__))
    parser.add_argument('-c', '--configFile', help='Path of config file for analysis', type=str,metavar='configFile',default='config/configcuts.yaml', required=True)
    parser.add_argument('-f', '--inputFile', action='store',
                      type=str, metavar='inputFile',
                      default='AnalysisResults.root',
                      help='Path of ROOT file containing TTrees')
    parser.add_argument('-o', '--outputDir', action='store',
                      type=str, metavar='outputDir',
                      default='./TestOutput',
                      help='Output directory for QA plots to be written to')
    parser.add_argument('-ct', '--configType', action='store',
                      type=str, metavar='configType',
                      default='default',
                      help='Configuration type of cuts to be chosen from config file')
    args = parser.parse_args()
    
    print('Configuring...')
    print('inputFile: \'{0}\''.format(args.inputFile))
    print('configFile: \'{0}\''.format(args.configFile))
    print('ouputDir: \'{0}\''.format(args.outputDir))
    print('configType: \'{0}\''.format(args.configType))
    print('----------------------------------------------------------------')
  
    if not os.path.exists(args.configFile):
        print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)

    hfaio = hfdio.HFAnalysisIO()

    hfa = HFAnalysisInvMass(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir, config_type=args.configType, hfio=hfaio)

    print("Set event selection cuts here")
    hfa.event_selection.add_selection_range_abs('z_vtx_reco', 10)
    hfa.event_selection.add_selection_equal('is_ev_rej', 0)

    print("Set D0 pt independent selection cuts here")
    #topomatic cut suggested by D2H.
    hfa.d0_selection.add_selection_range_abs('max_norm_d0d0exp',2)
    hfa.d0_selection.add_selection_range('pt_cand', 2, 1e3)
    hfa.d0_selection.add_selection_range_abs('eta_cand', 0.8)
    hfa.d0_selection.add_selection_nsig('nsigTPC_Pi_0', 'nsigTOF_Pi_0', 'nsigTPC_K_1', 'nsigTOF_K_1', 'nsigTPC_Pi_1', 'nsigTOF_Pi_1', 'nsigTPC_K_0', 'nsigTOF_K_0', 3, -900)

    hfa.d0_gen_selection.add_selection_range('dau_in_acc',0, 1.1)
    hfa.d0_gen_selection.add_selection_range_abs('eta_cand',0.8)
    hfa.d0_gen_selection.add_selection_range('pt_cand',2,1e3)


    hfaio.add_analysis(hfa)
    print("Start reading input output and perform analysis")
    hfa.HFAnalysis()
  

