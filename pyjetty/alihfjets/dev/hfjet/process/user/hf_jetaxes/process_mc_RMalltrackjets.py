#!/usr/bin/env python3

'''
This class process the MC data and
performs the Djet reconstruction generated and reconstructed level
perform jet matching and fill the histograms for histograms for D-jet reconstruction efficiency and Response Matrices.
Author: Preeti Dhankher(pdhankher@berkeley.edu)
'''
import numpy as np
import argparse
import os
from array import array
import pandas as pd


#for reference, see rey's code
#https://github.com/reynier0611/pyjetty/tree/master/pyjetty/alice_analysis/process/user/rey

#command for test job:
#./process_mc_RMalltrackjets.py -c ../../../config/hf_jetaxes/configcutsRMpreeti_ptbin.yaml -f /rstorage/alice/data/LHC18b8/569/LHC18b8_fast/20/282306/0003/AnalysisResults.root

from pyjetty.alice_analysis.process.base import process_io, process_utils, jet_info, process_base
import pyjetty.alihfjets.dev.hfjet.process.user.hf_jetaxes.process_io_mc_hf as hfdio
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
        self.initialize_config()
        print(self.config_file)
        self.overflow_gen = False #check
        self.overflow_det = False #check

    def HFAnalysis(self):
    
        hfaio = hfdio.HFAnalysisIO(input_file=self.input_file)
    
        #apply event and D0 selection cuts on tree
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
    
        start_time = time.time()
        print('--- {} seconds ---'.format(time.time() - start_time))
        
        # Initialize configuration and histograms
        self.initialize_config()
        self.initialize_user_output_objects()
        
        print("dataframe will be created jet analysis")
        
        self.df_fj_particles_gen = hfaio.load_data(m=self.track_mass, random_mass=self.track_random_mass, Is_treff_sys=self.Is_treff_sys, IsGen=True )
        # Find jets and fill histograms
        print("print generated D0 and track group")
        print(self.df_fj_particles_gen)
        
        self.df_fj_particles_rec = hfaio.load_data(m=self.track_mass, random_mass=self.track_random_mass, Is_treff_sys=self.Is_treff_sys, IsGen = False )
        print("print reconstructed D0 and track group")
        print(self.df_fj_particles_rec)
        
        print("Merge det-level and truth-level into a single dataframe grouped by event...")
        
        if len(self.df_fj_particles_gen)>0 and len(self.df_fj_particles_rec)>0:
            self.df_fj_merge_particles = pd.merge(self.df_fj_particles_gen, self.df_fj_particles_rec , on=self.unique_identifier)
            print(self.df_fj_merge_particles)
        else:
            print("No merge dataframe available")
            print("return empty dataframe")
            self.df_fj_merge_particles = pd.DataFrame()
            
        
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
        self.Is_treff_sys = self.config["Is_treff_sys"]
        # Mass assumption for jet reconstruction
        self.track_mass = self.config["track_mass"]
        self.track_random_mass = self.config["track_random_mass"]   # randomly assign K or p mass to tracks
        self.observable_list = self.config["observable_list"]
        self.sd_zcut = self.config["sd_zcut"]
        self.sd_beta = self.config["sd_beta"]
       

        # Detector kinematics
        self.eta_max = 0.9
    
    #---------------------------------------------------------------
    # Initialize histograms
    #---------------------------------------------------------------
    def initialize_user_output_objects(self):
    
        self.hTrackEtaPhi = ROOT.TH2F('hTrackEtaPhi', 'hTrackEtaPhi_det', 200, -1., 1., 628, 0., 6.28)
        self.hTrackPt = ROOT.TH1F('hTrackPt', 'hTrackPt_det', 300, 0., 300.)
        
    
        # Response Matrix THnSparse
                
        title = [ '#it{p}_{T,det}^{ch jet}', '#it{p}_{T,truth}^{ch jet}', '#it{p}_{T,det}^{D^{0}}', '#it{p}_{T,truth}^{D^{0}}' , '#it{#Delta R}_{det}', '#it{#Delta R}_{truth}']
        
        dim = 6
        #below is wrong!
        #delta_R = np.linspace(0, 0.8, 161)
        #below is correct!
        delta_Rbins = np.linspace(-0.1, 0.8, 181) #for bins of 0.005
        #delta_Rbins = np.linspace(-0.05, 1.5, 311)
        #below is for STD-WTA, WTA-SD, STD-D, SD-D and WTA-D (from Rey's AN)
        #delta_R = np.array([-0.02,0.,0.01,0.02,0.03,0.04,0.06,0.08,0.12,0.2])
        #delta_Rbins = np.array([-0.05,-0.04, -0.02, 0.0, 0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.12, 0.20])
        #below is for STD-SD
        #delta_R = [-0.003,0.0,0.0025,0.005,0.0075,0.01,0.0125,0.015,0.0175,0.02,0.025,0.0275,0.035]
        print("delta_Rbins ==")
        print(delta_Rbins)
        #print("delta_Rog ==")
        #print(delta_Rog)

        #originally at linspace(0,60,61)
        pt_bins_jet = np.linspace(0, 60, 61)
        pt_bins_dmeson = np.linspace(0, 60, 61)
        print("pt_bins_jet ==")
        print(pt_bins_jet)
        dmeson_mass_bin = np.linspace(1.7, 2.07, 371)
        
        
        nbins  = [len(pt_bins_jet)-1, len(pt_bins_jet)-1, len(pt_bins_dmeson)-1, len(pt_bins_dmeson)-1, len(delta_Rbins)-1, len(delta_Rbins)-1]
        min_li = [pt_bins_jet[0], pt_bins_jet[0],  pt_bins_dmeson[0], pt_bins_dmeson[0], delta_Rbins[0], delta_Rbins[0]]
        max_li = [pt_bins_jet[-1],  pt_bins_jet[0],   pt_bins_dmeson[-1], pt_bins_dmeson[-1], delta_Rbins[-1], delta_Rbins[-1]]

        nbins = (nbins)
        xmin = (min_li)
        xmax = (max_li)

        nbins_array = array('i', nbins)
        xmin_array = array('d', xmin)
        xmax_array = array('d', xmax)
        
      
        
        for observable in self.observable_list:
        
            name = 'hsparse_R{}_{}'.format(self.jetR_list[0], observable)
            h = ROOT.THnSparseD( name, name, dim,  nbins_array, xmin_array, xmax_array)
            h.Sumw2()
          
            for axis in range(0,dim):
                h.GetAxis(axis).SetTitle(title[axis])
                if axis == 0 or axis == 1:
                    h.SetBinEdges(axis, pt_bins_jet)
                if axis == 2 or axis == 3:
                    h.SetBinEdges(axis, pt_bins_dmeson)
                if axis == 4 or axis == 5:
                    h.SetBinEdges(axis, delta_Rbins)
                    
            setattr(self, name, h)
            
            name2 = 'hdeltaR_fd_gen_{}_{}'.format(self.jetR_list[0], observable)
            h2 = ROOT.TH1F(name2, name2, len(delta_Rbins)-1, delta_Rbins)
            h2.Sumw2()
            h2.AddBinContent(0)
            
            setattr(self, name2, h2) 

            name3 = 'hdeltaR_fd_det_{}_{}'.format(self.jetR_list[0], observable)
            h3 = ROOT.TH1F(name3, name3, len(delta_Rbins)-1, delta_Rbins)
            h3.Sumw2()
            h3.AddBinContent(0)

            setattr(self, name3, h3)

            name4 = 'hdeltaR_prompt_gen_{}_{}'.format(self.jetR_list[0], observable)
            h4 = ROOT.TH1F(name4, name4, len(delta_Rbins)-1, delta_Rbins)
            h4.Sumw2()
            h4.AddBinContent(0)

            setattr(self, name4, h4)

            name5 = 'hdeltaR_prompt_det_{}_{}'.format(self.jetR_list[0], observable)
            h5 = ROOT.TH1F(name5, name5, len(delta_Rbins)-1, delta_Rbins)
            h5.Sumw2()
            h5.AddBinContent(0)

            setattr(self, name5, h5)
        
        self.hjetvsDpt_fd_gen_STD_D = ROOT.TH2F('hjetvsDpt_fd_gen', 'hjetvsDpt_fd_gen', 60, pt_bins_jet, 60, pt_bins_dmeson)
        self.hjetvsDpt_prompt_gen_STD_D = ROOT.TH2F('hjetvsDpt_prompt_gen', 'hjetvsDpt_prompt_gen', 60, pt_bins_jet, 60, pt_bins_dmeson) 
         
        self.hjetvsDpt_fd_det_STD_D = ROOT.TH2F('hjetvsDpt_fd_det', 'hjetvsDpt_fd_det', 60, pt_bins_jet, 60, pt_bins_dmeson)
        self.hjetvsDpt_prompt_det_STD_D = ROOT.TH2F('hjetvsDpt_prompt_det', 'hjetvsDpt_prompt_det', 60, pt_bins_jet, 60, pt_bins_dmeson)

        #self.hjetvsDpt_fd_gen_WTA_D = ROOT.TH2F('hjetvsDpt_fd_gen', 'hjetvsDpt_fd_gen', 60, pt_bins_jet, 60, pt_bins_dmeson)
        #self.hjetvsDpt_prompt_gen_WTA_D = ROOT.TH2F('hjetvsDpt_prompt_gen', 'hjetvsDpt_prompt_gen', 60, pt_bins_jet, 60, pt_bins_dmeson)

        #self.hjetvsDpt_fd_det_WTA_D = ROOT.TH2F('hjetvsDpt_fd_det', 'hjetvsDpt_fd_det', 60, pt_bins_jet, 60, pt_bins_dmeson)
        #self.hjetvsDpt_prompt_det_WTA_D = ROOT.TH2F('hjetvsDpt_prompt_det', 'hjetvsDpt_prompt_det', 60, pt_bins_jet, 60, pt_bins_dmeson)

        #self.hjetvsDpt_fd_gen_SD_D = ROOT.TH2F('hjetvsDpt_fd_gen', 'hjetvsDpt_fd_gen', 60, pt_bins_jet, 60, pt_bins_dmeson)
        #self.hjetvsDpt_prompt_gen_SD_D = ROOT.TH2F('hjetvsDpt_prompt_gen', 'hjetvsDpt_prompt_gen', 60, pt_bins_jet, 60, pt_bins_dmeson)

        #self.hjetvsDpt_fd_det_SD_D = ROOT.TH2F('hjetvsDpt_fd_det', 'hjetvsDpt_fd_det', 60, pt_bins_jet, 60, pt_bins_dmeson)
        #self.hjetvsDpt_prompt_det_SD_D = ROOT.TH2F('hjetvsDpt_prompt_det', 'hjetvsDpt_prompt_det', 60, pt_bins_jet, 60, pt_bins_dmeson)

        #self.hjetvsDpt_fd_gen_STD_WTA = ROOT.TH2F('hjetvsDpt_fd_gen', 'hjetvsDpt_fd_gen', 60, pt_bins_jet, 60, pt_bins_dmeson)
        #self.hjetvsDpt_prompt_gen_STD_WTA = ROOT.TH2F('hjetvsDpt_prompt_gen', 'hjetvsDpt_prompt_gen', 60, pt_bins_jet, 60, pt_bins_dmeson)

        #self.hjetvsDpt_fd_det_STD_WTA = ROOT.TH2F('hjetvsDpt_fd_det', 'hjetvsDpt_fd_det', 60, pt_bins_jet, 60, pt_bins_dmeson)
        #self.hjetvsDpt_prompt_det_STD_WTA = ROOT.TH2F('hjetvsDpt_prompt_det', 'hjetvsDpt_prompt_det', 60, pt_bins_jet, 60, pt_bins_dmeson)

        #self.hjetvsDpt_fd_gen_STD_SD = ROOT.TH2F('hjetvsDpt_fd_gen', 'hjetvsDpt_fd_gen', 60, pt_bins_jet, 60, pt_bins_dmeson)
        #self.hjetvsDpt_prompt_gen_STD_SD = ROOT.TH2F('hjetvsDpt_prompt_gen', 'hjetvsDpt_prompt_gen', 60, pt_bins_jet, 60, pt_bins_dmeson)

        #self.hjetvsDpt_fd_det_STD_SD = ROOT.TH2F('hjetvsDpt_fd_det', 'hjetvsDpt_fd_det', 60, pt_bins_jet, 60, pt_bins_dmeson)
        #self.hjetvsDpt_prompt_det_STD_SD = ROOT.TH2F('hjetvsDpt_prompt_det', 'hjetvsDpt_prompt_det', 60, pt_bins_jet, 60, pt_bins_dmeson)

        #self.hjetvsDpt_fd_gen_WTA_SD = ROOT.TH2F('hjetvsDpt_fd_gen', 'hjetvsDpt_fd_gen', 60, pt_bins_jet, 60, pt_bins_dmeson)
        #self.hjetvsDpt_prompt_gen_WTA_SD = ROOT.TH2F('hjetvsDpt_prompt_gen', 'hjetvsDpt_prompt_gen', 60, pt_bins_jet, 60, pt_bins_dmeson)

        #self.hjetvsDpt_fd_det_WTA_SD = ROOT.TH2F('hjetvsDpt_fd_det', 'hjetvsDpt_fd_det', 60, pt_bins_jet, 60, pt_bins_dmeson)
        #self.hjetvsDpt_prompt_det_WTA_SD = ROOT.TH2F('hjetvsDpt_prompt_det', 'hjetvsDpt_prompt_det', 60, pt_bins_jet, 60, pt_bins_dmeson)

        self.hjetvsDptvsm_refl_STD_D = ROOT.TH3F('hjetvsDptvsm_refl', 'hjetvsDptvsm_refl', 60, pt_bins_jet, 60, pt_bins_dmeson, 370, dmeson_mass_bin)
        #self.hjetvsDptvsm_refl_WTA_D = ROOT.TH3F('hjetvsDptvsm_refl', 'hjetvsDptvsm_refl', 60, pt_bins_jet, 60, pt_bins_dmeson, 370, dmeson_mass_bin)
        #self.hjetvsDptvsm_refl_SD_D = ROOT.TH3F('hjetvsDptvsm_refl', 'hjetvsDptvsm_refl', 60, pt_bins_jet, 60, pt_bins_dmeson, 370, dmeson_mass_bin)
        #self.hjetvsDptvsm_refl_STD_WTA = ROOT.TH3F('hjetvsDptvsm_refl', 'hjetvsDptvsm_refl', 60, pt_bins_jet, 60, pt_bins_dmeson, 370, dmeson_mass_bin)
        #self.hjetvsDptvsm_refl_STD_SD = ROOT.TH3F('hjetvsDptvsm_refl', 'hjetvsDptvsm_refl', 60, pt_bins_jet, 60, pt_bins_dmeson, 370, dmeson_mass_bin)
        #self.hjetvsDptvsm_refl_WTA_SD = ROOT.TH3F('hjetvsDptvsm_refl', 'hjetvsDptvsm_refl', 60, pt_bins_jet, 60, pt_bins_dmeson, 370, dmeson_mass_bin)

    #---------------------------------------------------------------
    # Fill track histograms.
    #---------------------------------------------------------------
    def fillTrackHistograms(self, fj_particles_det):
        # Check that the entries exist appropriately
        # (need to check how this can happen -- but it is only a tiny fraction of events)
    
        for track in fj_particles_det:
            self.hTrackEtaPhi.Fill(track.eta(), track.phi())
            self.hTrackPt.Fill(track.pt())


    def analyzeEvents(self):
    
        if len(self.df_fj_particles_rec)>0:
            #filling some QA histograms
            #dont forget to seperate single track from multi track
            for fj_particles_det in self.df_fj_particles_rec['fj_D0_particles']:
                self.fillTrackHistograms(fj_particles_det)
    
        fj.ClusterSequence.print_banner()
        print()
                
        for jetR in self.jetR_list:
            print("radius of jet", jetR)
            #jet_def for Djets
            jet_def = jet_analysis.JetAnalysis(jet_R = jetR,  jet_RecombinationScheme=fj.E_scheme, particle_eta_max=self.eta_max, jet_pt_min=0.0, explicit_ghosts=False)
            #need to seperate out the single track jets from the multi-track jets (so new RMs for each)  
            
            if len(self.df_fj_particles_gen)>0:
                # Fill information all generated jets
                for fj_particles_truth, Isfd, Isprompt in \
                    zip(self.df_fj_particles_gen['fj_D0_particles_gen'], self.df_fj_particles_gen['ismcfd_gen'], self.df_fj_particles_gen['ismcprompt_gen']):
                    self.overflow_gen = False
                    self.overflow_det = False
                    jets_truth = self.analyze_jets(fj_particles_truth, jet_def, False, False, False, jetR, True)
                    Dcand_truth = self.analyze_jets(fj_particles_truth, jet_def, True, False, False, jetR, True)
                    
                    if jets_truth:
                        if Isfd:
                            #generated info used to calculate reconstruction efficiency of fd
                            self.hjetvsDpt_fd_gen_STD_D.Fill(jets_truth.pt(), Dcand_truth.pt())
                        if Isprompt:
                            #generated info used to calculate reconstruction efficiency of prompt
                            self.hjetvsDpt_prompt_gen_STD_D.Fill(jets_truth.pt(), Dcand_truth.pt()) 

                    '''
                    if jets_SD_truth:
                        if Isfd:
                            #generated info used to calculate reconstruction efficiency of fd
                            self.hjetvsDpt_fd_gen_WTA_D.Fill(jets_WTA_truth.pt(), Dcand_truth.pt())
                        if Isprompt:
                            #generated info used to calculate reconstruction efficiency of prompt
                            self.hjetvsDpt_prompt_gen_WTA_D.Fill(jets_WTA_truth.pt(), Dcand_truth.pt())

                    if jets_SD_truth:
                        if Isfd:
                            #generated info used to calculate reconstruction efficiency of fd
                            self.hjetvsDpt_fd_gen_SD_D.Fill(jets_SD_truth.pt(), Dcand_truth.pt())
                        if Isprompt:
                            #generated info used to calculate reconstruction efficiency of prompt
                            self.hjetvsDpt_prompt_gen_SD_D.Fill(jets_SD_truth.pt(), Dcand_truth.pt())
                    '''        
                            
            if not self.df_fj_merge_particles.empty:
                # jet reconstruction and matching
                for fj_particles_truth, Isfd_gen, Isprompt_gen, fj_particles_det, Isfd, Isprompt, Isrefl  in \
                zip(self.df_fj_merge_particles['fj_D0_particles_gen'], self.df_fj_merge_particles['ismcfd_gen'], self.df_fj_merge_particles['ismcprompt_gen'], self.df_fj_merge_particles['fj_D0_particles'], self.df_fj_merge_particles['ismcfd'], self.df_fj_merge_particles['ismcprompt'],  self.df_fj_merge_particles['ismcrefl']):
                    #need to seperate single track jets from multitrack jets 
                    #>1 is multitrack and =1 is singletrack
                    self.overflow_gen = False
                    self.overflow_det = False
                    #initially set the overflow to false
                    jets_det = self.analyze_jets(fj_particles_det, jet_def, False, False, False, jetR, False)
                    jets_truth = self.analyze_jets(fj_particles_truth, jet_def, False, False, False, jetR, True)
                    
                    Dcand_det = self.analyze_jets(fj_particles_det, jet_def, True, False, False, jetR, False)
                    Dcand_truth = self.analyze_jets(fj_particles_truth, jet_def, True, False, False, jetR, True)
                    
                    jets_WTA_det = self.analyze_jets(fj_particles_det, jet_def, False, True, False, jetR, False)
                    jets_WTA_truth = self.analyze_jets(fj_particles_truth, jet_def, False, True, False, jetR, True)
                    
                    #print('checkpoint, what is oveflow? ', self.overflow)
                    jets_SD_det = self.analyze_jets(fj_particles_det, jet_def, False, False, True, jetR, False)
                    
                    jets_SD_truth = self.analyze_jets(fj_particles_truth, jet_def, False, False, True, jetR, True)
                    
                    #now do jet matching and then fill histograms!

                    STD_D_matching = self.jet_matching(jets_truth,Dcand_truth,jets_det,Dcand_det,jetR,True,False)
                    
                    WTA_D_matching = self.jet_matching(jets_WTA_truth,Dcand_truth,jets_WTA_det,Dcand_det,jetR,True,False)
                    
                    SD_D_matching = self.jet_matching(jets_SD_truth,Dcand_truth,jets_SD_det,Dcand_det,jetR,True,False)
                    
                    STD_WTA_matching = self.jet_matching(jets_truth,jets_WTA_truth,jets_det,jets_WTA_det,jetR,False,True)
                    
                    STD_SD_matching = self.jet_matching(jets_truth,jets_SD_truth,jets_det,jets_SD_det,jetR,False,True)
                    
                    WTA_SD_matching = self.jet_matching(jets_WTA_truth,jets_SD_truth,jets_WTA_det,jets_SD_det,jetR,False,True)
                    
                    
                    if STD_D_matching == True:
                        print("matched")

                        print(jets_det.delta_R(Dcand_det))
                        if Isrefl:
                            print("matched candidate is reflection")
                            self.hjetvsDptvsm_refl_STD_D.Fill(jets_det.pt(), Dcand_det.pt(), Dcand_det.m()) #info used to create reflection templates

                        elif Isfd:
                            print("STD_D matched candidate is feeddown")
                            #reconstructed info used to calculate reconstruction efficiency of fd
                            print('overflow? ', self.overflow_gen, ', ', self.overflow_det)
                            print('STD_D fd deltaR = ', jets_truth.delta_R(Dcand_truth))
                            self.hjetvsDpt_fd_det_STD_D.Fill(jets_det.pt(), Dcand_det.pt())
                            #self.hjetvsDpt_fd_gen_STD_D.Fill(jets_truth.pt(), Dcand_truth.pt())
                            #STD-D0
                            x_array = array( 'd', ( jets_det.pt(), jets_truth.pt(), Dcand_det.pt(), Dcand_truth.pt(), jets_det.delta_R(Dcand_det), jets_truth.delta_R(Dcand_truth)))
                            getattr(self, 'hsparse_R{}_{}'.format(jetR, self.observable_list[1])).Fill(x_array) #also make sure here the observable_list is calling the appropriate histogram
                            getattr(self, 'hdeltaR_fd_gen_{}_{}'.format(jetR, self.observable_list[1])).Fill(jets_truth.delta_R(Dcand_truth))
                            getattr(self, 'hdeltaR_fd_det_{}_{}'.format(jetR, self.observable_list[1])).Fill(jets_det.delta_R(Dcand_det))    
         
                        elif Isprompt:
                            print("STD_D matched candidate is prompt")
                            #reconstructed info used to calculate reconstruction efficiency of prompt
                            print('overflow? ', self.overflow_gen, ', ', self.overflow_det)
                            self.hjetvsDpt_prompt_det_STD_D.Fill(jets_det.pt(), Dcand_det.pt())
                            #self.hjetvsDpt_prompt_gen_STD_D.Fill(jets_truth.pt(), Dcand_truth.pt())
                            #STD-D0
                            x_array = array( 'd', ( jets_det.pt(), jets_truth.pt(), Dcand_det.pt(), Dcand_truth.pt(), jets_det.delta_R(Dcand_det), jets_truth.delta_R(Dcand_truth)))
                            getattr(self, 'hsparse_R{}_{}'.format(jetR, self.observable_list[0])).Fill(x_array)
                            getattr(self, 'hdeltaR_prompt_gen_{}_{}'.format(jetR, self.observable_list[1])).Fill(jets_truth.delta_R(Dcand_truth))
                            getattr(self, 'hdeltaR_prompt_det_{}_{}'.format(jetR, self.observable_list[1])).Fill(jets_det.delta_R(Dcand_det))

                    if WTA_D_matching == True:
                        print("matched") 
                        if Isrefl:
                            print("WTA_D matched candidate is reflection")
                            #self.hjetvsDptvsm_refl_WTA_D.Fill(jets_WTA_det.pt(), Dcand_det.pt(), Dcand_det.m()) #info used to create reflection templates
                        elif Isfd:
                            print("WTA_D matched candidate is feeddown")
                            #reconstructed info used to calculate reconstruction efficiency of fd
                            #self.hjetvsDpt_fd_det_WTA_D.Fill(jets_WTA_det.pt(), Dcand_det.pt())
                            #self.hjetvsDpt_fd_gen_WTA_D.Fill(jets_WTA_truth.pt(), Dcand_truth.pt())
                            #WTA-D0
                            
                            print('WTA_D fd deltaR = ', jets_WTA_truth.delta_R(Dcand_truth))
                            x_array = array( 'd', ( jets_det.pt(), jets_truth.pt(), Dcand_det.pt(), Dcand_truth.pt(), jets_WTA_det.delta_R(Dcand_det), jets_WTA_truth.delta_R(Dcand_truth)))
                            getattr(self, 'hsparse_R{}_{}'.format(jetR, self.observable_list[3])).Fill(x_array)
                            getattr(self, 'hdeltaR_fd_gen_{}_{}'.format(jetR, self.observable_list[3])).Fill(jets_WTA_truth.delta_R(Dcand_truth))
                            getattr(self, 'hdeltaR_fd_det_{}_{}'.format(jetR, self.observable_list[3])).Fill(jets_WTA_det.delta_R(Dcand_det))
                        elif Isprompt:  
                            print("matched candidate is prompt")
                            #reconstructed info used to calculate reconstruction efficiency of prompt
                            #self.hjetvsDpt_prompt_det_WTA_D.Fill(jets_WTA_det.pt(), Dcand_det.pt())
                            #self.hjetvsDpt_prompt_gen_WTA_D.Fill(jets_WTA_truth.pt(), Dcand_truth.pt())
                            #WTA-D0
                            x_array = array( 'd', ( jets_det.pt(), jets_truth.pt(), Dcand_det.pt(), Dcand_truth.pt(), jets_WTA_det.delta_R(Dcand_det), jets_WTA_truth.delta_R(Dcand_truth)))
                            getattr(self, 'hsparse_R{}_{}'.format(jetR, self.observable_list[2])).Fill(x_array)
                            getattr(self, 'hdeltaR_prompt_gen_{}_{}'.format(jetR, self.observable_list[3])).Fill(jets_WTA_truth.delta_R(Dcand_truth))
                            getattr(self, 'hdeltaR_prompt_det_{}_{}'.format(jetR, self.observable_list[3])).Fill(jets_WTA_det.delta_R(Dcand_det))

                    if SD_D_matching == True:
                        print("matched")
                        if Isrefl:
                            print("SD_D matched candidate is reflection")
                            #self.hjetvsDptvsm_refl_SD_D.Fill(jets_SD_det.pt(), Dcand_det.pt(), Dcand_det.m()) #info used to create reflection templates
                        elif Isfd:
                            print("SD_D matched candidate is feeddown")
                            #reconstructed info used to calculate reconstruction efficiency of fd
                            #self.hjetvsDpt_fd_det_SD_D.Fill(jets_SD_det.pt(), Dcand_det.pt())
                            #self.hjetvsDpt_fd_gen_SD_D.Fill(jets_SD_truth.pt(), Dcand_truth.pt())
                            #SD-D0

                            deltaR_SD_D_gen = jets_SD_truth.delta_R(Dcand_truth)
                            deltaR_SD_D_det = jets_SD_det.delta_R(Dcand_det)
                            print('overflow? ', self.overflow_gen, ', ', self.overflow_det)
                            if self.overflow_gen==True:
                                deltaR_SD_D_gen = -0.1
                            if self.overflow_det==True:
                                deltaR_SD_D_det = -0.1
                            
                            print("deltaR_SD_D_gen = ",deltaR_SD_D_gen)
                            print("deltaR_SD_D_det = ",deltaR_SD_D_det)

                            x_array = array( 'd', ( jets_det.pt(), jets_truth.pt(), Dcand_det.pt(), Dcand_truth.pt(), deltaR_SD_D_det, deltaR_SD_D_gen))
                            getattr(self, 'hsparse_R{}_{}'.format(jetR, self.observable_list[5])).Fill(x_array)
                            getattr(self, 'hdeltaR_fd_gen_{}_{}'.format(jetR, self.observable_list[5])).Fill(deltaR_SD_D_gen)
                            getattr(self, 'hdeltaR_fd_det_{}_{}'.format(jetR, self.observable_list[5])).Fill(deltaR_SD_D_det)
                        elif Isprompt:
                            print("SD_D matched candidate is prompt")
                            #reconstructed info used to calculate reconstruction efficiency of prompt
                            #self.hjetvsDpt_prompt_det_SD_D.Fill(jets_SD_det.pt(), Dcand_det.pt())
                            #self.hjetvsDpt_prompt_gen_SD_D.Fill(jets_SD_truth.pt(), Dcand_truth.pt())
                            #SD-D0
                            deltaR_SD_D_gen = jets_SD_truth.delta_R(Dcand_truth)
                            deltaR_SD_D_det = jets_SD_det.delta_R(Dcand_det)
                            print('overflow? ', self.overflow_gen, ', ', self.overflow_det)
                            if self.overflow_gen==True:
                                deltaR_SD_D_gen = -0.1
                            if self.overflow_det==True:
                                deltaR_SD_D_det = -0.1
                            
                            print("deltaR_SD_D_gen = ",deltaR_SD_D_gen)
                            print("deltaR_SD_D_det = ",deltaR_SD_D_det)
                            x_array = array( 'd', ( jets_det.pt(), jets_truth.pt(), Dcand_det.pt(), Dcand_truth.pt(), deltaR_SD_D_det, deltaR_SD_D_gen))
                            getattr(self, 'hsparse_R{}_{}'.format(jetR, self.observable_list[4])).Fill(x_array)
                            getattr(self, 'hdeltaR_prompt_gen_{}_{}'.format(jetR, self.observable_list[5])).Fill(deltaR_SD_D_gen)
                            getattr(self, 'hdeltaR_prompt_det_{}_{}'.format(jetR, self.observable_list[5])).Fill(deltaR_SD_D_det)

                    if STD_WTA_matching == True:
                        print("matched")
                        if Isrefl:
                            print("matched candidate is reflection")
                            #self.hjetvsDptvsm_refl_STD_WTA.Fill(jets_det.pt(), Dcand_det.m()) #info used to create reflection templates
                        elif Isfd:
                            print("matched candidate is feeddown")
                            #reconstructed info used to calculate reconstruction efficiency of fd
                            #self.hjetvsDpt_fd_det_STD_WTA.Fill(jets_WTA_det.pt(), Dcand_det.pt())
                            #self.hjetvsDpt_fd_gen_STD_WTA.Fill(jets_WTA_truth.pt(), Dcand_truth.pt())
                            #STD-WTA
                            print('STD_WTA fd deltaR = ', jets_truth.delta_R(jets_WTA_truth))
                            x_array = array( 'd', ( jets_det.pt(), jets_truth.pt(), Dcand_det.pt(), Dcand_truth.pt(), jets_det.delta_R(jets_WTA_det), jets_truth.delta_R(jets_WTA_truth)))
                            getattr(self, 'hsparse_R{}_{}'.format(jetR, self.observable_list[7])).Fill(x_array)
                            getattr(self, 'hdeltaR_fd_gen_{}_{}'.format(jetR, self.observable_list[7])).Fill(jets_truth.delta_R(jets_WTA_truth))
                            getattr(self, 'hdeltaR_fd_det_{}_{}'.format(jetR, self.observable_list[7])).Fill(jets_det.delta_R(jets_WTA_det))
                        elif Isprompt:
                            print("matched candidate is prompt")
                            #reconstructed info used to calculate reconstruction efficiency of prompt
                            #self.hjetvsDpt_prompt_det_STD_WTA.Fill(jets_WTA_det.pt(), Dcand_det.pt())
                            #self.hjetvsDpt_prompt_gen_STD_WTA.Fill(jets_WTA_truth.pt(), Dcand_truth.pt())
                            #STD-WTA
                            x_array = array( 'd', ( jets_det.pt(), jets_truth.pt(), Dcand_det.pt(), Dcand_truth.pt(), jets_det.delta_R(jets_WTA_det), jets_truth.delta_R(jets_WTA_truth)))
                            getattr(self, 'hsparse_R{}_{}'.format(jetR, self.observable_list[6])).Fill(x_array)
                            getattr(self, 'hdeltaR_prompt_gen_{}_{}'.format(jetR, self.observable_list[7])).Fill(jets_truth.delta_R(jets_WTA_truth))
                            getattr(self, 'hdeltaR_prompt_det_{}_{}'.format(jetR, self.observable_list[7])).Fill(jets_det.delta_R(jets_WTA_det))

                    if STD_SD_matching == True:
                        print("matched")
                        if Isrefl:
                            print("STD_SD matched candidate is reflection")
                            #self.hjetvsDptvsm_refl_STD_SD.Fill(jets_WTA_det.pt(), Dcand_det.pt(), Dcand_det.m()) #info used to create reflection templates
                        elif Isfd:
                            print("STD_SD matched candidate is feeddown")
                            #reconstructed info used to calculate reconstruction efficiency of fd
                            #self.hjetvsDpt_fd_det_STD_SD.Fill(jets_WTA_det.pt(), Dcand_det.pt())
                            #self.hjetvsDpt_fd_gen_STD_SD.Fill(jets_WTA_truth.pt(), Dcand_truth.pt())
                            #STD-SD
                            deltaR_STD_SD_gen = jets_truth.delta_R(jets_SD_truth)
                            deltaR_STD_SD_det = jets_det.delta_R(jets_SD_det)
                            print('overflow? ', self.overflow_gen, ', ', self.overflow_det)
                            if self.overflow_gen==True:
                                deltaR_STD_SD_gen = -0.1
                            if self.overflow_det==True:
                                deltaR_STD_SD_det = -0.1
                            print("deltaR_STD_SD_gen = ",deltaR_STD_SD_gen)
                            print("deltaR_STD_SD_det = ",deltaR_STD_SD_det)
                            x_array = array( 'd', ( jets_det.pt(), jets_truth.pt(), Dcand_det.pt(), Dcand_truth.pt(), deltaR_STD_SD_det, deltaR_STD_SD_gen))
                            getattr(self, 'hsparse_R{}_{}'.format(jetR, self.observable_list[9])).Fill(x_array)
                            getattr(self, 'hdeltaR_fd_gen_{}_{}'.format(jetR, self.observable_list[9])).Fill(deltaR_STD_SD_gen)
                            getattr(self, 'hdeltaR_fd_det_{}_{}'.format(jetR, self.observable_list[9])).Fill(deltaR_STD_SD_det)
                        elif Isprompt:
                            print("STD_SD matched candidate is prompt")
                            #reconstructed info used to calculate reconstruction efficiency of prompt
                            #self.hjetvsDpt_prompt_det_STD_SD.Fill(jets_WTA_det.pt(), Dcand_det.pt())
                            #self.hjetvsDpt_prompt_gen_STD_SD.Fill(jets_WTA_truth.pt(), Dcand_truth.pt())
                            #STD-SD
                            deltaR_STD_SD_gen = jets_truth.delta_R(jets_SD_truth)
                            deltaR_STD_SD_det = jets_det.delta_R(jets_SD_det)
                            print('overflow? ', self.overflow_gen, ', ', self.overflow_det)
                            if self.overflow_gen==True:
                                deltaR_STD_SD_gen = -0.1
                            if self.overflow_det==True:
                                deltaR_STD_SD_det = -0.1
                            print("deltaR_STD_SD_gen = ",deltaR_STD_SD_gen)
                            print("deltaR_STD_SD_det = ",deltaR_STD_SD_det)
                            x_array = array( 'd', ( jets_det.pt(), jets_truth.pt(), Dcand_det.pt(), Dcand_truth.pt(), deltaR_STD_SD_det, deltaR_STD_SD_gen))
                            getattr(self, 'hsparse_R{}_{}'.format(jetR, self.observable_list[8])).Fill(x_array)
                            getattr(self, 'hdeltaR_prompt_gen_{}_{}'.format(jetR, self.observable_list[9])).Fill(deltaR_STD_SD_gen)
                            getattr(self, 'hdeltaR_prompt_det_{}_{}'.format(jetR, self.observable_list[9])).Fill(deltaR_STD_SD_det)

                    if WTA_SD_matching == True:
                        print("matched")
                        if Isrefl:
                            print("WTA_SD matched candidate is reflection")
                            #self.hjetvsDptvsm_refl_WTA_SD.Fill(jets_WTA_det.pt(), Dcand_det.pt(), Dcand_det.m()) #info used to create reflection templates
                        elif Isfd:
                            print("WTA_SD matched candidate is feeddown")
                            #reconstructed info used to calculate reconstruction efficiency of fd
                            #self.hjetvsDpt_fd_det_WTA_SD.Fill(jets_WTA_det.pt(), Dcand_det.pt())
                            #self.hjetvsDpt_fd_gen_WTA_SD.Fill(jets_WTA_truth.pt(), Dcand_truth.pt())
                            #WTA-SD
                            deltaR_WTA_SD_gen = jets_WTA_truth.delta_R(jets_SD_truth)
                            deltaR_WTA_SD_det = jets_WTA_det.delta_R(jets_SD_det)
                            print('overflow? ', self.overflow_gen, ', ', self.overflow_det)
                            if self.overflow_gen==True:
                                deltaR_WTA_SD_gen = -0.1
                            if self.overflow_det==True:
                                deltaR_WTA_SD_det = -0.1
                            print('WTA_SD fd deltaR gen = ', deltaR_WTA_SD_gen)
                            print("deltaR_WTA_SD_det = ",deltaR_WTA_SD_det)
                            x_array = array( 'd', ( jets_det.pt(), jets_truth.pt(), Dcand_det.pt(), Dcand_truth.pt(), deltaR_WTA_SD_det, deltaR_WTA_SD_gen))
                            getattr(self, 'hsparse_R{}_{}'.format(jetR, self.observable_list[11])).Fill(x_array)
                            getattr(self, 'hdeltaR_fd_gen_{}_{}'.format(jetR, self.observable_list[11])).Fill(deltaR_WTA_SD_gen)
                            getattr(self, 'hdeltaR_fd_det_{}_{}'.format(jetR, self.observable_list[11])).Fill(deltaR_WTA_SD_det)
                        elif Isprompt:
                            print("WTA_SD matched candidate is prompt")
                            #reconstructed info used to calculate reconstruction efficiency of prompt
                            #self.hjetvsDpt_prompt_det_WTA_SD.Fill(jets_WTA_det.pt(), Dcand_det.pt())
                            #self.hjetvsDpt_prompt_gen_WTA_SD.Fill(jets_WTA_truth.pt(), Dcand_truth.pt())
                            #WTA-SD
                            deltaR_WTA_SD_gen = jets_WTA_truth.delta_R(jets_SD_truth)
                            deltaR_WTA_SD_det = jets_WTA_det.delta_R(jets_SD_det)
                            print('overflow? ', self.overflow_gen, ', ', self.overflow_det)
                            if self.overflow_gen==True:
                                deltaR_WTA_SD_gen = -0.1
                            if self.overflow_det==True:
                                deltaR_WTA_SD_det = -0.1
                            
                            print("deltaR_WTA_SD_gen = ",deltaR_WTA_SD_gen)
                            print("deltaR_WTA_SD_det = ",deltaR_WTA_SD_det)
                            x_array = array( 'd', ( jets_det.pt(), jets_truth.pt(), Dcand_det.pt(), Dcand_truth.pt(), deltaR_WTA_SD_det, deltaR_WTA_SD_gen))
                            getattr(self, 'hsparse_R{}_{}'.format(jetR, self.observable_list[10])).Fill(x_array)
                            getattr(self, 'hdeltaR_prompt_gen_{}_{}'.format(jetR, self.observable_list[11])).Fill(deltaR_WTA_SD_gen)
                            getattr(self, 'hdeltaR_prompt_det_{}_{}'.format(jetR, self.observable_list[11])).Fill(deltaR_WTA_SD_det)

    def jet_matching(self,jet1truth,jet2truth,jet1det,jet2det,jetR,IsDcand,IsJet):
        #note that jet1 is always a jet, but jet2 could be a jet or Dmeson. KEEP THIS ORDER!!!
        #check that first jet matches
        if jet1truth and jet1det:
            if jet1truth.delta_R(jet1det) < (0.6*jetR):
                #check for Dmeson matching if IsDcand is True
                if IsDcand and jet2truth and jet2det:
                    if jet2truth.delta_R(jet2det) < (0.1): #this should be 0.6*jetR - nevermind
                        return True
                    else:
                        return False
                #check for jet matching if IsJet is True - #check that i did this b4
                elif IsJet and jet2truth and jet2det:
                    if jet2truth.delta_R(jet2det) < (0.6*jetR):
                        return True
                    else:
                        return False
        else:
            return False


    def analyze_jets(self, fj_D0_particles, jet_def, IsDzero, IsWTA, IsSD, jetR, isTruth):
        # perform Djet reconstruction
        djmm = fjtools.DJetMatchMaker()
        _parts_and_ds = fj_D0_particles

        jet_def.analyze_event(_parts_and_ds)
        
        if len(jet_def.jets) < 1:
            return False
            
        jets = jet_def.jets_as_psj_vector()
        
        djets = djmm.filter_D0_jets(jets)

        #initialize overflow
        #self.overflow_gen = False #check
        #self.overflow_det = False

        if len(djets) > 0:
            j = djets[0]
            dcand = djmm.get_Dcand_in_jet(j)
            if len(j.constituents()) > 0:
                #print("number of jet constituents = ", len(j.constituents()))
                if IsDzero:
                    print("Dzero candidate!")
                    return dcand[0]
                
                elif IsWTA:
                    #print("################################ IS WTA #######################################")                
                    #jets with the winner take all axis################################
                    #print("##############is wta")
                    jet_def_wta = fj.JetDefinition(fj.cambridge_algorithm, 2*jetR)#originally 2*self.jetR_list[0])
                    #print("##############cambridge-achean applied")
                    jet_def_wta.set_recombination_scheme(fj.WTA_pt_scheme)
                    #print("##############recombination scheme applied")
                    reclusterer_wta =  fjcontrib.Recluster(jet_def_wta)
                    #print("##############reclustered")
                    jet_wta = reclusterer_wta.result(j)
                    #################################
                    print("Winner Take All with jet = ",jet_wta)
                    return jet_wta
                elif IsSD:
                    print("############################## IS SD ##############################")
                    ##jets with the SD axis################################
                    gshop = fjcontrib.GroomerShop(j, jetR, fj.cambridge_algorithm)
                    #gshop.soft_drop(beta, zcut, jetR)
                    jet_groomed_lund = gshop.soft_drop(self.sd_beta,self.sd_zcut,jetR)
                    jet_groomed = jet_groomed_lund.pair()
                    if jet_groomed_lund.Delta() < 0:
                        print("no groomed jet was returned :( ")
                        if isTruth == True:
                            self.overflow_gen = True
                        elif isTruth == False:
                            self.overflow_det = True
                        print('overflow (inside analyze_jets)= ', self.overflow_gen, ', ', self.overflow_det)
                        print('GROOMED AWAY JET INFO ****************************************')
                        #print('jet_groomed.constituents() = ', jet_groomed.constituents())
                        print('is Truth? ', isTruth)

                        return jet_groomed
                    else:
                        print("groomed jet returned :) ")
                        print("jet_groomed = ", jet_groomed)
                        #print('self.overflow = ', self.overflow)
                        return jet_groomed
                
                else:
                    print("Standard with jet = ",j)
                    return j
            else:
                print("number of jets is 0 or 1")
        if len(djets) > 1:
            perror("more than one jet per D candidate?")
            return False
            
        

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
    args = parser.parse_args()
    
    print('Configuring...')
    print('inputFile: \'{0}\''.format(args.inputFile))
    print('configFile: \'{0}\''.format(args.configFile))
    print('ouputDir: \'{0}\"'.format(args.outputDir))
    print('----------------------------------------------------------------')
  
    if not os.path.exists(args.configFile):
        print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)

    hfa = HFAnalysisInvMass(input_file = args.inputFile, config_file = args.configFile, output_dir = args.outputDir)
    print("Start reading input output and perform analysis")
    hfa.HFAnalysis()
  


