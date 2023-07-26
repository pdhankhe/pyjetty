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

    def HFAnalysis(self):

        hfaio = hfdio.HFAnalysisIO(input_file=self.input_file)

        #apply event and D0 selection cuts on tree
        print("Set event selection cuts here")
        hfa.event_selection.add_selection_range_abs('z_vtx_reco', 10)
        hfa.event_selection.add_selection_equal('is_ev_rej', 0)
        print("Set D0 pt independent selection cuts here")
        #topomatic cut suggested by D2H.
        #hfa.d0_selection.add_selection_range_abs('max_norm_d0d0exp',2)
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
        if len(self.df_fj_particles_gen) > 0 and len(self.df_fj_particles_rec) > 0:
            self.df_fj_merge_particles = pd.merge(self.df_fj_particles_gen, self.df_fj_particles_rec , on=self.unique_identifier)
            print(self.df_fj_merge_particles)

        #print('Merge det-level and truth-level into a single dataframe grouped by event...')
        #self.df_fj_particles = pd.concat([self.df_fj_particles_rec , self.df_fj_particles_gen], axis=1)
        #self.df_fj_particles = self.df_fj_particles.dropna()

        #for index, row in self.df_fj_particles.iterrows():
        #    print(row)
        #print('Find jets...')

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
        self.observable_RM_list = self.config["observable_RM_list"]
        self.ENC_pair_cut = self.config["ENC_pair_cut"]
        self.trk_thrd = self.config["trk_thrd"]
        # Detector kinematics
        self.eta_max = 0.9

    #---------------------------------------------------------------
    # Initialize histograms
    #---------------------------------------------------------------
    def initialize_user_output_objects(self):

                                        
        self.hTrackEtaPhi = ROOT.TH2F('hTrackEtaPhi', 'hTrackEtaPhi_det', 200, -1., 1., 628, 0., 6.28)
        self.hTrackPt = ROOT.TH1F('hTrackPt', 'hTrackPt_det', 300, 0., 300.)
        
        self.hNevents = ROOT.TH1F("hNevents", "hNevents", 2, array('f', [-0.5, 0.5, 1.5]) )
        
        self.hTrackEtaPhi_prompt = ROOT.TH2F("hTrackEtaPhi_prompt", "#eta#phi distribution of matched D0 candidate;#eta;#phi",600,-0.3,0.3,600,-0.3,0.3)
        self.hTrackEtaPhi_nonprompt = ROOT.TH2F("hTrackEtaPhi_nonprompt", "#eta#phi distribution of matched non prompt D0 candidate;#eta;#phi",600,-0.3,0.3,600,-0.3,0.3)
        
        self.hDeltaR_prompt = ROOT.TH1F('hDeltaR_prompt', 'Angular distance between D0', 20, 0., 0.1)
        self.hDeltaR_nonprompt = ROOT.TH1F('hDeltaR_nonprompt', 'Angular distance between non prompt D0', 20, 0., 0.1)

        # Eenergy energy correlator THnSparse
        title = [ '#it{p}_{T}^{ch jet}', '#it{p}_{T}^{D^{0}}', 'Invmass', '#it{R}_{L}']

        RL_bins  = np.logspace(np.log10(1E-4), np.log10(1), 51)
        print("print RL bins")
        #print(RL_bins)
        pt_bins_jet = np.linspace(0, 60, 61)
        pt_bins_dmeson = np.linspace(0, 60, 61)
        dmeson_mass_bin = np.linspace(1.7, 2.07, 371)

        dim = 4
        nbins  = [len(pt_bins_jet)-1, len(pt_bins_dmeson)-1, len(dmeson_mass_bin)-1, 50]
        min_li = [pt_bins_jet[0],     pt_bins_dmeson[0],     dmeson_mass_bin[0],     RL_bins[0]    ]
        max_li = [pt_bins_jet[-1],    pt_bins_dmeson[-1],    dmeson_mass_bin[-1],    RL_bins[-1]   ]

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
                if axis == 2:
                    h.SetBinEdges(axis, dmeson_mass_bin)
                if axis == 3:
                    h.SetBinEdges(axis, RL_bins)

            setattr(self, name, h)

        # Response Matrix THnSparse

        title = [ '#it{p}_{T,det}^{ch jet}', '#it{p}_{T,truth}^{ch jet}', '#it{p}_{T,det}^{D^{0}}', '#it{p}_{T,truth}^{D^{0}}']

        nbins  = [len(pt_bins_jet)-1, len(pt_bins_jet)-1, len(pt_bins_dmeson)-1, len(pt_bins_dmeson)-1]
        min_li = [pt_bins_jet[0], pt_bins_jet[0],  pt_bins_dmeson[0], pt_bins_dmeson[0] ]
        max_li = [pt_bins_jet[-1],  pt_bins_jet[0],   pt_bins_dmeson[-1], pt_bins_dmeson[-1] ]

        nbins = (nbins)
        xmin = (min_li)
        xmax = (max_li)

        nbins_array = array('i', nbins)
        xmin_array = array('d', xmin)
        xmax_array = array('d', xmax)



        for observable in self.observable_RM_list:

            name = 'hsparse_R{}_{}'.format(self.jetR_list[0], observable)
            h = ROOT.THnSparseD( name, name, dim,  nbins_array, xmin_array, xmax_array)
            h.Sumw2()

            for axis in range(0,dim):
                h.GetAxis(axis).SetTitle(title[axis])
                if axis == 0 or axis == 1:
                    h.SetBinEdges(axis, pt_bins_jet)
                if axis == 2 or axis == 3:
                    h.SetBinEdges(axis, pt_bins_dmeson)

            setattr(self, name, h)



        self.hjetvsDpt_fd_gen = ROOT.TH2F('hjetvsDpt_fd_gen', 'hjetvsDpt_fd_gen', 60, pt_bins_jet, 60, pt_bins_dmeson)
        self.hjetvsDpt_prompt_gen = ROOT.TH2F('hjetvsDpt_prompt_gen', 'hjetvsDpt_prompt_gen', 60, pt_bins_jet, 60, pt_bins_dmeson)


        self.hjetvsDpt_fd_det = ROOT.TH2F('hjetvsDpt_fd_det', 'hjetvsDpt_fd_det', 60, pt_bins_jet, 60, pt_bins_dmeson)
        self.hjetvsDpt_prompt_det = ROOT.TH2F('hjetvsDpt_prompt_det', 'hjetvsDpt_prompt_det', 60, pt_bins_jet, 60, pt_bins_dmeson)

        self.hjetvsDptvsm_refl = ROOT.TH3F('hjetvsDptvsm_refl', 'hjetvsDptvsm_refl', 60, pt_bins_jet, 60, pt_bins_dmeson, 370, dmeson_mass_bin)


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
            for fj_particles_det in self.df_fj_particles_rec['fj_D0_particles']:
                self.fillTrackHistograms(fj_particles_det)

        fj.ClusterSequence.print_banner()
        print()

        for jetR in self.jetR_list:
            print("radius of jet", jetR)
            jet_def = jet_analysis.JetAnalysis(jet_R = jetR,  jet_RecombinationScheme=fj.E_scheme, particle_eta_max=self.eta_max, jet_pt_min=0.0, explicit_ghosts=False)
            
            if len(self.df_fj_particles_gen)>0:
                for fj_particles_truth, Isfd, Isprompt in \
                    zip(self.df_fj_particles_gen['fj_D0_particles_gen'], self.df_fj_particles_gen['ismcfd_gen'], self.df_fj_particles_gen['ismcprompt_gen']):

                    jets_truth = self.analyze_jets(fj_particles_truth, jet_def, False, False)
                    Dcand_truth = self.analyze_jets(fj_particles_truth, jet_def, True, False)

                    if jets_truth:
                        if Isfd:
                            self.hjetvsDpt_fd_gen.Fill(jets_truth.pt(), Dcand_truth.pt())
                        if Isprompt:
                            self.hjetvsDpt_prompt_gen.Fill(jets_truth.pt(), Dcand_truth.pt())
                            
                            
            if len(self.df_fj_merge_particles)>0:
                for fj_particles_truth, Isfd_gen, Isprompt_gen, fj_particles_det, Isfd, Isprompt, Isrefl  in \
                    zip(self.df_fj_merge_particles['fj_D0_particles_gen'], self.df_fj_merge_particles['ismcfd_gen'], self.df_fj_merge_particles['ismcprompt_gen'], self.df_fj_merge_particles['fj_D0_particles'], self.df_fj_merge_particles['ismcfd'], self.df_fj_merge_particles['ismcprompt'],  self.df_fj_merge_particles['ismcrefl']):

                    jets_det = self.analyze_jets(fj_particles_det, jet_def, False, False)
                    jets_truth = self.analyze_jets(fj_particles_truth, jet_def, False, False)

                    Dcand_det = self.analyze_jets(fj_particles_det, jet_def, True, False)
                    Dcand_truth = self.analyze_jets(fj_particles_truth, jet_def, True, False)

                    Constituents_det = self.analyze_jets(fj_particles_det, jet_def, False, True)
                    Constituents_truth = self.analyze_jets(fj_particles_truth, jet_def, False, True)

                    
                    if jets_truth:
                        if jets_det:
                            if jets_truth.delta_R(jets_det) < (0.6*0.4): #jet matching
                                if Dcand_truth.delta_R(Dcand_det) < (0.1): # Dcandidate inside jet matching
                                    print("matched")
                                    etadiff = Dcand_det.eta() - Dcand_truth.eta()
                                    phidiff = Dcand_det.phi() - Dcand_truth.phi()
                                   
                                    if Isrefl:
                                        print("matched candidate is reflection")
                                        self.hjetvsDptvsm_refl.Fill(jets_det.pt(), Dcand_det.pt(), Dcand_det.m())
                                        self.fill_EEC_jet_histograms(jets_truth, Constituents_truth, Dcand_truth, jetR,  self.observable_list[8], False)

                                    elif Isfd:
                                        print("matched candidate is feeddown")

                                        self.hjetvsDpt_fd_det.Fill(jets_det.pt(), Dcand_det.pt())
                                        self.hTrackEtaPhi_nonprompt.Fill(etadiff, phidiff)
                                        self.hDeltaR_nonprompt.Fill(Dcand_det.delta_R(Dcand_truth))
                                      
                                        x_array = array( 'd', ( jets_det.perp(), jets_truth.perp(), Dcand_det.perp(), Dcand_truth.perp()))
                                        getattr(self, 'hsparse_R{}_{}'.format(jetR, self.observable_RM_list[1])).Fill(x_array)

                                        self.fill_EEC_jet_histograms(jets_truth, Constituents_truth, Dcand_truth, jetR,  self.observable_list[1], False)
                                        self.fill_EEC_jet_histograms(jets_det, Constituents_det, Dcand_det, jetR,  self.observable_list[3], False)

                                        self.fill_EEC_jet_histograms(jets_truth, Constituents_truth, Dcand_truth, jetR,  self.observable_list[5], True)
                                        self.fill_EEC_jet_histograms(jets_det, Constituents_det, Dcand_det, jetR,  self.observable_list[7], True)

                                    elif Isprompt:
                                        print("matched candidate is prompt")

                                        self.hjetvsDpt_prompt_det.Fill(jets_det.pt(), Dcand_det.pt())
                                        self.hTrackEtaPhi_prompt.Fill(etadiff, phidiff)
                                        self.hDeltaR_prompt.Fill(Dcand_det.delta_R(Dcand_truth))

                                        x_array = array( 'd', ( jets_det.perp(), jets_truth.perp(), Dcand_det.perp(), Dcand_truth.perp()))
                                        getattr(self, 'hsparse_R{}_{}'.format(jetR, self.observable_RM_list[0])).Fill(x_array)


                                        self.fill_EEC_jet_histograms(jets_truth, Constituents_truth, Dcand_truth, jetR,  self.observable_list[0], False)
                                        self.fill_EEC_jet_histograms(jets_det, Constituents_det, Dcand_det, jetR,  self.observable_list[2], False)


                                        self.fill_EEC_jet_histograms(jets_truth, Constituents_truth, Dcand_truth, jetR,  self.observable_list[4], True)
                                        self.fill_EEC_jet_histograms(jets_det, Constituents_det, Dcand_det, jetR,  self.observable_list[6], True)
                                    #also have smal fraction of cases where prompt, fd and ref are all zero





    def analyze_jets(self, fj_D0_particles, jet_def, IsDzero, IsConstituent):

        djmm = fjtools.DJetMatchMaker()
        _parts_and_ds = fj_D0_particles

        jet_def.analyze_event(_parts_and_ds)

        if len(jet_def.jets) < 1:
            return False

        jets = jet_def.jets_as_psj_vector()
        djets = djmm.filter_D0_jets(jets)


        if len(djets) > 0:
            j = djets[0]
            dcand = djmm.get_Dcand_in_jet(j)

            if IsDzero:
                return dcand[0]

            elif IsConstituent:
                return j.constituents()

            else:
                return j

        if len(djets) > 1:
            perror("more than one jet per D candidate?")
            return False


    def fill_EEC_jet_histograms(self, jet, Constituents, dcand, jetR, obs_label, IsPair):
        constituents = fj.sorted_by_pt(Constituents)
        c_select = fj.vectorPJ()
        
        if jet.pt() < 15.0:
            trk_thrd = self.trk_thrd[0]
        else:
            trk_thrd = self.trk_thrd[1]

        for c in constituents:
            if c.pt() < trk_thrd:
                break
            c_select.append(c) # NB: use the break statement since constituents are already sorted

        if self.ENC_pair_cut:
            dphi_cut = -9999 # means no dphi cut
            deta_cut = 0.008
        else:
            dphi_cut = -9999
            deta_cut = -9999
 
        hname = 'hsparse_R{}_{}'
        new_corr = ecorrel.CorrelatorBuilder(c_select, jet.perp(), 2, 1, dphi_cut, deta_cut)


        for index in range(new_corr.correlator(2).rs().size()):
            x_EEC_array = array( 'd', ( jet.perp(), dcand.perp(), dcand.m(), new_corr.correlator(2).rs()[index]))
            
            if IsPair:
                getattr(self, 'hsparse_R{}_{}'.format(jetR, obs_label)).Fill(x_EEC_array)
            else:
                getattr(self, 'hsparse_R{}_{}'.format(jetR, obs_label)).Fill(x_EEC_array, new_corr.correlator(2).weights()[index])
        
        
 

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
  

