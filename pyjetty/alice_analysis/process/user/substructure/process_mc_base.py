#!/usr/bin/env python3

"""
Base class to read a ROOT TTree of track information
and do jet-finding, and save basic histograms.
  
To use this class, the following should be done:

  - Implement a user analysis class inheriting from this one, such as in user/james/process_mc_XX.py
    You should implement the following functions:
      - initialize_user_output_objects_R()
      - fill_observable_histograms()
      - fill_matched_jet_histograms()
    
  - You should include the following histograms:
      - Response matrix: hResponse_JetPt_[obs]_R[R]_[subobs]_[grooming setting]
      - Residual distribution: hResidual_JetPt_[obs]_R[R]_[subobs]_[grooming setting]

  - You also should modify observable-specific functions at the top of common_utils.py
  
Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import time

# Data analysis and plotting
import pandas
import numpy as np
from array import *
import ROOT
import yaml
import random
import math

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjtools
import pythiafjext

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io
from pyjetty.alice_analysis.process.base import process_io_emb
from pyjetty.alice_analysis.process.base import process_base
from pyjetty.alice_analysis.process.base import thermal_generator
from pyjetty.alice_analysis.process.base import jet_info
from pyjetty.mputils.csubtractor import CEventSubtractor

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class ProcessMCBase(process_base.ProcessBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', event_start_offset=0, dstar=0, debug_level=0, **kwargs):
    # Initialize base class
    super(ProcessMCBase, self).__init__(input_file, config_file, output_dir, event_start_offset, dstar, debug_level, **kwargs)  
    
    # Initialize configuration
    self.initialize_config()
    
    # find pt_hat for set of events in input_file, assumes all events in input_file are in the same pt_hat bin
    if self.mcprod == True:
      slashes_from_end = 4
    else:
      slashes_from_end = 3
    self.pt_hat_bin = int(input_file.split('/')[len(input_file.split('/')) - slashes_from_end]) # depends on exact format of input_file name
    
    # need to specify comp_system and generator in config file!
    file_basepath = ''
    if self.compsystem == 'perlmutter':
      file_basepath = '/global/cfs/cdirs/alice/alicepro/hiccup'

    if self.generator == 'pythia':
      if self.mcprod:
        with open("{}/rstorage/alice/data/LHC18b8_charge/scaleFactors.yaml".format(file_basepath), 'r') as stream:
          pt_hat_yaml = yaml.safe_load(stream)
      else: #pythia fastsim
        print("FILE!", "{}/rstorage/generators/pythia_alice/tree_fastsim/scaleFactors.yaml".format(file_basepath))
        with open("{}/rstorage/generators/pythia_alice/tree_fastsim/scaleFactors.yaml".format(file_basepath), 'r') as stream:
          pt_hat_yaml = yaml.safe_load(stream)
    elif self.generator == 'herwig' and self.mcprod == False: #no anchored mc for herwig
      with open("{}/rstorage/generators/herwig_alice/tree_fastsim/scaleFactors.yaml".format(file_basepath), 'r') as stream:
        pt_hat_yaml = yaml.safe_load(stream)

    print("DEBUGGING", self.compsystem, "AND", self.generator, "AND", self.mcprod)

    print("FILE OUTPUT!", pt_hat_yaml)
    self.pt_hat = pt_hat_yaml[self.pt_hat_bin]
    print("pt hat bin : " + str(self.pt_hat_bin))
    print("pt hat weight : " + str(self.pt_hat))

    
  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_config(self):
    
    # Call base class initialization
    process_base.ProcessBase.initialize_config(self)
    
    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)
      
    self.fast_simulation = config['fast_simulation']
    if self.fast_simulation == True:
      if 'ENC_fastsim' in config:
          self.ENC_fastsim = config['ENC_fastsim']
      else:
          self.ENC_fastsim = False
    else: # if not fast simulation, set ENC_fastsim flag to False
      self.ENC_fastsim = False  
    if self.ENC_fastsim == True:
      self.pair_eff_file = config['pair_eff_file'] # load pair efficiency input for fastsim

    if 'mc_prod' in config:
      self.mcprod = config['mc_prod']
    else:
      self.mcprod = False
    
    if 'ENC_pair_cut' in config:
        self.ENC_pair_cut = config['ENC_pair_cut']
    else:
        self.ENC_pair_cut = False
    if 'ENC_pair_like' in config:
        self.ENC_pair_like = config['ENC_pair_like']
    else:
        self.ENC_pair_like = False
    if 'ENC_pair_unlike' in config:
        self.ENC_pair_unlike = config['ENC_pair_unlike']
    else:
        self.ENC_pair_unlike = False
    if 'jetscape' in config:
        self.jetscape = config['jetscape']
    else:
        self.jetscape = False
    if 'event_plane_angle' in config:
      self.event_plane_range = config['event_plane_angle']
    else:
      self.event_plane_range = None
    if 'matching_systematic' in config:
      self.matching_systematic = config['matching_systematic']
    else:
      self.matching_systematic = False
    if 'study_D0' in config:
      self.use_D0_info = config['study_D0']
    else:
      self.use_D0_info = False
    self.dry_run = config['dry_run']
    self.skip_deltapt_RC_histograms = True
    self.fill_RM_histograms = True

    if 'comp_system' in config:
      self.compsystem = config['comp_system']
    else:
      self.compsystem = '' #'perlmutter'

    if 'generator' in config:
      self.generator = config['generator']
    else:
      self.generator = '' #'pythia'

    if 'leadingtrack_pt_cut' in config:
      self.leading_parton_pt_cut = config['leadingtrack_pt_cut']
    else:
      self.leading_parton_pt_cut = 0.
    
    self.jet_matching_distance = config['jet_matching_distance']
    self.reject_tracks_fraction = config['reject_tracks_fraction']
    if 'mc_fraction_threshold' in config:
      self.mc_fraction_threshold = config['mc_fraction_threshold']
    if 'do_rho_subtraction' in config:
      self.do_rho_subtraction = config['do_rho_subtraction']
    else:
      self.do_rho_subtraction = False
    if 'do_jetcone' in config:
      self.do_jetcone = config['do_jetcone']
    else:
      self.do_jetcone = False
    if self.do_jetcone and 'jetcone_R_list' in config:
      self.jetcone_R_list = config['jetcone_R_list']
    else:
      self.jetcone_R_list = [0.4] # NB: set default value to 0.4
    if 'leading_pt' in config:
        self.leading_pt = config['leading_pt']
    else:
        self.leading_pt = -1 # negative means no leading track cut
    
    if self.do_constituent_subtraction:
        self.is_pp = False
        self.emb_file_list = config['emb_file_list']
        self.main_R_max = config['constituent_subtractor']['main_R_max']
    else:
        self.is_pp = True
        
    if 'thermal_model' in config:
      self.thermal_model = True
      beta = config['thermal_model']['beta']
      N_avg = config['thermal_model']['N_avg']
      sigma_N = config['thermal_model']['sigma_N']
      self.thermal_generator = thermal_generator.ThermalGenerator(N_avg, sigma_N, beta)
    else:
      self.thermal_model = False

    # Create dictionaries to store grooming settings and observable settings for each observable
    # Each dictionary entry stores a list of subconfiguration parameters
    #   The observable list stores the observable setting, e.g. subjetR
    #   The grooming list stores a list of grooming settings {'sd': [zcut, beta]} or {'dg': [a]}
    self.observable_list = config['process_observables']
    self.obs_settings = {}
    self.obs_grooming_settings = {}
    for observable in self.observable_list:
    
      obs_config_dict = config[observable]
      obs_config_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]
      
      obs_subconfig_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]
      self.obs_settings[observable] = self.utils.obs_settings(observable, obs_config_dict, obs_subconfig_list)
      self.obs_grooming_settings[observable] = self.utils.grooming_settings(obs_config_dict)
      
    # Construct set of unique grooming settings
    self.grooming_settings = []
    lists_grooming = [self.obs_grooming_settings[obs] for obs in self.observable_list]
    for observable in lists_grooming:
      for setting in observable:
        if setting not in self.grooming_settings and setting != None:
          self.grooming_settings.append(setting)


  # Check if two four momenta are (approximately) equal
  def four_mom_equal(self, a, b, tol=0.1):
    return (abs(a.px() - b.px()) < tol and
            abs(a.py() - b.py()) < tol and
            abs(a.pz() - b.pz()) < tol)

    # return (abs(a.px() - b.px()) < tol and
    #           abs(a.py() - b.py()) < tol and
    #           abs(a.pz() - b.pz()) < tol and
    #           abs(a.e()  - b.e())  < tol)


  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def process_mc(self):
    
    self.start_time = time.time()
    
    # ------------------------------------------------------------------------
    
    # Use IO helper class to convert detector-level ROOT TTree into
    # a SeriesGroupBy object of fastjet particles per event
    print('--- {} seconds ---'.format(time.time() - self.start_time))
    if self.fast_simulation:
      tree_dir = ''
    else:
      tree_dir = 'PWGHF_TreeCreator'

    print("Starting to load data for detector level tracks")
    io_det = process_io.ProcessIO(input_file=self.input_file, tree_dir=tree_dir,
                                  track_tree_name='tree_Particle', use_ev_id_ext=False, use_D0_info=self.use_D0_info,
                                  is_jetscape=self.jetscape, event_plane_range=self.event_plane_range, is_ENC=self.ENC_fastsim, is_det_level=True, is_mcprod=self.mcprod)
    df_fjparticles_det = io_det.load_data(m=self.m, reject_tracks_fraction=self.reject_tracks_fraction)
    self.nEvents_det = len(df_fjparticles_det.index)
    self.nTracks_det = len(io_det.track_df.index)
    print('--- {} seconds ---'.format(time.time() - self.start_time))
    
    # If jetscape, store also the negative status particles (holes)
    if self.jetscape:
      io_det_holes = process_io.ProcessIO(input_file=self.input_file, tree_dir=tree_dir,
                                          track_tree_name='tree_Particle', use_ev_id_ext=False,
                                          is_jetscape=self.jetscape, holes=True,
                                          event_plane_range=self.event_plane_range)
      df_fjparticles_det_holes = io_det_holes.load_data(m=self.m, reject_tracks_fraction=self.reject_tracks_fraction)
      self.nEvents_det_holes = len(df_fjparticles_det_holes.index)
      self.nTracks_det_holes = len(io_det_holes.track_df.index)
      print('--- {} seconds ---'.format(time.time() - self.start_time))
    
    
    # ------------------------------------------------------------------------

    # Use IO helper class to convert truth-level ROOT TTree into
    # a SeriesGroupBy object of fastjet particles per event
    
    print("Starting to load data for truth level tracks", self.use_D0_info)
    io_truth = process_io.ProcessIO(input_file=self.input_file, tree_dir=tree_dir,
                                    track_tree_name='tree_Particle_gen', use_ev_id_ext=False, use_D0_info=self.use_D0_info,
                                    is_jetscape=self.jetscape, event_plane_range=self.event_plane_range, is_ENC=self.ENC_fastsim, is_det_level=False, is_mcprod=self.mcprod)
    df_fjparticles_truth = io_truth.load_data(m=self.m) # no dropping of tracks at truth level (important for the det-truth association because the index of the truth particle is used)
    self.nEvents_truth = len(df_fjparticles_truth.index)
    self.nTracks_truth = len(io_truth.track_df.index)
    print('--- {} seconds ---'.format(time.time() - self.start_time))

    # print('Input truth Data Frame',df_fjparticles_truth)
    
    # If jetscape, store also the negative status particles (holes)
    if self.jetscape:
      io_truth_holes = process_io.ProcessIO(input_file=self.input_file, tree_dir=tree_dir,
                                            track_tree_name='tree_Particle_gen', use_ev_id_ext=False,
                                            is_jetscape=self.jetscape, holes=True,
                                            event_plane_range=self.event_plane_range)
      df_fjparticles_truth_holes = io_truth_holes.load_data(m=self.m, reject_tracks_fraction=self.reject_tracks_fraction)
      self.nEvents_truth_holes = len(df_fjparticles_truth_holes.index)
      self.nTracks_truth_holes = len(io_truth_holes.track_df.index)
      print('--- {} seconds ---'.format(time.time() - self.start_time))

    #if D0, replace all kaon/pion pairs with the D0 here! //TODO: save d0 rapidity!
    # start by getting the D0s
    if self.use_D0_info:
      # # get the det level D0s
      # print("Starting to load data for detector level D0s")
      # io_D0_det = process_io.ProcessIO(input_file=self.input_file, tree_dir=tree_dir,
      #                               track_tree_name='tree_D0', use_ev_id_ext=False, use_D0_info=True, using_dstar=self.dstar,
      #                               is_jetscape=self.jetscape, event_plane_range=self.event_plane_range, is_ENC=self.ENC_fastsim, is_det_level=True) #for D0 herwig case, mcprod is set to false
      # df_D0particles_det = io_D0_det.load_data(m=self.m, reject_tracks_fraction=self.reject_tracks_fraction) # for D0 herwig case, reject_tracks_fraction is set to 0
      # self.nEvents_det = len(df_D0particles_det.index)
      # self.nTracks_det = len(io_D0_det.track_df.index)
      # print('--- {} seconds ---'.format(time.time() - self.start_time))

      # get the truth level D0s
      print("Starting to load data for truth level D0s")
      io_D0_truth = process_io.ProcessIO(input_file=self.input_file, tree_dir=tree_dir,
                                    track_tree_name='tree_D0_gen', use_ev_id_ext=False, use_D0_info=True, using_dstar=self.dstar,
                                    is_jetscape=self.jetscape, event_plane_range=self.event_plane_range, is_ENC=self.ENC_fastsim, is_det_level=False)
      df_D0particles_truth = io_D0_truth.load_data(m=self.m) # no dropping of tracks at truth level (important for the det-truth association because the index of the truth particle is used)
      self.nEvents_truth = len(df_D0particles_truth.index)
      self.nTracks_truth = len(io_D0_truth.track_df.index)
      print('--- {} seconds ---'.format(time.time() - self.start_time))


      print("COLS!", df_D0particles_truth.columns)

      # remove kaons and pions with D0 mother
      D0_PIDs = [421, -421]
      print("LOOK HERE", df_fjparticles_truth.columns)
      print("length", len(df_fjparticles_truth))
      # print("keep for debug, len of entry 21 with a D0", len(df_fjparticles_truth['MotherPID'].values[21]))


      #print("index", df_fjparticles_truth.index)    #gets the row names!! - these rows are indexed by MultiIndex!
      run_num = df_fjparticles_truth.index[0][0]

      # loop over the particles in the DET level and reconstruct any D0s -- this does not separate the D*'s
      mother_pids_det = df_fjparticles_det['MotherPID'].values
      track_pids_det = df_fjparticles_det['ParticlePID'].values
      track_mcindices_det = df_fjparticles_det['ParticleMCIndex'].values
      track_4vec_det = df_fjparticles_det['fj_particle'].values
      track_ev_id_det = df_fjparticles_det['ev_id'].values

      track_4vec_truth = df_fjparticles_truth['fj_particle'].values

      d0_ev_id_truth = df_D0particles_truth['ev_id'].values
      d0_4vec_truth = df_D0particles_truth['fj_particle'].values
      d0_pid_truth = df_D0particles_truth['ParticlePID'].values
      d0_mid_truth = df_D0particles_truth['MotherPID'].values
      d0_rap_truth = df_D0particles_truth['ParticleRapidity'].values

      #d0_ev_id_truth is in the form: [array([6]) array([25]) array([37]) ... array([3920,3920]) ... array([19942]) array([19959]) where each index gives you an array(a,b) for one event number.
      # which means we want it in the form [6 25 etc] -- needs to be flattened
      d0_ev_id_truth_flat = [item for sublist in d0_ev_id_truth for item in sublist]
      d0_4vec_truth_flat = [item for sublist in d0_4vec_truth for item in sublist]
      d0_pid_truth_flat = [item for sublist in d0_pid_truth for item in sublist]
      d0_mid_truth_flat = [item for sublist in d0_mid_truth for item in sublist]
      d0_rap_truth_flat = [item for sublist in d0_rap_truth for item in sublist]

      d0_ev_id_det = df_D0particles_det['ev_id'].values #todo: flatten this here too???
      d0_mcindices_det = df_D0particles_det['ParticleMCIndex'].values

      # Start event loop here! -- detector level edits
      for ind_ev, (iev_arr, event_mpids) in enumerate(zip(track_ev_id_det, mother_pids_det)):
        iev = iev_arr[0]
        ievent_adj = self.event_start_offset+iev # this is the event number (used in case any events are skipped)
        ind_ev_adj = self.event_start_offset+ind_ev # this is the enumerate -- indexing for accessing arrays
        # print("self.event_start_offset", self.event_start_offset)
        # print("iev", iev, "ievent_adj", ievent_adj, "ind_ev_adj", ind_ev_adj, "iev_arr", iev_arr)
        dau_inds = []

        # ind_ev_adj is the indexer for det level events -- gives the actual event ID! and skips any events that don't have det level particles
        # ievent_adj is the event counter for truth level events! (just counts in order)
        if D0_PIDs[0] in event_mpids or D0_PIDs[1] in event_mpids: #found a particle with a D0 mother
          print()
          print("iev", iev, "ievent_adj", ievent_adj, "ind_ev_adj", ind_ev_adj) 
          # print("EVEBT MPID", event_mpids)
          # now look for the other daughter of the D0
          remaining_dau_count = sum(1 for mpid in event_mpids if abs(mpid) == 421) #the number of particles that have a D0 mother in the det-level track list
          # print(remaining_dau_count)

          if remaining_dau_count <= 1: #if there are 1 or less particles with D0 mothers left, then there are no more D0s to reconstruct
            continue
          
          possible_daughter_indices = [ind for ind, mpid in enumerate(event_mpids) if abs(mpid) == 421]
          # print("possible_daughter_indices", possible_daughter_indices)
          rev_possible_daughter_indices = list(reversed(possible_daughter_indices))
          possible_daughter_mother_ids = [event_mpids[ind] for ind in rev_possible_daughter_indices]
          # print("possible_daughter_mother_ids", possible_daughter_mother_ids)
          # print("track_pids_det[ind_ev_adj]", track_pids_det[ind_ev_adj])

          possible_daughter_pids = [track_pids_det[ind_ev_adj][ind] for ind in rev_possible_daughter_indices]
          # print("possible_daughter_pids", possible_daughter_pids)
          possible_daughter_mcindices = [track_mcindices_det[ind_ev_adj][ind] for ind in rev_possible_daughter_indices]
          possible_daughter_4vec_det = [track_4vec_det[ind_ev_adj][ind] for ind in rev_possible_daughter_indices]

          # print("REVERESED poss dau ind", rev_possible_daughter_indices)
          # print("possible_daughter_mother_ids", possible_daughter_mother_ids)
          # print("possible_daughter_mcindices", possible_daughter_mcindices)

          possible_daughter_4vec_truth = [track_4vec_truth[ievent_adj][int(ind)] for ind in possible_daughter_mcindices] #rev_possible_daughter_indices]                    
          
          list_of_D0_indices_in_event_truth = [ind for ind, evid in enumerate(d0_ev_id_truth_flat) if evid == ievent_adj]
          rev_list_of_D0_indices_in_event_truth = list(reversed(list_of_D0_indices_in_event_truth))
          list_of_D0s_4vec_truth = [d0_4vec_truth_flat[ind] for ind in rev_list_of_D0_indices_in_event_truth] #need 0'th index because it's the form of an array of size 1
          list_of_D0s_pid_truth = [d0_pid_truth_flat[ind] for ind in rev_list_of_D0_indices_in_event_truth]
          list_of_D0s_mid_truth = [d0_mid_truth_flat[ind] for ind in rev_list_of_D0_indices_in_event_truth]
          list_of_D0s_rap_truth = [d0_rap_truth_flat[ind] for ind in rev_list_of_D0_indices_in_event_truth]

          [print("index of D0 in D0 tree:", ind, "// and corresponding evid", evid) for ind, evid in enumerate(d0_ev_id_truth_flat) if evid == ievent_adj]

          list_of_D0_indices_in_event_det = [ind for ind, evid in enumerate(d0_ev_id_det) if evid[0] == ievent_adj]
          list_of_D0s_mcindices_det = [d0_mcindices_det[ind] for ind in list_of_D0_indices_in_event_det]
          # print("rev_list_of_D0_indices_in_event_truth", rev_list_of_D0_indices_in_event_truth)

          if ievent_adj < 1000:
            print("keep for debug, length before adjustments", len(df_fjparticles_det['fj_particle'].values[ind_ev_adj])) #, 
                # len(df_fjparticles_det['ParticlePID'].values[ind_ev_adj]), len(df_fjparticles_det['MotherPID'].values[ind_ev_adj]))
          

          # Loop through daughter particles - start matching and replacing with D0s!
          while remaining_dau_count > 1: #if there are 1 or less particles with D0 mothers left, then there are no more D0s to reconstruct
            # checks: that both particles come from D0 or D0-bar // that the particles are opp. signed k+pi // that the 4-mom of k+pi == D0 [but need to do this at truth level]
            # look at the first track in the list and compare to the others
            # print(ievent_adj, " // remaining DAUGHTER count:", remaining_dau_count)

            match_made = False # select this when trying to break out of loop

            for j in range(1, len(rev_possible_daughter_indices)):

              if possible_daughter_mother_ids[0] == possible_daughter_mother_ids[j]:

                if abs(possible_daughter_pids[0] + possible_daughter_pids[j]) == 110:

                  # do one more check here -- check if truth level k+pi == D0 in 4-mom
                  # make sure to get the mcid from detector level to access the correct truth level particle
                  for k,D0_4vec in enumerate(list_of_D0s_4vec_truth):
                    print("DAU 1: ", possible_daughter_4vec_truth[0], possible_daughter_4vec_truth[0].pt(), possible_daughter_4vec_truth[0].eta(), possible_daughter_4vec_truth[0].phi())
                    print("DAU 2: ", possible_daughter_4vec_truth[int(j)], possible_daughter_4vec_truth[int(j)].pt(), possible_daughter_4vec_truth[int(j)].eta(), possible_daughter_4vec_truth[int(j)].phi())
                    print("D0: ", D0_4vec, D0_4vec.pt(), D0_4vec.eta(), D0_4vec.phi())
                    if ( self.four_mom_equal(possible_daughter_4vec_truth[0] + possible_daughter_4vec_truth[int(j)], D0_4vec) ):

                      # then do replacement
                      for col in df_fjparticles_det.loc[(run_num, ievent_adj)].index: #lists the columns 
                        if col == 'fj_particle':
                          # remove kaon and pion
                          df_fjparticles_det.loc[(run_num, ievent_adj)][col] = pythiafjext.removeByIndex(df_fjparticles_det.loc[(run_num, ievent_adj)][col], rev_possible_daughter_indices[0]) 
                          df_fjparticles_det.loc[(run_num, ievent_adj)][col] = pythiafjext.removeByIndex(df_fjparticles_det.loc[(run_num, ievent_adj)][col], rev_possible_daughter_indices[j]) 
                        else:
                          # print("col", col, df_fjparticles_det.loc[(run_num, ievent_adj),col])
                          df_fjparticles_det.loc[(run_num, ievent_adj)][col] = np.delete(df_fjparticles_det.loc[(run_num, ievent_adj)][col], rev_possible_daughter_indices[0])
                          # print("col", col, df_fjparticles_det.loc[(run_num, ievent_adj),col])
                          df_fjparticles_det.loc[(run_num, ievent_adj)][col] = np.delete(df_fjparticles_det.loc[(run_num, ievent_adj)][col], rev_possible_daughter_indices[j])
                      if ievent_adj < 1000:
                        print("keep for debug, length after kpi removal", len(df_fjparticles_det['fj_particle'].values[ind_ev_adj])) #, 
                            # len(df_fjparticles_det['ParticlePID'].values[ind_ev_adj]), len(df_fjparticles_det['MotherPID'].values[ind_ev_adj]))

                      # and add the D0

                      col = 'fj_particle'
                      df_fjparticles_det.loc[(run_num, ievent_adj)][col] = pythiafjext.addByIndex(df_fjparticles_det.loc[(run_num, ievent_adj)][col], 
                                                                             rev_possible_daughter_indices[j], possible_daughter_4vec_det[0]+possible_daughter_4vec_det[j])
                      col = 'ParticlePID'
                      df_fjparticles_det.loc[(run_num, ievent_adj),col] = np.insert(df_fjparticles_det.loc[(run_num, ievent_adj),col],
                                                                            rev_possible_daughter_indices[j], list_of_D0s_pid_truth[k])
                      # print("innnn into this array:", df_fjparticles_det.loc[(run_num, ievent_adj)][col])
                      
                      col = 'MotherPID'
                      # print("inserting", list_of_D0s_mid_truth[k], "into", col)
                      # print("into this array:", df_fjparticles_det.loc[(run_num, ievent_adj)][col])
                      df_fjparticles_det.loc[(run_num, ievent_adj)][col] = np.insert(df_fjparticles_det.loc[(run_num, ievent_adj)][col], 
                                                                             rev_possible_daughter_indices[j], list_of_D0s_mid_truth[k]) #-1) 
                      # print("into this array:", df_fjparticles_det.loc[(run_num, ievent_adj)][col])
                    
                      col = 'ParticleRapidity'
                      df_fjparticles_det.loc[(run_num, ievent_adj)][col] = np.insert(df_fjparticles_det.loc[(run_num, ievent_adj)][col], rev_possible_daughter_indices[j], list_of_D0s_rap_truth[k])
                      # print("just added d0 particles here what is happening?")
                      
                      # these are placeholder values right now. replace these indices properly in next loop
                      col = 'ParticleMCIndex'
                      # print("rev_list_of_D0_indices_in_event_truth[k]: ", rev_list_of_D0_indices_in_event_truth[k], " // list_of_D0s_mcindices_det:", list_of_D0s_mcindices_det)
                      # if ( len(list_of_D0s_mcindices_det) > 0 ): # and ( int(len(list_of_D0s_4vec_truth) - 1 - k) in list_of_D0s_mcindices_det[0] ):
                        # print("if passed if", list_of_D0s_mcindices_det[0], len(list_of_D0s_4vec_truth), k, len(list_of_D0s_4vec_truth) - 1 - k)
                      D0_mcid = int(len(list_of_D0s_4vec_truth) - 1 - k) #rev_list_of_D0_indices_in_event_truth[k]
                      print("D0 mcid is ", D0_mcid)

                        # print("inserting", type(d0_rap_truth[d0_event_counter][i_d0]), d0_rap_truth[d0_event_counter][i_d0], "of", d0_rap_truth[d0_event_counter], "into", col)
                      df_fjparticles_det.loc[(run_num, ievent_adj)][col] = np.insert(df_fjparticles_det.loc[(run_num, ievent_adj)][col], rev_possible_daughter_indices[j], D0_mcid)
                      
                      '''
                      if ( len(list_of_D0s_mcindices_det) > 0 ) and ( int(len(list_of_D0s_4vec_truth) - 1 - k) in list_of_D0s_mcindices_det[0] ):
                        print("if passed if", list_of_D0s_mcindices_det[0], len(list_of_D0s_4vec_truth), k, len(list_of_D0s_4vec_truth) - 1 - k)
                        D0_mcid = int(len(list_of_D0s_4vec_truth) - 1 - k) #rev_list_of_D0_indices_in_event_truth[k]
                        print("D0 mcid is ", D0_mcid)
                      else: # negative means the D0 was not in the detector level tree
                        D0_mcid = -1 * int(len(list_of_D0s_4vec_truth) - 1 - k)
                      '''
                      
                      if ievent_adj < 1000:
                        # print("keep for debug, length at end", len(df_fjparticles_det['fj_particle'].values[ind_ev_adj]), 
                        #     len(df_fjparticles_det['ParticlePID'].values[ind_ev_adj]), len(df_fjparticles_det['MotherPID'].values[ind_ev_adj]))
                        print("keep for debug, length at end", len(df_fjparticles_det['fj_particle'].values[ind_ev_adj]))


                      # clear the track arrays
                      rev_possible_daughter_indices.pop(j)
                      rev_possible_daughter_indices.pop(0)
                      possible_daughter_mother_ids.pop(j)
                      possible_daughter_mother_ids.pop(0)
                      possible_daughter_pids.pop(j) 
                      possible_daughter_pids.pop(0) 
                      possible_daughter_4vec_det.pop(j)
                      possible_daughter_4vec_det.pop(0)

                      possible_daughter_4vec_truth.pop(j)
                      possible_daughter_4vec_truth.pop(0)

                      # clear the D0 arrays
                      list_of_D0s_4vec_truth.pop(k)
                      list_of_D0_indices_in_event_truth.pop(k)
                      list_of_D0s_pid_truth.pop(k)
                      list_of_D0s_mid_truth.pop(k)
                      list_of_D0s_rap_truth.pop(k)
                                
                      match_made = True
                      remaining_dau_count -= 2
                      break
            
              
                  # print(ievent_adj, " // check: DAUGHTER count:", remaining_dau_count, "match_made", match_made)
                  if match_made:
                    break  # stop inner loop after modifying arr - get out of for loop of matching daughters, and go back to the while loop
                  

            
            # if no pairs made, the daughter particles stays and the needs to be removed from the list for the while loop (not from the actual dataframe) 
            # clear the track arrays
            if not match_made: #remaining_dau_count >= 1:
              rev_possible_daughter_indices.pop(0)
              possible_daughter_mother_ids.pop(0)
              possible_daughter_pids.pop(0) 
              possible_daughter_4vec_det.pop(0)

              possible_daughter_4vec_truth.pop(0) 
              remaining_dau_count -= 1
      
      # print("CHECK: len of df_fjparticles_det", len(df_fjparticles_det))


      
      # loop over the D0s in the TRUTH level and reconstruct
      mother_pids = df_fjparticles_truth['MotherPID'].values

      d0_event_counter = 0
      num_d0s_in_event = 0
      self.alld0counter_truth = 0
      self.d0nodstar_counter_truth = 0
      self.alld0counter_det = 0
      self.d0nodstar_counter_det = 0
      self.alld0counter_truthmatched = 0
      self.d0nodstar_counter_truthmatched = 0
      self.alld0counter_detmatched = 0
      self.d0nodstar_counter_detmatched = 0

      # Start event loop here! -- truth level edits
      for iev, event_mpids in enumerate(mother_pids):
        ievent_adj = self.event_start_offset+iev
        # if iev == 500:
        #   break
        dau_inds = []
        if D0_PIDs[0] in event_mpids or D0_PIDs[1] in event_mpids:

          print()
          print("iev", iev, "ievent_adj", ievent_adj)

          if iev < 1000:
            print("keep for debug, length to start", len(df_fjparticles_truth['fj_particle'].values[iev]), 
                len(df_fjparticles_truth['ParticlePID'].values[iev]), len(df_fjparticles_truth['MotherPID'].values[iev]))

          dau_inds = [i for i, mpid in enumerate(event_mpids) if abs(mpid) == 421] # List comprehension
          rev_dau_inds = list(reversed(dau_inds))
          # print("debug index", iev, "with D0 daughter indices", dau_inds)

          num_d0s_in_event = len(d0_ev_id_truth[d0_event_counter]) # get rid of this

          list_of_D0_indices_in_event_truth = [ind for ind, evid in enumerate(d0_ev_id_truth) if evid[0] == ievent_adj]
          # list_of_D0_indices_in_event_truth = [ind for ind, evid in enumerate(d0_ev_id_truth_flat) if evid == ievent_adj]
          rev_list_of_D0_indices_in_event_truth = list(reversed(list_of_D0_indices_in_event_truth))
          list_of_D0s_4vec_truth = [d0_4vec_truth[ind][0] for ind in rev_list_of_D0_indices_in_event_truth] 
          # list_of_D0s_4vec_truth = [d0_4vec_truth_flat[ind][0] for ind in rev_list_of_D0_indices_in_event_truth] 

        
          for i_dau,dau_ind in enumerate(rev_dau_inds):
            
            # DELETING the rows that have the kaon and pion
            for col in df_fjparticles_truth.loc[(run_num, ievent_adj)].index: #lists the columns
              
              if col == 'fj_particle':
                # print("dau ind", dau_ind, type(dau_ind))
                # print("Check: what is being removed", i_dau, ":", df_fjparticles_truth.loc[(run_num, ievent_adj)]['ParticlePID'][dau_ind], df_fjparticles_truth.loc[(run_num, ievent_adj)]['fj_particle'][dau_ind].pt(), df_fjparticles_truth.loc[(run_num, ievent_adj)]['fj_particle'][dau_ind].eta())
                df_fjparticles_truth.loc[(run_num, ievent_adj)][col] = pythiafjext.removeByIndex(df_fjparticles_truth.loc[(run_num, ievent_adj)][col], dau_ind) 
              else:
                df_fjparticles_truth.loc[(run_num, ievent_adj)][col] = np.delete(df_fjparticles_truth.loc[(run_num, ievent_adj)][col], dau_ind)

            if iev < 1000:
              print("keep for debug, length after kpi removal", len(df_fjparticles_truth['fj_particle'].values[iev]), 
                  len(df_fjparticles_truth['ParticlePID'].values[iev]), len(df_fjparticles_truth['MotherPID'].values[iev]))

            # and ADDING row for the D0
            if i_dau%2 == 1:
              
              D0_index = int((len(dau_inds)-1-i_dau)/2) #int(i_dau/2) #rev_list_of_D0_indices_in_event_truth[int(i_dau/2)]
              # print(" i_dau:", i_dau, "int(i_dau/2):", int(i_dau/2), "rev_list_of_D0_indices_in_event_truth:", rev_list_of_D0_indices_in_event_truth)
              # print(" rev_dau_inds:", rev_dau_inds)
              print(" D0_index:", D0_index, "d0_event_counter:", d0_event_counter)

              col = 'fj_particle'
              # print("inserting", d0_4vec_truth[d0_event_counter][D0_index], "of", len(d0_4vec_truth[d0_event_counter]), " items into", col)
              df_fjparticles_truth.loc[(run_num, ievent_adj)][col] = pythiafjext.addByIndex(df_fjparticles_truth.loc[(run_num, ievent_adj)][col], 
                                                                                          rev_dau_inds[i_dau], d0_4vec_truth[d0_event_counter][D0_index])
                      
              col = 'ParticlePID'
              # print("inserting", d0_pid_truth[d0_event_counter][D0_index], "of", d0_pid_truth[d0_event_counter], "into", col)
              df_fjparticles_truth.loc[(run_num, ievent_adj)][col] = np.insert(df_fjparticles_truth.loc[(run_num, ievent_adj)][col],
                                                                            rev_dau_inds[i_dau], d0_pid_truth[d0_event_counter][D0_index])
              col = 'MotherPID'
              # print("inserting MID", d0_mid_truth[d0_event_counter][D0_index], "into", col)
              df_fjparticles_truth.loc[(run_num, ievent_adj)][col] = np.insert(df_fjparticles_truth.loc[(run_num, ievent_adj)][col], 
                                                                              rev_dau_inds[i_dau], d0_mid_truth[d0_event_counter][D0_index]) #-1)

              col = 'ParticleRapidity'
              # print("inserting", d0_rap_truth[d0_event_counter][D0_index], "of", d0_rap_truth[d0_event_counter], "into", col)
              df_fjparticles_truth.loc[(run_num, ievent_adj)][col] = np.insert(df_fjparticles_truth.loc[(run_num, ievent_adj)][col], rev_dau_inds[i_dau], d0_rap_truth[d0_event_counter][D0_index])


              # D0_part_that_was_added = df_fjparticles_truth.loc[(run_num, ievent_adj)]['fj_particle'][rev_dau_inds[i_dau]]
              # print("Check: D0 being added to spot:", rev_dau_inds[i_dau], " pt:", D0_part_that_was_added.pt(), "eta:", D0_part_that_was_added.eta(), "pid:", d0_pid_truth[d0_event_counter][D0_index], "mid:", d0_mid_truth[d0_event_counter][D0_index], "rap:", d0_rap_truth[d0_event_counter][D0_index])


              # Adjusting the ParticleMCIndex of the det-level tree since the gen-level df is changing here with replacement of kpi->D0
              # note: all mcids point to the index in truth particle array for each event (including for the D0s)
              if (run_num, ievent_adj) in df_fjparticles_det.index: # check if det-level particles exist in this event first!
                # print("lenth det, before mc index adjustment:", len(df_fjparticles_det.loc[(run_num, ievent_adj)]['ParticleMCIndex']),df_fjparticles_det.loc[(run_num, ievent_adj)]['ParticleMCIndex'])
                # print(type(df_fjparticles_det.loc[(run_num, ievent_adj)]['ParticleMCIndex']))
                # print(df_fjparticles_det.loc[(run_num, ievent_adj)]['ParticleMCIndex'].ndim)
                if df_fjparticles_det.loc[(run_num, ievent_adj)]['ParticleMCIndex'].ndim == 0: # this is a 0-d array but it looks like its saved as a "scalar" and it's the D0!
                  # print("this is a scalar!!!!")
                  df_fjparticles_det.loc[(run_num, ievent_adj)]['ParticleMCIndex'] = np.array([rev_dau_inds[i_dau]])
                  # print(df_fjparticles_det.loc[(run_num, ievent_adj)]['ParticleMCIndex'].ndim)
                  # print(df_fjparticles_det.loc[(run_num, ievent_adj)]['ParticleMCIndex'])
                else:
                  for i_item,mcid in enumerate(df_fjparticles_det.loc[(run_num, ievent_adj)]['ParticleMCIndex']): #loop over items in array of ParticleMCIndex
                    if abs(df_fjparticles_det.loc[(run_num, ievent_adj)]['ParticlePID'][i_item]) == 421 and D0_index == mcid:
                      df_fjparticles_det.loc[(run_num, ievent_adj)]['ParticleMCIndex'][i_item] = rev_dau_inds[i_dau]
                      # print("D0 found at", i_item)
                    if mcid > rev_dau_inds[i_dau]:
                      df_fjparticles_det.loc[(run_num, ievent_adj)]['ParticleMCIndex'][i_item] = mcid - 1

                  # if abs(mcid) > rev_dau_inds[i_dau]: #rev_possible_daughter_indices[0]:
                  #   if mcid < 0: #this is a negative D0 mc index to point to not being there at det-level... need to add 1 instead of subtract
                  #     df_fjparticles_det.loc[(run_num, ievent_adj)]['ParticleMCIndex'][i_item] = mcid + 1
                  #   else:
                  #     df_fjparticles_det.loc[(run_num, ievent_adj)]['ParticleMCIndex'][i_item] = mcid - 1 #note: mcid points to the index in truth particle array for each event
                  # else:
                  #   print("D0 found at", i_item) #here mcid also points to the index in D0 tree amongst each event
                # print("Check: after mc index adjustment        :", len(df_fjparticles_det.loc[(run_num, ievent_adj)]['ParticleMCIndex']),df_fjparticles_det.loc[(run_num, ievent_adj)]['ParticleMCIndex'])

              if iev < 1000:
                print("keep for debug, length after D0 removal", len(df_fjparticles_truth['fj_particle'].values[iev]), 
                    len(df_fjparticles_truth['ParticlePID'].values[iev]), len(df_fjparticles_truth['MotherPID'].values[iev]))
              
          d0_event_counter+=1
              
          
        
    # print("COLS!", df_D0particles_truth.columns) #COLS! Index(['fj_particle', 'ParticlePID', 'MotherPID'], dtype='object')


    
    # ------------------------------------------------------------------------

    # Now merge the two SeriesGroupBy to create a groupby df with [ev_id, run_number, fj_1, fj_2]
    # (Need a structure such that we can iterate event-by-event through both fj_1, fj_2 simultaneously)
    # In the case of jetscape, we merge also the hole collections fj_3, fj_4
    print('Merge det-level and truth-level into a single dataframe grouped by event...')
    print("df fj particles det")
    print(df_fjparticles_det)
    print("df fj particles truth")
    print(df_fjparticles_truth)

    if self.jetscape:
      self.df_fjparticles = pandas.concat([df_fjparticles_det, df_fjparticles_truth, df_fjparticles_det_holes, df_fjparticles_truth_holes], axis=1)
      self.df_fjparticles.columns = ['fj_particles_det', 'fj_particles_truth', 'fj_particles_det_holes', 'fj_particles_truth_holes']
    elif self.ENC_fastsim:

      if self.use_D0_info:


        self.df_fjparticles = pandas.concat([df_fjparticles_truth, df_fjparticles_det, df_D0particles_truth], axis=1) #, df_D0particles_det], axis=1)
        self.df_fjparticles.columns = ['fj_particles_truth', 'ParticlePID_truth', 'ParticleRapidity_truth', 'MotherPID_truth', 
                                       'fj_particles_det', 'ev_id_det', 'ParticleMCIndex_det', 'ParticlePID_det', 'ParticleRapidity_det', 'MotherPID_det', 
                                       'fj_D0_truth', 'ev_id_D0_truth', 'D0Rapidity_truth', 'D0PID_truth', 'D0MotherPID_truth']
                                      #  'fj_D0_det', 'ev_id_D0_det', 'D0Rapidity_det', 'D0MCIndex_det', 'D0PID_det', 'D0MotherPID_det'] #, "ev_id_corr"]
        
        # Combine repeat columns and drop unnecessary
        # self.df_fjparticles["D0PID"] = self.df_fjparticles["D0PID_truth"].combine_first(self.df_fjparticles["D0PID_det"]) # combine these columns into one
        # self.df_fjparticles["D0MotherPID"] = self.df_fjparticles["D0MotherPID_truth"].combine_first(self.df_fjparticles["D0MotherPID_det"]) # combine these columns into one
        # self.df_fjparticles = self.df_fjparticles.drop(columns=["ev_id_det", "ev_id_D0_det", "ev_id_D0_truth", "D0PID_truth", "D0PID_det", "D0MotherPID_truth", "D0MotherPID_det"]) # drop these columns
        self.df_fjparticles = self.df_fjparticles.drop(columns=["ev_id_det", "ev_id_D0_truth"]) # drop these columns
        
        # By the end, the updated columns will be:
        # self.df_fjparticles.columns = ['fj_particles_truth', 'ParticlePID_truth', 'ParticleRapidity_truth', 'MotherPID_truth', 
        #                                'fj_particles_det', 'ParticleMCIndex_det', 'ParticlePID_det', 'ParticleRapidity_det', 'MotherPID_det', 
        #                                'fj_D0_truth', 'D0Rapidity_truth', 'D0PID_truth', 'D0MotherPID_truth'] 

      else:
        self.df_fjparticles = pandas.concat([df_fjparticles_det, df_fjparticles_truth], axis=1)
        self.df_fjparticles.columns = ['fj_particles_det', 'ParticleMCIndex', 'fj_particles_truth', 'ParticlePID']
      print('Merged output',self.df_fjparticles.columns)
      print(self.df_fjparticles)
    elif self.mcprod:
      self.df_fjparticles = pandas.concat([df_fjparticles_det, df_fjparticles_truth], axis=1)
      self.df_fjparticles.columns = ['fj_particles_det', 'ParticleCharge_det', 'ParticleMCid_det', 'fj_particles_truth', 'ParticleCharge_truth', 'ParticleMCid_truth']
      print('Merged output',self.df_fjparticles.columns)
      print(self.df_fjparticles)
      # print("last 20 rows")
      # print(self.df_fjparticles[-20:])
    else:
      self.df_fjparticles = pandas.concat([df_fjparticles_det, df_fjparticles_truth], axis=1)
      self.df_fjparticles.columns = ['fj_particles_det', 'fj_particles_truth']
    print('--- {} seconds ---'.format(time.time() - self.start_time))

    # ------------------------------------------------------------------------
    
    # Set up the Pb-Pb embedding object
    if not self.is_pp and not self.thermal_model:
      self.process_io_emb = process_io_emb.ProcessIO_Emb(self.emb_file_list, track_tree_name='tree_Particle', m=self.m)
    
    # ------------------------------------------------------------------------

    # Initialize histograms
    if not self.dry_run:
      self.initialize_output_objects()
    
    # Create constituent subtractor, if configured
    if self.do_constituent_subtraction:
      self.constituent_subtractor = [CEventSubtractor(max_distance=R_max, alpha=self.alpha, max_eta=self.max_eta, bge_rho_grid_size=self.bge_rho_grid_size, max_pt_correct=self.max_pt_correct, ghost_area=self.ghost_area, distance_type=fjcontrib.ConstituentSubtractor.deltaR) for R_max in self.max_distance]
    
    print(self)
    
    # Find jets and fill histograms
    print('Find jets...')
    self.analyze_events()

    if self.use_D0_info:
      print("There were", self.alld0counter_truth, " truth D0's")
      print("There were", self.d0nodstar_counter_truth, "truth D0's that did not come from charged D*")
      print("There were", self.alld0counter_det, "det D0's")
      print("There were", self.d0nodstar_counter_det, "det D0's that did not come from charged D*")
      print("On the truth matched side, there were", self.alld0counter_truthmatched, "D0's")
      print("On the truth matched side, there were", self.d0nodstar_counter_truthmatched, "D0's that did not come from charged D*")
      print("On the det matched side, there were", self.alld0counter_detmatched, "D0's")
      print("On the det matched side, there were", self.d0nodstar_counter_detmatched, "D0's that did not come from charged D*")
    
    # Plot histograms
    print('Save histograms...')
    process_base.ProcessBase.save_output_objects(self)
    
    print('--- {} seconds ---'.format(time.time() - self.start_time))

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_output_objects(self):
    
    self.hNevents = ROOT.TH1F('hNevents', 'hNevents', 2, -0.5, 1.5)
    # self.hNevents.Fill(1, self.nEvents_det)
    self.hNevents.Fill(1, self.nEvents_truth)
    
    self.hTrackEtaPhi = ROOT.TH2F('hTrackEtaPhi', 'hTrackEtaPhi', 200, -1., 1., 628, 0., 6.28)
    self.hTrackPt = ROOT.TH1F('hTrackPt', 'hTrackPt', 300, 0., 300.)
    
    if not self.is_pp:
      self.hRho =  ROOT.TH1F('hRho', 'hRho', 1000, 0., 1000.)
      
    if not self.skip_deltapt_RC_histograms:
      name = 'hN_MeanPt'
      h = ROOT.TH2F(name, name, 200, 0, 5000, 200, 0., 2.)
      setattr(self, name, h)

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_output_objects_R(self, jetR):
  
      # Call user-specific initialization
      self.initialize_user_output_objects_R(jetR)
      
      # Base histograms
      if self.is_pp:
      
          name = 'hJES_R{}'.format(jetR)
          h = ROOT.TH2F(name, name, 300, 0, 300, 200, -1., 1.)
          setattr(self, name, h)
      
          name = 'hDeltaR_All_R{}'.format(jetR)
          h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 2.)
          setattr(self, name, h)
          
      else:
      
          for R_max in self.max_distance:
          
            name = 'hJES_R{}_Rmax{}'.format(jetR, R_max)
            h = ROOT.TH2F(name, name, 300, 0, 300, 200, -1., 1.)
            setattr(self, name, h)
          
            name = 'hDeltaPt_emb_R{}_Rmax{}'.format(jetR, R_max)
            h = ROOT.TH2F(name, name, 300, 0, 300, 400, -200., 200.)
            setattr(self, name, h)
            
            if not self.skip_deltapt_RC_histograms:
              name = 'hDeltaPt_RC_beforeCS_R{}_Rmax{}'.format(jetR, R_max)
              h = ROOT.TH1F(name, name, 400, -200., 200.)
              setattr(self, name, h)
              
              name = 'hDeltaPt_RC_afterCS_R{}_Rmax{}'.format(jetR, R_max)
              h = ROOT.TH1F(name, name, 400, -200., 200.)
              setattr(self, name, h)
      
            name = 'hDeltaR_ppdet_pptrue_R{}_Rmax{}'.format(jetR, R_max)
            h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 2.)
            setattr(self, name, h)
            
            name = 'hDeltaR_combined_ppdet_R{}_Rmax{}'.format(jetR, R_max)
            h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 2.)
            setattr(self, name, h)
              
      name = 'hZ_Truth_R{}'.format(jetR)
      h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 1.)
      setattr(self, name, h)
      
      name = 'hZ_Det_R{}'.format(jetR)
      h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 1.)
      setattr(self, name, h)

      # new histograms for debugging
      name = 'hEtaRap_Truth_R{}'.format(jetR)
      h = ROOT.TH2F(name, name, 100, -1.5, 1.5, 100, -1.5, 1.5)
      setattr(self, name, h)

      name = 'hEtaRap_Det_R{}'.format(jetR)
      h = ROOT.TH2F(name, name, 100, -1.5, 1.5, 100, -1.5, 1.5)
      setattr(self, name, h)

  #---------------------------------------------------------------
  # Main function to loop through and analyze events
  #---------------------------------------------------------------
  def analyze_events(self):
    
    # # Fill track histograms
    # if not self.dry_run:
    #   [self.fill_track_histograms(fj_particles_det) for fj_particles_det in self.df_fjparticles['fj_particles_det']]
    
    fj.ClusterSequence.print_banner()
    print()
        
    self.event_number = 0
    self.jet_number = -1 # so that jet counting starts at 0
    
    for jetR in self.jetR_list:
      if not self.dry_run:
        self.initialize_output_objects_R(jetR)
    
    # Then can use list comprehension to iterate over the groupby and do jet-finding
    # simultaneously for fj_1 and fj_2 per event, so that I can match jets -- and fill histograms
    if self.jetscape:
      result = [self.analyze_event(fj_particles_det, fj_particles_truth, fj_particles_det_holes, fj_particles_truth_holes) for fj_particles_det, fj_particles_truth, fj_particles_det_holes, fj_particles_truth_holes in zip(self.df_fjparticles['fj_particles_det'], self.df_fjparticles['fj_particles_truth'], self.df_fjparticles['fj_particles_det_holes'], self.df_fjparticles['fj_particles_truth_holes'])]
    elif self.ENC_fastsim:
      if self.use_D0_info:
        # result = [self.analyze_event_nodet(fj_particles_truth=fj_particles_truth, particles_pid_truth=particles_pid_truth, particles_rap_truth=particles_rap_truth, particles_mid_truth=particles_mid_truth) for fj_particles_truth, particles_pid_truth, particles_mid_truth, particles_rap_truth in zip(self.df_fjparticles['fj_particles_truth'], self.df_fjparticles['ParticlePID'], self.df_fjparticles['MotherPID'],self.df_fjparticles['ParticleRapidity'])]
        result = [self.analyze_event(fj_particles_det=fj_particles_det, fj_particles_truth=fj_particles_truth, 
                                   particles_mcid_det=particles_mcid_det, particles_pid_truth=particles_pid_truth, 
                                   particles_pid_det=particles_pid_det, particles_mid_truth=particles_mid_truth, 
                                   particles_rap_truth=particles_rap_truth, D0s_truth=D0s_truth) 
                                   for fj_particles_det, fj_particles_truth, particles_mcid_det, 
                                   particles_pid_truth, particles_pid_det, particles_mid_truth, 
                                   particles_rap_truth, D0s_truth 
                                   in zip(self.df_fjparticles['fj_particles_det'], self.df_fjparticles['fj_particles_truth'], 
                                          self.df_fjparticles['ParticleMCIndex_det'], self.df_fjparticles['ParticlePID_truth'], 
                                          self.df_fjparticles['ParticlePID_det'], self.df_fjparticles['MotherPID_truth'], 
                                          self.df_fjparticles['ParticleRapidity_truth'], self.df_fjparticles['fj_D0_truth'])]
      else:
        #don't use no det for now...
        result = [self.analyze_event(fj_particles_det=fj_particles_det, fj_particles_truth=fj_particles_truth, particles_mcid_det=particles_mcid_det, particles_pid_truth=particles_pid_truth) for fj_particles_det, fj_particles_truth, particles_mcid_det, particles_pid_truth in zip(self.df_fjparticles['fj_particles_det'], self.df_fjparticles['fj_particles_truth'], self.df_fjparticles['ParticleMCIndex'], self.df_fjparticles['ParticlePID'])]
        # result = [self.analyze_event_nodet(fj_particles_truth=fj_particles_truth, particles_pid_truth=particles_pid_truth) for fj_particles_truth, particles_pid_truth in zip(self.df_fjparticles['fj_particles_truth'], self.df_fjparticles['ParticlePID'])]

      
    elif self.mcprod:
      self.crazycounter = 0
      result = [self.analyze_event(fj_particles_det=fj_particles_det, fj_particles_truth=fj_particles_truth, particles_mcid_det=particles_mcid_det, particles_mcid_truth=particles_mcid_truth, particles_charge_det=particles_charge_det, particles_charge_truth=particles_charge_truth) for fj_particles_det, fj_particles_truth, particles_mcid_det, particles_mcid_truth, particles_charge_det, particles_charge_truth in zip(self.df_fjparticles['fj_particles_det'], self.df_fjparticles['fj_particles_truth'], self.df_fjparticles['ParticleMCid_det'], self.df_fjparticles['ParticleMCid_truth'], self.df_fjparticles['ParticleCharge_det'], self.df_fjparticles['ParticleCharge_truth'])]
    else:
      result = [self.analyze_event(fj_particles_det, fj_particles_truth) for fj_particles_det, fj_particles_truth in zip(self.df_fjparticles['fj_particles_det'], self.df_fjparticles['fj_particles_truth'])]
    
    if self.debug_level > 0:
      for attr in dir(self):
        obj = getattr(self, attr)
        print('size of {}: {}'.format(attr, sys.getsizeof(obj)))
        
    print('Save thn...')
    process_base.ProcessBase.save_thn_th3_objects(self)
    
  #---------------------------------------------------------------
  # Fill track histograms.
  #---------------------------------------------------------------
  def fill_track_histograms(self, fj_particles_det):

    # Check that the entries exist appropriately
    # (need to check how this can happen -- but it is only a tiny fraction of events)
    if type(fj_particles_det) != fj.vectorPJ:
      return
    
    for track in fj_particles_det:
      self.hTrackEtaPhi.Fill(track.eta(), track.phi())
      self.hTrackPt.Fill(track.pt())
      
  #---------------------------------------------------------------
  # Analyze jets of a given event.
  # fj_particles is the list of fastjet pseudojets for a single fixed event.
  #---------------------------------------------------------------
  def analyze_event_nodet(self, fj_particles_truth, particles_pid_truth=None, particles_rap_truth=None, particles_mid_truth=None):
  
    self.event_number += 1
    if self.event_number > self.event_number_max:
      return
    if self.debug_level > 1:
      print('-------------------------------------------------')
      print('event {}'.format(self.event_number))

    # print('debug5 det parts',fj_particles_det)
    # print('debug5 mcid',particles_mcid_det)
    # print('debug5 truth parts',fj_particles_truth)
    # print('debug5 pid',particles_pid_truth)

    count_charged = 0
    count_neutral = 0
    count_d0 = 0

    if self.ENC_fastsim:
      # make charge array from pid info, needed for pair efficiency determination
        particles_charge_truth = np.array([])
        for pid in particles_pid_truth:
            # charged hadrons
            if abs(pid)==211 or abs(pid)==321 or abs(pid)==2212 or abs(pid)==3222:
                if pid>0:
                    particles_charge_truth = np.append(particles_charge_truth, 1)
                else:
                    particles_charge_truth = np.append(particles_charge_truth, -1)
                count_charged+=1
            # electrons and muons
            elif abs(pid)==11 or abs(pid)==13 or abs(pid)==3112 or abs(pid)==3312 or abs(pid)==3334:
                if pid>0:
                    particles_charge_truth = np.append(particles_charge_truth, -1)
                else:
                    particles_charge_truth = np.append(particles_charge_truth, 1)
                count_charged+=1
            # long lived weak decay particles (<2% of the total number of charged particles)
            # for now mark as charge 0 and later NOT applying pair efficiency for 0-charged or 0-0 pairs
            # NB: this can be avoided by decaying these paritcles within the generation step
            else:
                particles_charge_truth = np.append(particles_charge_truth, 0)
                count_neutral+=1
                if abs(pid)==421:
                  count_d0+=1

    if self.debug_level > 2:
      print("there are ", count_charged, "charged particles")
      print("there are ", count_neutral, "neutrl particles")
      print("there are ", count_d0, "d0 particles")
      if count_neutral != count_d0:
        print("There are neutrals that are not D0's!") # checked this - there aren't (for my herwig generation)


    # Check that the entries exist appropriately
    # (need to check how this can happen -- but it is only a tiny fraction of events)
    if type(fj_particles_truth) != fj.vectorPJ:
      print('fj_particles type mismatch -- skipping event')
      return
    else:
      # Todo
      ## for full simulation, match det-level and truth level particles
      ## sort both list by pT, phi and eta first before matching

      # add associated truth info and charge info in fj_particles_det using the JetInfo object
      if self.ENC_fastsim:

        for index in range( len(fj_particles_truth) ):
          if fj_particles_truth[index].has_user_info():
            ecorr_user_info = fj_particles_truth[index].python_info()
          else:
            # note: goes into here!!
            ecorr_user_info = jet_info.JetInfo()
          ecorr_user_info.particle_truth = fj_particles_truth[index]
          ecorr_user_info.charge = particles_charge_truth[index]
          if (self.use_D0_info):
            ecorr_user_info.particle_pid = particles_pid_truth[index]
            ecorr_user_info.particle_rap = particles_rap_truth[index]
            ecorr_user_info.particle_mid = particles_mid_truth[index]
            # if abs(particles_pid_truth[index]) == 421:
              # print("CROSS CHECKING", particles_rap_truth)
              # print("Analyzing event and found a D0 in event", self.event_number, "with pid", particles_pid_truth[index], "and rapidity", particles_rap_truth[index])
            # print(ecorr_user_info)
            
          else:
            ecorr_user_info.particle_pid = particles_pid_truth[index]
            ecorr_user_info.particle_rap = -99
            ecorr_user_info.particle_mid = -99

          fj_particles_truth[index].set_python_info(ecorr_user_info)
          # if rapidity needs to be saved, maybe here??
          # fj_particles_truth[index].set_user_index(int(index))

      
    if len(fj_particles_truth) > 1:
      if np.abs(fj_particles_truth[0].pt() - fj_particles_truth[1].pt()) <  1e-10:
        print('WARNING: Duplicate particles may be present')
        print([p.user_index() for p in fj_particles_truth])
        print([p.pt() for p in fj_particles_truth])


    if self.dry_run:
      return

    # Loop through jetR, and process event for each R
    for jetR in self.jetR_list:        
    
      # Keep track of whether to fill R-independent histograms
      self.fill_R_indep_hists = (jetR == self.jetR_list[0])

      # Set jet definition and a jet selector
      jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
      jet_selector_det = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9 - jetR)
      jet_selector_truth_matched = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9)
      if self.debug_level > 2:
        print('')
        print('jet definition is:', jet_def)
        print('jet selector for det-level is:', jet_selector_det)
        print('jet selector for truth-level matches is:', jet_selector_truth_matched)
      
      # Analyze

      # Find pp det jets
      # if self.ENC_fastsim:
      #   # FIX ME: should treat long lived charged particle differently (check how the existing fast herwig and pythia handles it)
      #   fj_particles_det_ch = fj.vectorPJ()
      #   if not (isinstance(fj_particles_det, float) and np.isnan(fj_particles_det)) and not self.use_D0_info: # check that det-level particles exist, and that we are not using D0 case! (bc D0 needs to be added to jet)
      #     for part in fj_particles_det:
      #       if part.python_info().charge!=0: # only use charged particles
      #         fj_particles_det_ch.append(part)
      #     cs_det = fj.ClusterSequence(fj_particles_det_ch, jet_def)
      #   else: # D0 case!
      #     if isinstance(fj_particles_det, float) and np.isnan(fj_particles_det): #make no detector level particles event into an empty PJ array
      #       fj_particles_det = fj.vectorPJ() 
      #     cs_det = fj.ClusterSequence(fj_particles_det, jet_def)
      # else:
      #   if isinstance(fj_particles_det, float) and np.isnan(fj_particles_det): #make no detector level particles event into an empty PJ array
      #     fj_particles_det = fj.vectorPJ() 
      #   cs_det = fj.ClusterSequence(fj_particles_det, jet_def)

      # jets_det_pp = fj.sorted_by_pt(cs_det.inclusive_jets())
      # # make sure the user info (on the jet side) for jets are all empty right after the jet-clustering 
      # for jet in jets_det_pp:
      #   if jet.has_user_info():
      #     jet.python_info().clear_jet_info()
      #   jets_det_pp_selected = jet_selector_det(jets_det_pp)

      # Find pp truth jets
      if self.ENC_fastsim:
        # FIX ME: should treat long lived charged particle differently (check how the existing fast herwig and pythia handles it)
        fj_particles_truth_ch = fj.vectorPJ()
        for part in fj_particles_truth:
          if part.python_info().charge!=0: # only use charged particles
            fj_particles_truth_ch.append(part)
        cs_truth = fj.ClusterSequence(fj_particles_truth_ch, jet_def)
      else:
        cs_truth = fj.ClusterSequence(fj_particles_truth, jet_def)

      jets_truth = fj.sorted_by_pt(cs_truth.inclusive_jets())
      # make sure the user info (on the jet side) for jets are all empty right after the jet-clustering  
      for jet in jets_truth:
        if jet.has_user_info():
          jet.python_info().clear_jet_info()
      jets_truth_selected = jet_selector_det(jets_truth) # has jet cuts of min pt = 5 gev, max eta = 0.9  - jetR --> USING THIS ONE
      jets_truth_selected_matched = jet_selector_truth_matched(jets_truth) # has jet cuts of min pt = 5 gev, max eta = 0.9
    
      # self.analyze_jets(jets_det_pp_selected, jets_truth_selected, jets_truth_selected_matched, jetR)
      self.analyze_jets(jets_truth_selected, jetR)

        
      
  def analyze_event(self, fj_particles_det, fj_particles_truth, fj_particles_det_holes=None, fj_particles_truth_holes=None, particles_mcid_det=None, particles_mcid_truth=None,
                    particles_charge_det=None, particles_charge_truth=None, particles_pid_truth=None, particles_pid_det=None, particles_mid_truth=None, particles_rap_truth=None, D0s_truth=None):
  
    # add condition to skip events that have no detector level particles in them?? -- maybe get rid of
    if self.event_number % 1000 == 0:
      print("EVENT", self.event_number)
    
    self.event_number += 1
    if self.event_number > self.event_number_max:
      return
    if self.debug_level > 1:
      print('-------------------------------------------------')
      print('event {}'.format(self.event_number))


    if self.ENC_fastsim:
      # make charge array from pid info, needed for pair efficiency determination
        particles_charge_truth = np.array([])
        for pid in particles_pid_truth:
            # charged hadrons
            if abs(pid)==211 or abs(pid)==321 or abs(pid)==2212 or abs(pid)==3222:
                if pid>0:
                    particles_charge_truth = np.append(particles_charge_truth, 1)
                else:
                    particles_charge_truth = np.append(particles_charge_truth, -1)
            # electrons and muons
            elif abs(pid)==11 or abs(pid)==13 or abs(pid)==3112 or abs(pid)==3312 or abs(pid)==3334:
                if pid>0:
                    particles_charge_truth = np.append(particles_charge_truth, -1)
                else:
                    particles_charge_truth = np.append(particles_charge_truth, 1)
            # long lived weak decay particles (<2% of the total number of charged particles)
            # for now mark as charge 0 and later NOT applying pair efficiency for 0-charged or 0-0 pairs
            # NB: this can be avoided by decaying these paritcles within the generation step
            else:
                particles_charge_truth = np.append(particles_charge_truth, 0)
      
    # Check that the entries exist appropriately
    # (need to check how this can happen -- but it is only a tiny fraction of events)
    if type(fj_particles_truth) != fj.vectorPJ:
      print('fj_particles type mismatch -- skipping event')
      return
    else:
      # Todo
      ## for full simulation, match det-level and truth level particles
      ## sort both list by pT, phi and eta first before matching

      # add associated truth info and charge info in fj_particles_det using the JetInfo object
      if self.ENC_fastsim:

        if isinstance(particles_mcid_det, float) and np.isnan(particles_mcid_det):
          print("Nan value - no detector level particles in this event.")
        else:
          if particles_mcid_det.ndim == 0: # this is a 0-d array but it looks like its saved as a "scalar" and it's the D0!
            print("particles_mcid_det", particles_mcid_det)
            print("this is a scalar -- only one detector-level particle in event!!!!")
            print("fj_particles_det[0]", fj_particles_det[0])
            if fj_particles_det[0].has_user_info():
              ecorr_user_info = fj_particles_det[0].python_info()
            else:
              ecorr_user_info = jet_info.JetInfo()
            
            # Get the PID of the corresponding truth particle
            if self.use_D0_info: # we have the particle pid saved here!!
              corresponding_truth_pid = particles_pid_det
            else:
              corresponding_truth_pid = particles_pid_truth[int(particles_mcid_det)]
            print("corresponding truth pid", corresponding_truth_pid)

            # Get the charge of the particle - there shouldn't be any neutrals (except D0 case)
            det_charge = int(corresponding_truth_pid / abs(corresponding_truth_pid))
            if abs(corresponding_truth_pid) == 421:
              det_charge = 0
            
            # Get the fj particle of the corresponding truth particle
            corresponding_truth_fj_particle = fj_particles_truth[int(abs(particles_mcid_det))]
            if abs(corresponding_truth_pid) == 421: # it's a D0 so the mcid will take it back to the D0 tree!
              print("D0 fj particle here!")

            ecorr_user_info.particle_mcid = int(particles_mcid_det)
            ecorr_user_info.particle_truth = corresponding_truth_fj_particle
            ecorr_user_info.particle_pid = corresponding_truth_pid
            ecorr_user_info.charge = det_charge 

            if (self.use_D0_info):
              ecorr_user_info.particle_rap = particles_rap_truth[int(abs(particles_mcid_det))]
              ecorr_user_info.particle_mid = particles_mid_truth[int(abs(particles_mcid_det))]
            
            fj_particles_det[0].set_python_info(ecorr_user_info)

          else:
            for index, mcid in enumerate(particles_mcid_det):
              if fj_particles_det[index].has_user_info():
                ecorr_user_info = fj_particles_det[index].python_info()
              else:
                ecorr_user_info = jet_info.JetInfo()

              # Get the PID of the corresponding truth particle
              if self.use_D0_info: # we have the particle pid saved here!!
                corresponding_truth_pid = particles_pid_det[int(index)]
              else:
                corresponding_truth_pid = particles_pid_truth[int(mcid)]

              # Get the charge of the particle - there shouldn't be any neutrals (except D0 case)
              det_charge = int(corresponding_truth_pid / abs(corresponding_truth_pid))
              if abs(corresponding_truth_pid) == 421:
                det_charge = 0
              
              # Get the fj particle of the corresponding truth particle
              corresponding_truth_fj_particle = fj_particles_truth[int(abs(mcid))]
              if abs(corresponding_truth_pid) == 421: # it's a D0 so the mcid will take it back to the D0 tree!
                print("D0 fj particle here!")

              ecorr_user_info.particle_mcid = int(mcid)
              ecorr_user_info.particle_truth = corresponding_truth_fj_particle
              ecorr_user_info.particle_pid = corresponding_truth_pid
              ecorr_user_info.charge = det_charge #int(particles_charge_det[index])

              if (self.use_D0_info):
                ecorr_user_info.particle_rap = particles_rap_truth[int(abs(mcid))]
                ecorr_user_info.particle_mid = particles_mid_truth[int(abs(mcid))]
              
              fj_particles_det[index].set_python_info(ecorr_user_info)
              # fj_particles_det[index].set_user_index(int(mcid))


        # now add appropriate information for truth particles
        for index, pid in enumerate(particles_pid_truth): #for index in range( len(fj_particles_truth) ):
          if fj_particles_truth[index].has_user_info():
            ecorr_user_info = fj_particles_truth[index].python_info()
          else:
            # note: goes into here!!
            ecorr_user_info = jet_info.JetInfo()

          # Get the charge of the particle - there shouldn't be any neutrals (except D0 case)
          truth_charge = int(pid / abs(pid))
          if abs(pid) == 421:
            truth_charge = 0

          ecorr_user_info.particle_mcid = index  #int(mcid)
          ecorr_user_info.particle_truth = fj_particles_truth[index]
          ecorr_user_info.particle_pid = pid
          ecorr_user_info.charge = truth_charge #particles_charge_truth[index] #int(particles_charge_truth[index])

          if (self.use_D0_info):
            ecorr_user_info.particle_rap = particles_rap_truth[index]
            ecorr_user_info.particle_mid = particles_mid_truth[index]
          # else:
          #   ecorr_user_info.particle_rap = -99
          #   ecorr_user_info.particle_mid = -99

          fj_particles_truth[index].set_python_info(ecorr_user_info)
          # if rapidity needs to be saved, maybe here??

          # fj_particles_truth[index].set_user_index(int(index))


      if self.mcprod:
        # save det level information first
        # print("particles mc id det", self.crazycounter, "and", particles_mcid_det)
        if isinstance(particles_mcid_det, float) and np.isnan(particles_mcid_det):
          print("Nan value - no detector level particles in this event.")
        else:
          # len(particles_mcid_det) >= 1 #there might be some events with no detector level particles, so this is to account for those
          for index, mcid in enumerate(particles_mcid_det):
            if fj_particles_det[index].has_user_info():
              ecorr_user_info = fj_particles_det[index].python_info()
            else:
              ecorr_user_info = jet_info.JetInfo()
            ecorr_user_info.particle_mcid = int(mcid)
            ecorr_user_info.charge = int(particles_charge_det[index])
            fj_particles_det[index].set_python_info(ecorr_user_info)
            # fj_particles_det[index].set_user_index(int(mcid))
            # self.crazycounter += 1

        # now save truth level information
        for index, mcid in enumerate(particles_mcid_truth):
          if fj_particles_truth[index].has_user_info():
            ecorr_user_info = fj_particles_truth[index].python_info()
          else:
            ecorr_user_info = jet_info.JetInfo()
          ecorr_user_info.particle_mcid = int(mcid)
          ecorr_user_info.charge = int(particles_charge_truth[index] / 3) #for some reason charge is saved as a multiple of 3?
          fj_particles_truth[index].set_python_info(ecorr_user_info)
          # fj_particles_truth[index].set_user_index(int(index))

    if self.jetscape:
      if type(fj_particles_det_holes) != fj.vectorPJ or type(fj_particles_truth_holes) != fj.vectorPJ:
        print('fj_particles_holes type mismatch -- skipping event')
        return
    
    if len(fj_particles_truth) > 1:
      if np.abs(fj_particles_truth[0].pt() - fj_particles_truth[1].pt()) <  1e-10:
        print('WARNING: Duplicate particles may be present')
        print([p.user_index() for p in fj_particles_truth])
        print([p.pt() for p in fj_particles_truth])

    # If Pb-Pb, construct embedded event (do this once, for all jetR)
    if not self.is_pp:
        
        # If thermal model, generate a thermal event and add it to the det-level particle list
        if self.thermal_model:
          fj_particles_combined_beforeCS = self.thermal_generator.load_event()
          
          # Form the combined det-level event
          # The pp-det tracks are each stored with a unique user_index >= 0
          #   (same index in fj_particles_combined and fj_particles_det -- which will be used in prong-matching)
          # The thermal tracks are each stored with a unique user_index < 0
          [fj_particles_combined_beforeCS.push_back(p) for p in fj_particles_det]

        # Main case: Get Pb-Pb event and embed it into the det-level particle list
        else:
          fj_particles_combined_beforeCS = self.process_io_emb.load_event()
              
          # Form the combined det-level event
          # The pp-det tracks are each stored with a unique user_index >= 0
          #   (same index in fj_particles_combined and fj_particles_det -- which will be used in prong-matching)
          # The Pb-Pb tracks are each stored with a unique user_index < 0
          [fj_particles_combined_beforeCS.push_back(p) for p in fj_particles_det]
         
        # Perform constituent subtraction for each R_max
        fj_particles_combined = [self.constituent_subtractor[i].process_event(fj_particles_combined_beforeCS) for i, R_max in enumerate(self.max_distance)]
        # for i, R_max in enumerate(self.max_distance):
        #   rho = self.constituent_subtractor[i].bge_rho.rho()
        #   print('rho is ',rho)
        # print('**************Before CS subtraction*****************')
        # n_sig_before = 0
        # n_bkg_before = 0
        # for part in fj_particles_combined_beforeCS:
        #   if part.user_index() < 0:
        #     n_bkg_before += 1
        #   else:
        #     n_sig_before += 1
        #     print('index checking:',part.user_index(),'pt',part.perp(),'phi',part.phi(),'eta',part.eta())
        # print('n_sig',n_sig_before,'n_bkg',n_bkg_before)
        # print('**************After CS subtraction*****************')
        # n_sig_after = 0
        # n_bkg_after = 0
        # for part in fj_particles_combined[0]:
        #   if part.user_index() < 0:
        #     n_bkg_after += 1
        #   else:
        #     n_sig_after += 1
        #     print('index checking:',part.user_index(),'pt',part.perp(),'phi',part.phi(),'eta',part.eta())
        # print('n_sig',n_sig_after,'n_bkg',n_bkg_after)
        # print('**************After CS subtraction*****************')
        
        if self.debug_level > 3:
          print([p.user_index() for p in fj_particles_truth])
          print([p.pt() for p in fj_particles_truth])
          print([p.user_index() for p in fj_particles_det])
          print([p.pt() for p in fj_particles_det])
          print([p.user_index() for p in fj_particles_combined_beforeCS])
          print([p.pt() for p in fj_particles_combined_beforeCS])
          
    if self.dry_run:
      return

    # Loop through jetR, and process event for each R
    for jetR in self.jetR_list:  

      # Keep track of whether to fill R-independent histograms
      self.fill_R_indep_hists = (jetR == self.jetR_list[0])

      # Set jet definition and a jet selector
      jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
      jet_selector_det = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9 - jetR)
      jet_selector_truth_matched = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9)
      if self.debug_level > 2:
        print('')
        print('jet definition is:', jet_def)
        print('jet selector for det-level is:', jet_selector_det)
        print('jet selector for truth-level matches is:', jet_selector_truth_matched)
      
      # Analyze
      if self.is_pp:

        # Find pp det and truth jets
        if self.ENC_fastsim:
          # FIX ME: should treat long lived charged particle differently (check how the existing fast herwig and pythia handles it)
          fj_particles_det_ch = fj.vectorPJ()
          if not (isinstance(fj_particles_det, float) and np.isnan(fj_particles_det)) and not self.use_D0_info: # check that det-level particles exist, and that we are not using D0 case! (bc D0 needs to be added to jet)
            for part in fj_particles_det:
              if part.python_info().charge!=0: # only use charged particles
                fj_particles_det_ch.append(part)
            cs_det = fj.ClusterSequence(fj_particles_det_ch, jet_def)
          else: # D0 case!
            if isinstance(fj_particles_det, float) and np.isnan(fj_particles_det): #make no detector level particles event into an empty PJ array
              fj_particles_det = fj.vectorPJ() 
            cs_det = fj.ClusterSequence(fj_particles_det, jet_def)
        else:
          if isinstance(fj_particles_det, float) and np.isnan(fj_particles_det): #make no detector level particles event into an empty PJ array
            fj_particles_det = fj.vectorPJ() 
          cs_det = fj.ClusterSequence(fj_particles_det, jet_def)
        
        # print("here")
        # [print("here charge wrong", p.python_info().charge) for p in fj_particles_truth if np.abs(p.python_info().charge) != 1 ] #TODO: is this relevant??
        # print("here", [p.user_index() for p in fj_particles_truth])
        # print([p.user_index() for p in fj_particles_det])
        # print([p.pt() for p in fj_particles_truth])
        # print([p.pt() for p in fj_particles_det])
        # print([p.python_info().particle_mcid for p in fj_particles_truth])
        # print([p.python_info().particle_mcid for p in fj_particles_det])
        # print([p.python_info().charge for p in fj_particles_truth])
        # print([p.python_info().charge for p in fj_particles_det])
        
        jets_det_pp = fj.sorted_by_pt(cs_det.inclusive_jets())
        # make sure the user info (on the jet side) for jets are all empty right after the jet-clustering 
        for jet in jets_det_pp:
          if jet.has_user_info():
            jet.python_info().clear_jet_info()
        jets_det_pp_selected = jet_selector_det(jets_det_pp)
        
        if self.ENC_fastsim:
          # FIX ME: should treat long lived charged particle differently (check how the existing fast herwig and pythia handles it)
          fj_particles_truth_ch = fj.vectorPJ()
          if not self.use_D0_info:
            for part in fj_particles_truth:
              if part.python_info().charge!=0: # only use charged particles
                fj_particles_truth_ch.append(part)
            cs_truth = fj.ClusterSequence(fj_particles_truth_ch, jet_def)
          else: # D0 case!
            cs_truth = fj.ClusterSequence(fj_particles_truth, jet_def)
        else:
          cs_truth = fj.ClusterSequence(fj_particles_truth, jet_def)

        jets_truth = fj.sorted_by_pt(cs_truth.inclusive_jets())
        # make sure the user info (on the jet side) for jets are all empty right after the jet-clustering  
        for jet in jets_truth:
          if jet.has_user_info():
            jet.python_info().clear_jet_info()
        jets_truth_selected = jet_selector_det(jets_truth)
        jets_truth_selected_matched = jet_selector_truth_matched(jets_truth)
      
        self.analyze_jets(jets_det_pp_selected, jets_truth_selected, jets_truth_selected_matched, jetR)
        # self.analyze_jets(jets_truth_selected, jetR)
        
        
      else:
        print("haha skipping byeee")
        ''' don't want to do any background stuff
        for i, R_max in enumerate(self.max_distance):
            
          if self.debug_level > 1:
            print('')
            print('R_max: {}'.format(R_max))
            print('Total number of combined particles: {}'.format(len([p.pt() for p in fj_particles_combined_beforeCS])))
            print('After constituent subtraction {}: {}'.format(i, len([p.pt() for p in fj_particles_combined[i])))
            
          # Keep track of whether to fill R_max-independent histograms
          self.fill_Rmax_indep_hists = (i == 0)
          
          # Perform constituent subtraction on det-level, if applicable
          self.fill_background_histograms(fj_particles_combined_beforeCS, fj_particles_combined[i], jetR, i)
          rho = self.constituent_subtractor[i].bge_rho.rho() 
      
          # Do jet finding (re-do each time, to make sure matching info gets reset)
          cs_det = fj.ClusterSequence(fj_particles_det, jet_def)
          jets_det_pp = fj.sorted_by_pt(cs_det.inclusive_jets())
          jets_det_pp_selected = jet_selector_det(jets_det_pp)
          
          cs_truth = fj.ClusterSequence(fj_particles_truth, jet_def)
          jets_truth = fj.sorted_by_pt(cs_truth.inclusive_jets())
          jets_truth_selected = jet_selector_det(jets_truth)
          jets_truth_selected_matched = jet_selector_truth_matched(jets_truth)
          
          cs_combined = fj.ClusterSequence(fj_particles_combined[i], jet_def)
          jets_combined = fj.sorted_by_pt(cs_combined.inclusive_jets())
          jets_combined_selected = jet_selector_det(jets_combined)

          
          if self.do_rho_subtraction:
            cs_combined_beforeCS = fj.ClusterSequenceArea(fj_particles_combined_beforeCS, jet_def, fj.AreaDefinition(fj.active_area_explicit_ghosts))
            jets_combined_beforeCS = fj.sorted_by_pt(cs_combined_beforeCS.inclusive_jets())
            jets_combined_selected_beforeCS = jet_selector_det(jets_combined_beforeCS)

            jets_combined_reselected_beforeCS = self.reselect_jets(jets_combined_selected_beforeCS, jetR, rho_bge = rho)

            if self.do_jetcone:
              self.analyze_jets(jets_combined_reselected_beforeCS, jets_truth_selected, jets_truth_selected_matched, jetR,
                            jets_det_pp_selected = jets_det_pp_selected, R_max = R_max,
                            fj_particles_det_holes = fj_particles_det_holes,
                            fj_particles_truth_holes = fj_particles_truth_holes, rho_bge = rho, fj_particles_det_cones = fj_particles_combined_beforeCS, fj_particles_truth_cones = fj_particles_truth)
            else:
              self.analyze_jets(jets_combined_reselected_beforeCS, jets_truth_selected, jets_truth_selected_matched, jetR,
                            jets_det_pp_selected = jets_det_pp_selected, R_max = R_max,
                            fj_particles_det_holes = fj_particles_det_holes,
                            fj_particles_truth_holes = fj_particles_truth_holes, rho_bge = rho)
          else:
            if self.do_jetcone:
              self.analyze_jets(jets_combined_selected, jets_truth_selected, jets_truth_selected_matched, jetR,
                            jets_det_pp_selected = jets_det_pp_selected, R_max = R_max,
                            fj_particles_det_holes = fj_particles_det_holes,
                            fj_particles_truth_holes = fj_particles_truth_holes, rho_bge = 0, fj_particles_det_cones = fj_particles_combined_beforeCS, fj_particles_truth_cones = fj_particles_truth) # NB: feed all particles for cone around the CS subtracted jet. An alternate way is to use CS subtracted particles
            else:
              self.analyze_jets(jets_combined_selected, jets_truth_selected, jets_truth_selected_matched, jetR,
                            jets_det_pp_selected = jets_det_pp_selected, R_max = R_max,
                            fj_particles_det_holes = fj_particles_det_holes,
                            fj_particles_truth_holes = fj_particles_truth_holes, rho_bge = 0)
        '''

  

  #---------------------------------------------------------------
  # Jet selection cuts.
  #---------------------------------------------------------------
  def reselect_jets(self, jets_selected, jetR, rho_bge = 0):
    # re-apply jet pt > 5GeV cut after rho subtraction and leading track pt cut if there is any. NB: need to be applied inside apply_events to make sure matching work properly
    jets_reselected = []
    for jet in jets_selected:
      is_jet_selected = True
      
      # leading track selection
      if self.leading_pt > 0:
        constituents = fj.sorted_by_pt(jet.constituents())
        if constituents[0].perp() < self.leading_pt:
          is_jet_selected = False
      
      # if rho subtraction, require jet pt > 5 after subtration
      if self.do_rho_subtraction and rho_bge > 0:
        if jet.perp()-rho_bge*jet.area() < 5:
          # FIX ME: not sure whether to apply the area selection or not yet. jet.area() > 0.6*np.pi*jetR*jetR
          is_jet_selected = False

      if is_jet_selected:
        jets_reselected.append(jet)

    return jets_reselected

  #---------------------------------------------------------------
  # Analyze jets of a given event.
  #---------------------------------------------------------------
  def analyze_jets(self, jets_det_selected, jets_truth_selected, jets_truth_selected_matched, jetR,
                   jets_det_pp_selected = None, R_max = None,
                   fj_particles_det_holes = None, fj_particles_truth_holes = None, rho_bge = 0, fj_particles_det_cones = None, fj_particles_truth_cones = None):
  # def analyze_jets(self, jets_truth_selected, jetR,
  #                  fj_particles_truth_holes = None, rho_bge = 0, fj_particles_truth_cones = None):
  
    if self.debug_level > 1 and self.debug_level != 3:
      print('Number of det-level jets: {}'.format(len(jets_det_selected)))

    # print("-- Number of det jets in this event: ", len(jets_det_selected), "--")
    
    # i dont care about det-level
    # Fill det-level jet histograms (before matching)
    for ijet,jet_det in enumerate(jets_det_selected):

      self.ijet = ijet
      
      # Check additional acceptance criteria
      # skip event if not satisfied -- since first jet in event is highest pt
      if not self.utils.is_det_jet_accepted(jet_det):
        if self.fill_R_indep_hists:
          self.hNevents.Fill(0)
        if self.debug_level > 1:
          print('event rejected due to jet acceptance')
        return
      
      self.fill_det_before_matching(jet_det, jetR, R_max, rho_bge)
    
    # print("-- Number of truth jets in this event: ", len(jets_truth_selected), "--")
  
    # Fill truth-level jet histograms (before matching)
    for ijet,jet_truth in enumerate(jets_truth_selected):

      self.jet_number += 1 #starts counting jets at 0
      self.ijet = ijet

      leading_parton = fj.sorted_by_pt(jet_truth.constituents())[0]
      leading_parton_pt = leading_parton.pt()
      if (leading_parton_pt < self.leading_parton_pt_cut):
        # print("leading parton pt cut!!!, skipping jet", leading_parton_pt)
        continue
    
      if self.is_pp or self.fill_Rmax_indep_hists:
        self.fill_truth_before_matching(jet_truth, jetR)


    # NOW DO MATCHING!
    # Loop through jets and set jet matching candidates for each jet in user_info
    if self.is_pp:
        [[self.set_matching_candidates(jet_det, jet_truth, jetR, 'hDeltaR_All_R{}'.format(jetR)) for jet_truth in jets_truth_selected_matched] for jet_det in jets_det_selected]
    else:
        # First fill the combined-to-pp matches, then the pp-to-pp matches
        [[self.set_matching_candidates(jet_det_combined, jet_det_pp, jetR, 'hDeltaR_combined_ppdet_R{{}}_Rmax{}'.format(R_max), fill_jet1_matches_only=True) for jet_det_pp in jets_det_pp_selected] for jet_det_combined in jets_det_selected]
        [[self.set_matching_candidates(jet_det_pp, jet_truth, jetR, 'hDeltaR_ppdet_pptrue_R{{}}_Rmax{}'.format(R_max)) for jet_truth in jets_truth_selected_matched] for jet_det_pp in jets_det_pp_selected]

    # # debug
    # for jet_det_combined in jets_det_selected:
    #   print('debug7.1--jet_det',jet_det_combined.pt(),'user_info',jet_det_combined.has_user_info())
    #   if jet_det_combined.has_user_info() and jet_det_combined.python_info().closest_jet:
    #     print('matches to',jet_det_combined.python_info().closest_jet.pt())
    #     print('debug7.1--jet_det',len(jet_det_combined.constituents()))
    #     print('matches to',len(jet_det_combined.python_info().closest_jet.constituents()))
        
    # Loop through jets and set accepted matches
    if self.is_pp:
        hname = 'hJetMatchingQA_R{}'.format(jetR)
        [self.set_matches_pp(jet_det, hname, self.use_D0_info) for jet_det in jets_det_selected]
    else:
        hname = 'hJetMatchingQA_R{}_Rmax{}'.format(jetR, R_max)
        [self.set_matches_AA(jet_det_combined, jetR, hname) for jet_det_combined in jets_det_selected]
          
    # Loop through jets and fill response histograms if both det and truth jets are unique match
    result = [self.fill_jet_matches(jet_det, jetR, R_max, fj_particles_det_holes, fj_particles_truth_holes, rho_bge, fj_particles_det_cones, fj_particles_truth_cones) for jet_det in jets_det_selected]
    

  #---------------------------------------------------------------
  # Fill some background histograms
  #---------------------------------------------------------------
  def fill_background_histograms(self, fj_particles_combined_beforeCS, fj_particles_combined, jetR, i):

    # Fill rho
    rho = self.constituent_subtractor[i].bge_rho.rho()
    if self.fill_R_indep_hists and self.fill_Rmax_indep_hists:
      getattr(self, 'hRho').Fill(rho)
    
    # Fill random cone delta-pt before constituent subtraction
    if not self.skip_deltapt_RC_histograms:
      R_max = self.max_distance[i]
      self.fill_deltapt_RC_histogram(fj_particles_combined_beforeCS, rho, jetR, R_max, before_CS=True)
          
      # Fill random cone delta-pt after constituent subtraction
      self.fill_deltapt_RC_histogram(fj_particles_combined, rho, jetR, R_max, before_CS=False)
    
  #---------------------------------------------------------------
  # Fill delta-pt histogram
  #---------------------------------------------------------------
  def fill_deltapt_RC_histogram(self, fj_particles, rho, jetR, R_max, before_CS=False):
  
    # Choose a random eta-phi in the fiducial acceptance
    phi = random.uniform(0., 2*np.pi)
    eta = random.uniform(-0.9+jetR, 0.9-jetR)
    
    # Loop through tracks and sum pt inside the cone
    pt_sum = 0.
    pt_sum_global = 0.
    for track in fj_particles:
        if self.utils.delta_R(track, eta, phi) < jetR:
            pt_sum += track.pt()
        pt_sum_global += track.pt()
            
    if before_CS:
        delta_pt = pt_sum - rho * np.pi * jetR * jetR
        getattr(self, 'hDeltaPt_RC_beforeCS_R{}_Rmax{}'.format(jetR, R_max)).Fill(delta_pt)
    else:
        delta_pt = pt_sum
        getattr(self, 'hDeltaPt_RC_afterCS_R{}_Rmax{}'.format(jetR, R_max)).Fill(delta_pt)
        
    # Fill mean pt
    if before_CS and self.fill_R_indep_hists and self.fill_Rmax_indep_hists:
      N_tracks = len(fj_particles)
      mean_pt = pt_sum_global/N_tracks
      getattr(self, 'hN_MeanPt').Fill(N_tracks, mean_pt)

  #---------------------------------------------------------------
  # Fill truth jet histograms
  #---------------------------------------------------------------
  def fill_truth_before_matching(self, jet, jetR):
    
    jet_pt = jet.pt()
    for constituent in jet.constituents():
      z = constituent.pt() / jet.pt()
      getattr(self, 'hZ_Truth_R{}'.format(jetR)).Fill(jet.pt(), z)
    
    # Fill 2D histogram of truth (pt, obs)
    hname = 'h_{{}}_JetPt_Truth_R{}_{{}}'.format(jetR)
    self.fill_unmatched_jet_histograms(jet, jetR, hname)

    getattr(self, 'hEtaRap_Truth_R{}'.format(jetR)).Fill(jet.eta(), jet.rap())

  #---------------------------------------------------------------
  # Fill det jet histograms
  #---------------------------------------------------------------
  def fill_det_before_matching(self, jet, jetR, R_max, rho_bge = 0):
    
    if self.is_pp or self.fill_Rmax_indep_hists:
      jet_pt = jet.pt()
      if self.do_rho_subtraction:
        jet_pt = jet.pt()-rho_bge*jet.area()
      for constituent in jet.constituents():
        z = constituent.pt() / jet_pt
        getattr(self, 'hZ_Det_R{}'.format(jetR)).Fill(jet_pt, z)
      
      getattr(self, 'hEtaRap_Det_R{}'.format(jetR)).Fill(jet.eta(), jet.rap())
      
    # for const in jet.constituents():
    #   if const.perp()>0.15:
    #     print('index',const.user_index(),'pt',const.perp())
    # print('fill det hist')
    
    # Fill groomed histograms
    if self.thermal_model:
      # hname = 'h_{{}}_JetPt_R{}_{{}}_Rmax{}'.format(jetR, R_max)
      hname = 'h_{{}}_JetPt_Det_R{}_{{}}_Rmax{}'.format(jetR, R_max)
      self.fill_unmatched_jet_histograms(jet, jetR, hname, rho_bge)

    if self.is_pp:
      hname = 'h_{{}}_JetPt_Det_R{}_{{}}'.format(jetR)
      self.fill_unmatched_jet_histograms(jet, jetR, hname, rho_bge)

    if self.do_rho_subtraction:
      # hname = 'h_{{}}_JetPt_R{}_{{}}'.format(jetR)
      hname = 'h_{{}}_JetPt_Det_R{}_{{}}'.format(jetR)
      self.fill_unmatched_jet_histograms(jet, jetR, hname, rho_bge)
  
  #---------------------------------------------------------------
  # This function is called once for each jet
  #---------------------------------------------------------------
  def fill_unmatched_jet_histograms(self, jet, jetR, hname, rho_bge = 0):

    # Loop through each jet subconfiguration (i.e. subobservable / grooming setting)
    observable = self.observable_list[0]
    for i in range(len(self.obs_settings[observable])):
      if i==0:
        self.firsttimejet = True
      else:
        self.firsttimejet = False

      obs_setting = self.obs_settings[observable][i]
      grooming_setting = self.obs_grooming_settings[observable][i]
      obs_label = self.utils.obs_label(obs_setting, grooming_setting)

      # Groom jet, if applicable
      if grooming_setting:
        gshop = fjcontrib.GroomerShop(jet, jetR, self.reclustering_algorithm)
        jet_groomed_lund = self.utils.groom(gshop, grooming_setting, jetR)
        if not jet_groomed_lund:
          continue
      else:
        jet_groomed_lund = None
        
      if self.do_rho_subtraction and rho_bge > 0:
        jet_pt = jet.perp()-rho_bge*jet.area() # use subtracted jet pt for energy weight calculation and pt selection for there is a non-zero UE energy density
      else:
        jet_pt = jet.perp()

      # Call user function to fill histograms
      # print("filling here!")
      self.fill_observable_histograms(hname, jet, jet_groomed_lund, jetR, obs_setting,
                                      grooming_setting, obs_label, jet_pt)
  
  def find_parts_around_jet(self, parts, jet, cone_R):
    # select particles around jet axis
    cone_parts = fj.vectorPJ()
    for part in parts:
      if jet.delta_R(part) <= cone_R:
        cone_parts.push_back(part)
    
    return cone_parts

  #---------------------------------------------------------------
  # Loop through jets and call user function to fill matched
  # histos if both det and truth jets are unique match.
  #---------------------------------------------------------------
  def fill_jet_matches(self, jet_det, jetR, R_max, fj_particles_det_holes, fj_particles_truth_holes, rho_bge = 0, fj_particles_det_cones = None, fj_particles_truth_cones = None):
  
    # Set suffix for filling histograms
    if R_max:
      suffix = '_Rmax{}'.format(R_max)
    else:
      suffix = ''
    
    # Get matched truth jet
    if jet_det.has_user_info():
      jet_truth = jet_det.python_info().match
      if self.do_rho_subtraction and rho_bge > 0:
        jet_det_pt = jet_det.perp()-rho_bge*jet_det.area() # use subtracted jet pt for energy weight calculation and pt selection for there is a non-zero UE energy density
      else:
        jet_det_pt = jet_det.perp()

      if jet_truth:

        # # debug
        # print('debug8--jet det', jet_det_pt, 'size', len(jet_det.constituents()))
        # print('debug8--jet_truth', jet_truth.pt(), 'size', len(jet_truth.constituents()))
        
        jet_pt_det_ungroomed = jet_det_pt
        jet_pt_truth_ungroomed = jet_truth.pt()
        JES = (jet_pt_det_ungroomed - jet_pt_truth_ungroomed) / jet_pt_truth_ungroomed
        getattr(self, 'hJES_R{}{}'.format(jetR, suffix)).Fill(jet_pt_truth_ungroomed, JES)
        
        # If Pb-Pb case, we need to keep jet_det, jet_truth, jet_pp_det
        jet_pp_det = None
        if not self.is_pp:
        
          # Get pp-det jet
          jet_pp_det = jet_truth.python_info().match
            
          # Fill delta-pt histogram
          if jet_pp_det:
            jet_pp_det_pt = jet_pp_det.pt()
            delta_pt = (jet_pt_det_ungroomed - jet_pp_det_pt)
            getattr(self, 'hDeltaPt_emb_R{}_Rmax{}'.format(jetR, R_max)).Fill(jet_pt_truth_ungroomed, delta_pt)
            
        # Loop through each jet subconfiguration (i.e. subobservable / grooming setting)
        observable = self.observable_list[0]
        for i in range(len(self.obs_settings[observable])):
          
          if i==0:
            self.firsttimejet = True
          else:
            self.firsttimejet = False

          obs_setting = self.obs_settings[observable][i]
          grooming_setting = self.obs_grooming_settings[observable][i]
          obs_label = self.utils.obs_label(obs_setting, grooming_setting)
          
          if self.debug_level > 3:
            print('obs_label: {}'.format(obs_label))
          
          # Groom jets, if applicable
          if grooming_setting:
                    
            # Groom det jet
            gshop_det = fjcontrib.GroomerShop(jet_det, jetR, self.reclustering_algorithm)
            jet_det_groomed_lund = self.utils.groom(gshop_det, grooming_setting, jetR)
            if not jet_det_groomed_lund:
              continue

            # Groom truth jet
            gshop_truth = fjcontrib.GroomerShop(jet_truth, jetR, self.reclustering_algorithm)
            jet_truth_groomed_lund = self.utils.groom(gshop_truth, grooming_setting, jetR)
            if not jet_truth_groomed_lund:
              continue
              
          else:
          
            jet_det_groomed_lund = None
            jet_truth_groomed_lund = None
            
          # If jetscape, pass the list of holes within R of the jet to the user
          holes_in_det_jet = None
          holes_in_truth_jet = None
          if self.jetscape:
            holes_in_det_jet = [hadron for hadron in fj_particles_det_holes if jet_det.delta_R(hadron) < jetR]
            holes_in_truth_jet = [hadron for hadron in fj_particles_truth_holes if jet_truth.delta_R(hadron) < jetR]
            
            # Get the corrected jet pt by subtracting the negative recoils within R
            for hadron in holes_in_det_jet:
                jet_pt_det_ungroomed -= hadron.pt()
                
            for hadron in holes_in_truth_jet:
                jet_pt_truth_ungroomed -= hadron.pt()
          
          # Call user function to fill histos
          self.fill_matched_jet_histograms(jet_det, jet_det_groomed_lund, jet_truth,
                               jet_truth_groomed_lund, jet_pp_det, jetR,
                               obs_setting, grooming_setting, obs_label,
                               jet_pt_det_ungroomed, jet_pt_truth_ungroomed,
                               R_max, suffix, holes_in_det_jet=holes_in_det_jet,
                               holes_in_truth_jet=holes_in_truth_jet, cone_parts_in_det_jet=None, cone_parts_in_truth_jet=None, cone_R=0)
          
          # If check cone, pass the list of cone particles
          if self.do_jetcone:
            for jetcone_R in self.jetcone_R_list:
              
              cone_parts_in_det_jet = self.find_parts_around_jet(fj_particles_det_cones, jet_det, jetcone_R)
              cone_parts_in_truth_jet = self.find_parts_around_jet(fj_particles_truth_cones, jet_truth, jetcone_R)

              # Call user function to fill histos
              self.fill_matched_jet_histograms(jet_det, jet_det_groomed_lund, jet_truth,
                                 jet_truth_groomed_lund, jet_pp_det, jetR,
                                 obs_setting, grooming_setting, obs_label,
                                 jet_pt_det_ungroomed, jet_pt_truth_ungroomed,
                                 R_max, suffix, holes_in_det_jet=holes_in_det_jet,
                                 holes_in_truth_jet=holes_in_truth_jet, cone_parts_in_det_jet=cone_parts_in_det_jet, cone_parts_in_truth_jet=cone_parts_in_truth_jet, cone_R=jetcone_R)

  #---------------------------------------------------------------
  # Fill response histograms -- common utility function
  #---------------------------------------------------------------
  def fill_response(self, observable, jetR, jet_pt_det_ungroomed, jet_pt_truth_ungroomed,
                    obs_det, obs_truth, obs_label, R_max, prong_match = False):

    if self.fill_RM_histograms:
      x = ([jet_pt_det_ungroomed, jet_pt_truth_ungroomed, obs_det, obs_truth])
      x_array = array('d', x)
      name = 'hResponse_JetPt_{}_R{}_{}'.format(observable, jetR, obs_label)
      if not self.is_pp:
        name += '_Rmax{}'.format(R_max)
      getattr(self, name).Fill(x_array)
      
    if obs_truth > 1e-5:
      obs_resolution = (obs_det - obs_truth) / obs_truth
      name = 'hResidual_JetPt_{}_R{}_{}'.format(observable, jetR, obs_label)
      if not self.is_pp:
        name += '_Rmax{}'.format(R_max)
      getattr(self, name).Fill(jet_pt_truth_ungroomed, obs_truth, obs_resolution)
    
    # Fill prong-matched response
    if not self.is_pp and R_max == self.main_R_max:
      if prong_match:
      
        name = 'hResponse_JetPt_{}_R{}_{}_Rmax{}_matched'.format(observable, jetR, obs_label, R_max)
        getattr(self, name).Fill(x_array)
        
        if obs_truth > 1e-5:
          name = 'hResidual_JetPt_{}_R{}_{}_Rmax{}_matched'.format(observable, jetR, obs_label, R_max)
          getattr(self, name).Fill(jet_pt_truth_ungroomed, obs_truth, obs_resolution)

  #---------------------------------------------------------------
  # This function is called once for each jetR
  # You must implement this
  #---------------------------------------------------------------
  def initialize_user_output_objects_R(self, jetR):
      
    raise NotImplementedError('You must implement initialize_user_output_objects_R()!')

  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  # You must implement this
  #---------------------------------------------------------------
  def fill_observable_histograms(self, hname, jet, jet_groomed_lund, jetR, obs_setting,
                                 grooming_setting, obs_label, jet_pt_ungroomed):

    raise NotImplementedError('You must implement fill_observable_histograms()!')

  #---------------------------------------------------------------
  # This function is called once for each matched jet subconfiguration
  # You must implement this
  #---------------------------------------------------------------
  def fill_matched_jet_histograms(self, jet_det, jet_det_groomed_lund, jet_truth,
                                  jet_truth_groomed_lund, jet_pp_det, jetR,
                                  obs_setting, grooming_setting, obs_label,
                                  jet_pt_det_ungroomed, jet_pt_truth_ungroomed,
                                  R_max, suffix,
                                  **kwargs):

    raise NotImplementedError('You must implement fill_matched_jet_histograms()!')