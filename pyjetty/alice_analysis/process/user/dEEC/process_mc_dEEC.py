#!/usr/bin/env python3

"""
  Analysis class to read a ROOT TTree of MC track information
  and do jet-finding, and save response histograms.
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import argparse

# Data analysis and plotting
import numpy as np
import ROOT
import yaml
import array
import math
# from array import *

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjtools
import ecorrel
import othercorrel

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io
from pyjetty.alice_analysis.process.base import process_io_emb
from pyjetty.alice_analysis.process.base import jet_info
from pyjetty.alice_analysis.process.user.substructure import process_mc_base
from pyjetty.alice_analysis.process.base import thermal_generator
from pyjetty.mputils.csubtractor import CEventSubtractor

def linbins(xmin, xmax, nbins):
  lspace = np.linspace(xmin, xmax, nbins+1)
  arr = array.array('d', lspace) #'f', lspace) # needs to be d to be bins for histograms??
  return arr

def logbins(xmin, xmax, nbins):
  lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
  arr = array.array('d', lspace) #'f', lspace)
  return arr

################################################################
class ProcessMC_dEEC(process_mc_base.ProcessMCBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', save_tuples=1, debug_level=0, **kwargs):
  
    # Initialize base class
    super(ProcessMC_dEEC, self).__init__(input_file, config_file, output_dir, save_tuples, debug_level, **kwargs)
    
    self.observable = self.observable_list[0]

    if self.ENC_fastsim:
      self.pair_eff_file = ROOT.TFile.Open(self.pair_eff_file,"READ")
      # self.dpbin = 5
      # self.dp_lo = [0, 0.1, 0.2, 0.4, 1]
      # self.dp_hi = [0.1, 0.2, 0.4, 1, 2]
      self.dpbin = 10
      self.dp_lo = [0, 0.02, 0.06, 0.1, 0.14, 0.2, 0.3, 0.4, 0.6, 1]
      self.dp_hi = [0.02, 0.06, 0.1, 0.14, 0.2, 0.3, 0.4, 0.6, 1, 2]
      self.h1d_eff_vs_dR_in_dq_over_p = []
      for idp in range(self.dpbin):
          hname = 'h1d_eff_vs_dR_in_dq_over_p_{}'.format(idp)
          self.h1d_eff_vs_dR_in_dq_over_p.append( ROOT.TH1D(self.pair_eff_file.Get(hname)) )

  #---------------------------------------------------------------
  # Determine pair efficiency with the pair
  # property and input histograms
  #---------------------------------------------------------------
  def get_pair_eff(self, dist, dq_over_p):
    # return pair efficiency (from 0 to 1)
    idpbin = -9999
    for idp in range(self.dpbin):
        if math.fabs(dq_over_p)>=self.dp_lo[idp] and math.fabs(dq_over_p)<self.dp_hi[idp]:
            idpbin = idp
    
    pair_eff = 1 # set pair efficeincy to 1 if dq_over_p>=2
    if idpbin>=0:
      if math.log10(dist)<0 and math.log10(dist)>-3:
          ibin = self.h1d_eff_vs_dR_in_dq_over_p[idpbin].FindBin(math.log10(dist))
          pair_eff = self.h1d_eff_vs_dR_in_dq_over_p[idpbin].GetBinContent(ibin)
      elif math.log10(dist)>=0:
          pair_eff = 1 # overflow
      else:
          pair_eff = 0 # NB: underflow set to 0 efficiency. Maybe too aggressive but should be fine since we plan to measure down to dist ~1E-2
    
    return pair_eff

  #---------------------------------------------------------------
  # Calculate pair distance of two fastjet particles
  #---------------------------------------------------------------
  def calculate_distance(self, p0, p1):   
    dphiabs = math.fabs(p0.phi() - p1.phi())
    dphi = dphiabs

    if dphiabs > math.pi:
      dphi = 2*math.pi - dphiabs

    deta = p0.eta() - p1.eta()
    return math.sqrt(deta*deta + dphi*dphi)

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_user_output_objects_R(self, jetR):

    self.fout.cd()

    for observable in self.observable_list:
      print("HISTOGRAM OBSERVABLE", observable)

      for trk_thrd in self.obs_settings[observable]:
        # print("CP TRK THRD", trk_thrd)

        obs_label = self.utils.obs_label(trk_thrd, None) 
        # print("OBS_LABEL OG", obs_label)

        self.pair_type_labels = ['']
        if self.do_rho_subtraction or self.do_constituent_subtraction:
          self.pair_type_labels = ['_bb','_sb','_ss']

        # Init ENC histograms (both det and truth level)
        # print("pair type label", self.pair_type_labels)
        for pair_type_label in self.pair_type_labels:
            if 'ENC' in observable:
              for ipoint in range(2, 3):

                pt_bins = linbins(0,200,200)
                RL_bins = logbins(1E-4,1,50)
                ptRL_bins = logbins(1E-3,1E2,60)

                # Truth histograms
                name = 'h_{}{}{}_JetPt_Truth_R{}_{}'.format(observable, ipoint, pair_type_label, jetR, obs_label) #pair_type_label is blank for me
                h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
                h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
                h.GetYaxis().SetTitle('#it{R}_{L}')
                setattr(self, name, h)

                name = 'h_{}{}{}Pt_JetPt_Truth_R{}_{}'.format(observable, ipoint, pair_type_label, jetR, obs_label) # pt scaled histograms (currently only for unmatched jets)
                h = ROOT.TH2D(name, name, 200, pt_bins, 60, ptRL_bins)
                h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
                h.GetYaxis().SetTitle('#it{p}_{T,ch jet}#it{R}_{L}') # NB: y axis scaled by jet pt (applied jet by jet)
                setattr(self, name, h)

                # Det-level histograms
                name = 'h_{}{}{}_JetPt_Det_R{}_{}'.format(observable, ipoint, pair_type_label, jetR, obs_label) #pair_type_label is blank for me
                RL_bins = logbins(1E-4,1,50)
                h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
                h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
                h.GetYaxis().SetTitle('#it{R}_{L}')
                setattr(self, name, h)

                name = 'h_{}{}{}Pt_JetPt_Det_R{}_{}'.format(observable, ipoint, pair_type_label, jetR, obs_label) # pt scaled histograms (currently only for unmatched jets)
                h = ROOT.TH2D(name, name, 200, pt_bins, 60, ptRL_bins)
                h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
                h.GetYaxis().SetTitle('#it{p}_{T,ch jet}#it{R}_{L}') # NB: y axis scaled by jet pt (applied jet by jet)
                setattr(self, name, h)


                # Matched det histograms
                name = 'h_matched_{}{}{}_JetPt_R{}_{}'.format(observable, ipoint, pair_type_label, jetR, obs_label)
                h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
                h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
                h.GetYaxis().SetTitle('#it{R}_{L}')
                setattr(self, name, h)

                # Matched det histograms (with matched truth jet pT filled to the other axis)
                name = 'h_matched_extra_{}{}{}_JetPt_R{}_{}'.format(observable, ipoint, pair_type_label, jetR, obs_label)
                h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
                h.GetXaxis().SetTitle('#it{p}_{T,ch jet}^{truth}')
                h.GetYaxis().SetTitle('#it{R}_{L}')
                setattr(self, name, h)

                # Matched truth histograms
                name = 'h_matched_{}{}{}_JetPt_Truth_R{}_{}'.format(observable, ipoint, pair_type_label, jetR, obs_label)
                h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
                h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
                h.GetYaxis().SetTitle('#it{R}_{L}')
                setattr(self, name, h)

                if self.do_jetcone:
                  for jetcone_R in self.jetcone_R_list:
                    # Matched det histograms
                    name = 'h_jetcone{}_matched_{}{}{}_JetPt_R{}_{}'.format(jetcone_R, observable, ipoint, pair_type_label, jetR, obs_label)
                    h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
                    h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
                    h.GetYaxis().SetTitle('#it{R}_{L}')
                    setattr(self, name, h)

                    # Matched det histograms (with matched truth jet pT filled to the other axis)
                    name = 'h_jetcone{}_matched_extra_{}{}{}_JetPt_R{}_{}'.format(jetcone_R, observable, ipoint, pair_type_label, jetR, obs_label)
                    h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
                    h.GetXaxis().SetTitle('#it{p}_{T,ch jet}^{truth}')
                    h.GetYaxis().SetTitle('#it{R}_{L}')
                    setattr(self, name, h)

                    # Matched truth histograms
                    name = 'h_jetcone{}_matched_{}{}{}_JetPt_Truth_R{}_{}'.format(jetcone_R, observable, ipoint, pair_type_label, jetR, obs_label)
                    h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
                    h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
                    h.GetYaxis().SetTitle('#it{R}_{L}')
                    setattr(self, name, h)
                
                if self.thermal_model:
                  for R_max in self.max_distance:
                    name = 'h_{}{}{}_JetPt_R{}_{}_Rmax{}'.format(observable, ipoint, pair_type_label, jetR, obs_label, R_max)
                    h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
                    h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
                    h.GetYaxis().SetTitle('#it{R}_{L}')
                    setattr(self, name, h)

                
            if 'EEC_noweight' in observable:
              pt_bins = linbins(0,200,200)
              RL_bins = logbins(1E-4,1,50)
              
              # Truth histograms
              name = 'h_{}{}_JetPt_Truth_R{}_{}'.format(observable, pair_type_label, jetR, obs_label)  
              h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
              h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
              h.GetYaxis().SetTitle('#it{R}_{L}')
              setattr(self, name, h)

              # Det histogram
              name = 'h_{}{}_JetPt_Det_R{}_{}'.format(observable, pair_type_label, jetR, obs_label)
              h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
              h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
              h.GetYaxis().SetTitle('#it{R}_{L}')
              setattr(self, name, h)

              # Matched det histograms
              name = 'h_matched_{}{}_JetPt_R{}_{}'.format(observable, pair_type_label, jetR, obs_label)
              h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
              h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
              h.GetYaxis().SetTitle('#it{R}_{L}')
              setattr(self, name, h)

              # Matched det histograms (with matched truth jet pT filled to the other axis)
              name = 'h_matched_extra_{}{}_JetPt_R{}_{}'.format(observable, pair_type_label, jetR, obs_label)
              h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
              h.GetXaxis().SetTitle('#it{p}_{T,ch jet}^{truth}')
              h.GetYaxis().SetTitle('#it{R}_{L}')
              setattr(self, name, h)

              # Matched truth histograms
              name = 'h_matched_{}{}_JetPt_Truth_R{}_{}'.format(observable, pair_type_label, jetR, obs_label)
              h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
              h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
              h.GetYaxis().SetTitle('#it{R}_{L}')
              setattr(self, name, h)

              if self.do_jetcone:
                for jetcone_R in self.jetcone_R_list:
                  # Matched det histograms
                  name = 'h_jetcone{}_matched_{}{}_JetPt_R{}_{}'.format(jetcone_R, observable, pair_type_label, jetR, obs_label)
                  h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
                  h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
                  h.GetYaxis().SetTitle('#it{R}_{L}')
                  setattr(self, name, h)

                  # Matched det histograms (with matched truth jet pT filled to the other axis)
                  name = 'h_jetcone{}_matched_extra_{}{}_JetPt_R{}_{}'.format(jetcone_R, observable, pair_type_label, jetR, obs_label)
                  h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
                  h.GetXaxis().SetTitle('#it{p}_{T,ch jet}^{truth}')
                  h.GetYaxis().SetTitle('#it{R}_{L}')
                  setattr(self, name, h)

                  # Matched truth histograms
                  name = 'h_jetcone{}_matched_{}{}_JetPt_Truth_R{}_{}'.format(jetcone_R, observable, pair_type_label, jetR, obs_label)
                  h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
                  h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
                  h.GetYaxis().SetTitle('#it{R}_{L}')
                  setattr(self, name, h)
              
              if self.thermal_model:
                for R_max in self.max_distance:
                  name = 'h_{}{}_JetPt_R{}_{}_Rmax{}'.format(observable, pair_type_label, jetR, obs_label, R_max)
                  h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
                  h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
                  h.GetYaxis().SetTitle('#it{R}_{L}')
                  setattr(self, name, h)
        
              

        if 'jet_pt' in observable:

          pt_bins = linbins(0,200,200)
          Nconst_bins = linbins(0,50,50)

          # Truth histograms
          name = 'h_1D{}_JetPt_Truth_R{}_{}'.format(observable, jetR, obs_label)
          h = ROOT.TH1D(name, name, 200, pt_bins)
          h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
          h.GetYaxis().SetTitle('Counts')
          setattr(self, name, h)

          name = 'h_Nconst_JetPt_Truth_R{}_{}'.format(jetR, obs_label)
          h = ROOT.TH2D(name, name, 200, pt_bins, 50, Nconst_bins)
          h.GetXaxis().SetTitle('p_{T,ch jet}')
          h.GetYaxis().SetTitle('N_{const}')
          setattr(self, name, h)
          
          name = 'tn_JETINFO{}_Truth_R{}_{}'.format(observable, jetR, obs_label)
          tn = ROOT.TNtuple(name, name, "event_id:jet_id:jet_num_in_ev:jet_pt:total_num_const:num_const_aftercut:corr_rc:leading_q:subleading_q:total_num_baryons:num_baryons_aftercut:total_num_mesons:num_mesons_aftercut")
          setattr(self, name, tn)

          # Det histograms
          name = 'h_1D{}_JetPt_R{}_{}'.format(observable, jetR, obs_label)
          h = ROOT.TH1D(name, name, 200, pt_bins)
          h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
          h.GetYaxis().SetTitle('Counts')
          setattr(self, name, h)

          name = 'h_Nconst_JetPt_R{}_{}'.format(jetR, obs_label)
          h = ROOT.TH2D(name, name, 200, pt_bins, 50, Nconst_bins)
          h.GetXaxis().SetTitle('p_{T,ch jet}')
          h.GetYaxis().SetTitle('N_{const}')
          setattr(self, name, h)

          name = 'tn_JETINFO{}_R{}_{}'.format(observable, jetR, obs_label)
          tn = ROOT.TNtuple(name, name, "event_id:jet_id:jet_num_in_ev:jet_pt:total_num_const:num_const_aftercut:corr_rc:leading_q:subleading_q:total_num_baryons:num_baryons_aftercut:total_num_mesons:num_mesons_aftercut")
          setattr(self, name, tn)


          # Matched det histograms
          name = 'h_matched_1D{}_JetPt_R{}_{}'.format(observable, jetR, obs_label)
          pt_bins = linbins(0,200,200)
          h = ROOT.TH1D(name, name, 200, pt_bins)
          h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
          h.GetYaxis().SetTitle('Counts')
          setattr(self, name, h)

          # Matched truth histograms
          name = 'h_matched_1D{}_JetPt_Truth_R{}_{}'.format(observable, jetR, obs_label)
          pt_bins = linbins(0,200,200)
          h = ROOT.TH1D(name, name, 200, pt_bins)
          h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
          h.GetYaxis().SetTitle('Counts')
          setattr(self, name, h)

          # Correlation between matched det and truth
          name = 'h_matched_1D{}_JetPt_Truth_vs_Det_R{}_{}'.format(observable, jetR, obs_label)
          pt_bins = linbins(0,200,200)
          h = ROOT.TH2D(name, name, 200, pt_bins, 200, pt_bins)
          h.GetXaxis().SetTitle('#it{p}_{T,ch jet}^{det}')
          h.GetYaxis().SetTitle('#it{p}_{T,ch jet}^{truth}')
          setattr(self, name, h)



          #Make some blank arrays to be filled if thnsparse - TODO: idk if this is in the right place...
          self.fsparsepartonJetvalue = array.array( 'd', ( 0, 0, 0, 0, 0 ))
          self.fsparsejetlevelJetvalue = array.array( 'd', ( 0, 0, 0, 0 ))

          

        
        

        # Init pair distance histograms (both det and truth level)
        # average track pt bins
        self.trk_pt_lo = [0, 1, 2, 3, 5, 7, 10]
        self.trk_pt_hi = [1, 2, 3, 5, 7, 10, 100]
        # track pt asymmetry bins: (pt_trk1-pt_trk2)/(pt_trk1+pt_trk2)
        self.trk_alpha_lo = [0, 0.2, 0.4, 0.6, 0.8]
        self.trk_alpha_hi = [0.2, 0.4, 0.6, 0.8, 1]


        '''
        # Residuals and responses (currently not filled or used)
        for trk_thrd in self.obs_settings[observable]:
        
          for ipoint in range(2, 3):
            if not self.is_pp:
              for R_max in self.max_distance:
                self.create_response_histograms(observable, ipoint, jetR, trk_thrd, R_max)          
            else:
              self.create_response_histograms(observable, ipoint, jetR, trk_thrd)
        '''

        #correlation histograms
        # self.new_observables = ["deltap", "deltapt", "charge", "unweightedRL"]
        # print("OBS_LABEL", obs_label)
        if "corr" in observable:
          if observable == "corr_beg":
            self.tuple_obs_string = "event_id:jet_id:jet_num_in_ev:jet_pt:RL:weights"   
          elif observable == "corr_end": #only purpose of this observable is to signify the end
            name = 'tn_pairlevel_Truth_R{}_{}'.format(jetR, obs_label)
            tn = ROOT.TNtuple(name, name, self.tuple_obs_string)
            setattr(self, name, tn)
            colon_count = self.tuple_obs_string.count(':')
            self.fsparsepartonJetvalue_tuple = array.array( 'd', np.zeros(colon_count+1)) # >=18 to match the number of axes
          
            name = 'tn_pairlevel_R{}_{}'.format(jetR, obs_label)
            tn = ROOT.TNtuple(name, name, self.tuple_obs_string)
            setattr(self, name, tn)
            self.fsparsepartonJetvalue_tuple = array.array( 'd', np.zeros(colon_count+1))
          else:
            self.create_corr_tuples(observable, jetR, obs_label)

          

  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  # Fill 2D histogram of (pt, obs)
  #---------------------------------------------------------------
  def create_response_histograms(self, observable, ipoint, jetR, trk_thrd, R_max = None):
  
    if R_max:
      suffix = '_Rmax{}'.format(R_max)
    else:
      suffix = ''

    # # Create THn of response for ENC
    # dim = 4;
    # title = ['#it{p}_{T,det}', '#it{p}_{T,truth}', '#it{R}_{L,det}', '#it{R}_{L,truth}']
    # nbins = [30, 20, 100, 100]
    # min = [0., 0., 0., 0.]
    # max = [150., 200., 1., 1.]
    # name = 'hResponse_JetPt_{}{}_R{}_{}{}'.format(observable, ipoint, jetR, trk_thrd, suffix)
    # self.create_thn(name, title, dim, nbins, min, max)
    
    # name = 'hResidual_JetPt_{}{}_R{}_{}{}'.format(observable, ipoint, jetR, trk_thrd, suffix)
    # h = ROOT.TH3F(name, name, 20, 0, 200, 100, 0., 1., 200, -2., 2.)
    # h.GetXaxis().SetTitle('#it{p}_{T,truth}')
    # h.GetYaxis().SetTitle('#it{R}_{L}')
    # h.GetZaxis().SetTitle('#frac{#it{R}_{L,det}-#it{R}_{L,truth}}{#it{R}_{L,truth}}')
    # setattr(self, name, h)
  

  
  def create_corr_tuples(self, observable, jetR, obs_label):

    if observable == "corr_energyweights" or observable == "corr_rc":
      return
    
    if observable.startswith("corr_"):
      obs_string = observable[5:] #cut out the "corr_"
      if observable == "corr_charge":
        obs_string = "q1q2"
    

    # tuple_obs_string = getattr(self, 'tuple_obs_string') 
    self.tuple_obs_string += ":" + obs_string

    if observable == "corr_deltap":
      self.tuple_obs_string += ":p1"
      self.tuple_obs_string += ":p2"
    elif observable == "corr_deltapt":
      self.tuple_obs_string += ":pt1"
      self.tuple_obs_string += ":pt2"
    elif observable == "corr_deltapl":
      self.tuple_obs_string += ":pl1"
      self.tuple_obs_string += ":pl2"
    elif observable == "corr_charge":
      self.tuple_obs_string += ":q1"
      self.tuple_obs_string += ":q2"
    elif observable == "corr_baryonmeson":
      self.tuple_obs_string += ":pid1"
      self.tuple_obs_string += ":pid2"


    
    if ("baryonmeson" in observable):
      name = 'tn_{}_JetPt_Truth_R{}_{}'.format("baryon", jetR, obs_label)
      tn_b = ROOT.TNtuple(name, name, "jet_pt:baryon_pt")
      setattr(self, name, tn_b)
      name = 'tn_{}_JetPt_Truth_R{}_{}'.format("meson", jetR, obs_label)
      tn_m = ROOT.TNtuple(name, name, "jet_pt:meson_pt")
      setattr(self, name, tn_m)



  def get_pair_eff_weights(self, corr_builder, ipoint, constituents):
    # NB: currently applying the pair eff weight to both 2 point correlator and higher point correlators. Need to check if the same pair efficiency effect still work well for higher point correlators
    weights_pair = []
    for index in range(corr_builder.correlator(ipoint).rs().size()):
      part1 = int(corr_builder.correlator(ipoint).indices1()[index])
      part2 = int(corr_builder.correlator(ipoint).indices2()[index])
      if part1!=part2: # FIX ME: not sure, but for now only apply pair efficiency for non auto-correlations
        # Need to find the associated truth information for each pair (charge and momentum)
        part1_truth = constituents[part1].python_info().particle_truth
        part2_truth = constituents[part2].python_info().particle_truth
        q1 = int(constituents[part1].python_info().charge)
        q2 = int(constituents[part2].python_info().charge)
        dist = corr_builder.correlator(ipoint).rs()[index] # NB: use reconstructed distance since it's faster and should be equivalent to true distance because there is no angular smearing on the track momentum. To switch back to the true distance, use: self.calculate_distance(part1_truth, part2_truth)
        dq_over_p = q1/part1_truth.pt()-q2/part2_truth.pt()
        # calculate pair efficeincy and apply it as an additional weight
        weights_pair.append( self.get_pair_eff(dist, dq_over_p) )
      else:
        weights_pair.append( 1 )
    return weights_pair

  def is_same_charge(self, corr_builder, ipoint, constituents, index):
    part1 = int(corr_builder.correlator(ipoint).indices1()[index])
    part2 = int(corr_builder.correlator(ipoint).indices2()[index])
    q1 = int(constituents[part1].python_info().charge)
    q2 = int(constituents[part2].python_info().charge)

    if q1*q2 > 0:
      return True
    else:
      return False
    
  def charge_p1p2(self, corr_builder, ipoint, constituents, index):
    part1 = int(corr_builder.correlator(ipoint).indices1()[index])
    part2 = int(corr_builder.correlator(ipoint).indices2()[index])
    q1 = int(constituents[part1].python_info().charge)
    q2 = int(constituents[part2].python_info().charge)

    return q1, q2
    

  def mom_p1p2(self, corr_builder, ipoint, constituents, index):
    part1 = int(corr_builder.correlator(ipoint).indices1()[index])
    part2 = int(corr_builder.correlator(ipoint).indices2()[index])
    p1x = constituents[part1].px()
    p1y = constituents[part1].py()
    p1z = constituents[part1].pz()
    p2x = constituents[part2].px()
    p2y = constituents[part2].py()
    p2z = constituents[part2].pz()
    p1 = np.sqrt(p1x*p1x + p1y*p1y + p1z*p1z)
    p2 = np.sqrt(p2x*p2x + p2y*p2y + p2z*p2z)

    return p1, p2
  
  def transmom_p1p2(self, corr_builder, ipoint, constituents, index):
    part1 = int(corr_builder.correlator(ipoint).indices1()[index])
    part2 = int(corr_builder.correlator(ipoint).indices2()[index])
    pt1 = constituents[part1].pt()
    pt2 = constituents[part2].pt()

    return pt1, pt2
  
  def longmom_p1p2(self, corr_builder, ipoint, constituents, index):
    part1 = int(corr_builder.correlator(ipoint).indices1()[index])
    part2 = int(corr_builder.correlator(ipoint).indices2()[index])
    pl1 = constituents[part1].pz()
    pl2 = constituents[part2].pz()

    return pl1, pl2

  # returns:
  # 1 if baryon+baryon (proton + proton)
  # 0 if baryon+meson (proton + pion)
  # -1 if meson+meson (pion + pion)
  # 2 if any other particles
  def is_pair_baryonmeson(self, corr_builder, ipoint, constituents, index):
    part1 = int(corr_builder.correlator(ipoint).indices1()[index])
    part2 = int(corr_builder.correlator(ipoint).indices2()[index])
    pid1 = int(constituents[part1].python_info().particle_pid)
    pid2 = int(constituents[part2].python_info().particle_pid)

    if abs(pid1) == 2212 and abs(pid2) == 2212:
      return 1
    elif abs(pid1) == 211 and abs(pid2) == 211:
      return -1
    elif (abs(pid1) == 2212 and abs(pid2) == 211) or (abs(pid1) == 211 and abs(pid2) == 2212):
      return 0
    else:
      return 2
    
    
  def charge_bm1bm2(self, corr_builder, ipoint, constituents, index):
    part1 = int(corr_builder.correlator(ipoint).indices1()[index])
    part2 = int(corr_builder.correlator(ipoint).indices2()[index])
    pid1 = int(constituents[part1].python_info().particle_pid)
    pid2 = int(constituents[part2].python_info().particle_pid)

    return pid1, pid2
  
  def get_leadsublead_q1q2(self, constituents):
    lead_q = constituents[0].python_info().charge
    if (len(constituents) >= 2):
      sublead_q = constituents[1].python_info().charge
      return lead_q * sublead_q
    else:
      return -99
    
  
  def check_pair_type(self, corr_builder, ipoint, constituents, index):
    part1 = int(corr_builder.correlator(ipoint).indices1()[index])
    part2 = int(corr_builder.correlator(ipoint).indices2()[index])
    type1 = constituents[part1].user_index()
    type2 = constituents[part2].user_index()

    # NB: match the strings in self.pair_type_label = ['bb','sb','ss']
    if type1*type2 >= 0:
      if type1 < 0 or type2 < 0:
        # print('bkg-bkg (',type1,type2,') pt1',constituents[part1].perp(),'pt2',constituents[part2].perp())
        return 0 # means bkg-bkg
      else:
        # print('sig-sig (',type1,type2,') pt1',constituents[part1].perp(),'pt2',constituents[part2].perp())
        return 2 # means sig-sig
    else:
      # print('sig-bkg (',type1,type2,') pt1',constituents[part1].perp(),'pt2',constituents[part2].perp())
      return 1 # means sig-bkg

  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  # Fill 2D histogram of (pt, obs)
  #---------------------------------------------------------------
  def fill_observable_histograms(self, hname, jet, jet_groomed_lund, jetR, obs_setting,
                                 grooming_setting, obs_label, jet_pt_ungroomed):
    
    # For ENC in PbPb, jet_pt_ungroomed stores the corrected jet pT
    constituents = fj.sorted_by_pt(jet.constituents())
    c_select = fj.vectorPJ()
    trk_thrd = obs_setting

    num_baryons_tot = 0
    num_mesons_tot = 0
    num_baryons_aftercut = 0
    num_mesons_aftercut = 0

    for c in constituents:
      pid = c.python_info().particle_pid
      if abs(pid) == 2212: #proton
        num_baryons_tot+=1
      elif abs(pid) == 211: #pion
        num_mesons_tot+=1

      if c.pt() < trk_thrd:
        break
      c_select.append(c) # NB: use the break statement since constituents are already sorted

      if abs(pid) == 2212: #proton
        num_baryons_aftercut+=1
      elif abs(pid) == 211: #pion
        num_mesons_aftercut+=1
    
    nconst_jet = len(c_select)


    if (self.debug_level == 3):
      print("CHARGE and PID HERE")
      print([c.python_info().charge for c in c_select])
      print([c.python_info().particle_pid for c in c_select])

    if self.ENC_pair_cut and (not 'Truth' in hname):
      dphi_cut = -9999 # means no dphi cut
      deta_cut = 0.008
    else:
      dphi_cut = -9999
      deta_cut = -9999

    # NB: use jet_pt_ungroomed instead of jet.perp() for PbPb, which include the UE subtraction
    if self.do_rho_subtraction:
      jet_pt = jet_pt_ungroomed
    else:
      jet_pt = jet.perp()

    # Fill 1D jet pt histogram regardless
    getattr(self, hname.format("1Djet_pt",obs_label)).Fill(jet_pt)


    new_corr = ecorrel.CorrelatorBuilder(c_select, jet_pt, 2, 1, dphi_cut, deta_cut)
    deltap_obs_corr = othercorrel.OtherCorrelatorBuilder(c_select, jet_pt, 2, 1, dphi_cut, deta_cut, "deltap")
    deltapt_obs_corr = othercorrel.OtherCorrelatorBuilder(c_select, jet_pt, 2, 1, dphi_cut, deta_cut, "deltapt")
    deltapl_obs_corr = othercorrel.OtherCorrelatorBuilder(c_select, jet_pt, 2, 1, dphi_cut, deta_cut, "deltapl")
    # print("THERE ARE ", new_corr.correlator(2).rs().size(), "NUM OF ENTRIES IN NEW CORR")
    
    # save jet pt - because it is a jet quantity (not pair)
    self.fsparsepartonJetvalue_tuple[0] = self.event_number
    self.fsparsepartonJetvalue_tuple[1] = self.jet_number   
    self.fsparsepartonJetvalue_tuple[2] = self.ijet    
    self.fsparsepartonJetvalue_tuple[3] = jet_pt
        

    for observable in self.observable_list:
      # print("CP OBSERVABLE", observable)
      
      if 'ENC' in observable or 'EEC_noweight' in observable:
        for ipoint in range(2, 3):

          if self.ENC_fastsim and (not 'Truth' in hname): # NB: only apply pair efficiency effect for fast sim and det level distributions
            weights_pair = self.get_pair_eff_weights(new_corr, ipoint, c_select)

          for index in range(new_corr.correlator(ipoint).rs().size()):

            # processing only like-sign pairs when self.ENC_pair_like is on
            if self.ENC_pair_like and (not self.is_same_charge(new_corr, ipoint, c_select, index)):
              continue

            # processing only unlike-sign pairs when self.ENC_pair_unlike is on
            if self.ENC_pair_unlike and self.is_same_charge(new_corr, ipoint, c_select, index):
              continue
            
            # separate out sig-sig, sig-bkg, bkg-bkg correlations for EEC pairs
            pair_type_label = ''
            if self.do_rho_subtraction or self.do_constituent_subtraction:
              pair_type = self.check_pair_type(new_corr, ipoint, c_select, index)
              pair_type_label = self.pair_type_labels[pair_type]

            if 'ENC' in observable:
              # Det histograms
              if self.ENC_fastsim and (not 'Truth' in hname):
                getattr(self, hname.format(observable + str(ipoint) + pair_type_label,obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index]*weights_pair[index])
                getattr(self, hname.format(observable + str(ipoint) + pair_type_label + 'Pt',obs_label)).Fill(jet_pt, jet_pt*new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index]*weights_pair[index]) # NB: fill pt*RL

              else: # Truth histogram
                getattr(self, hname.format(observable + str(ipoint) + pair_type_label,obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index])
                
                getattr(self, hname.format(observable + str(ipoint) + pair_type_label + 'Pt',obs_label)).Fill(jet_pt, jet_pt*new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index])


            if ipoint==2 and 'EEC_noweight' in observable:
              # Det level
              if self.ENC_fastsim and (not 'Truth' in hname):
                getattr(self, hname.format(observable + pair_type_label,obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], weights_pair[index])
              else: # Truth level
                getattr(self, hname.format(observable + pair_type_label,obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index])


      if 'jet_pt' in observable:

        leadq_subleadq = self.get_leadsublead_q1q2(c_select)
        leading_q = constituents[0].python_info().charge
        if (len(constituents) >= 2):
          subleading_q = constituents[1].python_info().charge
        else:
          subleading_q = -99

        # Det level
        if self.ENC_fastsim and (not 'Truth' in hname):
          getattr(self, hname.format("1D"+observable,obs_label)).Fill(jet_pt, weights_pair[index])
          getattr(self, hname.format('Nconst', obs_label)).Fill(jet_pt, nconst_jet, weights_pair[index]) 
          
          # TODO: don't do this until you've figured out how...
          # getattr(self, 'tn_JETINFOjet_pt_R{}_{}'.format(jetR, obs_label)).Fill(self.event_number, self.jet_number, self.ijet, jet_pt*weights_pair[index], len(constituents), nconst_jet, leadq_subleadq, leading_q, subleading_q, num_baryons_tot, num_baryons_aftercut, num_mesons_tot, num_mesons_aftercut)

        else: # Truth level
          getattr(self, hname.format("1D"+observable,obs_label)).Fill(jet_pt)
          getattr(self, hname.format('Nconst', obs_label)).Fill(jet_pt, nconst_jet) 

          getattr(self, 'tn_JETINFOjet_pt_Truth_R{}_{}'.format(jetR, obs_label)).Fill(self.event_number, self.jet_number, self.ijet, jet_pt, len(constituents), nconst_jet, leadq_subleadq, leading_q, subleading_q, num_baryons_tot, num_baryons_aftercut, num_mesons_tot, num_mesons_aftercut)


        

    # now save pair level dEEC information
    for index in range(new_corr.correlator(ipoint).rs().size()):
      for observable in self.observable_list:
        # print("CP OBSERVABLE", observable)

        # Fill pair level tuple here
        if (observable == "jet_EEC_noweight_RL"):
          continue
        if (observable == "jet_pt"):
          continue
        if (observable == "corr_beg" or observable == "corr_rc"):
          continue

        elif 'ENC' in observable:
          # Fill the RL values here
          self.fsparsepartonJetvalue_tuple[4] = new_corr.correlator(ipoint).rs()[index]
          # print("indices 1:",self.fsparsepartonJetvalue_tuple[1] )
        
        elif "energyweights" in observable:
          self.fsparsepartonJetvalue_tuple[5] = new_corr.correlator(ipoint).weights()[index]
          # print("indices 2:",self.fsparsepartonJetvalue_tuple[2] )

        elif observable == 'corr_deltap':
          self.fsparsepartonJetvalue_tuple[6] = deltap_obs_corr.correlator(ipoint).rs()[index]
          p1, p2 = self.mom_p1p2(new_corr, ipoint, c_select, index)
          self.fsparsepartonJetvalue_tuple[7] = p1
          self.fsparsepartonJetvalue_tuple[8] = p2
          # print("indices 3, 4, 5:",self.fsparsepartonJetvalue_tuple[3], self.fsparsepartonJetvalue_tuple[4], self.fsparsepartonJetvalue_tuple[5] )

        elif observable == 'corr_deltapt':
          self.fsparsepartonJetvalue_tuple[9] = deltapt_obs_corr.correlator(ipoint).rs()[index]
          pt1, pt2 = self.transmom_p1p2(new_corr, ipoint, c_select, index)
          self.fsparsepartonJetvalue_tuple[10] = pt1
          self.fsparsepartonJetvalue_tuple[11] = pt2
          # print("indices 6, 7, 8:",self.fsparsepartonJetvalue_tuple[6], self.fsparsepartonJetvalue_tuple[7], self.fsparsepartonJetvalue_tuple[8] )
          
        elif observable == 'corr_deltapl':
          self.fsparsepartonJetvalue_tuple[12] = deltapl_obs_corr.correlator(ipoint).rs()[index]
          pl1, pl2 = self.longmom_p1p2(new_corr, ipoint, c_select, index)
          self.fsparsepartonJetvalue_tuple[13] = pl1
          self.fsparsepartonJetvalue_tuple[14] = pl2
          # print("indices 9, 10, 11:",self.fsparsepartonJetvalue_tuple[9], self.fsparsepartonJetvalue_tuple[10], self.fsparsepartonJetvalue_tuple[11] )

        elif observable == 'corr_charge':
          samecharge_boolean = self.is_same_charge(new_corr, ipoint, c_select, index)
          self.fsparsepartonJetvalue_tuple[15] = 1 if samecharge_boolean else -1
          q1, q2 = self.charge_p1p2(new_corr, ipoint, c_select, index)
          self.fsparsepartonJetvalue_tuple[16] = q1
          self.fsparsepartonJetvalue_tuple[17] = q2
          # print("indices 12, 13, 14:",self.fsparsepartonJetvalue_tuple[12], self.fsparsepartonJetvalue_tuple[13], self.fsparsepartonJetvalue_tuple[14] )

        elif ("baryonmeson" in observable):
          baryonmeson_quantity = self.is_pair_baryonmeson(new_corr, ipoint, c_select, index)
          pid1, pid2 = self.charge_bm1bm2(new_corr, ipoint, c_select, index)
          self.fsparsepartonJetvalue[18] = baryonmeson_quantity
          self.fparsepartonJetvalue[19] = pid1
          self.fparsepartonJetvalue[20] = pid2

          baryon_tn_name = hname.format("baryon",obs_label)
          baryon_tn_name = "tn" + baryon_tn_name[1:]
          meson_tn_name = hname.format("meson",obs_label)
          meson_tn_name = "tn" + meson_tn_name[1:]
          for c in constituents:
            pid = c.python_info().particle_pid
            if abs(pid) == 2212: #proton
              getattr(self, baryon_tn_name).Fill(jet_pt, c.pt())
            elif abs(pid) == 211: #pion
              getattr(self, meson_tn_name).Fill(jet_pt, c.pt())
            
        elif 'end' in observable:
          # print("IN CORR_END WRITING", self.fsparsepartonJetvalue_tuple)
          # print("IN CORR_END WRITING", self.fsparsepartonJetvalue_tuple[0])
          #TODO: again, don't do det level until you know what you're doing
          # if self.ENC_fastsim and (not 'Truth' in hname):
          #   getattr(self, 'tn_pairlevel_R{}_{}'.format(jetR, obs_label)).Fill(np.array(list(self.fsparsepartonJetvalue_tuple), dtype='float32'))
          # else:
          # print("IN END!!!!", self.fsparsepartonJetvalue_tuple )
          getattr(self, 'tn_pairlevel_Truth_R{}_{}'.format(jetR, obs_label)).Fill(np.array(list(self.fsparsepartonJetvalue_tuple), dtype='float32'))

        


  #---------------------------------------------------------------
  # This function is called per observable per jet subconfigration 
  # used in fill_matched_jet_histograms
  # This function is created because we cannot use fill_observable_histograms 
  # directly because observable list loop inside that function
  #---------------------------------------------------------------
  #TODO: fix this later!

  def fill_matched_observable_histograms(self, hname, observable, jet, jet_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_ungroomed, jet_pt_matched, cone_parts = None):
    
    constituents = fj.sorted_by_pt(jet.constituents())
    if cone_parts!=None:
      constituents = fj.sorted_by_pt(cone_parts)

    # if cone_parts!=None and 'Truth' in hname:
    #   print('Nconst in cone',len(constituents))
    # if cone_parts==None and 'Truth' in hname:
    #   print('Nconst in jet',len(constituents))

    c_select = fj.vectorPJ()
    trk_thrd = obs_setting

    for c in constituents:
      if c.pt() < trk_thrd:
        break
      c_select.append(c) # NB: use the break statement since constituents are already sorted
    
    if self.ENC_pair_cut and (not 'Truth' in hname):
      dphi_cut = -9999 # means no dphi cut
      deta_cut = 0.008
    else:
      dphi_cut = -9999
      deta_cut = -9999

    # Only need rho subtraction for det-level jets 
    if self.do_rho_subtraction and (not 'Truth' in hname):
      jet_pt = jet_pt_ungroomed
    else:
      jet_pt = jet.perp()

    new_corr = ecorrel.CorrelatorBuilder(c_select, jet_pt, 2, 1, dphi_cut, deta_cut)
    if 'ENC' in observable or 'EEC_noweight' in observable or 'EEC_weight2' in observable:
      for ipoint in range(2, 3):
        if self.ENC_fastsim and (not 'Truth' in hname): # NB: only apply pair efficiency effect for fast sim and det level distributions
          weights_pair = self.get_pair_eff_weights(new_corr, ipoint, c_select)

        for index in range(new_corr.correlator(ipoint).rs().size()):
          pair_type_label = ''
          if self.do_rho_subtraction or self.do_constituent_subtraction:
            pair_type = self.check_pair_type(new_corr, ipoint, c_select, index)
            pair_type_label = self.pair_type_labels[pair_type]
          
          if 'ENC' in observable:
            if self.ENC_fastsim and (not 'Truth' in hname):
              getattr(self, hname.format(observable + str(ipoint) + pair_type_label,obs_label)).Fill(jet_pt_matched, new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index]*weights_pair[index])
            else:
              getattr(self, hname.format(observable + str(ipoint) + pair_type_label,obs_label)).Fill(jet_pt_matched, new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index]) # NB: use jet_pt_matched instead of jet_pt so if jet_pt_matched is different from jet_pt, it will be used. This is mainly for matched jets study

          if ipoint==2 and 'EEC_noweight' in observable:
            if self.ENC_fastsim and (not 'Truth' in hname):
              getattr(self, hname.format(observable + pair_type_label,obs_label)).Fill(jet_pt_matched, new_corr.correlator(ipoint).rs()[index], weights_pair[index])
            else:
              getattr(self, hname.format(observable + pair_type_label,obs_label)).Fill(jet_pt_matched, new_corr.correlator(ipoint).rs()[index])

          if ipoint==2 and 'EEC_weight2' in observable:
            if self.ENC_fastsim and (not 'Truth' in hname):
              getattr(self, hname.format(observable + pair_type_label,obs_label)).Fill(jet_pt_matched, new_corr.correlator(ipoint).rs()[index], pow(new_corr.correlator(ipoint).weights()[index]*weights_pair[index],2))
            else:
              getattr(self, hname.format(observable + pair_type_label,obs_label)).Fill(jet_pt_matched, new_corr.correlator(ipoint).rs()[index], pow(new_corr.correlator(ipoint).weights()[index],2))

    if 'jet_pt' in observable:
      getattr(self, hname.format(observable,obs_label)).Fill(jet_pt)

  #---------------------------------------------------------------
  # This function is called per jet subconfigration 
  # Fill matched jet histograms
  #---------------------------------------------------------------
  def fill_matched_jet_histograms(self, jet_det, jet_det_groomed_lund, jet_truth,
                                  jet_truth_groomed_lund, jet_pp_det, jetR,
                                  obs_setting, grooming_setting, obs_label,
                                  jet_pt_det_ungroomed, jet_pt_truth_ungroomed, R_max, suffix, **kwargs):
    # If jetscape, we will need to correct substructure observable for holes (pt is corrected in base class)
    # For ENC in PbPb, jet_pt_det_ungroomed stores the corrected jet pT
    if self.jetscape:
      holes_in_det_jet = kwargs['holes_in_det_jet']
      holes_in_truth_jet = kwargs['holes_in_truth_jet']

    cone_parts_in_det_jet = kwargs['cone_parts_in_det_jet']
    cone_parts_in_truth_jet = kwargs['cone_parts_in_truth_jet']
    cone_R = kwargs['cone_R']

    # Todo: add additonal weight for jet pT spectrum
    # if self.rewight_pt:
    #   w_pt = 1+pow(jet_truth,0.2)
    # else:
    #   w_pt = 1
    
    if self.do_rho_subtraction:
      # print('evt #',self.event_number)
      jet_pt_det = jet_pt_det_ungroomed
      # print('Det: pT',jet_det.perp(),'(',jet_pt_det,')','phi',jet_det.phi(),'eta',jet_det.eta())
      # print('Truth: pT',jet_truth.perp(),'phi',jet_truth.phi(),'eta',jet_truth.eta())
      # print('Difference pT (truth-det)',jet_truth.perp()-jet_pt_det_ungroomed)
    else:
      jet_pt_det = jet_det.perp()

    for observable in self.observable_list:
      
      if cone_R == 0: # fill for jet constituents
        hname = 'h_matched_{{}}_JetPt_R{}_{{}}'.format(jetR)
        self.fill_matched_observable_histograms(hname, observable, jet_det, jet_det_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_det, jet_pt_det)

        hname = 'h_matched_{{}}_JetPt_Truth_R{}_{{}}'.format(jetR)
        self.fill_matched_observable_histograms(hname, observable, jet_truth, jet_truth_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_det, jet_truth.pt())

        # fill RL vs matched truth jet pT for det jets (only fill these extra histograms for ENC or pair distributions)
        if 'ENC' in observable or 'EEC_noweight' in observable or 'EEC_weight2' in observable:
          hname = 'h_matched_extra_{{}}_JetPt_R{}_{{}}'.format(jetR)
          self.fill_matched_observable_histograms(hname, observable, jet_det, jet_det_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_det, jet_truth.pt()) # NB: use the truth jet pt so the reco jets histograms are comparable to matched truth jets. However this also means that two identical histograms will be filled fot jet_pt observable

        # Fill correlation between matched det and truth jets
        if 'jet_pt' in observable:
          hname = 'h_matched_{}_JetPt_Truth_vs_Det_R{}_{}'.format(observable, jetR, obs_label)
          getattr(self, hname).Fill(jet_pt_det, jet_truth.pt())

      else: # fill for cone parts around jet
        if 'ENC' in observable or 'EEC_noweight' in observable or 'EEC_weight2' in observable:
          hname = 'h_jetcone{}_matched_{{}}_JetPt_R{}_{{}}'.format(cone_R, jetR)
          self.fill_matched_observable_histograms(hname, observable, jet_det, jet_det_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_det, jet_pt_det, cone_parts_in_det_jet)

          hname = 'h_jetcone{}_matched_{{}}_JetPt_Truth_R{}_{{}}'.format(cone_R, jetR)
          self.fill_matched_observable_histograms(hname, observable, jet_truth, jet_truth_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_det, jet_truth.pt(), cone_parts_in_truth_jet)

          hname = 'h_jetcone{}_matched_extra_{{}}_JetPt_R{}_{{}}'.format(cone_R, jetR)
          self.fill_matched_observable_histograms(hname, observable, jet_det, jet_det_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_det, jet_truth.pt(), cone_parts_in_det_jet)          


    # # Find all subjets
    # trk_thrd = obs_setting
    # cs_subjet_det = fj.ClusterSequence(jet_det.constituents(), self.subjet_def[trk_thrd])
    # subjets_det = fj.sorted_by_pt(cs_subjet_det.inclusive_jets())

    # cs_subjet_truth = fj.ClusterSequence(jet_truth.constituents(), self.subjet_def[trk_thrd])
    # subjets_truth = fj.sorted_by_pt(cs_subjet_truth.inclusive_jets())
    
    # if not self.is_pp:
    #   cs_subjet_det_pp = fj.ClusterSequence(jet_pp_det.constituents(), self.subjet_def[trk_thrd])
    #   subjets_det_pp = fj.sorted_by_pt(cs_subjet_det_pp.inclusive_jets())

    # # Loop through subjets and set subjet matching candidates for each subjet in user_info
    # if self.is_pp:
    #     [[self.set_matching_candidates(subjet_det, subjet_truth, subjetR, 'hDeltaR_ppdet_pptrue_ENC_R{}_{}'.format(jetR, subjetR)) for subjet_truth in subjets_truth] for subjet_det in subjets_det]
    # else:
    #     # First fill the combined-to-pp matches, then the pp-to-pp matches
    #     [[self.set_matching_candidates(subjet_det_combined, subjet_det_pp, subjetR, 'hDeltaR_combined_ppdet_ENC_R{}_{}_Rmax{}'.format(jetR, subjetR, R_max), fill_jet1_matches_only=True) for subjet_det_pp in subjets_det_pp] for subjet_det_combined in subjets_det]
    #     [[self.set_matching_candidates(subjet_det_pp, subjet_truth, subjetR, 'hDeltaR_ppdet_pptrue_ENC_R{}_{}_Rmax{}'.format(jetR, subjetR, R_max)) for subjet_truth in subjets_truth] for subjet_det_pp in subjets_det_pp]
      
    # # Loop through subjets and set accepted matches
    # if self.is_pp:
    #     [self.set_matches_pp(subjet_det, 'hSubjetMatchingQA_R{}_{}'.format(jetR, subjetR)) for subjet_det in subjets_det]
    # else:
    #     [self.set_matches_AA(subjet_det_combined, subjetR, 'hSubjetMatchingQA_R{}_{}'.format(jetR, subjetR)) for subjet_det_combined in subjets_det]

    # # Loop through matches and fill histograms
    # for observable in self.observable_list:
    
    #   # Fill inclusive subjets
    #   if 'inclusive' in observable:

    #     for subjet_det in subjets_det:
        
    #       z_det = subjet_det.pt() / jet_det.pt()
          
    #       # If z=1, it will be default be placed in overflow bin -- prevent this
    #       if np.isclose(z_det, 1.):
    #         z_det = 0.999
          
    #       successful_match = False

    #       if subjet_det.has_user_info():
    #         subjet_truth = subjet_det.python_info().match
          
    #         if subjet_truth:
            
    #           successful_match = True

    #           # For subjet matching radius systematic, check distance between subjets
    #           if self.matching_systematic:
    #             if subjet_det.delta_R(subjet_truth) > 0.5 * self.jet_matching_distance * subjetR:
    #               continue
              
    #           z_truth = subjet_truth.pt() / jet_truth.pt()
              
    #           # If z=1, it will be default be placed in overflow bin -- prevent this
    #           if np.isclose(z_truth, 1.):
    #             z_truth = 0.999
              
    #           # In Pb-Pb case, fill matched pt fraction
    #           if not self.is_pp:
    #             self.fill_subjet_matched_pt_histograms(observable,
    #                                                    subjet_det, subjet_truth,
    #                                                    z_det, z_truth,
    #                                                    jet_truth.pt(), jetR, subjetR, R_max)
              
    #           # Fill histograms
    #           # Note that we don't fill 'matched' histograms here, since that is only
    #           # meaningful for leading subjets
    #           self.fill_response(observable, jetR, jet_pt_det_ungroomed, jet_pt_truth_ungroomed,
    #                              z_det, z_truth, obs_label, R_max, prong_match=False)
                                 
    #       # Fill number of subjets with/without unique match, as a function of zr
    #       if self.is_pp:
    #         name = 'h_match_fraction_{}_R{}_{}'.format(observable, jetR, subjetR)
    #         getattr(self, name).Fill(jet_truth.pt(), z_det, successful_match)
                               
    #   # Get leading subjet and fill histograms
    #   if 'leading' in observable:
      
    #     leading_subjet_det = self.utils.leading_jet(subjets_det)
    #     leading_subjet_truth = self.utils.leading_jet(subjets_truth)
        
    #     # Note that we don't want to check whether they are geometrical matches
    #     # We rather want to correct the measured leading subjet to the true leading subjet
    #     if leading_subjet_det and leading_subjet_truth:
          
    #       z_leading_det = leading_subjet_det.pt() / jet_det.pt()
    #       z_leading_truth = leading_subjet_truth.pt() / jet_truth.pt()
          
    #       # If z=1, it will be default be placed in overflow bin -- prevent this
    #       if np.isclose(z_leading_det, 1.):
    #         z_leading_det = 0.999
    #       if np.isclose(z_leading_truth, 1.):
    #         z_leading_truth = 0.999
          
    #       # In Pb-Pb case, fill matched pt fraction
    #       if not self.is_pp:
    #         match = self.fill_subjet_matched_pt_histograms(observable,
    #                                                        leading_subjet_det, leading_subjet_truth,
    #                                                        z_leading_det, z_leading_truth,
    #                                                        jet_truth.pt(), jetR, subjetR, R_max)
    #       else:
    #         match = False
          
    #       # Fill histograms
    #       self.fill_response(observable, jetR, jet_pt_det_ungroomed, jet_pt_truth_ungroomed,
    #                          z_leading_det, z_leading_truth, obs_label, R_max, prong_match=match)
          
    #       # Plot deltaR distribution between the detector and truth leading subjets
    #       # (since they are not matched geometrically, the true leading may not be the measured leading
    #       deltaR = leading_subjet_det.delta_R(leading_subjet_truth)
    #       name = 'hDeltaR_det_truth_{}_R{}_{}'.format(observable, jetR, subjetR)
    #       if not self.is_pp:
    #         name += '_Rmax{}'.format(R_max)
    #       getattr(self, name).Fill(jet_truth.pt(), z_leading_truth, deltaR)
        
  #---------------------------------------------------------------
  # Do prong-matching
  #---------------------------------------------------------------
  # def fill_subjet_matched_pt_histograms(self, observable, subjet_det, subjet_truth,
  #                                       z_det, z_truth, jet_pt_truth, jetR, subjetR, R_max):
    
    # # Get pp det-level subjet
    # # Inclusive case: This is matched to the combined subjet (and its pp truth-level subjet)
    # # Leading case: This is matched only to the pp truth-level leading subjet
    # subjet_pp_det = None
    # if subjet_truth.has_user_info():
    #   subjet_pp_det = subjet_truth.python_info().match
    # if not subjet_pp_det:
    #   return
                                     
    # matched_pt = fjtools.matched_pt(subjet_det, subjet_pp_det)
    # name = 'h_{}_matched_pt_JetPt_R{}_{}_Rmax{}'.format(observable, jetR, subjetR, R_max)
    # getattr(self, name).Fill(jet_pt_truth, z_det, matched_pt)
    
    # # Plot dz between det and truth subjets
    # deltaZ = z_det - z_truth
    # name = 'h_{}_matched_pt_deltaZ_JetPt_R{}_{}_Rmax{}'.format(observable, jetR, subjetR, R_max)
    # getattr(self, name).Fill(jet_pt_truth, matched_pt, deltaZ)

    # # Plot dR between det and truth subjets
    # deltaR = subjet_det.delta_R(subjet_truth)
    # name = 'h_{}_matched_pt_deltaR_JetPt_R{}_{}_Rmax{}'.format(observable, jetR, subjetR, R_max)
    # getattr(self, name).Fill(jet_pt_truth, matched_pt, deltaR)

    # match = (matched_pt > 0.5)
    # return match

##################################################################
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description='Process MC')
  parser.add_argument('-f', '--inputFile', action='store',
                      type=str, metavar='inputFile',
                      default='AnalysisResults.root',
                      help='Path of ROOT file containing TTrees')
  parser.add_argument('-c', '--configFile', action='store',
                      type=str, metavar='configFile',
                      default='config/analysis_config.yaml',
                      help="Path of config file for analysis")
  parser.add_argument('-o', '--outputDir', action='store',
                      type=str, metavar='outputDir',
                      default='./TestOutput',
                      help='Output directory for output to be written to')
  parser.add_argument('-t', '--saveTuples', action='store',
                      type=int, metavar='saveTuples',
                      default=1,
                      help='Save correlations as tuples if 1 (otherwise save as histograms if 0)')
  
  # Parse the arguments
  args = parser.parse_args()
  
  print('Configuring...')
  print('inputFile: \'{0}\''.format(args.inputFile))
  print('configFile: \'{0}\''.format(args.configFile))
  print('ouputDir: \'{0}\"'.format(args.outputDir))
  print('saveTuples: \'{0}\"'.format(args.saveTuples))

  # If invalid inputFile is given, exit
  if not os.path.exists(args.inputFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
    sys.exit(0)
  
  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = ProcessMC_dEEC(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir, save_tuples=args.saveTuples)
  analysis.process_mc()