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

import pythiafjext

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io
from pyjetty.alice_analysis.process.base import process_io_emb
from pyjetty.alice_analysis.process.base import jet_info
from pyjetty.alice_analysis.process.user.substructure import process_mc_base
from pyjetty.alice_analysis.process.base import thermal_generator
from pyjetty.mputils.csubtractor import CEventSubtractor

def linbins(xmin, xmax, nbins):
  lspace = np.linspace(xmin, xmax, nbins+1)
  arr = array.array('f', lspace)
  return arr

def logbins(xmin, xmax, nbins):
  lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
  arr = array.array('f', lspace)
  return arr

################################################################
class ProcessMC_ENC_HF(process_mc_base.ProcessMCBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', event_start_offset=0, debug_level=0, **kwargs):
  
    # Initialize base class
    super(ProcessMC_ENC_HF, self).__init__(input_file, config_file, output_dir, event_start_offset, debug_level, **kwargs)
    
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
                
                dim = 4
                pt_bins = linbins(0,200,200)
                ptRL_bins = logbins(1E-3,1E2,60)
                rapi_bins = np.linspace(-5,5,201)
                RL_bins = logbins(1E-4,1,50)

                # Truth histograms
                name = 'h_{}{}{}_JetPt_Truth_R{}_{}'.format(observable, ipoint, pair_type_label, jetR, obs_label) #pair_type_label is blank for me
                title_truth = [ '#it{p}_{T}^{ch jet}', '#it{p}_{T}^{D^{0}}', 'D^{0} y', '#it{R}_{L}' ]
                binnings = [pt_bins, rapi_bins, RL_bins]
                self.create_thn(name, title_truth, dim, binnings, obs='rl')

                name = 'h_{}{}{}Pt_JetPt_Truth_R{}_{}'.format(observable, ipoint, pair_type_label, jetR, obs_label) # pt scaled histograms (currently only for unmatched jets)
                title = [ '#it{p}_{T}^{ch jet}', '#it{p}_{T}^{D^{0}}', 'D^{0} y', '#it{p}_{T}#it{R}_{L}' ]
                binnings = (pt_bins, rapi_bins, ptRL_bins)
                self.create_thn(name, title_truth, dim, binnings, obs='ptrl') 
                
            if 'EEC_noweight' in observable or 'EEC_weight2' in observable:
              
              dim = 4
              pt_bins = linbins(0,200,200)
              rapi_bins = np.linspace(-5,5,201)
              RL_bins = logbins(1E-4,1,50)

              # Truth histograms
              name = 'h_{}{}_JetPt_Truth_R{}_{}'.format(observable, pair_type_label, jetR, obs_label)
              title = [ '#it{p}_{T}^{ch jet}', '#it{p}_{T}^{D^{0}}', 'D^{0} y', '#it{R}_{L}' ]
              binnings = (pt_bins, rapi_bins, RL_bins)
              self.create_thn(name, title_truth, dim, binnings, obs='rl') # this will need to be fixed later! - add more dimensions to also include jet pt?

        
        if 'jet_pt' in observable:

          dim = 3
          pt_bins = linbins(0,200,200)
          rapi_bins = np.linspace(-5,5,201)

          # Truth histograms
          name = 'h_{}_JetPt_Truth_R{}_{}'.format(observable, jetR, obs_label)
          title = [ '#it{p}_{T}^{ch jet}', '#it{p}_{T}^{D^{0}}', 'D^{0} y' ]
          binnings = (pt_bins, rapi_bins)
          self.create_thn(name, title_truth, dim, binnings) # this will need to be fixed later! - add more dimensions to also include jet pt?

        #Make some blank arrays to be filled if thnsparse
        self.fsparsepartonJetvalue = array.array( 'd', ( 0, 0, 0 ,0 ))
        self.fsparsejetlevelJetvalue = array.array( 'd', ( 0, 0, 0 ))

        # Init pair distance histograms (both det and truth level)
        # average track pt bins
        self.trk_pt_lo = [0, 1, 2, 3, 5, 7, 10]
        self.trk_pt_hi = [1, 2, 3, 5, 7, 10, 100]
        # track pt asymmetry bins: (pt_trk1-pt_trk2)/(pt_trk1+pt_trk2)
        self.trk_alpha_lo = [0, 0.2, 0.4, 0.6, 0.8]
        self.trk_alpha_hi = [0.2, 0.4, 0.6, 0.8, 1]
        if 'EEC_detail' in observable: # not currently being used
          # inclusive
          name = 'h_{}_JetPt_R{}_{}'.format(observable, jetR, obs_label)
          pt_bins = linbins(0,200,200)
          RL_bins = logbins(1E-4,1,50)
          h = ROOT.TH2D(name, name, 50, pt_bins, 50, RL_bins)
          h.GetXaxis().SetTitle('p_{T,ch jet}')
          h.GetYaxis().SetTitle('R_{L}')
          setattr(self, name, h)

          name = 'h_{}_JetPt_Truth_R{}_{}'.format(observable, jetR, obs_label)
          pt_bins = linbins(0,200,200)
          RL_bins = logbins(1E-4,1,50)
          h = ROOT.TH2D(name, name, 50, pt_bins, 50, RL_bins)
          h.GetXaxis().SetTitle('p_{T,ch jet}')
          h.GetYaxis().SetTitle('R_{L}')
          setattr(self, name, h)

          # fine bins
          for ipt in range( len(self.trk_pt_lo) ):
            for ialpha in range( len(self.trk_alpha_lo) ):
              name = 'h_{}{}{}_{:.1f}{:.1f}_JetPt_R{}_{}'.format(observable, self.trk_pt_lo[ipt], self.trk_pt_hi[ipt], self.trk_alpha_lo[ialpha], self.trk_alpha_hi[ialpha], jetR, obs_label)
              pt_bins = linbins(0,200,200)
              RL_bins = logbins(1E-4,1,50)
              h = ROOT.TH2D(name, name, 50, pt_bins, 50, RL_bins)
              h.GetXaxis().SetTitle('p_{T,ch jet}')
              h.GetYaxis().SetTitle('R_{L}')
              setattr(self, name, h)

              name = 'h_{}{}{}_{:.1f}{:.1f}_JetPt_Truth_R{}_{}'.format(observable, self.trk_pt_lo[ipt], self.trk_pt_hi[ipt], self.trk_alpha_lo[ialpha], self.trk_alpha_hi[ialpha], jetR, obs_label)
              pt_bins = linbins(0,200,200)
              RL_bins = logbins(1E-4,1,50)
              h = ROOT.TH2D(name, name, 50, pt_bins, 50, RL_bins)
              h.GetXaxis().SetTitle('p_{T,ch jet}')
              h.GetYaxis().SetTitle('R_{L}')
              setattr(self, name, h)



        #correlation histograms
        # self.new_observables = ["deltap", "deltapt", "charge", "unweightedRL"]
        # print("OBS_LABEL", obs_label)
        if "corr" in observable: # not being used in herwig d0 analysis
          self.create_corr_histograms(observable, ipoint, jetR, obs_label)
          

  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  # Fill 2D histogram of (pt, obs)
  #---------------------------------------------------------------
  def create_corr_histograms(self, observable, ipoint, jetR, obs_label):
  
    pt_bins = linbins(0,200,200)
    # RL_bins = logbins(1E-4,1,50)
    ptRL_bins = logbins(1E-3,1E2,60)
    deltap_bins = linbins(0, 1., 100)
    charge_bins = linbins(-1.5, 1.5, 3)

    a1 = np.logspace(np.log10(1E-4),np.log10(1E-2),26)
    arr1 = array.array('f', a1)
    arr2 = logbins(1E-2,2E-2,5)
    arr3 = logbins(2E-2,3E-2,5)
    arr4 = logbins(3E-2,7E-2,5)
    arr5 = logbins(7E-2,2E-1,5)
    arr6 = logbins(2E-1,4E-1,5)
    arr7 = logbins(4E-1,1,5)
    RL_bins = arr1+arr2+arr3+arr4+arr5+arr6+arr7


    # Create histograms
    # delta p, truth
    dim = 3;
    if (observable == "corr_deltap"):
      title_truth = ['p_{T,ch jet,truth}', 'R_{L,truth}', '#Deltap_{truth}']
      title = ['p_{T,ch jet,det}', 'R_{L,det}', '#Deltap_{det}']
      obs_bins = deltap_bins

    # delta pt, truth
    if (observable == "corr_deltapt"):
      title_truth = ['p_{T,ch jet,truth}', 'R_{L,truth}', '#Deltap_{T, truth}']
      title = ['p_{T,ch jet,det}', 'R_{L,det}', '#Deltap_{T, det}']
      obs_bins = deltap_bins

    # charge
    if (observable == "corr_samecharge"):
      title_truth = ['p_{T,ch jet,truth}', 'R_{L,truth}', 'same charge_{truth}']
      title = ['p_{T,ch jet,det}', 'R_{L,det}', 'same charge_{det}']
      obs_bins = charge_bins
    if (observable == "corr_oppcharge"):
      title_truth = ['p_{T,ch jet,truth}', 'R_{L,truth}', 'opp charge_{truth}']
      title = ['p_{T,ch jet,det}', 'R_{L,det}', 'opp charge_{det}']
      obs_bins = charge_bins

    # unweighted RL
    if (observable == "corr_unweightedRL"):
      title_truth = ['p_{T,ch jet,truth}', 'R_{L,truth}', 'unweighted R_{L,truth}']
      title = ['p_{T,ch jet,det}', 'R_{L,det}', 'unweighted R_{L,det}']
      obs_bins = RL_bins


    # name = 'h_JetPt_{}{}_R{}_{}'.format(observable, ipoint, jetR, obs_label)
    name = 'h_{}{}_JetPt_Truth_R{}_{}'.format(observable, ipoint, jetR, obs_label)
    # print("NAME", name)
    nbins  = [len(pt_bins)-1, len(RL_bins)-1, len(obs_bins)-1]
    min = [pt_bins[0],      RL_bins[0],     obs_bins[0]]
    max = [pt_bins[-1],     RL_bins[-1],    obs_bins[-1]]
    self.create_thn(name, title_truth, dim, nbins, min, max)

    name = 'h_{}{}_JetPt_R{}_{}'.format(observable, ipoint, jetR, obs_label)
    self.create_thn(name, title, dim, nbins, min, max)
    



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

    # NB: use jet_pt_ungroomed instead of jet.perp() for PbPb, which include the UE subtraction
    if self.do_rho_subtraction:
      jet_pt = jet_pt_ungroomed
    else:
      jet_pt = jet.perp()

    new_corr = ecorrel.CorrelatorBuilder(c_select, jet_pt, 2, 1, dphi_cut, deta_cut)

    # save jet level information to thnsparse arrays
    self.fsparsejetlevelJetvalue[0] = jet_pt
    self.fsparsepartonJetvalue[0] = jet_pt
    
    # find the D0's - assuming one D0 per jet
    d0injetfound_bool = False
    for part in c_select:
      if self.findD0(part): # if self.D0particleinfo != None:
        d0injetfound_bool = True
        D0_px = self.D0particleinfo.px()
        D0_py = self.D0particleinfo.py()
        self.fsparsejetlevelJetvalue[1] = math.sqrt(D0_px*D0_px + D0_py*D0_py)
        self.fsparsepartonJetvalue[1] = math.sqrt(D0_px*D0_px + D0_py*D0_py)
        self.fsparsejetlevelJetvalue[2] = self.D0particleinfo.python_info().particle_rap #self.D0particleinfo.rap()
        self.fsparsepartonJetvalue[2] = self.D0particleinfo.python_info().particle_rap #self.D0particleinfo.rap()
        if self.firsttimejet:
          print('event {}'.format(self.event_number))
          print("D0 pt is ", math.sqrt(D0_px*D0_px + D0_py*D0_py), "and D0 rapidity is", self.D0particleinfo.python_info().particle_rap) #self.D0particleinfo.rap())
          self.alld0counter+=1          
        break # break after one D0 found in a jet
    if d0injetfound_bool == False:
      # print("in else NOT a D0 jet IS THIS WRONG")
      return # don't want to look at jets that don't have a D0

    # now check if the D0 comes from a D*. if yes, skip. if no, move on
    # print("checking D*")
    if self.firsttimejet:
      print("D0 mother is", self.D0particleinfo.python_info().particle_mid)
    # print("here", hname)
    if abs(self.D0particleinfo.python_info().particle_mid) == 413: # D*
      return

    if self.firsttimejet:
      self.d0nodstar_counter+=1    
      print("------- saving to hists -------")

    for observable in self.observable_list:
      # print("CP OBSERVABLE", observable)
      if 'ENC' in observable or 'EEC_noweight' in observable or 'EEC_weight2' in observable:
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
              if self.ENC_fastsim and (not 'Truth' in hname):
                getattr(self, hname.format(observable + str(ipoint) + pair_type_label,obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index]*weights_pair[index])
                getattr(self, hname.format(observable + str(ipoint) + pair_type_label + 'Pt',obs_label)).Fill(jet_pt, jet_pt*new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index]*weights_pair[index]) # NB: fill pt*RL

              else:
                # only change truth level hists to thnsparse
                # save in thnsparse
                self.fsparsepartonJetvalue[3] = new_corr.correlator(ipoint).rs()[index]
                getattr(self, hname.format(observable + str(ipoint) + pair_type_label,obs_label)).Fill(self.fsparsepartonJetvalue, new_corr.correlator(ipoint).weights()[index])
                
                self.fsparsepartonJetvalue[3] = jet_pt*new_corr.correlator(ipoint).rs()[index]
                getattr(self, hname.format(observable + str(ipoint) + pair_type_label + 'Pt',obs_label)).Fill(self.fsparsepartonJetvalue, new_corr.correlator(ipoint).weights()[index])

                # getattr(self, hname.format(observable + str(ipoint) + pair_type_label,obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index])
                # getattr(self, hname.format(observable + str(ipoint) + pair_type_label + 'Pt',obs_label)).Fill(jet_pt, jet_pt*new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index])

            if ipoint==2 and 'EEC_noweight' in observable:
              if self.ENC_fastsim and (not 'Truth' in hname):
                getattr(self, hname.format(observable + pair_type_label,obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], weights_pair[index])
              else:
                #change truth level hists to thnsparse
                self.fsparsepartonJetvalue[3] = new_corr.correlator(ipoint).rs()[index]
                getattr(self, hname.format(observable + pair_type_label,obs_label)).Fill(self.fsparsepartonJetvalue)

            if ipoint==2 and 'EEC_weight2' in observable:
              if self.ENC_fastsim and (not 'Truth' in hname):
                getattr(self, hname.format(observable + pair_type_label,obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], pow(new_corr.correlator(ipoint).weights()[index]*weights_pair[index],2))
              else:
                getattr(self, hname.format(observable + pair_type_label,obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], pow(new_corr.correlator(ipoint).weights()[index],2))

      if 'jet_pt' in observable:
        getattr(self, hname.format(observable,obs_label)).Fill(self.fsparsejetlevelJetvalue)
      
      # NB: for now, only perform this check on data and full sim
      if 'EEC_detail' in observable and self.ENC_fastsim==False: 
        ipoint = 2 # EEC is 2 point correlator
        for index in range(new_corr.correlator(ipoint).rs().size()):
          part1 = new_corr.correlator(ipoint).indices1()[index]
          part2 = new_corr.correlator(ipoint).indices2()[index]
          pt1 = c_select[part1].perp()
          pt2 = c_select[part2].perp()
          pt_avg = (pt1+pt2)/2
          alpha = math.fabs(pt1-pt2)
          
          for _, (pt_lo, pt_hi) in enumerate(zip(self.trk_pt_lo,self.trk_pt_hi)):
            for _, (alpha_lo, alpha_hi) in enumerate(zip(self.trk_alpha_lo,self.trk_alpha_hi)):
              if pt_avg >= pt_lo and pt_avg < pt_hi and alpha >= alpha_lo and alpha < alpha_hi:
                if 'noweight' in observable:
                  getattr(self, hname.format(observable + str(pt_lo) + str(pt_hi) + '_' + '{:.1f}'.format(alpha_lo) + '{:.1f}'.format(alpha_hi),obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index])
                else:
                  getattr(self, hname.format(observable + str(pt_lo) + str(pt_hi) + '_' + '{:.1f}'.format(alpha_lo) + '{:.1f}'.format(alpha_hi),obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index])
                break

          # fill inclusively
          if 'noweight' in observable:
            getattr(self, hname.format(observable,obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index])
          else:
            getattr(self, hname.format(observable,obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index])

      # correlation histograms 
      if 'corr' in observable :
        # print("OBS! HERE!", observable)

        if (observable == "corr_deltap"):
          # print("what is this name1", hname.format(observable + pair_type_label,obs_label))
          # getattr(self, hname.format(observable,obs_label)).Fill(jet_pt, ??)

          # ecorrel.CorrelatorBuilder
          new_obs_corr = othercorrel.OtherCorrelatorBuilder(c_select, jet_pt, 2, 1, dphi_cut, deta_cut, "deltap")
          
        # delta pt, truth
        if (observable == "corr_deltapt"):
          new_obs_corr = othercorrel.OtherCorrelatorBuilder(c_select, jet_pt, 2, 1, dphi_cut, deta_cut, "deltapt")

        # # charge
        # if ("charge" in observable):
        #   # new_obs_corr = othercorrel.OtherCorrelatorBuilder(c_select, jet_pt, 2, 1, dphi_cut, deta_cut, "charge")
        #   self.is_same_charge(new_corr, ipoint, c_select, index)

        # unweighted RL
        if (observable == "corr_unweightedRL"): # this is WRONG, but leave for now
          new_obs_corr = ecorrel.CorrelatorBuilder(c_select, jet_pt, 2, 1, dphi_cut, deta_cut)

        # assuming the length of new_corr is the same as new_obs_corr
        fsparsejetlevelJetvalue = array.array( 'd', ( 0, 0, 0 ))
        fsparsejetlevelJetvalue[0] = jet_pt

        # print("OBS! HERE!", observable)
        for index in range(new_corr.correlator(ipoint).rs().size()):
          fsparsejetlevelJetvalue[1] = new_corr.correlator(ipoint).rs()[index]
          if ("charge" in observable):
            samecharge_boolean = self.is_same_charge(new_corr, ipoint, c_select, index)
            fsparsejetlevelJetvalue[2] = 1 if samecharge_boolean else -1
            if samecharge_boolean and observable == "corr_samecharge":
              getattr(self, hname.format(observable+str(ipoint),obs_label)).Fill(fsparsejetlevelJetvalue, new_corr.correlator(ipoint).weights()[index])
            elif not samecharge_boolean and observable == "corr_oppcharge":
              getattr(self, hname.format(observable+str(ipoint),obs_label)).Fill(fsparsejetlevelJetvalue, new_corr.correlator(ipoint).weights()[index])
          else:
            fsparsejetlevelJetvalue[2] = new_obs_corr.correlator(ipoint).rs()[index]
            getattr(self, hname.format(observable+str(ipoint),obs_label)).Fill(fsparsejetlevelJetvalue, new_corr.correlator(ipoint).weights()[index])


  #---------------------------------------------------------------
  # This function is called per observable per jet subconfigration 
  # used in fill_matched_jet_histograms
  # This function is created because we cannot use fill_observable_histograms 
  # directly because observable list loop inside that function
  #---------------------------------------------------------------
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
  # function to find D0's
  def findD0(self, particle):
    # print("user index", particle.user_index()) #idk what this number is
    # print("user info", particle.python_info()) # this is super useful! is of the form "JetInfo" from jetinfo.py
    part_info = particle.python_info()
    if abs(part_info.particle_pid) == 421:
      if self.firsttimejet:
        print("D0 in jet!", part_info.particle_pid)
      self.D0particleinfo = particle
      return True
    else:
      self.D0particleinfo = None
      return False
  #---------------------------------------------------------------

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
  parser.add_argument('-e', '--eventStartOffset', action='store',
                      type=int, metavar='eventStartOffset',
                      default=0,
                      help='If input file does not start with event 1, offset it by this amount')
  
  # Parse the arguments
  args = parser.parse_args()
  
  print('Configuring...')
  print('inputFile: \'{0}\''.format(args.inputFile))
  print('configFile: \'{0}\''.format(args.configFile))
  print('ouputDir: \'{0}\"'.format(args.outputDir))
  print('eventOffset: \'{0}\"'.format(args.eventStartOffset))

  # If invalid inputFile is given, exit
  if not os.path.exists(args.inputFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
    sys.exit(0)
  
  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = ProcessMC_ENC_HF(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir, event_start_offset=args.eventStartOffset)
  analysis.process_mc()

