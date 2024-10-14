#!/usr/bin/env python3

"""
  Analysis class to read a ROOT TTree of track information
  and do jet-finding, and save basic histograms.
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import argparse
import sys

# Data analysis and plotting
import ROOT
import yaml
import numpy as np
import array 
import math
import random

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import ecorrel
import othercorrel

# Base class
from pyjetty.alice_analysis.process.user.substructure import process_data_base

def linbins(xmin, xmax, nbins):
  lspace = np.linspace(xmin, xmax, nbins+1)
  arr = array.array('f', lspace)
  return arr

def logbins(xmin, xmax, nbins):
  lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
  arr = array.array('f', lspace)
  return arr

################################################################
class ProcessData_ENC(process_data_base.ProcessDataBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
  
    # print(sys.path)
    # print(sys.modules)

    # Initialize base class
    super(ProcessData_ENC, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)
    
    self.observable = self.observable_list[0]


  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_user_output_objects(self):

    self.fout.cd()

    for jetR in self.jetR_list:
      perpcone_R_list = []
      if self.do_jetcone:
        if self.do_only_jetcone:
          for jetcone_R in self.jetcone_R_list:
            perpcone_R_list.append(jetcone_R)
        else:
          perpcone_R_list.append(jetR)
          for jetcone_R in self.jetcone_R_list:
            if jetcone_R != jetR: # just a safeguard since jetR is already added in the list
              perpcone_R_list.append(jetcone_R)
      else:
        perpcone_R_list.append(jetR)

      for observable in self.observable_list:
        print("HISTOGRAM OBSERVABLE", observable)

        for trk_thrd in self.obs_settings[observable]:

          obs_label = self.utils.obs_label(trk_thrd, None) 

          jet_type_labels = []
          if self.is_dijet:
            xjbin = int(0.99/self.xj_interval + 1)
            for ixj in range(xjbin):
              jet_type_labels.append( '_xj{:.1f}{:.1f}'.format(self.xj_interval*ixj, self.xj_interval*(ixj+1)) )
          else:
            jet_type_labels.append( '' )
          
          # only apply jet type labels for jets (not perpcone nor jetcone). FIX ME: not sure if needed for perpcone (sig-bkg correlations)
          for jet_type_label in jet_type_labels:
            if 'ENC' in observable:
              for ipoint in range(2, 3):
                name = 'h_{}_JetPt_R{}_{}{}'.format(observable + str(ipoint), jetR, trk_thrd, jet_type_label)
                pt_bins = linbins(0,200,200)
                RL_bins = logbins(1E-4,1,50)
                h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
                h.GetXaxis().SetTitle('p_{T,ch jet}')
                h.GetYaxis().SetTitle('R_{L}')
                setattr(self, name, h)

                name = 'h_{}Pt_JetPt_R{}_{}{}'.format(observable + str(ipoint), jetR, trk_thrd, jet_type_label)
                pt_bins = linbins(0,200,200)
                ptRL_bins = logbins(1E-3,1E2,60)
                h = ROOT.TH2D(name, name, 200, pt_bins, 60, ptRL_bins)
                h.GetXaxis().SetTitle('p_{T,ch jet}')
                h.GetYaxis().SetTitle('p_{T,ch jet}R_{L}') # NB: y axis scaled by jet pt (applied jet by jet)
                setattr(self, name, h)

            if 'EEC_noweight' in observable:
              name = 'h_{}_JetPt_R{}_{}{}'.format(observable, jetR, obs_label, jet_type_label)
              pt_bins = linbins(0,200,200)
              RL_bins = logbins(1E-4,1,50)
              h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
              h.GetXaxis().SetTitle('p_{T,ch jet}')
              h.GetYaxis().SetTitle('R_{L}')
              setattr(self, name, h)

            if 'jet_pt' in observable:
              name = 'h_1D{}_JetPt_R{}_{}'.format(observable, jetR, obs_label)
              pt_bins = linbins(0,200,200)
              h = ROOT.TH1D(name, name, 200, pt_bins)
              h.GetXaxis().SetTitle('p_{T,ch jet}')
              h.GetYaxis().SetTitle('Counts')
              setattr(self, name, h)

              name = 'h_Nconst_JetPt_R{}_{}{}'.format(jetR, trk_thrd, jet_type_label)
              pt_bins = linbins(0,200,200)
              Nconst_bins = linbins(0,50,50)
              h = ROOT.TH2D(name, name, 200, pt_bins, 50, Nconst_bins)
              h.GetXaxis().SetTitle('p_{T,ch jet}')
              h.GetYaxis().SetTitle('N_{const}')
              setattr(self, name, h)

              name = 'tn_JETINFO_R{}_{}'.format(jetR, obs_label)
              tn = ROOT.TNtuple(name, name, "event_id:jet_id:jet_num_in_ev:jet_pt:total_num_const:num_const_aftercut:corr_rc:leading_q:subleading_q") #total_num_baryons:num_baryons_aftercut:total_num_mesons:num_mesons_aftercut")
              setattr(self, name, tn)
              
            
            if "corr" in observable:
              if observable == "corr_beg":
                self.tuple_obs_string = "event_id:jet_id:jet_num_in_ev:jet_pt:RL:weights"   
              elif observable == "corr_end": #only purpose of this observable is to signify the end
                name = 'tn_pairlevel_R{}_{}'.format(jetR, obs_label)
                tn = ROOT.TNtuple(name, name, self.tuple_obs_string) #don't need weights for each observable?? I'll worry about this another time
                setattr(self, name, tn)
                colon_count = self.tuple_obs_string.count(':')
                self.fsparsepartonJetvalue_tuple = array.array( 'd', np.zeros(colon_count+1)) #18 to match the number of axes
              else:
                self.create_corr_tuples(observable, jetR, obs_label)

          # fill perp cone histograms
          self.pair_type_labels = ['']

          if self.do_rho_subtraction:
            self.pair_type_labels = ['_bb','_sb','_ss']
          
          if self.do_perpcone:
            
            for perpcone_R in perpcone_R_list:
              if 'jet_pt' in observable:
                name = 'h_perpcone{}_{}_JetPt_R{}_{}'.format(perpcone_R, 'pt', jetR, trk_thrd)
                pt_bins = linbins(-200,200,400)
                h = ROOT.TH1D(name, name, 400, pt_bins)
                h.GetXaxis().SetTitle('p_{T,perp cone}')
                h.GetYaxis().SetTitle('Counts')
                setattr(self, name, h)

                name = 'h_perpcone{}_Nconst_JetPt_R{}_{}'.format(perpcone_R, jetR, trk_thrd)
                pt_bins = linbins(0,200,200)
                Nconst_bins = linbins(0,50,50)
                h = ROOT.TH2D(name, name, 200, pt_bins, 50, Nconst_bins)
                h.GetXaxis().SetTitle('p_{T,ch jet}')
                h.GetYaxis().SetTitle('N_{const}')
                setattr(self, name, h)

                name = 'h_perpcone{}_rho_local_JetPt_R{}_{}'.format(perpcone_R, jetR, trk_thrd)
                pt_bins = linbins(0,200,200)
                rho_bins = linbins(0,500,100)
                h = ROOT.TH2D(name, name, 200, pt_bins, 100, rho_bins)
                h.GetXaxis().SetTitle('p_{T,ch jet}')
                h.GetYaxis().SetTitle('local rho')
                setattr(self, name, h)

              for pair_type_label in self.pair_type_labels:

                if 'ENC' in observable:
                  for ipoint in range(2, 3):
                    name = 'h_perpcone{}_{}_JetPt_R{}_{}'.format(perpcone_R, observable + str(ipoint) + pair_type_label, jetR, trk_thrd)
                    pt_bins = linbins(0,200,200)
                    RL_bins = logbins(1E-4,1,50)
                    h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
                    h.GetXaxis().SetTitle('p_{T,ch jet}')
                    h.GetYaxis().SetTitle('R_{L}')
                    setattr(self, name, h)

                    name = 'h_perpcone{}_{}Pt_JetPt_R{}_{}'.format(perpcone_R, observable + str(ipoint) + pair_type_label, jetR, trk_thrd)
                    pt_bins = linbins(0,200,200)
                    ptRL_bins = logbins(1E-3,1E2,60)
                    h = ROOT.TH2D(name, name, 200, pt_bins, 60, ptRL_bins)
                    h.GetYaxis().SetTitle('p_{T,ch jet}R_{L}') # NB: y axis scaled by jet pt (applied jet by jet)
                    setattr(self, name, h)

                if 'EEC_noweight' in observable or 'EEC_weight2' in observable:
                  name = 'h_perpcone{}_{}_JetPt_R{}_{}'.format(perpcone_R, observable + pair_type_label, jetR, obs_label)
                  pt_bins = linbins(0,200,200)
                  RL_bins = logbins(1E-4,1,50)
                  h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
                  h.GetXaxis().SetTitle('p_{T,ch jet}')
                  h.GetYaxis().SetTitle('R_{L}')
                  setattr(self, name, h)

          if self.do_jetcone:

            for jetcone_R in self.jetcone_R_list:
              if 'ENC' in observable:
                for ipoint in range(2, 3):
                  name = 'h_jetcone{}_{}_JetPt_R{}_{}'.format(jetcone_R, observable + str(ipoint), jetR, trk_thrd)
                  pt_bins = linbins(0,200,200)
                  RL_bins = logbins(1E-4,1,50)
                  h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
                  h.GetXaxis().SetTitle('p_{T,ch jet}')
                  h.GetYaxis().SetTitle('R_{L}')
                  setattr(self, name, h)

                  name = 'h_jetcone{}_{}Pt_JetPt_R{}_{}'.format(jetcone_R, observable + str(ipoint), jetR, trk_thrd)
                  pt_bins = linbins(0,200,200)
                  ptRL_bins = logbins(1E-3,1E2,60)
                  h = ROOT.TH2D(name, name, 200, pt_bins, 60, ptRL_bins)
                  h.GetYaxis().SetTitle('p_{T,ch jet}R_{L}') # NB: y axis scaled by jet pt (applied jet by jet)
                  setattr(self, name, h)

              if 'EEC_noweight' in observable or 'EEC_weight2' in observable:
                name = 'h_jetcone{}_{}_JetPt_R{}_{}'.format(jetcone_R, observable, jetR, obs_label)
                pt_bins = linbins(0,200,200)
                RL_bins = logbins(1E-4,1,50)
                h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
                h.GetXaxis().SetTitle('p_{T,ch jet}')
                h.GetYaxis().SetTitle('R_{L}')
                setattr(self, name, h)
                    
  
  #---------------------------------------------------------------
  # Add branches to dEEC tuple
  #---------------------------------------------------------------
  def create_corr_tuples(self, observable, jetR, obs_label):

    if observable == "corr_energyweights" or observable == "corr_rc":
      return
    # if observable == "corr_oppcharge":
    #   return
    
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
    if type1 < 0 and type2 < 0:
      # print('bkg-bkg (',type1,type2,') pt1',constituents[part1].perp()
      return 0 # means bkg-bkg
    if type1 < 0 and type2 >= 0:
      # print('sig-bkg (',type1,type2,') pt1',constituents[part1].perp(),'pt2',constituents[part2].perp())
      return 1 # means sig-bkg
    if type1 >= 0 and type2 < 0:
      # print('sig-bkg (',type1,type2,') pt1',constituents[part1].perp(),'pt2',constituents[part2].perp())
      return 1 # means sig-bkg
    if type1 >= 0 and type2 >= 0:
      # print('sig-sig (',type1,type2,') pt1',constituents[part1].perp()
      return 2 # means sig-sig

  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  #---------------------------------------------------------------
  def fill_jet_histograms(self, jet, jet_groomed_lund, jetR, obs_setting, grooming_setting,
                          obs_label, jet_pt_ungroomed, suffix):

    constituents = fj.sorted_by_pt(jet.constituents())
    c_select = fj.vectorPJ()
    trk_thrd = obs_setting

    for c in constituents:
      if c.pt() < trk_thrd:
        break
      c_select.append(c) # NB: use the break statement since constituents are already sorted
    nconst_jet = len(c_select)

    if (self.debug_level == 3):
      print("CHARGE and PID HERE")
      print([c.python_info().charge for c in c_select])

    if self.ENC_pair_cut:
      dphi_cut = -9999 # means no dphi cut
      deta_cut = 0.008
    else:
      dphi_cut = -9999
      deta_cut = -9999

    
    # print("WHAT IS THIS SUFFIX?", suffix)
    hname = 'h_{}_JetPt_R{}_{}{}'
    if self.do_rho_subtraction:
      jet_pt = jet_pt_ungroomed # jet_pt_ungroomed stores subtracted jet pt for energy weight calculation and pt selection for there is a non-zero UE energy density
      if jet.area() == 0:
        return # NB: skip the zero area jets for now (also skip the perp-cone and jet-cone w.r.t. the zero area jets)
    else:
      jet_pt = jet.perp()
    # print('unsubtracted pt',jet.perp(),'subtracted',jet_pt,'# of constituents >',trk_thrd,'is',len(c_select))
    
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
    
    # save ENC and jet level information
    for observable in self.observable_list:
      # if analyze jet cones and only analyze jet cones, then only fill jet pt histograms for standard jets (just to speed things up)
      if self.do_jetcone and self.do_only_jetcone and not ('jet_pt' in observable):
        pass
      else:
        if 'ENC' in observable or 'EEC_noweight' in observable:
          for ipoint in range(2, 3):
            for index in range(new_corr.correlator(ipoint).rs().size()):

              # processing only like-sign pairs when self.ENC_pair_like is on
              if self.ENC_pair_like and (not self.is_same_charge(new_corr, ipoint, c_select, index)):
                continue

              # processing only unlike-sign pairs when self.ENC_pair_unlike is on
              if self.ENC_pair_unlike and self.is_same_charge(new_corr, ipoint, c_select, index):
                continue

              if 'ENC' in observable:
                getattr(self, hname.format(observable + str(ipoint), jetR, obs_label, suffix)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index])
                getattr(self, hname.format(observable + str(ipoint) + 'Pt', jetR, obs_label, suffix)).Fill(jet_pt, jet_pt*new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index]) # NB: fill pt*RL

              if ipoint==2 and 'EEC_noweight' in observable:
                getattr(self, hname.format(observable, jetR, obs_label, suffix)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index])

      if 'jet_pt' in observable:
        getattr(self, hname.format("1D"+observable, jetR, obs_label, suffix)).Fill(jet_pt) 
        getattr(self, hname.format('Nconst', jetR, obs_label, suffix)).Fill(jet_pt, nconst_jet) 
        
        # if observable == "corr_rc": # doesn't need this
        leadq_subleadq = self.get_leadsublead_q1q2(c_select)
        leading_q = constituents[0].python_info().charge
        if (len(constituents) >= 2):
          subleading_q = constituents[1].python_info().charge
        else:
          subleading_q = -99
        # print("indices 15, 16, 17:",self.fsparsepartonJetvalue_tuple[15], self.fsparsepartonJetvalue_tuple[16], self.fsparsepartonJetvalue_tuple[17] )   
        # below line assumes suffix = ""
        getattr(self, 'tn_JETINFO_R{}_{}'.format(jetR, obs_label)).Fill(self.event_number, self.jet_number, self.ijet, jet_pt, len(constituents), len(c_select), leadq_subleadq, leading_q, subleading_q)
    
      

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

        elif 'end' in observable:
          # print("IN CORR_END WRITING", self.fsparsepartonJetvalue_tuple)
          # print("IN CORR_END WRITING", self.fsparsepartonJetvalue_tuple[0])
          getattr(self, 'tn_pairlevel_R{}_{}'.format(jetR, obs_label)).Fill(np.array(list(self.fsparsepartonJetvalue_tuple), dtype='float32'))



          
  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  #---------------------------------------------------------------
  def fill_perp_cone_histograms(self, cone_parts, cone_R, jet, jet_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_ungroomed, suffix, rho_bge = 0):

    # calculate perp cone pt after subtraction. Notice that the perp cone already contain the particles from signal "jet". Signal and background can be identified using user_index()
    cone_px = 0
    cone_py = 0
    cone_npart = 0
    for part in cone_parts:
      if part.user_index() < 0:
        cone_px = cone_px + part.px()
        cone_py = cone_py + part.py()
        cone_npart = cone_npart + 1
    cone_pt = math.sqrt(cone_px*cone_px + cone_py*cone_py)
    # print('cone pt', cone_pt-rho_bge*jet.area(), '(', cone_pt, ')')
    cone_pt = cone_pt-rho_bge*jet.area() # ideally this should fluctuate around 0
    # print('jet pt', jet_pt_ungroomed, '(', jet.perp(), ')')

    # combine sig jet and perp cone with trk threshold cut
    trk_thrd = obs_setting
    c_select = fj.vectorPJ()
    c_select_perp = fj.vectorPJ()

    cone_parts_sorted = fj.sorted_by_pt(cone_parts)
    # print('perp cone nconst:',len(cone_parts_sorted))
    for part in cone_parts_sorted:
      if part.pt() < trk_thrd:
        break
      c_select.append(part) # NB: use the break statement since constituents are already sorted
      if part.user_index() < 0:
        c_select_perp.append(part)

    nconst_perp = len(c_select_perp)
    # print('cone R',cone_R)
    # print('total cone nconst (with thrd cut):',len(c_select))
    # print('perp cone nconst (with thrd cut):',nconst_perp)
    pt_sum_perp = 0.
    for c in c_select_perp:
      pt_sum_perp += c.pt()

    if cone_R != jetR or self.do_only_jetcone:
      rho_local_perp = pt_sum_perp / (np.pi * cone_R * cone_R)
    else:
      if self.static_perpcone == True:
        rho_local_perp = pt_sum_perp / (np.pi * jetR * jetR)
      elif jet.has_area():
        if jet.area() == 0:
        # NB: this type of hets are currently skipped in analyze_accepted_jets()
          rho_local_perp = -1
        else:
          rho_local_perp = pt_sum_perp / jet.area()
      else:
        # NB: this type of hets are currently skipped in analyze_accepted_jets()
        rho_local_perp = -1

    if self.ENC_pair_cut:
      dphi_cut = -9999 # means no dphi cut
      deta_cut = 0.008
    else:
      dphi_cut = -9999
      deta_cut = -9999

    hname = 'h_perpcone{}_{}_JetPt_R{}_{}{}'
    if self.do_rho_subtraction:
      jet_pt = jet_pt_ungroomed # jet_pt_ungroomed stores subtracted jet pt for energy weight calculation and pt selection for there is a non-zero UE energy density
      if jet.area() == 0:
        return # NB: skip the zero area jets for now (also skip the perp-cone and jet-cone w.r.t. the zero area jets)
    else:
      jet_pt = jet.perp()

    new_corr = ecorrel.CorrelatorBuilder(c_select, jet_pt, 2, 1, dphi_cut, deta_cut)
    for observable in self.observable_list:

      if 'jet_pt' in observable:
        getattr(self, hname.format(cone_R, 'pt', jetR, obs_label, suffix)).Fill(cone_pt)
        getattr(self, hname.format(cone_R, 'Nconst', jetR, obs_label, suffix)).Fill(jet_pt, nconst_perp)
        getattr(self, hname.format(cone_R, 'rho_local', jetR, obs_label, suffix)).Fill(jet_pt, rho_local_perp)

      if 'ENC' in observable or 'EEC_noweight' in observable or 'EEC_weight2' in observable:
        for ipoint in range(2, 3):
          for index in range(new_corr.correlator(ipoint).rs().size()):

            # processing only like-sign pairs when self.ENC_pair_like is on
            if self.ENC_pair_like and (not self.is_same_charge(new_corr, ipoint, c_select, index)):
              continue

            # processing only unlike-sign pairs when self.ENC_pair_unlike is on
            if self.ENC_pair_unlike and self.is_same_charge(new_corr, ipoint, c_select, index):
              continue

            # separate out sig-sig, sig-bkg, bkg-bkg correlations for EEC pairs
            pair_type_label = ''
            if self.do_rho_subtraction:
              pair_type = self.check_pair_type(new_corr, ipoint, c_select, index)
              pair_type_label = self.pair_type_labels[pair_type]

            if 'ENC' in observable:
              # print('hname is',hname.format(cone_R, observable + str(ipoint) + pair_type_label, jetR, obs_label))
              getattr(self, hname.format(cone_R, observable + str(ipoint) + pair_type_label, jetR, obs_label, suffix)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index])
              getattr(self, hname.format(cone_R, observable + str(ipoint) + pair_type_label + 'Pt', jetR, obs_label, suffix)).Fill(jet_pt, jet_pt*new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index]) # NB: fill pt*RL

            if ipoint==2 and 'EEC_noweight' in observable:
              getattr(self, hname.format(cone_R, observable + pair_type_label, jetR, obs_label, suffix)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index])

            if ipoint==2 and 'EEC_weight2' in observable:
              getattr(self, hname.format(cone_R, observable + pair_type_label, jetR, obs_label, suffix)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], pow(new_corr.correlator(ipoint).weights()[index],2))

  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  #---------------------------------------------------------------
  def fill_jet_cone_histograms(self, cone_parts, cone_R, jet, jet_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_ungroomed, suffix, rho_bge = 0):

    # combine sig jet and perp cone with trk threshold cut
    trk_thrd = obs_setting
    c_select = fj.vectorPJ()

    cone_parts_sorted = fj.sorted_by_pt(cone_parts)
    # print('perp cone nconst:',len(cone_parts_sorted))
    for part in cone_parts_sorted:
      if part.pt() < trk_thrd:
        break
      c_select.append(part) # NB: use the break statement since constituents are already sorted

    if self.ENC_pair_cut:
      dphi_cut = -9999 # means no dphi cut
      deta_cut = 0.008
    else:
      dphi_cut = -9999
      deta_cut = -9999

    hname = 'h_jetcone{}_{}_JetPt_R{}_{}{}'
    if self.do_rho_subtraction:
      jet_pt = jet_pt_ungroomed # jet_pt_ungroomed stores subtracted jet pt for energy weight calculation and pt selection for there is a non-zero UE energy density
      if jet.area() == 0:
        return # NB: skip the zero area jets for now (also skip the perp-cone and jet-cone w.r.t. the zero area jets)
    else:
      jet_pt = jet.perp()

    new_corr = ecorrel.CorrelatorBuilder(c_select, jet_pt, 2, 1, dphi_cut, deta_cut)
    for observable in self.observable_list:

      if 'ENC' in observable or 'EEC_noweight' in observable or 'EEC_weight2' in observable:
        for ipoint in range(2, 3):
          for index in range(new_corr.correlator(ipoint).rs().size()):

            # processing only like-sign pairs when self.ENC_pair_like is on
            if self.ENC_pair_like and (not self.is_same_charge(new_corr, ipoint, c_select, index)):
              continue

            # processing only unlike-sign pairs when self.ENC_pair_unlike is on
            if self.ENC_pair_unlike and self.is_same_charge(new_corr, ipoint, c_select, index):
              continue

            if 'ENC' in observable:
              # print('hname is',hname.format(cone_R, observable + str(ipoint), jetR, obs_label))
              getattr(self, hname.format(cone_R, observable + str(ipoint), jetR, obs_label, suffix)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index])
              getattr(self, hname.format(cone_R, observable + str(ipoint) + 'Pt', jetR, obs_label, suffix)).Fill(jet_pt, jet_pt*new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index]) # NB: fill pt*RL

            if ipoint==2 and 'EEC_noweight' in observable:
              getattr(self, hname.format(cone_R, observable, jetR, obs_label, suffix)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index])

            if ipoint==2 and 'EEC_weight2' in observable:
              getattr(self, hname.format(cone_R, observable, jetR, obs_label, suffix)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], pow(new_corr.correlator(ipoint).weights()[index],2))

##################################################################
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description='Process data')
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
  
  # Parse the arguments
  args = parser.parse_args()
  
  print('Configuring...')
  print('inputFile: \'{0}\''.format(args.inputFile))
  print('configFile: \'{0}\''.format(args.configFile))
  print('ouputDir: \'{0}\"'.format(args.outputDir))
  print('----------------------------------------------------------------')
  
  # If invalid inputFile is given, exit
  if not os.path.exists(args.inputFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
    sys.exit(0)
  
  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = ProcessData_ENC(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_data()