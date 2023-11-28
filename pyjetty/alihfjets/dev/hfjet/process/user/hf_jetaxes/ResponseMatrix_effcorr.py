#! /usr/bin/env python

"""
HF post processing analysis code
Preeti Dhankher, 2020 (pdhankher@berkeley.edu)
"""
import pandas as pd
import sys
import os
import argparse
from array import *
import numpy as np
import math   # for exp()
import ROOT
import uproot
import time
import pyjetty.alihfjets.dev.hfjet.process.base.process_io_data_hf as hfdio
#ROOT.gSystem.Load("$HEPPY_DIR/external/roounfold/roounfold-current/lib/libRooUnfold.so")
import yaml
# For log y plots which ROOT just decides not to work for
import matplotlib
#matplotlib.rcParams['text.usetex'] = True   # LaTeX labels
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams["yaxis.labellocation"] = 'top'
plt.rcParams["xaxis.labellocation"] = 'right'

from pyjetty.alice_analysis.analysis.user.substructure import run_analysis

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io, process_utils, jet_info, process_base

# Load pyjetty ROOT utils
ROOT.gSystem.Load('libpyjetty_rutil')
#ROOT.gSystem.Load('libpyjetty_rutil.dylib')

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(False)


#----------------------------------------------------------------------
# Extrapolate y-values for values in xlist_new given points (x,y) in xlist and ylist
# Use power=1 for linear, or power=2 for quadratic extrapolation
#----------------------------------------------------------------------
    
def list_interpolate(xlist, ylist, xlist_new, power=1, require_positive=False):

  if len(xlist) < (power + 1):
    raise ValueError("list_interpolate() requires at least %i points!" % (power + 1))

  ylist_new = []
  ix = 0
  for xval in xlist_new:

    while (ix + power) < len(xlist) and xlist[ix+power] <= xval:
      ix += 1

    x1 = xlist[ix]; y1 = ylist[ix]

    # Check if data point is identical
    if xval == x1:
      if require_positive and y1 < 0:
        ylist_new.append(0)
        continue
      ylist_new.append(y1)
      continue

    # Set value to 0 if out-of-range for extrapolation
    if x1 > xval or (ix + power) >= len(xlist):
      ylist_new.append(0)
      continue

    x2 = xlist[ix+1]; y2 = ylist[ix+1]

    yval = None
    if power == 1:  # linear
      yval = linear_extrapolate(x1, y1, x2, y2, xval)
    elif power == 2:  # quadratic
      x3 = xlist[ix+2]; y3 = ylist[ix+2]
      yval = quadratic_extrapolate(x1, y1, x2, y2, x3, y3, xval)
    else:
      raise ValueError("Unrecognized power", power, "/ please use either 1 or 2")

    # Require positive values
    if require_positive and yval < 0:
      ylist_new.append(0)
      continue

    ylist_new.append(yval)

  return ylist_new


#---------------------------------------------------------------
# Given two data points, find linear fit and y-value for x
#---------------------------------------------------------------
def linear_extrapolate(x1, y1, x2, y2, x):

  return (y2 - y1) / (x2 - x1) * x + (y1 - (y2 - y1) / (x2 - x1) * x1)


#---------------------------------------------------------------
# Given three data points, find quadratic fit and y-value for x
#---------------------------------------------------------------
def quadratic_extrapolate(x1, y1, x2, y2, x3, y3, x):

  a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2))
  b = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3))
  c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3)

  return (a * x * x + b * x + c) / ((x1 - x2) * (x1 - x3) * (x2 - x3))


#---------------------------------------------------------------
# Set LHS of distributions to 0 if crosses to 0 at some point (prevents multiple peaks)
#---------------------------------------------------------------
def set_zero_range(yvals):

  found_nonzero_val = False

  # Step through list backwards
  for i in range(len(yvals)-1, -1, -1):
    if yvals[i] <= 0:
      if found_nonzero_val:
        for j in range(0, i+1):
          yvals[j] = 0
        break
      yvals[i] = 0
      continue
    else:
      found_nonzero_val = True
      continue

  return yvals


# Where there are single values pos/neg between two neg/pos, interpolate point
def fix_fluctuations(yvals):

  for i in range(1, len(yvals) - 1):
    if yvals[i] > 0:
      if yvals[i+1] < 0 and yvals[i-1] < 0:
        yvals[i] = (yvals[i+1] + yvals[i-1]) / 2
    else:  # yvals[i] <= 0
      if yvals[i+1] > 0 and yvals[i-1] > 0:
        yvals[i] = (yvals[i+1] + yvals[i-1]) / 2

  return yvals

def process_canvas(Canvas):
    ROOT.gStyle.SetOptStat(0)
    Canvas.SetHighLightColor(2)
    Canvas.SetFillColor(0)
    Canvas.SetBorderMode(0)
    Canvas.SetBorderSize(2)
    Canvas.SetTickx(1)
    Canvas.SetTicky(1)
    Canvas.SetFrameBorderMode(0)
    Canvas.SetFrameLineWidth(2)
    Canvas.SetFrameBorderMode(1)

################################################################
#######################  RUN ANALYSIS  #########################
################################################################
class RunAnalysisAng(run_analysis.RunAnalysis):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, config_file='', **kwargs):
    super(RunAnalysisAng, self).__init__(config_file, **kwargs)
    
    start_time = time.time()
    # Initialize yaml config
    self.initialize_user_config()
    #after initializing, loop over each axis type
    for i in range(len(self.hsparse_jet)):
        isparse = self.hsparse_jet[i]
        print("Beginning analysis for ",isparse)
        self.initializeHistograms()
        self.file_simulation = ROOT.TFile(self.datafile, 'READ')
        self.fill_truth_before_matching()
        self.fill_response_histograms(isparse)
        # Plot histograms
        print('Save histograms...')
        process_base.ProcessBase.save_output_objects(self)
    
    #print(self)
  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_user_config(self):
    
    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)
      
    self.figure_approval_status = config['figure_approval_status']
    self.plot_overlay_list = \
      self.obs_config_dict['common_settings']['plot_overlay_list']
    
    self.jet_matching_distance = config['jet_matching_distance']
    
    if 'constituent_subtractor' in config:
        self.is_pp = False
    else:
        self.is_pp = True
    print('is_pp: {}'.format(self.is_pp))

    self.histutils = ROOT.RUtil.HistUtils()
    self.jetrange = config['jet_range']
    self.fptbinsJetFinalA = config['jetptbinsFinalA']
    self.fptbinsJetMeasN = config['jetptbinsMeasN']
    self.fptbinsJetMeasA = config['jetptbinsMeasA']
    self.jetR_list = config['jetR']
    self.beta_list = config['betas']
    self.fobs_bins = config['obs_bins']
    self.IsUseRefl = config['fUseRefl']
    self.fptbinsDN = config['ptbinsDN']
    self.fptbinsDA = config['ptbinsDA']
    self.RebinMass = config['fRebinMass']
    self.massrange = config['fMassRange']
    self.pt_dmeson_minpt_cut = config['pt_dmeson_cut']
    #self.ptjet_cut = config['ptjet_cut']
    self.Dmass = config['fDmass']
    self.Dsigma = config['fDsigma']
    self.Bkgsigma = config ['fsigmaBkg']
    self.Signalsigma = config ['fsigmaSignal']
    self.bkgtype = config['fbkgtype']
    self.IsUseRefl = config['fUseRefl']
    self.xsection_inel = config['xsection_inel']
    self.p_nevents = config['p_nevents']
    self.branching_ratio = config['branching_ratio']

    self.datafile = config['main_data_simulation']
    print(self.datafile)

    self.efffile = config['main_efffile']
    print(self.efffile)
 
    self.hsparse_jet = config['hsparse_jet']
    print(self.hsparse_jet)

    self.isRefSys = False
   
    

    # Whether or not to use the previous preliminary result in final plots
    self.use_prev_prelim = config['use_prev_prelim']

    self.histutils = ROOT.RUtil.HistUtils()

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  
  def initializeHistograms(self):
    obs_bin_array = array('d', self.fobs_bins)
    jet_pt_bin_array = array('d', self.fptbinsJetFinalA)
    cand_pt_bin_array = array('d', self.fptbinsDA)
    
        
  #---------------------------------------------------------------
  # efficiency function
  #---------------------------------------------------------------

  def efficiency(self):
    file_eff = ROOT.TFile(self.efffile, 'READ')
    print(self.efffile)
    self.efficiency_prompt = file_eff.Get("PromptEfficiency")
    self.efficiency_prompt.SetDirectory(0)
    self.efficiency_nonprompt = file_eff.Get("NonPromptEfficiency")
    self.efficiency_nonprompt.SetDirectory(0)
    return True
    
  #---------------------------------------------------------------
  # Fill truth jet histograms
  #---------------------------------------------------------------
  def fill_truth_before_matching(self):
   
    obs_bin_array = array('d', self.fobs_bins)
    jet_pt_bin_array = array('d', self.fptbinsJetFinalA)
    cand_pt_bin_array = array('d', self.fptbinsDA)
  
    hjetptvsdmesonpt_truth_prompt = self.file_simulation.Get("hjetvsDpt_prompt_gen")
    # Rebin to correct binning
    hjetptvsdmesonpt_truth_prompt_rebinned = self.histutils.rebin_th2(
                  hjetptvsdmesonpt_truth_prompt, hjetptvsdmesonpt_truth_prompt.GetName()+"_rebinned", jet_pt_bin_array,
                      len(self.fptbinsJetFinalA)-1, cand_pt_bin_array,
              len(self.fptbinsDA)-1)
    hjetptvsdmesonpt_truth_prompt_rebinned_proj =hjetptvsdmesonpt_truth_prompt_rebinned.ProjectionX()
    hjetptvsdmesonpt_truth_prompt_rebinned_proj.SetName("hJetPt_Truth_R0.4")
   
    hjetptvsdmesonpt_truth_nonprompt = self.file_simulation.Get("hjetvsDpt_fd_gen")
    #hjetptvsdmesonpt_det_prompt = self.file_simulation.Get("particle_jet_prompt")
    #hjetptvsdmesonpt_det_nonprompt = self.file_simulation.Get("particle_jet_fd")

 
  #---------------------------------------------------------------
  # Fill response histograms
  #---------------------------------------------------------------
  def fill_response_histograms(self, isparse):
    self.efficiency()
    obs_bin_array = array('d', self.fobs_bins)
    jet_pt_bin_array = array('d', self.fptbinsJetFinalA)
    cand_pt_bin_array = array('d', self.fptbinsDA)
    
    sparseprompt = self.file_simulation.Get('hsparse_R0.4_prompt_'+isparse)
    
    sparsenonprompt = self.file_simulation.Get('hsparse_R0.4_nonprompt_'+isparse) #"hsparse_R0.4_nonprompt_STD_D0" 
    print("prompt = ",sparseprompt,", nonprompt = ",sparsenonprompt)

    ndim = 4
    dim = [array('i',(0, 1, 4, 5))] #1,0,6,5 or 
    
    sparseprompt_proj_jetpt_truth = sparseprompt.Projection(1)
    sparseprompt_proj_jetpt_det = sparseprompt.Projection(0)
    sparseprompt_proj_obs_truth = sparseprompt.Projection(5)
    sparseprompt_proj_obs_det = sparseprompt.Projection(4)
    
    #get bin array for all bins
    pt_det_bin = []
    pt_truth_bin = []
    obs_det_bin = []
    obs_truth_bin = []
    
    for i in range(1,sparseprompt_proj_jetpt_truth.GetNbinsX()):
        pt_truth_bin.append(sparseprompt_proj_jetpt_truth.GetBinLowEdge(i))
        
    for i in range(1,sparseprompt_proj_jetpt_det.GetNbinsX()):
        pt_det_bin.append(sparseprompt_proj_jetpt_det.GetBinLowEdge(i))
        
    for i in range(1,sparseprompt_proj_obs_truth.GetNbinsX()):
        obs_truth_bin.append(sparseprompt_proj_obs_truth.GetBinLowEdge(i))
        
    for i in range(1,sparseprompt_proj_obs_det.GetNbinsX()):
        obs_det_bin.append(sparseprompt_proj_obs_det.GetBinLowEdge(i))

    
    pt_det_bin_array = array('d', pt_det_bin)
    pt_truth_bin_array = array('d', pt_truth_bin)
    obs_truth_bin_array = array('d', obs_truth_bin)
    obs_det_bin_array = array('d', obs_det_bin)
    
    for beta in self.beta_list:
        print("beta value is ", beta)
        label = "R%s_%s" % (str(0.4), str(beta))

        prompt_response_matrix_list = []
        for ibinpt in range(len(self.fptbinsJetFinalA)-1):
            first_fit = 0
            for idptbin in range(self.fptbinsDN):
                if self.fptbinsDA[idptbin] < self.pt_dmeson_minpt_cut[ibinpt]:
                    print("analysis exit because of low dmeson")
                    continue
                if self.fptbinsDA[idptbin+1] > self.fptbinsJetFinalA[ibinpt+1]:
                    print("analysis exit because of jet threshold")
                    break
                ptbin = (self.fptbinsDA[idptbin]+self.fptbinsDA[idptbin+1]) / 2.
                sparse_response_prompt = sparseprompt.Clone("sparse_prompt_{}".format(idptbin))
                sparse_response_prompt.GetAxis(0).SetRangeUser(self.fptbinsJetFinalA[ibinpt+1], self.fptbinsJetFinalA[ibinpt+1]) #dmesonpt
                sparse_response_prompt.GetAxis(1).SetRangeUser(self.fptbinsJetFinalA[ibinpt+1], self.fptbinsJetFinalA[ibinpt+1]) #dmesonpt
                sparse_response_prompt.GetAxis(3).SetRangeUser(self.fptbinsDA[idptbin], self.fptbinsDA[idptbin+1]) #dmesonpt
                #sparse_response_prompt.GetAxis(2).SetRangeUser(self.fptbinsDA[idptbin], self.fptbinsDA[idptbin+1])
                print(sparse_response_prompt)
                
                prompt_response = sparse_response_prompt.Projection(ndim, dim[0], "E")
                #prompt_response.Scale(1. /(self.efficiency_prompt.GetBinContent(self.efficiency_prompt.GetXaxis().FindBin(ptbin))))
                name_ch = "sparse_prompt_{}".format(idptbin)
                print(prompt_response)
                prompt_response_rebin=self.histutils.rebin_thn(prompt_response, '%s_Rebinned' % (name_ch), ndim, len(pt_det_bin_array)-1, pt_det_bin_array, len(obs_det_bin_array)-1, obs_det_bin_array, len(pt_truth_bin_array)-1, pt_truth_bin_array, len(obs_truth_bin_array)-1, obs_truth_bin_array,label)

                # efficiency correction after rebinning otherwise the error set won't be correct
                prompt_response_rebin.Scale(1. /(self.efficiency_prompt.GetBinContent(self.efficiency_prompt.GetXaxis().FindBin(ptbin))))
                
                if (first_fit == 0):
                    prompt_response_matrix =  prompt_response_rebin.Clone("hResponse_JetPt_jet_axis_%sScaled" % label)
                    first_fit = 1
                else:
                    prompt_response_matrix.Add(prompt_response_rebin)
                    
            prompt_response_matrix_list.append(prompt_response_matrix)
            
            if ibinpt==0:
                prompt_response_perbin= prompt_response_matrix_list[ibinpt].Clone("prompt_response_perbin_{}".format(self.fptbinsJetFinalA[ibinpt]))
            else:
                prompt_response_perbin.Add(prompt_response_matrix_list[ibinpt])
                
        
        title = ['#it{p}_{T,det}^{ch jet}', '#it{p}_{T,truth}^{ch jet}',
                 '#it{#DeltaR}_{det}', '#it{#DeltaR}_{truth}']
    
        for i in range(0, 4):
            prompt_response_perbin.GetAxis(i).SetTitle(title[i])
        name = "hResponse_JetPt_jet_axis_%sScaled" % (label)
        prompt_response_perbin.SetName(name)
        prompt_response_perbin.Sumw2()
        setattr(self, name, prompt_response_perbin)
        
    
        nonprompt_response_matrix_list = []
        for ibinpt in range(len(self.fptbinsJetFinalA)-1):
            first_fit = 0
            for idptbin in range(self.fptbinsDN):
                if self.fptbinsDA[idptbin] < self.pt_dmeson_minpt_cut[ibinpt]:
                    print("analysis exit because of low dmeson")
                    continue
                if self.fptbinsDA[idptbin+1] > self.fptbinsJetFinalA[ibinpt+1]:
                    print("analysis exit because of jet threshold")
                    break
                ptbin = (self.fptbinsDA[idptbin]+self.fptbinsDA[idptbin+1]) / 2.
                sparse_response_nonprompt = sparsenonprompt.Clone("sparse_nonprompt_{}".format(idptbin))
                sparse_response_nonprompt.GetAxis(0).SetRangeUser(self.fptbinsJetFinalA[ibinpt+1], self.fptbinsJetFinalA[ibinpt+1]) #dmesonpt
                sparse_response_nonprompt.GetAxis(1).SetRangeUser(self.fptbinsJetFinalA[ibinpt+1], self.fptbinsJetFinalA[ibinpt+1]) #dmesonpt
                sparse_response_nonprompt.GetAxis(3).SetRangeUser(self.fptbinsDA[idptbin], self.fptbinsDA[idptbin+1]) #dmesonpt
                nonprompt_response = sparse_response_nonprompt.Projection(ndim, dim[0], "E")
                
                name_ch = "sparse_nonprompt_{}".format(idptbin)
                nonprompt_response_rebin=self.histutils.rebin_thn(nonprompt_response, '%s_Rebinned' % (name_ch), ndim, len(pt_det_bin_array)-1, pt_det_bin_array, len(obs_det_bin_array)-1, obs_det_bin_array, len(pt_truth_bin_array)-1, pt_truth_bin_array, len(obs_truth_bin_array)-1, obs_truth_bin_array,label)
                # efficiency correction after rebinning otherwise the error set won't be correct
                nonprompt_response_rebin.Scale(1. /(self.efficiency_prompt.GetBinContent(self.efficiency_prompt.GetXaxis().FindBin(ptbin))))
                if (first_fit == 0):
                    nonprompt_response_matrix =  nonprompt_response_rebin.Clone("hResponse_JetPt_jet_axis_%sScaledfd" % label)
                    first_fit = 1
                else:
                    nonprompt_response_matrix.Add(nonprompt_response_rebin)
                    
            nonprompt_response_matrix_list.append(nonprompt_response_matrix)
            
            if ibinpt==0:
                nonprompt_response_perbin= nonprompt_response_matrix_list[ibinpt].Clone("nonprompt_response_perbin_{}".format(self.fptbinsJetFinalA[ibinpt]))
            else:
                nonprompt_response_perbin.Add(nonprompt_response_matrix_list[ibinpt])
                
        
        title = ['#it{p}_{T,det}^{ch jet}', '#it{p}_{T,truth}^{ch jet}',
                 '#it{#DeltaR}_{det}', '#it{#DeltaR}_{truth}']
    
        for i in range(0, 4):
            nonprompt_response_perbin.GetAxis(i).SetTitle(title[i])
        name = "hResponse_JetPt_jet_axis_%sScaledfd" % label
        nonprompt_response_perbin.SetName(name)
        nonprompt_response_perbin.Sumw2()
        setattr(self, name, nonprompt_response_perbin)

      
  #---------------------------------------------------------------
  # This function is called once for each subconfiguration
  #---------------------------------------------------------------
  def plot_single_result(self, jetR, obs_label, obs_setting, grooming_setting):
    #print('Plotting each individual result...')

    # Plot final result for each 1D substructure distribution (with PYTHIA)
    self.plot_final_result(jetR, obs_label, obs_setting, grooming_setting)



#----------------------------------------------------------------------
if __name__ == '__main__':

  # Define arguments
  parser = argparse.ArgumentParser(description='Jet substructure analysis')
  parser.add_argument('-c', '--configFile', action='store',
                      type=str, metavar='configFile',
                      default='analysis_config.yaml',
                      help='Path of config file for analysis')

  # Parse the arguments
  args = parser.parse_args()
  
  print('Configuring...')
  print('configFile: \'{0}\''.format(args.configFile))
  
  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = RunAnalysisAng(config_file = args.configFile)


