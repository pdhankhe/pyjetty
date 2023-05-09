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
import pyjetty.alihfjets.dev.hfjet.process_io_data_hf as hfdio
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
    self.initializeHistograms()
    self.file_simulation = ROOT.TFile(self.datafile, 'READ')
    self.fill_truth_before_matching()
    self.fill_matching_histograms()
    self.fill_response_histograms()
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
    self.ptjet_cut = config['ptjet_cut']
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
  
    hjetptvsdmesonpt_truth_prompt = self.file_simulation.Get("truth_jet_prompt")
    # Rebin to correct binning
    hjetptvsdmesonpt_truth_prompt_rebinned = self.histutils.rebin_th2(
                  hjetptvsdmesonpt_truth_prompt, hjetptvsdmesonpt_truth_prompt.GetName()+"_rebinned", jet_pt_bin_array,
                      len(self.fptbinsJetFinalA)-1, cand_pt_bin_array,
              len(self.fptbinsDA)-1)
    hjetptvsdmesonpt_truth_prompt_rebinned_proj =hjetptvsdmesonpt_truth_prompt_rebinned.ProjectionX()
    hjetptvsdmesonpt_truth_prompt_rebinned_proj.SetName("hJetPt_Truth_R0.4")
   
    hjetptvsdmesonpt_truth_nonprompt = self.file_simulation.Get("truth_jet_fd")
    hjetptvsdmesonpt_det_prompt = self.file_simulation.Get("particle_jet_prompt")
    hjetptvsdmesonpt_det_nonprompt = self.file_simulation.Get("particle_jet_fd")

  #---------------------------------------------------------------
  # Loop through jets and fill matching histos
  #---------------------------------------------------------------
  def fill_matching_histograms(self):
    obs_bin_array = array('d', self.fobs_bins)
    jet_pt_bin_array = array('d', self.fptbinsJetFinalA)
    cand_pt_bin_array = array('d', self.fptbinsDA)
    
    sparseprompt = self.file_simulation.Get("THnSparse_prompt_signal_alpha_1")
    sparsenonprompt = self.file_simulation.Get("THnSparse_nonprompt_signal_alpha_1")
    
    sparse_JES = sparseprompt.Clone("sparse_JES")
    sparse_JES.Add(sparsenonprompt)
    #sparse_JES.GetAxis(0).SetRangeUser(self.jetrange[0], self.jetrange[1])
    #sparse_JES.GetAxis(1).SetRangeUser(self.jetrange[0], self.jetrange[1])

    
    jet_pt_prompt_info = sparse_JES.Projection(0,1)
    jet_pt_prompt_info.SetName("Jet_pt_prompt_info")
    
    #jet_pt_info_rebinned = self.histutils.rebin_th2(
     #                     jet_pt_prompt_info, jet_pt_prompt_info.GetName()+"_rebinned", jet_pt_bin_array,
      #                len(self.fptbinsJetFinalA)-1, jet_pt_bin_array, len(self.fptbinsJetFinalA)-1)
    
    jet_pt_truth = jet_pt_prompt_info.ProjectionY()
    jet_pt_det = jet_pt_prompt_info.ProjectionX()
    
    JES = jet_pt_det.Clone("JES_prompt")
    JES.Add(jet_pt_truth, -1)
    JES.Divide(JES, jet_pt_truth, 1 , 1, "B")
    hJES = ROOT.TH1F("hJES", "hJES" ,120, -0.6, 0.6)
    
    for ibin in range(JES.GetNbinsX()):
        JES_content = JES.GetBinContent(ibin+1)
        print(JES_content)
        hJES.Fill(JES_content)
        
    canvas_jes = ROOT.TCanvas()
    hJES.Draw("same,ep")
    canvas_jes.SaveAs("JES_prompt.pdf")
            
    


    #getattr(self, 'hResponse_JetPt_R%s' %
    #        str(jetR).replace('.', '')).Fill(jet_pt_det_ungroomed, jet_pt_truth_ungroomed)

  #---------------------------------------------------------------
  # Fill response histograms
  #---------------------------------------------------------------
  def fill_response_histograms(self):
    self.efficiency()
    obs_bin_array = array('d', self.fobs_bins)
    jet_pt_bin_array = array('d', self.fptbinsJetFinalA)
    cand_pt_bin_array = array('d', self.fptbinsDA)
    
    sparseprompt = self.file_simulation.Get("THnSparse_prompt_signal_alpha_1")
    sparsenonprompt = self.file_simulation.Get("THnSparse_nonprompt_signal_alpha_1")
    ndim = 4
    dim = [array('i',(1, 0, 6, 5))]
    
    for beta in self.beta_list:
        print("beta value is ", beta)
        #label = "R%s_%s" % (str(0.4).replace('.', ''), str(beta).replace('.', ''))
        label = "R%s_%s" % (str(0.4), str(beta))


        for ipt in range(self.fptbinsDN):
    
            ptbin = (self.fptbinsDA[ipt]+self.fptbinsDA[ipt+1]) / 2.
            sparse_response_prompt = sparseprompt.Clone("sparse_prompt_{}".format(ipt))
            sparse_response_prompt.GetAxis(0).SetRangeUser(2, 60)
            sparse_response_prompt.GetAxis(1).SetRangeUser(2, 60)
            sparse_response_prompt.GetAxis(3).SetRangeUser(self.fptbinsDA[ipt], self.fptbinsDA[ipt+1]) #dmesonpt
            prompt_response = sparse_response_prompt.Projection(ndim, dim[0], "E")
            prompt_response.Scale(1. /(self.efficiency_prompt.GetBinContent(self.efficiency_prompt.GetXaxis().FindBin(ptbin))))
            if (ipt == 0):
                prompt_response_matrix =  prompt_response.Clone("hResponse_JetPt_ang_%sScaled" % label)
            else:
                prompt_response_matrix.Add(prompt_response)
        
        title = ['#it{p}_{T,det}^{ch jet}', '#it{p}_{T,truth}^{ch jet}',
                 '#it{#lambda}_{#it{#alpha=1},det}', '#it{#lambda}_{#it{#alpha=1},truth}']
    
        for i in range(0, 4):
            prompt_response_matrix.GetAxis(i).SetTitle(title[i])
        name = "hResponse_JetPt_ang_%sScaled" % label
        prompt_response_matrix.Sumw2()
        setattr(self, name, prompt_response_matrix)

        for ipt in range(self.fptbinsDN):
            ptbin = (self.fptbinsDA[ipt]+self.fptbinsDA[ipt+1]) / 2.
            sparse_response_nonprompt = sparsenonprompt.Clone("sparse_nonprompt_{}".format(ipt))
            sparse_response_nonprompt.GetAxis(0).SetRangeUser(2, 60)
            sparse_response_nonprompt.GetAxis(1).SetRangeUser(2, 60)
            sparse_response_nonprompt.GetAxis(3).SetRangeUser(self.fptbinsDA[ipt], self.fptbinsDA[ipt + 1]) #dmesonpt
            nonprompt_response = sparse_response_nonprompt.Projection(ndim, dim[0], "E")
            nonprompt_response.Scale(1. /(self.efficiency_prompt.GetBinContent(self.efficiency_nonprompt.GetXaxis().FindBin(ptbin))))
        
            if (ipt == 0):
                nonprompt_response_matrix =  nonprompt_response.Clone("hResponse_nonprompt_JetPt_ang_%sScaled" % label)
            else:
                nonprompt_response_matrix.Add(nonprompt_response)
    
        title = ['#it{p}_{T,det}^{ch jet}', '#it{p}_{T,truth}^{ch jet}',
                 '#it{#lambda}_{#it{#alpha=1},det}', '#it{#lambda}_{#it{#alpha=1},truth}']
                 
        for i in range(0, 4):
            nonprompt_response_matrix.GetAxis(i).SetTitle(title[i])
    
        nonprompt_response_matrix.Sumw2()
        name = "hResponse_nonprompt_JetPt_ang_%sScaled" % label
        setattr(self, name, nonprompt_response_matrix)

      
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

