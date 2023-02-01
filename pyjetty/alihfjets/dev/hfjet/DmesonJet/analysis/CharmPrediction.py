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

# Load pyjetty ROOT utils
ROOT.gSystem.Load('libpyjetty_rutil')
#ROOT.gSystem.Load('libpyjetty_rutil.dylib')

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(False)


################################################################
# Helper functions
################################################################


def folding(h_input, response_matrix, h_output):
    h_folded = h_output.Clone("h_folded")
    for a in range(h_output.GetNbinsX()):
        for b in range(h_output.GetNbinsY()):
            val = 0.0
            val_err = 0.0
            for k in range(h_input.GetNbinsX()):
                for l in range(h_input.GetNbinsY()):
                    index_x_out = a+ h_output.GetNbinsX()*b
                    index_x_in = k + h_input.GetNbinsX()*l
                    val = val + h_input.GetBinContent(k+1, l+1) * \
                        response_matrix(index_x_out, index_x_in)
                    val_err = val_err + h_input.GetBinError(k+1, l+1) * \
                        h_input.GetBinError(k+1, l+1)* \
                        response_matrix(index_x_out, index_x_in) * \
                        response_matrix(index_x_out, index_x_in)
            h_folded.SetBinContent(a+1, b+1, val)
            h_folded.SetBinError(a+1, b+1, math.sqrt(val_err))
    return h_folded
    
def seldf_singlevar(dataframe, var, minval, maxval):
    """
    Make projection on variable using [X,Y), e.g. pT or multiplicity
    """
    dataframe = dataframe.loc[(dataframe[var] >= minval) & (dataframe[var] < maxval)]
    return dataframe

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



################################################################
#######################  RUN ANALYSIS  #########################
################################################################
class RunAnalysisAng(run_analysis.RunAnalysis):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, config_file='', **kwargs):
    super(RunAnalysisAng, self).__init__(config_file, **kwargs)
    
    # Initialize yaml config
    self.initialize_user_config()
    self.charmprediction()
    
    
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
    self.Dmass = config['fDmass']
    self.Dsigma = config['fDsigma']
    self.Bkgsigma = config ['fsigmaBkg']
    self.Signalsigma = config ['fsigmaSignal']
    self.bkgtype = config['fbkgtype']
    self.IsUseRefl = config['fUseRefl']
    self.xsection_inel = config['xsection_inel']
    self.p_nevents = config['p_nevents']
    self.branching_ratio = config['branching_ratio']

    self.datafile = config['main_data']
    print(self.datafile)
    self.simulationfile = config['main_charm_simulation']
    print(self.simulationfile)
    self.efffile = config['main_efffile']
    print(self.efffile)
    self.reflection_template = config['reflection_template']
    self.isRefSys = False
    print(self.reflection_template[0])
    

    # Whether or not to use the previous preliminary result in final plots
    self.use_prev_prelim = config['use_prev_prelim']

    self.histutils = ROOT.RUtil.HistUtils()

  

  def get_simulated_yields(self, file_path: str, dim: int, prompt: bool, xsec = False):
    """Create a histogram from a simulation tree.
    file_path - input file path
    dim - dimension of the output histogram: 2, 3
    prompt - prompt or non-prompt: True, False"""
    print("Starting the histogram extraction from an MC tree\nInput file: {}".format(file_path))
    
    if dim not in (2, 3):
        print("Error: {} is not a supported dimension.".format(dim))
    
    file_sim = ROOT.TFile.Open(file_path)
    if not file_sim:
        print("Error: file not found")
    pr_xsec = file_sim.Get("fHistXsection")
    if not pr_xsec:
        print("Error: fHistXsection histogram not found ")
        
    self.scale_factor_pr_xsec = pr_xsec.GetBinContent(1) / pr_xsec.GetEntries()
    print("scale factor ",self.scale_factor_pr_xsec )
    file_sim.Close()
    
    # Load the tree.
    
    
    tree_name = "tree_D0"
        
    tree_sim = uproot.open(file_path)[tree_name]
    if not tree_sim:
        print("Error: Tree not found")
        
    print("Converting")
    list_branches = ["pt_cand", "eta_cand", "phi_cand", "y_cand", "pdg_parton", "pt_jet", \
            "eta_jet", "phi_jet", "delta_r_jet", "LB1K1_jet", "LB15K1_jet", "LB2K1_jet", "LB3K1_jet", "delta_r_jet_std", "delta_r_jet_wta", "delta_r_jet_sd", "delta_r_jet_std_sd", "delta_r_jet_std_wta", "delta_r_jet_wta_sd" ]
            
    if self.obs_axis == 3:
        self.var = "LB1K1_jet"
    if self.obs_axis == 4:
        self.var = "LB15K1_jet"
    if self.obs_axis == 5:
        self.var = "LB2K1_jet"
    if self.obs_axis == 6:
        self.var = "LB3K1_jet"
        
        
    if self.obs_axis == 7:
        self.var = "delta_r_jet_std"
    if self.obs_axis == 8:
        self.var = "delta_r_jet_wta"
    if self.obs_axis == 9:
        self.var = "delta_r_jet_sd"
        
    if self.obs_axis == 10:
        self.var = "delta_r_jet_std_sd"
    if self.obs_axis == 11:
        self.var = "delta_r_jet_std_wta"
    if self.obs_axis == 12:
        self.var = "delta_r_jet_wta_sd"
    
    df_sim = uproot.concatenate(tree_sim, branches=list_branches, library="pd")
    #df_sim = tree_sim.pandas.df(branches=list_branches)
    hfaio = hfdio.HFAnalysis()
    sel_cand_array = hfaio.apply_cut_fiducial_acceptance(df_sim)
    pdg_parton_good = 4 if prompt else 5
    df_sim = df_sim[df_sim["pdg_parton"] == pdg_parton_good]

    # cut on jet pt
    df_sim = seldf_singlevar(df_sim, "pt_jet", 5, 50)
    # cut on Dmeson pt
    df_sim = seldf_singlevar(df_sim, "pt_cand", 2, 36)
    # cut on shape
    df_sim = seldf_singlevar(df_sim, self.var, 0, 0.8)
    # acceptance cut
    df_sim = seldf_singlevar(df_sim, "eta_jet", -0.5, 0.5)
    
    print("Entries after filtering:", len(df_sim))
    print("Filling a {}D histogram".format(dim))
    
    
    obs_bin_array = array('d', self.fobs_bins)
    jet_pt_bin_array = array('d', self.fptbinsJetFinalA)
    cand_pt_bin_array = array('d', self.fptbinsDA)
    
    if dim == 2:
        # Binning: x - shape, y - jet pt
        his2 = ROOT.TH2D("simulation","simulation", len(self.fobs_bins)-1, obs_bin_array, len(self.fptbinsJetFinalA)-1, jet_pt_bin_array)
        
        for index, row in df_sim.iterrows():
            his2.Fill(row[self.var], row['pt_jet'])
            his2.Scale(scale_factor)
        
        return his2
        
    if dim == 3:
        # Binning: x - shape, y - jet pt, z - dmeson pt
        his3 = ROOT.TH3D("simulation","simulation", len(self.fobs_bins)-1, obs_bin_array, len(self.fptbinsJetFinalA)-1, jet_pt_bin_array, len(self.fptbinsDA)-1, cand_pt_bin_array)
        
        for index, row in df_sim.iterrows():
            his3.Fill(row[self.var], row['pt_jet'], row['pt_cand'])
        return his3
        
    return None
  
  def charmprediction(self):
    print("===============Evaluating feed-down contribution")
    self.obs_axis = 3
    fout = ROOT.TFile.Open("AnalysisResults_charmprediction.root", 'recreate')
    
    for i in range(10):
    
        input_data = self.get_simulated_yields(self.simulationfile, 3, True, False)
        input_data_projectionx = input_data.ProjectionX()
        
        input_data_projectionx.Scale(1 / input_data_projectionx.Integral(), "width")
        input_data_projectionx.SetName("angvsobs{}".format(self.var))
        
        input_data.Scale(1 / input_data.Integral(), "width")
        input_data.SetName("angvsjetpt_obs{}".format(self.var))
        
        fout.cd()
        input_data.Write()
        input_data_projectionx.Write()
        
        self.obs_axis+=1
        
    fout.Close()
    print("end of prediction")
    
    
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

