#! /usr/bin/env python

import sys
import os
import argparse
from array import *
import numpy as np
import ROOT
import yaml

from pyjetty.alice_analysis.analysis.user.substructure import run_analysis
# from pyjetty.alice_analysis.analysis.user.james import run_analysis_james_base
# from pyjetty.alice_analysis.analysis.user.james import plotting_utils_theta_g

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
# class RunAnalysisThetaG(run_analysis_james_base.RunAnalysisJamesBase):
class RunAnalysisdEEC(run_analysis.RunAnalysis):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, config_file='', **kwargs):
    super(RunAnalysisdEEC, self).__init__(config_file, **kwargs)
    
    print(self)

  #---------------------------------------------------------------
  # This function is called once for each subconfiguration
  #---------------------------------------------------------------
  def plot_single_result(self, jetR, obs_label, obs_setting, grooming_setting):
    print('Plotting each individual result...')
  
    # Plot final result for each 1D substructure distribution (with PYTHIA)
    self.plot_final_result(jetR, obs_label, obs_setting, grooming_setting)
    
    
  
  #---------------------------------------------------------------
  # This function is called once after all subconfigurations have been looped over, for each R
  #---------------------------------------------------------------
  def plot_all_results(self, jetR):
    print('Plotting overlay of all results...')
    
    for i_config, overlay_list in enumerate(self.plot_overlay_list):
    
      if len(overlay_list) > 1:
      
        self.plot_final_result_overlay(i_config, jetR, overlay_list)

        self.plot_NPcorrections(i_config, jetR, overlay_list)
  

  
  #----------------------------------------------------------------------
  def plot_final_result(self, jetR, obs_label, obs_setting, grooming_setting):
    print('Plot final results for {}: R = {}, {}'.format(self.observable, jetR, obs_label))

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()
    
    # Construct histogram of tagging fraction, to write to file
    if 'sd' in grooming_setting:
      name = 'h_tagging_fraction_R{}_{}'.format(jetR, obs_label)
      h_tagging = ROOT.TH1D(name, name, len(self.pt_bins_reported) - 1, array('d', self.pt_bins_reported))

    # Loop through pt slices, and plot final result for each 1D theta_g distribution
    for bin in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[bin]
      max_pt_truth = self.pt_bins_reported[bin+1]
      
      self.plot_observable(jetR, obs_label, obs_setting, grooming_setting, min_pt_truth, max_pt_truth, plot_pythia=True)
      
      # Fill tagging fraction
      if 'sd' in grooming_setting:
        fraction_tagged = getattr(self, 'tagging_fraction_R{}_{}_{}-{}'.format(jetR, obs_label, min_pt_truth, max_pt_truth))
        pt = (min_pt_truth + max_pt_truth)/2
        h_tagging.Fill(pt, fraction_tagged)
      
    # Write tagging fraction to ROOT file
    if 'sd' in grooming_setting:
      output_dir = getattr(self, 'output_dir_final_results')
      final_result_root_filename = os.path.join(output_dir, 'fFinalResults.root')
      fFinalResults = ROOT.TFile(final_result_root_filename, 'UPDATE')
      h_tagging.Write()
      fFinalResults.Close()
   
  #----------------------------------------------------------------------




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

  analysis = RunAnalysisdEEC(config_file = args.configFile)
  analysis.run_analysis()
