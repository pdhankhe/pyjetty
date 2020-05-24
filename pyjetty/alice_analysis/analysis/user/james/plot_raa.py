#! /usr/bin/env python

import sys
import os
import argparse
from array import *
import numpy as np
import ROOT
import yaml

from pyjetty.alice_analysis.process.base import common_base
from pyjetty.alice_analysis.analysis.user.substructure import analysis_utils_obs

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class PlotRAA(common_base.CommonBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, config_file='', **kwargs):
    super(PlotRAA, self).__init__(**kwargs)
    
    self.utils = analysis_utils_obs.AnalysisUtils_Obs()
    
    self.initialize()

    print(self)
  
  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize(self):
      
    self.output_dir = '/Users/jamesmulligan/Analysis_theta_g/TheoryPredictions/'
    self.observables = ['theta_g']
    self.jetR_list = [0.2, 0.4]
    self.plot_data = False
    self.plot_theory = True
    
    self.colors = [600-6, 632-4]
    self.markers = [20, 21, 33]
    self.theory_colors  = [ROOT.kViolet-8, ROOT.kTeal-8, ROOT.kOrange+7, ROOT.kPink+1, 416-2, 1]
    
    self.obs_label = 'SD_zcut02_B0'
    self.formatted_grooming_label = 'SD #it{z}_{cut}=0.2, #beta=0'
    self.figure_approval_status = 'Preliminary'
    
    if not os.path.exists(self.output_dir):
      os.makedirs(self.output_dir)

  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def init_observable(self, observable, jetR):
  
    if jetR == 0.2:
      self.filedir_pp_name = '/Users/jamesmulligan/Analysis_theta_g/roounfold_output/pp_ref/67940-02/'
      self.filedir_AA_name = '/Users/jamesmulligan/Analysis_theta_g/roounfold_output/PbPb/78488-02/'
      self.min_pt = 60
      self.max_pt = 80
    elif jetR == 0.4:
      self.filedir_pp_name = '/Users/jamesmulligan/Analysis_theta_g/roounfold_output/pp_ref/67940-04/'
      self.filedir_AA_name = '/Users/jamesmulligan/Analysis_theta_g/roounfold_output/PbPb/78488-04/'
      self.min_pt = 80
      self.max_pt = 100
  
    self.file_pp_name = os.path.join(self.filedir_pp_name, '{}/final_results/fFinalResults.root'.format(observable))
    self.file_AA_name = os.path.join(self.filedir_AA_name, '{}/final_results/fFinalResults.root'.format(observable))

    if observable == 'theta_g':
      self.xtitle = '#theta_{g}'
    if observable == 'zg':
      self.xtitle = '#it{z}_{g}'
    self.ytitle = '#frac{{1}}{{#sigma_{{jet, inc}}}} #frac{{d#sigma}}{{d{}}}'.format(self.xtitle)
 
    self.main_result_name = 'hmain_{}_R{}_{}_{}-{}'.format(observable, jetR, self.obs_label, self.min_pt, self.max_pt)
    self.sys_total_name = 'hResult_{}_systotal_R{}_{}_{}-{}'.format(observable, jetR, self.obs_label, self.min_pt, self.max_pt)
    self.h_tagging_fraction_name = 'h_tagging_fraction_R{}_{}'.format(jetR, self.obs_label)
 
    self.file_pp = ROOT.TFile(self.file_pp_name, 'READ')
    self.file_AA = ROOT.TFile(self.file_AA_name, 'READ')

    self.h_main_pp = self.file_pp.Get(self.main_result_name)
    self.h_sys_pp = self.file_pp.Get(self.sys_total_name)
    self.h_tagging_fraction_pp = self.file_pp.Get(self.h_tagging_fraction_name)
    self.h_main_AA = self.file_AA.Get(self.main_result_name)
    self.h_sys_AA = self.file_AA.Get(self.sys_total_name)
    self.h_tagging_fraction_AA = self.file_AA.Get(self.h_tagging_fraction_name)
    
    integral_inclusive = self.h_main_pp.Integral(1, self.h_main_pp.GetNbinsX()+1, 'width')
    integral_reported = self.h_main_pp.Integral(2, self.h_main_pp.GetNbinsX()+1, 'width')
    tagging_fraction = integral_reported / integral_inclusive
    print('tagging fraction pp re-computed: {}'.format(tagging_fraction))
    
    integral_inclusive = self.h_main_AA.Integral(1, self.h_main_pp.GetNbinsX()+1, 'width')
    integral_reported = self.h_main_AA.Integral(2, self.h_main_pp.GetNbinsX()+1, 'width')
    tagging_fraction = integral_reported / integral_inclusive
    print('tagging fraction AA re-computed: {}'.format(tagging_fraction))
    
    pt = (self.min_pt + self.max_pt)/2
    bin = self.h_tagging_fraction_pp.GetXaxis().FindBin(pt)
    self.f_tagging_pp = self.h_tagging_fraction_pp.GetBinContent(bin)
    pt = (self.min_pt + self.max_pt)/2
    bin = self.h_tagging_fraction_AA.GetXaxis().FindBin(pt)
    self.f_tagging_AA = self.h_tagging_fraction_AA.GetBinContent(bin)

    self.ymax = self.h_main_pp.GetMaximum()
    
    if self.plot_theory:
      name = 'hRAA_{}_R{}_Theory.pdf'.format(observable, self.utils.remove_periods(jetR))
    else:
      name = 'hRAA_{}_R{}.pdf'.format(observable, self.utils.remove_periods(jetR))
    self.output_filename = os.path.join(self.output_dir, name)
    
    # Get theory predictions
    if self.plot_theory:
      if observable == 'theta_g':
        config_file = 'theory_predictions_rg.yaml'
      if observable == 'zg':
        config_file = 'theory_predictions_zg.yaml'
        
      # Get binning, for re-binning
      bin_edges = np.array(self.h_main_pp.GetXaxis().GetXbins())[1:]
      self.bin_array = array('d', bin_edges)

      with open(config_file, 'r') as stream:
        config = yaml.safe_load(stream)
              
      self.label_list = []
      self.observable_name_list = []
      self.sublabel_list = []
      self.prediction_g_list = []
      
      self.prediction_types = [name for name in list(config.keys())]
      for type in self.prediction_types:
      
        theory = config[type]
        configs = [name for name in list(theory.keys()) if 'config' in name ]
        for prediction in configs:
          theory_prediction = theory[prediction]
          config_jetR = theory_prediction['jetR']
          if config_jetR == jetR:
            
            plot_list = theory['plot_list']
            
            if type == 'lbnl':
            
              x = np.array(theory_prediction['x'])
              y_pp = np.array(theory_prediction['y_pp'])
              
              n = len(x)
              xerr = np.zeros(n)
              
              if 'f_q' in theory_prediction:
                f_q = theory_prediction['f_q']
                f_q_min = theory_prediction['f_q_min']
                f_q_max = theory_prediction['f_q_max']
                y_AA_quark = np.array(theory_prediction['y_AA_quark'])
                y_AA_gluon = np.array(theory_prediction['y_AA_gluon'])
                
                y_AA = f_q*y_AA_quark + (1-f_q)*y_AA_gluon
                y_AA_min = f_q_min*y_AA_quark + (1-f_q_min)*y_AA_gluon
                y_AA_max = f_q_max*y_AA_quark + (1-f_q_max)*y_AA_gluon
                
                ratio = np.divide(y_AA, y_pp)
                ratio_lower = np.divide(y_AA_min, y_pp)
                ratio_upper = np.divide(y_AA_max, y_pp)
                g = ROOT.TGraphAsymmErrors(n, x, ratio, xerr, xerr, ratio-ratio_lower, ratio_upper-ratio)

              else:
                y_AA = np.array(theory_prediction['y_AA'])
                ratio = np.divide(y_AA, y_pp)
                g = ROOT.TGraph(n, x, ratio)

            elif type == 'caucal':
                        
              x = np.array(theory_prediction['x'])
              ratio = np.array(theory_prediction['ratio'])
              ratio_neg_unc_tot = np.array(theory_prediction['ratio_neg_unc_tot'])
              ratio_pos_unc_tot = np.array(theory_prediction['ratio_pos_unc_tot'])
              
              n = len(x)
              xerr = np.zeros(n)
              g = ROOT.TGraphAsymmErrors(n, x, ratio, xerr, xerr, ratio_neg_unc_tot, ratio_pos_unc_tot)
            
            elif type == 'hybrid_model':
            
              xbins = np.array(theory_prediction['xbins'])
              
              # Get distributions
              y_pp = np.array(theory_prediction['y_pp'])
              y_pp_err = np.array(theory_prediction['y_pp_err'])
              y_AA_lower = np.array(theory_prediction['y_AA_lower'])
              y_AA_upper = np.array(theory_prediction['y_AA_upper'])
              
              # Rebin distributions
              x, y_pp = self.rebin_arrays(xbins, y_pp)
              _, y_pp_err = self.rebin_arrays(xbins, y_pp_err)
              _, y_AA_lower = self.rebin_arrays(xbins, y_AA_lower)
              _, y_AA_upper = self.rebin_arrays(xbins, y_AA_upper)

              # Form ratio and propagate uncertainty
              y_AA = (y_AA_lower + y_AA_upper) / 2.
              y_AA_lower_fraction = np.divide(y_AA - y_AA_lower, y_AA)
              y_AA_upper_fraction = np.divide(y_AA_upper - y_AA, y_AA)
              y_pp_err_fraction = np.divide(y_pp_err, y_pp)
              
              ratio = np.divide(y_AA, y_pp)
              ratio_lower_fraction = np.sqrt(np.square(y_AA_lower_fraction) + np.square(y_pp_err_fraction))
              ratio_upper_fraction = np.sqrt(np.square(y_AA_upper_fraction) + np.square(y_pp_err_fraction))
            
              ratio_lower = np.multiply(ratio, ratio_lower_fraction)
              ratio_upper = np.multiply(ratio, ratio_upper_fraction)
    
              print(x)
              n = len(x)
              xerr = np.zeros(n)
              g = ROOT.TGraphAsymmErrors(n, x, ratio, xerr, xerr, ratio_lower, ratio_upper)

            if prediction in plot_list:
              self.prediction_g_list.append(g)
              self.label_list.append(theory['label'])
              self.sublabel_list.append(theory_prediction['sublabel'])
              self.observable_name_list.append(theory['observable'])
  
  #---------------------------------------------------------------
  # Rebin numpy arrays (xbins,y) representing a histogram,
  # where x represents the bin edges, and y the bin content
  #
  # Return (x_rebinned, y_rebinned) where x_rebinned is the bin centers
  #
  # I don't find any easy way to do this...so I
  #        construct TH1, rebin, and then extract array)
  #----------------------------------------------------------------------
  def rebin_arrays(self, x, y):

    x_array = array('d', x)
    h = ROOT.TH1F('h', 'h', len(x)-1, x_array)
    for i, content in enumerate(y):
      h.SetBinContent(i+1, y[i])
    h_rebinned = h.Rebin(len(self.bin_array)-1, 'h_rebinned', self.bin_array)
    x_rebinned = []
    y_rebinned = []
    for i in range(len(self.bin_array)-1):
      x_rebinned.append(h_rebinned.GetBinCenter(i+1))
      y_rebinned.append(h_rebinned.GetBinContent(i+1))
      
    return (np.array(x_rebinned), np.array(y_rebinned))
    
  #---------------------------------------------------------------
  # This function is called once for each subconfiguration
  #----------------------------------------------------------------------
  def plot_raa(self):
  
    for observable in self.observables:
      for jetR in self.jetR_list:
  
        self.init_observable(observable, jetR)
        self.plot_final_result(observable, jetR)

  #---------------------------------------------------------------
  # This function is called once for each subconfiguration
  #----------------------------------------------------------------------
  def plot_final_result(self, observable, jetR):
    print('Plot final results for {}: R = {}, {}'.format(observable, jetR, self.obs_label))

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()

    c = ROOT.TCanvas('c', 'c', 600, 650)
    c.Draw()
    
    c.cd()
    pad1 = ROOT.TPad('myPad', 'The pad',0,0.4,1,1)
    pad1.SetLeftMargin(0.2)
    pad1.SetTopMargin(0.08)
    pad1.SetRightMargin(0.04)
    pad1.SetBottomMargin(0.)
    #pad1.SetLogy()
    pad1.Draw()
    pad1.cd()

    myLegend = ROOT.TLegend(0.25,0.65,0.45,0.85)
    self.utils.setup_legend(myLegend,0.05)
    
    self.h_main_pp.SetMarkerSize(1.5)
    self.h_main_pp.SetMarkerStyle(self.markers[0])
    self.h_main_pp.SetMarkerColor(self.colors[0])
    self.h_main_pp.SetLineStyle(1)
    self.h_main_pp.SetLineWidth(2)
    self.h_main_pp.SetLineColor(self.colors[0])
      
    self.h_sys_pp.SetLineColor(0)
    self.h_sys_pp.SetFillColor(self.colors[0])
    self.h_sys_pp.SetFillColorAlpha(self.colors[0], 0.3)
    self.h_sys_pp.SetFillStyle(1001)
    self.h_sys_pp.SetLineWidth(0)
    
    self.h_main_AA.SetMarkerSize(1.5)
    self.h_main_AA.SetMarkerStyle(self.markers[1])
    self.h_main_AA.SetMarkerColor(self.colors[1])
    self.h_main_AA.SetLineStyle(1)
    self.h_main_AA.SetLineWidth(2)
    self.h_main_AA.SetLineColor(self.colors[1])
      
    self.h_sys_AA.SetLineColor(0)
    self.h_sys_AA.SetFillColor(self.colors[1])
    self.h_sys_AA.SetFillColorAlpha(self.colors[1], 0.3)
    self.h_sys_AA.SetFillStyle(1001)
    self.h_sys_AA.SetLineWidth(0)
      
    xmin = self.h_sys_pp.GetBinLowEdge(2)
    xmax = self.h_sys_pp.GetBinLowEdge(self.h_sys_pp.GetNbinsX()+1)
    #if jetR == 0.4 and observable == 'theta_g':
    #  xmax = self.h_sys_pp.GetBinLowEdge(self.h_sys_pp.GetNbinsX()-1)
        
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, xmin, xmax)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetXTitle(self.xtitle)
    myBlankHisto.SetYTitle(self.ytitle)
    if observable == 'theta_g':
      if jetR == 0.2:
        ymax = 1.9*self.ymax
      else:
        ymax = 2.5*self.ymax
    else:
      ymax = 3.*self.ymax
    myBlankHisto.SetMaximum(ymax)
    myBlankHisto.SetMinimum(2e-4) # Don't draw 0 on top panel
    myBlankHisto.GetYaxis().SetTitleSize(0.065)
    myBlankHisto.GetYaxis().SetTitleSize(0.08)
    myBlankHisto.GetYaxis().SetTitleOffset(1.15)
    myBlankHisto.GetYaxis().SetLabelSize(0.06)
    myBlankHisto.Draw('E')
    
    if observable == 'theta_g':
      rg_axis_tf1 = ROOT.TF1('rg_axis_tf1', 'x', 0, jetR*xmax-0.01)
      rg_axis = ROOT.TGaxis(xmin, ymax, xmax, ymax, 'rg_axis_tf1', 505, '- S')
      rg_axis.SetTitle('#it{R}_{g}')
      rg_axis.SetTitleSize(25)
      rg_axis.SetTitleFont(43)
      rg_axis.SetTitleOffset(0.62)
      rg_axis.SetLabelFont(43)
      rg_axis.SetLabelSize(25)
      rg_axis.SetTickSize(0.025)
      rg_axis.SetLabelOffset(0.016)
      rg_axis.Draw()
 
    c.cd()
    pad2 = ROOT.TPad('pad2', 'pad2', 0, 0.02, 1, 0.4)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.3)
    pad2.SetLeftMargin(0.2)
    pad2.SetRightMargin(0.04)
    pad2.SetTicks(0,1)
    pad2.Draw()
    pad2.cd()
    
    ratio_legend = ROOT.TLegend(0.45,0.65,0.75,0.95)
    self.utils.setup_legend(ratio_legend,0.05)
          
    myBlankHisto2 = myBlankHisto.Clone('myBlankHisto_C')
    myBlankHisto2.SetYTitle('#frac{Pb-Pb}{pp}')
    myBlankHisto2.SetXTitle(self.xtitle)
    myBlankHisto2.GetXaxis().SetTitleSize(30)
    myBlankHisto2.GetXaxis().SetTitleFont(43)
    myBlankHisto2.GetXaxis().SetTitleOffset(3.)
    myBlankHisto2.GetXaxis().SetLabelFont(43)
    myBlankHisto2.GetXaxis().SetLabelSize(25)
    myBlankHisto2.GetYaxis().SetTitleSize(25)
    myBlankHisto2.GetYaxis().SetTitleFont(43)
    myBlankHisto2.GetYaxis().SetTitleOffset(2.)
    myBlankHisto2.GetYaxis().SetLabelFont(43)
    myBlankHisto2.GetYaxis().SetLabelSize(25)
    myBlankHisto2.GetYaxis().SetNdivisions(505)
    if observable == 'theta_g':
      if jetR == 0.2:
        myBlankHisto2.GetYaxis().SetRangeUser(0., 1.99)
      else:
        myBlankHisto2.GetYaxis().SetRangeUser(0., 2.1)
    else:
      myBlankHisto2.GetYaxis().SetRangeUser(0.4, 1.6)

    myBlankHisto2.Draw('')
  
    line = ROOT.TLine(xmin,1,xmax,1)
    line.SetLineColor(920+2)
    line.SetLineStyle(2)
    line.Draw()
      
    h_ratio = self.h_main_AA.Clone()
    h_ratio.Divide(self.h_main_pp)
      
    h_ratio_sys = self.h_sys_AA.Clone()
    h_ratio_sys.Divide(self.h_sys_pp)

    pad1.cd()
    if self.plot_data:
      self.h_sys_pp.DrawCopy('E2 same')
      self.h_sys_AA.DrawCopy('E2 same')
      self.h_main_pp.DrawCopy('PE X0 same')
      self.h_main_AA.DrawCopy('PE X0 same')
      
    pad2.cd()

    if self.plot_theory:
    
      for i, g in enumerate(self.prediction_g_list):

        color = self.theory_colors[i]
        g.SetLineColor(color)
        g.SetLineWidth(2)
        g.SetFillColor(color)
        #g.SetFillColorAlpha(color, 0.65)
        if type(g) == ROOT.TGraphAsymmErrors:
          g.Draw("3 same")
        elif type(g) == ROOT.TGraph:
          g.Draw("same")
        
        label = self.label_list[i]
        sublabel = self.sublabel_list[i]
        ratio_legend.AddEntry(g, '{}, {}'.format(label, sublabel), 'L')
     
    ratio_legend.Draw()
     
    if self.plot_data:
      h_ratio_sys.DrawCopy('E2 same')
      h_ratio.DrawCopy('PE X0 same')
   
    pad1.cd()
    myLegend.AddEntry(self.h_main_pp, 'pp', 'PE')
    myLegend.AddEntry(self.h_main_AA, 'Pb-Pb 0-10%', 'PE')
    myLegend.AddEntry(self.h_sys_pp, 'Sys. uncertainty', 'f')
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    
    text_latex.SetTextSize(0.05)
    text = '#it{{f}}_{{tagged}}^{{pp}} = {:.2f}, #it{{f}}_{{tagged}}^{{AA}} = {:.2f}'.format(self.f_tagging_pp, self.f_tagging_AA)
    text_latex.DrawLatex(0.56, 0.39, text)

    text_latex.SetTextSize(0.06)
    text = 'ALICE {}'.format(self.figure_approval_status)
    text_latex.DrawLatex(0.56, 0.84, text)
    
    text = '#sqrt{#it{s_{#it{NN}}}} = 5.02 TeV'
    text_latex.DrawLatex(0.56, 0.77, text)

    text = 'Charged jets   anti-#it{k}_{T}'
    text_latex.DrawLatex(0.56, 0.7, text)
    
    text = '#it{R} = ' + str(jetR) + '   | #eta_{jet}| < 0.5'
    text_latex.DrawLatex(0.56, 0.63, text)
    
    text = str(self.min_pt) + ' < #it{p}_{T, ch jet} < ' + str(self.max_pt) + ' GeV/#it{c}'
    text_latex.DrawLatex(0.56, 0.56, text)
    
    text = self.formatted_grooming_label
    text_latex.DrawLatex(0.56, 0.47, text)
    
    myLegend.Draw()

    c.SaveAs(self.output_filename)
    c.Close()
    
    # Write result to ROOT file
    #final_result_root_filename = os.path.join(output_dir, 'fFinalResults.root')
    #fFinalResults = ROOT.TFile(final_result_root_filename, 'UPDATE')
    #h.Write()
    #h_sys.Write()
    #hPythia.Write()
    #fFinalResults.Close()

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
  
  analysis = PlotRAA(config_file = args.configFile)
  analysis.plot_raa()
