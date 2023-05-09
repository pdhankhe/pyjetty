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


################################################################
# Helper functions
################################################################

def process_histo(h, size, col, style):
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    h.SetMarkerSize(size)
    h.SetMarkerColor(col)
    h.SetLineWidth(2)
    h.SetLineColor(col)
    h.SetMarkerStyle(style)
    h.GetYaxis().SetTitleOffset(1.)
    h.GetYaxis().SetTitleSize(0.042)
    h.GetYaxis().SetLabelSize(0.042)
    h.GetYaxis().SetLabelFont(42)
    h.GetXaxis().SetLabelFont(42)
    h.GetYaxis().SetTitleFont(42)
    h.GetXaxis().SetTitleFont(42)
    h.GetXaxis().SetTitleOffset(1.0)
    h.GetXaxis().SetTitleSize(0.042)
    h.GetXaxis().SetLabelSize(0.042)

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


def ErrorCalUnCorr(x,  y,  dx,  dy):
    return (x/y)*math.sqrt((dx/x)*(dx/x)+(dy/y)*(dy/y))

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
    self.hf_jet_shape_extraction()
    
    process_base.ProcessBase.save_output_objects(self)
    print(self)
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
    self.use_MC_gauss_info = config['use_MC_gauss_info']
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
    self.rebin_theory_response = config['rebin_theory_response']

    self.datafile = config['main_rawdatafile']
    print(self.datafile)
    self.simulationfile = config['main_feeddown_simulation']
    print(self.simulationfile)
    self.responsefile = config['main_response']
    self.gauss_info = config['main_gauss_info']
    self.output_dir_fd = config['output_dir_fd']
    self.efffile = config['main_efffile']
    self.nonprompt_kinematic_efffile = config['main_nonprompt_kinematic_efffile']
    print(self.efffile)
    self.reflection_template = config['reflection_template']
    self.isRefSys = config['isRefSys']
    self.refScale = config['refScale']
    print(self.reflection_template[3])
    

    # Whether or not to use the previous preliminary result in final plots
    self.use_prev_prelim = config['use_prev_prelim']

    self.histutils = ROOT.RUtil.HistUtils()

    # ---------------------------------------------------------------
    # This function is used for jet signal extraction
    # ---------------------------------------------------------------

  def hf_jet_shape_extraction(self):
    file_data = ROOT.TFile(self.datafile, 'READ')
    jetR = self.jetR_list[0]
    sparse = file_data.Get("hsparsejet")
 
    okSignalExt = self.rawJetSpectra(sparse, jetR)
    if not okSignalExt:
        print("!!!!!! Something wrong in the raw signal extraction !!!!")
        
            
    # ---------------------------------------------------------------
    # This function is used for raw jet signal extraction
    # ---------------------------------------------------------------

  def rawJetSpectra(self, sparse, jetR):

    obs_bin_array = array('d', self.fobs_bins)
    jet_pt_bin_array = array('d', self.fptbinsJetFinalA)
    
    ndim = 4
    
    hangvsjetpt_sb_sub = ROOT.TH2F("hangvsjetpt_sb_sub", "hangvsjetpt_sb_sub",  len(self.fptbinsJetFinalA)-1, jet_pt_bin_array, len(self.fobs_bins)-1, obs_bin_array)
    hangvsjetpt_eff_corr = ROOT.TH2F("hangvsjetpt_eff_corr", "hangvsjetpt_eff_corr",  len(self.fptbinsJetFinalA)-1, jet_pt_bin_array, len(self.fobs_bins)-1, obs_bin_array)
    hangvsjetpt_woeff_corr = ROOT.TH2F("hangvsjetpt_woeff_corr", "hangvsjetpt_woeff_corr",  len(self.fptbinsJetFinalA)-1, jet_pt_bin_array, len(self.fobs_bins)-1, obs_bin_array)
    
    hangvsjetpt_sb_sub.Sumw2()
    hangvsjetpt_eff_corr.Sumw2()
    hangvsjetpt_woeff_corr.Sumw2()
    
    file_gauss_info = ROOT.TFile(self.gauss_info, 'READ')
    print(self.gauss_info)
    
    for beta in self.beta_list:
        
        if beta == 1:
            obs_axis = 3
        if beta == 1.5:
            obs_axis = 4
        if beta == 2:
            obs_axis = 5
        if beta == 3:
            obs_axis = 6
            
        print("calculating shape for the observable ", obs_axis)
        #loop over jet pt
        first_jetpt = 0
        
        hang_sb_sub_rebinned_list = []
        hang_eff_corr_rebinned_list = []
        hang_woeff_corr_rebinned_list = []
        
        for ibin2 in range(len(self.fptbinsJetFinalA)-1):
            
            dim = array('i',(2, obs_axis, 1, 0))
            #Clone sparse info
            suffix_jet = str(self.fptbinsJetFinalA[ibin2])+" < p_{T,jet} < "+str(self.fptbinsJetFinalA[ibin2+1])
            first_fit = 0
            first_fit_effcor = 0
            sparsehadron = sparse.Clone("sparsehadron"+suffix_jet)
    
            print("Analysis in jet pt bins = "+str(self.fptbinsJetFinalA[ibin2])+"-"+str(self.fptbinsJetFinalA[ibin2+1]))
            sparsehadron.GetAxis(0).SetRangeUser(self.fptbinsJetFinalA[ibin2], self.fptbinsJetFinalA[ibin2+1])
            #0-D-mass, 1-shape, 2-Dmeson pt, 3-Jet pt
           
            h_invmass_shape_Dmesonpt = sparsehadron.Projection(4, dim,"e")
            h_invmass_shape_Dmesonpt.Sumw2()
            h_invmass_shape_Dmesonpt.GetAxis(0).SetTitle("Counts per 5 MeV/c^{2}")
            h_invmass_shape_Dmesonpt.GetAxis(0).SetTitle("m(K#pi)(GeV/c^{2})")
            h_invmass_shape_Dmesonpt.GetAxis(0).SetTitleSize(0.06)
            h_invmass_shape_Dmesonpt.GetAxis(0).SetTitleOffset(0.7)
            h_invmass_shape_Dmesonpt.GetAxis(0).SetTitleOffset(0.7)
            h_invmass_shape_Dmesonpt.SetTitle("")
    
            #for systematics get information of mean and sigma form MC
            hmean = file_gauss_info.Get("mean_jetpt_{}_{}".format(self.fptbinsJetFinalA[ibin2], self.fptbinsJetFinalA[ibin2+1]))
            hisgma = file_gauss_info.Get("sigma_jetpt_{}_{}".format(self.fptbinsJetFinalA[ibin2], self.fptbinsJetFinalA[ibin2+1]))
            
            RS = 0
            input_data_invmass_list = []
            hsignal_sb_sub_rebinned_list = []
            hmass_l_list = []
            hmass_u_list = []
            hmass_c_list = []
            bkg_func_list = []
            
            cinvmass = ROOT.TCanvas("cinvmass ", "cinvmass")
            cinvmass.Divide(4,3)
            
            #loop over cand. pt
            for ipt in range(self.fptbinsDN):
              #max pt cut on Dmeson in a jet bin
              if self.fptbinsDA[ipt+1] > self.fptbinsJetFinalA[ibin2+1]:
                break
                    
              print("min dmeson pt cut"+str(self.pt_dmeson_minpt_cut[ibin2]))
              #minimum pt cut on Dmeson in a jet bin
              if self.fptbinsDA[ipt] < self.pt_dmeson_minpt_cut[ibin2]:
                input_data_invmass_list.append(None)
                hsignal_sb_sub_rebinned_list.append(None)
                hmass_l_list.append(None)
                hmass_u_list.append(None)
                hmass_c_list.append(None)
                bkg_func_list.append(None)
                continue
                
              print("Analysis in Dmeson pt bins = "+str(self.fptbinsDA[ipt])+"-"+str(self.fptbinsDA[ipt+1]))
              text_for_invmass = str(self.fptbinsJetFinalA[ibin2])+" < p_{T,jet} < "+str(self.fptbinsJetFinalA[ibin2+1])+", "+\
                                 str(self.fptbinsDA[ipt])+" < #it{p}_{T}^{D^{0}} < "+str(self.fptbinsDA[ipt + 1])
              suffix = text_for_invmass
              ptbin = (self.fptbinsDA[ipt]+self.fptbinsDA[ipt+1]) / 2.
              h_shape_Dmesonpt_jetpt = h_invmass_shape_Dmesonpt.Clone("shapevsdptvsjetpt"+suffix)
              h_invmass_shape_Dmesonpt.GetAxis(2).SetRangeUser(self.fptbinsDA[ipt], self.fptbinsDA[ipt+1])
              h_invmass_shape_Dmesonpt.GetAxis(3).SetRangeUser(self.fptbinsJetFinalA[ibin2], self.fptbinsJetFinalA[ibin2+1])
              hinvmass = h_invmass_shape_Dmesonpt.Projection(0,"e")
              hinvmass.SetName("invmass"+suffix_jet)
        
              hinvmass.Sumw2()
              hinvmass.Rebin(self.RebinMass)
              hinvmass.SetTitle(text_for_invmass)
              hinvmass.GetYaxis().SetTitle("Counts per 5 MeV/c^{2}")
              hinvmass.GetYaxis().SetTitleSize(0.06)
              hinvmass.GetYaxis().SetTitleOffset(0.7)
              hinvmass.GetXaxis().SetRangeUser(self.massrange[0], self.massrange[1])
              

              hmassfit = hinvmass.Clone()
              hmassfit.SetMaximum(hmassfit.GetMaximum() * 1.3)
              hmin = ROOT.TMath.Max(self.massrange[0], hmassfit.GetBinLowEdge(2))
              hmax = ROOT.TMath.Min(self.massrange[1], hmassfit.GetBinLowEdge(hmassfit.GetNbinsX()))

              histo_invmassfit = ROOT.TH1F()
              histo_invmassfit.Sumw2()
              hmassfit.Copy(histo_invmassfit)

              fitterp = ROOT.RUtil.AliHFInvMassFitter(histo_invmassfit, hmin, hmax, self.bkgtype,0)
              if self.use_MC_gauss_info:
                #use fix pdg mass
                #fitterp.SetFixGaussianMean(self.Dmass)
                #use fix mean mass from MC
                fitterp.SetFixGaussianMean(hmean.GetBinContent(hmean.GetXaxis().FindBin(ptbin)))
                print("mean value used from MC ", hmean.GetBinContent(hmean.GetXaxis().FindBin(ptbin)))
              else:
                fitterp.SetInitialGaussianMean(self.Dmass)
                
              fitterp.SetInitialGaussianSigma(self.Dsigma)
              
              if self.IsUseRefl:
                self.SetReflection(fitterp, hmin, hmax, self.fptbinsJetFinalA[ibin2], self.fptbinsJetFinalA[ibin2+1], self.fptbinsDA[ipt], self.fptbinsDA[ipt + 1])
                print("RS value calculated"+str(self.RS))

              fitterp.MassFitter(0)
            

              h = fitterp.GetHistoClone()
              bkg_func = fitterp.GetBackgroundRecalcFunc()
              
              if not bkg_func: # if there is no background fit it continues
                print('FAIL: bkg. func. not found')
                print("exiting analysis in Dmeson pt bins = "+str(self.fptbinsDA[ipt])+"-"+str(self.fptbinsDA[ipt+1]))
                continue
                
              bkg_func.SetRange(hmin, hmax)
              bkg_func.SetName("bkground" + suffix)
              bkg_func.SetLineColor(633)
              
              bkg_func_list.append(bkg_func)
              input_data_invmass_list.append(histo_invmassfit)
              
              mass_func = fitterp.GetMassFunc()
              mass_func.SetRange(hmin, hmax)
              mass_func.SetLineColor(4)

              if self.IsUseRefl:
                ref_func = fitterp.GetReflFunc()
                ref_func.SetLineColor(418)
              else:
                None

              sig_func = fitterp.GetSignalFunc()
              sig_func.SetRange(hmin, hmax)

              fullfit = h.GetFunction("funcmass")
              hmass = histo_invmassfit.Clone()

              if self.IsUseRefl:
                bkgRfit = fitterp.GetBkgPlusReflFunc()
                bkgRfit.SetRange(hmin, hmax)
                bkgRfit.SetLineColor(15)
             
              if not fullfit:
                print("======= Fit failed for bin: ")
                continue
                

              Dsigma = fitterp.GetSigma()
              DsigmaUnc = fitterp.GetSigmaUncertainty()
              Dmean = fitterp.GetMean()
              DmeanUnc = fitterp.GetMeanUncertainty()


              signal_c_min = Dmean-self.Signalsigma*Dsigma
              signal_c_max = Dmean+self.Signalsigma*Dsigma
              signal_l_min = Dmean+self.Bkgsigma[0]*Dsigma
              signal_l_max = Dmean+self.Bkgsigma[1]*Dsigma
              signal_u_min = Dmean+self.Bkgsigma[2]*Dsigma
              signal_u_max = Dmean+self.Bkgsigma[3]*Dsigma
              
              if signal_l_min<hmin:
                signal_l_min = hmin
              if signal_u_max>hmax:
                signal_u_max = hmax
                
              binwidth = hmass.GetXaxis().GetBinWidth(1)*0.5
              binmin = hmass.GetXaxis().FindBin(signal_c_min)
              binmax = hmass.GetXaxis().FindBin(signal_c_max)
              min = hmass.GetXaxis().GetBinCenter(binmin)-binwidth
              max = hmass.GetXaxis().GetBinCenter(binmax-1)+binwidth
              
              binmin_sb1 = hmass.GetXaxis().FindBin(signal_l_min)
              binmax_sb1 = hmass.GetXaxis().FindBin(signal_l_max)
              min_sb1 = hmass.GetXaxis().GetBinCenter(binmin_sb1)-binwidth
              max_sb1 = hmass.GetXaxis().GetBinCenter(binmax_sb1-1)+binwidth
              
              binmin_sb2 = hmass.GetXaxis().FindBin(signal_u_min)
              binmax_sb2 = hmass.GetXaxis().FindBin(signal_u_max)
              min_sb2 = hmass.GetXaxis().GetBinCenter(binmin_sb2)-binwidth
              max_sb2 = hmass.GetXaxis().GetBinCenter(binmax_sb2-1)+binwidth
              
              nbkg = bkg_func.Integral(signal_c_min, signal_c_max) / binwidth
              nsignal = mass_func.Integral(signal_c_min, signal_c_max) / binwidth - nbkg
      
              srelerr = 0.0
              ref = 0.0
              ref1 = 0.0
              ref2 = 0.0

              signal_info = fitterp.Signal(min, max)
              bkg_info = fitterp.Background(min, max)

              signal= signal_info[0]
              signal_err = signal_info[1]

              bkg = bkg_info[0]
              bkg_err = bkg_info[1]

              sb1 = bkg_func.Integral(min_sb1, max_sb1) / hmass.GetBinWidth(1)
              sb2 = bkg_func.Integral(min_sb2, max_sb2) / hmass.GetBinWidth(1)
              
              if sb1+sb2 == 0:
                continue

              if self.IsUseRefl:
                bkgref = bkgRfit.Integral(min, max) / hmass.GetBinWidth(1)
                ref = bkgref - bkg
                ref1 = (bkgRfit.Integral(min_sb1, max_sb1) / hmass.GetBinWidth(1)) - sb1
                ref2 = (bkgRfit.Integral(min_sb2, max_sb2) / hmass.GetBinWidth(1)) - sb2

              signf = 0
              signferr = 0
              sob = 0
              soberr = 0

              significance_info =fitterp.Significance(min, max)
              signf = significance_info[0]
              signferr = significance_info[1]


              if signal:
                srelerr = signal_err/ signal
              if bkg:
                sob = signal / bkg
              else:
                sob = signal

              if bkg and bkg_err:
                soberr = ROOT.TMath.Sqrt((signal_err / bkg) * (signal_err / bkg) + (signal / bkg / bkg * bkg_err) * (signal / bkg / bkg * bkg_err))
              else:
                soberr = signal_err

              #---------------- side-band drawing
              hmass_l = hmass.Clone("hmass_l" + suffix)
              hmass_l.GetXaxis().SetRangeUser(signal_l_min,signal_l_max)
           
              hmass_u = hmass.Clone("hmass_u" + suffix)
              hmass_u.GetXaxis().SetRangeUser(signal_u_min,signal_u_max)
           
              hmass_c = hmass.Clone("hmass_c" + suffix)
              hmass_c.GetXaxis().SetRangeUser(signal_c_min,signal_c_max)
    

              hmass_l.SetFillColor(602)
              hmass_u.SetFillColor(602)
              hmass_c.SetFillColor(633)
              hmass_l.SetFillStyle(3004)
              hmass_u.SetFillStyle(3004)
              hmass_c.SetFillStyle(3005)

              hmass_l_list.append(hmass_l)
              hmass_u_list.append(hmass_u)
              hmass_c_list.append(hmass_c)
              
              cinvmass.cd(ipt+1)
              input_data_invmass_list[ipt].Draw("same,ep")
              bkg_func_list[ipt].Draw("same")
              hmass_l_list[ipt].Draw("hsame")
              hmass_u_list[ipt].Draw("hsame")
              hmass_c_list[ipt].Draw("hsame")
              mass_func.Draw("same")
              if self.IsUseRefl:
                bkgRfit.Draw("same")
                ref_func.Draw("same")
              
              # ---------------- clone sparse for signal shape extraction
              h_shape_Dmesonpt_jetptleft = h_shape_Dmesonpt_jetpt.Clone("signalleft"+suffix)
              h_shape_Dmesonpt_jetptright = h_shape_Dmesonpt_jetpt.Clone("signalright"+suffix)
            
              # ---------------- Extracting signal shape
              h_shape_Dmesonpt_jetpt.GetAxis(0).SetRangeUser(signal_c_min,signal_c_max)
              h_shape_Dmesonpt_jetpt.GetAxis(2).SetRangeUser(self.fptbinsDA[ipt], self.fptbinsDA[ipt + 1])
              h_shape_Dmesonpt_jetpt.GetAxis(3).SetRangeUser(self.fptbinsJetFinalA[ibin2],self.fptbinsJetFinalA[ibin2+1])
              hang_sig_raw = h_shape_Dmesonpt_jetpt.Projection(1,3,"e")
              hang_sig_raw.Sumw2()
              hang_sig = self.histutils.rebin_th2(hang_sig_raw,"hang_sig"+suffix, jet_pt_bin_array, len(self.fptbinsJetFinalA)-1, obs_bin_array, len(self.fobs_bins)-1)
              ROOT.RUtil.delete_h(hang_sig_raw)
              
              
              # ---------------- Extracting left SB shape
              h_shape_Dmesonpt_jetptleft.GetAxis(0).SetRangeUser(signal_l_min,signal_l_max)
              h_shape_Dmesonpt_jetptleft.GetAxis(2).SetRangeUser(self.fptbinsDA[ipt], self.fptbinsDA[ipt + 1])
              h_shape_Dmesonpt_jetptleft.GetAxis(3).SetRangeUser(self.fptbinsJetFinalA[ibin2],self.fptbinsJetFinalA[ibin2+1])
              hang_sig_sb1_raw = h_shape_Dmesonpt_jetptleft.Projection(1,3,"e")
              hang_sig_sb1_raw.Sumw2()
              hang_sig_sb1 = self.histutils.rebin_th2(hang_sig_sb1_raw,"hang_sig_sb1"+suffix, jet_pt_bin_array, len(self.fptbinsJetFinalA)-1, obs_bin_array, len(self.fobs_bins)-1)
              ROOT.RUtil.delete_h(hang_sig_sb1_raw)

              # ---------------- Extracting right SB shape
              h_shape_Dmesonpt_jetptright.GetAxis(0).SetRangeUser(signal_u_min,signal_u_max)
              h_shape_Dmesonpt_jetptright.GetAxis(2).SetRangeUser(self.fptbinsDA[ipt], self.fptbinsDA[ipt + 1])
              h_shape_Dmesonpt_jetptright.GetAxis(3).SetRangeUser(self.fptbinsJetFinalA[ibin2],self.fptbinsJetFinalA[ibin2+1])
              hang_sig_sb2_raw = h_shape_Dmesonpt_jetptright.Projection(1,3,"e")
              hang_sig_sb2_raw.Sumw2()
              hang_sig_sb2 = self.histutils.rebin_th2(hang_sig_sb2_raw,"hang_sig_sb2"+suffix, jet_pt_bin_array, len(self.fptbinsJetFinalA)-1, obs_bin_array, len(self.fobs_bins)-1)
              ROOT.RUtil.delete_h(hang_sig_sb2_raw)
              
              hang_bkg = hang_sig_sb1.Clone("hang_sb_bkg"+ suffix)
              hang_bkg.Add(hang_sig_sb2)

              # Scaling - - Reflection signal OFF
              scalingB = bkg / (sb1 + sb2)
              print("scaling B is"+str(scalingB))
              
              # Scaling - - Reflection signal ON
              if self.IsUseRefl:
                scalingS = signal / (signal + ref - ((ref1 + ref2) * bkg) / (sb1 + sb2))
                print("scaling S is"+str(scalingS))
                
              hang_bkg_scaled = hang_bkg.Clone("hbkg_scaled" + suffix)
              hang_bkg_scaled.Scale(scalingB)
              
              # ------- subtract background from signal jet
              hang_sig_bkgsub = hang_sig.Clone("hang_sig_bkgsub" + suffix)
              hang_sig_bkgsub.Sumw2()
              hang_sig_bkgsub.Add(hang_bkg_scaled, -1)
            
              if self.IsUseRefl:
                hang_sig_bkgsub.Scale(scalingS)
                
              if self.Signalsigma == 2:
                print("Scaling for limited 2σ width of the signal region")
                hang_sig_bkgsub.Scale(1/ 0.9545)
                
              hsignal_sb_sub_rebinned_list.append(hang_sig_bkgsub)
             
              #process_histo(hsignal_sb_sub_rebinned_list[ipt],0.8, 1 ,20)
              #hsignal_sb_sub_rebinned_list[ipt].Draw("same,text")
              
              # setting the negative bins to zero
              for ibinx in range(hang_sig_bkgsub.GetNbinsX()):
                for ibiny in range(hang_sig_bkgsub.GetNbinsY()):
                    content = hang_sig_bkgsub.GetBinContent(ibinx+1,ibiny+1)
                    if content < 0:
                        hang_sig_bkgsub.SetBinContent(ibinx+1,ibiny+1,0)
                        hang_sig_bkgsub.SetBinError(ibinx+1,ibiny+1,0)
              
            
              # ------- sum all Dmeson pt bins
              if first_fit == 0:
                hrawjetptspectrum = hang_sig_bkgsub.Clone("hrawjetptspectrum"+suffix)
                first_fit = 1
                hrawjetptspectrum.Sumw2()
              else:
                hrawjetptspectrum.Add(hang_sig_bkgsub)
                
             
              
              # ------- correct shape for prompt efficiency
              hang_sig_bkgsub_effcorr = hang_sig_bkgsub.Clone("hang_sig_bkgsub_effcorr" + suffix)
              self.efficiency()
              print("efficiency values" + str(self.efficiency_prompt.GetBinContent(self.efficiency_prompt.GetXaxis().FindBin(ptbin))))
              eff_value = self.efficiency_prompt.GetBinContent(self.efficiency_prompt.GetXaxis().FindBin(ptbin))
            
              hang_sig_bkgsub_effcorr.Scale(1/eff_value)
           
              # ------- sum all Dmeson pt bins
              if first_fit_effcor == 0:
                hang_corr = hang_sig_bkgsub_effcorr.Clone("hjetptspectrum"+suffix)
                hang_corr.Sumw2()
                first_fit_effcor = 1
              else:
                hang_corr.Add(hang_sig_bkgsub_effcorr)
                
             
            
            hang_sb_sub_rebinned_list.append(hrawjetptspectrum)
            hang_eff_corr_rebinned_list.append(hang_corr)
            cinvmass.SaveAs("invmass_{}_{}.pdf".format(self.fptbinsJetFinalA[ibin2],self.fptbinsJetFinalA[ibin2+1]))
            
            if ibin2 == 0:
                hangvsjetpt_sb_sub = hang_sb_sub_rebinned_list[ibin2].Clone("hangvsjet_sb_sub"+suffix_jet)
                hangvsjetpt_sb_sub.Sumw2()
                hangvsjetpt_eff_corr = hang_eff_corr_rebinned_list[ibin2].Clone("hangvsjet_sb_sub_eff"+suffix_jet)
                hangvsjetpt_eff_corr.Sumw2()
            else:
                hangvsjetpt_sb_sub.Add(hang_sb_sub_rebinned_list[ibin2])
                hangvsjetpt_eff_corr.Add(hang_eff_corr_rebinned_list[ibin2])
                
 
        hfeeddown = self.feeddown(obs_axis, jetR, beta)
        hangvsjetpt_feedsub = hangvsjetpt_eff_corr.Clone("h_ang_R%s_%s" % (jetR, beta))
        angvsjetpt_feedsub_proj_y = hangvsjetpt_feedsub.ProjectionX()
        hangvsjetpt_feedsub.Add(hfeeddown,-1)
        
        # setting the negative bins to zero
        for ibinx in range(hangvsjetpt_feedsub.GetNbinsX()):
            for ibiny in range(hangvsjetpt_feedsub.GetNbinsY()):
                content = hangvsjetpt_feedsub.GetBinContent(ibinx+1,ibiny+1)
                if content < 0:
                    hangvsjetpt_feedsub.SetBinContent(ibinx+1,ibiny+1,0)
                    hangvsjetpt_feedsub.SetBinError(ibinx+1,ibiny+1,0)
        
        name = "h_ang_JetPt_R%s_%s" % (jetR, beta)
        hangvsjetpt_feedsub.SetName(name)
        setattr(self, name, hangvsjetpt_feedsub)
        
        #fout = ROOT.TFile("AnalysisResults_output_{}.root".format(beta), 'recreate')
        #hangvsjetpt_eff_corr.Write()
        #hangvsjetpt_feedsub.Write()
        #fout.Close()


    #print("length of histogram list"+str(len(input_data_invmass_list)))
    
    
    return True

  def efficiency(self):
    file_eff = ROOT.TFile(self.efffile, 'READ')
    print(self.efffile)
    self.efficiency_prompt = file_eff.Get("PromptEfficiency")
    self.efficiency_prompt.SetDirectory(0)
    self.efficiency_nonprompt = file_eff.Get("NonPromptEfficiency")
    self.efficiency_nonprompt.SetDirectory(0)
    
    file_nonprompt_kinematic_eff = ROOT.TFile(self.nonprompt_kinematic_efffile,'READ')
    self.particle_level_kinematic_eff = file_nonprompt_kinematic_eff.Get("h_nonprompt_kinematic_particle_level_efficiency")
    self.particle_level_kinematic_eff.SetDirectory(0)
    self.gen_level_kinematic_eff = file_nonprompt_kinematic_eff.Get("h_nonprompt_kinematic_gen_level_efficiency")
    self.gen_level_kinematic_eff.SetDirectory(0)
    return True

  def SetReflection(self, fitter, fLeftFitRange, fRightFitRange, jetptmin, jetptmax, Dptmin, Dptmax):
    file_ref = ROOT.TFile(self.reflection_template[3], 'READ')
    if not file_ref:
        print("File "+str(reflection_template)+" (reflection templates) cannot be opened! check your file path!")
        return False
    
    histRefl = file_ref.Get("histRflFittedDoubleGaus_jetpt_{}_{}_Dmesonpt_{}_{}".format(jetptmin, jetptmax, Dptmin, Dptmax))
    histSign = file_ref.Get("histSgn_jetpt_{}_{}_Dmesonpt_{}_{}".format(jetptmin, jetptmax, Dptmin, Dptmax))
    
    
    if not histRefl or not histSign:
        print("Error in loading the template/signal histrograms! Exiting...")
        return False
        
    fitter.SetTemplateReflections(histRefl,"template",fLeftFitRange,fRightFitRange)
  
    
    histo_ref = ROOT.TH1F()
    histRefl.Copy(histo_ref)
    
    fitter.SetHistogramReflectionFit(histo_ref)
        
    RoverS=histRefl.Integral(histRefl.FindBin(fLeftFitRange),histRefl.FindBin(fRightFitRange))/histSign.Integral(histSign.FindBin(fLeftFitRange),histSign.FindBin(fRightFitRange))
    
    if self.isRefSys:
        RoverS*=self.refScale
 
        
    print("R/S ratio in fit range for bin jetpt {} {} and Dmesonpt {} {} = {}".format(jetptmin, jetptmax, Dptmin, Dptmax, RoverS))
    fitter.SetFixReflOverS(RoverS)
    self.RS = RoverS
    #RoverS
    print("value of reflection over signal"+str(self.RS))
    return True
    

  def get_simulated_yields(self, file_path, dim, prompt, xsec, obs_axis):
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
    
    if obs_axis == 3:
        self.var = "LB1K1_jet"
    if obs_axis == 4:
        self.var = "LB15K1_jet"
    if obs_axis == 5:
        self.var = "LB2K1_jet"
    if obs_axis == 6:
        self.var = "LB3K1_jet"

 
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
        # Binning: x - jet pt, y - shape
        his2 = ROOT.TH2D("simulation","simulation", len(self.fptbinsJetFinalA)-1, jet_pt_bin_array, len(self.fobs_bins)-1, obs_bin_array)
        
        for index, row in df_sim.iterrows():
            his2.Fill(row['pt_jet'], row[self.var])
            his2.Scale(scale_factor)
        
        return his2
        
    if dim == 3:
        # Binning: x - jet pt, y - shape , z - dmeson pt
        his3 = ROOT.TH3D("simulation","simulation",  len(self.fptbinsJetFinalA)-1, jet_pt_bin_array, len(self.fobs_bins)-1, obs_bin_array, len(self.fptbinsDA)-1, cand_pt_bin_array)
        
        for index, row in df_sim.iterrows():
            his3.Fill(row['pt_jet'], row[self.var], row['pt_cand'])
        return his3
        
    return None
  
  def feeddown(self, obs_axis, jetR, beta ):
    print("===============Evaluating feed-down contribution")
    
    input_data = self.get_simulated_yields(self.simulationfile, 3, False, False, obs_axis)
    
    # get prompt and non prompt efficiency
    self.efficiency()
    
    label = "R%s_%s" % (str(jetR).replace('.', ''), str(beta).replace('.', ''))
    print("label is ", label)
    obs_bin_array = array('d', self.fobs_bins)
    jet_pt_bin_array = array('d', self.fptbinsJetFinalA)
    cand_pt_bin_array = array('d', self.fptbinsDA)
    # get nonprompt response
    if not os.path.exists(self.output_dir_fd):
        os.makedirs(self.output_dir_fd)
    roounfold_filename = os.path.join(self.output_dir_fd, "fRoounfold_R%s_%s.root" % (jetR, beta))
    roounfold_exists = os.path.exists(roounfold_filename)
    
    response = ROOT.TFile(self.responsefile, 'READ')
    name_ch = "hResponse_nonprompt_JetPt_ang_R%s_%sScaled" % (jetR, beta)
    thn_ch = response.Get(name_ch)
    print(thn_ch)
    # create Roounfold object
    name_roounfold_ch = '%s_Roounfold' % (name_ch)
    
    if roounfold_exists and not self.rebin_theory_response:
        fRoo = ROOT.TFile(roounfold_filename, 'READ')
        roounfold_response_ch = fRoo.Get(name_roounfold_ch)
        fRoo.Close()
    else:
        det_pt_bin_array = jet_pt_bin_array
        tru_pt_bin_array = det_pt_bin_array
        det_obs_bin_array = obs_bin_array
        tru_obs_bin_array = det_obs_bin_array
        obs_bins_gr = None; det_obs_bin_array_gr = None; tru_obs_bin_array_gr = None;
    
        n_dim = 4
    
        self.histutils.rebin_thn(
              roounfold_filename, thn_ch, '%s_Rebinned' % (name_ch), name_roounfold_ch, n_dim,
              len(det_pt_bin_array)-1, det_pt_bin_array, len(det_obs_bin_array)-1, det_obs_bin_array,
              len(tru_pt_bin_array)-1, tru_pt_bin_array, len(tru_obs_bin_array)-1, tru_obs_bin_array,
              label)
    
    f_resp = ROOT.TFile(roounfold_filename, 'READ')
    roounfold_response_ch = f_resp.Get(name_roounfold_ch)
    f_resp.Close()
    setattr(self, name_roounfold_ch, roounfold_response_ch)

    print("roo unfold object is ", name_roounfold_ch)
    
    
    input_data_angvsjetpt_wdmesonptcut_list = []
    hangvsjetpt_fd = ROOT.TH2F("hangvsjetpt_fd", "hangvsjetpt_fd",  len(self.fptbinsJetFinalA)-1, jet_pt_bin_array, len(self.fobs_bins)-1, obs_bin_array)
    for ibin2 in range(len(self.fptbinsJetFinalA)-1):
        suffix_jet = str(self.fptbinsJetFinalA[ibin2])+" < p_{T,jet} < "+str(self.fptbinsJetFinalA[ibin2+1])
        print("jet pt bins" + suffix_jet)
        hjetptvsangvsdmesonpt = input_data.Clone("hjetptvsangvsdmesonpt_"+suffix_jet)
        hjetptvsangvsdmesonpt.GetXaxis().SetRangeUser(self.fptbinsJetFinalA[ibin2],self.fptbinsJetFinalA[ibin2+1])
        hjetptvsangvsdmesonpt_injetbin = hjetptvsangvsdmesonpt.Clone()
        
        first_fit = 0
        
        input_data_effscaled = ROOT.TH2F()
        input_data_angvsjetpt_in_dmesonbins_list = []
        for ipt in range(self.fptbinsDN):
            suffix = str(self.fptbinsJetFinalA[ibin2])+" < p_{T,jet} < "+str(self.fptbinsJetFinalA[ibin2+1])+", "+\
                                 str(self.fptbinsDA[ipt])+" < #it{p}_{T}^{D^{0}} < "+str(self.fptbinsDA[ipt + 1])
            #minimum pt cut on Dmeson in a jet bin
            if self.fptbinsDA[ipt] < self.pt_dmeson_minpt_cut[ibin2]:
                input_data_angvsjetpt_in_dmesonbins_list.append(None)
                continue
            #maximum pt cut on Dmeson in a jet bin
            if self.fptbinsDA[ipt+1] > self.fptbinsJetFinalA[ibin2+1]:
                break
            
            # scale by ratio of non prompt/ prompt efficiency
            hjetptvsangvsdmesonpt_injetbin_clonned = hjetptvsangvsdmesonpt_injetbin.Clone("shapevsjetpt_"+suffix)
            hjetptvsangvsdmesonpt_injetbin_clonned.GetZaxis().SetRange(ipt+1, ipt + 1)
            input_data_angvsjetpt_in_dmesonbins_list.append(hjetptvsangvsdmesonpt_injetbin_clonned.Project3D("yxe"))
            
    
            for ibinshape in range(len(self.fobs_bins)-1):
                ratio_efficiency = (self.efficiency_nonprompt.GetBinContent(ipt+1))/(self.efficiency_prompt.GetBinContent(ipt+1))
                input_data_angvsjetpt_in_dmesonbins_list[ipt].SetBinContent(ibin2 + 1, ibinshape + 1, input_data_angvsjetpt_in_dmesonbins_list[ipt].GetBinContent( ibin2 + 1, ibinshape + 1)*ratio_efficiency)
                input_data_angvsjetpt_in_dmesonbins_list[ipt].SetBinError(ibin2 + 1, ibinshape + 1, input_data_angvsjetpt_in_dmesonbins_list[ipt].GetBinError( ibin2 + 1, ibinshape + 1)*ratio_efficiency)
                
            #sum over all Dmeson pt bins
            if first_fit == 0:
                input_data_effscaled = input_data_angvsjetpt_in_dmesonbins_list[ipt].Clone("input_data_effscaled_{}_{}_{}_{}.pdf".format(self.fptbinsJetFinalA[ibin2],self.fptbinsJetFinalA[ibin2+1],self.fptbinsDA[ipt],self.fptbinsDA[ipt + 1]))
                first_fit = 1
            else:
                input_data_effscaled.Add(input_data_angvsjetpt_in_dmesonbins_list[ipt])
           
        #add all jet pt bin shape
        input_data_angvsjetpt_wdmesonptcut_list.append(input_data_effscaled)
   
    #extract shape vs jetpt
    for ibin2 in range(len(self.fptbinsJetFinalA)-1):
        for ibiny in range(len(self.fobs_bins)-1):
            hangvsjetpt_fd.SetBinContent(ibin2+1,ibiny+1,input_data_angvsjetpt_wdmesonptcut_list[ibin2].GetBinContent(1,ibiny+1))
            hangvsjetpt_fd.SetBinError(ibin2+1,ibiny+1,input_data_angvsjetpt_wdmesonptcut_list[ibin2].GetBinError(1,ibiny+1))

    hangvsjetpt_fd_uncorr =  hangvsjetpt_fd.Clone("hangvsjetpt_raw")
   
    # apply gen. level kinematic efficiency
    hangvsjetpt_fd_part_level_corr = hangvsjetpt_fd.Clone("hangvsjetpt_fd_part_level_corr")
    hangvsjetpt_fd_part_level_corr.Multiply(self.particle_level_kinematic_eff)
    
    # scale with real luminosity and branching ratio
    scale_factor = float(self.scale_factor_pr_xsec)* float(self.p_nevents) * float(self.branching_ratio) / float(self.xsection_inel)
    hangvsjetpt_fd_part_level_corr.Scale(scale_factor)
    hangvsjetpt_fd_uncorr.Scale(scale_factor)
    
    # fold with the non-prompt response matrix
    folded = roounfold_response_ch.ApplyToTruth(hangvsjetpt_fd_part_level_corr)
    folded.SetName("foldedSpectrum" + str (obs_axis))
    
    # apply rec. level kinematic efficiency
    folded.Divide(self.gen_level_kinematic_eff)
    
    fout = ROOT.TFile("AnalysisResults_outputfd_{}.root".format(beta), 'recreate')
    folded.Write()
    hangvsjetpt_fd_uncorr.Write()
    fout.Close()
    
    return folded
    print("end of feed-down")
    
  #---------------------------------------------------------------
  # This function is called once for each subconfiguration
  #---------------------------------------------------------------
  def plot_single_result(self, jetR, obs_label, obs_setting, grooming_setting):
    #print('Plotting each individual result...')

    # Plot final result for each 1D substructure distribution
    self.plot_final_result(jetR, obs_label, obs_setting, grooming_setting)


  #----------------------------------------------------------------------
  def plot_final_result(self, jetR, obs_label, obs_setting, grooming_setting):
    print('Plot final results for {}: R = {}, {}'.format(self.observable, jetR, obs_label))
    
    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()
    # Loop through pt slices, and plot final result for each 1D theta_g distribution
    for i in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[i]
      max_pt_truth = self.pt_bins_reported[i+1]
      maxbin = self.obs_max_bins(obs_label)[i]
      
      self.plot_observable(jetR, obs_label, obs_setting, grooming_setting,
                           min_pt_truth, max_pt_truth, maxbin, plot_MC=True)
                           
                           
  #----------------------------------------------------------------------
  def plot_observable(self, jetR, obs_label, obs_setting, grooming_setting,
                      min_pt_truth, max_pt_truth, maxbin, plot_MC=False,
                      plot_theory=False, plot_theory_Fnp=False):
    
    self.set_logy = False
    
    name = 'cResult_R{}_{}_{}-{}'.format(jetR, obs_label, min_pt_truth, max_pt_truth)
    c = ROOT.TCanvas(name, name, 600, 450)
    c.Draw()
    
    c.cd()
    myPad = ROOT.TPad('myPad', 'The pad',0,0,1,1)
    myPad.SetLeftMargin(0.2)
    myPad.SetTopMargin(0.07)
    myPad.SetRightMargin(0.04)
    myPad.SetBottomMargin(0.13)
    if self.set_logy:
      myPad.SetLogy()
    myPad.Draw()
    myPad.cd()
    
    xtitle = getattr(self, 'xtitle')
    ytitle = getattr(self, 'ytitle')
    color = 1   # black for data
    
    # Get histograms
    name = 'hmain_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label,
                                              min_pt_truth, max_pt_truth)
    if grooming_setting:
      fraction_tagged = getattr(self, 'tagging_fraction_R{}_{}_{}-{}'.format(
        jetR, obs_label, min_pt_truth, max_pt_truth))
      #fraction_tagged = getattr(self, '{}_fraction_tagged'.format(name))
      # maxbin+1 in grooming case to account for extra tagging bin
    new_name = name+'_trunc'
    
    if grooming_setting and maxbin:
      h = self.truncate_hist(getattr(self, name), None, maxbin+1, name+'_trunc')
    else:
      h = self.truncate_hist(getattr(self, name), None, maxbin, new_name)
    h.SetMarkerSize(1.5)
    h.SetMarkerStyle(20)
    h.SetMarkerColor(color)
    h.SetLineStyle(1)
    h.SetLineWidth(2)
    h.SetLineColor(color)
    
    h_sys = getattr(self, 'hResult_{}_systotal_R{}_{}_{}-{}'.format(
      self.observable, jetR, obs_label, min_pt_truth, max_pt_truth))
    h_sys.SetLineColor(0)
    h_sys.SetFillColor(color)
    h_sys.SetFillColorAlpha(color, 0.3)
    h_sys.SetFillStyle(1001)
    h_sys.SetLineWidth(0)

    n_obs_bins_truth = self.n_bins_truth(obs_label)
    truth_bin_array = self.truth_bin_array(obs_label)
    if maxbin:
      truth_bin_array = truth_bin_array[0:maxbin+1]
      n_obs_bins_truth = len(truth_bin_array)-1
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_obs_bins_truth, truth_bin_array)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetXTitle(xtitle)
    myBlankHisto.GetXaxis().SetTitleOffset(1.02)
    myBlankHisto.GetXaxis().SetTitleSize(0.055)
    myBlankHisto.SetYTitle(ytitle)
    myBlankHisto.GetYaxis().SetTitleOffset(1.1)
    myBlankHisto.GetYaxis().SetTitleSize(0.055)
    ymin = 1e-4 if self.set_logy else 0
    myBlankHisto.SetMinimum(ymin)
    if not plot_theory or not show_parton_theory:
      maxval = max(2.3*h.GetBinContent(int(0.4*h.GetNbinsX())), 1.7*h.GetMaximum())
      myBlankHisto.SetMaximum(maxval)
      myBlankHisto.Draw("E")
      
    output_dir = getattr(self, 'output_dir_final_results')
    output_dir_single = output_dir + '/single_results'
    if not os.path.exists(output_dir_single):
      os.mkdir(output_dir_single)
    outputFilename = os.path.join(output_dir_single, name)
    c.SaveAs(outputFilename)
    c.Close()
    
    # Write result to ROOT file
    final_result_root_filename = os.path.join(output_dir, 'fFinalResults.root')
    fFinalResults = ROOT.TFile(final_result_root_filename, 'UPDATE')
    h.Write()
    h_sys.Write()
  #---------------------------------------------------------------
  # This function is called once after all subconfigurations have been looped over, for each R
  #---------------------------------------------------------------
  def plot_all_results(self, jetR):

    print('Plotting overlay of all results...')
    
    for i_config, overlay_list in enumerate(self.plot_overlay_list):
    
      if len(overlay_list) > 1:
      
        self.plot_final_result_overlay(i_config, jetR, overlay_list)


  #----------------------------------------------------------------------
  def plot_final_result_overlay(self, i_config, jetR, overlay_list):
    print('Plotting overlay of', overlay_list)

    # Plot overlay of different subconfigs, for fixed pt bin
    for i in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[i]
      max_pt_truth = self.pt_bins_reported[i+1]
      maxbins = [self.obs_max_bins(obs_label)[i] for obs_label in self.obs_labels]

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
  analysis.run_analysis()
