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
    self.rebin_theory_response = config['rebin_theory_response']

    self.datafile = config['main_rawdatafile']
    print(self.datafile)
    self.simulationfile = config['main_feeddown_simulation']
    print(self.simulationfile)
    self.responsefile = config['main_response']
    self.output_dir_fd = config['output_dir_fd']
    self.efffile = config['main_efffile']
    print(self.efffile)
    self.reflection_template = config['reflection_template']
    self.isRefSys = False
    print(self.reflection_template[0])
    

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
    obs_axis = 3
    
    obs_bin_array = array('d', self.fobs_bins)
    jet_pt_bin_array = array('d', self.fptbinsJetFinalA)
    
    hangvsjetpt = ROOT.TH2F("hangvsjetpt", "hangvsjetpt",  len(self.fptbinsJetFinalA)-1, jet_pt_bin_array, len(self.fobs_bins)-1, obs_bin_array)
    hangvsjetpt.Sumw2()
    for beta in self.beta_list:
        print("calculating shape for the observable ", obs_axis)
        #loop over jet pt
        for ibin2 in range(len(self.fptbinsJetFinalA)-1):
            #Clone sparse info
            suffix_jet = str(self.fptbinsJetFinalA[ibin2])+" < p_{T,jet} < "+str(self.fptbinsJetFinalA[ibin2+1])
            first_fit = 0
            sparsehadron = sparse.Clone("sparsehadron")
            #apply jet selection range
            sparsehadron.GetAxis(0).SetRangeUser(self.jetrange[0], self.jetrange[1])
            #select jet pt bins
            print("Analysis in jet pt bins = "+str(self.fptbinsJetFinalA[ibin2])+"-"+str(self.fptbinsJetFinalA[ibin2+1]))
            sparsehadron.GetAxis(0).SetRangeUser(self.fptbinsJetFinalA[ibin2], self.fptbinsJetFinalA[ibin2+1])
            #x-D-mass, y-shape, z- Dmeson pt
            h_invmass_shape_Dmesonpt = sparsehadron.Projection(2, obs_axis, 1)
            h_invmass_shape_Dmesonpt.GetXaxis().SetTitle("Counts per 5 MeV/c^{2}")
            h_invmass_shape_Dmesonpt.GetXaxis().SetTitle("m(K#pi)(GeV/c^{2})")
            h_invmass_shape_Dmesonpt.GetXaxis().SetTitleSize(0.06)
            h_invmass_shape_Dmesonpt.GetXaxis().SetTitleOffset(0.7)
            h_invmass_shape_Dmesonpt.GetYaxis().SetTitleOffset(0.7)
            h_invmass_shape_Dmesonpt.SetTitle("")
    
            
            RS = 0
            #loop over cand. pt
            input_data_invmass_list = []
            hmass_l_list = []
            hmass_u_list = []
            hmass_c_list = []
            bkg_func_list = []
            cinvmass = ROOT.TCanvas("cinvmass ", "cinvmass")
            cinvmass.Divide(4,3)
        
            for ipt in range(self.fptbinsDN):
              print("min dmeson pt cut"+str(self.pt_dmeson_minpt_cut[ibin2]))
              if self.fptbinsDA[ipt] < self.pt_dmeson_minpt_cut[ibin2]:
                
                input_data_invmass_list.append(None)
                hmass_l_list.append(None)
                hmass_u_list.append(None)
                hmass_c_list.append(None)
                bkg_func_list.append(None)
                continue
                
              text_for_invmass = str(self.fptbinsJetFinalA[ibin2])+" < p_{T,jet} < "+str(self.fptbinsJetFinalA[ibin2+1])+", "+\
                                 str(self.fptbinsDA[ipt])+" < #it{p}_{T}^{D^{0}} < "+str(self.fptbinsDA[ipt + 1])
              suffix = text_for_invmass
              
              print("Analysis in Dmeson pt bins = "+str(self.fptbinsDA[ipt])+"-"+str(self.fptbinsDA[ipt+1]))
              hinvmass = h_invmass_shape_Dmesonpt.ProjectionX("hinvmass"+suffix, h_invmass_shape_Dmesonpt.GetYaxis().FindBin(0.001),
                                           h_invmass_shape_Dmesonpt.GetYaxis().FindBin(0.8) - 1,
                                           h_invmass_shape_Dmesonpt.GetZaxis().FindBin(self.fptbinsDA[ipt]),
                                           h_invmass_shape_Dmesonpt.GetZaxis().FindBin(self.fptbinsDA[ipt + 1]) - 1,"e")
              
              
                                 
              
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
              fitterp.SetInitialGaussianMean(self.Dmass)
              fitterp.SetInitialGaussianSigma(self.Dsigma)
              
              if self.IsUseRefl:
                self.SetReflection(fitterp, hmin, hmax, ipt)
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

              ref_func = fitterp.GetReflFunc() if self.IsUseRefl else None

              sig_func = fitterp.GetSignalFunc()
              sig_func.SetRange(hmin, hmax)

              fullfit = histo_invmassfit.GetFunction("funcmass")

              hmass = histo_invmassfit.Clone()

              if self.IsUseRefl:
                bkgRfit = fitterp.GetBkgPlusReflFunc()
                bkgRfit.SetRange(hmin, hmax)
                bkgRfit.SetLineColor(15)
                hReflRS.SetBinContent(ipt + 1, self.RS)
                hReflRS.SetBinError(ipt + 1, 0)



              if not fullfit:
                print("======= Fit failed for bin: ")
                

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
              process_histo(input_data_invmass_list[ipt],0.8, 1 ,20)
              input_data_invmass_list[ipt].Draw("same,ep")
              bkg_func_list[ipt].Draw("same")
              hmass_l_list[ipt].Draw("hsame")
              hmass_u_list[ipt].Draw("hsame")
              hmass_c_list[ipt].Draw("hsame")
             
             
              
              # ---------------- Extracting signal shape
              hang_sig = h_invmass_shape_Dmesonpt.ProjectionY("signal"+ suffix, h_invmass_shape_Dmesonpt.GetXaxis().FindBin(signal_c_min),
              h_invmass_shape_Dmesonpt.GetXaxis().FindBin(signal_c_max) - 1,
              h_invmass_shape_Dmesonpt.GetZaxis().FindBin(self.fptbinsDA[ipt]),
              h_invmass_shape_Dmesonpt.GetZaxis().FindBin(self.fptbinsDA[ipt + 1]) - 1,"e")


              # ---------------- Extracting left SB shape
              hang_sig_sb1 = h_invmass_shape_Dmesonpt.ProjectionY("signalleft" + suffix, h_invmass_shape_Dmesonpt.GetXaxis().FindBin(signal_l_min),
              h_invmass_shape_Dmesonpt.GetXaxis().FindBin(signal_l_max) - 1,
              h_invmass_shape_Dmesonpt.GetZaxis().FindBin(self.fptbinsDA[ipt]),
              h_invmass_shape_Dmesonpt.GetZaxis().FindBin(self.fptbinsDA[ipt + 1]) - 1,"e")

              # ---------------- Extracting right SB shape
              hang_sig_sb2 = h_invmass_shape_Dmesonpt.ProjectionY("signalright" + suffix , h_invmass_shape_Dmesonpt.GetXaxis().FindBin(signal_u_min),
              h_invmass_shape_Dmesonpt.GetXaxis().FindBin(signal_u_max) -1,
              h_invmass_shape_Dmesonpt.GetZaxis().FindBin(self.fptbinsDA[ipt]),
              h_invmass_shape_Dmesonpt.GetZaxis().FindBin(self.fptbinsDA[ipt + 1]) - 1,"e")

              hang_bkg = hang_sig_sb1.Clone("hang_bkg"+ suffix)
              hang_bkg.Add(hang_sig_sb2)


              # Scaling - - Reflection signal OFF
              scalingB = bkg / (sb1 + sb2)

              print("scaling B is"+str(scalingB))

              # Scaling - - Reflection signal ON
              if self.IsUseRefl:
                scalingS = signal / (signal + ref - ((ref1 + ref2) * bkg) / (sb1 + sb2))
                print("scaling S is"+str(scalingS))

              hang_bkg_scaled = hang_bkg.Clone("hzbkg_scaled" + suffix)
              hang_bkg_scaled.Scale(scalingB)

              # ------- subtract background from signal jet
              hang_sig_bkgsub = hang_sig.Clone("hang_sig_bkgsub" + suffix)
              hang_sig_bkgsub.Add(hang_bkg_scaled, -1)

              if self.IsUseRefl:
                hang_sig_bkgsub.Scale(scalingS)

              if self.Signalsigma == 2:
                print("Scaling for limited 2σ width of the signal region")
                hang_sig_bkgsub.Scale(1. / 0.9545)
                
              if first_fit == 0:
                hrawjetptspectrum = hang_sig_bkgsub.Clone("hrawjetptspectrum"+suffix)
                first_fit = 1
              else:
                hrawjetptspectrum.Add(hang_sig_bkgsub)


              # ------- correct shape for prompt efficiency
              hang_sig_bkgsub_effcorr = hang_sig_bkgsub.Clone("hang_sig_bkgsub_effcorr" + suffix)
              self.efficiency()
              ptbin = (self.fptbinsDA[ipt]+self.fptbinsDA[ipt+1]) / 2.
              print("efficiency values" + str(self.efficiency_prompt.GetBinContent(self.efficiency_prompt.GetXaxis().FindBin(ptbin))))
              hang_sig_bkgsub_effcorr.Scale(1. / (self.efficiency_prompt.GetBinContent(self.efficiency_prompt.GetXaxis().FindBin(ptbin))))

              #input_data_shape_list.append(hang_sig_bkgsub_effcorr)
              
              if first_fit-1 == 0:
                hang_corr = hang_sig_bkgsub_effcorr.Clone("hjetptspectrum"+suffix)
                first_fit = 2
              else:
                hang_corr.Add(hang_sig_bkgsub_effcorr)
                
                   
            cinvmass.SaveAs("invmass_{}_{}.pdf".format(self.fptbinsJetFinalA[ibin2],self.fptbinsJetFinalA[ibin2+1]))
                
            canvas_angdist = ROOT.TCanvas("canvas_angdist_{}_{}".format(self.fptbinsJetFinalA[ibin2], self.fptbinsJetFinalA[ibin2+1]),"canvas_angdist_{}_{}".format(self.fptbinsJetFinalA[ibin2], self.fptbinsJetFinalA[ibin2+1]),1200,1200)
            hAngdist = hang_corr.Clone("hAngdist" + suffix_jet +"_"+str(beta))
            hRebinAngdist = hAngdist.Rebin(len(self.fobs_bins)-1, "hRebinAngdist", obs_bin_array)
            print("the number of jets ", str(hRebinAngdist.Integral()))
            hRebinAngdist.Scale(1/hRebinAngdist.Integral(),"width")
            
            hRebinAngdist.GetXaxis().SetTitle("#lambda^{k=1}_{#alpha=1}")
            hRebinAngdist.GetYaxis().SetTitle("1/#it{N}_{jet} d#it{N}/d#lambda")
            name = "h_ang_JetPt%d_R%s_%s" % (self.fptbinsJetFinalA[ibin2],jetR, beta)
            hRebinAngdist.Draw("ep")
            canvas_angdist.SaveAs(name+".pdf")
            
            
            # fill the 2D histogram shape vs jet pt
            for angbins in range(len(self.fobs_bins)-1):
                    if (hRebinAngdist.GetBinContent(angbins + 1) < 0):
                        hangvsjetpt.SetBinContent( ibin2 + 1, angbins + 1, 0)
                        hangvsjetpt.SetBinError( ibin2 + 1, angbins + 1, 0)
                    
                    else:
                        hangvsjetpt.SetBinContent( ibin2 + 1, angbins + 1,hRebinAngdist.GetBinContent(angbins + 1))
                        hangvsjetpt.SetBinError( ibin2 + 1, angbins + 1, hRebinAngdist.GetBinError(angbins + 1))
         
         
        #hfeeddown = self.feeddown(obs_axis, jetR, beta)
        hangvsjetpt_feedsub = hangvsjetpt.Clone("h_ang_JetPt_R%s_%s" % (jetR, beta))
        #hangvsjetpt_feedsub.Add(hfeeddown,-1)
        #check the bins with negative yield after subtraction and set them to zero
        hangvsjetpt_feedsub_projy = hangvsjetpt_feedsub.ProjectionY("y")
        
        
        for ibin2 in range(len(self.fptbinsJetFinalA)-1):
            for angbins in range(len(self.fobs_bins)-1):
                if (hangvsjetpt_feedsub_projy.GetBinContent(angbins + 1) < 0):
                    hangvsjetpt_feedsub.SetBinContent( ibin2 + 1, angbins + 1, 0)
                    hangvsjetpt_feedsub.SetBinError( ibin2 + 1, angbins + 1, 0)
                    
    
        name = "h_ang_JetPt_R%s_%s" % (jetR, beta)
        setattr(self, name, hangvsjetpt_feedsub)
        
        
        fout = ROOT.TFile("AnalysisResults_output_{}.root".format(beta), 'recreate')
        hangvsjetpt_feedsub.Write()
        fout.Close()
        
        #move to next shape observable
        #obs_axis += 1

    print("length of histogram list"+str(len(input_data_invmass_list)))
    return True

  def efficiency(self):
    file_eff = ROOT.TFile(self.efffile, 'READ')
    print(self.efffile)
    self.efficiency_prompt = file_eff.Get("PromptEfficiency")
    self.efficiency_prompt.SetDirectory(0)
    self.efficiency_nonprompt = file_eff.Get("NonPromptEfficiency")
    self.efficiency_nonprompt.SetDirectory(0)
    return True
    
  def SetReflection(self, fitter, fLeftFitRange, fRightFitRange, iBin):
    file_ref = ROOT.TFile(self.reflection_template[0], 'READ')
    if not file_ref:
        print("File "+str(reflection_template)+" (reflection templates) cannot be opened! check your file path!")
        return False
    
    fHistnameRefl = "histRflFittedDoubleGaus_ptBin"
    fHistnameSign = "histSgn_"
    print("{}{}".format(fHistnameRefl,iBin-3))
    print("{}{}".format(fHistnameSign,iBin-3))
    histRefl = file_ref.Get("{}{}".format(fHistnameRefl,iBin-3))
    histSign = file_ref.Get("{}{}".format(fHistnameSign,iBin-3))
    
    if not histRefl or not histSign:
        print("Error in loading the template/signal histrograms! Exiting...")
        return False
        
    fitter.SetTemplateReflections(histRefl,"template",fLeftFitRange,fRightFitRange)
  
    
    histo_ref = ROOT.TH1F()
    histRefl.Copy(histo_ref)
    
    fitter.SetHistogramReflectionFit(histo_ref)
        
    RoverS=histRefl.Integral(histRefl.FindBin(fLeftFitRange),histRefl.FindBin(fRightFitRange))/histSign.Integral(histSign.FindBin(fLeftFitRange),histSign.FindBin(fRightFitRange))
    
    if self.isRefSys:
        RoverS*=refScale
        
    print("R/S ratio in fit range for bin {} = {}".format(iBin,RoverS))
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
    
    input_data_angvsjetpt_list = []
    input_data_scaled = ROOT.TH2F()
    
    # scale by ratio of non prompt/ prompt efficiency
    for ipt in range(self.fptbinsDN):
            input_data.GetZaxis().SetRange(ipt+1, ipt + 2)
            input_data_angvsjetpt_list.append(input_data.Project3D("yx"))
            for ibin2 in range(len(self.fptbinsJetFinalA)-1):
                for ibinshape in range(len(self.fobs_bins)-1):
                    ratio_efficiency = (self.efficiency_nonprompt.GetBinContent(ipt+1))/(self.efficiency_prompt.GetBinContent(ipt+1))
                    input_data_angvsjetpt_list[ipt].SetBinContent(  ibin2 + 1, ibinshape + 1, input_data_angvsjetpt_list[ipt].GetBinContent ( ibin2 + 1, ibinshape + 1)*ratio_efficiency)
                    
            if ipt == 0:
                input_data_scaled = input_data_angvsjetpt_list[ipt].Clone("input_data_scaled")
            else:
                input_data_scaled.Add(input_data_angvsjetpt_list[ipt])
        
    
    # apply gen. level kinematic efficiency
    #input_data_scaled.Multiply(1)
    # scale with real luminosity and branching ratio
    scale_factor = float(self.scale_factor_pr_xsec)* float(self.p_nevents) * float(self.branching_ratio) / float(self.xsection_inel)
    input_data_scaled.Scale(scale_factor)
    # fold with the non-prompt response matrix
    
    
    folded = roounfold_response_ch.ApplyToTruth(input_data_scaled)
    
    #folded = folding(input_data_scaled, response_matrix, output_template)
    folded.SetName("foldedSpectrum" + str (obs_axis))
    # apply rec. level kinematic efficiency
    
    fout = ROOT.TFile("AnalysisResults_outputfd_{}.root".format(beta), 'recreate')
    input_data_scaled.Write()
    folded.Write()
    roounfold_response_ch.Write()
    fout.Close()
    
    return folded
    print("end of feed-down")
    
    
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
  #analysis.run_analysis()
