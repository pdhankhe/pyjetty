#! /usr/bin/env python

""" theory_comp.py
Post processing code
Preeti Dhankher, 2020 (pdhankher@berkeley.edu)
"""

import sys
import os
import argparse
from array import *
import numpy as np
import math   # for exp()
import ROOT
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

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)


################################################################
# Helper functions
################################################################

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
# Reflection function
#---------------------------------------------------------------
def SetReflection(fitter, fLeftFitRange, fRightFitRange, RS, iBin):

  fileRefl = ROOT.TFile.Open(fReflFilename.Data())

  if not fileRefl:
    print("File "+str(fReflFilename)" (reflection templates) cannot be opened! check your file path!")
    return False

  fHistnameRefl = "histRflFittedDoubleGaus_ptBin"
  fHistnameSign = "histSgn_"
  print(Form("%s%d",fHistnameRefl.Data(),iBin))
  print(Form("%s%d",fHistnameSign.Data(),iBin))
  histRefl = fileRefl.Get(Form("%s%d",fHistnameRefl.Data(),iBin-3))
  histSign = fileRefl.Get(Form("%s%d",fHistnameSign.Data(),iBin-3))

  if not histRefl or not histSign:
    print("Error in loading the template/signal histrograms! Exiting...")
    return False

  fitter->SetTemplateReflections(histRefl,"template",fLeftFitRange,fRightFitRange);
  Double_t RoverS = histRefl->Integral(histRefl->FindBin(fLeftFitRange),histRefl->FindBin(fRightFitRange))/histSign->Integral(histSign->FindBin(fLeftFitRange),histSign->FindBin(fRightFitRange));
  if(isRefSys) RoverS*=refScale;
  printf("R/S ratio in fit range for bin %d = %1.3f\n",iBin,RoverS);
  fitter->SetFixReflOverS(RoverS);

  RS = RoverS
  return True

#---------------------------------------------------------------
# Given two data points, find linear fit and y-value for x
#---------------------------------------------------------------
def raw_jet_shape(outdir, prod, alpha):
  hmean = ROOT.TH1F('mean', 'mean', fptbinsJetMeasN, fptbinsJetMeasA)
  hsigma = ROOT.TH1F("hsigma", "hsigma", fptbinsJetMeasN, fptbinsJetMeasA)
  hrelErr = ROOT.TH1F("hrelErr", "hrelErr", fptbinsJetMeasN, fptbinsJetMeasA)
  hsign = ROOT.TH1F("hsign", "hsign", fptbinsJetMeasN, fptbinsJetMeasA)
  hsb = ROOT.TH1F("hsb", "hsb", fptbinsJetMeasN, fptbinsJetMeasA)
  hSignal = ROOT.TH1F("hSignal", "hSignal", fptbinsJetMeasN, fptbinsJetMeasA)
  hSignal.Sumw2()
  hReflRS = ROOT.TH1F("hReflRS", "hReflRS", fptbinsJetMeasN, fptbinsJetMeasA)

  hInvMassptD.GetYaxis().SetTitle("Counts per 5 MeV/c^{2}")
  hInvMassptD.GetXaxis().SetTitle(Form("m(%s)(GeV/c^{2})", "K#pi"))
  hInvMassptD.GetXaxis().SetTitleSize(0.06)
  hInvMassptD.GetXaxis().SetTitleOffset(0.7)
  hInvMassptD.GetYaxis().SetTitleOffset(0.7)
  hInvMassptD.SetTitle("")

  xnx = 3
  xny = 4

  if fptbinsDN > 4 and fptbinsDN < 7:
    xnx = 2
    xny = 3
  if fptbinsDN > 6 ND fptbinsDN < 10:
    xnx = 3
    xny = 3
  if fptbinsDN > 9 ND fptbinsDN < 13:
    xnx = 3
    xny = 4
  else:
    xnx = 4
    xny = 4

  c2 = ROOT.TCanvas("c2", "c2", 1200, 1200)
  c2.Divide(xnx, xny)

  c2jet = ROOT.TCanvas("c2jet", "c2jet", 1200, 1200)
  c2jet.Divide(xnx, xny)

  c2jetcorr = ROOT.TCanvas("c2jetcorr", "c2jetcorr", 1200, 1200)
  c2jetcorr.Divide(xnx, xny)

  firstPtBin = 0

  if fptbinsDA[0] == 2:
    firstPtBin = 3
  if fptbinsDA[0] == 3:
    firstPtBin = 4
  if fptbinsDA[0] == 4:
    firstPtBin = 5
  if fptbinsDA[0] == 5:
    firstPtBin = 6
  if fptbinsDA[0] == 6:
    firstPtBin = 7
  if fptbinsDA[0] == 7:
    firstPtBin = 8

  if not firstPtBin:
    print("==== Wrong first value of the D pT (should be 2,3 or 4) ===")


  RS = 0

  for (int i=0; i < 6; i++):
    hh = hInvMassptD.ProjectionX(Form("hh_%d", i),hInvMassptD.GetYaxis().FindBin(jetmin),
                               hInvMassptD.GetYaxis().FindBin(jetmax) - 1,hInvMassptD.GetZaxis().FindBin(fptbinsDA[i]),
                               hInvMassptD.GetZaxis().FindBin(fptbinsDA[i + 1]) - 1)

    hh.Rebin(fRebinMass)
    hh.GetXaxis().SetRangeUser(minf, maxf)
    hh.SetTitle(Form(" %.1lf < p_{T,jet} < %.1lf,  %.1lf < p_{T}^{%s} < %.1lf", fptbinsJetFinalA[0], fptbinsJetFinalA[1],
                  fptbinsDA[i], fDmesonS.Data(), fptbinsDA[i + 1]))
    hh.GetYaxis().SetTitle("Counts per 5 MeV/c^{2}")
    hh.GetYaxis().SetTitleSize(0.06)
    hh.GetYaxis().SetTitleOffset(0.7)

    hmassfit = hh.Clone()
    hmassfit.SetMaximum(hmassfit.GetMaximum()* 1.3)
    hmin = TMath::Max(minf, hmassfit.GetBinLowEdge(2))
    hmax = TMath::Min(maxf, hmassfit.GetBinLowEdge(hmassfit.GetNbinsX()));


    AliHFInvMassFitter * fitterp = new AliHFInvMassFitter((TH1F *)hmassfit, hmin, hmax, fbkgtype, 0)
    fitterp.SetInitialGaussianMean(fDmass)
    fitterp.SetInitialGaussianSigma(fDsigma)

    if fUseRefl:
      SetReflection(fitterp, hmin, hmax, RS, i+firstPtBin)

    fitterp.MassFitter(kFALSE)

    h = fitterp.GetHistoClone()
    massfit[i] = fitterp.GetMassFunc()
    massfit[i].SetRange(hmin, hmax)
    massfit[i].SetLineColor(4)
    fullfit[i] = h.GetFunction("funcmass")

    if fullfit[i]:
      fullfit[i].SetName(Form("fullfit_%d", i))

    hmass[i] = h.Clone(Form("hmass_%d", i))
    hmass[i].SetName(Form("hmass_%d", i))
    hmass[i].GetYaxis().SetTitle("Entries")
    hmass[i].GetXaxis().SetTitle("Invariant Mass (GeV/c^{2})")
    bkgfit[i] = fitterp.GetBackgroundRecalcFunc()
    bkgfit[i].SetRange(hmin, hmax)
    bkgfit[i].SetLineColor(2)
    bkgfit[i].SetName(Form("bkgFit_%d", i))\

    # bkg + reflection function

    if fUseRefl:
      bkgRfit[i] = fitterp.GetBkgPlusReflFunc()
      bkgRfit[i].SetName(Form("bkgFitWRef_%d", i))
      bkgRfit[i].SetRange(hmin, hmax)
      bkgRfit[i].SetLineColor(15)
      hReflRS.SetBinContent(i+1, RS)
      hReflRS.SetBinError(i+1, 0)

    TVirtualPad * pad = (TVirtualPad *)
    c2.GetPad(i + 1)
    fitterp.DrawHere(pad, 3, 0)
    Dsigma = 0
    Dmean = 0
    DmeanUnc = 0
    DsigmaUnc = 0

    if not fullfit[i]:
      print("======= Fit failed for bin: ")

    Dsigma = fitterp.GetSigma()
    DsigmaUnc = fitterp.GetSigmaUncertainty()
    Dmean = fitterp.GetMean()
    DmeanUnc = fitterp.GetMeanUncertainty()

    signal_c_min = Dmean - fsigmaSignal * Dsigma
    signal_c_max = Dmean + fsigmaSignal * Dsigma
    signal_l_min = Dmean + fsigmaBkg[0] * Dsigma
    signal_l_max = Dmean + fsigmaBkg[1] * Dsigma
    signal_u_min = Dmean + fsigmaBkg[2] * Dsigma
    signal_u_max = Dmean + fsigmaBkg[3] * Dsigma

    if signal_l_min < hmin:
      signal_l_min = hmin

    if signal_u_max > hmax:
      signal_u_max = hmax

    binwidth = hmass[i].GetXaxis().GetBinWidth(1) * 0.5
    # signal
    binmin = hmass[i].GetXaxis().FindBin(signal_c_min)
    binmax = hmass[i].GetXaxis().FindBin(signal_c_max)
    min = hmass[i].GetXaxis().GetBinCenter(binmin) - binwidth
    max = hmass[i].GetXaxis().GetBinCenter(binmax - 1) + binwidth

    # side - bands
    binmin_sb1 = hmass[i].GetXaxis().FindBin(signal_l_min)
    binmax_sb1 = hmass[i].GetXaxis().FindBin(signal_l_max)
    min_sb1 = hmass[i].GetXaxis().GetBinCenter(binmin_sb1) - binwidth
    max_sb1 = hmass[i].GetXaxis().GetBinCenter(binmax_sb1 - 1) + binwidth
    # side - bands
    binmin_sb2 = hmass[i]->GetXaxis()->FindBin(signal_u_min)
    binmax_sb2 = hmass[i]->GetXaxis()->FindBin(signal_u_max)
    min_sb2 = hmass[i]->GetXaxis()->GetBinCenter(binmin_sb2) - binwidth
    max_sb2 = hmass[i]->GetXaxis()->GetBinCenter(binmax_sb2 - 1) + binwidth

    s = 0, serr = 0, srelerr = 0, bkg = 0, bkgerr = 0, bkgref = 0, ref = 0
    ref1 = 0, ref2 = 0

    fitterp.Signal(fsigmaSignal, s, serr)
    fitterp.Background(min, max, bkg, bkgerr)
    sb1 = bkgfit[i].Integral(min_sb1, max_sb1) / (Double_t)hmass[i].GetBinWidth(1)
    sb2 = bkgfit[i].Integral(min_sb2, max_sb2) / (Double_t)hmass[i].GetBinWidth(1)

    if fUseRefl:
      bkgref = bkgRfit[i].Integral(min, max) / (Double_t)hmass[i].GetBinWidth(1)
      ref = bkgref - bkg
      ref1 = (bkgRfit[i].Integral(min_sb1, max_sb1) / hmass[i].GetBinWidth(1)) - sb1
      ref2 = (bkgRfit[i].Integral(min_sb2, max_sb2) / hmass[i].GetBinWidth(1)) - sb2

    signf = 0, signferr = 0, sob = 0, soberr = 0
    fitterp.Significance(fsigmaSignal, signf, signferr)

    if s:
      srelerr = serr / s
    if bkg:
      sob = s / bkg
    else:
      sob = s

    if bkg and bkgerr:
      soberr = TMath::Sqrt((serr / bkg) * (serr / bkg) + (s / bkg / bkg * bkgerr) * (s / bkg / bkg * bkgerr))
    else:
      soberr = serr

    #---------------- fitting results
    hmean.SetBinContent(i + 1, Dmean * 1000)
    hmean.SetBinError(i + 1, DmeanUnc * 1000)

    hsigma.SetBinContent(i + 1, Dsigma * 1000)
    hsigma.SetBinError(i + 1, DsigmaUnc * 1000)

    hrelErr.SetBinContent(i + 1, srelerr)
    hsign.SetBinContent(i + 1, signf)
    hsign.SetBinError(i + 1, signferr)
    hsb.SetBinContent(i + 1, sob)
    hsb.SetBinError(i + 1, soberr)
    hSignal.SetBinContent(i + 1, s)
    hSignal.SetBinError(i + 1, serr)

    #---------------- side - band drawing
    hmass_l[i] = hmass[i].Clone("hmass_l")
    hmass_l[i].GetXaxis().SetRangeUser(signal_l_min, signal_l_max)
    hmass_l[i].SetName(Form("hmass_l_%d", i))
    hmass_u[i] = hmass[i].Clone("hmass_u")
    hmass_u[i].GetXaxis().SetRangeUser(signal_u_min, signal_u_max)
    hmass_u[i].SetName(Form("hmass_u_%d", i))
    hmass_c[i] = hmass[i].Clone("hmass_c")
    hmass_c[i].GetXaxis().SetRangeUser(signal_c_min, signal_c_max)
    hmass_c[i].SetName(Form("hmass_c_%d", i))
    hmass_l[i].SetFillColor(kBlue + 2)
    hmass_u[i].SetFillColor(kBlue + 2)
    hmass_c[i].SetFillColor(kRed + 2)
    hmass_l[i].SetFillStyle(3004)
    hmass_u[i].SetFillStyle(3004)
    hmass_c[i].SetFillStyle(3005)
    hmass_l[i].Draw("hsame")
    hmass_u[i].Draw("hsame")
    hmass_c[i].Draw("hsame")

    #---------------- Extracting signal shape ---------
    hjetpt[i] = hInvMassptD.ProjectionY(Form("hjetpt_%d", i), hInvMassptD.GetXaxis().FindBin(signal_c_min),
                                      hInvMassptD.GetXaxis().FindBin(signal_c_max) - 1,
                                      hInvMassptD.GetZaxis().FindBin(fptbinsDA[i]),
                                      hInvMassptD.GetZaxis().FindBin(fptbinsDA[i + 1]) - 1)

    hjetpt[i].SetTitle(Form("%.1lf < p_{T}^{%s} < %.1lf", fptbinsDA[i], fDmesonS.Data(), fptbinsDA[i + 1]))
    hjetpt[i].SetMarkerColor(kBlue + 2)
    hjetpt[i].SetLineColor(kBlue + 2)

    # ---------------- Extracting left SB shape
    hjetpt_sb1 = hInvMassptD.ProjectionY(Form("hjetpt_sb1%d", i), hInvMassptD.GetXaxis().FindBin(signal_l_min),
                                          hInvMassptD.GetXaxis().FindBin(signal_l_max) - 1,
                                          hInvMassptD.GetZaxis().FindBin(fptbinsDA[i]),
                                          hInvMassptD.GetZaxis().FindBin(fptbinsDA[i + 1]) - 1)
    hjetpt_sb1.SetTitle(Form("%.1lf < p_{T}^{%s} < %.1lf", fptbinsDA[i], fDmesonS.Data(), fptbinsDA[i + 1]))
    hjetpt_sb1.SetMarkerColor(kGray + 2)
    hjetpt_sb1.SetLineColor(kGray + 2)\

    # ---------------- Extracting right SB shape
    hjetpt_sb2 = hInvMassptD.ProjectionY(Form("hjetpt_sb2%d", i), hInvMassptD.GetXaxis().FindBin(signal_u_min),
                                         hInvMassptD.GetXaxis().FindBin(signal_u_max) - 1,
                                         hInvMassptD.GetZaxis().FindBin(fptbinsDA[i]),
                                         hInvMassptD.GetZaxis().FindBin(fptbinsDA[i + 1]) - 1)

    hjetpt_sb2.SetTitle(Form("%.1lf < p_{T}^{%s} < %.1lf", fptbinsDA[i], fDmesonS.Data(), fptbinsDA[i + 1]));
    hjetpt_sb2.SetMarkerColor(kBlue + 2)
    hjetpt_sb2.SetLineColor(kBlue + 2)

    hjetpt_sb[i] = hjetpt_sb1.Clone(Form("hjetpt_sb_%d", i))
    hjetpt_sb[i].Add(hjetpt_sb2)

    # Scaling - - Reflection signal OFF
    scalingB = bkg / (sb1 + sb2)

    # Scaling - - Reflection signal ON
    if fUseRefl:
      scalingS = s / (s + ref - ((ref1 + ref2) * bkg) / (sb1 + sb2))

    hjetpt_sb[i].Scale(scalingB)

    # ------- subtract background from signal jet
    hjetptsub[i] = hjetpt[i].Clone(Form("hjetptsub_%d", i))
    hjetptsub[i].Add(hjetpt_sb[i], -1)
    if fUseRefl:
      hjetptsub[i].Scale(scalingS)

    if fsigmaSignal == 2:
      print("Scaling for limited 2σ width of the signal region")
      hjetptsub[i].Scale(1. / 0.9545)

    hjetptsub[i].SetMarkerColor(kGreen + 3)
    hjetptsub[i].SetLineColor(kGreen + 3)
    c2jet.cd(i + 1)
    gPad.SetLogy()

    hjetpt[i].Draw("ep")
    hjetpt[i].GetXaxis()->SetTitle("#lambda_{#alpha=3}")
    hjetpt[i].GetYaxis()->SetTitle("Counts")
    hjetpt[i].GetXaxis()->SetTitleSize(0.06)
    hjetpt[i].GetXaxis()->SetTitleOffset(0.7)
    hjetpt[i].GetYaxis()->SetTitleSize(0.06)
    hjetpt[i].GetYaxis()->SetTitleOffset(0.7)
    hjetpt_sb[i].Draw("epsame")
    hjetptsub[i].Draw("epsame")

    if not hrawjetptspectrum:
      hrawjetptspectrum = hjetptsub[i].Clone("hrawjetptspectrum")
    else:
      hrawjetptspectrum.Add(hjetptsub[i])

    # ------- correct for D * efficiency
    hjetptcorr[i] = (TH1F * )hjetptsub[i]->Clone(Form("hjetptcorr_%d", i))
    print("efficiency values"+str(efficiency[i]))
    if bEff:
      hjetptcorr[i].Scale(1. / efficiency[i])

    hjetptcorr[i].SetMarkerColor(kBlue+3)
    hjetptcorr[i].SetLineColor(kBlue+3)

    c2jetcorr.cd(i+1)
    hjetptcorr[i].Draw("ep")


    # ----- add corrected jet pt spectrum distributions in each D pt bin into a one distribution
    if not hjetptspectrum:
      hjetptspectrum = hjetptcorr[i].Clone("hjetptspectrum")
    else:
      hjetptspectrum.Add(hjetptcorr[i])

  c2->cd(fptbinsDN+1)
  c2jet->cd(fptbinsDN+1)


  canvas_angdist = ROOT.TCanvas("canvas_angdist", "canvas_angdist", 1200, 1200)
  fbinsJetAngularityN = 16
  fbinsJetAngularityA = array('d',(0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
                                   0.6, 0.65, 0.7, 0.75, 0.8))
  hAngdist= hjetptspectrum.Clone("hAngdist")
  hRebinAngdist= hAngdist.Rebin(15, "hRebinAngdist", fbinsJetAngularityA)
  hRebinAngdist.Scale(1 / hRebinAngdist.Integral(), "width")
  #ProcessHisto(hRebinAngdist, marker_size, colors[0], markers[0])
  hRebinAngdist.GetXaxis().SetTitle("#lambda^{k=1}_{#alpha=1}")
  hRebinAngdist.GetYaxis().SetTitle("1/#it{N}_{jet} d#it{N}/d#lambda")
  hRebinAngdist.Draw("ep")

  return True


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
    
    c = ROOT.TCanvas("", "", 600, 650)
    
    sparse = ROOT.TFile(self.datafile,'READ').Get("hsparsejet")
    if not sparse:
        print("reflection file missing")
    hInvMassptD =sparse.Projection(2)
    hInvMassptD.Draw("E")
    print("working here")
    #print(self)
    
    print(self.efffile)
   
    eff_prompt = ROOT.TFile("$HOME/Analysis/analysisHF/ang/reco_efficiency/D0JetRecoEff.root",'READ').Get("NonPromptEfficiency")
    if not eff_prompt:
        print("efficiency file missing")
    #eff_prompt.Draw("E")
    
    reflection = ROOT.TFile(self.reflection_template[0],'READ').Get("histSgn_0")
    if not reflection:
            print("reflection file missing")
  
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
    
    
    self.datafile = config['main_data']
    print(self.datafile)
    self.efffile = config['main_efffile']
    print(self.efffile)
    self.reflection_template = config['reflection_template']
    print(self.reflection_template[0])
    

    # Whether or not to use the previous preliminary result in final plots
    self.use_prev_prelim = config['use_prev_prelim']

    self.histutils = ROOT.RUtil.HistUtils()


    


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
