"""
Script containing all helper functions related to plotting with ROOT

Script also contains the "class Errors", used for systematic uncertainties (to
replace AliHFSystErr from AliPhysics).
"""
# pylint: disable=too-many-lines
from array import array
import math
import numpy as np

from ROOT import TH1F, TH2F, TH2, TFile, TH1, TH3F, TGraphAsymmErrors
from ROOT import TPad, TCanvas, TLegend, kBlack, kGreen, kRed, kBlue, kWhite
from ROOT import gStyle, gROOT, TMatrixD

def makefill2dhist(df_, titlehist, arrayx, arrayy, nvar1, nvar2):
    """
    Create a TH2F histogram and fill it with two variables from a dataframe.
    """
    histo = build2dhisto(titlehist, arrayx, arrayy)
    df_rd = df_[[nvar1, nvar2]]
    arr2 = df_rd.to_numpy()
    fill_hist(histo, arr2)
    return histo
