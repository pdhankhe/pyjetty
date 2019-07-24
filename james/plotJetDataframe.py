#!/usr/bin/env python3

"""
  Simple analysis script to read a pandas dataframe of track information
  and do jet-finding, and plot some basic histograms.
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import argparse
import math
import time

# Data analysis and plotting
import pandas as pd
import numpy as np
import ROOT

# Fastjet via python (from external library fjpydev)
import fastjet as fj
from recursivetools import pyrecursivetools as rt

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

# Set debug level (0 = no debug info, 1 = some debug info, 2 = all debug info)
debugLevel = 0

#---------------------------------------------------------------
def plotJetDataframe(inputFileTrack, outputDir, fileFormat):
  
  start_time = time.time()
  
  # Create output dir
  if not outputDir.endswith("/"):
    outputDir = outputDir + "/"
  if not os.path.exists(outputDir):
    os.makedirs(outputDir)

  # Load dataframe from .pkl file
  # dfTrack is a dataframe with one row per jet constituent: run_number, ev_id, ParticlePt, ParticleEta, ParticlePhi
  print('--- %s seconds ---' % (time.time() - start_time))
  print('Load track dataframe and unpickle...')
  dfTrack = pd.read_pickle(inputFileTrack)
  if debugLevel > 0:
    print(dfTrack.dtypes)
  if debugLevel > 1:
    print(dfTrack)

  # Transform the track dataframe into a SeriesGroupBy object of fastjet particles per event
  # This is currently the slowest part of the analysis, by factor ~10...
  print('--- %s seconds ---' % (time.time() - start_time))
  print('Transform the track dataframe into a series object of fastjet particles per event...')

  # (i) Group the track dataframe by event
  #     df_tracks is a DataFrameGroupBy object with one track dataframe per event
  df_tracks = dfTrack.groupby(['run_number','ev_id'])

  # (ii) Transform the DataFrameGroupBy object to a SeriesGroupBy of fastjet particles
  df_fjparticles = df_tracks.apply(get_fjparticles)
  if debugLevel > 0:
    print(df_fjparticles.dtypes)
    print(df_fjparticles)
  print('--- %s seconds ---' % (time.time() - start_time))

  # Print number of events
  nEvents = df_tracks.size().count()
  print('Number of events: {}'.format(nEvents))
  nTracks = len(dfTrack.index)
  print('Number of tracks: {}'.format(nTracks))

  # Initialize histogram dictionary
  hDict = initializeHistograms()
  
  # Find jets and fill histograms
  print('Find jets...')
  analyzeEvents(df_fjparticles, hDict, outputDir, fileFormat)

  # Plot histograms
  print('Plot histograms...')
  plotHistograms(hDict, outputDir, fileFormat)

  print('--- %s seconds ---' % (time.time() - start_time))

#---------------------------------------------------------------
def get_fjparticles(df_tracks):
  
  fj_particles = []
  for index, row in df_tracks.iterrows(): # Not sure how to do this without iterating...
    pt = row['ParticlePt']
    eta = row['ParticleEta']
    phi = row['ParticlePhi']
    px = pt * math.cos(phi)
    py = pt * math.sin(phi)
    pz = pt * math.sinh(eta)
    e = pt
    fj_particles.append(fj.PseudoJet(px, py, pz, e))
  
  return fj_particles

#---------------------------------------------------------------
def initializeHistograms():
  
  hDict = {}
  
  hJetPt = ROOT.TH1F('hJetPt', 'hJetPt', 200, 0, 200)
  hJetPt.GetXaxis().SetTitle('p_{T,jet}')
  hJetPt.GetYaxis().SetTitle('dN/dp_{T}')
  hDict['hJetPt'] = hJetPt

  hAJ = ROOT.TH2F('hAJ', 'hAJ', 100, 0, 1., 100, 0, 4.)
  hAJ.GetXaxis().SetTitle('A_{J}')
  hAJ.GetYaxis().SetTitle('#Delta #phi')
  hDict['hAJ'] = hAJ

  hZg = ROOT.TH1F('hZg', 'hZg', 100, 0, 1.)
  hZg.GetXaxis().SetTitle('z_{g}')
  hZg.GetYaxis().SetTitle('dN/dz_{g}')
  hDict['hZg'] = hZg

  hRg = ROOT.TH1F('hRg', 'hRg', 100, 0, 1.)
  hRg.GetXaxis().SetTitle('R_{g}')
  hRg.GetYaxis().SetTitle('dN/dR_{g}')
  hDict['hRg'] = hRg

  return hDict

#---------------------------------------------------------------
def analyzeEvents(df_fjparticles, hDict, outputDir, fileFormat):
  
  fj.ClusterSequence.print_banner()
  print()
  
  # Set jet definition and a jet selector
  jetR = 0.4
  jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
  jet_selector = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9 - jetR)
  print('jet definition is:', jet_def)
  print('jet selector is:', jet_selector,'\n')
  
  # Define SoftDrop settings
  beta = 0
  zcut = 0.1
  sd = rt.SoftDrop(beta, zcut, jetR)
  print('SoftDrop groomer is: {}'.format(sd.description()));

  # Use list comprehension to do jet-finding and fill histograms
  result = [analyzeJets(fj_particles, jet_def, jet_selector, sd, hDict) for fj_particles in df_fjparticles]

#---------------------------------------------------------------
def analyzeJets(fj_particles, jet_def, jet_selector, sd, hDict):
  
  # Do jet finding
  cs = fj.ClusterSequence(fj_particles, jet_def)
  jets = fj.sorted_by_pt(cs.inclusive_jets())
  jets_accepted = jet_selector(jets)

  fillJetHistograms(hDict, jets_accepted)

  # Loop through jets and perform SoftDrop grooming
  jets_sd = []
  for jet in jets_accepted:
    jets_sd.append(sd.result(jet))

  fillSoftDropHistograms(hDict, jets_sd)

#---------------------------------------------------------------
def fillJetHistograms(hDict, jets_accepted):

  # Loop through jets, and fill histograms
  for jet in jets_accepted:
    if debugLevel > 1:
      print('jet:')
      print(jet)
      
    hDict['hJetPt'].Fill(jet.pt())

  # Find di-jets and fill histograms
  if len(jets_accepted) > 1:

    pT1 = jets_accepted[0].pt()
    pT2 = jets_accepted[1].pt()
    phi1 = jets_accepted[0].phi()
    phi2 = jets_accepted[1].phi()
  
    AJ = (pT1 - pT2) / (pT1 + pT2)
    deltaPhi = abs(phi1 - phi2)
    
    hDict['hAJ'].Fill(AJ, deltaPhi)

#---------------------------------------------------------------
def fillSoftDropHistograms(hDict, jets_sd):
  
  for jet in jets_sd:
    
    sd_info = rt.get_SD_jet_info(jet)
    zg = sd_info.z
    Rg = sd_info.dR
    
    hDict['hZg'].Fill(zg)
    hDict['hRg'].Fill(Rg)

#---------------------------------------------------------------
def plotHistograms(hDict, outputDir, fileFormat):

  outputFilename = "{}hJetPt{}".format(outputDir, fileFormat)
  plotHist(hDict['hJetPt'], outputFilename, 'width', setLogy=True)

  outputFilename = "{}hAJ{}".format(outputDir, fileFormat)
  plotHist(hDict['hAJ'], outputFilename, 'colz')

  outputFilename = "{}hZg{}".format(outputDir, fileFormat)
  plotHist(hDict['hZg'], outputFilename)

  outputFilename = "{}hRg{}".format(outputDir, fileFormat)
  plotHist(hDict['hRg'], outputFilename)

#---------------------------------------------------------------
def plotHist(h, outputFilename, drawOptions = "", setLogy = False, setLogz = False):
  
  h.SetLineColor(1)
  h.SetLineWidth(1)
  h.SetLineStyle(1)
  
  c = ROOT.TCanvas("c","c: hist",600,450)
  c.cd()
  ROOT.gPad.SetLeftMargin(0.15)
  if setLogy:
    c.SetLogy()
  if setLogz:
    c.SetLogz()
  ROOT.gPad.SetLeftMargin(0.15)

  h.Draw(drawOptions)
  c.SaveAs(outputFilename)
  c.Close()

#----------------------------------------------------------------------
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description='Plot analysis histograms')
  parser.add_argument('-f', '--inputFileTrack', action='store',
                      type=str, metavar='inputFileTrack',
                      default='AnalysisResultsTrackMergedIndex0.pkl',
                      help='Path of pickled AnalysisResults file for track tree')
  parser.add_argument('-o', '--outputDir', action='store',
                      type=str, metavar='outputDir',
                      default='./myTestFigures',
                      help='Output directory for QA plots to be written to')
  parser.add_argument('-i', '--imageFormat', action='store',
                      type=str, metavar='imageFormat',
                      default='.pdf',
                      help='Image format to save plots in, e.g. \'.pdf\' or \'.png\'')
  
  # Parse the arguments
  args = parser.parse_args()
  
  print('Configuring...')
  print('inputFileTrack: \'{0}\''.format(args.inputFileTrack))
  print('ouputDir: \'{0}\"'.format(args.outputDir))
  print('imageFormat: \'{0}\''.format(args.imageFormat))
  print('----------------------------------------------------------------')
  
  # If invalid inputFile is given, exit
  if not os.path.exists(args.inputFileTrack):
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFileTrack))
    sys.exit(0)

  plotJetDataframe(inputFileTrack = args.inputFileTrack, outputDir = args.outputDir, fileFormat = args.imageFormat)