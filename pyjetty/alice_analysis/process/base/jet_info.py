#!/usr/bin/env python3

"""
    Class to store various jet info, used primarily for jet matching.
    
  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# Base class
from pyjetty.alice_analysis.process.base import common_base

################################################################
class JetInfo(common_base.CommonBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, **kwargs):
    super(JetInfo, self).__init__(**kwargs)

    # Store the matching candidates
    self.matching_candidates = []
    self.closest_jet = None
    self.closest_jet_deltaR = 1000.
    self.match = None

    # Store the associated truth info and particle charge (used primarily for applying pair efficeincy in the fast simulation)
    self.particle_truth = None
    self.charge= 1000.
    self.particle_pid = 0 # Store the particle PID, need for D0 identification
    self.particle_rap = 9999. #Store the particle rapidity, needed for D0 cuts
    self.particle_mid = 0 # Store the particle's mother's PID, need for D0 parent identification

  def clear_jet_info(self):
    self.matching_candidates.clear()
    self.closest_jet = None
    self.closest_jet_deltaR = 1000.
    self.match = None

  def clear_part_info(self):
    self.particle_truth = None
    self.charge= 1000.
    self.particle_pid = 0
    self.particle_rapi = 9999.
    self.particle_mid = 0

  def clear(self):
    self.matching_candidates.clear()
    self.closest_jet = None
    self.closest_jet_deltaR = 1000.
    self.match = None
    
    self.particle_truth = None
    self.charge= 1000.
    self.particle_pid = 0
    self.particle_rapi = 9999.
    self.particle_mid = 0