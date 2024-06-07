#!/usr/bin/env python

from __future__ import print_function

import os
import argparse
import sys

# import pyhepmc
import pyhepmc_ng

# import hepmc2antuple_base

# jit improves execution time by 18% - tested with jetty pythia8 events
# BUT produces invalid root file
# from numba import jit
# @jit


# Open the file in read mode
file_ahadic = "/rstorage/ploskon/eec_sherpa_charm/charm_jetpt15/list.txt"
file_lund = "/rstorage/ploskon/eec_sherpa_charm/charm_jetpt15_lund/list.txt"

with open(file_lund, 'r') as file:
  # Loop through each line in the file
  for ifile, input_file in enumerate(file):
    if ifile != 28:
      continue
    # Process each line as needed
    print("intput file: ", ifile, input_file.strip())  # strip() removes any leading/trailing whitespace
    input_file = input_file.strip()


    # input_file = "/rstorage/ploskon/eec_sherpa_charm/charm_jetpt15/999983/sherpa_LHC_jets_15.0.hepmc"
    # input_file = "/software/users/ploskon/eecmc/sherpa2x/charm_jetpt10_lund/sherpa_LHC_jets_10.0.hepmc"
    hepmc = 3
      
    if hepmc == 3:
      # input_hepmc = pyhepmc.io.ReaderAscii(input_file)
      input_hepmc = pyhepmc_ng.ReaderAscii(input_file)
    if hepmc == 2:
      # input_hepmc = pyhepmc.io.ReaderAsciiHepMC2(input_file)
      input_hepmc = pyhepmc_ng.ReaderAsciiHepMC2(input_file)

    if input_hepmc.failed():
      print ("[error] unable to read from {}".format(input_file))
      sys.exit(1)

    # event_hepmc = pyhepmc.GenEvent()
    event_hepmc = pyhepmc_ng.GenEvent()
    print("done!")

    evid = 1
    while not input_hepmc.failed():
      ev = input_hepmc.read_event(event_hepmc)
      # print("input_hepmc", evid, input_hepmc.failed())
      if input_hepmc.failed():
        break
      # self.fill_event(event_hepmc)
      # self.increment_event()
      if evid > 0 and 250000 < evid:
        break
      evid+=1

    print("Total number of events in hepmc file", ifile, "is ", evid)
      