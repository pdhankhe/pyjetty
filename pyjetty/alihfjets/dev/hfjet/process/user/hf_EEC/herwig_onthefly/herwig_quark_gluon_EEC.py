
#!/usr/bin/env python
'''
Script for looking at the quark vs gluon dependence in Herwig of substructure observables
Author: Beatrice Liang-Gilman
'''

from __future__ import print_function

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjext
import ecorrel

import ROOT

import tqdm
import yaml
import copy
import argparse
import os
import array
import numpy as np
from array import array
import math

from pyjetty.mputils import *
from pyjetty.mputils.mputils import pinfo, pwarning

from pyjetty.alice_analysis.process.base import process_base

from enum import Enum
import fjtools

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)
# Automatically set Sumw2 when creating new histograms
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

################################################################
class EMesonDecayChannel(Enum):
	kAnyDecay            = 0
	kUnknownDecay        = 1 #BIT(0)
	kDecayD0toKpi        = 2 #BIT(1)
	kDecayDStartoKpipi   = 3 #BIT(2)

class Promptness(Enum):
	kUnknown = 0
	kPrompt = 1
	kNonPrompt = 2

################################################################
class PythiaQuarkGluon(process_base.ProcessBase):
    #---------------------------------------------------------------
	# Constructor
	#---------------------------------------------------------------
	def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, args=None, **kwargs):

		super(PythiaQuarkGluon, self).__init__(
			input_file, config_file, output_dir, debug_level, **kwargs)

		# Call base class initialization
		process_base.ProcessBase.initialize_config(self)

		# Read config file
		with open(self.config_file, 'r') as stream:
			config = yaml.safe_load(stream)

		if not os.path.exists(self.output_dir):
			os.makedirs(self.output_dir)

		self.jetR_list = config["jetR"]

		self.user_seed = args.user_seed
		self.nev = args.nev

		self.noMPI = (bool)(1-args.MPIon)
		self.noISR = (bool)(1-args.ISRon)


################################################################
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='herwig fastjet on the fly',
									 prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	# Could use --py-seed
	parser.add_argument('--user-seed', help='PYTHIA starting seed', default=1111, type=int)
	parser.add_argument('-o', '--output-dir', action='store', type=str, default='./',
						help='Output directory for generated ROOT file(s)')
	parser.add_argument('--tree-output-fname', default="AnalysisResults.root", type=str,
						help="Filename for the (unscaled) generated particle ROOT TTree")
	parser.add_argument('--MPIon', action='store', type=int, default=1,
						help="MPI on or off")
	parser.add_argument('--ISRon', action='store', type=int, default=1,
						help="ISR on or off")
	parser.add_argument('-cf', '--config_file', action='store', type=str, default='config/angularity.yaml',
						help="Path of config file for observable configurations")


	args = parser.parse_args()
	pinfo("The arguments to run are: ", args)

	# If invalid configFile is given, exit
	if not os.path.exists(args.config_file):
		print('File \"{0}\" does not exist! Exiting!'.format(args.config_file))
		sys.exit(0)

	# Use PYTHIA seed for event generation
	if args.user_seed < 0:
		args.user_seed = 1111

	# Have at least 1 event
	if args.nev < 1:
		args.nev = 1

	print("args for charmdecay", args.nocharmdecay)

	process = PythiaQuarkGluon(config_file=args.config_file, output_dir=args.output_dir, args=args)
	process.pythia_quark_gluon(args)