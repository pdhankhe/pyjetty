
#!/usr/bin/env python
'''
Script for looking at the quark vs gluon dependence of substructure observables
Author: Beatrice Liang-Gilman, with most of the code from Ezra Lesser (elesser@berkeley.edu)
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

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext

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

		# self implemented variables to study
		self.charmdecaysOFF = (bool)(args.nocharmdecay) #charmdecaysOFF=True(F) when charmdecayon=1(0)
		pinfo("charm decay input value", args.nocharmdecay)
		self.weighted = (bool)(args.weightON) #weightON=True(F) means turn weights on(off)
		self.leading_parton_pt_cut = args.leadingptcut
		self.replaceKPpairs = (bool)(args.replaceKP) #replaceKP=True(F) means turn k/pi pairs are('nt) replaced
		self.Dstar = (bool)(args.DstarON) #Dstar=True means look at D* EEC, should be run with self.replaceKPpairs=True
		self.initscat = args.chinitscat #1=hard->ccbar, 2=gg->ccbar, 3=D0->Kpi channel, 4=hard->bbar w/ D0->Kpi
		self.D0wDstar = (bool)(args.D0withDstarON) #D0wDstar=True means looking at D-tagged jets including D0 from D*
		# self.softpion_action = args.softpion #1 = remove soft pion from D*, 2 = only pair soft pion with charged particles, 3 = only pair soft pion with D0, 4 = pair soft pion w everything
		# self.phimeson = (bool)(args.runphi) #1=don't let phi meson decay and look at its EEC


		# PDG ID values for quarks and gluons
		self.quark_pdg_ids = [1, 2, 3, 4, 5, 6, 7, 8, -1, -2, -3, -4, -5, -6, -7, -8]
		self.down_pdg_ids = [1, -1]
		self.up_pdg_ids = [2, -2]
		self.strange_pdg_ids = [3, -3]
		self.charm_pdg_ids = [4, -4]
		self.gluon_pdg_ids = [9, 21] 
		self.beauty_pdg_ids = [5, -5]

		# hadron level - ALICE tracking restriction
		self.max_eta_hadron = 0.9

		self.min_leading_track_pT = config["min_leading_track_pT"] if "min_leading_track_pT" in config else None

		self.pt_bins = array.array('d', list(range(5, 100, 5)) + list(range(100, 210, 10)))
		self.obs_bins_jet_axis = np.array([0, 0.01, 0.03, 0.05, 0.12, 0.25]) #should we cut 0.25 bin?

		# jet axis specific - TODO: FIX
		self.sd_zcut = config["sd_zcut"]
		self.sd_beta = config["sd_beta"]

		# self.obs_bins_STD__D = np.array([0, 0.01, 0.03, 0.05, 0.12, 0.25])
		# self.obs_bins_WTA__D = np.array([0, 0.01, 0.03, 0.05, 0.12, 0.25])
		# self.obs_bins_SD__D = np.array([0, 0.01, 0.03, 0.05, 0.12, 0.25])
		# self.obs_bins_STD__WTA = np.array([0, 0.01, 0.03, 0.05, 0.12, 0.25])
		# self.obs_bins_WTA__SD = np.array([0, 0.01, 0.03, 0.05, 0.12, 0.25])
		# self.obs_bins_STD__SD = np.array([0, 0.01, 0.03, 0.05, 0.12, 0.25])


		self.observable_list = config['process_observables']
		self.obs_settings = {}
		self.obs_grooming_settings = {}
		self.obs_names_xtitle = {}
		self.jetaxes_observables = {}
		self.obs_bins = {}

		
		for observable in self.observable_list:

			# jet axis specific
			self.jetaxes_observables[observable] = config["observable_list"]

			obs_config_dict = config[observable]
			print("obs_config_dict", obs_config_dict)
			obs_config_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]
			print("obs_config_list", obs_config_list)

			self.obs_settings[observable] = self.utils.obs_settings(observable, obs_config_dict, obs_config_list)
			pinfo("self.obs_settings[observable]", self.obs_settings[observable]) # all the different axes
			self.obs_grooming_settings[observable] = self.utils.grooming_settings(obs_config_dict)
			pinfo("self.obs_grooming_settings[observable]", self.obs_grooming_settings[observable])
			self.obs_bins[observable] = self.utils.binning_settings(self.obs_settings[observable], obs_config_dict)
			pinfo("self.obs_bins[observable]", self.obs_bins[observable])

			self.obs_names_xtitle[observable] = obs_config_dict["common_settings"]["xtitle"]

			

	#---------------------------------------------------------------
	# Main processing function
	#---------------------------------------------------------------
	def pythia_quark_gluon(self, args):

		# Create ROOT TTree file for storing raw PYTHIA particle information
		outf_path = os.path.join(self.output_dir, args.tree_output_fname)
		outf = ROOT.TFile(outf_path, 'recreate')
		outf.cd()

		# Initialize response histograms 
		self.initialize_hist()

		pinfo('user seed for pythia', self.user_seed) #TODO: what does this do?? it doesn't work...
#        print('user seed for pythia', self.user_seed)
		mycfg = ['Random:setSeed=on', 'Random:seed={}'.format(self.user_seed)]
		mycfg.append('HadronLevel:all=off')
		pinfo("charmdecays value", self.charmdecaysOFF)
		if (self.charmdecaysOFF == True and self.replaceKPpairs == False):
			pinfo("charm decays turning OFF")
			mycfg.append('411:mayDecay = no')
			mycfg.append('421:mayDecay = no')
			mycfg.append('10411:mayDecay = no')
			mycfg.append('10421:mayDecay = no')
			mycfg.append('413:mayDecay = no')
			mycfg.append('423:mayDecay = no')
			mycfg.append('10413:mayDecay = no')
			mycfg.append('10423:mayDecay = no')
			mycfg.append('20413:mayDecay = no')
			mycfg.append('20423:mayDecay = no')
			mycfg.append('415:mayDecay = no')
			mycfg.append('425:mayDecay = no')
			mycfg.append('431:mayDecay = no')
			mycfg.append('10431:mayDecay = no')
			mycfg.append('433:mayDecay = no')
			mycfg.append('10433:mayDecay = no')
			mycfg.append('20433:mayDecay = no')
			mycfg.append('435:mayDecay = no')

			mycfg.append('4122:mayDecay = no')
			mycfg.append('4222:mayDecay = no')
			mycfg.append('4212:mayDecay = no')
			mycfg.append('4112:mayDecay = no')
			mycfg.append('4224:mayDecay = no')
			mycfg.append('4214:mayDecay = no')
			mycfg.append('4114:mayDecay = no')
			mycfg.append('4232:mayDecay = no')
			mycfg.append('4132:mayDecay = no')
			mycfg.append('4322:mayDecay = no')
			mycfg.append('4312:mayDecay = no')
			mycfg.append('4324:mayDecay = no')
			mycfg.append('4314:mayDecay = no')
			mycfg.append('4332:mayDecay = no')
			mycfg.append('4334:mayDecay = no')
			mycfg.append('4412:mayDecay = no')
			mycfg.append('4422:mayDecay = no')
			mycfg.append('4414:mayDecay = no')
			mycfg.append('4424:mayDecay = no')
			mycfg.append('4432:mayDecay = no')
			mycfg.append('4434:mayDecay = no')
			mycfg.append('4444:mayDecay = no')

		if (self.initscat == 2): #if (self.gg2ccbar):
			mycfg.append('HardQCD:all = off')
			mycfg.append('HardQCD:gg2ccbar = on')

		if (self.initscat == 1): #if (self.hardccbar):
			mycfg.append('HardQCD:all = off')
			mycfg.append('HardQCD:hardccbar = on')

		if (self.initscat == 3): # just D0->Kpi
			mycfg.append('HardQCD:all = off')
			mycfg.append('HardQCD:hardccbar = on')

			mycfg.append('421:onMode = off')
			mycfg.append('421:onIfMatch = 321 211')

		if (self.initscat == 4): # hard->bbar with D0 -> (only) Kpi
			mycfg.append('HardQCD:all = off')
			mycfg.append('HardQCD:hardbbbar = on')

			mycfg.append('421:onMode = off')
			mycfg.append('421:onIfMatch = 321 211')

		if (self.replaceKPpairs):
			if (not (self.Dstar or self.D0wDstar)):
				pinfo("turning D*'s OFF")
				mycfg.append('10411:mayDecay = no')
				mycfg.append('10421:mayDecay = no')
				mycfg.append('413:mayDecay = no')
				mycfg.append('423:mayDecay = no')
				mycfg.append('10413:mayDecay = no')
				mycfg.append('10423:mayDecay = no')
				mycfg.append('20413:mayDecay = no')
				mycfg.append('20423:mayDecay = no')
				mycfg.append('415:mayDecay = no')
				mycfg.append('425:mayDecay = no')
				mycfg.append('431:mayDecay = no')
				mycfg.append('10431:mayDecay = no')
				mycfg.append('433:mayDecay = no')
				mycfg.append('10433:mayDecay = no')
				mycfg.append('20433:mayDecay = no')
				mycfg.append('435:mayDecay = no')

		# print the banner first
		fj.ClusterSequence.print_banner()
		print()

		# -------------------------------
		# Setting MPIs and ISRs
		print('Will run no MPI:',self.noMPI)
		print('Will run no ISR:',self.noISR)
		setattr(args, "py_noMPI", self.noMPI)
		setattr(args, "py_noISR", self.noISR)
		# -------------------------------

		pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
		# print("----------------- PARTICLE DATA INFO HERE -----------------")
		# pythia.particleData.listAll()
		# print("----------------- PARTICLE DATA INFO END -----------------")

		self.init_jet_tools()
		self.calculate_events(pythia)
		pythia.stat()
		print()

		self.scale_print_final_info(pythia)

		outf.Write()
		outf.Close()

		self.save_output_objects()

	#---------------------------------------------------------------
	# Initialize histograms
	#---------------------------------------------------------------
	def initialize_hist(self):

		self.hNevents = ROOT.TH1I("hNevents", 'Number accepted events (unscaled)', 2, -0.5, 1.5)
		self.hD0Nevents = ROOT.TH1I("hD0Nevents", "Total Number of D0 events (unscaled)", 2, -0.5, 1.5)
		self.hD0KpiNevents = ROOT.TH1I("hD0KpiNevents", "Number of D0->Kpi events (unscaled)", 2, -0.5, 1.5)
		self.hD0KpiNjets = ROOT.TH1I("hD0KpiNjets", "Number of D0->Kpi jets (unscaled)", 2, -0.5, 1.5) 
		self.hDstarNjets = ROOT.TH1I("hDstarNjets", "Number of D* jets (unscaled)", 2, -0.5, 1.5)
		# self.hsoftpionpT = ROOT.TH1D("hsoftpionpT", "pT of soft pion from D*", 50, 0, 50)
		
		for jetR in self.jetR_list:

			# Store a list of all the histograms just so that we can rescale them later
			hist_list_name = "hist_list_R%s" % str(jetR).replace('.', '')
			setattr(self, hist_list_name, [])

			# R_label = str(jetR).replace('.', '') + 'Scaled'

			for observable in self.observable_list:
				# Should only be one: observable == "jet_axis"
				pinfo("observable", observable)
				if observable != "jet_axis":
					raise ValueError("Observable %s is not implemented in this script" % observable)
				
				# Use more finely binned pT bins for TH2s than for the RMs
				pt_bins = array.array('d', list(range(0, 201, 1)))
				rapi_bins = np.linspace(-5,5,201)
				
				print("jetaxes obs", self.jetaxes_observables[observable])
				for i, axisname in enumerate(self.jetaxes_observables[observable]):
					# should have all the jet axes

					print("axisname", axisname)
					obs_name_xtitle = self.obs_names_xtitle[observable] #xaxis title
					# pinfo("obs_name_xtitle", obs_name_xtitle)

					#TODO: FIX OBS BINS
					# dashindex = axisname.find('-')
					# obs_adjustedname = axisname[:dashindex]+'__'+axisname[dashindex+1:]
					# obs_bins = getattr(self, "obs_bins_" + obs_adjustedname)
					pinfo("obs_bins", self.obs_bins[observable])
					obs_bins = self.obs_bins[observable][i].get(axisname)
					obs_bins = np.array(obs_bins)
					pinfo("obs_bins", type(obs_bins), obs_bins)


					dim = 4
					nbins  = [len(pt_bins)-1, len(pt_bins)-1, len(rapi_bins)-1, len(obs_bins)-1]
					min_li = [pt_bins[0],     pt_bins[0],      rapi_bins[0],      obs_bins[0]]
					max_li = [pt_bins[-1],    pt_bins[-1],     rapi_bins[-1],     obs_bins[-1]]

					nbins = (nbins)
					xmin = (min_li)
					xmax = (max_li)
					
					nbins_array = array.array('i', nbins)
					xmin_array = array.array('d', xmin)
					xmax_array = array.array('d', xmax)
					
					# Loop over subobservable (alpha value)
	#                for i in range(len(self.obs_settings[observable])):

					obs_setting = self.obs_settings[observable][i] #[observable] #['SD-D']
					grooming_setting = self.obs_grooming_settings[observable][i] if None != self.obs_grooming_settings[observable][i] else {} #{'sd': [0.2, 0]}
					obs_label = self.utils.obs_label(obs_setting, grooming_setting) #form of STD-WTA OR SD-D_SD_zcut02_B0
					pinfo("all the settings", obs_setting, grooming_setting, obs_label)

					self.fsparsepartonJetvalue = array.array( 'd', ( 0, 0, 0 ,0 ))
					
					partontypeslist = ["charm", "light", "gluon", "inclusive"] #got rid of quark
					if (self.initscat == 4):
						partontypeslist.append("beauty")
						
					for parton_type in partontypeslist:

						title = [ '#it{p}_{T}^{ch jet}', '#it{p}_{T}^{D^{0}}', 'y', obs_name_xtitle ]

						name = ('h_%s_JetPt_%s_R%s_%s' % (observable, parton_type, jetR, obs_label)) if \
							len(obs_label) else ('h_%s_JetPt_%s_R%s' % (observable, parton_type, jetR))
						hsparse = ROOT.THnSparseD(name,"%s-init_hsparsejet; #it{p}_{T,%s}^{ch jet}; #it{p}_{T}^{D^{0}}; y;%s^{%s}" %(parton_type[0], parton_type[0] + "-init", obs_name_xtitle, parton_type[0] + "-init"), dim,  nbins_array, xmin_array, xmax_array)
						hsparse.Sumw2()
						for i in range(0,dim):
							hsparse.GetAxis(i).SetTitle(title[i])
							if i == 0 or i == 1:
								hsparse.SetBinEdges(i, pt_bins)
							if i == 2:
								hsparse.SetBinEdges(i, rapi_bins)
							if i == 3:
								hsparse.SetBinEdges(i, obs_bins)
						setattr(self, name, hsparse)
						getattr(self, hist_list_name).append(hsparse)


						# h = ROOT.TH2F(name, name, len(pt_bins)-1, pt_bins, len(obs_bins)-1, obs_bins)
						# h.GetXaxis().SetTitle('#it{p}_{T,%s}^{ch jet}' % (parton_type[0] + "-init"))
						# h.GetYaxis().SetTitle(obs_name_xtitle + '^{%s}' % (parton_type[0] + "-init"))
						# h.GetYaxis().SetTitle("R_{L}" + '^{%s}' % (parton_type[0] + "-init"))
						# h.Sumw2()
						# setattr(self, name, h)
						# getattr(self, hist_list_name).append(h)

	#---------------------------------------------------------------
	# Initiate jet defs, selectors, and sd (if required)
	#---------------------------------------------------------------
	def init_jet_tools(self):

		for jetR in self.jetR_list:
			jetR_str = str(jetR).replace('.', '')

			# set up our jet definition and a jet selector
			jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
			setattr(self, "jet_def_R%s" % jetR_str, jet_def)
			print(jet_def)

		pwarning('max eta for particles after hadronization set to', self.max_eta_hadron)
#        print('max eta for particles after hadronization set to', self.max_eta_hadron)
		parts_selector_h = fj.SelectorAbsEtaMax(self.max_eta_hadron)
		track_selector_ch = fj.SelectorPtMin(0.15) & parts_selector_h #ALICE parameters
		setattr(self, "track_selector_ch", track_selector_ch)
		
		for jetR in self.jetR_list:
			jetR_str = str(jetR).replace('.', '')

			jet_selector = fj.SelectorPtMin(5.0) & fj.SelectorAbsEtaMax(self.max_eta_hadron - jetR)
			#jet_selector = fj.SelectorPtMin(0.) & fj.SelectorAbsEtaMax(self.max_eta_hadron - jetR)
			setattr(self, "jet_selector_R%s" % jetR_str, jet_selector)

			count1 = 0  # Number of partonic parents which match to >1 ch-jets
			setattr(self, "count1_R%s" % jetR_str, count1)
			count2 = 0  # Number of partonic parents which match to zero ch-jets
			setattr(self, "count2_R%s" % jetR_str, count2)

	#---------------------------------------------------------------
	# Calculate events and pass information on to jet finding
	#---------------------------------------------------------------
	def calculate_events(self, pythia):

		iev = 0  # Event loop count

		while iev < self.nev:
			if not pythia.next():
				continue
			
			if (iev%10000 == 0):
				print("Event", iev)
			# print("Event", iev)

			self.parents = []
			self.event = pythia.event
			#print(self.event)
			fs_parton_5 = fj.PseudoJet(pythia.event[5].px(), pythia.event[5].py(), pythia.event[5].pz(), pythia.event[5].e())
			fs_parton_6 = fj.PseudoJet(pythia.event[6].px(), pythia.event[6].py(), pythia.event[6].pz(), pythia.event[6].e())
			self.parents = [fs_parton_5, fs_parton_6] # parent partons in dijet

			# Save PDG code of the parent partons
			self.parent_ids = [pythia.event[5].id(), pythia.event[6].id()]

			# parton level
			#parts_pythia_p = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)

			hstatus = pythia.forceHadronLevel()
			if not hstatus:
				continue

			# full-hadron level
			#parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)

			# charged-hadron level
			if ( self.replaceKPpairs == False ):
				parts_pythia_hch = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged], 0, True)
			else: #replace D0->Kpi
				parts_pythia_hch = pythiafjext.vectorize_select_replaceD0(pythia, [pythiafjext.kFinal, pythiafjext.kCharged], 0, True, True)
			
			
			# Count D0 events
			D0found = False
			D0Kpidecayfound = False
			self.DstarKpipidecayfound = False
			for particle in self.event: #for event in pythia.event:
				if particle.id() == 421 or particle.id() == -421: #D0
					D0found = True
					if self.checkDecayChannel(particle, self.event) == EMesonDecayChannel.kDecayD0toKpi:
						D0Kpidecayfound = True

					if self.checkDecayChannel(particle, self.event) == EMesonDecayChannel.kDecayDStartoKpipi:
						self.DstarKpipidecayfound = True #can't fill histogram here because it will fill at particle level
			
			#if D0->Kpi found, count the events; if not, check that length of charged final state hadrons vector is 0
			if (D0Kpidecayfound):
				self.hD0KpiNevents.Fill(0)
			if (D0found):
				self.hD0Nevents.Fill(0)
			
			
			
			# Some "accepted" events don't survive hadronization step -- keep track here
			self.hNevents.Fill(0)
			self.find_jets_fill_histograms(parts_pythia_hch, iev, D0Kpidecayfound)

			iev += 1

	#---------------------------------------------------------------
	# Find primordial parent
	#---------------------------------------------------------------
	def primordial_parent(self,p):
		parent1 = parent2 = -10
		while p > 6:
			parent1 = self.event[p].mother1()
			parent2 = self.event[p].mother2()
			if parent1 != parent2:
				p = max(parent1,parent2)
			else:
				p = parent1
		return p
	
	# trk_thrd default set 0, meaning all tracks would pass
	def checkIfPartInJetConst(self, jet_const_arr, pythia_particle_index, trk_thrd=0):
		in_jet = False
		for c in jet_const_arr:
			# print("jet const user index", c.user_index(), pythiafjext.getPythia8Particle(c).name())
			if (c.user_index() == pythia_particle_index and c.pt() >= trk_thrd):
				in_jet = True
				# print("ifpartinjet", c.user_index(), pythia_particle_index)
				break
		return in_jet

	#---------------------------------------------------------------
	# Find jets, do matching between levels, and fill histograms
	#---------------------------------------------------------------
	def find_jets_fill_histograms(self, parts_pythia_hch, iev, D0Kpidecayfound):
		# Loop over jet radii
		for jetR in self.jetR_list:

			jetR_str = str(jetR).replace('.', '')
			jet_selector = getattr(self, "jet_selector_R%s" % jetR_str)
			jet_def = getattr(self, "jet_def_R%s" % jetR_str)
			track_selector_ch = getattr(self, "track_selector_ch")

			count1 = getattr(self, "count1_R%s" % jetR_str)
			count2 = getattr(self, "count2_R%s" % jetR_str)

			# Get the jets at different levels
			#jets_p  = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_p  ))) # parton level
			#jets_h  = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_h  ))) # full hadron level
			if (not self.replaceKPpairs):
				jets_ch = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_hch))) # charged hadron level
			else:
				jets_ch = fj.sorted_by_pt(jet_selector(jet_def(track_selector_ch(parts_pythia_hch)))) # charged hadron level

			# R_label = str(jetR).replace('.', '') + 'Scaled'

			# Find the charged jet closest to the axis of the original parton
			# Require that the match is within some small angle, and that it is unique
			jet_matching_distance = 0.6  # Match jets with deltaR < jet_matching_distance*jetR
			self.parent0match, self.parent1match = None, None
			for i_jch, jch in enumerate(jets_ch):
				# Do constituent pT cut
#                pinfo("self.min_leading_track_pT", self.min_leading_track_pT)
# 				if self.min_leading_track_pT and not \
# 					self.utils.is_truth_jet_accepted(jch):
# #                   self.utils.is_truth_jet_accepted(jch, self.min_leading_track_pT):
# 					continue
				for i_parent, parent in enumerate(self.parents):
					parentmatch_name = "parent%imatch" % i_parent
					if jch.delta_R(parent) < jet_matching_distance * jetR:
						match = getattr(self, parentmatch_name)
						if not match:
							setattr(self, parentmatch_name, jch)
						else:  # Already found a match
							# Set flag value so that we know to ignore this one
							setattr(self, parentmatch_name, 0)

			# If we have matches, fill histograms
			for i_parent, parent in enumerate(self.parents):
				jet = getattr(self, "parent%imatch" % i_parent)
				if not jet:
#                    pinfo("in not jet")
					if jet == 0:
						# More than one match -- take note and continue
						count1 += 1
						continue
					else:  # jet == None
						# No matches -- take note and continue
						count2 += 1
						continue

				# One unique match
				# Identify the histograms which need to be filled
#                pinfo("passed not jet")
				parton_id = self.parent_ids[i_parent]
				parton_types = []
				if parton_id in self.quark_pdg_ids:
					# parton_types += ["quark"]
					if parton_id in self.charm_pdg_ids:
						parton_types += ["charm"]
					elif parton_id in self.up_pdg_ids or parton_id in self.down_pdg_ids or parton_id in self.strange_pdg_ids:
						parton_types += ["light"]
					elif (parton_id in self.beauty_pdg_ids and self.initscat == 4):
						parton_types += ["beauty"]
				elif parton_id in self.gluon_pdg_ids:
					parton_types += ["gluon"]
				parton_types += ["inclusive"]

				# If parent parton not identified, skip for now
				if not len(parton_types):
					continue
				
				if D0Kpidecayfound:
					print("parton types", parton_types)


				# Select for just D0-tagged jets #TODO: check if this D0 goes to kaon pion??
				D0taggedjet = False
				Dstartaggedjet = False
				if ( self.replaceKPpairs ):
					# print("There are ", len(jet.constituents()), "constituents.")
					for c in jet.constituents():
						constituent_pdg_idabs = pythiafjext.getPythia8Particle(c).idAbs()
						# constituent_pdg_index = c.user_index()
						if (constituent_pdg_idabs == 421): #TODO: this is assuming there is only one D0 per jet!
							if (self.checkDecayChannel(pythiafjext.getPythia8Particle(c), self.event) == EMesonDecayChannel.kDecayD0toKpi): # or self.checkDecayChannel(pythiafjext.getPythia8Particle(c), self.event) == EMesonDecayChannel.kDecayDStartoKpipi ):
								self.getD0Info(pythiafjext.getPythia8Particle(c))
								D0taggedjet = True
								break
							if (self.checkDecayChannel(pythiafjext.getPythia8Particle(c), self.event) == EMesonDecayChannel.kDecayDStartoKpipi):
								self.getD0Info(pythiafjext.getPythia8Particle(c))
								Dstartaggedjet = True
								break

					# prevent jets that are not dtagged or dstar tagged from continuing
					if ( not D0taggedjet and not Dstartaggedjet):
						# print("Not a D0 or D* jet")
						continue
							

					# print("booleans are", self.Dstar, Dstartaggedjet)
					if ( not self.Dstar and not self.D0wDstar ):
						if ( not D0taggedjet ): #if not a D0 tagged jet, move to next jet
							# print("Dstar is false, D0wDstar is false, and this is not D0tagged jet")
							continue
					if ( self.Dstar and not Dstartaggedjet ): #if only looking at D*s and D* is not tagged, move to next jet
						# print("Dstar is true and Dstar is not tagged")
						continue
					


				

				# Fill histograms
				# pinfo("if theres no list....", self.observable_list)
				for observable in self.observable_list:
					for i, jetaxisname in enumerate(self.jetaxes_observables[observable]): #looping through different configurations?

						obs_setting = self.obs_settings[observable][i] 
						grooming_setting = self.obs_grooming_settings[observable][i] if None != self.obs_grooming_settings[observable][i] else {} #[{'sd': [0.2, 0]}]
						# pinfo("GROOMING SETTING", grooming_setting)
						obs_label = self.utils.obs_label(obs_setting, grooming_setting)


						# Apply cut on leading track pT - does this apply to groomed jets too??
						leading_parton = fj.sorted_by_pt(jet.constituents())[0]
						leading_parton_pt = leading_parton.pt()
						if (leading_parton_pt < self.leading_parton_pt_cut):
							continue


						# TODO HERE
						# Groom jet, if applicable
						# jet_groomed_lund = None
						# if grooming_setting:
						# 	# gshop = fjcontrib.GroomerShop(jet, jetR, self.reclustering_algorithm)
						# 	gshop = fjcontrib.GroomerShop(j, jetR, fj.cambridge_algorithm)
						# 	jet_groomed_lund = self.utils.groom(gshop, grooming_setting, jetR)
						# 	if not jet_groomed_lund:
						# 		continue



						if (D0taggedjet): #D0Kpidecayfound):
							self.hD0KpiNjets.Fill(0)
						if (Dstartaggedjet): #self.DstarKpipidecayfound):
							self.hDstarNjets.Fill(0) 

					

	#                        obs = self.calculate_observable(
	#                            observable, jet, jet_groomed_lund, jetR, obs_setting,
	#                            grooming_setting, obs_label, jet.pt())
						obs = self.calculate_observable(
							jetaxisname, jet, jetR, jet.pt())
						
	#                            
						for parton_type in parton_types:
							# move filling thnsparse stuff outside the parton_types loop?
							self.fsparsepartonJetvalue[0] = jet.pt()
							self.fsparsepartonJetvalue[3] = obs
							if ( self.replaceKPpairs ):
								D0_px = self.D0particleinfo.px()
								D0_py = self.D0particleinfo.py()
								self.fsparsepartonJetvalue[1] = math.sqrt(D0_px*D0_px + D0_py*D0_py)
								self.fsparsepartonJetvalue[2] = self.D0particleinfo.y()
							else:
								self.fsparsepartonJetvalue[1] = -1
								self.fsparsepartonJetvalue[2] = -99

							getattr(self, ('h_%s_JetPt_%s_R%s_%s' % (observable, parton_type, jetR, obs_label)) if \
							len(obs_label) else ('h_%s_JetPt_%s_R%s' % (observable, parton_type, jetR))).Fill(self.fsparsepartonJetvalue)
							# getattr(self, ('h_JetPt_%s_R%s_%s_jetlevel' % (parton_type, jetR, obs_label)) if \
							#     len(obs_label) else ('h_JetPt_%s_R%s_jetlevel' % (parton_type, jetR))).Fill(jet.pt())

			setattr(self, "count1_R%s" % jetR_str, count1)
			setattr(self, "count2_R%s" % jetR_str, count2)

	#---------------------------------------------------------------
	# Calculate the observable given a jet
	#---------------------------------------------------------------
#    def calculate_observable(self, observable, jet, jet_groomed_lund,
#        jetR, obs_setting, grooming_setting, obs_label, jet_pt_ungroomed):
	def calculate_observable(self, observable, jet,
		jetR, jet_pt_ungroomed):

		# print("The observable is ", observable)

		if observable == "STD-D": #"STD-D":
			leading_parton = fj.sorted_by_pt(jet.constituents())[0]
			jet_axis_dif_STD_D0 = jet.delta_R(leading_parton)

			return jet_axis_dif_STD_D0
		
		elif observable == "WTA-D":
		   
			leading_parton = fj.sorted_by_pt(jet.constituents())[0]
			jet_wta = self.getWTAjet(jet, jetR)
			jet_axis_dif_WTA_D0 = jet_wta.delta_R(leading_parton)

			return jet_axis_dif_WTA_D0
		
		elif observable == "SD-D":
			leading_parton = fj.sorted_by_pt(jet.constituents())[0]
			jet_groomed_lund = self.getSDjet(jet, jetR)
			if jet_groomed_lund.Delta() < 0:
				jet_axis_dif_SD_D0 = -0.05
			else:
				jet_sd = jet_groomed_lund.pair()
				jet_axis_dif_SD_D0 = jet_sd.delta_R(leading_parton)
			return jet_axis_dif_SD_D0

		elif observable == "STD-WTA":
			leading_parton = fj.sorted_by_pt(jet.constituents())[0]
			jet_wta = self.getWTAjet(jet, jetR)
			jet_axis_dif_STD_WTA = jet.delta_R(jet_wta)

			return jet_axis_dif_STD_WTA
		
		elif observable == "WTA-SD":
			jet_wta = self.getWTAjet(jet, jetR)
			jet_groomed_lund = self.getSDjet(jet, jetR)
			if jet_groomed_lund.Delta() < 0:
				jet_axis_dif_WTA_SD = -0.05
			else:
				jet_sd = jet_groomed_lund.pair()
				jet_axis_dif_WTA_SD = jet_wta.delta_R(jet_sd)
			return jet_axis_dif_WTA_SD
		
		elif observable == "STD-SD":
			jet_groomed_lund = self.getSDjet(jet, jetR)
			if jet_groomed_lund.Delta() < 0:
				jet_axis_dif_STD_SD = -0.05
			else:
				jet_sd = jet_groomed_lund.pair()
				jet_axis_dif_STD_SD = jet.delta_R(jet_sd)
			return jet_axis_dif_STD_SD

		#    if grooming_setting:
		#        j_groomed = jet_groomed_lund.pair()
		#        if not j_groomed.has_constituents():
		#            # Untagged jet -- record underflow value
		#            return -1
		#        else:
		#            return j_groomed.m()

		# Should not be any other observable
		raise ValueError("Observable %s not implemented" % observable)


	def checkDecayChannel(self, particle, event): #(part, mcArray): # what type is part
		
		if(not event):
			return EMesonDecayChannel.kUnknownDecay
	 
		decay = EMesonDecayChannel.kUnknownDecay 

		absPdgPart = particle.idAbs()
		
		if(len(particle.daughterList()) == 2):
			d1_index = particle.daughterList()[0] #don't use daughter1() and daughter(2)
			d2_index = particle.daughterList()[1]
			d1 = event[d1_index]
			d2 = event[d2_index]

			if(not d1 or not d2):
				return decay

			# print("checkpoint 3")
	

			absPdg1 = d1.idAbs()
			absPdg2 = d2.idAbs()

			if(absPdgPart == 421):  # D0 -> K pi
				if((absPdg1 == 211 and absPdg2 == 321) or (absPdg1 == 321 and absPdg2 == 211)): # pi K or K pi - QUESTION: does this account for k and pi being opposite signs?
					decay = EMesonDecayChannel.kDecayD0toKpi 

			# TODO: can insert if (self.Dstar) later
	  
			# Look at D0's mother particles
			# print("current particle ID is", absPdgPart)
			mother_indices = particle.motherList()
			if (len(mother_indices) != 1):
				return decay #just return D0->Kpi because D0 didn't come from a D*
			# print("MOTHERS", len(mother_indices)) # there's a lot of these...
			# print(mother_indices)
			for mother_index in mother_indices:
				mother = event[mother_index]
				absPdg_mother = mother.idAbs()

				if (absPdg_mother == 413): # if mother is D*+/-
					# if (len(mother_indices != 1)):
					#     print("There were", len(mother_indices), "mothers in this event!")
					# look at daughters of mother
					if(len(mother.daughterList()) == 2):
						d1_index = mother.daughterList()[0] #don't use daughter1() and daughter(2)
						d2_index = mother.daughterList()[1]
						d1 = event[d1_index]
						d2 = event[d2_index]
						if(not d1 or not d2):
							return decay                
						absPdg1 = d1.idAbs()
						absPdg2 = d2.idAbs()

						if((absPdg1 == 421 and absPdg2 == 211) or (absPdg1 == 211 and absPdg2 == 421)): # D0 pi or pi D0
							decay = EMesonDecayChannel.kDecayDStartoKpipi 
							break #TODO: should this break be earlier? is it possible to have multiple mothers that are D*?
					
			# print(event)

		return decay


	# save D0 particle info to save to THnSparse
	def getD0Info(self, particle): 
		self.D0particleinfo = particle
		return
	
	def getWTAjet(self, jet, jetR):

		jet_def_wta = fj.JetDefinition(fj.cambridge_algorithm, 2*jetR)#originally 2*self.jetR_list[0])
		jet_def_wta.set_recombination_scheme(fj.WTA_pt_scheme)
		reclusterer_wta =  fjcontrib.Recluster(jet_def_wta)
		jet_wta = reclusterer_wta.result(jet)
		
		return jet_wta
	
	def getSDjet(self, jet, jetR):
		 
		#  if grooming_setting:
		gshop = fjcontrib.GroomerShop(jet, jetR, fj.cambridge_algorithm)
		# jet_groomed_lund = self.utils.groom(gshop, grooming_setting, jetR)
		jet_groomed_lund = gshop.soft_drop(self.sd_beta,self.sd_zcut,jetR)
		
		
		# jet_groomed = jet_groomed_lund.pair()


		# if jet_groomed_lund.Delta() < 0:
		# 	self.overflow = True
		# 	return jet_groomed
		# else:
		# 	self.overflow = False
		# 	return jet_groomed
		return jet_groomed_lund



		# if not jet_groomed_lund:
		# 	continue
		 


	#---------------------------------------------------------------
	# Initiate scaling of all histograms and print final simulation info
	#---------------------------------------------------------------
	def scale_print_final_info(self, pythia):
		# Scale all jet histograms by the appropriate factor from generated cross section and the number of accepted events
		scale_f = pythia.info.sigmaGen() / self.hNevents.GetBinContent(1)

		for jetR in self.jetR_list:
			hist_list_name = "hist_list_R%s" % str(jetR).replace('.', '')
			for h in getattr(self, hist_list_name):
			#     if 'jetlevel' in h.GetTitle():
			#         continue
				h.Scale(scale_f)

		print("N total final events:", int(self.hNevents.GetBinContent(1)), "with",
			  int(pythia.info.nAccepted() - self.hNevents.GetBinContent(1)),
			  "events rejected at hadronization step")
		self.hNevents.SetBinError(1, 0)
		self.hD0Nevents.SetBinError(1, 0)
		self.hD0KpiNevents.SetBinError(1, 0)
		self.hD0KpiNjets.SetBinError(1, 0)
		self.hDstarNjets.SetBinError(1, 0)

################################################################
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly',
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
	parser.add_argument('--nocharmdecay', action='store', type=int, default=0, help="'1' turns charm decays off")
	parser.add_argument('--weightON', action='store', type=int, default=0, help="'1' turns weights on")
	parser.add_argument('--leadingptcut', action='store', type=float, default=0, help="leading track pt cut")
	parser.add_argument('--replaceKP', action='store', type=int, default=0, help="'1' replaces the K/pi pairs with D0")
	parser.add_argument('--DstarON', action='store', type=int, default=0, help="'1' looks at EEC for D* only")
	parser.add_argument('--chinitscat', action='store', type=int, default=0, help="'0' runs all events, \
						'1' runs only hard->ccbar events, '2' runs only gg->ccbar events, '3' runs only D0->Kpi events")
	parser.add_argument('--D0withDstarON', action='store', type=int, default=0, help="'1' looks at EEC for D0 and D0 from D*")
	

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
