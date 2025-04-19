
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

		# self implemented variables to study
		self.weighted = (bool)(args.weightON) #weightON=True(F) means turn weights on(off)
		self.leading_parton_pt_cut = args.leadingptcut
		
		self.replaceKPpairs = (bool)(args.replaceKP) #replaceKP=True(F) means turn k/pi pairs are('nt) replaced
		self.softqcd = args.softqcd # do all softqcd scatterings

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

		# self.pt_bins = array.array('d', list(range(5, 100, 5)) + list(range(100, 210, 10)))
#        self.obs_bins_ang = np.concatenate((np.linspace(0, 0.009, 10), np.linspace(0.01, 0.1, 19),
#                                            np.linspace(0.11, 0.8, 70)))
#        self.obs_bins_mass = np.concatenate(
#          (np.linspace(0, 0.9, 10), np.linspace(1, 9.8, 45), np.linspace(10, 14.5, 10),
#           np.linspace(15, 19, 5), np.linspace(20, 60, 9)))
		self.obs_bins_EEC = np.logspace(np.log10(1E-4), np.log10(1), 51)

		self.observable_list = config['process_observables']
		self.obs_settings = {}
		self.obs_grooming_settings = {}
		self.obs_names = {}
		for observable in self.observable_list:

			obs_config_dict = config[observable]
			obs_config_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]

			self.obs_settings[observable] = self.utils.obs_settings(observable, obs_config_dict, obs_config_list)
			pinfo("self.obs_settings[observable]", self.obs_settings[observable])
			self.obs_grooming_settings[observable] = self.utils.grooming_settings(obs_config_dict)
			pinfo("self.obs_grooming_settings[observable]", self.obs_grooming_settings[observable])

			self.obs_names[observable] = obs_config_dict["common_settings"]["xtitle"]

	#---------------------------------------------------------------
	# Main processing function
	#---------------------------------------------------------------
	def pythia_quark_gluon(self, args):

		''' # this should be done in process_base.py
		# Create ROOT TTree file for storing raw PYTHIA particle information
		outf_path = os.path.join(self.output_dir, args.tree_output_fname)
		fout = ROOT.TFile(outf_path, 'recreate') #restructured for new file save
		fout.cd() #restructured for new file save
		'''
		# Initialize response histograms
		self.initialize_hist()

		pinfo('user seed for pythia', self.user_seed) #TODO: what does this do?? it doesn't work...
#        print('user seed for pythia', self.user_seed)
		mycfg = ['Random:setSeed=on', 'Random:seed={}'.format(self.user_seed)]
		mycfg.append('HadronLevel:all=off')
		

		if (self.softqcd):
			mycfg.append('HardQCD:all = off')
			mycfg.append('SoftQCD:all = on')

		if (self.replaceKPpairs):
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

		# outf.Write()
		# outf.Close()
		# self.fout.Write() #restructued for new file structure

		self.save_output_objects()

	#---------------------------------------------------------------
	# Initialize histograms
	#---------------------------------------------------------------
	def initialize_hist(self):

		self.fout.cd() # new addition for new file save format

		self.hNevents = ROOT.TH1I("hNevents", 'Number accepted events (unscaled)', 2, -0.5, 1.5)
		self.hD0Nevents = ROOT.TH1I("hD0Nevents", "Total Number of D0 events (unscaled)", 2, -0.5, 1.5)
		self.hD0KpiNevents = ROOT.TH1I("hD0KpiNevents", "Number of D0->Kpi events (unscaled)", 2, -0.5, 1.5)
		self.hD0KpiNjets = ROOT.TH1I("hD0KpiNjets", "Number of D0->Kpi jets (unscaled)", 2, -0.5, 1.5) #accidentally called "hD0KpiNehD0KpiNjetsvents"
		self.hDeltaR = ROOT.TH1F("hDeltaR", 'Delta R between jet and each parent', 40, 0, 0.4)
		self.hnumconst = ROOT.TH1I("hnumconst", "Number of constituents per jet (unscaled)", 50, 0, 50)
		self.hnumconstwTrackcut = ROOT.TH1I("hnumconstwTrackcut", "Number of constituents per jet with the track cut (unscaled)", 50, 0, 50)
		self.hD0pT = ROOT.TH1F("hD0pT", 'pt of D0', 200, 0., 200.)
		self.hD0z = ROOT.TH2F("hD0z", 'z of D0;D0 z;jet pT', 100, 0., 1.01, 200, 0., 200.)
		# self.hNgluonsplitjets = ROOT.TH1I("hNgluonsplitjets", "Number of jets that come from gg->ccbar (unscaled)", 2, -0.5, 1.5)
		# self.hD0pT_gluonsplit = ROOT.TH1F("hD0pT_gluonsplit", 'pt of D0s that come from gluon splitting', 200, 0., 200.)
		# self.hD0z_gluonsplit = ROOT.TH2F("hD0z_gluonsplit", 'z of D0s that come from gluon splitting;D0 z;jet pT', 100, 0., 1.01, 200, 0., 200.)
		self.hNevents_nocollision = ROOT.TH1I("hNevents_nocollision", 'Number of events where pp did not collide;N_{events};Counts', 2, 0, 2)
		self.hpartons_from_2D0events = ROOT.TH2I("hpartons_from_2D0events", 'Initiating Parton of events with 2 D^{0}s;abs(Parton 1 PID);abs(Parton 2 PID)', 25, 0, 25, 25, 0, 25)
		self.hnumD0InEvent = ROOT.TH1I("hnumD0InEvent", 'Number of D^{0}s in each event;# of D^{0}s;Counts', 5, 0, 5)
		self.hnumD0_injets_InEvent = ROOT.TH1I("hnumD0_injets_InEvent", 'Number of D^{0}s in the saved jets in each event;# of D^{0}s;Counts', 5, 0, 5)

		for jetR in self.jetR_list:

			# Store a list of all the histograms just so that we can rescale them later
			hist_list_name = "hist_list_R%s" % str(jetR).replace('.', '')
			setattr(self, hist_list_name, [])

			R_label = str(jetR).replace('.', '') + 'Scaled'

			for observable in self.observable_list:
				# Should only be one: observable == "EEC"
				if observable != "EEC":
					raise ValueError("Observable %s is not implemented in this script" % observable)

				obs_bins = getattr(self, "obs_bins_" + observable)
				# Use more finely binned pT bins for TH2s than for the RMs
				pt_bins = array.array('d', list(range(0, 201, 1)))
				rapi_bins = np.linspace(-5,5,201)
				z_bins = np.linspace(0, 1.01, 102)


				dim = 5
				nbins  = [len(pt_bins)-1, len(pt_bins)-1, len(rapi_bins)-1, len(z_bins)-1, 50]
				min_li = [pt_bins[0],     pt_bins[0],      rapi_bins[0],      obs_bins[0],      z_bins[0]]
				max_li = [pt_bins[-1],    pt_bins[-1],     rapi_bins[-1],     obs_bins[-1],     z_bins[-1]]

				nbins = (nbins)
				xmin = (min_li)
				xmax = (max_li)
				
				nbins_array = array.array('i', nbins)
				xmin_array = array.array('d', xmin)
				xmax_array = array.array('d', xmax)

				# Loop over subobservable (alpha value)
#                for i in range(len(self.obs_settings[observable])):

				obs_setting = self.obs_settings[observable]
				grooming_setting = self.obs_grooming_settings[observable]
				obs_label = self.utils.obs_label(obs_setting, grooming_setting)
				pinfo("all the settings", obs_setting, grooming_setting, obs_label)



				self.fsparsepartonJetvalue = array.array( 'd', ( 0, 0, 0, 0, 0 ))
				self.fsparsejetlevelJetvalue = array.array( 'd', ( 0, 0, 0, 0 ))
		
				partontypeslist = ["charm", "light", "gluon", "gluon2c", "beauty", "inclusive", "anything2c"] #got rid of quark

				for parton_type in partontypeslist:

					title = [ '#it{p}_{T}^{ch jet}', '#it{p}_{T}^{#phi}', 'y', 'z', '#it{R}_{L}' ]

					# make THnSparse for parton EECs
					name = ('hsparse_%s_JetPt_%s_R%s_%s' % (observable, parton_type, jetR, obs_label)) if \
						len(obs_label) else ('h_%s_JetPt_%s_R%s' % (observable, parton_type, jetR))
					hsparse = ROOT.THnSparseD(name,"%s-init_hsparsejet; #it{p}_{T,%s}^{ch jet}; #it{p}_{T}^{D^{0}}; y;R_{L}^{%s}" %(parton_type[0], parton_type[0] + "-init", parton_type[0] + "-init"), dim,  nbins_array, xmin_array, xmax_array)
					# hsparse.GetXaxis().SetTitle('#it{p}_{T,%s}^{ch jet}' % (parton_type[0] + "-init"))
					# hsparse.GetYaxis().SetTitle("R_{L}" + '^{%s}' % (parton_type[0] + "-init"))
					hsparse.Sumw2()
					for i in range(0,dim):
						hsparse.GetAxis(i).SetTitle(title[i])
						if i == 0 or i == 1:
							hsparse.SetBinEdges(i, pt_bins)
						if i == 2:
							hsparse.SetBinEdges(i, rapi_bins)
						if i == 3:
							hsparse.SetBinEdges(i, z_bins)
						if i == 4:
							hsparse.SetBinEdges(i, obs_bins)
					setattr(self, name, hsparse)
					getattr(self, hist_list_name).append(hsparse)



					# make another of THnSparse for the jet level (above is pair level)
					name_jetpt = ('h_JetPt_%s_R%s_%s_jetlevel' % (parton_type, jetR, obs_label)) if \
						len(obs_label) else ('h_JetPt_%s_R%s_jetlevel' % (parton_type, jetR))
					hsparse_jetpt = ROOT.THnSparseD(name_jetpt,"%s-init_hsparsejet_jetlevel; #it{p}_{T,%s}^{ch jet}; #it{p}_{T}^{D^{0}}; y" %(parton_type[0], parton_type[0] + "-init"), dim-1,  nbins_array[:-1], xmin_array[:-1], xmax_array[:-1])
					# hsparse_jetpt.GetXaxis().SetTitle('#it{p}_{T,%s}^{ch jet}' % (parton_type[0] + "-init"))
					# hsparse_jetpt.GetYaxis().SetTitle('Counts')
					hsparse_jetpt.Sumw2()
					for i in range(0,dim-1):
						hsparse_jetpt.GetAxis(i).SetTitle(title[i])
						if i == 0 or i == 1:
							hsparse_jetpt.SetBinEdges(i, pt_bins)
						if i == 2:
							hsparse_jetpt.SetBinEdges(i, rapi_bins)
						if i == 3:
							hsparse_jetpt.SetBinEdges(i, z_bins)
					setattr(self, name_jetpt, hsparse_jetpt)
					getattr(self, hist_list_name).append(hsparse_jetpt)


	#---------------------------------------------------------------
	# Initiate jet defs, selectors, and sd (if required)
	#---------------------------------------------------------------
	def init_jet_tools(self):

		for jetR in self.jetR_list:
			jetR_str = str(jetR).replace('.', '')

			# set up our jet definition and a jet selector
			jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
			setattr(self, "jet_def_R%s" % jetR_str, jet_def)

		pwarning('max eta for particles after hadronization set to', self.max_eta_hadron)
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

			# note: iev only counts the number of "successful" events -- or those that find jets
			if (iev%1000 == 0): #iev%10000 == 0):
				print("Event", iev)
			# print("Event", iev)

			self.parents = []
			self.event = pythia.event
			# print("Event table!", type(pythia.event), "//", self.event) # to print out a table of the event information

			# make sure the protons have collided
			# if (iev == 3):
			# 	print("checking", pythia.event.size())
			if pythia.event.size()<=5:
				# print("Event", iev, ": Protons whiffed by -- no generated partons :(")
				self.hNevents_nocollision.Fill(1)
				continue

			fs_parton_5 = fj.PseudoJet(pythia.event[5].px(), pythia.event[5].py(), pythia.event[5].pz(), pythia.event[5].e())
			fs_parton_6 = fj.PseudoJet(pythia.event[6].px(), pythia.event[6].py(), pythia.event[6].pz(), pythia.event[6].e())
			self.parents = [fs_parton_5, fs_parton_6] # parent partons in dijet

			# Save PDG code of the parent partons
			self.parent_ids = [pythia.event[5].id(), pythia.event[6].id()]

			# parton level
			# if (self.partonlevel):
			# 	parts_pythia_p = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kVisible], 0, True) # parton level, full jets
			# else:
			parts_pythia_p = None

			# Force hadronization here
			hstatus = pythia.forceHadronLevel()
			if not hstatus:
				continue

			# full-hadron level
			parts_pythia_h = None
			

			# charged-hadron level
			if ( self.replaceKPpairs == False ):
				parts_pythia_hch = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kVisible, pythiafjext.kCharged], 0, True)		
			else: #replace D0->Kpi
				parts_pythia_hch = pythiafjext.vectorize_select_replaceD0(pythia, [pythiafjext.kFinal, pythiafjext.kVisible, pythiafjext.kCharged], 0, True, True)
				
			# look at events in charged hadron??
			# print("!! pythia hadron (after vectorization) event size is ", pythia.event.size())
			#TODO: move this block above choosing final state charged particles??
			self.particlecounter = 0
			self.D0_particle_list = []
			D0found = False
			D0Kpidecayfound = False
			# self.DstarKpipidecayfound = False
			for particle in self.event: #for event in pythia.event:
				if particle.id() == 421 or particle.id() == -421: #D0
					D0found = True
					self.D0_particle_list.append(particle)
				#     print(self.particlecounter, "D0 with particle id", particle.id())
					if self.checkDecayChannel(particle, self.event) == EMesonDecayChannel.kDecayD0toKpi:
						# print(self.particlecounter, "This is a D0->Kpi decay!", particle.id())
						D0Kpidecayfound = True

					# if self.checkDecayChannel(particle, self.event) == EMesonDecayChannel.kDecayDStartoKpipi:
					# 	self.DstarKpipidecayfound = True #can't fill histogram here because it will fill at particle level
					# 	# print("Dstar found!") 


				self.particlecounter+=1

			if not D0found:
				continue # at this point, if there is no D0, just move on and make a new event!!!


			#if D0->Kpi found, count the events; if not, check that length of charged final state hadrons vector is 0
			if (D0Kpidecayfound):
				self.hD0KpiNevents.Fill(0)
			if (D0found):
				self.hD0Nevents.Fill(0)

			self.hnumD0InEvent.Fill(len(self.D0_particle_list))


			# Some "accepted" events don't survive hadronization step -- keep track here
			self.hNevents.Fill(0)
			self.find_jets_fill_histograms(parts_pythia_hch, iev, D0Kpidecayfound, parts_pythia_h, parts_pythia_p)

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
	def find_jets_fill_histograms(self, parts_pythia_hch, iev, D0Kpidecayfound, parts_pythia_h, parts_pythia_p):
		# Loop over jet radii
		for jetR in self.jetR_list:

			jetR_str = str(jetR).replace('.', '')
			jet_selector = getattr(self, "jet_selector_R%s" % jetR_str)
			jet_def = getattr(self, "jet_def_R%s" % jetR_str)
			track_selector_ch = getattr(self, "track_selector_ch")

			count1 = getattr(self, "count1_R%s" % jetR_str)
			count2 = getattr(self, "count2_R%s" % jetR_str)

			# Get the jets at different levels
			jets_ch = fj.sorted_by_pt(jet_selector(jet_def(track_selector_ch(parts_pythia_hch)))) # charged hadron level


			# Find the charged jet closest to the axis of the original parton
			# Require that the match is within some small angle, and that it is unique
			jet_matching_distance = 0.6  # Match jets with deltaR < jet_matching_distance*jetR
			self.parent0match, self.parent1match = None, None
			# if len(jets_ch) > 0:
			# 	print("number of jets", len(jets_ch))
			for i_jch, jch in enumerate(jets_ch):
				for i_parent, parent in enumerate(self.parents):
					parentmatch_name = "parent%imatch" % i_parent
					self.hDeltaR.Fill(jch.delta_R(parent))
					if jch.delta_R(parent) < jet_matching_distance * jetR:
						print("match found!!!!")
						match = getattr(self, parentmatch_name)
						if not match:
							setattr(self, parentmatch_name, jch)
						else:  # Already found a match
							# Set flag value so that we know to ignore this one
							setattr(self, parentmatch_name, 0)


			# If we have matches, fill histograms
			self.D0injet_partonparent_list = []
			for i_parent, parent in enumerate(self.parents):
				jet = getattr(self, "parent%imatch" % i_parent)
				if not jet:
					if jet == 0: # More than one match -- take note and continue
						count1 += 1
						continue
					else:  # jet == None
						# No matches -- take note and continue
						count2 += 1
						continue


				# One unique match
				# Identify the histograms which need to be filled
				parton_id = self.parent_ids[i_parent]
				parton_types = []
				if parton_id in self.quark_pdg_ids:
					# parton_types += ["quark"]
					if parton_id in self.charm_pdg_ids:
						print("lookartmeeee, charm found!!")
						parton_types += ["charm"]
					elif parton_id in self.up_pdg_ids or parton_id in self.down_pdg_ids or parton_id in self.strange_pdg_ids:
						parton_types += ["light"]
					elif (parton_id in self.beauty_pdg_ids):
						parton_types += ["beauty"]
				elif parton_id in self.gluon_pdg_ids:
					parton_types += ["gluon"]
				# if not self.replaceKPpairs:
				parton_types += ["inclusive"]

				if self.particlecounter == 2:
					parton_types += ["anything2c"]
					if parton_id in self.gluon_pdg_ids:
						parton_types += ["gluon2c"]


				# If parent parton not identified, skip for now
				if not len(parton_types):
					continue
				
				# print(D0Kpidecayfound)
				if D0Kpidecayfound:
					print("parton types", parton_types)


				# Select for just D0-tagged jets #TODO: check if this D0 goes to kaon pion??
				D0taggedjet = False
				# Dstartaggedjet = False
				if ( self.replaceKPpairs ):
					print("There are ", len(jet.constituents()), "constituents.")
					for c in jet.constituents():
						constituent_pdg_idabs = pythiafjext.getPythia8Particle(c).idAbs()
						constituent_pdg_index = c.user_index()
						# print("D0 index from pythiafjext", pythiafjext.getPythia8Particle(c).index())
						# print("D0 user_index from pythiafjext", c.user_index())
						if (constituent_pdg_idabs == 421): #TODO: this is assuming there is only one D0 per jet!
							print("The decay channel is ", self.checkDecayChannel(pythiafjext.getPythia8Particle(c), self.event))
							if (self.checkDecayChannel(pythiafjext.getPythia8Particle(c), self.event) == EMesonDecayChannel.kDecayD0toKpi): # or self.checkDecayChannel(pythiafjext.getPythia8Particle(c), self.event) == EMesonDecayChannel.kDecayDStartoKpipi ):
								# print("Check the momentum!", pythiafjext.getPythia8Particle(c).px(), pythiafjext.getPythia8Particle(c).py())
								self.getD0Info(pythiafjext.getPythia8Particle(c))
								# print("D0 index from pythiafjext", pythiafjext.getPythia8Particle(c).index())
								# print("D0 user_index from pythiafjext", c.user_index())
								
								# save information on potential multiple D0s per event
								self.D0injet_partonparent_list.append(parton_types[0])
								print("ADDING PARTON HERE!", parton_types[0], self.D0injet_partonparent_list)

								D0taggedjet = True
								break
							# if (self.checkDecayChannel(pythiafjext.getPythia8Particle(c), self.event) == EMesonDecayChannel.kDecayDStartoKpipi):
							# 	self.getD0Info(pythiafjext.getPythia8Particle(c))
							# 	Dstartaggedjet = True

							# 	break

					# TODO: need to prevent jets that are not dtagged or dstar tagged from 
					if ( not D0taggedjet): # and not Dstartaggedjet):
						# print("Not a D0 or D* jet")
						continue


				# Fill histograms
				for observable in self.observable_list:
#                    pinfo("len(self.obs_settings[observable])", len(self.obs_settings[observable]))

					obs_setting = self.obs_settings[observable]
					grooming_setting = self.obs_grooming_settings[observable]
					obs_label = self.utils.obs_label(obs_setting, grooming_setting)

					# Groom jet, if applicable
					jet_groomed_lund = None
					if grooming_setting:
						gshop = fjcontrib.GroomerShop(jet, jetR, self.reclustering_algorithm)
						jet_groomed_lund = self.utils.groom(gshop, grooming_setting, jetR)
						if not jet_groomed_lund:
							continue


					# Apply cut on leading track pT
					# print("num jet constituents", len(jet.constituents()))
					leading_parton = fj.sorted_by_pt(jet.constituents())[0]
					leading_parton_pt = leading_parton.pt()
					if (leading_parton_pt < self.leading_parton_pt_cut):
						continue



					# count the number of D0-tagged jets. If the observable is not EEC, might have to change where this is
					if (D0taggedjet): #D0Kpidecayfound):
						self.hD0KpiNjets.Fill(0)
					# if (Dstartaggedjet): #self.DstarKpipidecayfound):
					# 	self.hDstarNjets.Fill(0) 


					# print("filling jet level thnsparse")
					# fill jet pt histogram to give the normalization
					self.fsparsejetlevelJetvalue[0] = jet.pt()
					if ( self.replaceKPpairs): # phimeson has bad naming convention but is properly filled here
						D0_px = self.D0particleinfo.px()
						D0_py = self.D0particleinfo.py()
						D0_pt = math.sqrt(D0_px*D0_px + D0_py*D0_py)
						# print("momentum confirmed", D0_px, D0_py)
						self.fsparsejetlevelJetvalue[1] = D0_pt
						self.fsparsejetlevelJetvalue[2] = self.D0particleinfo.y()
						self.fsparsejetlevelJetvalue[3] = D0_pt/jet.pt()

						# D0 information
						self.hD0pT.Fill(D0_pt)
						self.hD0z.Fill(D0_pt/jet.pt(), jet.pt())

						
					else:
						self.fsparsejetlevelJetvalue[1] = -1
						self.fsparsejetlevelJetvalue[2] = -99
						self.fsparsejetlevelJetvalue[3] = -99

					for parton_type in parton_types:
						getattr(self, ('h_JetPt_%s_R%s_%s_jetlevel' % (parton_type, jetR, obs_label)) if \
							len(obs_label) else ('h_JetPt_%s_R%s_jetlevel' % (parton_type, jetR))).Fill(self.fsparsejetlevelJetvalue)
					

					# Fill number of constituents per jet
					self.hnumconst.Fill(len(jet.constituents()))


					obs = self.calculate_observable(
						observable, jet, jet_groomed_lund, jetR, jet.pt())


					# print("filling pair level thnsparse")
					for index in range(obs.correlator(2).rs().size()):
						self.fsparsepartonJetvalue[0] = jet.pt()
						self.fsparsepartonJetvalue[4] = obs.correlator(2).rs()[index]
						if ( self.replaceKPpairs): # phimeson has bad naming convention but is properly filled here
							self.fsparsepartonJetvalue[1] = D0_pt
							self.fsparsepartonJetvalue[2] = self.D0particleinfo.y()
							self.fsparsepartonJetvalue[3] = D0_pt/jet.pt()
						else:
							self.fsparsepartonJetvalue[1] = -1
							self.fsparsepartonJetvalue[2] = -99
							self.fsparsepartonJetvalue[3] = -99

						for parton_type in parton_types:
							#fill parton hnsparse info
							if self.weighted:
								getattr(self, ('h_%s_JetPt_%s_R%s_%s' % (observable, parton_type, jetR, obs_label)) if \
									len(obs_label) else ('h_%s_JetPt_%s_R%s' % (observable, parton_type, jetR))).Fill(self.fsparsepartonJetvalue, obs.correlator(2).weights()[index])
							else:
								getattr(self, ('h_%s_JetPt_%s_R%s_%s' % (observable, parton_type, jetR, obs_label)) if \
									len(obs_label) else ('h_%s_JetPt_%s_R%s' % (observable, parton_type, jetR))).Fill(self.fsparsepartonJetvalue)
							
			setattr(self, "count1_R%s" % jetR_str, count1)
			setattr(self, "count2_R%s" % jetR_str, count2)

			# Fill information of potential multiple D0s per event
			self.hnumD0_injets_InEvent.Fill(len(self.D0injet_partonparent_list))
			if (len(self.D0injet_partonparent_list) == 2):
				self.hpartons_from_2D0events.Fill(self.D0injet_partonparent_list[0], self.D0injet_partonparent_list[1])
			
	#---------------------------------------------------------------
	# Calculate the observable given a jet
	#---------------------------------------------------------------
#    def calculate_observable(self, observable, jet, jet_groomed_lund,
#        jetR, obs_setting, grooming_setting, obs_label, jet_pt_ungroomed):
	def calculate_observable(self, observable, jet, jet_groomed_lund,
		jetR, jet_pt_ungroomed):

		if observable == "EEC":
			
			# Extract information for EEC
			constituents = fj.sorted_by_pt(jet.constituents())
			c_select = fj.vectorPJ()
			trk_thrd = 1 # track pt threshold
				
			# apply pT threshold on jet constituents
			for c in constituents:
				if c.pt() < trk_thrd:
					break
				#print("constituent used for pair =", c)
				c_select.append(c)
				
			dphi_cut = -9999
			deta_cut = -9999

			
			# Fill num const after track cut
			self.hnumconstwTrackcut.Fill(len(c_select))
			
			new_corr = ecorrel.CorrelatorBuilder(c_select, jet.perp(), 2, 1, dphi_cut, deta_cut)
			
			return new_corr #new_corr.correlator(2).rs()[index]
			

		# Should not be any other observable
		raise ValueError("Observable %s not implemented" % observable)
	

	def find_parts_around_jet(self, parts, jet, cone_R):
		# select particles around jet axis
		cone_parts = fj.vectorPJ()
		for part in parts:
			if jet.delta_R(part) <= cone_R:
				cone_parts.push_back(part)
		
		return cone_parts

	
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


	def getParticleAsPseudojet(self, particle):
		psjet = fj.PseudoJet(particle.px(), particle.py(), particle.pz(), particle.e())

		psjet.set_user_index(particle.index()) #should be + user_index_offset but that is 0
		# _pinfo = PythiaParticleInfo(pythia.event[particle.index()])
		# psjet.set_user_info(_pinfo)
		
		return psjet

	
	# check if D0 is prompt - should only send D0 that does not come from D* here
	# also assuming that D0's mother is c (direct mother, not with other generations in between)
	def checkPrompt(self, D0particle, event):

		promptness = Promptness.kUnknown

		absPdgPart = D0particle.idAbs()
		motherlist_indices = D0particle.motherList()
		# if (len(motherlist_indices) != 1):
		#     return  
		print("D0's mothers", motherlist_indices)
		for mother_index in motherlist_indices:
			mother = event[mother_index]
			absPdg_mother = mother.idAbs()
			print("D0 mother ID", absPdg_mother)

			if (absPdg_mother == 4): #charm
				# check if mother of charm is beauty
				charms_mother_indices = mother.motherList()
				print("charm's mothers", charms_mother_indices)

				# if there are no mothers???
				if len(charms_mother_indices) == 0:
					promptness = Promptness.kPrompt
					break

				for charms_mother_index in charms_mother_indices:
					charms_mother = event[charms_mother_index]
					absPdg_charms_mother = charms_mother.idAbs()
					print("charm mother ID", absPdg_charms_mother)

					if (absPdg_charms_mother == 4): #charm
						promptness = Promptness.kPrompt
						break
					if (absPdg_charms_mother == 5): #beauty
						promptness = Promptness.kNonPrompt
						break
					#else: would be unknown (if c's parentage is something but not a b...??)
				break

		return promptness


			

	
	def printD0mothers(self, particle, event, num):
		# if num == 10: #break statement
		#         print("Exited with num=10 ")
		#         return

		print("This is generation", num)

		motherlist_indices = particle.motherList()
		motherlist = [event[i].name() for i in motherlist_indices]
		motherlist_status = [event[i].status() for i in motherlist_indices]
		print("The indices are", motherlist_indices)
		print("The mothers are", motherlist)
		print("The statuss are", motherlist_status)

		if len(motherlist_indices) == 0:
			return

		for mother_index in motherlist_indices:

			# if mother_index < 5: #break statement
			#     print("Exited with mother_index of ", mother_index)
			#     break

			mother = event[mother_index]
			print("Following mother ", mother.name(), "with index", mother_index)
			self.printD0mothers(mother, event, num+1)



	# check if D0's mother is D*
	def checkD0motherIsDstar(self, D0particle, event):
		motherisDstar = False
		
		if (D0particle.idAbs() == 421): #D0
		
			mother_indices = D0particle.motherList()
			if len(mother_indices) == 1: # assuming D* is the only mother to D0
				mo1 = mother_indices[0]
				if event[mo1].idAbs() == 413: #D*
					motherisDstar = True
	
		# std::cout << "is mother a Dstar?  " << motherisDstar << std::endl;
		return motherisDstar
	

	

	#---------------------------------------------------------------
	# Initiate scaling of all histograms and print final simulation info
	#---------------------------------------------------------------
	def scale_print_final_info(self, pythia):
		# Scale all jet histograms by the appropriate factor from generated cross section and the number of accepted events
		scale_f = pythia.info.sigmaGen() / self.hNevents.GetBinContent(1)
		print("pythia.info.sigmaGen() is", pythia.info.sigmaGen())
		print("scale_f is", scale_f)
		print("int(pythia.info.nAccepted())", int(pythia.info.nAccepted()))

		for jetR in self.jetR_list:
			hist_list_name = "hist_list_R%s" % str(jetR).replace('.', '')
			# print(hist_list_name)
			for h in getattr(self, hist_list_name):
				h.Scale(scale_f)

		print("N total final events:", int(self.hNevents.GetBinContent(1)), "with",
			  int(pythia.info.nAccepted() - self.hNevents.GetBinContent(1)),
			  "events rejected at hadronization step")
		self.hNevents.SetBinError(1, 0)
		self.hD0Nevents.SetBinError(1, 0)
		self.hD0KpiNevents.SetBinError(1, 0)
		self.hD0KpiNjets.SetBinError(1, 0)


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
	parser.add_argument('--weightON', action='store', type=int, default=1, help="'1' turns weights on, '0' turns them off")
	parser.add_argument('--leadingptcut', action='store', type=float, default=0, help="leading track pt cut")
	parser.add_argument('--replaceKP', action='store', type=int, default=0, help="'1' replaces the K/pi pairs with D0")
	parser.add_argument('--softqcd', action='store', type=int, default=0, help="run all softqcd processes")

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

	process = PythiaQuarkGluon(config_file=args.config_file, output_dir=args.output_dir, args=args)
	process.pythia_quark_gluon(args)
