
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

from pyjetty.mputils import *
from pyjetty.mputils.mputils import pinfo, pwarning

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext

from pyjetty.alice_analysis.process.base import process_base

from enum import Enum

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

        # PDG ID values for quarks and gluons
        self.quark_pdg_ids = [1, 2, 3, 4, 5, 6, 7, 8]
        self.down_pdg_ids = [1]
        self.up_pdg_ids = [2]
        self.strange_pdg_ids = [3]
        self.charm_pdg_ids = [4]
        self.gluon_pdg_ids = [9, 21]

        # hadron level - ALICE tracking restriction
        self.max_eta_hadron = 0.9

        self.min_leading_track_pT = config["min_leading_track_pT"] if "min_leading_track_pT" in config else None

        self.pt_bins = array.array('d', list(range(5, 100, 5)) + list(range(100, 210, 10)))
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

            obs_subconfig_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]
            pinfo("obs_subconfig_list", obs_subconfig_list)
            self.obs_settings[observable] = self.utils.obs_settings(observable, obs_config_dict, obs_subconfig_list)
            pinfo("self.obs_settings[observable]", self.obs_settings[observable])
            self.obs_grooming_settings[observable] = self.utils.grooming_settings(obs_config_dict)
            pinfo("self.obs_grooming_settings[observable]", self.obs_grooming_settings[observable])

            self.obs_names[observable] = obs_config_dict["common_settings"]["xtitle"]

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

        for jetR in self.jetR_list:

            # Store a list of all the histograms just so that we can rescale them later
            hist_list_name = "hist_list_R%s" % str(jetR).replace('.', '')
            setattr(self, hist_list_name, [])

            R_label = str(jetR).replace('.', '') + 'Scaled'

            for observable in self.observable_list:
                # Should only be one: observable == "EEC"
                if observable != "EEC":
                    raise ValueError("Observable %s is not implemented in this script" % observable)

                obs_name = self.obs_names[observable]
                obs_bins = getattr(self, "obs_bins_" + observable)
                # Use more finely binned pT bins for TH2s than for the RMs
                pt_bins = array.array('d', list(range(0, 201, 1)))

                # Loop over subobservable (alpha value)
#                for i in range(len(self.obs_settings[observable])):

                obs_setting = self.obs_settings[observable]
                grooming_setting = self.obs_grooming_settings[observable]
                obs_label = self.utils.obs_label(obs_setting, grooming_setting)
                pinfo("all the settings", obs_setting, grooming_setting, obs_label)

                for parton_type in ["charm", "light", "gluon", "inclusive"]: #got rid of quark

                    name = ('h_%s_JetPt_%s_R%s_%s' % (observable, parton_type, jetR, obs_label)) if \
                        len(obs_label) else ('h_%s_JetPt_%s_R%s' % (observable, parton_type, jetR))
                    h = ROOT.TH2F(name, name, len(pt_bins)-1, pt_bins, len(obs_bins)-1, obs_bins)
                    h.GetXaxis().SetTitle('#it{p}_{T,%s}^{ch jet}' % (parton_type[0] + "-init"))
#                    h.GetYaxis().SetTitle(obs_name + '^{%s}' % (parton_type[0] + "-init"))
                    h.GetYaxis().SetTitle("R_{L}" + '^{%s}' % (parton_type[0] + "-init"))
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    # make another of histogram for the unweighted EEC
                    name = ('h_%s_JetPt_%s_R%s_%s_unweighted' % (observable, parton_type, jetR, obs_label)) if \
                        len(obs_label) else ('h_%s_JetPt_%s_R%s_unweighted' % (observable, parton_type, jetR))
                    h = ROOT.TH2F(name, name, len(pt_bins)-1, pt_bins, len(obs_bins)-1, obs_bins)
                    h.GetXaxis().SetTitle('#it{p}_{T,%s}^{ch jet}' % (parton_type[0] + "-init"))
#                    h.GetYaxis().SetTitle(obs_name + '^{%s}' % (parton_type[0] + "-init"))
                    h.GetYaxis().SetTitle("R_{L}" + '^{%s}' % (parton_type[0] + "-init"))
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)


                    # make another of histogram for the jet level (above is pair level)
                    name_jetpt = ('h_JetPt_%s_R%s_%s_jetlevel' % (parton_type, jetR, obs_label)) if \
                        len(obs_label) else ('h_JetPt_%s_R%s_jetlevel' % (parton_type, jetR))
                    h_jetpt = ROOT.TH1F(name_jetpt, name_jetpt, len(pt_bins)-1, pt_bins)
                    h_jetpt.GetXaxis().SetTitle('#it{p}_{T,%s}^{ch jet}' % (parton_type[0] + "-init"))
                    h_jetpt.GetYaxis().SetTitle('Counts')
                    h_jetpt.Sumw2()
                    setattr(self, name_jetpt, h_jetpt)
                    getattr(self, hist_list_name).append(h_jetpt)

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

            self.parents = []
            self.event = pythia.event
            # print(self.event) # to print out a table of the event information
            fs_parton_5 = fj.PseudoJet(pythia.event[5].px(), pythia.event[5].py(), pythia.event[5].pz(), pythia.event[5].e())
            fs_parton_6 = fj.PseudoJet(pythia.event[6].px(), pythia.event[6].py(), pythia.event[6].pz(), pythia.event[6].e())
            self.parents = [fs_parton_5, fs_parton_6] # parent partons in dijet

            # Save PDG code of the parent partons
            self.parent_ids = [pythia.event[5].id(), pythia.event[6].id()]

            # Print out some information
            # print("NEW EVENT", iev)
            # print("pythia.event size is ", pythia.event.size())
            # eventcounter = 0
            # for event in pythia.event:
                # '''
                # if eventcounter >=3 and eventcounter < 10:
                #     print(eventcounter, "parton! with event id", event.id())
                # if event.id() == 111 or event.id() == 211 or event.id() == -211: #pi0 or pi+ or pi-
                #     print(eventcounter, "pion with event id", event.id())
                # elif event.id() == 321 or event.id() == -321:
                #     print(eventcounter, "kaon with event id", event.id())
                # elif not(event.id() == 21 or (event.id() >=-8 and event.id() <= 8)): #gluon or quark
                #     print(eventcounter, "not q/g with event id", event.id())
                #     '''
                # print(event.id())
                # eventcounter+=1
            # print("event counter is ", eventcounter)


            # parton level
            #parts_pythia_p = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)

            hstatus = pythia.forceHadronLevel()
            if not hstatus:
                continue

            # full-hadron level
            #parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)

            # print("!! pythia hadron (before vectorization) event size is ", pythia.event.size())
            # eventcounter = 0
            # for event in pythia.event:
            #     # if event.id() == 111 or event.id() == 211 or event.id() == -211: #pi0 or pi+ or pi-
            #         # print(eventcounter, "pion with event id", event.id())
            #     eventcounter+=1


            # charged-hadron level
            if ( self.replaceKPpairs == False ):
                parts_pythia_hch = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged], 0, True)
            else: #replace D0->Kpi
                parts_pythia_hch = pythiafjext.vectorize_select_replaceD0(pythia, [pythiafjext.kFinal, pythiafjext.kCharged], 0, True)
            # print("Size of new vector", len(parts_pythia_hch))

            # look at events in charged hadron??
            # print("!! pythia hadron (after vectorization) event size is ", pythia.event.size())
            particlecounter = 0
            # eventcounter1 = 0
            # eventcounter2 = 0
            # indeximportant = -1
            D0found = False
            D0Kpidecayfound = False
            for particle in self.event: #for event in pythia.event:
            #     # if particle.id() == 111 or particle.id() == 211 or particle.id() == -211: #pi0 or pi+ or pi-
            #         # print(particlecounter, "pion with particle id", particle.id())
                if particle.id() == 421 or particle.id() == -421: #D0
                    D0found = True
                #     print(particlecounter, "D0 with particle id", particle.id())
                    if self.checkDecayChannel(particle, self.event) == EMesonDecayChannel.kDecayD0toKpi:
                        print(particlecounter, "This is a D0->Kpi decay!", particle.id())
                        print("Size of new vector", len(parts_pythia_hch))
                        D0Kpidecayfound = True

                particlecounter+=1
            #         # print("D0 daughter indices", pythia.event[event.daughter1()].id(), pythia.event[event.daughter2()].id())
            #     if event.isCharged() and event.isFinal():
            #         eventcounter1+=1
            #         if indeximportant == -1:
            #             indeximportant = eventcounter
            #     # if event.isNeutral():
            #     else:
            #         eventcounter2+=1
            #         # print(eventcounter2, "neutral hadron with event id", event.id())
            #     eventcounter+=1
            # print("eventcounter1", eventcounter1)
            # # print(eventcounter1+eventcounter2)

            #if D0->Kpi found, count the events; if not, check that length of charged final state hadrons vector is 0
            if (D0Kpidecayfound):
                self.hD0KpiNevents.Fill(0)
            if (D0found):
                self.hD0Nevents.Fill(0)
            # else:
                # if len(parts_pythia_hch) > 0:
                #         print(particlecounter, "There ARE particles in this jet (but shouldn't be)", len(parts_pythia_hch))

            
            


            # self.event = pythia.event
            # iev_to_print = [28, 31, 32, 58, 61, 89, 110, 113, 118, 123, 137, 144, 147, 150, 161, 173, 174, 175, 185, 194, 200, 202,210,211,216, 228, 231, 236, 240, 247, 248,251, 257, 259,263, 272, 273,284, 292,301,972]
            # if (iev in iev_to_print):
            #     print("event printed!")
            #     print(self.event)


            # if there is no D0 for charm jet, remove
            # if ( self.replaceKPpairs ):
            #     D0_indices = list(filter(lambda event: (event.idAbs() == 421), pythia.event))
            #     print("There are", len(D0_indices), "D0s in this list")
            #     if (len(D0_indices) == 0):
            #         parts_pythia_hch = []


            # print("!! length (parts_pythi_hch)", len(parts_pythia_hch))
            # print("ID COMPARISON", pythia.event[indeximportant].id(), parts_pythia_hch[0].description()) #user_index())


                

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

    #---------------------------------------------------------------
    # Find jets, do matching between levels, and fill histograms
    #---------------------------------------------------------------
    def find_jets_fill_histograms(self, parts_pythia_hch, iev, D0Kpidecayfound):
        # Loop over jet radii
        for jetR in self.jetR_list:

            jetR_str = str(jetR).replace('.', '')
            jet_selector = getattr(self, "jet_selector_R%s" % jetR_str)
            jet_def = getattr(self, "jet_def_R%s" % jetR_str)

            count1 = getattr(self, "count1_R%s" % jetR_str)
            count2 = getattr(self, "count2_R%s" % jetR_str)

            # Get the jets at different levels
            #jets_p  = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_p  ))) # parton level
            #jets_h  = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_h  ))) # full hadron level
            jets_ch = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_hch))) # charged hadron level
            # print("!! length of jets_ch", len(jets_ch))

            R_label = str(jetR).replace('.', '') + 'Scaled'

            # look at events in charged hadron??
            # print("pythia charged hadron event size is ", pythia.event.size())


            # # Replace k/pi pairs with D0
            # if (self.replaceKPpairs):
            #     print("hi")
            anothacounter=0

            # Find the charged jet closest to the axis of the original parton
            # Require that the match is within some small angle, and that it is unique
            jet_matching_distance = 0.6  # Match jets with deltaR < jet_matching_distance*jetR
            self.parent0match, self.parent1match = None, None
            # print("LOOPING OVER JETS")
            for i_jch, jch in enumerate(jets_ch):
                # print(i_jch)
                # Do constituent pT cut
#                pinfo("self.min_leading_track_pT", self.min_leading_track_pT)
#                 if self.min_leading_track_pT and not \
#                     self.utils.is_truth_jet_accepted(jch):
# #                   self.utils.is_truth_jet_accepted(jch, self.min_leading_track_pT):
#                     continue
                # print("PARENTS:",self.parents)
                for i_parent, parent in enumerate(self.parents):
                    anothacounter+=1
                    parentmatch_name = "parent%imatch" % i_parent
                    print("CHECKING PARENT", i_parent)
                    print("DELTA R TO JET:", jch.delta_R(parent))
                    if jch.delta_R(parent) < jet_matching_distance * jetR:
                        match = getattr(self, parentmatch_name)
                        print("MATCH FOR",i_parent,":",match)
                        if not match:
                            setattr(self, parentmatch_name, jch)
                            print("MATCH SET TO JET WITH pT", jch.pt())
                        else:  # Already found a match
                            # Set flag value so that we know to ignore this one
                            print("already found a match flagged")
                            setattr(self, parentmatch_name, 0)
                    # print(i_jch, "anothacounter", anothacounter)

            # print("event num", iev)
            # print("Cehckpoint 1")

            # If we have matches, fill histograms
            for i_parent, parent in enumerate(self.parents):
                # print("in nexr loop")
                jet = getattr(self, "parent%imatch" % i_parent)
                # print(jet)
                if not jet:
                    # pinfo("in not jet")
                    if jet == 0: # More than one match -- take note and continue
                        # print("case  1")
                        count1 += 1
                        continue
                    else:  # jet == None
                        # No matches -- take note and continue
                        # print("case  2")
                        count2 += 1
                        continue

                # print("CHECKPOINT 2")

                # One unique match
                # Identify the histograms which need to be filled
#                pinfo("passed not jet")
                parton_id = self.parent_ids[i_parent]
                # print("parton_id is ", parton_id)
                parton_types = []
                if parton_id in self.quark_pdg_ids:
                    # parton_types += ["quark"]
                    if parton_id in self.charm_pdg_ids:
                        parton_types += ["charm"]
                    elif parton_id in self.up_pdg_ids or parton_id in self.down_pdg_ids or parton_id in self.strange_pdg_ids:
                        parton_types += ["light"]
                elif parton_id in self.gluon_pdg_ids:
                    parton_types += ["gluon"]
                parton_types += ["inclusive"]

                # If parent parton not identified, skip for now
                if not len(parton_types):
                    continue
                
                print(D0Kpidecayfound)
                if D0Kpidecayfound:
                    print("parton types")

                # Fill histograms
                for observable in self.observable_list:
#                    pinfo("len(self.obs_settings[observable])", len(self.obs_settings[observable]))
#                    for i in range(len(self.obs_settings[observable])): #looping through different configurations?

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

#                        obs = self.calculate_observable(
#                            observable, jet, jet_groomed_lund, jetR, obs_setting,
#                            grooming_setting, obs_label, jet.pt())
                    obs = self.calculate_observable(
                        observable, jet, jet_groomed_lund, jetR, jet.pt())
                    

#                        for parton_type in parton_types:
#                            getattr(self, ('h_%s_JetPt_%s_R%s_%s' % (observable, parton_type, jetR, obs_label)) if \
#                                len(obs_label) else ('h_%s_JetPt_%s_R%s' % (observable, parton_type, jetR))).Fill(
#                                jet.pt(), obs)
                    for parton_type in parton_types:
                            getattr(self, ('h_JetPt_%s_R%s_%s_jetlevel' % (parton_type, jetR, obs_label)) if \
                                len(obs_label) else ('h_JetPt_%s_R%s_jetlevel' % (parton_type, jetR))).Fill(jet.pt())
                    for index in range(obs.correlator(2).rs().size()):
                        for parton_type in parton_types:
                            # if (self.weighted == True):
                            getattr(self, ('h_%s_JetPt_%s_R%s_%s' % (observable, parton_type, jetR, obs_label)) if \
                                len(obs_label) else ('h_%s_JetPt_%s_R%s' % (observable, parton_type, jetR))).Fill(jet.pt(), obs.correlator(2).rs()[index], obs.correlator(2).weights()[index])
                            # getattr(self, ('h_%s_JetPt_%s_R%s_%s_unweighted' % (observable, parton_type, jetR, obs_label)) if \
                            #     len(obs_label) else ('h_%s_JetPt_%s_R%s_unweighted' % (observable, parton_type, jetR))).Fill(jet.pt(), obs.correlator(2).rs()[index])

            setattr(self, "count1_R%s" % jetR_str, count1)
            setattr(self, "count2_R%s" % jetR_str, count2)

    #---------------------------------------------------------------
    # Calculate the observable given a jet
    #---------------------------------------------------------------
#    def calculate_observable(self, observable, jet, jet_groomed_lund,
#        jetR, obs_setting, grooming_setting, obs_label, jet_pt_ungroomed):
    def calculate_observable(self, observable, jet, jet_groomed_lund,
        jetR, jet_pt_ungroomed):

        if observable == "EEC":

            # idk about this grooming setting....
            # WRITE STUFF HERE
            
            #TODO: get the dcand from the jet! (need the djmm??)
#            j = jet[0]
#            dcand = djmm.get_Dcand_in_jet(jet)
            
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
            
#            new_corr = ecorrel.CorrelatorBuilder(c_select, dcand, jet_pt_ungroomed, 2, 1, dphi_cut, deta_cut) #jet, D, scale, max, power, dphicut, detacut
            new_corr = ecorrel.CorrelatorBuilder(c_select, jet.perp(), 2, 1, dphi_cut, deta_cut)
            
            #print out weights here
            # for index in range(new_corr.correlator(2).rs().size()):
            #     print("weight", new_corr.correlator(2).weights()[index])

            return new_corr #new_corr.correlator(2).rs()[index]
            
#            return fjext.lambda_beta_kappa(jet, jet_groomed_lund.pair(), obs_setting, 1, jetR) \
#                   if grooming_setting else fjext.lambda_beta_kappa(jet, obs_setting, 1, jetR)

#        elif observable == "mass":
#
#            if grooming_setting:
#                j_groomed = jet_groomed_lund.pair()
#                if not j_groomed.has_constituents():
#                    # Untagged jet -- record underflow value
#                    return -1
#                else:
#                    return j_groomed.m()
#
#            return jet.m()

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
    

            absPdg1 = d1.idAbs()
            absPdg2 = d2.idAbs()

            if(absPdgPart == 421):  # D0 -> K pi
                if((absPdg1 == 211 and absPdg2 == 321) or (absPdg1 == 321 and absPdg2 == 211)): # pi K or K pi - QUESTION: does this account for k and pi being opposite signs?
                    decay = EMesonDecayChannel.kDecayD0toKpi 
      
            if(absPdgPart == 413):  # D* -> D0 pi
                if(absPdg1 == 421 and absPdg2 == 211):   # D0 pi
                    D0decay = self.checkDecayChannel(d1, event)
                    if(D0decay == EMesonDecayChannel.kDecayD0toKpi):
                        decay = EMesonDecayChannel.kDecayDStartoKpipi
                elif(absPdg1 == 211 and absPdg2 == 421):   # pi D0
                    D0decay = self.checkDecayChannel(d2, event)
                    if(D0decay == EMesonDecayChannel.kDecayD0toKpi):
                        decay = EMesonDecayChannel.kDecayDStartoKpipi

        return decay
    
    
    
    

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

    args = parser.parse_args()

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
