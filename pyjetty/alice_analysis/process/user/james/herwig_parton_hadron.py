#!/usr/bin/env python

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext

import ROOT

import tqdm
import yaml
import copy
import argparse
import os

from pyjetty.mputils import *

from pyjetty.alice_analysis.process.base import process_base
from pyjetty.alice_analysis.process.user.ang_pp.helpers import lambda_alpha_kappa_i

from array import array
import numpy as np

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)
# Automatically set Sumw2 when creating new histograms
ROOT.TH1.SetDefaultSumw2()

################################################################
class herwig_parton_hadron(process_base.ProcessBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, input_file='', config_file='', output_dir='',
                 debug_level=0, args=None, **kwargs):

        super(herwig_parton_hadron, self).__init__(
            input_file, config_file, output_dir, debug_level, **kwargs)

        self.initialize_config(args)


    #---------------------------------------------------------------
    # Initialize config file into class members
    #---------------------------------------------------------------
    def initialize_config(self, args):

        # Call base class initialization
        process_base.ProcessBase.initialize_config(self)

        # Read config file
        with open(self.config_file, 'r') as stream:
            config = yaml.safe_load(stream)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        self.herwig_file = args.input_file
        self.herwig_file_MPI = args.input_file_mpi

        # Defaults to None if not in use
        self.level = args.no_match_level

        self.jetR_list = config["jetR"]

        # Formatted LaTeX names for plotting
        self.obs_names = ["#it{#theta}_{g}", "#it{z}_{g}"]
        self.observables = config['process_observables']
        self.obs_settings = {}
        self.obs_grooming_settings = {}
        for observable in self.observables:
            obs_config_dict = config[observable]
            obs_config_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]

            obs_subconfig_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]
            self.obs_settings[observable] = self.utils.obs_settings(
                observable, obs_config_dict, obs_subconfig_list)
            self.obs_grooming_settings[observable] = self.utils.grooming_settings(obs_config_dict)

        # Construct set of unique grooming settings
        self.grooming_settings = []
        lists_grooming = [self.obs_grooming_settings[obs] for obs in self.observables]
        for observable in lists_grooming:
            for setting in observable:
                if setting not in self.grooming_settings and setting != None:
                    self.grooming_settings.append(setting)
        self.grooming_labels = [self.utils.grooming_label(gs) for gs in self.grooming_settings]

        # Observable binnings for theta_g and zg
        self.obs_bins_theta_g = array('d', np.arange(0, 1.01, 0.01))
        self.obs_bins_zg = array('d', np.arange(0, 0.505, 0.005))

        # We are not reporting zg theory so save time/memory by skipping these histograms
        self.skip_zg = True
        if self.skip_zg:
            self.obs_names = self.obs_names[:-1]
            self.observables = self.observables[:-1]

        # Manually added binnings for RM and scaling histograms
        if 'theory_pt_bins' in config:
            self.pt_bins = array('d', config['theory_pt_bins'])

        # hadron level - ALICE tracking restriction
        self.max_eta_hadron = 0.9

        # Whether or not to rescale final jet histograms based on sigma/N
        self.no_scale = args.no_scale

        # Whether or not to save particle info in raw tree structure
        self.no_tree = args.no_tree

        # Initialize variables for final cross sections from event generator
        self.xsec = None
        self.xsec_MPI = None


    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def herwig_parton_hadron(self, args):

        # Create ROOT TTree file for storing raw PYTHIA particle information
        outf_path = os.path.join(self.output_dir, args.tree_output_fname)
        outf = ROOT.TFile(outf_path, 'recreate')
        outf.cd()

        # Initialize response histograms
        self.initialize_hist()

        # Print the banner first
        fj.ClusterSequence.print_banner()
        print()

        self.init_jet_tools()
        self.parse_events()
        if self.herwig_file_MPI:
            self.parse_events(MPIon=True)

        if not self.no_tree:
            for jetR in self.jetR_list:
                getattr(self, "tw_R%s" % str(jetR).replace('.', '')).fill_tree()

        # Scale histograms
        self.scale_print_final_info()

        outf.Write()
        outf.Close()

        self.save_output_objects()


    #---------------------------------------------------------------
    # Initialize histograms
    #---------------------------------------------------------------
    def initialize_hist(self):

        self.hNevents = ROOT.TH1I("hNevents", 'Number accepted events (unscaled)', 2, -0.5, 1.5)
        self.hNeventsMPI = ROOT.TH1I("hNeventsMPI", 'Number accepted events (unscaled)', 2, -0.5, 1.5)

        for jetR in self.jetR_list:

            # Store a list of all the histograms just so that we can rescale them later
            hist_list_name = "hist_list_R%s" % str(jetR).replace('.', '')
            setattr(self, hist_list_name, [])
            hist_list_name_MPIon = "hist_list_MPIon_R%s" % str(jetR).replace('.', '')
            setattr(self, hist_list_name_MPIon, [])

            R_label = str(jetR) + 'Scaled'

            if self.level in [None, 'ch']:
                name = 'hJetPt_ch_R%s' % R_label
                h = ROOT.TH1F(name, name+';p_{T}^{ch jet};#frac{dN}{dp_{T}^{ch jet}};', 300, 0, 300)
                h.Sumw2()  # enables calculation of errors
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

                name = 'hNconstit_Pt_ch_R%s' % R_label
                h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                h.GetXaxis().SetTitle('#it{p}_{T}^{ch jet}')
                h.GetYaxis().SetTitle('#it{N}_{constit}^{ch jet}')
                h.Sumw2()
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

            if self.level in [None, 'h']:
                name = 'hJetPt_h_R%s' % R_label
                h = ROOT.TH1F(name, name+';p_{T}^{jet, h};#frac{dN}{dp_{T}^{jet, h}};', 300, 0, 300)
                h.Sumw2()
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

                name = 'hNconstit_Pt_h_R%s' % R_label
                h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                h.GetXaxis().SetTitle('#it{p}_{T}^{h jet}')
                h.GetYaxis().SetTitle('#it{N}_{constit}^{h jet}')
                h.Sumw2()
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

            if self.level in [None, 'p']:
                name = 'hJetPt_p_R%s' % R_label
                h = ROOT.TH1F(name, name+';p_{T}^{jet, parton};#frac{dN}{dp_{T}^{jet, parton}};',
                              300, 0, 300)
                h.Sumw2()
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

                name = 'hNconstit_Pt_p_R%s' % R_label
                h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                h.GetXaxis().SetTitle('#it{p}_{T}^{p jet}')
                h.GetYaxis().SetTitle('#it{N}_{constit}^{p jet}')
                h.Sumw2()
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

            if self.level == None:
                name = 'hJetPtRes_R%s' % R_label
                h = ROOT.TH2F(name, name, 300, 0, 300, 200, -1., 1.)
                h.GetXaxis().SetTitle('#it{p}_{T}^{parton jet}')
                h.GetYaxis().SetTitle(
                    '#frac{#it{p}_{T}^{parton jet}-#it{p}_{T}^{ch jet}}{#it{p}_{T}^{parton jet}}')
                h.Sumw2()
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

                name = 'hResponse_JetPt_R%s' % R_label
                h = ROOT.TH2F(name, name, 200, 0, 200, 200, 0, 200)
                h.GetXaxis().SetTitle('#it{p}_{T}^{parton jet}')
                h.GetYaxis().SetTitle('#it{p}_{T}^{ch jet}')
                h.Sumw2()
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

            for i_obs, obs in enumerate(self.observables):
                obs_bins = getattr(self, "obs_bins_%s" % obs)
                obs_name = self.obs_names[i_obs]

                for i, gs in enumerate(self.grooming_settings):
                    gl = self.grooming_labels[i]
                    label = "R%s_%s" % (str(jetR), gl)

                    if self.level in [None, 'ch']:
                        name = 'h_JetPt_%s_ch_MPIoff_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, len(self.pt_bins)-1, self.pt_bins,
                                      len(obs_bins)-1, obs_bins)
                        h.GetXaxis().SetTitle('p_{T}^{ch jet}')
                        h.GetYaxis().SetTitle('#frac{dN}{d%s^{ch}}' % obs_name)
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        name = 'h_JetPt_%s_ch_MPIon_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, len(self.pt_bins)-1, self.pt_bins,
                                      len(obs_bins)-1, obs_bins)
                        h.GetXaxis().SetTitle('p_{T}^{ch jet}')
                        h.GetYaxis().SetTitle('#frac{dN}{d%s^{ch}}' % obs_name)
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name_MPIon).append(h)

                    if self.level in [None, 'h']:
                        name = 'h_JetPt_%s_h_MPIoff_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, len(self.pt_bins)-1, self.pt_bins,
                                  len(obs_bins)-1, obs_bins)
                        h.GetXaxis().SetTitle('p_{T}^{jet, h}')
                        h.GetYaxis().SetTitle('#frac{dN}{d%s^{h}}' % obs_name)
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                    if self.level in [None, 'p']:
                        name = 'h_JetPt_%s_p_MPIoff_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, len(self.pt_bins)-1, self.pt_bins,
                                      len(obs_bins)-1, obs_bins)
                        h.GetXaxis().SetTitle('p_{T}^{jet, parton}')
                        h.GetYaxis().SetTitle('#frac{dN}{d%s^{parton}}' % obs_name)
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                    if self.level == None:
                        name = 'hResponse_%s_p_ch_MPIoff_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, 100, 0, 1, 100, 0, 1)
                        h.GetXaxis().SetTitle('%s^{parton}' % obs_name)
                        h.GetYaxis().SetTitle('%s^{ch}' % obs_name)
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        name = "hResidual_JetPt_%s_MPIoff_%sScaled" % (obs, label)
                        h = ROOT.TH2F(name, name, 300, 0, 300, 200, -3., 1.)
                        h.GetXaxis().SetTitle('p_{T}^{p jet}')
                        h.GetYaxis().SetTitle('#frac{%s^{p}-%s^{ch}}{%s^{p}}' % (obs_name, obs_name, obs_name))
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        name = "hDiff_JetPt_%s_MPIoff_%sScaled" % (obs, label)
                        h = ROOT.TH2F(name, name, 300, 0, 300, 200, -2., 2.)
                        h.GetXaxis().SetTitle('#it{p}_{T}^{ch jet}')
                        h.GetYaxis().SetTitle('%s^{p}-%s^{ch}' % (obs_name, obs_name))
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        # Create THn of response
                        dim = 4
                        title = ['p_{T}^{ch jet}', 'p_{T}^{parton jet}',
                                 obs_name + '^{ch}', obs_name + '^{parton}']
                        nbins  = [len(self.pt_bins)-1, len(self.pt_bins)-1,
                                  len(obs_bins)-1,     len(obs_bins)-1]
                        min_li = [self.pt_bins[0],     self.pt_bins[0],
                                  obs_bins[0],         obs_bins[0]    ]
                        max_li = [self.pt_bins[-1],    self.pt_bins[-1],
                                  obs_bins[-1],        obs_bins[-1]   ]

                        name = 'hResponse_JetPt_%s_p_ch_MPIoff_%sScaled' % (obs, label)
                        nbins = (nbins)
                        xmin = (min_li)
                        xmax = (max_li)
                        nbins_array = array('i', nbins)
                        xmin_array = array('d', xmin)
                        xmax_array = array('d', xmax)
                        h = ROOT.THnF(name, name, dim, nbins_array, xmin_array, xmax_array)
                        for i in range(0, dim):
                            h.GetAxis(i).SetTitle(title[i])
                            if i == 0 or i == 1:
                                h.SetBinEdges(i, self.pt_bins)
                            else:  # i == 2 or i == 3
                                h.SetBinEdges(i, obs_bins)
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        # Another set of THn for full hadron folding
                        title = ['p_{T}^{h jet}', 'p_{T}^{parton jet}',
                                 obs_name + '^{h}', obs_name + '^{parton}']

                        name = 'hResponse_JetPt_%s_p_h_MPIoff_%sScaled' % (obs, label)
                        h = ROOT.THnF(name, name, dim, nbins_array, xmin_array, xmax_array)
                        for i in range(0, dim):
                            h.GetAxis(i).SetTitle(title[i])
                            if i == 0 or i == 1:
                                h.SetBinEdges(i, self.pt_bins)
                            else:  # i == 2 or i == 3
                                h.SetBinEdges(i, obs_bins)
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        # Finally, a set of THn for folding H --> CH (with MPI on)
                        title = ['p_{T}^{ch jet}', 'p_{T}^{h jet}',
                                 obs_name + '^{ch}', obs_name + '^{h}']

                        name = 'hResponse_JetPt_%s_h_ch_MPIon_%sScaled' % (obs, label)
                        h = ROOT.THnF(name, name, dim, nbins_array, xmin_array, xmax_array)
                        for i in range(0, dim):
                            h.GetAxis(i).SetTitle(title[i])
                            if i == 0 or i == 1:
                                h.SetBinEdges(i, self.pt_bins)
                            else:  # i == 2 or i == 3
                                h.SetBinEdges(i, obs_bins)
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name_MPIon).append(h)


    #---------------------------------------------------------------
    # Initiate jet defs, selectors, and sd (if required)
    #---------------------------------------------------------------
    def init_jet_tools(self):

        for jetR in self.jetR_list:
            jetR_str = str(jetR).replace('.', '')

            if not self.no_tree:
                # Initialize tree writer
                name = 'particle_unscaled_R%s' % jetR_str
                t = ROOT.TTree(name, name)
                setattr(self, "t_R%s" % jetR_str, t)
                tw = RTreeWriter(tree=t)
                setattr(self, "tw_R%s" % jetR_str, tw)

            # set up our jet definition and a jet selector
            jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
            setattr(self, "jet_def_R%s" % jetR_str, jet_def)
            print(jet_def)

        pwarning('max eta for particles after hadronization set to', self.max_eta_hadron)
        parts_selector_h = fj.SelectorAbsEtaMax(self.max_eta_hadron)

        for jetR in self.jetR_list:
            jetR_str = str(jetR).replace('.', '')

            jet_selector = fj.SelectorPtMin(5.0) & \
                           fj.SelectorAbsEtaMax(self.max_eta_hadron - jetR)
            setattr(self, "jet_selector_R%s" % jetR_str, jet_selector)

            #max_eta_parton = self.max_eta_hadron + 2. * jetR
            #setattr(self, "max_eta_parton_R%s" % jetR_str, max_eta_parton)
            #pwarning("Max eta for partons with jet R =", jetR, "set to", max_eta_parton)
            #parts_selector_p = fj.SelectorAbsEtaMax(max_eta_parton)
            #setattr(self, "parts_selector_p_R%s" % jetR_str, parts_selector_p)

            count1 = 0  # Number of jets rejected from ch-h matching
            setattr(self, "count1_R%s" % jetR_str, count1)
            count2 = 0  # Number of jets rejected from h-p matching
            setattr(self, "count2_R%s" % jetR_str, count2)

    #---------------------------------------------------------------
    # Read events from output, find jets, and fill histograms
    #---------------------------------------------------------------
    def parse_events(self, MPIon=False):

        if MPIon:
            hNevents = self.hNeventsMPI
            infile = self.herwig_file_MPI
        else:
            hNevents = self.hNevents
            infile = self.herwig_file

        print("Reading events from %s..." % infile)

        with open(infile, 'r') as f:
            ev_num = 0

            # Flags to assist with keeping track of place within file
            reading_ev = False
            parton = False
            parton_final = False
            parton_finished = False
            hadron = False
            hadron_final = False

            partons_px = []
            partons_py = []
            partons_pz = []
            partons_e = []
            #partons_q = []
            hadrons_px = []
            hadrons_py = []
            hadrons_pz = []
            hadrons_e = []
            #hadrons_q = []
            ch_hadrons_px = []
            ch_hadrons_py = []
            ch_hadrons_pz = []
            ch_hadrons_e = []
            #ch_hadrons_q = []

            for line in f:

                # Waiting to start reading event
                if not reading_ev:
                    if "Event number" in line:
                        reading_ev = True
                        ev_num = int(line.split()[2])
                        if not ev_num % 1000:
                            print("Event number", ev_num, end="\r")
                        hNevents.Fill(0)
                    elif "Total integrated xsec:" in line:
                        if MPIon:
                            self.xsec_MPI = float(line.split()[3])
                        else:
                            self.xsec = float(line.split()[3])
                    continue

                # Reading event
                # First step is to read the parton info
                elif not parton:
                    if "ShowerHandler" in line:
                        parton = True
                    continue
                elif not parton_final:
                    if "final" in line:
                        parton_final = True
                    continue

                # Get showered partons
                elif not MPIon and not parton_finished:
                    # Read parton information
                    vals = line.split()
                    if line[0] == '-':
                        parton_finished = True
                    elif len(vals) == 5 and line[2] == ' ':
                        partons_px.append(vals[0])
                        partons_py.append(vals[1])
                        partons_pz.append(vals[2])
                        partons_e.append(vals[3])
                        #partons_q.append(vals[4])
                    continue

                # Get final hadrons
                elif not hadron:
                    if "DecayHandler" in line:
                        hadron = True
                    continue
                elif not hadron_final:
                    if "final" in line:
                        hadron_final = True
                    continue

                # Check if event is over
                elif line[0] == '-':
                    # Finished reading hadron info
                    reading_ev = False
                    parton = False
                    parton_final = False
                    parton_finished = False
                    hadron = False
                    hadron_final = False

                    # Get correct structure for finding jets
                    partons = None
                    hadrons = None
                    if not MPIon:
                        partons = fjext.vectorize_px_py_pz_e(
                            partons_px, partons_py, partons_pz, partons_e)

                        hadrons = fjext.vectorize_px_py_pz_e(
                            hadrons_px, hadrons_py, hadrons_pz, hadrons_e)

                    ch_hadrons = fjext.vectorize_px_py_pz_e(
                        ch_hadrons_px, ch_hadrons_py, ch_hadrons_pz, ch_hadrons_e)

                    self.find_jets_fill_hist(partons, hadrons, ch_hadrons, ev_num, MPIon)

                    partons_px = []
                    partons_py = []
                    partons_pz = []
                    partons_e = []
                    #partons_q = []
                    hadrons_px = []
                    hadrons_py = []
                    hadrons_pz = []
                    hadrons_e = []
                    #hadrons_q = []
                    ch_hadrons_px = []
                    ch_hadrons_py = []
                    ch_hadrons_pz = []
                    ch_hadrons_e = []
                    #ch_hadrons_q = []

                    continue

                elif line[2].isnumeric():
                    # Save the hadron name for charge descrimination
                    i = 1
                    while line[i].isnumeric():
                        i += 1
                    hadron_type = line[i:].split()[0]
                    continue

                elif line[2] == ' ':  # and len(line.split()) == 5:
                    # Reading hadron information
                    vals = line.split()
                    if not MPIon:
                        hadrons_px.append(vals[0])
                        hadrons_py.append(vals[1])
                        hadrons_pz.append(vals[2])
                        hadrons_e.append(vals[3])
                        #hadrons_q.append(vals[4])

                    if '+' in hadron_type or '-' in hadron_type:
                        ch_hadrons_px.append(vals[0])
                        ch_hadrons_py.append(vals[1])
                        ch_hadrons_pz.append(vals[2])
                        ch_hadrons_e.append(vals[3])
                        #ch_hadrons_q.append(vals[4])

        return partons, hadrons, ch_hadrons

    #---------------------------------------------------------------
    # Read events from output, find jets, and fill histograms
    #---------------------------------------------------------------
    def find_jets_fill_hist(self, partons, hadrons,
                            ch_hadrons, iev, MPIon=False):

        for jetR in self.jetR_list:
            #print("Filling jet histograms for R = %s..." % str(jetR))

            jetR_str = str(jetR).replace('.', '')
            jet_selector = getattr(self, "jet_selector_R%s" % jetR_str)
            jet_def = getattr(self, "jet_def_R%s" % jetR_str)
            t = None; tw = None;
            if not self.no_tree:
                t = getattr(self, "t_R%s" % jetR_str)
                tw = getattr(self, "tw_R%s" % jetR_str)
            count1 = getattr(self, "count1_R%s" % jetR_str)
            count2 = getattr(self, "count2_R%s" % jetR_str)

            #if not (iev+1) % 1000:
            #    print("Event number %s" % str(iev+1), end='\r')

            jets_ch = fj.sorted_by_pt(jet_selector(jet_def(ch_hadrons)))
            jets_p = fj.sorted_by_pt(jet_selector(jet_def(partons)))
            jets_h = fj.sorted_by_pt(jet_selector(jet_def(hadrons)))

            if MPIon:
                for jet in jets_ch:
                    self.fill_MPI_histograms(jetR, jet)

            if self.level and not MPIon:  # Only save info at one level w/o matching
                if not self.no_tree:
                    jets = locals()["jets_%s" % self.level]
                    for jet in jets:
                        self.fill_unmatched_jet_tree(tw, jetR, iev, jet)
                continue

            for i,jchh in enumerate(jets_ch):

                # match hadron (full) jet
                drhh_list = []
                for j, jh in enumerate(jets_h):
                    drhh = jchh.delta_R(jh)
                    if drhh < jetR / 2.:
                        drhh_list.append((j,jh))
                if len(drhh_list) != 1:
                    count1 += 1
                else:  # Require unique match
                    j, jh = drhh_list[0]

                    # match parton level jet
                    dr_list = []
                    for k, jp in enumerate(jets_p):
                        dr = jh.delta_R(jp)
                        if dr < jetR / 2.:
                            dr_list.append((k, jp))
                    if len(dr_list) != 1:
                        count2 += 1
                    else:
                        k, jp = dr_list[0]

                        if self.debug_level > 0:
                            pwarning('event', iev)
                            pinfo('matched jets: ch.h:', jchh.pt(), 'h:', jh.pt(),
                                  'p:', jp.pt(), 'dr:', dr)

                        if not MPIon:
                            self.fill_jet_histograms(jetR, jp, jh, jchh)
                            if not self.no_tree:
                                self.fill_matched_jet_tree(tw, jetR, iev, jp, jh, jchh)
                        else:
                            self.fill_jet_histograms_MPI(jetR, jp, jh, jchh)

                #print("  |-> SD jet params z={0:10.3f} dR={1:10.3f} mu={2:10.3f}".format(
                #    sd_info.z, sd_info.dR, sd_info.mu))

            if MPIon:
                setattr(self, "count1_R%s_MPIon" % jetR_str, count1)
                setattr(self, "count2_R%s_MPIon" % jetR_str, count2)
            else:
                setattr(self, "count1_R%s" % jetR_str, count1)
                setattr(self, "count2_R%s" % jetR_str, count2)


    #---------------------------------------------------------------
    # Fill jet tree with (unscaled/raw) matched parton/hadron tracks
    #---------------------------------------------------------------
    def fill_matched_jet_tree(self, tw, jetR, iev, jp, jh, jchh):

        tw.fill_branch('iev', iev)
        tw.fill_branch('ch', jchh)
        tw.fill_branch('h', jh)
        tw.fill_branch('p', jp)

        for i, gs in enumerate(self.grooming_settings):
            gl = self.grooming_labels[i]

            # Groomed jets
            gshop_chh = fjcontrib.GroomerShop(jchh, jetR, self.reclustering_algorithm)
            jet_ch_groomed_lund = self.utils.groom(gshop_chh, gs, jetR)
            if not jet_ch_groomed_lund:
                continue

            gshop_h = fjcontrib.GroomerShop(jh, jetR, self.reclustering_algorithm)
            jet_h_groomed_lund = self.utils.groom(gshop_h, gs, jetR)
            if not jet_h_groomed_lund:
                continue

            gshop_p = fjcontrib.GroomerShop(jp, jetR, self.reclustering_algorithm)
            jet_p_groomed_lund = self.utils.groom(gshop_p, gs, jetR)
            if not jet_p_groomed_lund:
                continue

            obs_dict = None
            for obs in self.observables:
                if obs == "theta_g":
                    obs_dict = {
                        "p" : jet_p_groomed_lund.Delta() / jetR,
                        "h" : jet_h_groomed_lund.Delta() / jetR,
                        "ch": jet_ch_groomed_lund.Delta() / jetR }
                elif obs == "zg":
                    if self.skip_zg:
                        continue
                    obs_dict = {
                        "p" : jet_p_groomed_lund.z(),
                        "h" : jet_h_groomed_lund.z(),
                        "ch": jet_ch_groomed_lund.z() }
                else:
                    raise ValueError("Unrecognized observable " + obs)

                for level in ["p", "h", "ch"]:
                    tw.fill_branch("%s_%s_%s" % (obs, level, gl), obs_dict[level])


    #---------------------------------------------------------------
    # Fill jet tree with (unscaled/raw) unmatched parton/hadron tracks
    #---------------------------------------------------------------
    def fill_unmatched_jet_tree(self, tw, jetR, iev, jet):

        tw.fill_branch('iev', iev)
        tw.fill_branch(self.level, jet)

        for i, gs in enumerate(self.grooming_settings):
            gl = self.grooming_labels[i]

            # Groomed jet
            gshop = fjcontrib.GroomerShop(jet, jetR, self.reclustering_algorithm)
            jet_groomed_lund = self.utils.groom(gshop, gs, jetR)
            if not jet_groomed_lund:
                continue

            obs_val = None
            for obs in self.observables:
                if obs == "theta_g":
                    obs_val = jet_p_groomed_lund.Delta() / jetR
                elif obs == "zg":
                    if self.skip_zg:
                        continue
                    obs_val = jet_p_groomed_lund.z()

                tw.fill_branch("%s_%s_%s" % (obs, self.level, gl), obs_val)


    #---------------------------------------------------------------
    # Fill jet histograms for MPI-on PYTHIA run-through
    #---------------------------------------------------------------
    def fill_MPI_histograms(self, jetR, jet):

        for i, gs in enumerate(self.grooming_settings):
            gl = self.grooming_labels[i]
            label = "R" + str(jetR) + '_' + gl

            # Groomed jet
            gshop = fjcontrib.GroomerShop(jet, jetR, self.reclustering_algorithm)
            jet_groomed_lund = self.utils.groom(gshop, gs, jetR)
            if not jet_groomed_lund:
                continue

            obs_val = None
            for obs in self.observables:
                if obs == "theta_g":
                    obs_val = jet_groomed_lund.Delta() / jetR
                elif obs == "zg":
                    if self.skip_zg:
                        continue
                    obs_val = jet_groomed_lund.z()
                else:
                    raise ValueError("Unrecognized observable " + obs)

                getattr(self, 'h_JetPt_%s_ch_MPIon_%sScaled' % (obs, label)).Fill(jet.pt(), obs_val)


    #---------------------------------------------------------------
    # Fill jet histograms
    #---------------------------------------------------------------
    def fill_jet_histograms(self, jetR, jp, jh, jch):

        R_label = str(jetR) + 'Scaled'

        # Fill jet histograms which are not dependant on angualrity
        if self.level in [None, 'ch']:
            getattr(self, 'hJetPt_ch_R%s' % R_label).Fill(jch.pt())
            getattr(self, 'hNconstit_Pt_ch_R%s' % R_label).Fill(jch.pt(), len(jch.constituents()))
        if self.level in [None, 'h']:
            getattr(self, 'hJetPt_h_R%s' % R_label).Fill(jh.pt())
            getattr(self, 'hNconstit_Pt_h_R%s' % R_label).Fill(jh.pt(), len(jh.constituents()))
        if self.level in [None, 'p']:
            getattr(self, 'hJetPt_p_R%s' % R_label).Fill(jp.pt())
            getattr(self, 'hNconstit_Pt_p_R%s' % R_label).Fill(jp.pt(), len(jp.constituents()))

        if self.level == None:
            if jp.pt():  # prevent divide by 0
                getattr(self, 'hJetPtRes_R%s' % R_label).Fill(jp.pt(), (jp.pt() - jch.pt()) / jp.pt())
            getattr(self, 'hResponse_JetPt_R%s' % R_label).Fill(jp.pt(), jch.pt())

        # Fill angularity histograms and response matrices
        for i, gs in enumerate(self.grooming_settings):
            gl = self.grooming_labels[i]
            self.fill_RMs(jetR, gs, gl, jp, jh, jch)


    #---------------------------------------------------------------
    # Fill jet histograms
    #---------------------------------------------------------------
    def fill_RMs(self, jetR, gs, gl, jp, jh, jch):

        # Groomed jets
        gshop_chh = fjcontrib.GroomerShop(jch, jetR, self.reclustering_algorithm)
        jet_ch_groomed_lund = self.utils.groom(gshop_chh, gs, jetR)
        if not jet_ch_groomed_lund:
            return

        gshop_h = fjcontrib.GroomerShop(jh, jetR, self.reclustering_algorithm)
        jet_h_groomed_lund = self.utils.groom(gshop_h, gs, jetR)
        if not jet_h_groomed_lund:
            return

        gshop_p = fjcontrib.GroomerShop(jp, jetR, self.reclustering_algorithm)
        jet_p_groomed_lund = self.utils.groom(gshop_p, gs, jetR)
        if not jet_p_groomed_lund:
            return

        label = "R%s_%s" % (jetR, gl)
        obs_dict = None
        for obs in self.observables:
            if obs == "theta_g":
                obs_dict = {
                    "p" : jet_p_groomed_lund.Delta() / jetR,
                    "h" : jet_h_groomed_lund.Delta() / jetR,
                    "ch": jet_ch_groomed_lund.Delta() / jetR }
            elif obs == "zg":
                if self.skip_zg:
                    continue
                obs_dict = {
                    "p" : jet_p_groomed_lund.z(),
                    "h" : jet_h_groomed_lund.z(),
                    "ch": jet_ch_groomed_lund.z() }
            else:
                raise ValueError("Unrecognized observable " + obs)

            if self.level in [None, 'ch']:
                getattr(self, 'h_JetPt_%s_ch_MPIoff_%sScaled' % (obs, label)).Fill(jch.pt(), obs_dict['ch'])

            if self.level in [None, 'h']:
                getattr(self, 'h_JetPt_%s_h_MPIoff_%sScaled' % (obs, label)).Fill(jh.pt(), obs_dict['h'])

            if self.level in [None, 'p']:
                getattr(self, 'h_JetPt_%s_p_MPIoff_%sScaled' % (obs, label)).Fill(jp.pt(), obs_dict['p'])

            if self.level == None:
                getattr(self, 'hResponse_%s_p_ch_MPIoff_%sScaled' % (obs, label)).Fill(obs_dict['p'], obs_dict['ch'])

                # Residual plots (with and without divisor in y-axis)
                getattr(self, "hDiff_JetPt_%s_MPIoff_%sScaled" % (obs, label)).Fill(
                    jch.pt(), obs_dict['p'] - obs_dict['ch'])
                if obs_dict['p']:  # prevent divide by 0
                    getattr(self, "hResidual_JetPt_%s_MPIoff_%sScaled" % (obs, label)).Fill(
                        jp.pt(), (obs_dict['p'] - obs_dict['ch']) / obs_dict['p'])

                # 4D response matrices for "forward folding" to ch level
                x = ([jch.pt(), jp.pt(), obs_dict['ch'], obs_dict['p']])
                x_array = array('d', x)
                getattr(self, 'hResponse_JetPt_%s_p_ch_MPIoff_%sScaled' % (obs, label)).Fill(x_array)

                x = ([jh.pt(), jp.pt(), obs_dict['h'], obs_dict['p']])
                x_array = array('d', x)
                getattr(self, 'hResponse_JetPt_%s_p_h_MPIoff_%sScaled' % (obs, label)).Fill(x_array)


    #---------------------------------------------------------------
    # Fill jet histograms for MPI (which are just the H-->CH RMs)
    #---------------------------------------------------------------
    def fill_jet_histograms_MPI(self, jetR, jp, jh, jch):

        for i, gs in enumerate(self.grooming_settings):
            gl = self.grooming_labels[i]

            # Groomed jets
            gshop_chh = fjcontrib.GroomerShop(jch, jetR, self.reclustering_algorithm)
            jet_ch_groomed_lund = self.utils.groom(gshop_chh, gs, jetR)
            if not jet_ch_groomed_lund:
                continue

            gshop_h = fjcontrib.GroomerShop(jh, jetR, self.reclustering_algorithm)
            jet_h_groomed_lund = self.utils.groom(gshop_h, gs, jetR)
            if not jet_h_groomed_lund:
                continue

            gshop_p = fjcontrib.GroomerShop(jp, jetR, self.reclustering_algorithm)
            jet_p_groomed_lund = self.utils.groom(gshop_p, gs, jetR)
            if not jet_p_groomed_lund:
                continue

            label = "R%s_%s" % (jetR, gl)
            obs_dict = None
            for obs in self.observables:
                if obs == "theta_g":
                    obs_dict = {
                        "p" : jet_p_groomed_lund.Delta() / jetR,
                        "h" : jet_h_groomed_lund.Delta() / jetR,
                        "ch": jet_ch_groomed_lund.Delta() / jetR }
                elif obs == "zg":
                    if self.skip_zg:
                        continue
                    obs_dict = {
                        "p" : jet_p_groomed_lund.z(),
                        "h" : jet_h_groomed_lund.z(),
                        "ch": jet_ch_groomed_lund.z() }
                else:
                    raise ValueError("Unrecognized observable " + obs)

                # 4D response matrices for "forward folding" from h to ch level
                x = ([jch.pt(), jh.pt(), obs_dict['ch'], obs_dict['h']])
                x_array = array('d', x)
                getattr(self, 'hResponse_JetPt_%s_h_ch_MPIon_%sScaled' % (obs, label)).Fill(x_array)


    #---------------------------------------------------------------
    # Initiate scaling of all histograms and print final simulation info
    #---------------------------------------------------------------
    def scale_print_final_info(self):

        # Scale all jet histograms by the appropriate factor from generated cross section
        # and the number of accepted events
        if not self.no_scale:
            scale_f = self.xsec / self.hNevents.GetBinContent(1)
            print("Weight MPIoff histograms by (cross section)/(N events) =", scale_f)

            MPI_scale_f = None
            if self.herwig_file_MPI:
                MPI_scale_f = self.xsec_MPI / self.hNeventsMPI.GetBinContent(1)
                print("Weight MPIon histograms by (cross section)/(N events) =", MPI_scale_f)

            self.scale_jet_histograms(scale_f, MPI_scale_f)
        print()

        self.hNevents.SetBinError(1, 0)
        if self.herwig_file_MPI:
            self.hNeventsMPI.SetBinError(1, 0)

        for jetR in self.jetR_list:
            jetR_str = str(jetR).replace('.', '')
            count1 = getattr(self, "count1_R%s" % jetR_str)
            count2 = getattr(self, "count2_R%s" % jetR_str)
            print(("For R=%s:  %i jets cut at first match criteria; " + \
                  "%i jets cut at second match criteria.") %
                  (str(jetR), count1, count2))
        print()


    #---------------------------------------------------------------
    # Scale all jet histograms by sigma/N
    #---------------------------------------------------------------
    def scale_jet_histograms(self, scale_f, MPI_scale_f):

        for jetR in self.jetR_list:
            hist_list_name = "hist_list_R%s" % str(jetR).replace('.', '')
            for h in getattr(self, hist_list_name):
                h.Scale(scale_f)

            hist_list_MPIon_name = "hist_list_MPIon_R%s" % str(jetR).replace('.', '')
            for h in getattr(self, hist_list_MPIon_name):
                h.Scale(MPI_scale_f)


################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Herwig7 debug level-1 output parser',
                                     prog=os.path.basename(__file__))
    parser.add_argument('-i', '--input-file', action='store', type=str, default='LHC.log',
                        help='Input .log file from Herwig7 analysis')
    parser.add_argument('-m', '--input-file-mpi', action='store', type=str, default=None,
                        help='Input .log file with MPI on from Herwig7 analysis')
    parser.add_argument('-o', '--output-dir', action='store', type=str, default='./', 
                        help='Output directory for generated ROOT file(s)')
    parser.add_argument('--tree-output-fname', default="AnalysisResults.root", type=str,
                        help="Filename for the (unscaled) generated particle ROOT TTree")
    parser.add_argument('--no-tree', default=False, action='store_true', 
                        help="Do not save tree of particle information, only create histograms")
    parser.add_argument('--no-match-level', help="Save simulation for only one level with " + \
                        "no matching. Options: 'p', 'h', 'ch'", default=None, type=str)
    parser.add_argument('--no-scale', help="Turn off rescaling all histograms by cross section / N",
                        action='store_true', default=False)
    parser.add_argument('-c', '--config-file', action='store', type=str,
                        default='config/angularity.yaml',
                        help="Path of config file for observable configurations")
    args = parser.parse_args()

    if args.no_match_level not in [None, 'p', 'h', 'ch']:
        print("ERROR: Unrecognized type %s. Please use 'p', 'h', or 'ch'" % args.type_only)
        exit(1)

    # If invalid configFile is given, exit
    if not os.path.exists(args.config_file):
        print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)

    process = herwig_parton_hadron(
        config_file=args.config_file, output_dir=args.output_dir, args=args)
    process.herwig_parton_hadron(args)
