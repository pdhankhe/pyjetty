#!/usr/bin/env python3

import numpy as np
import argparse
import os
from array import array
import pandas as pd

from pyjetty.alice_analysis.process.base import process_io, process_utils, jet_info, process_base
from pyjetty.alice_analysis.process.user.substructure import process_data_base
import pyjetty.alihfjets.dev.hfjet.process.base.process_io_data_hf as hfdio
from pyjetty.mputils.mputils import perror, pinfo, pwarning
from pyjetty.mputils import treewriter, jet_analysis

import fastjet as fj
import fjext
import fjtools
import tqdm
import fjcontrib
import math as ma
import yaml
import ROOT
import time
import ecorrel
ROOT.gROOT.SetBatch(True)


import ROOT
ROOT.gROOT.SetBatch(True)

def linbins(xmin, xmax, nbins):
	lspace = np.linspace(xmin, xmax, nbins+1)
	arr = array('f', lspace)
	return arr

def logbins(xmin, xmax, nbins):
	lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
	arr = array('f', lspace)
	return arr


class HFAnalysisInvMass(hfdio.HFAnalysis, process_base.ProcessBase):
    def __init__(self, **kwargs):
        super(HFAnalysisInvMass, self).__init__(**kwargs)
        self.initialize_config()
        print(self.config_file)
    
    def HFAnalysis(self):
        print("dataframe will be created jet analysis")
        start_time = time.time()
        print('--- {} seconds ---'.format(time.time() - start_time))
        
        # Initialize configuration and histograms
        self.initialize_config()
        self.initialize_user_output_objects()
        
        # Find jets and fill histograms
        print('Find jets...')
        self.analyzeEvents()
        
        
        # Plot histograms
        print('Save histograms...')
        print('--- {} seconds ---'.format(time.time() - start_time))
        process_base.ProcessBase.save_output_objects(self)
            
    def initialize_config(self):
        print("initializing the configuration")
        with open(self.config_file, 'r') as stream:
            self.config = yaml.load(stream,Loader=yaml.FullLoader)

        # Set configuration for analysis
        self.jetR_list = self.config["jetR"]
        self.JetR = self.jetR_list[0]
        # Mass assumption for jet reconstruction
        self.track_mass = self.config["track_mass"]
        self.track_random_mass = self.config["track_random_mass"]   # randomly assign K or p mass to tracks
        
        self.sd_zcut = self.config["sd_zcut"]
        self.sd_beta = self.config["sd_beta"]
        self.observable_list = self.config["observable_list"]
        # Detector kinematics
        self.eta_max = 0.9


    #---------------------------------------------------------------
    # Initialize histograms
    #---------------------------------------------------------------
    def initialize_user_output_objects(self):
        
        name = "hNevents"
        h = ROOT.TH1F("hNevents", "hNevents", 2, array('f', [-0.5, 0.5, 1.5]) )
        h.GetYaxis().SetTitle('counts')
        setattr(self, name, h)
        
        
        title = [ '#it{p}_{T}^{ch jet}', '#it{p}_{T}^{D^{0}}', 'Invmass', '#it{R}_{L}']
        
        RL_bins  = np.logspace(np.log10(1E-4), np.log10(1), 51)
        pt_bins_jet = np.linspace(0, 60, 61)
        pt_bins_dmeson = np.linspace(0, 60, 61)
        dmeson_mass_bin = np.linspace(1.7, 2.07, 371)

        dim = 4
        nbins  = [len(pt_bins_jet)-1, len(pt_bins_dmeson)-1, len(dmeson_mass_bin)-1, len(RL_bins)-1]
        min_li = [pt_bins_jet[0],     pt_bins_dmeson[0],     dmeson_mass_bin[0],     RL_bins[0]    ]
        max_li = [pt_bins_jet[-1],    pt_bins_dmeson[-1],    dmeson_mass_bin[-1],    RL_bins[-1]   ]

        nbins = (nbins)
        xmin = (min_li)
        xmax = (max_li)

        nbins_array = array('i', nbins)
        xmin_array = array('d', xmin)
        xmax_array = array('d', xmax)
        
        for observable in self.observable_list:
            name = 'hsparse_R{}_{}'.format(self.JetR, observable)
            h = ROOT.THnSparseD( name, name, dim,  nbins_array, xmin_array, xmax_array)
            h.Sumw2()
            
            for axis in range(0,dim):
                h.GetAxis(axis).SetTitle(title[axis])
                if axis == 0 or axis == 1:
                    h.SetBinEdges(axis, pt_bins_jet)
                if axis == 2:
                    h.SetBinEdges(axis, dmeson_mass_bin)
                if axis == 3:
                    h.SetBinEdges(axis, RL_bins)
                    
            setattr(self, name, h)
        
    def analyzeEvents(self):
        self.hfio.execute_analyses_on_inputfile(self.input_file)
        
    def event_by_event_D_jet_reconstruction(self, df):
        #print(df)
        m_array = np.full((self.df_tracks['ParticlePt'].values.size), 0.1396)
        djmm = fjtools.DJetMatchMaker()
        djmm.set_ch_pt_eta_phi_m(self.df_tracks['ParticlePt'].values, self.df_tracks['ParticleEta'].values, self.df_tracks['ParticlePhi'].values, m_array)
        djmm.set_Ds_pt_eta_phi_m(df['pt_cand'].values, df['eta_cand'].values, df['phi_cand'].values, df['inv_mass'].values)
        djmm.set_daughters0_pt_eta_phi(df['pt_prong0'].values, df['eta_prong0'].values, df['phi_prong0'].values)
        djmm.set_daughters1_pt_eta_phi(df['pt_prong1'].values, df['eta_prong1'].values, df['phi_prong1'].values)


        #run for each D candidate to build jet
        for id0, d0 in enumerate(djmm.Ds):
            #daughter tracks matching
            _parts_and_ds = djmm.match(0.005, id0)
            #replacing daughter tracks with matched D0 candidate
            #including D0
            _parts_and_ds.push_back(d0)
    
            #jet reconstruction with D0 and charged particle
            ja = jet_analysis.JetAnalysis(jet_R = self.JetR, jet_RecombinationScheme=fj.E_scheme, particle_eta_max=0.9, jet_pt_min=5.0)
            ja.analyze_event(_parts_and_ds)
            
            if len(ja.jets) < 1:
                continue

            jets = ja.jets_as_psj_vector()

            #filtering D0 jets
            djets = djmm.filter_D0_jets(jets)
            
            if len(djets) > 0:
                j = djets[0]
                dcand = djmm.get_Dcand_in_jet(j)
                # Extract information for EEC
                constituents = fj.sorted_by_pt(j.constituents())
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
                
                x_array = array( 'd', ( j.perp(), dcand[0].perp(), dcand[0].m(), j.delta_R(dcand[0])))
                getattr(self, 'hsparse_R{}_{}'.format(self.JetR, self.observable_list[0])).Fill(x_array)
                
                new_corr = ecorrel.CorrelatorBuilder(c_select, dcand, j.perp(), 2, 1, dphi_cut, deta_cut)
                
                for index in range(new_corr.correlator(2).rs().size()):
                    x_EEC_array = array( 'd', ( j.perp(), dcand[0].perp(), dcand[0].m(), new_corr.correlator(2).rs()[index]))
                    getattr(self, 'hsparse_R{}_{}'.format(self.JetR, self.observable_list[1])).Fill(x_array, new_corr.correlator(2).weights()[index])
                    
            if len(djets) > 1:
                perror("more than one jet per D candidate?")
                
        return True
	

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='D0 analysis on alice data', prog=os.path.basename(__file__))
    parser.add_argument('-c', '--configFile', help='Path of config file for analysis', type=str,metavar='configFile',default='config/configcuts.yaml', required=True)
    parser.add_argument('-f', '--inputFile', action='store',
                      type=str, metavar='inputFile',
                      default='AnalysisResults.root',
                      help='Path of ROOT file containing TTrees')
    parser.add_argument('-o', '--outputDir', action='store',
                      type=str, metavar='outputDir',
                      default='./TestOutput',
                      help='Output directory for QA plots to be written to')
    args = parser.parse_args()
    
    print('Configuring...')
    print('inputFile: \'{0}\''.format(args.inputFile))
    print('configFile: \'{0}\''.format(args.configFile))
    print('ouputDir: \'{0}\"'.format(args.outputDir))
    print('----------------------------------------------------------------')
  
    if not os.path.exists(args.configFile):
        print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)
        
    hfaio = hfdio.HFAnalysisIO()
    
    hfa = HFAnalysisInvMass(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir, hfio=hfaio)
    print("Set event selection cuts here")
    hfa.event_selection.add_selection_range_abs('z_vtx_reco', 10)
    hfa.event_selection.add_selection_equal('is_ev_rej', 0)

    print("Set D0 pt independent selection cuts here")
    #topomatic cut suggested by D2H.
    hfa.d0_selection.add_selection_range_abs('max_norm_d0d0exp',2)
    hfa.d0_selection.add_selection_range('pt_cand', 2, 1e3)
    hfa.d0_selection.add_selection_range_abs('eta_cand', 0.8)
    hfa.d0_selection.add_selection_nsig('nsigTPC_Pi_0', 'nsigTOF_Pi_0', 'nsigTPC_K_1', 'nsigTOF_K_1', 'nsigTPC_Pi_1', 'nsigTOF_Pi_1', 'nsigTPC_K_0', 'nsigTOF_K_0', 3, -900)

    hfa.d0_gen_selection.add_selection_range('dau_in_acc',0, 1.1)
    hfa.d0_gen_selection.add_selection_range_abs('eta_cand',0.8)
    hfa.d0_gen_selection.add_selection_range('pt_cand',2,1e3)


    hfaio.add_analysis(hfa)
    print("Start reading input output and perform analysis")
    hfa.HFAnalysis()
