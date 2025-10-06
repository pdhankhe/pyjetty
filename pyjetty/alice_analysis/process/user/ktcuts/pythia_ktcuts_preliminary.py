#!/usr/bin/env python3
# Code taken from Mateusz's https://github.com/blianggilman/alian/blob/main/alian/sandbox/jse/pythia_jse.py
# This is pythia on-the-fly to look at EECs if kt cut is applied
# can use --nev 100 or whatever to change the number of events
# method1 = take all EEC pairs, and calculate kT of each.
# python pythia_ktcuts_preliminary.py --nev 10000 --jet-pt-min 20 -o pythia_jse_output_ptmin20.root

from __future__ import print_function
import tqdm
import yaml
import argparse
import os
import numpy as np
import sys
import ROOT
import math
from array import array
# import yasp

import fastjet as fj
import fjcontrib
import pythiafjext
# import heppyy.util.fastjet_cppyy
# import heppyy.util.pythia8_cppyy
# import heppyy.util.heppyy_cppyy

# from cppyy.gbl import fastjet as fj
# from cppyy.gbl import Pythia8
# from cppyy.gbl.std import vector

# Analysis utilities
# from heppyy.util.mputils import logbins
# from heppyy.pythia_util import configuration as pyconf
from heppy.pythiautils import configuration as pyconf
from pyjetty.alice_analysis.process.base import jet_info

import ecorrel
import othercorrel

ROOT.TH1.SetDefaultSumw2(True)
ROOT.TH2.SetDefaultSumw2(True)
ROOT.TH3.SetDefaultSumw2(True)
# ROOT.THnSparseD.SetDefaultSumw2(True) # Apparently this doesn't exist

# from yasp import GenericObject
# from alian.sandbox.root_output import SingleRootFile

# we measured for k=1 and a>0
def angularity(jet, a, k, jetR):
    ang = 0.0
    for p in jet.constituents():
        dr = (np.sqrt((jet.delta_phi_to(p))**2 + (jet.eta() - p.eta())**2)) / jetR #jet.delta_R(p) / jetR
        pt = p.perp() / jet.perp()
        ang += ((dr)**a) * ((pt)**k)
    return ang

def mass(jet):
    m2 = jet.e()**2 - jet.px()**2 - jet.py()**2 - jet.pz()**2
    if m2 > 0:
        return math.sqrt(m2)
    return 0.0


def create_thn(name, title, dim, binnings=[]): # note: binnings could be given as (ptbins, rapi bins, obs bins)

    nbins = [len(x)-1 for x in binnings]
    print("NBINS", nbins)

    xmin = [x[0] for x in binnings]
    xmax = [x[-1] for x in binnings] 

    nbins_array = array('i', nbins) #('i', nbins_arr)
    xmin_array = array('d', xmin) #('d', xmin_arr)
    xmax_array = array('d', xmax) #('d', xmax_arr)
    h = ROOT.THnSparseD(name, name, dim, nbins_array, xmin_array, xmax_array)
    h.Sumw2()
    for i in range(0, dim):
      h.GetAxis(i).SetTitle(title[i])
      h.SetBinEdges(i, binnings[i])
    
    return h

def main():
    parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
    pyconf.add_standard_pythia_args(parser)
    parser.add_argument('-v', '--verbose', help="be verbose", default=False, action='store_true')
    parser.add_argument('--ncorrel', help='max n correlator', type=int, default=2)
    parser.add_argument('-o','--output', help='root output filename', default='pythia_jse_output.root', type=str)
    parser.add_argument('--jet-pt-min', help='jet pt min', default=20.0, type=float)
    parser.add_argument('--jet-pt-max', help='jet pt max', default=-1, type=float)
    parser.add_argument('--etadet', help='detector eta', default=2.5, type=float)
    parser.add_argument('--shape', help='fill the jet shape histograms', action='store_true', default=False)
    args = parser.parse_args()

    # pythia = Pythia8.Pythia()
    mycfg = [] #['Random:setSeed=on', 'Random:seed={}'.format(self.user_seed)]
    # mycfg.append('HadronLevel:all=off')
    pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
    if not pythia:
        print("[e] pythia initialization failed.")
        return
    # if args.nev < 10:
    #     args.nev = 10

    # jet finder
    # print the banner first
    fj.ClusterSequence.print_banner()
    print()

    Rs = [0.4] #[0.2, 0.4, 0.6]
    pt_bin_final = [20,40,60,80,100,120,150,200]
    pt_bin_collection_sum = np.zeros(len(pt_bin_final)-1)
    pt_bin_collection_counts = np.zeros(len(pt_bin_final)-1)
    pt_bin_avgs = []
    jet_defs = {}
    jet_selectors = {}
    for R in Rs:
        jet_defs[R] = fj.JetDefinition(fj.antikt_algorithm, R)
        jet_selectors[R] = fj.SelectorPtMin(args.jet_pt_min) * fj.SelectorAbsEtaMax(args.etadet - R * 1.05) #TODO: is this right?
        if args.jet_pt_max > 0:
            jet_selectors[R] *= fj.SelectorPtMax(args.jet_pt_max) * fj.SelectorPtMin(args.jet_pt_min) * fj.SelectorAbsEtaMax(args.etadet - R * 1.05)
        

    dphi_cut = -9999
    deta_cut = -9999

    # Make output file
    # fout = SingleRootFile(args.output)
    output_filepath = "/software/users/blianggi/mypyjetty/storage/ktcuts/rootfiles/" + args.output
    fout = ROOT.TFile(output_filepath, 'recreate')
    fout.cd()

    # bins
    pt_bins = np.linspace(0, 200, 201)
    RL_bins = np.logspace(np.log10(5E-3),np.log10(1),51) #(np.log10(1E-4),np.log10(1),51)
    ptRL_bins = np.logspace(np.log10(2E-1),np.log10(200),51)
    kt_bins = np.linspace(0, 7.5, 101) # (min, max, len of arr)
    kappa_bins = np.linspace(0, 0.25, 101)
    inverse_kt_bins = np.linspace(-10, 10, 201)
    inverse_Delta_bins = np.linspace(0, 10, 101)

    # Make histograms
    h_jet_pt = ROOT.TH1D("h_jet_pt", "h_jet_pt", 200, 0, 200)
    h_num_jets_per_event = ROOT.TH1I("h_num_jets_per_event", "h_num_jets_per_event", 5, 0, 5)
    h_lund_plane = ROOT.TH2D("h_lund_plane", "Lund Plane", 100, inverse_Delta_bins, 200, inverse_kt_bins)

    # h3D_pt_vs_RL_vs_kT = ROOT.TH3D("h3D_pt_vs_RL_vs_kT", "h3D_pt_vs_RL_vs_kT", 200, pt_bins, 50, RL_bins, 100, kt_bins)
    # h3D_pt_vs_ptRL_vs_kT = ROOT.TH3D("h3D_pt_vs_ptRL_vs_kT", "h3D_pt_vs_ptRL_vs_kT", 200, pt_bins, 50, ptRL_bins, 100, kt_bins)
    # h3D_pt_vs_ptRL_vs_kappa = ROOT.TH3D("h3D_pt_vs_ptRL_vs_kappa", "h3D_pt_vs_ptRL_vs_kappa", 200, pt_bins, 50, ptRL_bins, 100, kappa_bins)
    
    binnings = (pt_bins, RL_bins, ptRL_bins, kt_bins, kappa_bins)
    # print("BINNINGS", binnings)
    hND_all_pair_info = create_thn("hND_all_pair_info", "hND_all_pair_info", 5, binnings)
    # hND_all_pair_info = ROOT.THnSparseD("hND_all_pair_info", "hND_all_pair_info", 200, pt_bins, 50, ptRL_bins, 100, kappa_bins,)
    # h = ROOT.THnSparseD(name, name, dim, nbins_array, xmin_array, xmax_array)

    h_EEC_cutkTbelow05 = ROOT.TH1D("h_EEC_cutkTbelow05", "h_EEC: k_{T} #LT 0.5", 50, RL_bins)
    h_EEC_cutkTbelow1 = ROOT.TH1D("h_EEC_cutkTbelow1", "h_EEC: k_{T} #LT 1.0", 50, RL_bins)
    h_EEC_cutkTabove1 = ROOT.TH1D("h_EEC_cutkTabove1", "h_EEC: k_{T} #GEQ 1.0", 50, RL_bins)
    
    kt_of_all_pairs = ROOT.TH1D("kt_all_pairs", "kt_all_pairs", 100, kt_bins)
    kt_of_all_pairs_pt2040 = ROOT.TH1D("kt_all_pairs_pt2040", "kt_all_pairs_pt2040", 100, kt_bins)
    kt_of_all_pairs_pt4060 = ROOT.TH1D("kt_all_pairs_pt4060", "kt_all_pairs_pt4060", 100, kt_bins)
    kt_of_all_pairs_pt6080 = ROOT.TH1D("kt_all_pairs_pt6080", "kt_all_pairs_pt6080", 100, kt_bins)



    # loop through the events
    print("NUMBER OF EVENTS:", args.nev)
    if args.nev < 10:
        args.nev = 10
    
    count_jets = 0
    pbar = tqdm.tqdm(total=args.nev) #progress bar
    while pbar.n < args.nev: # looping through events here!!
        if not pythia.next():
            continue
        # print("Event #", pbar.n)

        parts = fj.vectorPJ([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal() and p.isCharged()])
        pythia_particles_chfinal = [p for p in pythia.event if p.isFinal() and p.isCharged()]
        # now assign pid to each particle PJ
        for ipj, (p,pj) in enumerate(zip(pythia_particles_chfinal,parts)):
            jetinfo = jet_info.JetInfo() # even though this is particle level info lol
            jetinfo.particle_pid = p.id()
            parts[ipj].set_python_info(jetinfo)
        # [print(p.id()) for p in pythia.event if p.isFinal() and p.isCharged()]
        # [print(pj.python_info().particle_pid) for pj in parts]
        
        # parts = vector[fj.PseudoJet]([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal() and p.isCharged()])
        # parts = pythiafjext.vectorize(pythia, True, -1, 1, False)

        for R in Rs:
            # print("analyzing jet R =", R)
            #get the pairs
            jets_ch = fj.sorted_by_pt(jet_selectors[R](jet_defs[R](parts)))
            if len(jets_ch) > 0:
                pbar.update(1)
                # print(len(jets_ch), "jets found")
                h_num_jets_per_event.Fill(len(jets_ch))
            
            '''
            # get the parents
            fs_parton_5 = fj.PseudoJet(pythia.event[5].px(), pythia.event[5].py(), pythia.event[5].pz(), pythia.event[5].e())
            fs_parton_6 = fj.PseudoJet(pythia.event[6].px(), pythia.event[6].py(), pythia.event[6].pz(), pythia.event[6].e())
            self.parents = [fs_parton_5, fs_parton_6] # parent partons in dijet
            self.parent_ids = [pythia.event[5].id(), pythia.event[6].id()]
            '''

            for i_jch, jch in enumerate(jets_ch):
 
                for i in range(len(pt_bin_final)-1):
                    if pt_bin_final[i] <= jch.perp() < pt_bin_final[i+1]:
                        pt_bin_collection_sum[i] += jch.perp()
                        pt_bin_collection_counts[i] += 1
                        break
                '''        
                for i_parent, parent in enumerate(self.parents):
                    parentmatch_name = "parent%imatch" % i_parent
                    if jch.delta_R(parent) < jet_matching_distance * jetR:
                        print("match found!!!!")
                        match = getattr(self, parentmatch_name)
                        if not match:
                            setattr(self, parentmatch_name, jch)
                        else:  # Already found a match
                            # Set flag value so that we know to ignore this one
                            setattr(self, parentmatch_name, 0)
                '''
            pt_bin_avgs = np.divide(pt_bin_collection_sum, pt_bin_collection_counts, out=np.zeros_like(pt_bin_collection_sum), where=pt_bin_collection_counts>0)

            '''
            # If we have matches, fill histograms
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
                    elif (parton_id in self.beauty_pdg_ids and (self.initscat == 4 or self.initscat == 5)):
                        parton_types += ["beauty"]
                elif parton_id in self.gluon_pdg_ids:
                    parton_types += ["gluon"]
                # if not self.replaceKPpairs:
                parton_types += ["inclusive"]
            '''

            for ijet, jet in enumerate(jets_ch):
                # print("analyzing jet #", ijet)

                jet_pt = jet.perp()
                jet_pt_avg = -1
                if jet_pt >= 200:
                    jet_pt_avg = 999
                for i in range(len(pt_bin_final)-1):
                    if jet_pt >= pt_bin_final[i] and jet_pt < pt_bin_final[i+1]:
                        jet_pt_avg = pt_bin_avgs[i]
                        break  

                h_jet_pt.Fill(jet_pt)
                constituents = fj.sorted_by_pt(jet.constituents())
                c_select = fj.vectorPJ()
                trk_thrd = 1 
                for c in constituents:
                    if c.pt() < trk_thrd:
                        break
                    c_select.append(c)
                # print("Event #", pbar.n, ": getting eec pairs for jet #", ijet)
                new_corr = ecorrel.CorrelatorBuilder(c_select, jet.perp(), 2, 1, dphi_cut, deta_cut)
        
                # method 1!
                EEC_indicies1 = new_corr.correlator(2).indices1() # a list of indices referring to the particles in pythia event
                EEC_indicies2 = new_corr.correlator(2).indices2()
                for index in range(new_corr.correlator(2).rs().size()):
                    part1_index = int(EEC_indicies1[index]) # getting the specific index of the particle in pythia event
                    part1 = c_select[part1_index] # getting the particle
                    part2_index = int(EEC_indicies2[index])
                    part2 = c_select[part2_index]

                    softer_pt = min(part1.pt(), part2.pt())
                    rl = new_corr.correlator(2).rs()[index]
                    kt = softer_pt * rl
                    weight = part1.pt() * part2.pt() / (jet_pt*jet_pt)
                    tuple_to_fill = array('d', [jet_pt, rl, jet_pt_avg*rl, kt, kt/jet_pt])
                    hND_all_pair_info.Fill(tuple_to_fill, weight) # all jet pt that are not between 20-200 will not be accurate
                    # h3D_pt_vs_RL_vs_kT.Fill(jet_pt, rl, kt, weight)
                    # if jet_pt_avg != -1:
                    #     h3D_pt_vs_ptRL_vs_kT.Fill(jet_pt, jet_pt_avg*rl, kt, weight)
                    #     h3D_pt_vs_ptRL_vs_kappa.Fill(jet_pt, jet_pt_avg*rl, kt/jet_pt, weight)
                
                    # make cuts on kT
                    if kt < 0.5:
                        h_EEC_cutkTbelow05.Fill(rl)
                    if kt < 1:
                        h_EEC_cutkTbelow1.Fill(rl)
                    else:
                        h_EEC_cutkTabove1.Fill(rl)
                # print("Number of final state particles:", len(c_select))
                # [print(c) for c in c_select]
                    

                # new:
                #first i need to get all the pairs
                # then I need to get the lund plane
                # then I need to go through every pair of eec, and find the most common last splitting in the lund plane
                # then I need to put that splitting into a histogram of kt_splittings

                # in lund plane: l.pair() returns the pair aka sum of two subjets
                #                l.harder() returs subjet w/ larger pT
                #                l.softer() returs subjet w/ smaller pT
                # https://github.com/fdreyer/LundPlane/blob/03414c2b126e2440155612eb4bf792bd229133d7/LundGenerator.hh#L52
                # print("Event #", pbar.n, ": getting lund plane for jet #", ijet)
                jet_def_ca = fj.JetDefinition(fj.cambridge_algorithm, R)
                lund_gen = fjcontrib.LundGenerator(jet_def_ca)
                lunds = lund_gen.result(jet)
                
                # print("lunds:")
                for i,l in enumerate(lunds):
                    # print(i, l.pair(), "--<", l.harder(), "&", l.softer())
                    h_lund_plane.Fill(np.log(1/l.Delta()), np.log(l.kt()))
                    kt_of_all_pairs.Fill(l.kt())
                    if jet_pt >= 20 and jet_pt < 40:
                        kt_of_all_pairs_pt2040.Fill(l.kt())
                    elif jet_pt >= 40 and jet_pt < 60:
                        kt_of_all_pairs_pt4060.Fill(l.kt())
                    elif jet_pt >= 60 and jet_pt < 80:
                        kt_of_all_pairs_pt6080.Fill(l.kt())

                

                # for i in range(len(EEC_rs)):
                '''
                EEC_new_corr = new_corr.correlator(2)

                for index in range(EEC_new_corr.rs().size()):
                    EEC_indicies1 = EEC_new_corr.indices1() # a list of indices referring to the particles in pythia event
                    EEC_indicies2 = EEC_new_corr.indices2()
                    # event_index1 = c_select[EEC_indicies1[index]].user_index()
                    # event_index2 = c_select[EEC_indicies2[index]].user_index()

                    # print("testing", EEC_indicies1, "//", EEC_indicies2)
                    # [print("ONE", e) for e in EEC_indicies1]
                    # [print("TWO", e) for e in EEC_indicies2]
                    part1_index = int(EEC_indicies1[index]) # getting the specific index of the particle in pythia event
                    part1 = c_select[part1_index] # getting the particle
                    part2_index = int(EEC_indicies2[index])
                    part2 = c_select[part2_index]
                    # print("testing", part1_index, part1, part2_index, part2)
                print("Number of final state particles:", len(c_select))
                [print(c) for c in c_select]

                    
                    # check if two particles are from the same latest splitting

                    
                    #get the two partivles
                    # how to find their common denominator?
                    # can loop through the lund plane by doing:
                    # for splitting in lunds:
                '''

                    



    '''
    jet_def_wta = fj.JetDefinition(fj.cambridge_algorithm, 1.0)
    jet_def_wta.set_recombination_scheme(fj.WTA_pt_scheme)
    reclusterer_wta =  fj.contrib.Recluster(jet_def_wta)

    sd01 = fj.contrib.SoftDrop(0, 0.1, 1.0)
    sd02 = fj.contrib.SoftDrop(0, 0.2, 1.0)
    '''

    h_lund_plane.GetXaxis().SetTitle("ln(1/#DeltaR)")
    h_lund_plane.GetYaxis().SetTitle("ln(k_{T})")
    hND_all_pair_info.GetAxis(0).SetTitle("p_{T,jet}")
    hND_all_pair_info.GetAxis(1).SetTitle("R_{L}")
    hND_all_pair_info.GetAxis(2).SetTitle("#LTp_{T,jet}#GTR_{L}")
    hND_all_pair_info.GetAxis(3).SetTitle("k_{T}")
    hND_all_pair_info.GetAxis(4).SetTitle("#kappa")

    h_EEC_cutkTbelow05.GetXaxis().SetTitle("R_{L}")
    h_EEC_cutkTbelow1.GetXaxis().SetTitle("R_{L}")
    h_EEC_cutkTabove1.GetXaxis().SetTitle("R_{L}")


    # Save histograms
    h_jet_pt.Write()
    h_num_jets_per_event.Write()
    h_lund_plane.Write()

    # h3D_pt_vs_RL_vs_kT.Write()
    # h3D_pt_vs_ptRL_vs_kT.Write()
    # h3D_pt_vs_ptRL_vs_kappa.Write()
    hND_all_pair_info.Write()
    h_EEC_cutkTbelow05.Write()
    h_EEC_cutkTbelow1.Write()
    h_EEC_cutkTabove1.Write()

    kt_of_all_pairs.Write()
    kt_of_all_pairs_pt2040.Write()
    kt_of_all_pairs_pt4060.Write()
    kt_of_all_pairs_pt6080.Write()
    fout.Close()
    

if __name__ == '__main__':
    main()