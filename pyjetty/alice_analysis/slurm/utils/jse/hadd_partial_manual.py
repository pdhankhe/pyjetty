import ROOT
import argparse
import os


# python3 -u hadd_manual.py --dir /global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/57259276

is_final_merge = False
rap = False

suffix = "" #"LHC18f3/803"
subdirs = ['0', '1', '2'] # ppb
runlist = ['559348', '559361', '559362', '559385', '559387', '559408', '559409', '559410', '559437', '559443', '559444', '559456']



def finishing_merge(dir, i):

    print("LOADING ROOUNFOLD")
    roo = False
    while roo==False:
        roo = ROOT.gSystem.Load("libRooUnfold")
    print(roo)
    roo2 = False
    while roo2==False:
        roo2 = ROOT.gSystem.Load("libRooUnfold.so")
    print(roo2)

    # fin1 = ROOT.TFile.Open("{}/merged_{}_{}.root".format(dir, subdirs[0], runlist[0]))
    fin1 = ROOT.TFile.Open("{}/{}/{}/0/response.root".format(dir, subdirs[0], runlist[0]))

    resp_jetpt_merged = fin1.resp_jetpt
    resp_groomed_jetpt_merged = fin1.resp_groomed_jetpt



    jet_match_gen_eff_num_ungroomed_merged = fin1.jet_match_gen_eff_num_ungroomed
    jet_all_gen_eff_den_ungroomed_merged = fin1.jet_all_gen_eff_den_ungroomed
    jet_match_rec_pur_num_ungroomed_merged = fin1.jet_match_rec_pur_num_ungroomed
    jet_all_rec_pur_den_ungroomed_merged = fin1.jet_all_rec_pur_den_ungroomed
    jet_efficiency_ungroomed_merged = fin1.jet_efficiency_ungroomed
    jet_purity_ungroomed_merged = fin1.jet_purity_ungroomed

    jet_match_gen_eff_num_groomed_merged = fin1.jet_match_gen_eff_num_groomed
    jet_all_gen_eff_den_groomed_merged = fin1.jet_all_gen_eff_den_groomed
    jet_match_rec_pur_num_groomed_merged = fin1.jet_match_rec_pur_num_groomed
    jet_all_rec_pur_den_groomed_merged = fin1.jet_all_rec_pur_den_groomed
    jet_efficiency_groomed_merged = fin1.jet_efficiency_groomed
    jet_purity_groomed_merged = fin1.jet_purity_groomed

    pair_match_gen_eff_num_AA_merged = fin1.pair_match_gen_eff_num_AA
    pair_all_gen_eff_den_AA_merged = fin1.pair_all_gen_eff_den_AA
    pair_match_rec_pur_num_AA_merged = fin1.pair_match_rec_pur_num_AA
    pair_all_rec_pur_den_AA_merged = fin1.pair_all_rec_pur_den_AA
    pair_efficiency_AA_merged = fin1.pair_efficiency_AA
    pair_purity_AA_merged = fin1.pair_purity_AA

    pair_match_gen_eff_num_AB_merged = fin1.pair_match_gen_eff_num_AB
    pair_all_gen_eff_den_AB_merged = fin1.pair_all_gen_eff_den_AB
    pair_match_rec_pur_num_AB_merged = fin1.pair_match_rec_pur_num_AB
    pair_all_rec_pur_den_AB_merged = fin1.pair_all_rec_pur_den_AB
    pair_efficiency_AB_merged = fin1.pair_efficiency_AB
    pair_purity_AB_merged = fin1.pair_purity_AB

    pair_match_gen_eff_num_BB_merged = fin1.pair_match_gen_eff_num_BB
    pair_all_gen_eff_den_BB_merged = fin1.pair_all_gen_eff_den_BB
    pair_match_rec_pur_num_BB_merged = fin1.pair_match_rec_pur_num_BB
    pair_all_rec_pur_den_BB_merged = fin1.pair_all_rec_pur_den_BB
    pair_efficiency_BB_merged = fin1.pair_efficiency_BB
    pair_purity_BB_merged = fin1.pair_purity_BB

    pair_match_gen_eff_num_rad_merged = fin1.pair_match_gen_eff_num_rad
    pair_all_gen_eff_den_rad_merged = fin1.pair_all_gen_eff_den_rad
    pair_match_rec_pur_num_rad_merged = fin1.pair_match_rec_pur_num_rad
    pair_all_rec_pur_den_rad_merged = fin1.pair_all_rec_pur_den_rad
    pair_efficiency_rad_merged = fin1.pair_efficiency_rad
    pair_purity_rad_merged = fin1.pair_purity_rad


    lund_matched_gen_merged = fin1.lund_matched_gen
    lund_matched_rec_merged = fin1.lund_matched_rec
    lund_all_gen_merged = fin1.lund_all_gen
    lund_all_rec_merged = fin1.lund_all_rec
    lund_split_efficiency_merged = fin1.lund_split_efficiency
    lund_split_purity_merged = fin1.lund_split_purity
    summary_efficiencies_merged = fin1.summary_efficiencies

    
    
    for subdir in subdirs:
        for run in runlist:

            for k in range(400): #10 # FIX THIS

                if subdir == subdirs[0] and run == runlist[0] and k == 0:
                    continue
        
                print("files from {}".format(run))

                try:
                    # fin = ROOT.TFile.Open("{}/merged_{}.root".format(dir, run))
                    # fin = ROOT.TFile.Open("{}/merged_{}_{}.root".format(dir, subdirs[0], run))
                    fin = ROOT.TFile.Open("{}/{}/{}/{}/response.root".format(dir, subdir, run, k))

                    resp_jetpt_merged.Add(fin.resp_jetpt)
                    resp_groomed_jetpt_merged.Add(fin.resp_groomed_jetpt)


                    jet_match_gen_eff_num_ungroomed_merged.Add(fin.jet_match_gen_eff_num_ungroomed)
                    jet_all_gen_eff_den_ungroomed_merged.Add(fin.jet_all_gen_eff_den_ungroomed)
                    jet_match_rec_pur_num_ungroomed_merged.Add(fin.jet_match_rec_pur_num_ungroomed)
                    jet_all_rec_pur_den_ungroomed_merged.Add(fin.jet_all_rec_pur_den_ungroomed)
                    jet_efficiency_ungroomed_merged.Add(fin.jet_efficiency_ungroomed)
                    jet_purity_ungroomed_merged.Add(fin.jet_purity_ungroomed)

                    jet_match_gen_eff_num_groomed_merged.Add(fin.jet_match_gen_eff_num_groomed)
                    jet_all_gen_eff_den_groomed_merged.Add(fin.jet_all_gen_eff_den_groomed)
                    jet_match_rec_pur_num_groomed_merged.Add(fin.jet_match_rec_pur_num_groomed)
                    jet_all_rec_pur_den_groomed_merged.Add(fin.jet_all_rec_pur_den_groomed)
                    jet_efficiency_groomed_merged.Add(fin.jet_efficiency_groomed)
                    jet_purity_groomed_merged.Add(fin.jet_purity_groomed)

                    pair_match_gen_eff_num_AA_merged.Add(fin.pair_match_gen_eff_num_AA)
                    pair_all_gen_eff_den_AA_merged.Add(fin.pair_all_gen_eff_den_AA)
                    pair_match_rec_pur_num_AA_merged.Add(fin.pair_match_rec_pur_num_AA)
                    pair_all_rec_pur_den_AA_merged.Add(fin.pair_all_rec_pur_den_AA)
                    pair_efficiency_AA_merged.Add(fin.pair_efficiency_AA)
                    pair_purity_AA_merged.Add(fin.pair_purity_AA)

                    pair_match_gen_eff_num_AB_merged.Add(fin.pair_match_gen_eff_num_AB)
                    pair_all_gen_eff_den_AB_merged.Add(fin.pair_all_gen_eff_den_AB)
                    pair_match_rec_pur_num_AB_merged.Add(fin.pair_match_rec_pur_num_AB)
                    pair_all_rec_pur_den_AB_merged.Add(fin.pair_all_rec_pur_den_AB)
                    pair_efficiency_AB_merged.Add(fin.pair_efficiency_AB)
                    pair_purity_AB_merged.Add(fin.pair_purity_AB)

                    pair_match_gen_eff_num_BB_merged.Add(fin.pair_match_gen_eff_num_BB)
                    pair_all_gen_eff_den_BB_merged.Add(fin.pair_all_gen_eff_den_BB)
                    pair_match_rec_pur_num_BB_merged.Add(fin.pair_match_rec_pur_num_BB)
                    pair_all_rec_pur_den_BB_merged.Add(fin.pair_all_rec_pur_den_BB)
                    pair_efficiency_BB_merged.Add(fin.pair_efficiency_BB)
                    pair_purity_BB_merged.Add(fin.pair_purity_BB)

                    pair_match_gen_eff_num_rad_merged.Add(fin.pair_match_gen_eff_num_rad)
                    pair_all_gen_eff_den_rad_merged.Add(fin.pair_all_gen_eff_den_rad)
                    pair_match_rec_pur_num_rad_merged.Add(fin.pair_match_rec_pur_num_rad)
                    pair_all_rec_pur_den_rad_merged.Add(fin.pair_all_rec_pur_den_rad)
                    pair_efficiency_rad_merged.Add(fin.pair_efficiency_rad)
                    pair_purity_rad_merged.Add(fin.pair_purity_rad)


                    lund_matched_gen_merged.Add(fin.lund_matched_gen)
                    lund_matched_rec_merged.Add(fin.lund_matched_rec)
                    lund_all_gen_merged.Add(fin.lund_all_gen)
                    lund_all_rec_merged.Add(fin.lund_all_rec)
                    lund_split_efficiency_merged.Add(fin.lund_split_efficiency)
                    lund_split_purity_merged.Add(fin.lund_split_purity)
                    summary_efficiencies_merged.Add(fin.summary_efficiencies)

                    

                except:
                    print("skipping files from {}, {}".format(subdirs[0], run))
                    # print("skipping files from {}, {}".format(i, run))
                    nonexistent.append((subdirs[0],run))
                    continue

    fout = ROOT.TFile("{}/response_merged_partial.root".format(dir), "RECREATE")

    # RECALCULATE EFF AND PURITY
    jet_eff_ungroomed_new = jet_match_gen_eff_num_ungroomed_merged.Clone("jet_efficiency_ungroomed_new")
    jet_eff_ungroomed_new.Divide(jet_all_gen_eff_den_ungroomed_merged)
    jet_pur_ungroomed_new = jet_match_rec_pur_num_ungroomed_merged.Clone("jet_purity_ungroomed_new")
    jet_pur_ungroomed_new.Divide(jet_all_rec_pur_den_ungroomed_merged)

    jet_eff_groomed_new = jet_match_gen_eff_num_groomed_merged.Clone("jet_efficiency_groomed_new")
    jet_eff_groomed_new.Divide(jet_all_gen_eff_den_groomed_merged)
    jet_pur_groomed_new = jet_match_rec_pur_num_groomed_merged.Clone("jet_purity_groomed_new")
    jet_pur_groomed_new.Divide(jet_all_rec_pur_den_groomed_merged)

    pair_eff_AA_new = pair_match_gen_eff_num_AA_merged.Clone("pair_efficiency_AA_new")
    pair_eff_AA_new.Divide(pair_all_gen_eff_den_AA_merged)
    pair_pur_AA_new = pair_match_rec_pur_num_AA_merged.Clone("pair_purity_AA_new")
    pair_pur_AA_new.Divide(pair_all_rec_pur_den_AA_merged)

    pair_eff_AB_new = pair_match_gen_eff_num_AB_merged.Clone("pair_efficiency_AB_new")
    pair_eff_AB_new.Divide(pair_all_gen_eff_den_AB_merged)
    pair_pur_AB_new = pair_match_rec_pur_num_AB_merged.Clone("pair_purity_AB_new")
    pair_pur_AB_new.Divide(pair_all_rec_pur_den_AB_merged)

    pair_eff_BB_new = pair_match_gen_eff_num_BB_merged.Clone("pair_efficiency_BB_new")
    pair_eff_BB_new.Divide(pair_all_gen_eff_den_BB_merged)
    pair_pur_BB_new = pair_match_rec_pur_num_BB_merged.Clone("pair_purity_BB_new")
    pair_pur_BB_new.Divide(pair_all_rec_pur_den_BB_merged)

    pair_eff_rad_new = pair_match_gen_eff_num_rad_merged.Clone("pair_efficiency_rad_new")
    pair_eff_rad_new.Divide(pair_all_gen_eff_den_rad_merged)
    pair_pur_rad_new = pair_match_rec_pur_num_rad_merged.Clone("pair_purity_rad_new")
    pair_pur_rad_new.Divide(pair_all_rec_pur_den_rad_merged)

    lund_split_efficiency_new = lund_matched_gen_merged.Clone("lund_split_efficiency_new")
    lund_split_efficiency_new.Divide(lund_all_gen_merged)
    lund_split_purity_new = lund_matched_rec_merged.Clone("lund_split_purity_new")
    lund_split_purity_new.Divide(lund_all_rec_merged)


    jet_eff_ungroomed_new.Write()
    jet_pur_ungroomed_new.Write()
    jet_eff_groomed_new.Write()
    jet_pur_groomed_new.Write()

    pair_eff_AA_new.Write()
    pair_pur_AA_new.Write()
    pair_eff_AB_new.Write()
    pair_pur_AB_new.Write()
    pair_eff_BB_new.Write()
    pair_pur_BB_new.Write()
    pair_eff_rad_new.Write()
    pair_pur_rad_new.Write()

    lund_split_efficiency_new.Write()
    lund_split_purity_new.Write()


    # WRITE FILES
    resp_jetpt_merged.Write()
    resp_groomed_jetpt_merged.Write()


    jet_match_gen_eff_num_ungroomed_merged.Write()
    jet_all_gen_eff_den_ungroomed_merged.Write()
    jet_match_rec_pur_num_ungroomed_merged.Write()
    jet_all_rec_pur_den_ungroomed_merged.Write()
    jet_efficiency_ungroomed_merged.Write()
    jet_purity_ungroomed_merged.Write()

    jet_match_gen_eff_num_groomed_merged.Write()
    jet_all_gen_eff_den_groomed_merged.Write()
    jet_match_rec_pur_num_groomed_merged.Write()
    jet_all_rec_pur_den_groomed_merged.Write()
    jet_efficiency_groomed_merged.Write()
    jet_purity_groomed_merged.Write()

    pair_match_gen_eff_num_AA_merged.Write()
    pair_all_gen_eff_den_AA_merged.Write()
    pair_match_rec_pur_num_AA_merged.Write()
    pair_all_rec_pur_den_AA_merged.Write()
    pair_efficiency_AA_merged.Write()
    pair_purity_AA_merged.Write()

    pair_match_gen_eff_num_AB_merged.Write()
    pair_all_gen_eff_den_AB_merged.Write()
    pair_match_rec_pur_num_AB_merged.Write()
    pair_all_rec_pur_den_AB_merged.Write()
    pair_efficiency_AB_merged.Write()
    pair_purity_AB_merged.Write()

    pair_match_gen_eff_num_BB_merged.Write()
    pair_all_gen_eff_den_BB_merged.Write()
    pair_match_rec_pur_num_BB_merged.Write()
    pair_all_rec_pur_den_BB_merged.Write()
    pair_efficiency_BB_merged.Write()
    pair_purity_BB_merged.Write()

    pair_match_gen_eff_num_rad_merged.Write()
    pair_all_gen_eff_den_rad_merged.Write()
    pair_match_rec_pur_num_rad_merged.Write()
    pair_all_rec_pur_den_rad_merged.Write()
    pair_efficiency_rad_merged.Write()
    pair_purity_rad_merged.Write()


    lund_matched_gen_merged.Write()
    lund_matched_rec_merged.Write()
    lund_all_gen_merged.Write()
    lund_all_rec_merged.Write()
    lund_split_efficiency_merged.Write()
    lund_split_purity_merged.Write()
    summary_efficiencies_merged.Write()

    fout.Write()
    fout.Close()

    print("written to {}/response_merged_partial.root".format(dir))
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--dir', default=None)
    parser.add_argument('--i', type=int, default=0)
    parser.add_argument('--j', type=int, default=0)
    parser.add_argument('--final', type=bool, default=False)
    flags = parser.parse_args()

    nonexistent = []
    finishing_merge(flags.dir, flags.i)
    print(nonexistent)
        