import ROOT
import argparse
import os


# python3 -u hadd_manual.py --dir /global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/57259276 --i 27 --final True

is_final_merge = False
rap = False

suffix = "" #"LHC18f3/803"
subdirs = ['0', '1', '2'] # ppb
runlist = ['559348', '559361', '559362', '559385', '559387', '559408', '559409', '559410', '559437', '559443', '559444', '559456']
# 559348 559361 559362 559385 559387 559408 559409 559410 559423 559437 559443 559444 559456

def process(dir, i, j):
    print("LOADING ROOUNFOLD")
    roo = False
    while roo==False:
        roo = ROOT.gSystem.Load("libRooUnfold")
    print(roo)
    roo2 = False
    while roo2==False:
        roo2 = ROOT.gSystem.Load("libRooUnfold.so")
    print(roo2)

    # j, run_i = failed[i]
    run_i = runlist[i]

    fin1 = None
    k = 0
    while fin1==None:
        k += 1
        if k > 42: #6 #change this???
            return
        try:
            fin1 = ROOT.TFile.Open("{}/{}/{}/000{}/AnalysisResults.root".format(dir, subdirs[j], run_i, k)) #subdirs[0]
            print(subdirs[j], run_i, k)
        except:
            pass

    resp_jetpt_merged = fin1.resp_jetpt
    resp_groomed_jetpt_merged = fin1.resp_groomed_jetpt

    resp6_AA_merged = fin1.resp6_AA
    AA_reco_merged = fin1.AA_reco
    AA_reco_unmatched_merged = fin1.AA_reco_unmatched
    AA_gen_merged = fin1.AA_gen
    AA_gen_unmatched_merged = fin1.AA_gen_unmatched

    resp6_AB_merged = fin1.resp6_AB
    AB_reco_merged = fin1.AB_reco
    AB_reco_unmatched_merged = fin1.AB_reco_unmatched
    AB_gen_merged = fin1.AB_gen
    AB_gen_unmatched_merged = fin1.AB_gen_unmatched

    resp6_BB_merged = fin1.resp6_BB
    BB_reco_merged = fin1.BB_reco
    BB_reco_unmatched_merged = fin1.BB_reco_unmatched
    BB_gen_merged = fin1.BB_gen
    BB_gen_unmatched_merged = fin1.BB_gen_unmatched

    resp6_rad_merged = fin1.resp6_rad
    rad_reco_merged = fin1.rad_reco
    rad_reco_unmatched_merged = fin1.rad_reco_unmatched
    rad_gen_merged = fin1.rad_gen
    rad_gen_unmatched_merged = fin1.rad_gen_unmatched



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

    # roounfold matrices
    roounfold_response_1D_merged = fin1.roounfold_response_1D
    AA_roounfold_response_merged = fin1.AA_roounfold_response
    AB_roounfold_response_merged = fin1.AB_roounfold_response
    BB_roounfold_response_merged = fin1.BB_roounfold_response
    rad_roounfold_response_merged = fin1.rad_roounfold_response
    

    print(rad_gen_merged)

    for subdir in subdirs:
    # if subdirs[j]:
        # subdir = subdirs[j]
        
        subruns =  [x[1] for x in os.walk("{}/{}/{}/{}".format(dir, suffix, subdir, run_i))][0] # list of subdirs in run dir

        print(subruns)

        for subrun in subruns:

            if subdir==subdirs[0] and subrun=='000{}'.format(k):
                continue

            try:
                
                fin = ROOT.TFile.Open("{}/{}/{}/{}/AnalysisResults.root".format(dir, subdir, run_i, subrun))

                print(resp_jetpt_merged)
                resp_jetpt_merged.Add(fin.resp_jetpt)
                resp_groomed_jetpt_merged.Add(fin.resp_groomed_jetpt)

                resp6_AA_merged.Add(fin.resp6_AA)
                AA_reco_merged.Add(fin.AA_reco)
                AA_reco_unmatched_merged.Add(fin.AA_reco_unmatched)
                AA_gen_merged.Add(fin.AA_gen)
                AA_gen_unmatched_merged.Add(fin.AA_gen_unmatched)

                resp6_AB_merged.Add(fin.resp6_AB)
                AB_reco_merged.Add(fin.AB_reco)
                AB_reco_unmatched_merged.Add(fin.AB_reco_unmatched)
                AB_gen_merged.Add(fin.AB_gen)
                AB_gen_unmatched_merged.Add(fin.AB_gen_unmatched)

                resp6_BB_merged.Add(fin.resp6_BB)
                BB_reco_merged.Add(fin.BB_reco)
                BB_reco_unmatched_merged.Add(fin.BB_reco_unmatched)
                BB_gen_merged.Add(fin.BB_gen)
                BB_gen_unmatched_merged.Add(fin.BB_gen_unmatched)

                resp6_rad_merged.Add(fin.resp6_rad)
                rad_reco_merged.Add(fin.rad_reco)
                rad_reco_unmatched_merged.Add(fin.rad_reco_unmatched)
                rad_gen_merged.Add(fin.rad_gen)
                rad_gen_unmatched_merged.Add(fin.rad_gen_unmatched)



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

                # roounfold matrices
                roounfold_response_1D_merged.Add(fin.roounfold_response_1D)
                AA_roounfold_response_merged.Add(fin.AA_roounfold_response)
                AB_roounfold_response_merged.Add(fin.AB_roounfold_response)
                BB_roounfold_response_merged.Add(fin.BB_roounfold_response)
                rad_roounfold_response_merged.Add(fin.rad_roounfold_response)

            except Exception as e:
                print("skipping file {}/{}/{}/{}/AnalysisResults.root".format(dir, subdir, run_i, subrun))
                print(f"Error: {e}")
                continue
            
            print("MERGED", subdir, run_i, subrun)

    # fout = ROOT.TFile("{}/merged_{}.root".format(dir, run_i), "RECREATE")
    fout = ROOT.TFile("{}/merged_{}_{}.root".format(dir, subdirs[j], run_i), "RECREATE")

    resp_jetpt_merged.Write()
    resp_groomed_jetpt_merged.Write()

    resp6_AA_merged.Write()
    AA_reco_merged.Write()
    AA_reco_unmatched_merged.Write()
    AA_gen_merged.Write()
    AA_gen_unmatched_merged.Write()

    resp6_AB_merged.Write()
    AB_reco_merged.Write()
    AB_reco_unmatched_merged.Write()
    AB_gen_merged.Write()
    AB_gen_unmatched_merged.Write()

    resp6_BB_merged.Write()
    BB_reco_merged.Write()
    BB_reco_unmatched_merged.Write()
    BB_gen_merged.Write()
    BB_gen_unmatched_merged.Write()

    resp6_rad_merged.Write()
    rad_reco_merged.Write()
    rad_reco_unmatched_merged.Write()
    rad_gen_merged.Write()
    rad_gen_unmatched_merged.Write()



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

    # roounfold matrices
    roounfold_response_1D_merged.Write()
    AA_roounfold_response_merged.Write()
    AB_roounfold_response_merged.Write()
    BB_roounfold_response_merged.Write()
    rad_roounfold_response_merged.Write()
    

    fout.Write()
    fout.Close()

    # print("written to {}/merged_{}.root".format(dir, run_i))
    print("written to {}/merged_{}_{}.root".format(dir, subdirs[j], run_i))
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


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

    resp6_AA_merged = fin1.resp6_AA
    AA_reco_merged = fin1.AA_reco
    AA_reco_unmatched_merged = fin1.AA_reco_unmatched
    AA_gen_merged = fin1.AA_gen
    AA_gen_unmatched_merged = fin1.AA_gen_unmatched

    resp6_AB_merged = fin1.resp6_AB
    AB_reco_merged = fin1.AB_reco
    AB_reco_unmatched_merged = fin1.AB_reco_unmatched
    AB_gen_merged = fin1.AB_gen
    AB_gen_unmatched_merged = fin1.AB_gen_unmatched

    resp6_BB_merged = fin1.resp6_BB
    BB_reco_merged = fin1.BB_reco
    BB_reco_unmatched_merged = fin1.BB_reco_unmatched
    BB_gen_merged = fin1.BB_gen
    BB_gen_unmatched_merged = fin1.BB_gen_unmatched

    resp6_rad_merged = fin1.resp6_rad
    rad_reco_merged = fin1.rad_reco
    rad_reco_unmatched_merged = fin1.rad_reco_unmatched
    rad_gen_merged = fin1.rad_gen
    rad_gen_unmatched_merged = fin1.rad_gen_unmatched



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

    # roounfold matrices
    roounfold_response_1D_merged = fin1.roounfold_response_1D
    AA_roounfold_response_merged = fin1.AA_roounfold_response
    AB_roounfold_response_merged = fin1.AB_roounfold_response
    BB_roounfold_response_merged = fin1.BB_roounfold_response
    rad_roounfold_response_merged = fin1.rad_roounfold_response
    
    
    if not is_final_merge:
        for run in runlist:
            
            # for subdir in subdirs:

                print("files from {}".format(run))

                try:
                    # fin = ROOT.TFile.Open("{}/merged_{}.root".format(dir, run))
                    # fin = ROOT.TFile.Open("{}/merged_{}_{}.root".format(dir, subdirs[0], run))
                    fin = ROOT.TFile.Open("{}/{}/{}/response.root".format(dir, subdirs[0], run))

                    resp_jetpt_merged.Add(fin.resp_jetpt)
                    resp_groomed_jetpt_merged.Add(fin.resp_groomed_jetpt)

                    resp6_AA_merged.Add(fin.resp6_AA)
                    AA_reco_merged.Add(fin.AA_reco)
                    AA_reco_unmatched_merged.Add(fin.AA_reco_unmatched)
                    AA_gen_merged.Add(fin.AA_gen)
                    AA_gen_unmatched_merged.Add(fin.AA_gen_unmatched)

                    resp6_AB_merged.Add(fin.resp6_AB)
                    AB_reco_merged.Add(fin.AB_reco)
                    AB_reco_unmatched_merged.Add(fin.AB_reco_unmatched)
                    AB_gen_merged.Add(fin.AB_gen)
                    AB_gen_unmatched_merged.Add(fin.AB_gen_unmatched)

                    resp6_BB_merged.Add(fin.resp6_BB)
                    BB_reco_merged.Add(fin.BB_reco)
                    BB_reco_unmatched_merged.Add(fin.BB_reco_unmatched)
                    BB_gen_merged.Add(fin.BB_gen)
                    BB_gen_unmatched_merged.Add(fin.BB_gen_unmatched)

                    resp6_rad_merged.Add(fin.resp6_rad)
                    rad_reco_merged.Add(fin.rad_reco)
                    rad_reco_unmatched_merged.Add(fin.rad_reco_unmatched)
                    rad_gen_merged.Add(fin.rad_gen)
                    rad_gen_unmatched_merged.Add(fin.rad_gen_unmatched)



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

                    # roounfold matrices
                    roounfold_response_1D_merged.Add(fin.roounfold_response_1D)
                    AA_roounfold_response_merged.Add(fin.AA_roounfold_response)
                    AB_roounfold_response_merged.Add(fin.AB_roounfold_response)
                    BB_roounfold_response_merged.Add(fin.BB_roounfold_response)
                    rad_roounfold_response_merged.Add(fin.rad_roounfold_response)

                except:
                    print("skipping files from {}, {}".format(subdirs[0], run))
                    # print("skipping files from {}, {}".format(i, run))
                    nonexistent.append((subdirs[0],run))
                    continue
    else:                
        for pthatbin in range(2,21):

            try:
                fin = ROOT.TFile.Open("{}/Stage0/{}/AnalysisResult.root".format(dir,pthatbin))

                resp_jetpt_merged.Add(fin.resp_jetpt)
                resp_groomed_jetpt_merged.Add(fin.resp_groomed_jetpt)

                resp6_AA_merged.Add(fin.resp6_AA)
                AA_reco_merged.Add(fin.AA_reco)
                AA_reco_unmatched_merged.Add(fin.AA_reco_unmatched)
                AA_gen_merged.Add(fin.AA_gen)
                AA_gen_unmatched_merged.Add(fin.AA_gen_unmatched)

                resp6_AB_merged.Add(fin.resp6_AB)
                AB_reco_merged.Add(fin.AB_reco)
                AB_reco_unmatched_merged.Add(fin.AB_reco_unmatched)
                AB_gen_merged.Add(fin.AB_gen)
                AB_gen_unmatched_merged.Add(fin.AB_gen_unmatched)

                resp6_BB_merged.Add(fin.resp6_BB)
                BB_reco_merged.Add(fin.BB_reco)
                BB_reco_unmatched_merged.Add(fin.BB_reco_unmatched)
                BB_gen_merged.Add(fin.BB_gen)
                BB_gen_unmatched_merged.Add(fin.BB_gen_unmatched)

                resp6_rad_merged.Add(fin.resp6_rad)
                rad_reco_merged.Add(fin.rad_reco)
                rad_reco_unmatched_merged.Add(fin.rad_reco_unmatched)
                rad_gen_merged.Add(fin.rad_gen)
                rad_gen_unmatched_merged.Add(fin.rad_gen_unmatched)



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

                # roounfold matrices
                roounfold_response_1D_merged.Add(fin.roounfold_response_1D)
                AA_roounfold_response_merged.Add(fin.AA_roounfold_response)
                AB_roounfold_response_merged.Add(fin.AB_roounfold_response)
                BB_roounfold_response_merged.Add(fin.BB_roounfold_response)
                rad_roounfold_response_merged.Add(fin.rad_roounfold_response)

            except:
                print("skipping files from {}".format(pthatbin))
                continue    

    fout = ROOT.TFile("{}/AnalysisResultFinal.root".format(dir), "RECREATE")

    resp_jetpt_merged.Write()
    resp_groomed_jetpt_merged.Write()

    resp6_AA_merged.Write()
    AA_reco_merged.Write()
    AA_reco_unmatched_merged.Write()
    AA_gen_merged.Write()
    AA_gen_unmatched_merged.Write()

    resp6_AB_merged.Write()
    AB_reco_merged.Write()
    AB_reco_unmatched_merged.Write()
    AB_gen_merged.Write()
    AB_gen_unmatched_merged.Write()

    resp6_BB_merged.Write()
    BB_reco_merged.Write()
    BB_reco_unmatched_merged.Write()
    BB_gen_merged.Write()
    BB_gen_unmatched_merged.Write()

    resp6_rad_merged.Write()
    rad_reco_merged.Write()
    rad_reco_unmatched_merged.Write()
    rad_gen_merged.Write()
    rad_gen_unmatched_merged.Write()



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

    # roounfold matrices
    roounfold_response_1D_merged.Write()
    AA_roounfold_response_merged.Write()
    AB_roounfold_response_merged.Write()
    BB_roounfold_response_merged.Write()
    rad_roounfold_response_merged.Write()

    fout.Write()
    fout.Close()

    print("written to {}/AnalysisResultFinal.root".format(dir))
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--dir', default=None)
    parser.add_argument('--i', type=int, default=0)
    parser.add_argument('--j', type=int, default=0)
    parser.add_argument('--final', type=bool, default=False)
    flags = parser.parse_args()

    if flags.final:
        nonexistent = []
        finishing_merge(flags.dir, flags.i)
        print(nonexistent)
    else:
        process(flags.dir, flags.i, flags.j)
        