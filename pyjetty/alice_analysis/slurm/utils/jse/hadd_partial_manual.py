import ROOT
import argparse
import os


# python3 -u hadd_partial_manual.py --dir /global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/57485566

is_final_merge = False
rap = False

suffix = "" #"LHC18f3/803"
subdirs = ['0', '1', '2'] # ppb
runlist = ['559348', '559361', '559362', '559385', '559387', '559408', '559409', '559410', '559437', '559443', '559444', '559456']

# hist_list = ['resp_jetpt', 'resp_groomed_jetpt', 
#              'jet_match_gen_eff_num_ungroomed', 'jet_all_gen_eff_den_ungroomed', 'jet_match_rec_pur_num_ungroomed', 'jet_all_rec_pur_den_ungroomed', 
#              'jet_match_gen_eff_num_groomed', 'jet_all_gen_eff_den_groomed', 'jet_match_rec_pur_num_groomed', 'jet_all_rec_pur_den_groomed', 
#              'resp6_AA', 'resp6_AB', 'resp6_BB', 'resp6_rad',
#              'pair_match_gen_eff_num_AA', 'pair_all_gen_eff_den_AA', 'pair_match_rec_pur_num_AA', 'pair_all_rec_pur_den_AA', 
#              'pair_match_gen_eff_num_AB', 'pair_all_gen_eff_den_AB', 'pair_match_rec_pur_num_AB', 'pair_all_rec_pur_den_AB', 
#              'pair_match_gen_eff_num_BB', 'pair_all_gen_eff_den_BB', 'pair_match_rec_pur_num_BB', 'pair_all_rec_pur_den_BB', 
#              'pair_match_gen_eff_num_rad', 'pair_all_gen_eff_den_rad', 'pair_match_rec_pur_num_rad', 'pair_all_rec_pur_den_rad', 
#              'lund_matched_gen', 'lund_matched_rec', 'lund_all_gen', 'lund_all_rec', 'summary_efficiencies']
#groomed_full_ungroomed is just the full ungroomed in bins of groomed jet pt
hist_list = ['resp_jetpt_{}', 'resp6_{}_full_ungroomed',
             'resp6_{}_AA', 'resp6_{}_AB', 'resp6_{}_BB', 'resp6_{}_rad',
             'jet_match_gen_eff_num_{}', 'jet_all_gen_eff_den_{}', 'jet_match_rec_pur_num_{}', 'jet_all_rec_pur_den_{}', 
             'lund_matched_gen', 'lund_matched_rec', 'lund_all_gen', 'lund_all_rec', 
             'pair_match_gen_eff_num_full_ungroomed', 'pair_all_gen_eff_den_full_ungroomed', 'pair_match_rec_pur_num_full_ungroomed', 'pair_all_rec_pur_den_full_ungroomed', 
             'pair_match_gen_eff_num_AA', 'pair_all_gen_eff_den_AA', 'pair_match_rec_pur_num_AA', 'pair_all_rec_pur_den_AA', 
             'pair_match_gen_eff_num_AB', 'pair_all_gen_eff_den_AB', 'pair_match_rec_pur_num_AB', 'pair_all_rec_pur_den_AB', 
             'pair_match_gen_eff_num_BB', 'pair_all_gen_eff_den_BB', 'pair_match_rec_pur_num_BB', 'pair_all_rec_pur_den_BB', 
             'pair_match_gen_eff_num_rad', 'pair_all_gen_eff_den_rad', 'pair_match_rec_pur_num_rad', 'pair_all_rec_pur_den_rad', 
             'summary_efficiencies']

default_subpath_file = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/anchmc_subpath_filelist.txt"


def read_subpaths(filelist):
    out = []
    with open(filelist) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                out.append(line.strip('/'))
    return out
    

def finishing_merge(dir, grooming_str):

    print("LOADING ROOUNFOLD")
    roo = False
    while roo==False:
        roo = ROOT.gSystem.Load("libRooUnfold")
    print(roo)
    roo2 = False
    while roo2==False:
        roo2 = ROOT.gSystem.Load("libRooUnfold.so")
    print(roo2)

    # Open subpath file
    subpaths = read_subpaths(default_subpath_file)

    # fin1 = ROOT.TFile.Open("{}/{}/{}/0/response.root".format(dir, subdirs[0], runlist[0]))
    fin1 = ROOT.TFile.Open("{}/{}/response.root".format(dir, subpaths[0]))

    if grooming_str == "ungroomed":
        hist_list = ['resp_jetpt_{}', 'resp6_{}_full_ungroomed',
                     'resp6_{}_AA', 'resp6_{}_AB', 'resp6_{}_BB', 'resp6_{}_rad',
                     'jet_match_gen_eff_num_{}', 'jet_all_gen_eff_den_{}', 'jet_match_rec_pur_num_{}', 'jet_all_rec_pur_den_{}', 
                     'summary_efficiencies']

    # merged_hists = {h: getattr(fin1, h) for h in hist_list}
    merged_hists = {(h.format(grooming_str) if "{}" in h else h): getattr(fin1, h.format(grooming_str) if "{}" in h else h) for h in hist_list}

    
    for ifile, sp in enumerate(subpaths):
        # Skip first file because already used to initialize merged_hists
        if sp == subpaths[0]:
            continue
        if ifile % 500 == 0:
            print(f"File {ifile}: {sp}")
        try:
            # fin = ROOT.TFile.Open("{}/merged_{}.root".format(dir, run))
            # fin = ROOT.TFile.Open("{}/merged_{}_{}.root".format(dir, subdirs[0], run))
            # fin = ROOT.TFile.Open("{}/{}/{}/{}/response.root".format(dir, subdir, run, k))
            fin = ROOT.TFile.Open("{}/{}/response.root".format(dir, sp))
            # current = [os.path.join(indir, sp, "response.root") for sp in subpaths]

            for h_name in hist_list:
                hist = getattr(fin, h_name, None)
                if hist:
                    merged_hists[h_name].Add(hist)

        except:
            print("skipping files from {}".format(sp))
            # print("skipping files from {}".format(sp))
            nonexistent.append(sp)
            continue

    fout = ROOT.TFile("{}/response_{}_merged_partial.root".format(dir, grooming_str), "RECREATE")

    # RECALCULATE EFF AND PURITY
    jet_eff_new = merged_hists[f'jet_match_gen_eff_num_{grooming_str}'].Clone(f"jet_efficiency_{grooming_str}_new")
    jet_eff_new.Divide(merged_hists[f'jet_all_gen_eff_den_{grooming_str}'])
    jet_pur_new = merged_hists[f'jet_match_rec_pur_num_{grooming_str}'].Clone(f"jet_purity_{grooming_str}_new")
    jet_pur_new.Divide(merged_hists[f'jet_all_rec_pur_den_{grooming_str}'])

    if grooming_str == "groomed":
        lund_split_efficiency_new = merged_hists['lund_matched_gen'].Clone("lund_split_efficiency_new")
        lund_split_efficiency_new.Divide(merged_hists['lund_all_gen'])
        lund_split_purity_new = merged_hists['lund_matched_rec'].Clone("lund_split_purity_new")
        lund_split_purity_new.Divide(merged_hists['lund_all_rec'])

        pair_eff_full_ungroomed_new = merged_hists['pair_match_gen_eff_num_full_ungroomed'].Clone("pair_efficiency_full_ungroomed_new")
        pair_eff_full_ungroomed_new.Divide(merged_hists['pair_all_gen_eff_den_full_ungroomed'])
        pair_pur_full_ungroomed_new = merged_hists['pair_match_rec_pur_num_full_ungroomed'].Clone("pair_purity_full_ungroomed_new")
        pair_pur_full_ungroomed_new.Divide(merged_hists['pair_all_rec_pur_den_full_ungroomed'])
        
        pair_eff_AA_new = merged_hists['pair_match_gen_eff_num_AA'].Clone("pair_efficiency_AA_new")
        pair_eff_AA_new.Divide(merged_hists['pair_all_gen_eff_den_AA'])
        pair_pur_AA_new = merged_hists['pair_match_rec_pur_num_AA'].Clone("pair_purity_AA_new")
        pair_pur_AA_new.Divide(merged_hists['pair_all_rec_pur_den_AA'])

        pair_eff_AB_new = merged_hists['pair_match_gen_eff_num_AB'].Clone("pair_efficiency_AB_new")
        pair_eff_AB_new.Divide(merged_hists['pair_all_gen_eff_den_AB'])
        pair_pur_AB_new = merged_hists['pair_match_rec_pur_num_AB'].Clone("pair_purity_AB_new")
        pair_pur_AB_new.Divide(merged_hists['pair_all_rec_pur_den_AB'])

        pair_eff_BB_new = merged_hists['pair_match_gen_eff_num_BB'].Clone("pair_efficiency_BB_new")
        pair_eff_BB_new.Divide(merged_hists['pair_all_gen_eff_den_BB'])
        pair_pur_BB_new = merged_hists['pair_match_rec_pur_num_BB'].Clone("pair_purity_BB_new")
        pair_pur_BB_new.Divide(merged_hists['pair_all_rec_pur_den_BB'])

        pair_eff_rad_new = merged_hists['pair_match_gen_eff_num_rad'].Clone("pair_efficiency_rad_new")
        pair_eff_rad_new.Divide(merged_hists['pair_all_gen_eff_den_rad'])
        pair_pur_rad_new = merged_hists['pair_match_rec_pur_num_rad'].Clone("pair_purity_rad_new")
        pair_pur_rad_new.Divide(merged_hists['pair_all_rec_pur_den_rad'])


    jet_eff_new.Write()
    jet_pur_new.Write()

    if grooming_str == "groomed":
        lund_split_efficiency_new.Write()
        lund_split_purity_new.Write()

        pair_eff_AA_new.Write()
        pair_pur_AA_new.Write()
        pair_eff_AB_new.Write()
        pair_pur_AB_new.Write()
        pair_eff_BB_new.Write()
        pair_pur_BB_new.Write()
        pair_eff_rad_new.Write()
        pair_pur_rad_new.Write()



    # WRITE FILES
    for h_name in hist_list:
        full_h_name = h_name.format(grooming_str) if "{}" in h_name else h_name
        merged_hists[full_h_name].Write()

    fout.Write()
    fout.Close()

    print("written to {}/response_{}_merged_partial.root".format(dir, grooming_str))
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--dir', default=None)
    parser.add_argument('--final', type=bool, default=False)
    flags = parser.parse_args()

    nonexistent = []
    finishing_merge(flags.dir, "groomed")
    # print(nonexistent)
    finishing_merge(flags.dir, "ungroomed")
        