import ROOT
import argparse
import os


# python3 -u hadd_manual_ptdifferential.py --dir /global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/57626568

is_final_merge = False
rap = False

# suffix = "" #"LHC18f3/803"
# subdirs = ['0', '1', '2'] # ppb
# runlist = ['559348', '559361', '559362', '559385', '559387', '559408', '559409', '559410', '559437', '559443', '559444', '559456']


default_subpath_file = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/anchmc_subpath_filelist.txt"


def read_subpaths(filelist):
    out = []
    with open(filelist) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                out.append(line.strip('/'))
    return out
    

def finishing_merge(dir):

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
    fin1 = ROOT.TFile.Open("{}/{}/response_ptdifferential.root".format(dir, subpaths[0]))

    # Dynamically find all histogram names in the file
    hist_list = [key.GetName() for key in fin1.GetListOfKeys() if "TH" in key.GetClassName()]
    print(f"Found {len(hist_list)} histograms to merge.")

    merged_hists = {h: getattr(fin1, h) for h in hist_list}


    for ifile, sp in enumerate(subpaths):
        # Skip first file because already used to initialize merged_hists
        if sp == subpaths[0]:
            continue
        if ifile % 500 == 0:
            print(f"File {ifile}: {sp}")
        try:
            fin = ROOT.TFile.Open("{}/{}/response_ptdifferential.root".format(dir, sp))

            for h_name in hist_list:
                hist = getattr(fin, h_name, None)
                if hist:
                    merged_hists[h_name].Add(hist)

        except:
            print("skipping files from {}".format(sp))
            # print("skipping files from {}".format(sp))
            nonexistent.append(sp)
            continue

    fout = ROOT.TFile("{}/response_ptdifferential_merged.root".format(dir), "RECREATE")

    # RECALCULATE EFF AND PURITY
    # Define the pairs for recalculation: (num_key, den_key, output_name)
    recalc_pairs = [
        ('h_lund_matched_gen_{}', 'h_lund_all_gen_{}', 'h_lund_split_efficiency_{}_new'),
        ('h_lund_matched_rec_{}', 'h_lund_all_rec_{}', 'h_lund_split_purity_{}_new'),
        ('pair_match_gen_eff_num_{}_full_ungroomed', 'pair_all_gen_eff_den_{}_full_ungroomed', 'pair_efficiency_full_ungroomed_{}_new'),
        ('pair_match_rec_pur_num_{}_full_ungroomed', 'pair_all_rec_pur_den_{}_full_ungroomed', 'pair_purity_full_ungroomed_{}_new'),
        ('pair_match_gen_eff_num_{}_AA', 'pair_all_gen_eff_den_{}_AA', 'pair_efficiency_AA_{}_new'),
        ('pair_match_rec_pur_num_{}_AA', 'pair_all_rec_pur_den_{}_AA', 'pair_purity_AA_{}_new'),
        ('pair_match_gen_eff_num_{}_AB', 'pair_all_gen_eff_den_{}_AB', 'pair_efficiency_AB_{}_new'),
        ('pair_match_rec_pur_num_{}_AB', 'pair_all_rec_pur_den_{}_AB', 'pair_purity_AB_{}_new'),
        ('pair_match_gen_eff_num_{}_BB', 'pair_all_gen_eff_den_{}_BB', 'pair_efficiency_BB_{}_new'),
        ('pair_match_rec_pur_num_{}_BB', 'pair_all_rec_pur_den_{}_BB', 'pair_purity_BB_{}_new'),
        ('pair_match_gen_eff_num_{}_rad', 'pair_all_gen_eff_den_{}_rad', 'pair_efficiency_rad_{}_new'),
        ('pair_match_rec_pur_num_{}_rad', 'pair_all_rec_pur_den_{}_rad', 'pair_purity_rad_{}_new'),
    ]

    PTBINS = [10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 150.0, 200.0, 500.0]
    for i in range(len(PTBINS)-1):
        ptstr = f"pt{PTBINS[i]}-{PTBINS[i+1]}"

        for num_k, den_k, out_name in recalc_pairs:
            full_num_k = num_k.format(ptstr)
            full_den_k = den_k.format(ptstr)
            full_out_name = out_name.format(ptstr)
            if full_num_k in merged_hists and full_den_k in merged_hists:
                h_new = merged_hists[full_num_k].Clone(full_out_name)
                h_new.Divide(merged_hists[full_den_k])
                h_new.Write()
            else:
                print(f"Warning: Missing histograms for {full_out_name} ({full_num_k} or {full_den_k}). Skipping.")


    # WRITE FILES
    for h_name in hist_list:
        merged_hists[h_name].Write()

    fout.Write()
    fout.Close()

    print("written to {}/response_ptdifferential_merged.root".format(dir))
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--dir', default=None)
    parser.add_argument('--final', type=bool, default=False)
    flags = parser.parse_args()

    nonexistent = []
    finishing_merge(flags.dir)
    # print(nonexistent)
        