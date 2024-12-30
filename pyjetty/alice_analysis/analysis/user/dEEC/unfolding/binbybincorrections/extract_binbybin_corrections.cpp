// This file finds the bin by bin corrections and applies them to data.
// This is to be run on perlmutter.
// outputted plots:
    // truth_and_det_[obs].pdf
    // ratio_det_truth_[obs].pdf
    // corr_data_[obs].pdf
    // corr_data_and_raw_data_[obs].pdf
// Beatrice Liang-Gilman, beatrice_lg@berkeley.edu


std::string attempt_dir = "binbybincorrections";


//===========================================================================
//================================ PLOTTING =================================
//===========================================================================

// filetype 1 = ratio_det_truth_[obs].pdf
// filetype 2 = corr_data_[obs].pdf
void plot_and_save_one_histogram(TCanvas * can, TFile * file, std::string obsname,
                                            TH1D * hist1, int filetype, std::string addname) {

    can->cd();
    // gPad->SetLogy();
    hist1->Draw("same");

    std::string filename = "";
    if (filetype == 1) filename = "ratio_det_truth_";
    else if (filetype == 2) filename = "corr_data_";
    std::string outputbase = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/plots/";
    std::string outputname = outputbase + attempt_dir + "/" + obsname + "/" + filename + obsname + addname + ".pdf";
    can->SaveAs(outputname.c_str());

    file->cd();
    hist1->Write();

    delete can;

}


// filetype 1 = truth_and_det_[obs].pdf
// filetype 2 = corr_data_and_raw_data_[obs].pdf
void plot_and_save_two_histograms_overlayed(TCanvas * can, TFile * file, std::string obsname,
                                            TH1D * hist1, TH1D * hist2, int filetype, std::string addname) {

    can->cd();
    gPad->SetLogy();
    hist1->Draw("same");
    hist2->Draw("same");

    std::string filename = "";
    if (filetype == 1) filename = "truth_and_det_";
    else if (filetype == 2) filename = "corr_data_and_raw_data_";
    std::string outputbase = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/plots/";
    std::string outputname = outputbase + attempt_dir + "/" + obsname + "/" + filename + obsname + addname + ".pdf";
    can->SaveAs(outputname.c_str());

    file->cd();
    hist1->Write();
    hist2->Write();

    delete can;

}


//===========================================================================
//============================= OTHER FUNCTIONS =============================
//===========================================================================

//===========================================================================
// This function gets the reco (det) level observable and the gen (truth) level 
// observable, and takes the ratio of the two.
// It plots this ratio and saves it as a separate histogram.
//===========================================================================
TH1D * get_binbybin_corrfactors(TFile *fin_mc, TFile *fout, std::string observable, 
                                int ptbin, int pt_min, int pt_max, int rlbin) {

    fin_mc->cd();

    std::string histname = Form("hResponse_JetPt_corr_%s_PTBIN%d_RLBIN%d_R0.4_1.0", observable.c_str(), ptbin, rlbin);
    cout << " HISTNAME: " << histname << endl;
    THnSparse *hsparse = (THnSparse *) fin_mc->Get(histname.c_str());

    // extract 1D histograms
    hsparse->GetAxis(0)->SetRangeUser(pt_min, pt_max);
    TH1D * hist_det = (TH1D *) hsparse->Projection(2);

    hsparse->GetAxis(1)->SetRangeUser(pt_min, pt_max);
    TH1D * hist_truth = (TH1D *) hsparse->Projection(3);

    // add name string
    std::string addname = Form("_PTBIN%d_RLBIN%d", ptbin, rlbin);

    // plot and save those histograms
    TCanvas *can_truth_det = new TCanvas();
    plot_and_save_two_histograms_overlayed(can_truth_det, fout, observable, hist_det, hist_truth, 1, addname);


    // find the bin by bin corrections
    TH1D * hratio = (TH1D *) hist_det->Clone(Form("hratio_%s_PTBIN%d_RLBIN%d", observable.c_str(), ptbin, rlbin));
    hratio->Divide(hist_truth);

    // plot and save ratio
    TCanvas *can_ratio = new TCanvas();
    plot_and_save_one_histogram(can_ratio, fout, observable, hratio, 1, addname);

    return hratio;

}

// Get the raw data
TH1D * get_rawdata(TFile *fin_data, std::string observable, int pt_min, int pt_max, double RL_min, double RL_max) {

    fin_data->cd();
    std::string histname = Form("h_%s_R0.4_t1.0_pt%d-%d_RL%.3f-%.3f_norm_by_jets", observable.c_str(), pt_min, pt_max, RL_min, RL_max);
    TH1D * hist = (TH1D *)gDirectory->Get(histname.c_str());

    return hist;

}

// Correct the raw data here
// Important!! Right now the bins need to be the same; otherwise need to rebin to make it match
void apply_corrfactor(TFile *fout, TH1D * h_raw_data, TH1D * h_binbybin_fcorr) {

    TH1D *h_corr_data = (TH1D*) h_raw_data->Clone(h_raw_data->GetName());
    h_corr_data->Divide(h_binbybin_fcorr);
}

// General analysis function
void analyze(TFile *f_in_data, TFile *f_in_mc, TFile *f_out, const int pt_bins[], int n_bins, const double RL_bins[][8],
             int n_RLbins, std::string weightstr, std::string jetRname, std::string thrname, bool include_RL0, bool include_RL1) {
    
    // do i need this?
    // save_noncorrected_hists(f_in, f_out, weightstr, jetRname, thrname);
    
    for (int i = 0; i < n_bins; i++) {
        cout << "in pt bin" << i << endl;
        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];


        for ( int j = 0; j < n_RLbins; j++ ) {
            int k = j;
            if (!include_RL0) {
                k = j-1;
                if (j == 0) continue; // can add something here to change the filename for ALL
            }
            if (!include_RL1 && j == n_RLbins-1) continue;

            double RL_min = RL_bins[i][j];
            double RL_max = RL_bins[i][j+1];

            
            // double fcorr_forRLbin = extract_corrfactor_forRLbin(fcorr_hist, f_EEC_hist, RL_min, RL_max);
            // apply_corrfactor(f_in, f_out, weightstr, jetRname, thrname, pt_min, pt_max, RL_min, RL_max, fcorr_forRLbin);
        
            TH1D * h_binbybin_fcorr = get_binbybin_corrfactors(f_in_mc, f_out, "deltap", i, pt_min, pt_max, k);
            TH1D * h_raw_data = get_rawdata(f_in_data, "deltap", pt_min, pt_max, RL_min, RL_max);
            // apply_corrfactor();
        
        }



    }


}

// //===========================================================================
// // This function extracts the thnsparse response matrix, changes the name, 
// // and saves it to its own file in the correct directory.
// //===========================================================================
// void extract_2DRM(TFile *fin, TFile *fout, std::string observable, int ptbin, int RLbin) {

//     std::string original_histname = Form("hResponse_JetPt_%s_PTBIN%d_RLBIN%d_R0.4_1.0", observable, ptbin, RLbin);

//     fin->cd();
//     THnSparse * hist = (THnSparse *)gDirectory->Get(original_histname.c_str());

//     // change the name - { hResponse_JetPt_[obs]_R[R]_[subobs]_[grooming setting] }
//     std::string new_histname = Form("hResponse_JetPt_%s_R0.4_1.0", observable);
//     hist->SetNameTitle();

//     fout->cd();
//     hist->Write();

// }

//===========================================================================
// This file takes in the following arguments:
// [[nothing]]
// THIS IS GOOD FOR PERLMUTTER
//===========================================================================
void extract_binbybin_corrections() { 

    bool include_RL0 = false;
    bool include_RL1 = false;

    std::string weightstr = ""; 
    std::string jetRname = "_R0.4"; 
    std::string thrname = "_t1.0";

    // analysis variables
    const int pt_bins[] = { 20, 40, 60, 80 };
    const int n_bins = sizeof(pt_bins) / sizeof(pt_bins[0]) - 1; //3;
    
    const double RL_bins[3][8] = { { 0, 1e-2, 3e-2, 7e-2, 1.5e-1, 3e-1, 4e-1, 1 },
                            { 0, 1e-2, 2.5e-2, 4e-2, 8e-2, 2.5e-1, 4e-1, 1 },
                            { 0, 1e-2, 2.5e-2, 3e-2, 4.5e-2, 2e-1, 4e-1, 1 } };
    const int n_RLbins = sizeof(RL_bins[0]) / sizeof(RL_bins[0][0]) - 1; //gets the columns //6; //7; //5;

    // filenames
    // file that needs correcting:
    // TString input_histograms_filename = "/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/data_firstattempt/DataHists.root"; //"/software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/analysis/user/dEEC/datahists/DataHists.root";
    TString input_histograms_filename = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/rootfiles/data_firstattempt/DataHists.root"; //"/software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/analysis/user/dEEC/datahists/DataHists.root";
    TFile* root_data_file = new TFile(input_histograms_filename, "READ");

    // file with anchored mc - truth vs det level information
    // TString input_mc_filename = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/34225534/scaling/AnalysisResultsFinal.root";
    TString input_mc_filename = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/34225534/scaling/20/AnalysisResults.root";
    TFile* root_mc_file = new TFile(input_mc_filename, "READ");

    // Output file with corrected results
    // std::string outfile = "/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/" + attempt_dir + "/DataHists.root";
    std::string outfile = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/rootfiles/" + attempt_dir + "/DataHists_BinByBinCorr.root"; 
    TFile* root_outfile = new TFile(outfile.c_str(), "RECREATE");

    analyze(root_data_file, root_mc_file, root_outfile, pt_bins, n_bins, RL_bins, n_RLbins, weightstr, jetRname, thrname, include_RL0, include_RL1);
    

    root_data_file->Close();
    root_mc_file->Close();
    root_outfile->Close();

    delete root_data_file;
    delete root_mc_file;
    delete root_outfile;


}