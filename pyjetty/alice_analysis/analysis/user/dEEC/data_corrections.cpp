// ROOT macro to make corrections to the dEEC data
// Beatrice Liang-Gilman (beatrice_lg@berkeley.edu)

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

// global variables
Double_t colors[16] = {kGray, kMagenta, kGreen+2, kBlue, kOrange+1, kViolet+1, kRed, kYellow+1, kCyan+1};
Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
Double_t marker_size = 1.5;

std::string attempt_dir = "data_correctionfactors";
std::string outdir = "/software/users/blianggi/mypyjetty/storage/dEEC/plots/" + attempt_dir;




void draw_save_del_hists(TFile *fout, TObject* obj, std::string obs_name, 
                         std::string ptname, std::string norm_string, std::string hist_addname,
                         bool logx, bool logy, bool logz=false) {
    
    TCanvas *can = new TCanvas(); 
    can->cd();
    if (logx) gPad->SetLogx();
    if (logy) gPad->SetLogy();

    if (TH2* hist2D = dynamic_cast<TH2*>(obj)) { // put this first bc TH2 is a subclass of TH1!! (and it will go into the other loop :( )
        gPad->SetRightMargin(0.12);
        if (logz) gPad->SetLogz();
        can->SetFillColor(kWhite);
        hist2D->Draw("COLZ");
    } else if (TH1* hist = dynamic_cast<TH1*>(obj)) {
        hist->Draw();
    } else if (TGraphErrors* graph = dynamic_cast<TGraphErrors*>(obj)) {
        graph->Draw("ALP");
    } else {
        cout << "Error: Unsupported object type. Only TH1, TGraphErrors, and TH2 are supported." << endl;
    }
    
    fout->cd();
    obj->Write(); //TODO: this might not be right! Might have to use the casted type


    std::string add_dir = "";
    if (obs_name != "jet_pt" && obs_name != "total_num_const" && obs_name != "num_const_aftercut") {
        if (obs_name == "rc") add_dir = "/" + ptname + "/" + norm_string + "/" + obs_name;
        else add_dir = "/" + ptname + "/" + norm_string + "/" + obs_name + "/individuals";
    }
    std::string fname_out = outdir + add_dir + "/corrhist_" + obs_name + hist_addname + ".pdf";
    can->SaveAs(fname_out.c_str());

    // delete hist;
    delete can;
}


void save_noncorrected_hists(TFile * f_in, TFile * f_out, std::string weightstr, std::string jetRname, std::string thrname) {
    // Loop through the keys in the file
    TIter nextKey(f_in->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)nextKey())) {
        // Read object from key
        TObject *obj = key->ReadObj();
        // if (obj->InheritsFrom("TH1")) { // Check if object is a histogram
        //     std::cout << "Histogram found: " << obj->GetName() << std::endl;
        // }

        // std::string name = obj->GetName();
        std::string objName = key->GetName();
        // cout << "obj name " << objName << endl;

        // Check if "deltap" is in the name
        if (objName.find("deltap") == std::string::npos && objName.find("weights") == std::string::npos) { //npos = no position found
            std::cout << "Object with 'deltap' and 'weights' in name not found: " << objName << std::endl;
            
            std::string suffix = "_hist";
            if (objName.size() >= suffix.size() && objName.compare(objName.size() - suffix.size(), suffix.size(), suffix) == 0) {
                // Remove the "_hist" part
                std::string obs = objName.substr(0, objName.size() - suffix.size());
                // std::cout << "Extracted string: " << result << std::endl;
                draw_save_del_hists(f_out, obj, obs, "", "", weightstr + jetRname + thrname, false, true);
            
            } else { // graphs
                draw_save_del_hists(f_out, obj, objName, "", "", weightstr + jetRname + thrname, false, false);
            
            }
            delete obj;
        }

    }

    return;
}

TH1D * get_EEC_hist(TFile * f_EEC, int pt_min, int pt_max) {

    std::string EEC_name = "h_jet_ENC_RL2_JetPt_R0.4_1.0"; //"h_jet_EEC_noweight_RL_JetPt_R0.4_1.0";
    std::string jetpt_name = "h_1Djet_pt_JetPt_R0.4_1.0";

    // get the histograms
    TH2D * EEC_jetpt_hist = (TH2D*) f_EEC->Get(EEC_name.c_str());
    TH1D * jetpt_hist = (TH1D*) f_EEC->Get(jetpt_name.c_str());

    // get 1D EEC after making jet pt cuts
    EEC_jetpt_hist->GetXaxis()->SetRangeUser(pt_min, pt_max);
    TH1D * EEC_hist = EEC_jetpt_hist->ProjectionY();
    
    // now normalize by the number of jets
    double num_jets = jetpt_hist->Integral();
    EEC_hist->Scale(num_jets, "width");

    return EEC_hist;

}

double extract_corrfactor_forRLbin(TH1D *f_corr_hist, TH1D * EEC_hist, double RL_min, double RL_max) {
    
    double numerator = 0; 
    double denominator = 0;

    // first get the bins of the histogram that make the RL range.
    int x_left_bin = f_corr_hist->FindBin(RL_min);
    int x_right_bin = f_corr_hist->FindBin(RL_max);

    // then get the scale factors for the RL bins, and weight it by unweighted EEC
    for (int i = x_left_bin; i <= x_right_bin; i++) {
        double f_corr = f_corr_hist->GetBinContent(i);
        double RL = EEC_hist->GetBinContent(i);

        double ind_weight_fcorr = f_corr * RL;
        numerator += ind_weight_fcorr;
        denominator += RL;
    }

    double weighted_fcorr = numerator/denominator;
    cout << "The weighted fcorr for RL=" << RL_min << "-" << RL_max << " is " << weighted_fcorr << endl;

    return weighted_fcorr;

}


void apply_corrfactor(TFile * f_in, TFile * f_out, std::string weightstr, std::string jetRname, std::string thrname,
                      int pt_min, int pt_max, double RL_min, double RL_max, double fcorr_forRLbin) {
    
    std::string ptname = to_string(pt_min) + "-" + to_string(pt_max);
    std::string RLname = Form("_RL%.3f-%.3f", RL_min, RL_max);
    
    std::string norms[] = { "unnormalized", "self_normalized", "norm_by_jets" };
    std::string obss[] = { "deltap", "deltapt", "deltapl", "weights", "weights_vs_deltapt" };

    for (std::string norm : norms) {
        
        std::string hist_addname = weightstr + jetRname + thrname + "_pt" + ptname + RLname + "_" + norm;

        for (std::string obs : obss) {
            std::string hist_name = Form("h_%s_R0.4_t1.0_pt%d-%d_RL%.3f-%.3f_%s",obs.c_str(), pt_min, pt_max, RL_min, RL_max, norm.c_str());
            // cout << "hist name is " << hist_name << endl;

            // get the appropriate histogram, scale it, and write it to file
            if (obs != "weights_vs_deltapt") {
                f_in->cd();
                TH1D* hist1D = (TH1D*)gDirectory->Get(hist_name.c_str());
                hist1D->Scale(fcorr_forRLbin);

                draw_save_del_hists(f_out, hist1D, obs, ptname, norm, hist_addname, false, true);
    
            } else {
                f_in->cd();
                TH2D* hist2D = (TH2D*)gDirectory->Get(hist_name.c_str());
                hist2D->Scale(fcorr_forRLbin);

                draw_save_del_hists(f_out, hist2D, obs, ptname, norm, hist_addname, false, false, true);
    
            }

        }
    }

    return;
    
}

void scale_with_corr_factors(TFile *f_in, TFile *f_out, TFile *f_fcorr, TFile *f_EEC, const int pt_bins[], int n_bins, const double RL_bins[][8],
                             int n_RLbins, std::string weightstr, std::string jetRname, std::string thrname, bool include_RL0, bool include_RL1) {
    
    save_noncorrected_hists(f_in, f_out, weightstr, jetRname, thrname);
    
    for (int i = 0; i < n_bins; i++) {
        cout << "in pt bin" << i << endl;
        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];

        std::string f_corr_hist_name = Form("Ratio_hJetE2C_norm_ch_jetR4_%d%d_trk10", pt_min, pt_max);

        TH1D * fcorr_hist = (TH1D*) f_fcorr->Get(f_corr_hist_name.c_str());
        TH1D * f_EEC_hist = get_EEC_hist(f_EEC, pt_min, pt_max);

        for ( int j = 0; j < n_RLbins; j++ ) {
            int k = j;
            if (!include_RL0) {
                k = j-1;
                if (j == 0) continue; // can add something here to change the filename for ALL
            }
            if (!include_RL1 && j == n_RLbins-1) continue;

            double RL_min = RL_bins[i][j];
            double RL_max = RL_bins[i][j+1];

            double fcorr_forRLbin = extract_corrfactor_forRLbin(fcorr_hist, f_EEC_hist, RL_min, RL_max);
            apply_corrfactor(f_in, f_out, weightstr, jetRname, thrname, pt_min, pt_max, RL_min, RL_max, fcorr_forRLbin);
        }


    }


}



// ========================================================================


void data_corrections() {

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
    TString wenqing_corr_factors_filename = "/software/users/blianggi/mypyjetty/analysis/corr_factors/pythiaMC_corr_hists.root";
    TFile* root_corr_factors_file = new TFile(wenqing_corr_factors_filename, "READ");

    TString input_EEC_filename = "/rstorage/alice/AnalysisResults/blianggi/dEEC/442528/MergedAnalysisResults.root";
    TFile * root_EEC_file = new TFile(input_EEC_filename, "READ");

    // file that needs correcting:
    TString input_histograms_filename = "/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/data_firstattempt/DataHists.root"; //"/software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/analysis/user/dEEC/datahists/DataHists.root";
    TFile* root_infile = new TFile(input_histograms_filename, "READ");

    // Output file with corrected results
    std::string outfile = "/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/" + attempt_dir + "/DataHists.root"; //plots/ntuples/DataHists.root"; //FinalDataHists.root
    TFile* root_outfile = new TFile(outfile.c_str(), "RECREATE");

    scale_with_corr_factors(root_infile, root_outfile, root_corr_factors_file, root_EEC_file, pt_bins, n_bins, RL_bins, n_RLbins, weightstr, jetRname, thrname, include_RL0, include_RL1);
    

    root_corr_factors_file->Close();
    root_EEC_file->Close();
    root_infile->Close();
    root_outfile->Close();

    delete root_corr_factors_file;
    delete root_EEC_file;
    delete root_infile;
    delete root_outfile;


}