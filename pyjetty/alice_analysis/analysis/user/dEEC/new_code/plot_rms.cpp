// This code reads in a file of the dEEC anchored MC and will plot some response matrices.
// specifically: ∆p, ∆jT, EW, RL, jet pt
// From 14-D thnsparse with axes: pt_det, pt_truth, pTRL_det, pTRL_truth, EW_det, EW_truth, delta_p_det, delta_p_truth, delta_pt_det, delta_pt_truth, delta_jt_det, delta_jt_truth, q1q2_det, q1q2_truth

std::string attempt_dir = "RMs/LHC23a3_secondattempt/";

const int pt_bins[] = { 20, 40, 60, 80 };
const int n_bins = sizeof(pt_bins) / sizeof(pt_bins[0]) - 1;
const double ptRL_bins[5] = { 2e-1, 8e-1, 5.0, 10.0, 30.0 };
const int n_ptRLbins = 4;

// Define class Observable
class Observable {
public:
    std::string name;

    int obs_det_axis;
    int obs_tr_axis;

    bool make_pt_cuts;
    bool make_ptrl_cuts;

    std::string filepath_plots;

    Observable(std::string name_val, int obs_det_axis_val, int obs_tr_axis_val, bool make_pt_cuts_val, bool make_ptrl_cuts_val) {
        name = name_val;
        
        obs_det_axis = obs_det_axis_val;
        obs_tr_axis = obs_tr_axis_val; 

        make_pt_cuts = make_pt_cuts_val;
        make_ptrl_cuts = make_ptrl_cuts_val; 

        filepath_plots = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/plots/" + attempt_dir + name;
    }
};


// Selecting the pairs from the TRUTH jet pt and TRUTH ptRL
TH2D * make_cuts_and_project(THnSparse * hsparse, Observable obs, int pt_min, int pt_max, double ptrl_min, double ptrl_max) {
    // Cut on pT and pTRL
    if (obs.make_pt_cuts) hsparse->GetAxis(1)->SetRangeUser(pt_min, pt_max);
    if (obs.make_ptrl_cuts) hsparse->GetAxis(3)->SetRangeUser(ptrl_min, ptrl_max);

    // Project onto the appropriate axes
    TH2D * hRM = hsparse->Projection(obs.obs_det_axis, obs.obs_tr_axis);
    return hRM;
}

// plot the 2D histogram on a canvas
void plot_2D_hist(TH2D * hist2D, Observable obs, std::string kin_name, TFile * fout) {
    TCanvas * can = new TCanvas();
    can->cd();
    gPad->SetLogz();
    if ( obs.name == "pTRL" ) {
        gPad->SetLogx();
        gPad->SetLogy();
    }
    hist2D->Draw("COLZ");

    // TODO: label the plot here!!

    can->SaveAs(Form("%s/%s_RM%s.pdf", obs.filepath_plots.c_str(), obs.name.c_str(), kin_name.c_str()));

    // save 
    fout->cd();
    hist2D->Write();
}


void get_observable_RMs(THnSparse * hsparse, Observable obs, TFile * fout) {

    for ( int i = 0; i < n_bins; i++ ) {
        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];

        for ( int j = 0; j < n_ptRLbins; j++ ) {
            double ptrl_min = ptRL_bins[j];
            double ptrl_max = ptRL_bins[j+1];
            std::string kin_name = "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_ptrlbin" + std::to_string(j);
            if ( obs.name == "jet_pt" ) kin_name = "";
            if ( obs.name == "pTRL" ) kin_name = "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);

            TH2D * hist2D = make_cuts_and_project(hsparse, obs, pt_min, pt_max, ptrl_min, ptrl_max);
            std::string hname = Form("hRM_%s_pt%d-%d_ptrlbin%d", obs.name.c_str(), pt_min, pt_max, j);
            if ( obs.name == "jet_pt" && j == 0 ) hname = Form("hRM_%s", obs.name.c_str());
            if ( obs.name == "pTRL" && j == 0 ) hname = Form("hRM_%s_pt%d-%d", obs.name.c_str(), pt_min, pt_max);
            hist2D->SetNameTitle(hname.c_str(), hname.c_str());
            plot_2D_hist(hist2D, obs, kin_name, fout);

            if ( obs.name == "jet_pt" && j == 0 ) return; // just do jet pt once
            if ( obs.name == "pTRL" && j == 0 ) break; // just do pTRL once for each jet pt bin
        }
    }

}

void plot_rms() {

    gStyle->SetOptStat(0);

    // Open the file
    TString anchmc_filename = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/50801801/scaling/AnalysisResultsFinal.root"; 
    // TString anchmc_filename = "/global/cfs/cdirs/alice/blianggi/mypyjetty/analysis/testing/AnalysisResults.root"; //testing file
    TFile* anchmc_file = new TFile(anchmc_filename, "READ");

    // Get the 14-D histogram ahahaha
    THnSparse * hsparse = (THnSparse * ) anchmc_file->Get("hResponse_14D_ALLOBS_JetPt_R0.4_1.0");

    // Define observables
    Observable obs_jetpt("jet_pt", 0, 1, false, false);
    Observable obs_pTRL("pTRL", 2, 3, true, false);
    Observable obs_deltap("deltap", 6, 7, true, true);
    Observable obs_deltajt("deltajt", 10, 11, true, true);
    Observable obs_weights("weights", 4, 5, true, true);

    vector<Observable> obs_list = { obs_jetpt, obs_pTRL, obs_deltap, obs_deltajt, obs_weights };

    // Define output filename
    std::string filename_out = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/rootfiles/" + attempt_dir + "RMHists.root";
    TFile * out_file = new TFile(filename_out.c_str(), "RECREATE");
    
    // Loop through observables
    for ( Observable obs : obs_list) {
        get_observable_RMs(hsparse, obs, out_file);
    }

    out_file->Close();
    
}