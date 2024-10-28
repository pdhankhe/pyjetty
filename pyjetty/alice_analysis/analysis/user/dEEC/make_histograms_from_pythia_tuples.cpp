// ROOT macro to take correlation tuples and turn them into histograms
// Beatrice Liang-Gilman (beatrice_lg@berkeley.edu)


using namespace std;

void SetStyle(Bool_t graypalette=true) {
    cout << "Setting style!" << endl;
  
    gStyle->Reset("Plain");
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    if(graypalette) gStyle->SetPalette(8,0);
    else gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetPadColor(10);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetHistLineWidth(1);
    gStyle->SetHistLineColor(kRed);
    gStyle->SetFuncWidth(2);
    gStyle->SetFuncColor(kGreen);
    gStyle->SetLineWidth(1);
    gStyle->SetLabelSize(0.045,"xyz");
    gStyle->SetLabelOffset(0.005,"y"); //(0.01,"y");
    gStyle->SetLabelOffset(0.005,"x"); //(0.01,"x");
    gStyle->SetLabelColor(kBlack,"xyz");
    gStyle->SetTitleSize(0.05,"xyz");
    gStyle->SetTitleOffset(1.25,"y");
    gStyle->SetTitleOffset(1.2,"x");
    gStyle->SetTitleFillColor(kWhite);
    gStyle->SetTextSizePixels(26);
    gStyle->SetTextFont(42);
    //gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y");
    
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(kWhite);
    //gStyle->SetFillColor(kWhite);
    gStyle->SetLegendFont(42);

}

TH1D * getObs1DHistFromTChain(TChain *chain, std::string branch_name, int num_bins, double hist_xmin, double hist_xmax,
                              int pt_min, int pt_max, double RL_min, double RL_max) {
    
    TH1D * hist1D = new TH1D(Form("%s_hist", branch_name.c_str()), Form("%s_hist", branch_name.c_str()), num_bins, hist_xmin, hist_xmax);
    chain->Draw(Form("%s>>%s_hist", branch_name.c_str(), branch_name.c_str()), Form("jet_pt >= %d && jet_pt < %d && RL >= %f && RL < %f", pt_min, pt_max, RL_min, RL_max), "e");
//    chain->Draw(Form("%s>>%s_hist", branch_name.c_str(), branch_name.c_str()), Form("jet_pt >= 20 && jet_pt < 40 && RL > 0.01 && RL < 0.4"), "e");

    return hist1D;
}

/* Save and delete histograms */
void draw_save_del_hists(TCanvas *can, TH1D* hist, std::string obs_name, std::string hist_addname) {
    can->cd();
    hist->Draw();
    
    std::string outdir = ""; //plots/fourthattempt/"; // + ptbin_name + "/";//"plots/test/";
    std::string fname_out = outdir + "testing/corrhist_" + obs_name + hist_addname + ".pdf";
    can->SaveAs(fname_out.c_str());

    delete hist;
    delete can;
}




// ======================================================= //
//                     MAIN FUNCTION
// ======================================================= //

void analyze_data_tuples() {
    gStyle->SetOptStat(0);
    SetStyle();
    
    // setup variables
    bool debug = true;
    bool debug2 = false;
    
    // ntuple names
    std::string JETINFO_name = "h_JETINFOjet_pt_JetPt_R0.4_1.0";
    std::string PAIRINFO_name = "tn_pairlevel_R0.4_1.0";
    
    // other names
    std::string jet1D_name = "h_1Djet_pt_JetPt_R0.4_1.0"; // this one is a histogram
    std::string deltap_str = "deltap";
    
    // filenames
    std::string filename = Form("~/Documents/research/othercorrelations/data_ntuples/AnalysisResults_0001.root");
    // TODO: FIX THIS! (below)
    std::string base_filepath = Form("/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/31843529");
        
    
    // Output file for binned results
    std::string root_outfile = "testing/DataHists.root"; //plots/ntuples/DataHists.root"; //FinalDataHists.root
    TFile* f_out = new TFile(root_outfile.c_str(), "RECREATE");
    std::string add_name = ""; //"_othercorrel";

    
    /*------------------------------------------------------------
    //----------------------- NTUPLE INFO ------------------------
    JETINFO: jet_pt; total_num_const; num_const_aftercut; total_num_baryons; num_baryons_aftercut; total_num_mesons; num_mesons_aftercut
    PAIRINFO: jet_pt; RL; weights; deltap; p1; p2; deltapt; pt1; pt2; deltapl; pl1; pl2; q1q2; q1; q2; baryonmeson; pid1; pid2
    baryon: jet_pt; baryon_pt
    meson: jet_pt; meson_pt
    //----------------------------------------------------------*/
    
    // analysis variables
    const int pt_bins[] = { 20, 40, 60, 80 };
    const int n_bins = sizeof(pt_bins) / sizeof(pt_bins[0]) - 1; //3;
//    std::string pt_bin_names[n_bins];
//    for (int ipt = 0; ipt < n_bins; ipt++) {
//        pt_bin_names[ipt] = std::to_string(pt_bins[ipt]) + "-" + std::to_string(pt_bins[ipt+1]);
//    }
    
    const double RL_bins[3][8] = { { 0, 1e-2, 3e-2, 7e-2, 1.5e-1, 3e-1, 4e-1, 1 },
                            { 0, 1e-2, 2.5e-2, 4e-2, 8e-2, 2.5e-1, 4e-1, 1 },
                            { 0, 1e-2, 2.5e-2, 3e-2, 4.5e-2, 2e-1, 4e-1, 1 } };
    const int n_RLbins = sizeof(RL_bins[0]) / sizeof(RL_bins[0][0]) - 1; //gets the columns //6; //7; //5;
    
    if (debug2) cout << "pt_bins " << n_bins << " n_RLbins " << n_RLbins << endl;
    
    std::string jetRname = "_R0.4"; // + jetR;
    std::string thrname = "_t1.0"; // + threshold;
    std::string weightstr = ""; //"_xx";
    
    // initializing objects
    TChain *JETINFO_tree = new TChain("JETINFO_tree");
    TChain *PAIRINFO_tree = new TChain("PAIRINFO_tree");
    
    
    
    // make TChains
    for ( int ifile = 1; ifile <= 1; ifile++ ) {
        
        // at this point, make a tchain
        //TODO: FIX THIS (below)
        //        std::string JETINFO_fulltreename = Form("%s/%d/%d/AnalysisResults.root/%s", base_filepath.c_str(), ifile, jfile, JETINFO_name.c_str()); //TODO: this needs to be fixed on perly
        std::string JETINFO_fulltreename = Form("%s/LHC17q_CENT_wSDD/000282365/0085/AnalysisResults.root/%s", base_filepath.c_str(), JETINFO_name.c_str());
        JETINFO_tree->Add(JETINFO_fulltreename.c_str());
        if (debug) cout << "num JETINFO tree entries " << JETINFO_tree->GetEntries() << endl;
        
        std::string PAIRINFO_fulltreename = Form("%s/LHC17q_CENT_wSDD/000282365/0085/AnalysisResults.root/%s", base_filepath.c_str(), PAIRINFO_name.c_str()); //TODO: this needs to be fixed on perly
        PAIRINFO_tree->Add(PAIRINFO_fulltreename.c_str());
        if (debug) cout << "num PAIRINFO tree entries " << PAIRINFO_tree->GetEntries() << endl;
    }
            
            

    // debug
    if (debug2) PAIRINFO_tree->Print();
    
    // now look at observables and make histograms
    // don't separate by pt or RL bin
    TH1D * jetpt_hist = getObs1DHistFromTChain(JETINFO_tree, "jet_pt", 250, 0, 200, 0, 200, 1e-4, 1.0);
    TCanvas *can_jetpt = new TCanvas();
    draw_save_del_hists(can_jetpt, jetpt_hist, "jetpt", weightstr + jetRname + thrname);
    
    // needs to be separated by pt and RL bin
    for ( int i = 0; i < n_bins; i++ ) {
        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];
        std::string ptname = "_pt" + to_string(pt_bins[i]) + "-" + to_string(pt_bins[i+1]);
        if (debug) cout << " in pt bin" << i << " with " << pt_min << " - " << pt_max << endl;
        
        
        for ( int j = 0; j < n_RLbins; j++ ) {
            
            double RL_min = RL_bins[i][j];
            double RL_max = RL_bins[i][j+1];
            std::string RLname = Form("_RL%.3f-%.3f", RL_min, RL_max);
            std::string hist_addname = weightstr + jetRname + thrname + ptname + RLname;
            if (debug) cout << " in RL bin" << j << " with " << RL_min << " - " << RL_max << endl;
            
            // get histograms
            TH1D * deltap_hist = getObs1DHistFromTChain(PAIRINFO_tree, "deltap", 250, 0, pt_max+5, pt_min, pt_max, RL_min, RL_max);
            cout << "checkpoint 1 " << deltap_hist->GetEntries() << endl;
            TH1D * deltapt_hist = getObs1DHistFromTChain(PAIRINFO_tree, "deltapt", 250, 0, pt_max+5, pt_min, pt_max, RL_min, RL_max);
            TH1D * deltapl_hist = getObs1DHistFromTChain(PAIRINFO_tree, "deltapl", 250, 0, pt_max+5, pt_min, pt_max, RL_min, RL_max);
            TH1D * energyweight_hist = getObs1DHistFromTChain(PAIRINFO_tree, "weights", 250, 0, 1.0, pt_min, pt_max, RL_min, RL_max);
            


            // format histograms

            // draw, save, and delete histograms
            TCanvas *can_deltap = new TCanvas();
            TCanvas *can_deltapt = new TCanvas();
            TCanvas *can_deltapl = new TCanvas();

            draw_save_del_hists(can_deltap, deltap_hist, "deltap", hist_addname);
            draw_save_del_hists(can_deltapt, deltapt_hist, "deltapt", hist_addname);
            draw_save_del_hists(can_deltapl, deltapl_hist, "deltapl", hist_addname);
            
            
        }
        
        
        
        
    }
    
    // delete objects after saving for new pt-hat bin
    delete JETINFO_tree;
    delete PAIRINFO_tree;
    
  
}


