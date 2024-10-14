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

TH1D * getObs1DHistFromTChain(TChain *chain, std::string branch_name, int num_bins, int pt_min, int pt_max, double RL_min, double RL_max) {
    
    TH1D * hist1D = new TH1D(Form("%s_hist", branch_name.c_str()), Form("%s_hist", branch_name.c_str()), num_bins, 0, pt_max+5);
    chain->Draw(Form("%s>>%s_hist", branch_name.c_str(), branch_name.c_str()), Form("jet_pt >= %d && jet_pt < %d && RL > %f && RL < %f", pt_min, pt_max, RL_min, RL_max), "e");
    
    return hist1D;
}




// ======================================================= //
//                     MAIN FUNCTION
// ======================================================= //

void make_histograms_from_tuples() {
    gStyle->SetOptStat(0);
    SetStyle();
    
    // setup variables
    int N_pthatbins = 1; //20
    int N_files_per_pthatbin = 1; //100
    bool debug = true;
    bool debug2 = false;
    
    // ntuple names
    std::string JETINFO_name = "h_JETINFOjet_pt_JetPt_Truth_R0.4_1.0";
    std::string PAIRINFO_name = "tn_pairlevel_Truth_R0.4_1.0";
    std::string baryon_name = "tn_baryon_JetPt_Truth_R0.4_1.0";
    std::string meson_name = "tn_meson_JetPt_Truth_R0.4_1.0";
    
    // other names
    std::string jet1D_name = "h_1Djet_pt_JetPt_Truth_R0.4_1.0"; // this one is a histogram
    std::string deltap_str = "deltap";
    
    // filenames
    std::string filename = Form("~/Documents/research/othercorrelations/AnalysisResults_ntuplesnew.root");
    std::string base_filepath = Form("/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/pythiagen/histograms/30370699/ntuples/26652369");
        
    
    // Output file for binned results
    std::string root_outfile = "plots/ntuples/FinalHistsAfterScaling.root";
    TFile* f_out = new TFile(root_outfile.c_str(), "RECREATE");
    std::string add_name = "_othercorrel";

    
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
    std::string weightstr = "";
    
    
    
    
    // make TChains
    for ( int ifile = 10; ifile <= 10; ifile++ ) { //1; ifile <= N_pthatbins; ifile++ ) {
        
        // at this point, make a tchain
        // initializing objects
        TChain *JETINFO_tree = new TChain("JETINFO_tree");
        TChain *PAIRINFO_tree = new TChain("PAIRINFO_tree");
        TChain *baryon_tree = new TChain("baryon_tree");
        TChain *meson_tree = new TChain("meson_tree");
        
        
        for ( int jfile = 1; jfile <= N_files_per_pthatbin; jfile++ ) {
            
            std::string JETINFO_fulltreename = Form("%s/%d/%d/AnalysisResults.root/%s", base_filepath.c_str(), ifile, jfile, JETINFO_name.c_str()); //TODO: this needs to be fixed on perly
            JETINFO_tree->Add(JETINFO_fulltreename.c_str());
            if (debug) cout << "num JETINFO tree entries " << JETINFO_tree->GetEntries() << endl;

            std::string PAIRINFO_fulltreename = Form("%s/%d/%d/AnalysisResults.root/%s", base_filepath.c_str(), ifile, jfile, PAIRINFO_name.c_str()); //TODO: this needs to be fixed on perly
            PAIRINFO_tree->Add(PAIRINFO_fulltreename.c_str());
            if (debug) cout << "num PAIRINFO tree entries " << PAIRINFO_tree->GetEntries() << endl;
            
            
        }

        // debug
        if (debug2) PAIRINFO_tree->Print();
        
        // now look at observables and make histograms
        // don't separate by pt or RL bin
        TH1D * jetpt_hist = getObs1DHistFromTChain(PAIRINFO_tree, "jet_pt", 250, 0, 200, 1e-4, 1.0);
        if (debug) cout << "num entries in jetpt " << jetpt_hist->GetEntries() << endl;

        // needs to be separated by pt and RL bin
        for ( int i = 0; i < n_bins; i++ ) {
            int pt_min = pt_bins[i];
            int pt_max = pt_bins[i+1];
            std::string ptname = to_string(pt_bins[i]) + "-" + to_string(pt_bins[i+1]);
            
            for ( int j = 0; j < n_RLbins; j++ ) {
                
                double RL_min = RL_bins[i][j];
                double RL_max = RL_bins[i][j+1];
                std::string RLname = Form("_RL%.3f-%.3f", RL_min, RL_max);
                std::string hist_addname = jetRname + thrname + ptname + RLname;
                if (debug) cout << " in RL bin" << j << " with " << RL_min << " - " << RL_max << endl;
                
                double binwidth_for_p = 0.5;
                int numbins_for_p = floor((pt_max+5)/binwidth_for_p); 
                // TH1D * deltap_hist = new TH1D("deltap_hist", "deltap_hist", numbins_for_p, 0, pt_max+5);
                // cout << "PT MIN IS " << pt_min << endl;
                // PAIRINFO_tree->Draw(Form("%s>>%s_hist", "deltap", "deltap"), Form("jet_pt >= %d && jet_pt < %d && RL > 1e-4 && RL < 1", pt_min, 10), "e");
                
                TH1D * deltap_hist = getObs1DHistFromTChain(PAIRINFO_tree, "deltap", numbins_for_p, pt_min, pt_max, RL_min, RL_max);
                TH1D * deltapt_hist = getObs1DHistFromTChain(PAIRINFO_tree, "deltapt", numbins_for_p, pt_min, pt_max, RL_min, RL_max);
                TH1D * deltapl_hist = getObs1DHistFromTChain(PAIRINFO_tree, "deltapl", numbins_for_p, pt_min, pt_max, RL_min, RL_max);
                
                if (debug) cout << "num entries in deltap " << deltap_hist->GetEntries() << endl;

                TCanvas *can = new TCanvas();
                can->cd();
                deltap_hist->Draw();
                
                std::string outdir = "plots/fourthattempt/"; // + ptbin_name + "/";//"plots/test/";
                std::string fname_out = outdir + "testing/corrhist_deltap" + weightstr + "_pt" + ptname + RLname + "_R0.4" + add_name + ".pdf";
                can->SaveAs(fname_out.c_str());
                
                
                delete deltap_hist;
                delete deltapt_hist;
                delete deltapl_hist;
                delete can;
            }
            
            
            
            
        }
        
        
        
        
        // delete objects after saving for new pt-hat bin
        delete JETINFO_tree;
        delete PAIRINFO_tree;
        delete baryon_tree;
        delete meson_tree;
        
        
    }
    
  
}


