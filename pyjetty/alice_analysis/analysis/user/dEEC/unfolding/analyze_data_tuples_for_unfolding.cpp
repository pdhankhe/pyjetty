// ROOT macro to analyze and plot data tuples for unfolding purposes
// Remember that this inputs a bunch of files with NTuples (not one merged file)
// Needed for unfolding:
//   EW: 3D hist (weight, RL, pT), 1D hist jet pt 
//     ^x3 to get each configuration of RL bins
// Needed for bin-by-bin corrections (not in this file):
//   ∆p/∆pt/∆pl/rc: raw distributions 
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

std::string attempt_dir = "data_for_unfolding";
std::string outdir = "/software/users/blianggi/mypyjetty/storage/dEEC/plots/" + attempt_dir;

void SetStyle(Bool_t graypalette=true) {
    cout << "Setting style!" << endl;
  
    gStyle->Reset("Plain");
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    // if(graypalette) gStyle->SetPalette(8,0);
    // else gStyle->SetPalette(1);
    gStyle->SetPalette(kRainbow); //kBird
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

void ProcessCanvas(TCanvas *Canvas, bool moveright=false) { 
	gStyle->SetOptStat(0);
    if (moveright) gStyle->SetPadRightMargin(0.15);

	Canvas->SetHighLightColor(1);
	Canvas->SetFillColor(0);
	Canvas->SetBorderMode(0);
	Canvas->SetBorderSize(2);
	Canvas->SetTickx(1);
	Canvas->SetTicky(1);
	Canvas->SetFrameBorderMode(0);
	Canvas->SetFrameLineWidth(1);
 	Canvas->SetFrameBorderMode(1);
}

/* get a typical 1D histogram from the TChain */
TH1D * getObs1DHistFromTChain(TChain *chain, std::string branch_name, int num_bins, double hist_xmin, double hist_xmax,
                              int pt_min, int pt_max, double RL_min, double RL_max) {
    
    TH1D * hist1D;
    if ( branch_name == "jet_pt" ) {
        double binEdges[] = {5, 10, 20, 40, 60, 80, 100, 150};
        hist1D = new TH1D(Form("%s_hist", branch_name.c_str()), Form("%s_hist", branch_name.c_str()), num_bins, binEdges);
    
        chain->Draw(Form("%s>>%s_hist", branch_name.c_str(), branch_name.c_str()), Form("jet_pt >= %d && jet_pt < %d", pt_min, pt_max), "e");

        // TH1D *rebinnedHist = hist1D->Rebin(7, Form("%s_hist", branch_name.c_str()), binEdges);
        // delete hist1D;
        // TH1D * hist1D = rebinnedHist->Clone(rebinnedHist->GetName());

    } else {
        hist1D = new TH1D(Form("%s_hist", branch_name.c_str()), Form("%s_hist", branch_name.c_str()), num_bins, hist_xmin, hist_xmax);
        if ( branch_name == "total_num_const" || branch_name == "num_const_aftercut") {
        chain->Draw(Form("%s>>%s_hist", branch_name.c_str(), branch_name.c_str()), Form("jet_pt >= %d && jet_pt < %d", pt_min, pt_max), "e");
        } else {
            chain->Draw(Form("%s>>%s_hist", branch_name.c_str(), branch_name.c_str()), Form("jet_pt >= %d && jet_pt < %d && RL >= %f && RL < %f", pt_min, pt_max, RL_min, RL_max), "e");
        }
    }
//    chain->Draw(Form("%s>>%s_hist", branch_name.c_str(), branch_name.c_str()), Form("jet_pt >= 20 && jet_pt < 40 && RL > 0.01 && RL < 0.4"), "e");

    return hist1D;
}




/* get a typical 2D histogram from the TChain */
TH2D * getObs2DHistFromTChain(TChain *chain, std::string branch_name_x, std::string branch_name_y,
                              int num_bins_x, double hist_xmin, double hist_xmax,
                              int num_bins_y, double hist_ymin, double hist_ymax,
                              int pt_min, int pt_max, double RL_min, double RL_max) {
    
    TH2D * hist2D = new TH2D(Form("%s_vs_%s_hist", branch_name_y.c_str(), branch_name_x.c_str()), Form("%s_vs_%s_hist", branch_name_y.c_str(), branch_name_x.c_str()), num_bins_x, hist_xmin, hist_xmax, num_bins_y, hist_ymin, hist_ymax);
    chain->Draw(Form("%s:%s>>%s_vs_%s_hist", branch_name_y.c_str(), branch_name_x.c_str(), branch_name_y.c_str(), branch_name_x.c_str()), Form("jet_pt >= %d && jet_pt < %d && RL >= %f && RL < %f", pt_min, pt_max, RL_min, RL_max)); //, "colz");
    
    return hist2D;
}

/* get a typical 3D histogram from the TChain */
TH3D * getObs3DHistFromTChain(TChain *chain, std::string branch_name_x, std::string branch_name_y,
                              std::string branch_name_z,
                              int num_bins_x, double hist_xmin, double hist_xmax,
                              const double RL_bins[], int n_RLbins) {
    
    int n_ptbins = 7;
    double pt_binedges[] = {5.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 150.0};

    double EW_binedges[] = {0.   , 0.006, 0.012, 0.018, 0.024, 0.03 , 0.036, 0.042, 0.048,
                            0.054, 0.06 , 0.066, 0.072, 0.078, 0.084, 0.09 , 0.096, 0.102,
                            0.108, 0.114, 0.12 , 0.126, 0.132, 0.138, 0.144, 0.15 , 0.156,
                            0.162, 0.168, 0.174, 0.18 , 0.186, 0.192, 0.198, 0.204, 0.21 ,
                            0.216, 0.222, 0.228, 0.234, 0.24 , 0.246, 0.252, 0.258, 0.264,
                            0.27 , 0.276, 0.282, 0.288, 0.294, 0.3  };
    
    // TH3D * hist3D = new TH3D(Form("%s_vs_%s_hist", branch_name_y.c_str(), branch_name_x.c_str()), Form("%s_vs_%s_hist", branch_name_y.c_str(), branch_name_x.c_str()), num_bins_x, hist_xmin, hist_xmax, num_bins_y, hist_ymin, hist_ymax);
    std::string name = Form("%s_vs_%s_vs_%s_hist", branch_name_z.c_str(), branch_name_y.c_str(), branch_name_x.c_str());
    TH3D * hist3D = new TH3D(name.c_str(), name.c_str(), num_bins_x, EW_binedges, n_RLbins, RL_bins, n_ptbins, pt_binedges);
    // TH3D* hist3D = new TH3D(name.c_str(), name.c_str(),
    //                     num_bins_x, hist_xmin, hist_xmax,  // Uniform bins for x-axis
    //                     n_RLbins, RL_bins,               // Custom bins for y-axis
    //                     n_ptbins, pt_binedges); 
    chain->Draw(Form("%s:%s:%s>>%s", branch_name_z.c_str(), branch_name_y.c_str(), branch_name_x.c_str(), name.c_str())); //, "colz");
    
    return hist3D;
}


/* Format and adjust histograms */
void Format1DHist(TH1D *hist, TH1D *jetpt_hist, std::string norm_string, int markercolor, double markeralpha,
                  int markerstyle, std::string xtitle, std::string ytitle, TLegend& leg, TString leg_text, 
                  std::string obs_name="", std::string hist_addname="") {

    hist->SetTitle(Form("h_%s%s", obs_name.c_str(), hist_addname.c_str()));
    hist->SetName(Form("h_%s%s", obs_name.c_str(), hist_addname.c_str()));

    // normalization
    if ( norm_string == "self_normalized" ) {
        double selfnorm_value = hist->Integral();
        hist->Scale(1/selfnorm_value, "width");
    } else if ( norm_string == "norm_by_jets" ) {
        double numjets = jetpt_hist->Integral();
        cout << "Number of jets in " << leg_text << ": " << numjets << endl;
        hist->Scale(1/numjets, "width");
    }

    // stylization
    hist->SetLineColorAlpha(markercolor, markeralpha);
    hist->SetMarkerColorAlpha(markercolor, markeralpha);
    hist->SetMarkerStyle(markerstyle);
    hist->SetMarkerSize(1.5);

    // axes
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetTitleFont(42);
	if (obs_name == "weights") {
        hist->GetXaxis()->SetTitleSize(0.035);
        hist->GetXaxis()->SetTitleOffset(1.5);
    } else {
        hist->GetXaxis()->SetTitleSize(0.06); //(0.042);
        hist->GetXaxis()->SetTitleOffset(1.0);
    }
	hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetTitle(xtitle.c_str());

    hist->GetYaxis()->SetLabelFont(42);
	hist->GetYaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleOffset(1.05); 
	hist->GetYaxis()->SetTitleSize(0.06); //(0.042);
	hist->GetYaxis()->SetLabelSize(0.05); //(0.042);
    hist->GetYaxis()->SetTitle(ytitle.c_str());

    // legend
    leg.AddEntry(hist, leg_text, "pl");

}

void Format2DHist(TH2D *hist2D, TH1D *jetpt_hist, std::string norm_string, std::string xtitle, std::string ytitle, 
                  bool scalebyRLbinwidth, double RL_bin_width, std::string obs_name_x, std::string obs_name_y, std::string hist_addname="") {

    double ptvsew_normbounds[3][2] = { { 0, 500 }, { 1e-2, 5e2 }, { 1e-4, 2e2 }}; //TODO: make this less pt specific?
    // also TODO: potentially set all lower boundaries to 0 for data??

    hist2D->SetTitle(Form("h_%s_vs_%s%s", obs_name_y.c_str(), obs_name_x.c_str(), hist_addname.c_str()));
    hist2D->SetName(Form("h_%s_vs_%s%s", obs_name_y.c_str(), obs_name_x.c_str(), hist_addname.c_str()));

    // normalization
    int norm_index = 0;
    if ( scalebyRLbinwidth ) hist2D->Scale(RL_bin_width);
    if ( norm_string == "self_normalized" ) {
        double selfnorm_value = hist2D->Integral();
        hist2D->Scale(1/selfnorm_value, "width");
        norm_index = 1;
    } else if ( norm_string == "norm_by_jets" ) {
        double numjets = jetpt_hist->Integral();
        hist2D->Scale(1/numjets, "width");
        norm_index = 2;
    }

    // set z axis bounds
    hist2D->GetZaxis()->SetRangeUser(ptvsew_normbounds[norm_index][0], ptvsew_normbounds[norm_index][1]);

    // label axes
    hist2D->GetXaxis()->SetLabelFont(42);
    hist2D->GetXaxis()->SetTitleFont(42);
    hist2D->GetXaxis()->SetTitleSize(0.04); //(0.042);
    hist2D->GetXaxis()->SetTitleOffset(1.3);
	hist2D->GetXaxis()->SetLabelSize(0.05);
    hist2D->GetXaxis()->SetTitle(xtitle.c_str());

    hist2D->GetYaxis()->SetLabelFont(42);
	hist2D->GetYaxis()->SetTitleFont(42);
    // if (obs_name_y == "weights") {
    //     hist2D->GetYaxis()->SetTitleSize(0.035);
    //     hist2D->GetYaxis()->SetTitleOffset(1.5);
    // } else {
        hist2D->GetYaxis()->SetTitleSize(0.04); //0.06 //(0.042);
        hist2D->GetYaxis()->SetTitleOffset(1.3);
    // }
    hist2D->GetYaxis()->SetLabelSize(0.04); //(0.042);
    hist2D->GetYaxis()->SetTitle(ytitle.c_str());
}

void Format3DHist(TH3D *hist3D, int ptbin,
                  std::string xtitle, std::string ytitle, std::string ztitle,
                  std::string obs_name_x, std::string obs_name_y, std::string obs_name_z) {


    // name = Form("h_%s_vs_%s_vs_%s_PTBIN%d", "weights", "RL", "jet_pt", ptbin);
    std::string name = Form("h_%s_vs_%s_vs_%s_PTBIN%d", obs_name_z.c_str(), obs_name_y.c_str(), obs_name_x.c_str(), ptbin);
    hist3D->SetNameTitle(name.c_str(), name.c_str());



    // label axes
    hist3D->GetXaxis()->SetLabelFont(42);
    hist3D->GetXaxis()->SetTitleFont(42);
    hist3D->GetXaxis()->SetTitleSize(0.04); //(0.042);
    hist3D->GetXaxis()->SetTitleOffset(1.3);
	hist3D->GetXaxis()->SetLabelSize(0.05);
    hist3D->GetXaxis()->SetTitle(xtitle.c_str());

    hist3D->GetYaxis()->SetLabelFont(42);
	hist3D->GetYaxis()->SetTitleFont(42);
    // if (obs_name_y == "weights") {
    //     hist3D->GetYaxis()->SetTitleSize(0.035);
    //     hist3D->GetYaxis()->SetTitleOffset(1.5);
    // } else {
        hist3D->GetYaxis()->SetTitleSize(0.04); //0.06 //(0.042);
        hist3D->GetYaxis()->SetTitleOffset(1.3);
    // }
    hist3D->GetYaxis()->SetLabelSize(0.04); //(0.042);
    hist3D->GetYaxis()->SetTitle(ytitle.c_str());

    hist3D->GetZaxis()->SetLabelFont(42);
    hist3D->GetZaxis()->SetTitleFont(42);
    hist3D->GetZaxis()->SetTitleSize(0.04); //(0.042);
    hist3D->GetZaxis()->SetTitleOffset(1.3);
	hist3D->GetZaxis()->SetLabelSize(0.05);
    hist3D->GetZaxis()->SetTitle(ztitle.c_str());
}


/* Save and delete histograms */
// (TFile *fout, TCanvas *can, TH1D* hist, std::string obs_name, 
//                          std::string ptname, std::string norm_string, std::string hist_addname,
//                          bool logx, bool logy)
void draw_save_del_hists(TFile *fout, TCanvas *can, TObject* obj, std::string obs_name, 
                         std::string ptname, std::string norm_string, std::string hist_addname,
                         bool savepdf, bool logx, bool logy, bool logz=false) {
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
    // hist->Draw();

    fout->cd();
    // hist->Write();
    obj->Write(); //TODO: this might not be right! Might have to use the casted type

    // size_t length = hist_vec.size();
    // if ( length > 0 ) hist_vec.push_back(hist);

    

    // std::string outdir = "/software/users/blianggi/mypyjetty/storage/dEEC/plots/data_firstattempt"; // + ptbin_name + "/";//"plots/test/";
    std::string add_dir = "";
    if (obs_name != "jet_pt" && obs_name != "total_num_const" && obs_name != "num_const_aftercut") {
        if (obs_name == "rc") add_dir = "/" + ptname + "/" + norm_string + "/" + obs_name;
        else add_dir = "/" + ptname + "/" + norm_string + "/" + obs_name + "/individuals";
    }
    std::string fname_out = outdir + add_dir + "/corrhist_" + obs_name + hist_addname + ".pdf";
    if (savepdf) can->SaveAs(fname_out.c_str());

    // delete hist;
    delete can;
}


/* Delete a vector of histograms */
void deleteVecOfHists(std::vector<TH1D*>& histVector) {
    // Loop through the vector and delete each TH1D pointer
    for (TH1D* hist : histVector) {
        delete hist; // Free the memory allocated for the histogram
    }

    // Clear the vector to remove all the pointers
    histVector.clear();
}


/* plot all RL bins in one plot */
void plotandsave_combined_hists(TCanvas *can_all, vector<TH1D*> h_vec, TLegend *l, 
                          std::string obs_name, std::string ptname, 
                          std::string norm_string, std::string hist_addname,
                          int pt_max, bool scalebyRLbinwidth, double RL_bin_width[],
                          bool logx, bool logy, double pl_axis_cut=-1, bool debug=false) {

    // go into canvas
    can_all->cd();
    if (logx) gPad->SetLogx();
    if (logy) gPad->SetLogy();

    // if momentum axis, adjust x bounds accordingly
    size_t length = h_vec.size();
    for (int j=0; j<length; j++) {
        // cout << j << ": " << RL_bin_width[j] << endl;
        if (scalebyRLbinwidth) h_vec[j]->Scale(RL_bin_width[j]); // this needs to be done before normalization
        // if (mom_axis) {
        //     // h_vec[j]->Rebin(4);
        //     // h_vec[j]->GetXaxis()->SetRangeUser(0, pt_max+5);
        //     // cout << h_vec[j]->GetEntries() << endl;
        //     h_vec[j]->Scale(RL_bin_width[j]); // this needs to be done before normalization
        // }
        if (pl_axis_cut > 0) {
            h_vec[j]->GetXaxis()->SetRangeUser(0, pl_axis_cut);
        }
    }

    // set maximum based on maximum of all curves
    double max = 0;
    for (int j=0; j<length; j++) {
        double max_cand = h_vec[j]->GetMaximum();
        if (max_cand > max) max = max_cand;
    }
    if (debug) cout << "max is " << max << " which goes to " << max*1.5 << endl;
    h_vec[0]->SetMaximum( max * 1.5 );

    // draw!
    for (int j=0; j<length; j++) {
        h_vec[j]->Draw("same");
    }

    // // axes
    // hist->GetXaxis()->SetLabelFont(42);
    // hist->GetXaxis()->SetTitleFont(42);
	// hist->GetXaxis()->SetTitleOffset(1.0);
	// hist->GetXaxis()->SetTitleSize(0.06); //(0.042);
	// hist->GetXaxis()->SetLabelSize(0.05);
    // hist->GetXaxis()->SetTitle(xtitle.c_str());

    // hist->GetYaxis()->SetLabelFont(42);
	// hist->GetYaxis()->SetTitleFont(42);
    // hist->GetYaxis()->SetTitleOffset(1.05); 
	// hist->GetYaxis()->SetTitleSize(0.06); //(0.042);
	// hist->GetYaxis()->SetLabelSize(0.05); //(0.042);
    // hist->GetYaxis()->SetTitle(ytitle.c_str());

    // can_all->Modified();
    // can_all->Update();
    l->Draw("same");

    //save as PDF
    // std::string outdir = "/software/users/blianggi/mypyjetty/storage/dEEC/plots/data_firstattempt/"; // + ptbin_name + "/";//"plots/test/";
    std::string add_dir = "/" + ptname + "/" + norm_string + "/" + obs_name;
    std::string fname_out = outdir + add_dir + "/corrhist_" + obs_name + "_ALL" + hist_addname + ".pdf";
    can_all->SaveAs(fname_out.c_str());

    deleteVecOfHists(h_vec);
    delete can_all;

}

// ======================================================= //
//                   SPECIFIC FUNCTIONS
// ======================================================= //



// ======================================================= //
//                     SOME FUNCTIONS
// ======================================================= //

void analyze_ptbin(TChain * JETINFO_tree, TChain * PAIRINFO_tree, 
             TFile *f_out, int ptbin, const double RL_bins[], int n_RLbins, 
             bool debug, bool debug2) {
    
    
    TH3D * hist3D_EW = getObs3DHistFromTChain(PAIRINFO_tree, "weights", "RL", "jet_pt", 50, 0, 0.3, RL_bins, n_RLbins);
    
    // hist3D_EW->SetTitle(Form("h_%s_vs_%s_vs_%s_PTBIN%d", "weights", "RL", "jet_pt", ptbin));
    // hist3D_EW->SetName(Form("h_%s_vs_%s_vs_%s_PTBIND%d", "weights", "RL", "jet_pt", ptbin));
    Format3DHist(hist3D_EW, ptbin, "EW", "R_{L}", "p_{T, jet}", "weights", "RL", "jet_pt");

    TCanvas *can_weights = new TCanvas();

    draw_save_del_hists(f_out, can_weights, hist3D_EW, "", "", "", "", false, false, false, false);
    
 
}


//TODO: do something about norm_string!!
void analyze(TChain * JETINFO_tree, TChain * PAIRINFO_tree, 
             TFile *f_out, const int pt_bins[], int n_bins, const double RL_bins[][8], int n_RLbins,
             bool include_RL0, bool include_RL1, bool debug, bool debug2 ) {

    
    // now look at observables and make histograms
    // don't separate by pt or RL bin
    cout << "running jet pt hist now" << endl;
    TH1D * jetpt_hist = getObs1DHistFromTChain(JETINFO_tree, "jet_pt", 7, 0, 200, 0, 200, 0, 0);
    jetpt_hist->GetXaxis()->SetTitle("p_{T,jet}");
    TCanvas *can_jetpt = new TCanvas();
    draw_save_del_hists(f_out, can_jetpt, jetpt_hist, "", "", "", "", false, false, true);
    
    
    // return;



    // needs to be separated by pt and RL bin
    for ( int i = 0; i < n_bins; i++ ) {
        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];
           
        if (debug) cout << " in pt bin" << i << " with " << pt_min << " - " << pt_max << endl;
        analyze_ptbin(JETINFO_tree, PAIRINFO_tree, f_out, i, RL_bins[i], n_RLbins, debug, debug2);
    }



    


}



// ======================================================= //
//                     MAIN FUNCTION
// ======================================================= //

void analyze_data_tuples_for_unfolding() {
    gStyle->SetOptStat(0);
    SetStyle();
    
    // setup variables
    bool debug = true;
    bool debug2 = false;
    
    // ntuple/histogram names
    std::string JETINFO_name = "tn_JETINFO_R0.4_1.0";
    std::string PAIRINFO_name = "tn_pairlevel_R0.4_1.0";
    std::string jet1D_name = "h_1Djet_pt_JetPt_R0.4_1.0"; // this one is a histogram
        
    // filenames
    std::string base_filepath_hic = Form("/rstorage/alice/AnalysisResults/blianggi/dEEC/451107");
    
    // Output file for binned results
    std::string root_outfile = "/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/" + attempt_dir + "/DataForUnfolding.root"; //plots/ntuples/DataHists.root"; //FinalDataHists.root
    TFile* f_out = new TFile(root_outfile.c_str(), "RECREATE");
    std::string add_name = ""; //"_othercorrel";


    // analysis variables
    const int pt_bins[] = { 20, 40, 60, 80 };
    const int n_bins = sizeof(pt_bins) / sizeof(pt_bins[0]) - 1; //3;
    
    const double RL_bins[3][8] = { { 0, 1e-2, 3e-2, 7e-2, 1.5e-1, 3e-1, 4e-1, 1 },
                            { 0, 1e-2, 2.5e-2, 4e-2, 8e-2, 2.5e-1, 4e-1, 1 },
                            { 0, 1e-2, 2.5e-2, 3e-2, 4.5e-2, 2e-1, 4e-1, 1 } };
    const int n_RLbins = sizeof(RL_bins[0]) / sizeof(RL_bins[0][0]) - 1; //gets the columns //6; //7; //5;


    if (debug2) cout << "pt_bins " << n_bins << " n_RLbins " << n_RLbins << endl;
    

    std::string jetRname = "_R0.4"; // + jetR;
    std::string thrname = "_t1.0"; // + threshold;
    std::string weightstr = ""; //"_xx";
    std::string norm_string = "";
    bool include_RL0 = false;
    bool include_RL1 = false;
    
    // initializing objects
    TChain *JETINFO_tree = new TChain("JETINFO_tree");
    TChain *PAIRINFO_tree = new TChain("PAIRINFO_tree");

    // ====================================================================================
    /*------------------------------------------------------------
    //----------------------- NTUPLE INFO ------------------------
    JETINFO: jet_pt; total_num_const; num_const_aftercut; total_num_baryons; num_baryons_aftercut; total_num_mesons; num_mesons_aftercut
    PAIRINFO: jet_pt; RL; weights; deltap; p1; p2; deltapt; pt1; pt2; deltapl; pl1; pl2; q1q2; q1; q2; baryonmeson; pid1; pid2
    baryon: jet_pt; baryon_pt
    meson: jet_pt; meson_pt
    //----------------------------------------------------------*/
    
    
    // make TChains
    std::ifstream filelist("/software/users/blianggi/mypyjetty/dEEC/filelist_datatuples_451107_shortname.txt");
    if (!filelist.is_open()) {
        std::cerr << "Error: Could not open /software/users/blianggi/mypyjetty/dEEC/filelist_datatuples_451107_shortname.txt" << std::endl;
        return;
    }

    std::string ntuple_filename;
    int filecounter = 0;
    int filecounter_cutoff = -1; //total: 7601
    // Loop through each line in filelist
    while (std::getline(filelist, ntuple_filename)) {

        if (filecounter == filecounter_cutoff) break;

        std::string JETINFO_fulltreename = Form("%s/%s/%s", base_filepath_hic.c_str(), ntuple_filename.c_str(), JETINFO_name.c_str());
        JETINFO_tree->Add(JETINFO_fulltreename.c_str());
        
        std::string PAIRINFO_fulltreename = Form("%s/%s/%s", base_filepath_hic.c_str(), ntuple_filename.c_str(), PAIRINFO_name.c_str()); //TODO: this needs to be fixed on perly
        PAIRINFO_tree->Add(PAIRINFO_fulltreename.c_str());

        if (debug) {
            if (filecounter%100 == 0) {
                cout << "num JETINFO tree entries " << JETINFO_tree->GetEntries() << endl;
                cout << "num PAIRINFO tree entries " << PAIRINFO_tree->GetEntries() << endl;
            }
        }
        

        filecounter++;
    }

    // Close the filelist.txt file
    filelist.close();
            
    // ====================================================================================

    // debug
    if (debug2) PAIRINFO_tree->Print();
    
    // analyze to get histograms needed for unfolding
    analyze(JETINFO_tree, PAIRINFO_tree, f_out, pt_bins, n_bins, RL_bins, n_RLbins, include_RL0, include_RL1, debug, debug2);


    // delete objects after saving for new pt-hat bin
    delete JETINFO_tree;
    delete PAIRINFO_tree;

    f_out->Close();
    
  
}


