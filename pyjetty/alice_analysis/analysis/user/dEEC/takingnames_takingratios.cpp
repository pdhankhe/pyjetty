// This file takes some ratios of the data curves.
// This is to be run on hiccup.
// Right now, this file is only compatible with 4 RL bins.
// outputted plots:
    // ratios_of_dif_regions_[obs_PTBINX].pdf
    // ratios_of_dif_jetptbins_[obs_RLBINY].pdf
// Beatrice Liang-Gilman, beatrice_lg@berkeley.edu

#include "library/fitfunctions.h"

Double_t colors[16] = {kGray, kMagenta, kBlue, kOrange+1, kViolet+1, kGreen+2, kRed, kYellow+1, kCyan+1};
Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
Double_t marker_size = 1.5;

std::string attempt_dir = "data_fourthattempt_ptrlbins/rebinx4";
std::string outputbase = "/software/users/blianggi/mypyjetty/storage/dEEC/plots/" + attempt_dir + "/";

// double rebin = 4;
bool ptrl_bins = true;
bool jetpt_bool = true;

bool deltap_bool = true;
bool p_bool = false; //true;
bool deltajt_bool = true;
bool jt_bool = false; //true;
bool ew_bool = false; //true;
bool twoDhists_bool = false;
bool rc_bool = false;

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



//===========================================================================
//============================== CALCULATIONS ===============================
//===========================================================================

/* Get the r_c from TChain */
double getRc(TH1D *h_q1q2) 
{
    // do i need to scale by the RL bin width here?? - I think this would be redundant.
    // if both like sign bin and unlike sign get scaled by RL bin width, then the ratio still stays the same

    // get # of like sign and # of unlike sign
    double num_likesign = h_q1q2->GetBinContent(h_q1q2->FindBin(1));
    double num_unlikesign = h_q1q2->GetBinContent(h_q1q2->FindBin(-1));
    // if (debug) 
    cout << "num like sign " << num_likesign << " num unlike sign " << num_unlikesign << endl;

    // calculate the rc value for this pt & RL bin
    double rc = (double)(num_likesign - num_unlikesign) / (double)(num_likesign + num_unlikesign);
    //if (debug) 
    cout << "and that makes rc " << rc << endl;

    // if (num_likesign == num_unlikesign) rc = 0;

    return rc;
}

//===========================================================================
//================================ PLOTTING =================================
//===========================================================================


void formathist(TH1D * hist, std::string xtitle, std::string ytitle, int markercolor) {
    // cout << "xtitle: " << xtitle << endl;
    // cout << "ytitle: " << ytitle << endl;
    hist->GetXaxis()->SetTitle(xtitle.c_str());
    hist->GetYaxis()->SetTitle(ytitle.c_str());

    hist->SetMarkerColorAlpha(markercolor, 0.8);
    hist->SetLineColorAlpha(markercolor, 0.8);

    hist->GetXaxis()->SetLabelFont(42);
	hist->GetXaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetTitleOffset(1.05); 
	hist->GetXaxis()->SetTitleSize(0.042); //(0.042);
	hist->GetXaxis()->SetLabelSize(0.042); //(0.042);

    hist->GetYaxis()->SetLabelFont(42);
	hist->GetYaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleOffset(1.45); 
	hist->GetYaxis()->SetTitleSize(0.042); //(0.042);
	hist->GetYaxis()->SetLabelSize(0.042); //(0.042);
}

void plot_ratios(std::vector<TH1D*> ratio_vec, std::vector<std::string> label_vec,
                 std::vector<Double_t> color_vec, std::string filename, TFile * fout,
                 std::string xtitle, std::string ytitle) {

    TCanvas *can = new TCanvas();
    TLegend *l = new TLegend(0.65, 0.62, 0.85, 0.85);
    l->SetTextSize(0.037);

    can->cd();
    // if there are any log axes, set that here?
    // gPad->SetLogy();

    // do we want to get a max/min? if so, set that here
    double max_val = 0;
    for (int a=0; a<ratio_vec.size(); a++) {
        if (ratio_vec[a]->GetMaximum() > max_val) max_val = ratio_vec[a]->GetMaximum();
        ratio_vec[a]->GetXaxis()->SetTitle(xtitle.c_str());
        ratio_vec[a]->GetYaxis()->SetTitle(ytitle.c_str());
    }
    ratio_vec[0]->SetMaximum(max_val*1.1);
    
    for (int a=0; a<ratio_vec.size(); a++) {
        formathist(ratio_vec[a], xtitle, ytitle, color_vec[a]); // fill out this function
        l->AddEntry(ratio_vec[a], label_vec[a].c_str(), "pl");
        ratio_vec[a]->Draw("same");

        fout->cd();
        ratio_vec[a]->Write();
    }

    l->Draw("same");

    // save the plot
    can->SaveAs(filename.c_str());

    delete can;


}



//===========================================================================
//============================= OTHER FUNCTIONS =============================
//===========================================================================

//===========================================================================
// This function gets the bin centers of a histogram and returns it as a 
// vector of doubles.
//===========================================================================
std::vector<double> get_bin_centers(TH1D * hist) {

    std::vector<double> bincenters;

    // Get the number of bins
    int nBins = hist->GetNbinsX();

    // Loop through bins and get their centers
    for (int bin = 1; bin <= nBins; bin++) { // Bins start at 1 in ROOT
        double center = hist->GetBinCenter(bin);
        bincenters.push_back(center);
        // std::cout << "Bin " << bin << " center: " << center << std::endl;
    }

    return bincenters;
}

//===========================================================================
//============================= OTHER FUNCTIONS =============================
//===========================================================================

//===========================================================================
// This function gets the mean of the distributions and plots them against RL.
//===========================================================================

TH1D * getHist(TFile * file, std::string observable, int pt_min, int pt_max, 
             double RL_min, double RL_max, std::string norm_string) {
    
    // std::string histname = Form("h_%s_R0.4_t1.0_pt%d-%d_pTRL%.3f-%.3f_%s", observable.c_str(), pt_min, pt_max, RL_min, RL_max, norm_string.c_str());
    std::string histname = Form("h_%s_R0.4_t1.0_pt%d-%d_pTRL%.3f-%.3f_%s_corrected", observable.c_str(), pt_min, pt_max, RL_min, RL_max, norm_string.c_str());
    cout << "HIST: " << observable << " // " << pt_min << "-" << pt_max << " // " << RL_min << "-" << RL_max << " // " << histname << endl;
    TH1D * hist = (TH1D * )file->Get(histname.c_str());

    return hist;
}

//===========================================================================
// This function takes the ratios of the functions.
//===========================================================================
TH1D * takeratio(TH1D * hist1, TH1D * hist2, 
                bool debug = false) {
    
    TH1D * hratio = (TH1D*) hist1->Clone();
    hratio->Divide(hist2);
    return hratio;
}

// void analyze_ptbin(TFile *fin, std::string observable, int pt_min, int pt_max, 
//                    double RL_min, double RL_max, std::string norm_string) {
    
//     for (int j=0; j<n_RLbins; j++) {

//     } 

// }

void analyze(TFile *fin, TFile *fout, std::string observable, int n_bins, const int pt_bins[], 
             int n_RLbins, const double RL_bins[][7], std::string norm_string, std::string xtitle,
             bool include_RL0, bool include_RL1) {
    //histname = h_deltap_R0.4_t1.0_pt20-40_pTRL0.200-0.800_self_normalized


    // get all the histograms
    // within each pt bin, take the ratios hadronic/pert, hadronization/pert, other pert/pert
    // plot those ratios on one plot (0-1??), one plot for each pt bin
    //then within each RL bin, also take ratios of 40-60/20-40, and 60-80/20-40
    //plot those as well, one plot for each rl bin

    // get all histograms
    std::vector<std::vector<TH1D *>> obs_vec;
    for (int i=0; i<n_bins; i++) {
        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];

        std::vector<TH1D*> temp_vec;
        for (int j=0; j<n_RLbins; j++) {
            int k = j;
            if (!include_RL0) {
                k = j-1;
                if (j == 0) continue; // can add something here to change the filename for ALL
            }
            if (!include_RL1 && j == n_RLbins-1) continue;

            double RL_min = RL_bins[i][j];
            double RL_max = RL_bins[i][j+1];

            TH1D * obs_hist = getHist(fin, observable, pt_min, pt_max, RL_min, RL_max, norm_string);
            temp_vec.push_back(obs_hist);
        } 
        obs_vec.push_back(temp_vec);
        // delete temp_vec here?? deleteVecOfHists(temp_vec);
    }

    // take ratios of different regions
    for (int i=0; i<n_bins; i++) {
        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];

        std::vector<TH1D*> ratio_of_dif_regions_overnonpert_vec;
        std::vector<TH1D*> ratio_of_dif_regions_overpert_vec;
        std::vector<std::string> label_overnonpert_vec;
        std::vector<std::string> label_overpert_vec;
        std::vector<Double_t> color_overnonpert_vec;
        std::vector<Double_t> color_overpert_vec;
        // cout << "LOOK HERE⁄!!!!" << obs_vec[0].size() << endl;
        for (int k=0; k<obs_vec[0].size(); k++) {
            int j = k;
            if (!include_RL0) {
                j = k+1;
            }
            // cout << "j: " << j << " k: " << k << endl;
            
            double RL_min = RL_bins[i][j];
            double RL_max = RL_bins[i][j+1];
            // std::string RLname_leg = Form("p_{T}R_{L} = %.1f-%.1f", RL_min, RL_max);
            std::string RLname_leg = Form("A = %.1f, B = %.1f", RL_min, RL_max);

            if (k != 0) {
                TH1D * hratio = takeratio(obs_vec[i][k], obs_vec[i][0]);
                ratio_of_dif_regions_overnonpert_vec.push_back(hratio);

                label_overnonpert_vec.push_back(RLname_leg);
                color_overnonpert_vec.push_back(colors[j]);
            }
            
            if (k != obs_vec[0].size()-1) {
                TH1D * hratio = takeratio(obs_vec[i][obs_vec[0].size()-1], obs_vec[i][k]);
                ratio_of_dif_regions_overpert_vec.push_back(hratio);

                label_overpert_vec.push_back(RLname_leg);
                color_overpert_vec.push_back(colors[j]);
            }
            
            
        }

        // now plot
        std::string ptname = to_string(pt_min) + "-" + to_string(pt_max);
        
        std::string obs_filename = observable;
        if (observable == "p") obs_filename = "deltap";
        if (observable == "jt") obs_filename = "deltajt";
        std::string output_filepath = outputbase + ptname + "/" + norm_string + "/" + obs_filename + "/"; 
        std::string filename = Form("%sratios_of_dif_regions_overnonpert_%s_PTBIN%d.pdf", output_filepath.c_str(), observable.c_str(), i);
        plot_ratios(ratio_of_dif_regions_overnonpert_vec, label_overnonpert_vec, color_overnonpert_vec, filename, fout, xtitle, "#frac{A #leq p_{T}R_{L} < B}{0.2 #leq p_{T}R_{L} < 0.8}"); // this will be wrong if highest rl bin changes

        filename = Form("%sratios_of_dif_regions_overpert_%s_PTBIN%d.pdf", output_filepath.c_str(), observable.c_str(), i);
        plot_ratios(ratio_of_dif_regions_overpert_vec, label_overpert_vec, color_overpert_vec, filename, fout, xtitle, Form("#frac{10 #leq p_{T}R_{L} < 30}{A #leq p_{T}R_{L} < B}")); //, RL_bins[i][1], RL_bins[i][2])); // this will be wrong if highest rl bin changes
        

    }

    // take ratios of different jet pt bins
    for (int j=0; j<n_RLbins; j++) {
        int k = j;
        if (!include_RL0) {
            k = j-1;
            if (j == 0) continue; // can add something here to change the filename for ALL
        }
        if (!include_RL1 && j == n_RLbins-1) continue;

        std::vector<TH1D*> ratio_of_dif_jetptbins_vec;
        std::vector<std::string> label_vec;
        std::vector<Double_t> color_vec = {kTeal+2, kRed};
        for (int i=1; i<n_bins; i++) {
            int pt_min = pt_bins[i];
            int pt_max = pt_bins[i+1];

            TH1D * hratio = takeratio(obs_vec[i][k], obs_vec[0][k]);
            ratio_of_dif_jetptbins_vec.push_back(hratio);

            std::string ptname_leg = Form("p_{T, jet} = %d-%d", pt_min, pt_max);
            label_vec.push_back(ptname_leg);
        }

        // now plot        
        std::string output_filepath = outputbase + "jetpt_ratios" + "/" + norm_string + "/"; 
        std::string filename = Form("%sratios_of_dif_jetptbins_%s_RLBIN%d.pdf", output_filepath.c_str(), observable.c_str(), k);
        plot_ratios(ratio_of_dif_jetptbins_vec, label_vec, color_vec, filename, fout, xtitle, "data/20-40 jet p_{T}");

        
    }




}


//===========================================================================
// This file takes in the following arguments:
// [[nothing]]
// THIS IS GOOD FOR PERLMUTTER
//===========================================================================
void takingnames_takingratios() { 

    gStyle->SetOptStat(0); // hide stats panel
    SetStyle();

    // ???fitfunctions();

    bool include_RL0 = false;
    bool include_RL1 = false;


    
    const int pt_bins[] = { 20, 40, 60, 80 };
    const int n_bins = sizeof(pt_bins) / sizeof(pt_bins[0]) - 1; //3;
    
    double RL_bins[3][7] = { { 0, 1e-2, 7e-2, 1.5e-1, 3e-1, 4e-1, 1 },
                            { 0, 1e-2, 4e-2, 8e-2, 2.5e-1, 4e-1, 1 },
                            { 0, 1e-2, 3e-2, 4.5e-2, 2e-1, 4e-1, 1 } };
    int n_RLbins = sizeof(RL_bins[0]) / sizeof(RL_bins[0][0]) - 1; //gets the columns //6; //7; //5;

    if (ptrl_bins == true) {
        // RL_bins = { { 0, 2e-1, 8e-1, 5.0, 10, 30, 100 },
        //                 { 0, 2e-1, 8e-1, 5.0, 10, 30, 100 },
        //                 { 0, 2e-1, 8e-1, 5.0, 10, 30, 100 } };
        for (int a=0; a<3; a++) {
            RL_bins[a][1] = 2e-1;
            RL_bins[a][2] = 8e-1;
            RL_bins[a][3] = 5.0;
            RL_bins[a][4] = 10.0;
            RL_bins[a][5] = 30.0;
            RL_bins[a][6] = 100.0;
        }
    }
    
    for (int a=0; a<3; a++) {
        for (int b=0; b<7; b++) {
            cout << RL_bins[a][b] << " ";
        }
        cout << endl;
    }



    
    // filenames
    // TString input_histograms_filename = Form("/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/%s/DataHists.root",attempt_dir.c_str()); // raw data
    TString input_histograms_filename = Form("/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/%s/DataHists_BinByBinCorr.root",attempt_dir.c_str()); // corr data

    TFile* root_data_file = new TFile(input_histograms_filename, "READ");

    
    // Output file with corrected results
    // std::string outfile = "/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/" + attempt_dir + "/DataHists.root";
    std::string outfile = "/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/" + attempt_dir + "/DataRatios.root"; 
    TFile* root_outfile = new TFile(outfile.c_str(), "RECREATE");

    // analyze(root_data_file, root_mc_file, root_outfile, pt_bins, n_bins, RL_bins, n_RLbins, "deltap", weightstr, jetRname, thrname, include_RL0, include_RL1);
    std::string norm_string = "self_normalized";
    // if (self_normalized_bool)
    if (deltap_bool) analyze(root_data_file, root_outfile, "deltap", n_bins, pt_bins, n_RLbins, RL_bins, norm_string, "#Deltap", include_RL0, include_RL1);
    if (p_bool) analyze(root_data_file, root_outfile, "p", n_bins, pt_bins, n_RLbins, RL_bins, norm_string, "p", include_RL0, include_RL1);
    if (deltajt_bool) analyze(root_data_file, root_outfile, "deltajt", n_bins, pt_bins, n_RLbins, RL_bins, norm_string, "#Deltaj_{T}", include_RL0, include_RL1);
    if (jt_bool) analyze(root_data_file, root_outfile, "jt", n_bins, pt_bins, n_RLbins, RL_bins, norm_string, "j_{T}", include_RL0, include_RL1);
    if (ew_bool) analyze(root_data_file, root_outfile, "weights", n_bins, pt_bins, n_RLbins, RL_bins, norm_string, "#frac{p_{T,1}p_{T,2}}{p_{T, jet}^{2}}", include_RL0, include_RL1);

    root_data_file->Close();
    root_outfile->Close();

    delete root_data_file;
    delete root_outfile;


}
