

// This file finds the bin by bin corrections and applies them to data.
// This is to be run on perlmutter.
// outputted plots:
    // truth_and_det_[obs_PTBINX_RLBINY].pdf
    // ratio_det_truth_[obs_PTBINX_RLBINY].pdf
    // corr_data_[obs_PTBINX_RLBINY].pdf
    // corr_data_and_raw_data_[obs_PTBINX_RLBINY].pdf
// Beatrice Liang-Gilman, beatrice_lg@berkeley.edu


// double rebin = 4;
bool ptrl_bins = true;
bool less_RLbins = true; // set true for using 3 "RL" bins, set false for using 5
std::string attempt_dir = ""; //Form("binbybincorrections/rebinx%.0f",rebin);


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


void formathist(TH1D * hist, std::string xtitle, std::string ytitle) {
    hist->GetXaxis()->SetTitle(xtitle.c_str());
    hist->GetYaxis()->SetTitle(ytitle.c_str());
}

// filetype 1 = ratio_det_truth_[obs_PTBINX_RLBINY].pdf
// filetype 2 = corr_data_[obs_PTBINX_RLBINY].pdf
// filetype 3 = fcorr_from_fit_[obs_PTBINX_RLBINY].pdf
// filetype 4 = ratio_det_truth_LINFIT_[obs_PTBINX_RLBINY].pdf
// filetype 5 = ratio_det_truth_QUADFIT_[obs_PTBINX_RLBINY].pdf
// filetype 6 = ratio_det_truth_EXPOFIT_[obs_PTBINX_RLBINY].pdf
void plot_and_save_one_histogram(TCanvas * can, TFile * file, std::string obsname,
                                 TH1D * hist1, int filetype, std::string addname,
                                 TPaveText * textbox = nullptr) {

    hist1->SetMarkerColor(kBlack);
    if (filetype == 1) hist1->GetYaxis()->SetRangeUser(0,5); //10);
    if (filetype == 3) hist1->GetYaxis()->SetRangeUser(0,1.5); //10);
    
    can->cd();
    // gPad->SetLogy();
    hist1->Draw("same");
    if (textbox != nullptr) textbox->Draw("same");

    std::string filename = "";
    if (filetype == 1) filename = "ratio_det_truth_";
    else if (filetype == 2) filename = "corr_data_";
    else if (filetype == 3) filename = "fcorr_from_fit_";
    else if (filetype == 4) filename = "ratio_det_truth_LINFIT_";
    else if (filetype == 5) filename = "ratio_det_truth_QUADFIT_";
    else if (filetype == 6) filename = "ratio_det_truth_EXPOFIT_";
    std::string outputbase = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/plots/";
    std::string outputname = outputbase + attempt_dir + "/" + obsname + "/" + filename + obsname + addname + ".pdf";
    can->SaveAs(outputname.c_str());

    if (filetype < 3) {
        file->cd();
        hist1->Write();
    }

    delete can;

}


// filetype 1 = truth_and_det_[obs_PTBINX_RLBINY].pdf
// filetype 2 = corr_data_and_raw_data_[obs_PTBINX_RLBINY].pdf
void plot_and_save_two_histograms_overlayed(TCanvas * can, TFile * file, std::string obsname,
                                            TH1D * hist1, TH1D * hist2, int filetype, std::string addname,
                                            int markercolor1, int markercolor2, std::string label1, std::string label2) {

    if (obsname == "charge") {
        hist1->SetMinimum(0);
        double hist_max = std::max(hist1->GetMaximum(), hist2->GetMaximum()) * 1.1;
        hist1->SetMaximum(hist_max);
    } else {
        double hist_max = std::max(hist1->GetMaximum(), hist2->GetMaximum()) * 2;
        hist1->SetMaximum(hist_max);
    }
    
    
    hist1->SetMarkerColorAlpha(markercolor1, 0.8);
    hist1->SetLineColorAlpha(markercolor1, 0.8);
    hist2->SetMarkerColorAlpha(markercolor2, 0.8);
    hist2->SetLineColorAlpha(markercolor2, 0.8);

    TLegend *l = new TLegend();
    l->AddEntry(hist1, label1.c_str(), "pl");
    l->AddEntry(hist2, label2.c_str(), "pl");


    can->cd();
    if (obsname != "charge") gPad->SetLogy();
    hist1->Draw("same");
    hist2->Draw("same");
    l->Draw("same");

    std::string filename = "";
    if (filetype == 1) filename = "truth_and_det_";
    else if (filetype == 2) filename = "corr_data_and_raw_data_";
    std::string outputbase = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/plots/";
    std::string outputname = outputbase + attempt_dir + "/" + obsname + "/" + filename + obsname + addname + ".pdf";
    can->SaveAs(outputname.c_str());

    if (filetype != 2) {
        file->cd();
        hist1->Write();
        hist2->Write();
    }

    delete can;

}


// filetype 1 = ratio_det_truth_[obs_PTBINX_RLBINY].pdf
// filetype 2 = corr_data_[obs_PTBINX_RLBINY].pdf
// filetype 3 = ratio_det_truth_FIT_[obs_PTBINX_RLBINY].pdf
void plot_and_save_one_graph(TCanvas * can, TFile * file, std::string obsname,
                             TGraphErrors * graph1, int filetype, std::string addname) {

    graph1->SetMarkerColor(kBlack);
    // if (filetype == 1) {graph1->GetYaxis()->SetRangeUser(0,10);}
    
    can->cd();
    // gPad->SetLogy();
    graph1->Draw("ALP");

    std::string filename = "";
    if (filetype == 1) filename = "ratio_det_truth_";
    else if (filetype == 2) filename = "corr_data_";
    else if (filetype == 3) filename = "ratio_det_truth_FIT_";
    std::string outputbase = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/plots/";
    std::string outputname = outputbase + attempt_dir + "/" + obsname + "/" + filename + obsname + addname + ".pdf";
    can->SaveAs(outputname.c_str());

    if (filetype != 3) {
        file->cd();
        graph1->Write();
    }

    delete can;

}

// filetype 1 = truth_and_det_[obs_PTBINX_RLBINY].pdf
// filetype 2 = corr_data_and_raw_data_[obs_PTBINX].pdf
void plot_and_save_two_graphs_overlayed(TCanvas * can, TFile * file, std::string obsname,
                                        TGraphErrors * graph1, TGraphErrors * graph2, int filetype, std::string addname,
                                        int markercolor1, int markercolor2, std::string label1, std::string label2) {

    if (obsname == "rc") {
        // double hist_min = std::min(graph1->GetMinimum(), graph2->GetMinimum());
        double minY = graph1->GetY()[0];
        for (int i = 0; i < graph1->GetN(); ++i) {
            double currentY = graph1->GetPointY(i); // Get y-value at point i
            if (currentY < minY) minY = currentY;
            currentY = graph2->GetPointY(i);
            if (currentY < minY) minY = currentY;
        }
        graph1->SetMinimum(minY*1.1);
        double hist_max = std::max(graph1->GetMaximum(), graph2->GetMaximum()) > 0.1 ? std::max(graph1->GetMaximum(), graph2->GetMaximum()) : 0.1;
        graph1->SetMaximum(hist_max);
    }
    
    graph1->SetMarkerColorAlpha(markercolor1, 0.8);
    graph1->SetLineColorAlpha(markercolor1, 0.8);
    graph2->SetMarkerColorAlpha(markercolor2, 0.8);
    graph2->SetLineColorAlpha(markercolor2, 0.8);

    TLegend *l = new TLegend();
    l->AddEntry(graph1, label1.c_str(), "pl");
    l->AddEntry(graph2, label2.c_str(), "pl");


    can->cd();
    // gPad->SetLogy();
    graph1->Draw("ALP SAME");
    graph2->Draw("LP SAME");
    l->Draw("same");

    std::string filename = "";
    if (filetype == 1) filename = "truth_and_det_";
    else if (filetype == 2) filename = "corr_data_and_raw_data_";
    std::string outputbase = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/plots/";
    std::string outputname = outputbase + attempt_dir + "/" + obsname + "/" + filename + obsname + addname + ".pdf";
    can->SaveAs(outputname.c_str());

    if (filetype != 2) {
        file->cd();
        graph1->Write();
        graph2->Write();
    }

    delete can;

}


//===========================================================================
//================================= FITTING =================================
//===========================================================================
TF1 * fit_histogram_linearfit(TH1D * hist, std::string observable, std::string addname, double fitmax) {

    /*
    // could also do:
    TF1 *linearFit = new TF1("linearFit", "pol1", 2.0, 8.0); // Fit range: [2, 8]
    hist->Fit(linearFit, "R"); // "R" = ensures the fit uses the specified range
    */

    // Perform the linear fit
    // hist->Fit("pol1", "R"); //, "Q"); // "pol1" = linear function, "Q" = quiet mode
    // TF1 *fitFunction = hist->GetFunction("pol1");
    // auto fitResult = hist->Fit("pol1", "S");

    TF1 *fitFunction = new TF1("fitFunction", "pol1", 0.0, fitmax); // Fit range: [2, 8]
    fitFunction->SetParameter(1,0.005); // set initial slope parameter to 0.005
    auto fitResult = hist->Fit(fitFunction, "SR");
    // "W": Ignore weights - set the weights of all non-zero bins to 1
    // "E": Perform better error estimation.
    // "S": The full result of the fit is returned in the TFitResultPtr (incl cov matrix)
    // 7.1.1. in https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html
    
    // Retrieve fit parameters
    double slope = fitFunction->GetParameter(1);  // Slope of the line
    double intercept = fitFunction->GetParameter(0); // Intercept of the line
    double chi2 = fitResult->Chi2(); //or fitFunction->GetChisquare();
    int ndf = fitResult->Ndf(); //or fitFunction->GetNDF();
    double pvalue = fitFunction->GetProb();

    // std::cout << "Fit Results: slope = " << slope << ", intercept = " << intercept << std::endl;
    // std::cout << "Chi2/Ndf = " << chi2 / ndf << std::endl;

    

    // Create a TPaveText (x1, y1, x2, y2 in NDC coordinates)
    TPaveText* pavetext = new TPaveText(0.2, 0.6, 0.5, 0.8, "NDC"); // Coordinates in normalized device space
    // pavetext->SetFillColor(0);         // Set background color (0 for transparent)
    // pavetext->SetTextColor(1);         // Set text color
    // pavetext->SetTextSize(0.03);       // Set text size
    // pavetext->SetBorderSize(1);        // Set border size

    // Add multiple lines
    pavetext->AddText(Form("Fit results: y = %.3fx + %.3f", slope, intercept));
    pavetext->AddText(Form("#Chi^{2}/Ndf = %.4f", chi2 / ndf));
    pavetext->AddText(Form("p-value = %.4f", pvalue));

    // Draw the histogram and the fit
    TCanvas *c1 = new TCanvas("c1", "Linear Fit", 800, 600);
    TFile *dummy_file;
    plot_and_save_one_histogram(c1, dummy_file, observable, hist, 4, addname, pavetext);

    return fitFunction;
}

TF1 * fit_histogram_quadfit(TH1D * hist, std::string observable, std::string addname, double fitmax) {

    TF1 *fitFunction = new TF1("fitFunction", "pol2", 0.0, fitmax); // Fit range: [2, 8]
    auto fitResult = hist->Fit(fitFunction, "SR");

    // Retrieve fit parameters
    double p0 = fitFunction->GetParameter(0); // a in ax^2 + bx + c
    double p1 = fitFunction->GetParameter(1); // b in ax^2 + bx + c
    double p2 = fitFunction->GetParameter(2); // b in ax^2 + bx + c
    double chi2 = fitResult->Chi2(); //or fitFunction->GetChisquare();
    int ndf = fitResult->Ndf(); //or fitFunction->GetNDF();
    double pvalue = fitFunction->GetProb();

    // Create a TPaveText (x1, y1, x2, y2 in NDC coordinates)
    TPaveText* pavetext = new TPaveText(0.2, 0.6, 0.5, 0.8, "NDC"); // Coordinates in normalized device space

    // Add multiple lines
    pavetext->AddText(Form("Fit results: y = %.5fx^{2} + %.3fx + %.3f", p2, p1, p0));
    pavetext->AddText(Form("#Chi^{2}/Ndf = %.4f", chi2 / ndf));
    pavetext->AddText(Form("p-value = %.4f", pvalue));

    // Draw the histogram and the fit
    TCanvas *c1 = new TCanvas("c1", "Quadratic Fit", 800, 600);
    TFile *dummy_file;
    plot_and_save_one_histogram(c1, dummy_file, observable, hist, 5, addname, pavetext);

    return fitFunction;

}

TF1 * fit_histogram_expofit(TH1D * hist, std::string observable, std::string addname, double fitmax) {

    TF1 *fitFunction = new TF1("fitFunction", "expo", 0.0, fitmax); // expo: f(x) = exp(p0+p1*x)
    auto fitResult = hist->Fit(fitFunction, "SR");

    // Retrieve fit parameters
    double p0 = fitFunction->GetParameter(0);  // a in ax^2 + bx + c
    double p1 = fitFunction->GetParameter(1); // b in ax^2 + bx + c
    double chi2 = fitResult->Chi2(); //or fitFunction->GetChisquare();
    int ndf = fitResult->Ndf(); //or fitFunction->GetNDF();
    double pvalue = fitFunction->GetProb();

    // Create a TPaveText (x1, y1, x2, y2 in NDC coordinates)
    TPaveText* pavetext = new TPaveText(0.2, 0.6, 0.5, 0.8, "NDC"); // Coordinates in normalized device space

    // Add multiple lines
    pavetext->AddText(Form("Fit results: y = exp(%.3fx + %.3f)", p1, p0));
    pavetext->AddText(Form("#Chi^{2}/Ndf = %.4f", chi2 / ndf));
    pavetext->AddText(Form("p-value = %.4f", pvalue));

    // Draw the histogram and the fit
    TCanvas *c1 = new TCanvas("c1", "Exponential Fit", 800, 600);
    TFile *dummy_file;
    plot_and_save_one_histogram(c1, dummy_file, observable, hist, 6, addname, pavetext);

    return fitFunction;
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
void getthemeans() {

}

void getHist(TFile * file, std::string observable, int pt_min, int pt_max, 
             double RL_min, double RL_max) {
    
    std::string histname = Form("h_%s_R0.4_t1.0_pt%d-%d_pTRL%.3f-%.3f_self_normalized", observable, pt_min, pt_max, RL_min, RL_max);
    TH1D * hist = (TH1D * )file->Get(histname.c_str());

    return hist;
}

//===========================================================================
// This function takes the ratios of the functions.
//===========================================================================
void takeratios(std::string observable, 
                bool debug = false) {
    
    vector<TH1D *> hists_vec;
    vector<TH1D *> ratios_vec;

    for (int i=0; i<n_bins; i++) {
        for (int j=0; j<n_RLbins; j++) {
            TH1D * obs_hist = getHist(observable, pt_min, pt_max, RL_min, RL_max);
            hists_vec.push_back(obs_hist);

            if (j>0) { // check this!! should j>1?? or k>0?
                TH1D * hratio = hists_vec[0]->Clone();
                hratio->Divide(hists_vec[j]); // also check if this should be [j]
                ratios_vec.push_back(hratio);
            }
        }
    }

    //format the ratio plots, add a legend
    
    // plot the ratio plots
    TCanvas * can_ratios = new TCanvas();

}


void analyze(double RL_min, double RL_max) {
    //histname = h_deltap_R0.4_t1.0_pt20-40_pTRL0.200-0.800_self_normalized

    getthemeans();
    takeratios();

}


//===========================================================================
// This file takes in the following arguments:
// [[nothing]]
// THIS IS GOOD FOR PERLMUTTER
//===========================================================================
void motivatethis() { 

    gStyle->SetOptStat(0); // hide stats panel
    SetStyle();


    // attempt_dir name
    // std::string matched_str = "matched";
    // if (unmatched) matched_str = "unmatched";
    // if (anchmc) attempt_dir = Form("binbybincorrections/%s/rebinx%.0f", matched_str.c_str(), rebin);
    // else attempt_dir = Form("fastsim_binbybincorrections/%s/rebinx%.0f", matched_str.c_str(), rebin);
    
    
    const int pt_bins[] = { 20, 40, 60, 80 };
    const int n_bins = sizeof(pt_bins) / sizeof(pt_bins[0]) - 1; //3;
    
    const double RL_bins[3][8] = { { 0, 1e-2, 3e-2, 7e-2, 1.5e-1, 3e-1, 4e-1, 1 },
                            { 0, 1e-2, 2.5e-2, 4e-2, 8e-2, 2.5e-1, 4e-1, 1 },
                            { 0, 1e-2, 2.5e-2, 3e-2, 4.5e-2, 2e-1, 4e-1, 1 } };
    const int n_RLbins = sizeof(RL_bins[0]) / sizeof(RL_bins[0][0]) - 1; //gets the columns //6; //7; //5;

    if (less_RLbins == true) {
        n_RLbins = 5; // (assuming that first and last bins get skipped)
        // colors[16] = {kGray, kMagenta, kGreen+2, kBlue, kOrange+1, kViolet+1, kRed, kYellow+1, kCyan+1};
        colors[2] = kBlue; colors[3] = kViolet+1; colors[5] = kGreen+2;
        if (ptrl_bins == false) {
            RL_bins[0][2] = 4e-2; RL_bins[1][2] = 2e-2; RL_bins[2][2] = 2e-2;
            RL_bins[0][3] = 1.5e-1; RL_bins[1][3] = 8e-2; RL_bins[2][3] = 7e-2;
            for (int a=0; a<3; a++) {
                RL_bins[a][4] = 4e-1; 
                RL_bins[a][5] = 1; 
                RL_bins[a][6] = -1; 
                RL_bins[a][7] = -1;
            }
        } else { //3 ptrl bins
            for (int a=0; a<3; a++) {
                RL_bins[a][1] = 2e-1;
                RL_bins[a][2] = 8e-1;
                RL_bins[a][3] = 5.0;
                RL_bins[a][4] = 30.0;
                RL_bins[a][5] = 100.0;
                RL_bins[a][6] = -1;
                RL_bins[a][7] = -1;
            }
        }
    }
    
    for (int a=0; a<3; a++) {
        for (int b=0; b<8; b++) {
            cout << RL_bins[a][b] << " ";
        }
        cout << endl;
    }

    
    // filenames
    // NEED TO CORRECT THIS! (like apply corrections and then fix this path)
    TString input_histograms_filename = "/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/data_thirdattempt_ptrlbins/rebinx4";
    TFile* root_data_file = new TFile(input_histograms_filename, "READ");

    
    // Output file with corrected results
    // std::string outfile = "/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/" + attempt_dir + "/DataHists.root";
    std::string outfile = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/rootfiles/" + attempt_dir + "/DataHists_BinByBinCorr.root"; 
    TFile* root_outfile = new TFile(outfile.c_str(), "RECREATE");

    // analyze(root_data_file, root_mc_file, root_outfile, pt_bins, n_bins, RL_bins, n_RLbins, "deltap", weightstr, jetRname, thrname, include_RL0, include_RL1);
    
    


    root_data_file->Close();
    root_outfile->Close();

    delete root_data_file;
    delete root_outfile;


}
