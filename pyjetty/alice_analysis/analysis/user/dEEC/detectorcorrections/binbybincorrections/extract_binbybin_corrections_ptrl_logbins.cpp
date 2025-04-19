// This file finds the bin by bin corrections and applies them to data.
// This is to be run on perlmutter.
// outputted plots:
    // truth_and_det_[obs_PTBINX_RLBINY].pdf
    // ratio_det_truth_[obs_PTBINX_RLBINY].pdf
    // corr_data_[obs_PTBINX_RLBINY].pdf
    // corr_data_and_raw_data_[obs_PTBINX_RLBINY].pdf
// Beatrice Liang-Gilman, beatrice_lg@berkeley.edu


std::string compsystem = "perlmutter"; //"local"; //"perlmutter"
bool unmatched = true; // set true for unmatched, false for matched
double rebin = 4;
bool anchmc = true; // set true for anchored mc, set false for fastsim

bool self_normalize = true; // DO NOT CHANGE // -- this from any data_thirdattempt
bool ptrl_bins = true; // DO NOT CHANGE
bool less_RLbins = true; // DO NOT CHANGE // set true for using 3 "RL" bins, set false for using 5

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
    if (compsystem == "local") outputbase = "/Volumes/WORK USB/dEEC/storage/plots/";
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
    if (compsystem == "local") outputbase = "/Volumes/WORK USB/dEEC/storage/plots/";
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
    if (compsystem == "local") outputbase = "/Volumes/WORK USB/dEEC/storage/plots/";
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
    if (compsystem == "local") outputbase = "/Volumes/WORK USB/dEEC/storage/plots/";
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
// This function takes in a fit function, and outputs a histogram with the  
// correction factors in each bin.
//===========================================================================
TH1D * extract_fcorr_from_fitfunction(TF1 * fitfunc, std::vector<double> bincenter_values,
                                      std::string observable, std::string add_name) {

    // std::vector<double> fcorr_values;

    // this is assuming equal bin sizes:
    int numbins = bincenter_values.size();
    double binwidth = bincenter_values[1] - bincenter_values[0];
    double hist_lowx = bincenter_values[0] - binwidth/2;
    double hist_highx = bincenter_values[numbins-1] + binwidth/2;

    // TODO: rename this histogram!
    std::string histname = Form("hist_fcorr_fromfit_%s_%s", observable.c_str(), add_name.c_str());
    TH1D * hist_fcorr_fromfit = new TH1D(histname.c_str(), histname.c_str(), numbins, hist_lowx, hist_highx);

    // Loop through x-values, evaluate y, and store the results
    for (double bincenter : bincenter_values) {
        double y = fitfunc->Eval(bincenter);
        // fcorr_values.push_back(y);
        // std::cout << "bincenter = " << bincenter << ", y = " << y << std::endl;

        hist_fcorr_fromfit->Fill(bincenter, y);
    }
    return hist_fcorr_fromfit;
}

//===========================================================================
// This function gets the reco (det) level observable and the gen (truth) level 
// observable, and takes the ratio of the two.
// It plots this ratio and saves it as a separate histogram.
//===========================================================================
TH1D * get_binbybin_corrfactors(TFile *fin_mc, TFile *fout, std::string observable, 
                                int ptbin, int pt_min, int pt_max, int rlbin, double RL_min, double RL_max,
                                bool scalebyRLbinwidth,
                                vector<double>& RL_vals_mc_det = *(new std::vector<double>()), vector<double>& RL_vals_mc_truth = *(new std::vector<double>())) {

    fin_mc->cd();

    // find observable maximums
    int obs_max_val = 0;
    if (observable == "deltap" || observable == "deltapt") obs_max_val = pt_max+5;
    else if (observable == "deltapl") obs_max_val = pt_max/2;
    else if (observable == "deltajt") obs_max_val = 5;
    else if (observable == "deltajl") obs_max_val = pt_max+5;

    std::string histname;
    TH1D * hist_det;
    TH1D * hist_truth;
    
    if (unmatched) { // this is the unmatched version:
        //det level
        histname = Form("reco_%s_unmatched_PTBIN%dScaled", observable.c_str(), ptbin);
        cout << " HISTNAME: " << histname << endl;
        TH3D *h3D_reco = (TH3D *) fin_mc->Get(histname.c_str()); // axes: obs, RL, jet pt
        cout << "RL MIN MAX" << RL_min << " " << RL_max << endl;
        // cout << "x axis " << h3D_reco->GetXaxis()->GetXmin() << " and " << h3D_reco->GetXaxis()->GetXmax() << endl;
        // cout << "y axis " << h3D_reco->GetYaxis()->GetXmin() << " and " << h3D_reco->GetYaxis()->GetXmax() << endl;
        // cout << "z axis " << h3D_reco->GetZaxis()->GetXmin() << " and " << h3D_reco->GetZaxis()->GetXmax() << endl;
        
        // make a clone
        
        // extract 1D histograms
        h3D_reco->GetZaxis()->SetRangeUser(pt_min, pt_max);
        h3D_reco->GetYaxis()->SetRangeUser(RL_min, RL_max);
        h3D_reco->GetXaxis()->SetRangeUser(0, obs_max_val); // to limit the observable range displayed (needed for divide later)
        // hist_det = (TH1D *) h3D_reco->ProjectionX();
        // ("x", xMin, xMax, yMin, yMax, zMin, zMax)
        hist_det = (TH1D * ) h3D_reco->Project3D("x"); //, 0, obs_max_val, RL_min, RL_max, pt_min, pt_max);

        //truth level
        histname = Form("gen_%s_unmatched_PTBIN%dScaled", observable.c_str(), ptbin);
        cout << " HISTNAME: " << histname << endl;
        TH3D *h3D_gen = (TH3D *) fin_mc->Get(histname.c_str()); // axes: obs, RL, jet pt

        // extract 1D histograms
        h3D_gen->GetZaxis()->SetRangeUser(pt_min, pt_max);
        h3D_gen->GetYaxis()->SetRangeUser(RL_min, RL_max);
        h3D_gen->GetXaxis()->SetRangeUser(0, obs_max_val); // to limit the observable range displayed (needed for divide later)
        // hist_truth = (TH1D *) h3D_gen->ProjectionX();
        hist_truth = (TH1D * ) h3D_gen->Project3D("x"); //, 0, obs_max_val, RL_min, RL_max, pt_min, pt_max);


        delete h3D_reco;
        delete h3D_gen;

    } else { // this is the matched version:
        histname = Form("hResponse_JetPt_corr_%s_PTBIN%d_RLBIN%d_R0.4_1.0Scaled", observable.c_str(), ptbin, rlbin);
        cout << " HISTNAME: " << histname << endl;
        THnSparse *hsparse = (THnSparse *) fin_mc->Get(histname.c_str());

        // extract 1D histograms
        hsparse->GetAxis(0)->SetRangeUser(pt_min, pt_max);
        hsparse->GetAxis(2)->SetRangeUser(0, obs_max_val); // to limit the observable range displayed (needed for divide later)
        hist_det = (TH1D *) hsparse->Projection(2);

        hsparse->GetAxis(1)->SetRangeUser(pt_min, pt_max);
        hsparse->GetAxis(3)->SetRangeUser(0, obs_max_val); // to limit the observable range displayed
        hist_truth = (TH1D *) hsparse->Projection(3);
        
        // not doing but should maybe do for completeness: scale by the # of jets. Didn't do bc in matched # truth jets = # det jets. Not using matched anymore so didn't bother implementing.
    }

    // rebin
    if (observable != "charge" && ptrl_bins == false) {
        hist_det->Rebin(rebin);
        hist_truth->Rebin(rebin);
    }

    // scale by the number of jets here
    if (unmatched) {
        TH1D * jetpt_det = (TH1D *) fin_mc->Get("h_1Djet_pt_JetPt_Det_R0.4_1.0Scaled");
        TH1D * jetpt_truth = (TH1D *) fin_mc->Get("h_1Djet_pt_JetPt_Truth_R0.4_1.0Scaled");

        jetpt_det->GetXaxis()->SetRangeUser(pt_min, pt_max);
        jetpt_truth->GetXaxis()->SetRangeUser(pt_min, pt_max);

        double num_jets_det = jetpt_det->Integral();
        double num_jets_truth = jetpt_truth->Integral();

        if (self_normalize) {
            num_jets_det = hist_det->Integral(); // naming is bad, this is really num_pairs_det
            num_jets_truth = hist_truth->Integral(); // naming is bad, this is really num_pairs_truth
        }

        hist_det->Scale(num_jets_det, "width");
        hist_truth->Scale(num_jets_truth, "width");
    }

    // scale by the RL bin width
    double RL_bin_width = RL_max - RL_min;
    if ( scalebyRLbinwidth ) hist_det->Scale(RL_bin_width);
    if ( scalebyRLbinwidth ) hist_truth->Scale(RL_bin_width);
        
    // add name string
    std::string addname = Form("_PTBIN%d_RLBIN%d", ptbin, rlbin);

    // add axes titles
    std::string obs_axis_title = "";
    if (observable == "deltap") obs_axis_title = "#Deltap";
    else if (observable == "deltapt") obs_axis_title = "#Deltap_{T}";
    else if (observable == "deltapl") obs_axis_title = "#Deltap_{L}";
    else if (observable == "charge") obs_axis_title = "q_{1}q_{2}";
    else if (observable == "deltajt") obs_axis_title = "#Deltaj_{T}";
    else if (observable == "deltajl") obs_axis_title = "#Deltaj_{L}";
    std::string yaxis_title = Form("#frac{1}{N_{jet}#DeltaR_{L}} #frac{dN}{d%s}", obs_axis_title.c_str());
    formathist(hist_det, obs_axis_title, yaxis_title);
    formathist(hist_truth, obs_axis_title, yaxis_title);

    // plot and save those histograms
    TCanvas *can_truth_det = new TCanvas();
    plot_and_save_two_histograms_overlayed(can_truth_det, fout, observable, hist_det, hist_truth, 1, addname, kBlue, kRed, "det level", "truth level");

    // need to get rc
    if (observable == "charge") {
        double rc_det = getRc(hist_det);
        double rc_truth = getRc(hist_truth);
        RL_vals_mc_det.push_back(rc_det);
        RL_vals_mc_truth.push_back(rc_truth);

        double fcorr_rc = rc_det/rc_truth;

        TH1D * hfcorr_rc = new TH1D("hfcorr_rc", "hfcorr_rc", 1, 0.0, 1.0);
        hfcorr_rc->Fill(0.0,fcorr_rc); // Fill histogram at 0 with a weight that is the correction factor
        return hfcorr_rc;
    }

    // find the bin by bin corrections
    TH1D * hratio = (TH1D *) hist_det->Clone(Form("hratio_%s_PTBIN%d_RLBIN%d", observable.c_str(), ptbin, rlbin));
    hratio->Divide(hist_truth);
    formathist(hratio, obs_axis_title, "f_{corr}");
    
    // plot and save ratio
    TCanvas *can_ratio = new TCanvas();
    plot_and_save_one_histogram(can_ratio, fout, observable, hratio, 1, addname);
    

    // fit the ratio
    double fit_max = pt_max;
    if (observable == "deltapl") fit_max = pt_max/2;
    if (observable == "deltajt") fit_max = 5;
    TF1 * linfit = fit_histogram_linearfit(hratio, observable, addname, fit_max);
    TF1 * quadfit = fit_histogram_quadfit(hratio, observable, addname, fit_max);
    TF1 * expofit = fit_histogram_expofit(hratio, observable, addname, fit_max);

    // extract the fcorr from fitted function, then return this new histogram instead of hratio
    // std::vector<double> bincenters_vector = get_bin_centers(hratio);
    // TH1D * hfcorr_from_fit = extract_fcorr_from_fitfunction(--, bincenters_vector, observable, addname);
    // TCanvas *can_fcorr_from_fit = new TCanvas();
    // plot_and_save_one_histogram(can_fcorr_from_fit, fout, observable, hfcorr_from_fit, 3, addname);

    

    return hratio;

}

// Get the raw data
TH1D * get_rawdata(TFile *fin_data, std::string observable, int pt_min, int pt_max, double RL_min, double RL_max) {

    fin_data->cd();
    std::string histname = Form("h_%s_R0.4_t1.0_pt%d-%d_RL%.3f-%.3f_norm_by_jets", observable.c_str(), pt_min, pt_max, RL_min, RL_max);
    TH1D * hist = (TH1D *)gDirectory->Get(histname.c_str())->Clone(histname.c_str());

    // rebin -- now already done in original DataHists file
    // hist->Rebin(rebin);

    return hist;

}

// Correct the raw data here
// Important!! Right now the bins need to be the same; otherwise need to rebin to make it match
void apply_corrfactor(TFile *fout, TH1D * h_raw_data, TH1D * h_binbybin_fcorr, std::string observable,
                      int ptbin, int rlbin, bool debug = false) {

    // h_raw_data->Sumw2();
    // h_binbybin_fcorr->Sumw2();

    if (debug) {
        cout << "  NEW HIST" << endl;
        cout << "num bins in raw" << h_raw_data->GetNbinsX() << endl;
        cout << "num bins in binbybin" << h_binbybin_fcorr->GetNbinsX() << endl;
        for (int i=0; i<h_raw_data->GetNbinsX(); i++) {
            cout << h_raw_data->GetBinContent(i) << " ";
        }
        cout << endl;
    }

    std::string newhistname = Form("%s_corrected", h_raw_data->GetName());
    TH1D *h_corr_data = (TH1D*) h_raw_data->Clone(newhistname.c_str());
    h_corr_data->Divide(h_binbybin_fcorr);

    if (debug) {
        cout << endl;
        for (int i=0; i<h_raw_data->GetNbinsX(); i++) {
            cout << h_raw_data->GetBinContent(i) << "/" << h_binbybin_fcorr->GetBinContent(i) << " = " << h_corr_data->GetBinContent(i) << "   ";
        }
        cout << endl;
        cout << "=====" << endl;
    }

    // add name string
    std::string addname = Form("_PTBIN%d_RLBIN%d", ptbin, rlbin);

    // plot and save corrected data
    TCanvas *can_corrdata = new TCanvas();
    plot_and_save_one_histogram(can_corrdata, fout, observable, h_corr_data, 2, addname);

    // plot and save corr data with raw data 
    TCanvas *can_corr_vs_raw_data = new TCanvas();
    plot_and_save_two_histograms_overlayed(can_corr_vs_raw_data, fout, observable, h_raw_data, h_corr_data, 2, addname, kBlue, kRed, "raw data", "bin-by-bin corrected data");


}

// General analysis function
void analyze_charge(TFile *f_in_data, TFile *f_in_mc, TFile *f_out, const int pt_bins[], int n_bins, const double RL_bins[][8],
             int n_RLbins, std::string observable, std::string weightstr, std::string jetRname, 
             std::string thrname, bool include_RL0, bool include_RL1) {
    
    
    for (int i = 0; i < n_bins; i++) {
        cout << "in pt bin" << i << endl;
        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];

        std::vector<double> RL_bin_centers_vec;
        std::vector<double> fcorr_rc_vec;
        std::vector<double> rc_vec;
        std::vector<double> rc_corr_vec;

        vector<double> RL_vals_mc_det;
        vector<double> RL_vals_mc_truth;


        for ( int j = 0; j < n_RLbins; j++ ) {
            int k = j;
            if (!include_RL0) {
                k = j-1;
                if (j == 0) continue; // can add something here to change the filename for ALL
            }
            if (!include_RL1 && j == n_RLbins-1) continue;

            double RL_min = RL_bins[i][j];
            double RL_max = RL_bins[i][j+1];

            // // add name string
            // std::string addname = Form("_PTBIN%d_RLBIN%d", i, k);

            
        
            // need to do something different for charge! pt cuts not appropriate

            // first we want to plot the det vs truth level same sign and opp sign. And get the correction factor
            // the correction factor is saved at hist->getbincontent(hist->findbin(0))
            TH1D * h_binbybin_fcorr = get_binbybin_corrfactors(f_in_mc, f_out, observable, i, pt_min, pt_max, k, RL_min, RL_max, false, RL_vals_mc_det, RL_vals_mc_truth);
            double fcorr_rc = h_binbybin_fcorr->GetBinContent(h_binbybin_fcorr->FindBin(0));
            fcorr_rc_vec.push_back(fcorr_rc);
            RL_bin_centers_vec.push_back((RL_min+RL_max)/2);

            // then get the data rc
            // TH1D * h_raw_data = get_rawdata(f_in_data, observable, pt_min, pt_max, RL_min, RL_max);
            // double rc_data = getRc(h_raw_data);
            // rc_vec.push_back(rc_data);
            f_in_data->cd();
            TGraphErrors * graph_rc = (TGraphErrors *) f_in_data->Get(Form("rc__R0.4_t1.0_pt%d-%d", pt_min, pt_max));
            double rc_data_val = graph_rc->GetPointY(k);
            rc_vec.push_back(rc_data_val);
            
            // then apply correction to data
            double rc_corr_data_val = rc_data_val / fcorr_rc;
            rc_corr_vec.push_back(rc_corr_data_val);

            //TODO: add in the errors?
        
        }

        // add name string
        std::string addname = Form("_PTBIN%d", i);

        
        // also plot the corr vs raw data here
        TGraphErrors * gr_rc_mc_det = new TGraphErrors(RL_bin_centers_vec.size(), RL_bin_centers_vec.data(), RL_vals_mc_det.data()); //add errors at end
        TGraphErrors * gr_rc_mc_truth = new TGraphErrors(RL_bin_centers_vec.size(), RL_bin_centers_vec.data(), RL_vals_mc_truth.data()); //add errors at end
        TGraphErrors * gr_fcorr = new TGraphErrors(RL_bin_centers_vec.size(), RL_bin_centers_vec.data(), fcorr_rc_vec.data()); //add errors at end
        TGraphErrors * gr_rc_raw = new TGraphErrors(RL_bin_centers_vec.size(), RL_bin_centers_vec.data(), rc_vec.data()); //add errors at end
        TGraphErrors * gr_rc_corr = new TGraphErrors(RL_bin_centers_vec.size(), RL_bin_centers_vec.data(), rc_corr_vec.data()); //add errors at end

        gr_rc_mc_det->SetNameTitle(Form("gr_rc_mc_det_%s", addname.c_str()), Form("gr_rc_mc_det_%s", addname.c_str()));
        gr_rc_mc_truth->SetNameTitle(Form("gr_rc_mc_truth_%s", addname.c_str()), Form("gr_rc_mc_truth_%s", addname.c_str()));
        gr_fcorr->SetNameTitle(Form("gr_fcorr_%s", addname.c_str()), Form("gr_fcorr_%s", addname.c_str()));
        gr_rc_raw->SetNameTitle(Form("gr_rc_raw_%s", addname.c_str()), Form("gr_rc_raw_%s", addname.c_str()));
        gr_rc_corr->SetNameTitle(Form("gr_rc_corr_%s", addname.c_str()), Form("gr_rc_corr_%s", addname.c_str()));
        

        // plot the truth vs det level distributions
        TCanvas *can_truth_det_rc = new TCanvas();
        plot_and_save_two_graphs_overlayed(can_truth_det_rc, f_out, "rc", gr_rc_mc_det, gr_rc_mc_truth, 1, addname, kBlue, kRed, "det level", "truth level");


        // plot the correction factors here, as a function of RL?
        TCanvas *can_rc_fcorr = new TCanvas();
        plot_and_save_one_graph(can_rc_fcorr, f_out, "rc", gr_fcorr, 1, addname);

        // plot the corrected data
        TCanvas *can_rc_corr = new TCanvas();
        plot_and_save_one_graph(can_rc_corr, f_out, "rc", gr_rc_corr, 2, addname);
        
        // plot the corrected vs raw data
        TCanvas *can_rc = new TCanvas();
        plot_and_save_two_graphs_overlayed(can_rc, f_out, "rc", gr_rc_raw, gr_rc_corr, 2, addname, kBlue, kRed, "raw data", "corrected data");

        delete gr_rc_mc_det;
        delete gr_rc_mc_truth;
        delete gr_fcorr;
        delete gr_rc_raw;
        delete gr_rc_corr;
        

    }

    

}

// General analysis function
void analyze(TFile *f_in_data, TFile *f_in_mc, TFile *f_out, const int pt_bins[], int n_bins, const double RL_bins[][8],
             int n_RLbins, std::string observable, std::string weightstr, std::string jetRname, 
             std::string thrname, bool include_RL0, bool include_RL1) {
    
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

                    
            // need to do something different for charge! pt cuts not appropriate
            TH1D * h_binbybin_fcorr = get_binbybin_corrfactors(f_in_mc, f_out, observable, i, pt_min, pt_max, k, RL_min, RL_max, true);
            TH1D * h_raw_data = get_rawdata(f_in_data, observable, pt_min, pt_max, RL_min, RL_max);
            apply_corrfactor(f_out, h_raw_data, h_binbybin_fcorr, observable, i, k);
        
        }
        


    }


}


//===========================================================================
// This file takes in the following arguments:
// [[nothing]]
// THIS IS GOOD FOR PERLMUTTER
//===========================================================================
void extract_binbybin_corrections() { 

    gStyle->SetOptStat(0); // hide stats panel
    SetStyle();

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


    // attempt_dir name
    std::string matched_str = "matched";
    if (unmatched) matched_str = "unmatched";
    if (anchmc) attempt_dir = Form("binbybincorrections/%s/rebinx%.0f", matched_str.c_str(), rebin);
    else attempt_dir = Form("fastsim_binbybincorrections/%s/rebinx%.0f", matched_str.c_str(), rebin);
    if (ptrl_bins) Form("binbybincorrections_ptrlbins");
    
    
    // filenames
    // file that needs correcting:
    // TString input_histograms_filename = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/rootfiles/data_secondattempt/rebinx4/DataHists.root"; //"/software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/analysis/user/dEEC/datahists/DataHists.root";
    TString input_histograms_filename = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/rootfiles/data_thirdattempt_ptrlbins/logbins/DataHists.root"; //"/software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/analysis/user/dEEC/datahists/DataHists.root";
    // if (compsystem == "local") input_histograms_filename = "/Volumes/WORK USB/dEEC/storage/rootfiles/data_secondattempt/rebinx4/DataHists.root";
    TFile* root_data_file = new TFile(input_histograms_filename, "READ");

    // file with anchored mc - truth vs det level information
    TString input_mc_filename = "";
    if (anchmc) {
        input_mc_filename = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/34547495/scaling/AnalysisResultsFinal.root"; //LHC23a3
        if (compsystem == "local") input_mc_filename = "/Volumes/WORK USB/dEEC/storage/slurmfiles/perly/34547495/AnalysisResultsFinal.root";
    } else {
        input_mc_filename = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/35235011/scaling/AnalysisResultsFinal.root"; //pythia fastsim
        if (compsystem == "local") input_mc_filename = ""; //pythia fastsim
    }
    TFile* root_mc_file = new TFile(input_mc_filename, "READ");
    
    TFile* root_mc_jtjl_file;
    if (compsystem == "local") {
        TString input_mc_jtjl_filename = "/Volumes/WORK USB/dEEC/storage/slurmfiles/perly/35235011/AnalysisResultsFinal.root";
        root_mc_jtjl_file = new TFile(input_mc_jtjl_filename, "READ");
    }
    

    // Output file with corrected results
    // std::string outfile = "/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/" + attempt_dir + "/DataHists.root";
    std::string outfile = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/rootfiles/" + attempt_dir + "/DataHists_BinByBinCorr.root";
    if (compsystem == "local") outfile = "/Volumes/WORK USB/dEEC/storage/rootfiles/" + attempt_dir + "/DataHists_BinByBinCorr.root";
    TFile* root_outfile = new TFile(outfile.c_str(), "RECREATE");

    analyze(root_data_file, root_mc_file, root_outfile, pt_bins, n_bins, RL_bins, n_RLbins, "deltap", weightstr, jetRname, thrname, include_RL0, include_RL1);
    // analyze(root_data_file, root_mc_file, root_outfile, pt_bins, n_bins, RL_bins, n_RLbins, "deltapt", weightstr, jetRname, thrname, include_RL0, include_RL1);
    // analyze(root_data_file, root_mc_file, root_outfile, pt_bins, n_bins, RL_bins, n_RLbins, "deltapl", weightstr, jetRname, thrname, include_RL0, include_RL1);
    
    // analyze_charge(root_data_file, root_mc_file, root_outfile, pt_bins, n_bins, RL_bins, n_RLbins, "charge", weightstr, jetRname, thrname, include_RL0, include_RL1);
    

    analyze(root_data_file, root_mc_jtjl_file, root_outfile, pt_bins, n_bins, RL_bins, n_RLbins, "deltajt", weightstr, jetRname, thrname, include_RL0, include_RL1);
    // analyze(root_data_file, root_mc_jtjl_file, root_outfile, pt_bins, n_bins, RL_bins, n_RLbins, "deltajl", weightstr, jetRname, thrname, include_RL0, include_RL1);
    


    root_data_file->Close();
    root_mc_file->Close();
    root_mc_jtjl_file->Close();
    root_outfile->Close();

    delete root_data_file;
    delete root_mc_file;
    delete root_mc_jtjl_file;
    delete root_outfile;


}
