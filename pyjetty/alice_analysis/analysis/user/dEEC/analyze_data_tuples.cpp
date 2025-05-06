// ROOT macro to analyze and plot data tuples 
// This version will only accomodate the FOUR RL/pTRL bins
// Beatrice Liang-Gilman (beatrice_lg@berkeley.edu)

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

// global variables
Double_t colors[16] = {kGray, kMagenta, kBlue, kOrange+1, kViolet+1, kGreen+2, kRed, kYellow+1, kCyan+1};
Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
Double_t marker_size = 1.5;

int rebin = 4;
bool ptrl_bins = true;
bool logbins = false;
std::string attempt_dir; // = Form("data_secondattempt/rebinx%d", rebin);
std::string outdir; // = "/software/users/blianggi/mypyjetty/storage/dEEC/plots/" + attempt_dir;

bool jetpt_bool = true;
bool deltap_bool = true;
bool deltapt_bool = false;
bool deltapl_bool = false;
bool deltajt_bool = true;
bool deltajl_bool = false;
bool ew_bool = true;
bool twoDhists_bool = true;
bool rc_bool = true;

bool unnormalized_bool = true;
bool self_normalized_bool = true;
bool norm_by_jets_bool = true;

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


/* Make a vector of the logarithmic bins given the minimum and maximum values of the bins. */
std::vector<double> makeLogBins(double minVal, double maxVal, int numBins, bool debug=false) {  
    
    std::vector<double> bins(numBins + 1);
    double logMin = log10(minVal);
    double logMax = log10(maxVal);
    double logStep = (logMax - logMin) / numBins;
    
    // if (debug) cout << "minval: " << minVal << " maxval: " << maxVal << " numbins: " << numBins << endl;
    // if (debug) cout << "logMin: " << logMin << " logMax: " << logMax << " logStep: " << logStep << endl;
    
    for (int i = 0; i <= numBins; i++) {
        bins[i] = std::pow(10, logMin + i * logStep);
        if (debug) std::cout << "Bin " << i << ": " << bins[i] << std::endl;
    }
    
    return bins;
}

/* get a typical 1D histogram from the TChain */
TH1D * getObs1DHistFromTChain(TChain *chain, std::string branch_name, int num_bins, double hist_xmin, double hist_xmax,
                              int pt_min, int pt_max, double RL_min, double RL_max, double ptavg=0.0, bool logbinning=false,
                              bool debug=true) {
    
    TH1D * hist1D;
    if (logbinning == false) {
        hist1D = new TH1D(Form("%s_hist", branch_name.c_str()), Form("%s_hist", branch_name.c_str()), num_bins, hist_xmin, hist_xmax);
    } else {
        double new_hist_xmin = (hist_xmin == 0.0) ? 1e-1 : hist_xmin;
        if (branch_name == "p1" || branch_name == "pt1" || branch_name == "jl1") new_hist_xmin = 0.8;
        std::vector<double> bins = makeLogBins(new_hist_xmin, hist_xmax, num_bins);
        hist1D = new TH1D(Form("%s_hist", branch_name.c_str()), Form("%s_hist", branch_name.c_str()), num_bins, bins.data());
    }

    if ( branch_name == "jet_pt" || branch_name == "total_num_const" || branch_name == "num_const_aftercut") {
        chain->Draw(Form("%s>>%s_hist", branch_name.c_str(), branch_name.c_str()), Form("jet_pt >= %d && jet_pt < %d", pt_min, pt_max), "e");
    } else {
        if (ptrl_bins == false) {
            chain->Draw(Form("%s>>%s_hist", branch_name.c_str(), branch_name.c_str()), Form("jet_pt >= %d && jet_pt < %d && RL >= %f && RL < %f", pt_min, pt_max, RL_min, RL_max), "e");
        } else {
            chain->Draw(Form("%s>>%s_hist", branch_name.c_str(), branch_name.c_str()), Form("jet_pt >= %d && jet_pt < %d && %f*RL >= %f && %f*RL < %f", pt_min, pt_max, ptavg, RL_min, ptavg, RL_max), "e");
        }
    }
//    chain->Draw(Form("%s>>%s_hist", branch_name.c_str(), branch_name.c_str()), Form("jet_pt >= 20 && jet_pt < 40 && RL > 0.01 && RL < 0.4"), "e");

    if (debug) {
        cout << "IN GETOBS1DHISTFROMTCHAIN" << endl;
        int nBins = hist1D->GetNbinsX();
        std::cout << "Histogram: " << hist1D->GetName() << std::endl;
        for (int i = 1; i <= nBins; i++) {  // Bins are 1-indexed in ROOT
            double binLowEdge = hist1D->GetBinLowEdge(i);
            double binUpEdge  = hist1D->GetBinLowEdge(i + 1);
            double binContent = hist1D->GetBinContent(i);
            
            std::cout << "Bin " << i << ": [" << binLowEdge << ", " << binUpEdge 
                    << "] Content: " << binContent << std::endl;
        }
    }

    return hist1D;
}

/* Get the r_c from TChain */
double getRcFromTChain(TChain *chain, std::string branch_name, int num_bins, double hist_xmin, double hist_xmax,
                              int pt_min, int pt_max, double RL_min, double RL_max, double ptavg) 
{
    // draw regular charge histogram, where like sign = +1, and unlike sign = -1
    TH1D *hist_charge = new TH1D("hist_charge", "hist_charge", num_bins, hist_xmin, hist_xmax);
    if (ptrl_bins == false) {
        chain->Draw("q1q2>>hist_charge", Form("jet_pt >= %d && jet_pt < %d && RL >= %f && RL < %f", pt_min, pt_max, RL_min, RL_max), "e");
    } else {
        chain->Draw("q1q2>>hist_charge", Form("jet_pt >= %d && jet_pt < %d && %f*RL >= %f && %f*RL < %f", pt_min, pt_max, ptavg, RL_min, ptavg, RL_max), "e");
    }
    // do i need to scale by the RL bin width here?? - I think this would be redundant.
    // if both like sign bin and unlike sign get scaled by RL bin width, then the ratio still stays the same

    // get # of like sign and # of unlike sign
    int num_likesign = hist_charge->GetBinContent(hist_charge->FindBin(1));
    int num_unlikesign = hist_charge->GetBinContent(hist_charge->FindBin(-1));
    // if (debug) cout << "num like sign " << num_likesign << " num unlike sign " << num_unlikesign << endl;

    // calculate the rc value for this pt & RL bin
    double rc = (double)(num_likesign - num_unlikesign) / (double)(num_likesign + num_unlikesign);
    //if (debug) cout << "and that makes rc " << rc << endl;

    return rc;
}

/* Get the r_c from TChain */
double getRcErr(TChain *chain, std::string branch_name, int num_bins, double hist_xmin, double hist_xmax,
                              int pt_min, int pt_max, double RL_min, double RL_max, double ptavg) 
{
    // draw regular charge histogram, where like sign = +1, and unlike sign = -1
    TH1D *hist_charge = new TH1D("hist_charge", "hist_charge", num_bins, hist_xmin, hist_xmax);
    if (ptrl_bins == false) {
        chain->Draw("q1q2>>hist_charge", Form("jet_pt >= %d && jet_pt < %d && RL >= %f && RL < %f", pt_min, pt_max, RL_min, RL_max), "e");
    } else {
        chain->Draw("q1q2>>hist_charge", Form("jet_pt >= %d && jet_pt < %d && %f*RL >= %f && %f*RL < %f", pt_min, pt_max, ptavg, RL_min, ptavg, RL_max), "e");
    }
    
    // do i need to scale by the RL bin width here?? - I think this would be redundant.
    // if both like sign bin and unlike sign get scaled by RL bin width, then the ratio still stays the same

    //---====---====---====---====---====---====---====---====---====---====---====
    //---====---====---====---====---====---====---====---====---====---====---====

    // get # of like sign and # of unlike sign
    double num_likesign = hist_charge->GetBinContent(hist_charge->FindBin(1));
    double num_unlikesign = hist_charge->GetBinContent(hist_charge->FindBin(-1));

    cout << "NUM LIKE SIGN PAIRS IS " << num_likesign << " and NUM DISLIKE" << num_unlikesign << endl;

    // calculate the rc error for this pt & RL bin
    double num_totalpairs = num_likesign + num_unlikesign;
    double rc_err = ( 2 * sqrt( num_totalpairs * num_likesign * num_unlikesign ) ) / (num_totalpairs * num_totalpairs);

    cout << "RC ERR IN FUNC IS " << rc_err << endl;
        
    return rc_err;
}

/* get a typical 2D histogram from the TChain */
TH2D * getObs2DHistFromTChain(TChain *chain, std::string branch_name_x, std::string branch_name_y,
                              int num_bins_x, double hist_xmin, double hist_xmax,
                              int num_bins_y, double hist_ymin, double hist_ymax,
                              int pt_min, int pt_max, double RL_min, double RL_max, double ptavg) {

    TH2D * hist2D = new TH2D(Form("%s_vs_%s_hist", branch_name_y.c_str(), branch_name_x.c_str()), Form("%s_vs_%s_hist", branch_name_y.c_str(), branch_name_x.c_str()), num_bins_x, hist_xmin, hist_xmax, num_bins_y, hist_ymin, hist_ymax);
    if (branch_name_x == "zi") { // zi vs zj
        if (ptrl_bins == false) {
            chain->Draw(Form("pt2/jet_pt:pt1/jet_pt>>%s_vs_%s_hist", branch_name_y.c_str(), branch_name_x.c_str()), Form("jet_pt >= %d && jet_pt < %d && RL >= %f && RL < %f", pt_min, pt_max, RL_min, RL_max)); //, "colz");
        } else {
            chain->Draw(Form("pt2/jet_pt:pt1/jet_pt>>%s_vs_%s_hist", branch_name_y.c_str(), branch_name_x.c_str()), Form("jet_pt >= %d && jet_pt < %d && %f*RL >= %f && %f*RL < %f", pt_min, pt_max, ptavg, RL_min, ptavg, RL_max));
        }
    } else if (branch_name_x == "zsmol") { //zbig vs zsmol
        cout << "this isn't implemented yet! (idk how to )" << endl;
        return hist2D;
    } else if (branch_name_x == "maxpt") { //weights vs maxpt
        cout << "this isn't implemented yet! (idk how to )" << endl;
        return hist2D;
    } else {
        if (ptrl_bins == false) {
            chain->Draw(Form("%s:%s>>%s_vs_%s_hist", branch_name_y.c_str(), branch_name_x.c_str(), branch_name_y.c_str(), branch_name_x.c_str()), Form("jet_pt >= %d && jet_pt < %d && RL >= %f && RL < %f", pt_min, pt_max, RL_min, RL_max)); //, "colz");
        } else {
            chain->Draw(Form("%s:%s>>%s_vs_%s_hist", branch_name_y.c_str(), branch_name_x.c_str(), branch_name_y.c_str(), branch_name_x.c_str()), Form("jet_pt >= %d && jet_pt < %d && %f*RL >= %f && %f*RL < %f", pt_min, pt_max, ptavg, RL_min, ptavg, RL_max));
        }
    }
    return hist2D;
}


/* Format and adjust histograms */
void Format1DHist(TH1D *hist, TH1D *jetpt_hist, std::string norm_string, int markercolor, double markeralpha,
                  int markerstyle, std::string xtitle, std::string ytitle, TLegend& leg, TString leg_text, 
                  bool scalebyRLbinwidth, double RL_bin_width_val, std::string obs_name="", std::string hist_addname="") {

    hist->SetTitle(Form("h_%s%s", obs_name.c_str(), hist_addname.c_str()));
    hist->SetName(Form("h_%s%s", obs_name.c_str(), hist_addname.c_str()));

    // rebin before scaling!!!
    if (rebin != 0 && obs_name != "weights") hist->Rebin(rebin);

    // normalization
    if ( norm_string == "self_normalized" ) {
        double selfnorm_value = hist->Integral();
        hist->Scale(1/selfnorm_value, "width");
    } else if ( norm_string == "norm_by_jets" ) {
        double numjets = jetpt_hist->Integral();
        cout << "Number of jets in " << leg_text << ": " << numjets << endl;
        hist->Scale(1/numjets, "width");
    }

    // cout << "RL BIN WIDTH! " << RL_bin_width_val << endl;
    // if (scalebyRLbinwidth) hist->Scale(1/RL_bin_width_val);

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
                  bool scalebyRLbinwidth, double RL_bin_width, std::string obs_name_x, std::string obs_name_y, std::string hist_addname="",
                  bool restrictzaxis=false) {

    double ptvsew_normbounds[3][2] = { { 0, 500 }, { 1e-2, 5e2 }, { 1e-4, 2e2 }}; //TODO: make this less pt specific?
    // also TODO: potentially set all lower boundaries to 0 for data??

    hist2D->SetTitle(Form("h_%s_vs_%s%s", obs_name_y.c_str(), obs_name_x.c_str(), hist_addname.c_str()));
    hist2D->SetName(Form("h_%s_vs_%s%s", obs_name_y.c_str(), obs_name_x.c_str(), hist_addname.c_str()));

    // rebin before scaling!!!
    // if (rebin != 0) hist2D->RebinX(rebin); //TODO: put this back in later??

    // normalization
    int norm_index = 0;
    if ( norm_string == "self_normalized" ) {
        double selfnorm_value = hist2D->Integral();
        hist2D->Scale(1/selfnorm_value, "width");
        norm_index = 1;
    } else if ( norm_string == "norm_by_jets" ) {
        double numjets = jetpt_hist->Integral();
        hist2D->Scale(1/numjets, "width");
        norm_index = 2;
    }
    // if ( scalebyRLbinwidth ) hist2D->Scale(1/RL_bin_width);

    // set z axis bounds
    if (restrictzaxis) hist2D->GetZaxis()->SetRangeUser(ptvsew_normbounds[norm_index][0], ptvsew_normbounds[norm_index][1]);

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


TGraphErrors * MakeFormatGraph(vector<double> xvals, vector<double> yvals, int markercolor, double markeralpha,
                  int markerstyle, std::string xtitle, std::string ytitle, std::string obs_name, std::string hist_addname="") {

    TGraphErrors * graph = new TGraphErrors(xvals.size(), xvals.data(), yvals.data());
    graph->SetTitle(Form("Charge Ratio;%s;%s", xtitle.c_str(), ytitle.c_str())); // Set the title and axis labels
    graph->SetName(Form("rc_%s", hist_addname.c_str()));

    // Set graph styles
    graph->SetLineColorAlpha(markercolor, markeralpha);
    graph->SetMarkerColorAlpha(markercolor, markeralpha);
    graph->SetMarkerStyle(markerstyle);
    graph->SetMarkerSize(1.5);

    // axes
    graph->GetXaxis()->SetLabelFont(42);
    graph->GetXaxis()->SetTitleFont(42);
	graph->GetXaxis()->SetTitleSize(0.06); //(0.042);
    graph->GetXaxis()->SetTitleOffset(1.0);
	graph->GetXaxis()->SetLabelSize(0.05);
    // graph->GetXaxis()->SetTitle(xtitle.c_str());

    graph->GetYaxis()->SetLabelFont(42);
	graph->GetYaxis()->SetTitleFont(42);
    graph->GetYaxis()->SetTitleOffset(1.05); 
	graph->GetYaxis()->SetTitleSize(0.06); //(0.042);
	graph->GetYaxis()->SetLabelSize(0.05); //(0.042);
    // graph->GetYaxis()->SetTitle(ytitle.c_str());


    return graph;
}


/* Save and delete histograms */
// (TFile *fout, TCanvas *can, TH1D* hist, std::string obs_name, 
//                          std::string ptname, std::string norm_string, std::string hist_addname,
//                          bool logx, bool logy)
void draw_save_del_hists(TFile *fout, TCanvas *can, TObject* obj, std::string obs_name, 
                         std::string ptname, std::string norm_string, std::string hist_addname,
                         bool logx, bool logy, bool logz=false, std::string obs_filename="") {
    
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
    std::string usethisname = (obs_filename == "") ? obs_name : obs_filename;

    std::string fname_out = outdir + add_dir + "/corrhist_" + usethisname + hist_addname + ".pdf";
    can->SaveAs(fname_out.c_str());

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
                          int pt_max, bool logx, bool logy, double pl_axis_cut=-1, 
                          std::string obs_filename="", bool debug=false) {

    // go into canvas
    can_all->cd();
    if (logx) gPad->SetLogx();
    if (logy) gPad->SetLogy();

    // if momentum axis, adjust x bounds accordingly
    size_t length = h_vec.size();
    for (int j=0; j<length; j++) {
        // cout << j << ": " << RL_bin_width[j] << endl;
        // if (scalebyRLbinwidth) h_vec[j]->Scale(RL_bin_width[j]); // this needs to be done before normalization
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
    
    if (norm_string == "self_normalized") {
        if (obs_filename == "") {
            if (obs_name == "deltap") {
                h_vec[0]->SetMinimum( 1e-5 ); h_vec[0]->SetMaximum( 1.0 );
            } else if (obs_name == "deltajt") {
                h_vec[0]->SetMinimum( 1e-4 ); h_vec[0]->SetMaximum( 100. );
            }
        } else if (obs_filename == "p") {
            h_vec[0]->SetMinimum( 1e-5 ); h_vec[0]->SetMaximum( 1. );  
        } else if (obs_filename == "jt") {
            h_vec[0]->SetMinimum( 1e-4 ); h_vec[0]->SetMaximum( 70.0 );  
        }
    } else {
        double max = 0;
        for (int j=0; j<length; j++) {
            double max_cand = h_vec[j]->GetMaximum();
            if (max_cand > max) max = max_cand;
        }
        if (debug) cout << "max is " << max << " which goes to " << max*1.5 << endl;
        h_vec[0]->SetMaximum( max * 1.5 );
    }

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
    std::string usethisname = (obs_filename == "") ? obs_name : obs_filename;

    std::string fname_out = outdir + add_dir + "/corrhist_" + usethisname + "_ALL" + hist_addname + ".pdf";
    can_all->SaveAs(fname_out.c_str());

    deleteVecOfHists(h_vec);
    delete can_all;

}

// ======================================================= //
//                   SPECIFIC FUNCTIONS
// ======================================================= //

void plot_rc(vector<vector<double>>& RL_vals, vector<vector<double>>& rc_vals, vector<double>& ptcenter_bins,
             vector<vector<double>> rc_errors) {
            //  TLegend& leg_RLbins, TLegend& leg_ptbins) {
    
    vector<TGraphErrors *> rc_graphs_func_of_RL;
    vector<TGraphErrors *> rc_graphs_func_of_pT;
    vector<vector<TGraphErrors*>> rc_graphs_func_of_RL_ind;
    vector<vector<TGraphErrors*>> rc_graphs_func_of_pT_ind;

    vector<double> RL_err;
    vector<double> pt_err;
    for (int i=0; i<ptcenter_bins.size(); i++) pt_err.push_back(0);
    for (int j=0; j<RL_vals[0].size(); j++) RL_err.push_back(0);

    cout <<" RL_err size " << RL_err.size() << endl;
    cout <<" pt_err size " << pt_err.size() << endl;

    // std::string outdir = "/software/users/blianggi/mypyjetty/storage/dEEC/plots/data_firstattempt"; // + ptbin_name + "/";//"plots/test/";
    std::string fname_func_of_RL_out = outdir + "/corrhist_rc_func_of_RL.pdf"; // could add jetR and threshold info later??, maybe not needed tho 
    std::string fname_func_of_pT_out = outdir + "/corrhist_rc_func_of_pT.pdf";

    TLegend leg_RLbins(0.2, 0.6, 0.4, 0.85); 
    TLegend leg_ptbins(0.5, 0.7, 0.65, 0.85); 

    leg_RLbins.SetTextSize(0.037);
    leg_RLbins.SetBorderSize(0);
    leg_ptbins.SetTextSize(0.037);
    leg_ptbins.SetBorderSize(0);

    // delete hist;

    // get graphs of r_c as a function of RL
    // loop over pt bins
    for ( int i = 0; i < ptcenter_bins.size(); i++ ) { 

        // for graphs as a function of RL
        TGraphErrors *g = new TGraphErrors(RL_vals[i].size(), RL_vals[i].data(), rc_vals[i].data(), RL_err.data(), rc_errors[i].data());
        g->SetMarkerStyle(markers[i]);
        g->SetMarkerSize(1.5);
        g->SetMarkerColorAlpha(kBlack, 1.0);
        g->SetLineColorAlpha(kBlack, 1.0);
        vector<TGraphErrors*> ind_temp_vec;

        for (int j=0; j<RL_vals[i].size(); j++){
            int k=j+1;
            TGraphErrors *g_ind = new TGraphErrors(1, &RL_vals[i][j], &rc_vals[i][j], &RL_err[j], &rc_errors[i][j]);
            cout << "i: " << i << " j: " << j << " err: " << rc_errors[i][j] << endl;
            g_ind->SetMarkerColorAlpha(colors[k], 1.0);
            g_ind->SetLineColorAlpha(colors[k], 1.0);
            g_ind->SetMarkerSize(1.5);
            g_ind->SetMarkerStyle(markers[i]);
            ind_temp_vec.push_back(g_ind);
            
            if (i==0) {
                leg_RLbins.AddEntry(g_ind, Form("R_{L} bin %d", j+1), "P");
            }
        }
        rc_graphs_func_of_RL.push_back(g); 
        rc_graphs_func_of_RL_ind.push_back(ind_temp_vec); 
        
        // Should fix what gets subbed into %d so it is more flexible if the bins are not 20 GeV big?
        leg_ptbins.AddEntry(g, Form("p_{T} = %d-%d", (int)ptcenter_bins[i]-10, (int)ptcenter_bins[i]+10), "P");
         

    }

    // plot r_c as a function of RL
    TCanvas *can_func_of_RL = new TCanvas("can_func_of_RL", "can_func_of_RL", 750, 500);
    can_func_of_RL->cd();
    for ( int i = 0; i < ptcenter_bins.size(); i++ ) {
        // rc_graphs_func_of_RL[i]->SetMarkerSize(1.0);
        // rc_graphs_func_of_RL[i]->SetMarkerStyle(markers[i]);
        if ( i == 0 ) {
            rc_graphs_func_of_RL[i]->SetMinimum(-0.3);  // Lower y limit
            rc_graphs_func_of_RL[i]->SetMaximum(0.1); 
            rc_graphs_func_of_RL[i]->GetXaxis()->SetTitle("R_{L} bin center"); 
            rc_graphs_func_of_RL[i]->GetYaxis()->SetTitle("r_{c}"); 
            rc_graphs_func_of_RL[i]->Draw("AP");
        } // else {
        //     rc_graphs_func_of_RL[i]->Draw("P SAME");
        // }

        for (int j = 0; j < RL_vals[0].size(); j++) {
            rc_graphs_func_of_RL_ind[i][j]->Draw("P SAME");
        }
    }
    leg_RLbins.Draw("same");
    leg_ptbins.Draw("same");
    can_func_of_RL->SaveAs(fname_func_of_RL_out.c_str());
    delete can_func_of_RL;

    //========================================================

    // get graphs of r_c as a function of pT
    vector<vector<double>> rc_vals_func_of_pT;
    vector<vector<double>> rc_err_vals_func_of_pT;
    // loop over RL bins
    for ( int j = 0; j < RL_vals[0].size(); j++ ) {
        vector<double> temp_vec;
        vector<TGraphErrors*> ind_temp_vec;
        vector<double> err_temp_vec;

        // save values into appropriate vectors
        for ( int i = 0; i < ptcenter_bins.size(); i++ ) { 
            temp_vec.push_back(rc_vals[i][j]);
            err_temp_vec.push_back(rc_errors[i][j]);
            
            int k=j+1;
            TGraphErrors *g_ind = new TGraphErrors(1, &ptcenter_bins[i], &rc_vals[i][j], &pt_err[i], &rc_errors[i][j]);
            cout << "i: " << i << " j: " << j << " err: " << rc_errors[i][j] << endl;
            g_ind->SetMarkerColorAlpha(colors[k], 1.0);
            g_ind->SetLineColorAlpha(colors[k], 1.0);
            g_ind->SetMarkerSize(1.5);
            g_ind->SetMarkerStyle(markers[i]);
            ind_temp_vec.push_back(g_ind);
        }
        rc_vals_func_of_pT.push_back(temp_vec);
        rc_err_vals_func_of_pT.push_back(err_temp_vec);
        rc_graphs_func_of_pT_ind.push_back(ind_temp_vec); 

        // for graphs as a function of RL
        TGraphErrors *g = new TGraphErrors(ptcenter_bins.size(), ptcenter_bins.data(), rc_vals_func_of_pT[j].data(), pt_err.data(), rc_err_vals_func_of_pT[j].data());
        for (int aa = 0; aa < ptcenter_bins.size(); aa++) {
            // cout << "studying pt=" << ptcenter_bins[aa] << " // " << rc_vals_func_of_pT[j][aa] << endl;
        }
        rc_graphs_func_of_pT.push_back(g); 
    }
    
    // plot r_c as a function of pT
    TCanvas *can_func_of_pT = new TCanvas("can_func_of_pT", "can_func_of_pT", 750, 500);
    can_func_of_pT->cd();
    for ( int j = 0; j < RL_vals[0].size(); j++ ) {
        int k = j+1;
        // rc_graphs_func_of_pT[j]->SetMarkerSize(1.0);
        // rc_graphs_func_of_pT[j]->SetMarkerStyle(markers[0]);
        rc_graphs_func_of_pT[j]->SetMarkerColorAlpha(colors[k], 0.0);
        if ( j == 0 ) {
            rc_graphs_func_of_pT[j]->SetMinimum(-0.3);  // Lower y limit
            rc_graphs_func_of_pT[j]->SetMaximum(0.1); 
            rc_graphs_func_of_pT[j]->GetXaxis()->SetTitle("p_{T} bin center"); 
            rc_graphs_func_of_pT[j]->GetYaxis()->SetTitle("r_{c}"); 
            rc_graphs_func_of_pT[j]->Draw("AP");
        } // else { 
        //     rc_graphs_func_of_pT[j]->Draw("P SAME");
        // }

        for (int i = 0; i < ptcenter_bins.size(); i++) {
            rc_graphs_func_of_pT_ind[j][i]->Draw("P SAME");
        }
    }
    leg_RLbins.Draw("same");
    leg_ptbins.Draw("same");
    can_func_of_pT->SaveAs(fname_func_of_pT_out.c_str());
    delete can_func_of_pT;



}

// ======================================================= //
//                     SOME FUNCTIONS
// ======================================================= //

void analyze_ptbin(TChain * JETINFO_tree, TChain * PAIRINFO_tree, 
             TFile *f_out, std::string weightstr, std::string jetRname, std::string thrname,
             std::string norm_string, int pt_min, int pt_max, double pt_avg,  const double RL_bins[], int n_RLbins, 
             vector<vector<double>>& RL_vals, vector<vector<double>>& rc_vals, vector<vector<double>>& rc_errors,
             bool include_RL0, bool include_RL1, bool debug, bool debug2) {
    
    std::string ptname = to_string(pt_min) + "-" + to_string(pt_max);
    
    double RL_bin_width[8] = {0}; 
    double RL_bin_centers[8] = {0};
    for (int j = 0; j < n_RLbins; ++j) {
        RL_bin_width[j] = RL_bins[j+1] - RL_bins[j];
        RL_bin_centers[j] = (RL_bins[j+1] + RL_bins[j])/2;
        // cout << "RL BIN WIDTH HERE" << RL_bin_width[i][j] << endl;
        // cout << " AND CENTERS " << RL_bin_centers[j] << endl;
    }

    std::string ytitle_norm = ""; //"#frac{1}{#DeltaR_{L}} ";
    if (ptrl_bins == false) {
        if (norm_string == "self_normalized") ytitle_norm = "#frac{1}{N_{pair}} "; //#DeltaR_{L}} ";
        else if (norm_string == "norm_by_jets") ytitle_norm = "#frac{1}{N_{jet}} "; //#DeltaR_{L}} ";
    } else {
        if (norm_string == "unnormalized") ytitle_norm = ""; //"#frac{1}{#Delta(#LTp_{T}#GTR_{L})} ";
        else if (norm_string == "self_normalized") ytitle_norm = "#frac{1}{N_{pair}} "; //#Delta(#LTp_{T}#GTR_{L})} ";
        else if (norm_string == "norm_by_jets") ytitle_norm = "#frac{1}{N_{jet}} "; //#Delta(#LTp_{T}#GTR_{L})} ";
    }
    
    
    vector<TH1D*> deltap_vec;
    vector<TH1D*> deltapt_vec;
    vector<TH1D*> deltapl_vec;
    vector<TH1D*> weights_vec;
    // vector<TH1D*> q1q2_vec;
    vector<double> rc_vec;
    vector<double> rc_err_vec;
    vector<double> RLcenters_vec;

    vector<TH1D*> deltajt_vec;
    vector<TH1D*> deltajl_vec;

    vector<TH1D*> p_vec;
    vector<TH1D*> pt_vec;
    vector<TH1D*> pl_vec;
    vector<TH1D*> jt_vec;
    vector<TH1D*> jl_vec;

    // TLegend *leg = new TLegend(0.2, 0.16, 0.4, 0.39); //(0.6, 0.6, 0.85, 0.87);
    TLegend *leg = new TLegend(0.5, 0.62, 0.85, 0.85);
    leg->SetTextSize(0.037);
    leg->SetBorderSize(0);
    // leg->AddEntry(NULL, Form("%d #leq p_{T, jet} < %d", pt_min, pt_max)); //, "pl");
    cout << "ADDING TO LEGEND HERE!!!!" << pt_min << " " << pt_max << endl;
    leg->AddEntry("NULL",Form("%d #leq p_{T, jet} < %d", pt_min, pt_max),"h");
    TLegend *leg_dummy = new TLegend();
    
    for ( int j = 0; j < n_RLbins; j++ ) {
        int k = j;
        if (!include_RL0) {
            k = j-1;
            if (j == 0) continue; // can add something here to change the filename for ALL
        }
        if (!include_RL1 && j == n_RLbins-1) continue;
        
        double RL_min = RL_bins[j];
        double RL_max = RL_bins[j+1];
        std::string RLname = Form("_RL%.3f-%.3f", RL_min, RL_max);
        std::string RLname_leg = Form("R_{L} = %.3f-%.3f", RL_min, RL_max);
        if (ptrl_bins == true) {
            RLname = Form("_pTRL%.3f-%.3f", RL_min, RL_max);
            RLname_leg = Form("#LT p_{T} #GT R_{L} = %.3f-%.3f", RL_min, RL_max);
        }
        std::string hist_addname = weightstr + jetRname + thrname + "_pt" + ptname + RLname + "_" + norm_string;
        if (debug) cout << " in RL bin" << j << " with " << RL_min << " - " << RL_max << endl;
        
        // bin sizes
        double deltap_binsize = 0.625; // for 2.5 size bins later!(=2.5/4) //0.5;
        // int deltap_numbins = int((pt_max+5)/deltap_binsize);
        int deltap_numbins = int((85)/deltap_binsize);
        double deltapl_binsize = 0.5;
        int deltapl_numbins = int((pt_max/2)/deltap_binsize);
        int deltajt_numbins = 200; 
        int weights_numbins = int(0.3/0.005);
        if (logbins) {
            deltap_numbins /= 2; //45 bins for 20-40 jets, 85 bins for 60-80 jets
            deltajt_numbins = 50;
        }

        // get histograms
        TH1D * deltap_hist, * deltapt_hist, * deltapl_hist, * weights_hist;
        TH1D * deltajt_hist, * deltajl_hist;
        TH1D * p_hist, * pt_hist, * pl_hist, * jt_hist, * jl_hist;
        
        TH1D * jetpt_inptbin_hist = getObs1DHistFromTChain(JETINFO_tree, "jet_pt", 100, 0, 200, pt_min, pt_max, 0, 0, pt_avg);
        
        // if (deltap_bool) deltap_hist = getObs1DHistFromTChain(PAIRINFO_tree, "deltap", deltap_numbins, 0, pt_max+5, pt_min, pt_max, RL_min, RL_max, pt_avg, logbins);
        if (deltap_bool) deltap_hist = getObs1DHistFromTChain(PAIRINFO_tree, "deltap", deltap_numbins, 0, 85, pt_min, pt_max, RL_min, RL_max, pt_avg, logbins);
        cout << "checkpoint 1 " << deltap_hist->GetEntries() << endl;
        if (deltapt_bool) deltapt_hist = getObs1DHistFromTChain(PAIRINFO_tree, "deltapt", deltap_numbins, 0, pt_max+5, pt_min, pt_max, RL_min, RL_max, pt_avg, logbins);
        if (deltapl_bool) deltapl_hist = getObs1DHistFromTChain(PAIRINFO_tree, "deltapl", deltapl_numbins, 0, pt_max/2, pt_min, pt_max, RL_min, RL_max, pt_avg, logbins);
        if (ew_bool) weights_hist = getObs1DHistFromTChain(PAIRINFO_tree, "weights", weights_numbins, 0, 0.3, pt_min, pt_max, RL_min, RL_max, pt_avg, logbins);
        
        if (deltajt_bool) deltajt_hist = getObs1DHistFromTChain(PAIRINFO_tree, "deltajt", deltajt_numbins, 0, 5, pt_min, pt_max, RL_min, RL_max, pt_avg, logbins);
        if (deltajl_bool) deltajl_hist = getObs1DHistFromTChain(PAIRINFO_tree, "deltajl", deltap_numbins, 0, pt_max+5, pt_min, pt_max, RL_min, RL_max, pt_avg, logbins);
        
        // if (deltap_bool) p_hist = getObs1DHistFromTChain(PAIRINFO_tree, "p1", deltap_numbins, 0, pt_max+5, pt_min, pt_max, RL_min, RL_max, pt_avg, logbins);
        if (deltap_bool) p_hist = getObs1DHistFromTChain(PAIRINFO_tree, "p1", deltap_numbins, 0, 85, pt_min, pt_max, RL_min, RL_max, pt_avg, logbins);
        if (deltapt_bool) pt_hist = getObs1DHistFromTChain(PAIRINFO_tree, "pt1", deltap_numbins, 0, pt_max+5, pt_min, pt_max, RL_min, RL_max, pt_avg, logbins);
        if (deltapl_bool) pl_hist = getObs1DHistFromTChain(PAIRINFO_tree, "pl1", 200, -5, 5, pt_min, pt_max, RL_min, RL_max, pt_avg, logbins);
        if (deltajt_bool) jt_hist = getObs1DHistFromTChain(PAIRINFO_tree, "jt1", deltajt_numbins, 0, 5, pt_min, pt_max, RL_min, RL_max, pt_avg, logbins);
        if (deltajl_bool) jl_hist = getObs1DHistFromTChain(PAIRINFO_tree, "jl1", deltap_numbins, 0, pt_max+5, pt_min, pt_max, RL_min, RL_max, pt_avg, logbins);
        
        // cout << "BACK IN ANALYZEPTBIN" << endl;
        // int nBins = deltap_hist->GetNbinsX();
        // std::cout << "Histogram: " << deltap_hist->GetName() << std::endl;
        // for (int i = 1; i <= nBins; i++) {  // Bins are 1-indexed in ROOT
        //     double binLowEdge = deltap_hist->GetBinLowEdge(i);
        //     double binUpEdge  = deltap_hist->GetBinLowEdge(i + 1);
        //     double binContent = deltap_hist->GetBinContent(i);
            
        //     std::cout << "Bin " << i << ": [" << binLowEdge << ", " << binUpEdge 
        //             << "] Content: " << binContent << std::endl;
        // }

        double rc_value = 0.0;
        double rc_err = 0.0;
        if (norm_string == "unnormalized" && rc_bool) {
            TH1D * q1q2_hist = getObs1DHistFromTChain(PAIRINFO_tree, "q1q2", 6, -3, 3, pt_min, pt_max, RL_min, RL_max, pt_avg);
            rc_value = getRcFromTChain(PAIRINFO_tree, "rc", 6, -3, 3, pt_min, pt_max, RL_min, RL_max, pt_avg);
            rc_err = getRcErr(PAIRINFO_tree, "rc", 6, -3, 3, pt_min, pt_max, RL_min, RL_max, pt_avg);
            cout << "RC ERR IS " << rc_err << "(pt_min=" << pt_min << ", j=" << j << ")" <<endl;
        }

        int nbins_2D = 50;
        if (pt_min == 40 || pt_min == 60) nbins_2D = 30;
        TH2D * weights_vs_deltapt_hist2D, * weights_vs_deltajt_hist2D, * zj_vs_zi_hist2D;

        if (twoDhists_bool) {
            weights_vs_deltapt_hist2D = getObs2DHistFromTChain(PAIRINFO_tree, "deltapt", "weights", nbins_2D, 0, pt_max+5, nbins_2D, 0, 0.3, pt_min, pt_max, RL_min, RL_max, pt_avg);
            weights_vs_deltajt_hist2D = getObs2DHistFromTChain(PAIRINFO_tree, "deltajt", "weights", 100, 0, 5, nbins_2D, 0, 0.3, pt_min, pt_max, RL_min, RL_max, pt_avg);
            zj_vs_zi_hist2D = getObs2DHistFromTChain(PAIRINFO_tree, "zi", "zj", 50, 0, 1, 50, 0, 1, pt_min, pt_max, RL_min, RL_max, pt_avg);
            // // TH2D * zbig_vs_zsmol_hist2D = getObs2DHistFromTChain(PAIRINFO_tree, "zsmol", "zbig", 50, 0, 1, 50, 0, 1, pt_min, pt_max, RL_min, RL_max, pt_avg);
            // // TH2D * weights_vs_maxpt_hist2D = getObs2DHistFromTChain(PAIRINFO_tree, "maxpt", "weights", nbins_2D, 0, pt_max+5, nbins_2D, 0, 0.3, pt_min, pt_max, RL_min, RL_max, pt_avg);
        }

        // push to vectors
        if (deltap_bool) deltap_vec.push_back((TH1D*) deltap_hist->Clone(deltap_hist->GetName()));
        if (deltapt_bool) deltapt_vec.push_back((TH1D*) deltapt_hist->Clone(deltapt_hist->GetName()));
        if (deltapl_bool) deltapl_vec.push_back((TH1D*) deltapl_hist->Clone(deltapl_hist->GetName()));
        if (ew_bool) weights_vec.push_back((TH1D*) weights_hist->Clone(weights_hist->GetName()));

        if (deltajt_bool) deltajt_vec.push_back((TH1D*) deltajt_hist->Clone(deltajt_hist->GetName()));
        if (deltajl_bool) deltajl_vec.push_back((TH1D*) deltajl_hist->Clone(deltajl_hist->GetName()));

        if (deltap_bool) p_vec.push_back((TH1D*) p_hist->Clone(p_hist->GetName()));
        if (deltapt_bool) pt_vec.push_back((TH1D*) pt_hist->Clone(pt_hist->GetName()));
        if (deltapl_bool) pl_vec.push_back((TH1D*) pl_hist->Clone(pl_hist->GetName()));
        if (deltajt_bool) jt_vec.push_back((TH1D*) jt_hist->Clone(jt_hist->GetName()));
        if (deltajl_bool) jl_vec.push_back((TH1D*) jl_hist->Clone(jl_hist->GetName()));
        
        if (norm_string == "unnormalized" && rc_bool) {
            rc_vec.push_back(rc_value);
            rc_err_vec.push_back(rc_err);
            RLcenters_vec.push_back( (RL_min+RL_max)/2 );
        }

        // format histograms in vector
        if (deltap_bool) Format1DHist(deltap_vec[k], jetpt_inptbin_hist, norm_string, colors[j], 0.6, markers[0], "#Deltap", ytitle_norm + "#frac{dN}{d#Deltap}", *leg, RLname_leg, true, RL_bin_width[j], "deltap", hist_addname);
        if (deltapt_bool) Format1DHist(deltapt_vec[k], jetpt_inptbin_hist, norm_string, colors[j], 0.6, markers[0], "#Deltap_{T}", ytitle_norm + "#frac{dN}{d#Deltap_{T}}", *leg_dummy, RLname_leg, true, RL_bin_width[j], "deltapt", hist_addname);
        if (deltapl_bool) Format1DHist(deltapl_vec[k], jetpt_inptbin_hist, norm_string, colors[j], 0.6, markers[0], "#Deltap_{L}", ytitle_norm + "#frac{dN}{d#Deltap_{L}}", *leg_dummy, RLname_leg, true, RL_bin_width[j], "deltapl", hist_addname);
        if (ew_bool) Format1DHist(weights_vec[k], jetpt_inptbin_hist, norm_string, colors[j], 0.6, markers[0], "#frac{p_{T,1}p_{T,2}}{p_{T,jet}^{2}}", ytitle_norm + "#frac{dN}{d[EW]}", *leg_dummy, RLname_leg, true, RL_bin_width[j], "weights", hist_addname);
        
        if (deltajt_bool) Format1DHist(deltajt_vec[k], jetpt_inptbin_hist, norm_string, colors[j], 0.6, markers[0], "#Deltaj_{T}", ytitle_norm + "#frac{dN}{d#Deltaj_{T}}", *leg_dummy, RLname_leg, true, RL_bin_width[j], "deltajt", hist_addname);
        if (deltajl_bool) Format1DHist(deltajl_vec[k], jetpt_inptbin_hist, norm_string, colors[j], 0.6, markers[0], "#Deltaj_{L}", ytitle_norm + "#frac{dN}{d#Deltaj_{L}}", *leg_dummy, RLname_leg, true, RL_bin_width[j], "deltajl", hist_addname);
        
        if (deltap_bool) Format1DHist(p_vec[k], jetpt_inptbin_hist, norm_string, colors[j], 0.6, markers[0], "p", ytitle_norm + "#frac{dN}{dp}", *leg_dummy, RLname_leg, true, RL_bin_width[j], "p", hist_addname);
        if (deltapt_bool) Format1DHist(pt_vec[k], jetpt_inptbin_hist, norm_string, colors[j], 0.6, markers[0], "p_{T}", ytitle_norm + "#frac{dN}{dp_{T}}", *leg_dummy, RLname_leg, true, RL_bin_width[j], "pt", hist_addname);
        if (deltapl_bool) Format1DHist(pl_vec[k], jetpt_inptbin_hist, norm_string, colors[j], 0.6, markers[0], "p_{L}", ytitle_norm + "#frac{dN}{dp_{L}}", *leg_dummy, RLname_leg, true, RL_bin_width[j], "pl", hist_addname);
        if (deltajt_bool) Format1DHist(jt_vec[k], jetpt_inptbin_hist, norm_string, colors[j], 0.6, markers[0], "j_{T}", ytitle_norm + "#frac{dN}{dj_{T}}", *leg_dummy, RLname_leg, true, RL_bin_width[j], "jt", hist_addname);
        if (deltajl_bool) Format1DHist(jl_vec[k], jetpt_inptbin_hist, norm_string, colors[j], 0.6, markers[0], "j_{L}", ytitle_norm + "#frac{dN}{dj_{L}}", *leg_dummy, RLname_leg, true, RL_bin_width[j], "jl", hist_addname);
        
        if (twoDhists_bool) Format2DHist(weights_vs_deltapt_hist2D, jetpt_inptbin_hist, norm_string, ytitle_norm + "#Deltap_{T}", ytitle_norm + "p_{T,1}p_{T,2} / p_{T,jet}^{2}", true, RL_bin_width[j], "deltapt", "weights", hist_addname, true);
        if (twoDhists_bool) Format2DHist(weights_vs_deltajt_hist2D, jetpt_inptbin_hist, norm_string, ytitle_norm + "#Deltaj_{T}", ytitle_norm + "p_{T,1}p_{T,2} / p_{T,jet}^{2}", true, RL_bin_width[j], "deltajt", "weights", hist_addname);
        if (twoDhists_bool) Format2DHist(zj_vs_zi_hist2D, jetpt_inptbin_hist, norm_string, ytitle_norm + "z_{i}", ytitle_norm + "z_{j}", true, RL_bin_width[j], "zi", "zj", hist_addname);

        // draw, save, and delete histograms
        TCanvas *can_deltap = new TCanvas();
        TCanvas *can_deltapt = new TCanvas();
        TCanvas *can_deltapl = new TCanvas();
        TCanvas *can_weights = new TCanvas();
        
        TCanvas *can_deltajt = new TCanvas();
        TCanvas *can_deltajl = new TCanvas();

        TCanvas *can_p = new TCanvas();
        TCanvas *can_pt = new TCanvas();
        TCanvas *can_pl = new TCanvas();
        TCanvas *can_jt = new TCanvas();
        TCanvas *can_jl = new TCanvas();
        
        TCanvas *can_weights_vs_deltapt = new TCanvas("can_weights_vs_deltapt", "can_weights_vs_deltapt", 800, 500);
        TCanvas *can_weights_vs_deltajt = new TCanvas("can_weights_vs_deltajt", "can_weights_vs_deltajt", 800, 500);
        TCanvas *can_zj_vs_zi = new TCanvas("can_zj_vs_zi", "can_zj_vs_zi", 800, 500);

        if (deltap_bool) draw_save_del_hists(f_out, can_deltap, deltap_vec[k], "deltap", ptname, norm_string, hist_addname, false, true);
        if (deltapt_bool) draw_save_del_hists(f_out, can_deltapt, deltapt_vec[k], "deltapt", ptname, norm_string, hist_addname, false, true);
        if (deltapl_bool) draw_save_del_hists(f_out, can_deltapl, deltapl_vec[k], "deltapl", ptname, norm_string, hist_addname, false, true);
        if (ew_bool) draw_save_del_hists(f_out, can_weights, weights_vec[k], "weights", ptname, norm_string, hist_addname, false, true);
        
        if (deltajt_bool) draw_save_del_hists(f_out, can_deltajt, deltajt_vec[k], "deltajt", ptname, norm_string, hist_addname, false, true);
        if (deltajl_bool) draw_save_del_hists(f_out, can_deltajl, deltajl_vec[k], "deltajl", ptname, norm_string, hist_addname, false, true);
        
        if (deltap_bool) draw_save_del_hists(f_out, can_p, p_vec[k], "deltap", ptname, norm_string, hist_addname, false, true, false, "p");
        if (deltapt_bool) draw_save_del_hists(f_out, can_pt, pt_vec[k], "deltapt", ptname, norm_string, hist_addname, false, true, false, "pt");
        if (deltapl_bool) draw_save_del_hists(f_out, can_pl, pl_vec[k], "deltapl", ptname, norm_string, hist_addname, false, true, false, "pl");
        if (deltajt_bool) draw_save_del_hists(f_out, can_jt, jt_vec[k], "deltajt", ptname, norm_string, hist_addname, false, true, false, "jt");
        if (deltajl_bool) draw_save_del_hists(f_out, can_jl, jl_vec[k], "deltajl", ptname, norm_string, hist_addname, false, true, false, "jl");

        if (twoDhists_bool) draw_save_del_hists(f_out, can_weights_vs_deltapt, weights_vs_deltapt_hist2D, "weights_vs_deltapt", ptname, norm_string, hist_addname, false, false, true);
        if (twoDhists_bool) draw_save_del_hists(f_out, can_weights_vs_deltajt, weights_vs_deltajt_hist2D, "weights_vs_deltajt", ptname, norm_string, hist_addname, false, false, true);
        if (twoDhists_bool) draw_save_del_hists(f_out, can_zj_vs_zi, zj_vs_zi_hist2D, "zj_vs_zi", ptname, norm_string, hist_addname, false, false, true);
       
    }

    /* do pt bin stuff here */
    std::string hist_all_addname = weightstr + jetRname + thrname + "_pt" + ptname;

	// combine RL plots to get 1 plot per pt bin

    TCanvas *can_deltap_all = new TCanvas();
    TCanvas *can_deltapt_all = new TCanvas();
    TCanvas *can_deltapl_all = new TCanvas();
    TCanvas *can_weights_all = new TCanvas();

    TCanvas *can_deltajt_all = new TCanvas();
    TCanvas *can_deltajl_all = new TCanvas();

    TCanvas *can_p_all = new TCanvas();
    TCanvas *can_pt_all = new TCanvas();
    TCanvas *can_pl_all = new TCanvas();
    TCanvas *can_jt_all = new TCanvas();
    TCanvas *can_jl_all = new TCanvas();

    // // size_t length_deltap = deltap_vec.size();
    // // cout << " LENGTH DELTA P " << length_deltap << endl;
	
    if (deltap_bool) plotandsave_combined_hists(can_deltap_all, deltap_vec, leg, "deltap", ptname, norm_string, hist_all_addname, pt_max, logbins, true, -1);
    if (deltapt_bool) plotandsave_combined_hists(can_deltapt_all, deltapt_vec, leg, "deltapt", ptname, norm_string, hist_all_addname, pt_max, logbins, true, -1);
    if (deltapl_bool) plotandsave_combined_hists(can_deltapl_all, deltapl_vec, leg, "deltapl", ptname, norm_string, hist_all_addname, pt_max, logbins, true, -1);
    if (ew_bool) plotandsave_combined_hists(can_weights_all, weights_vec, leg, "weights", ptname, norm_string, hist_all_addname, pt_max, logbins, true, -1);

    if (deltajt_bool) plotandsave_combined_hists(can_deltajt_all, deltajt_vec, leg, "deltajt", ptname, norm_string, hist_all_addname, pt_max, logbins, true, -1);
    if (deltajl_bool) plotandsave_combined_hists(can_deltajl_all, deltajl_vec, leg, "deltajl", ptname, norm_string, hist_all_addname, pt_max, logbins, true, -1);
    
    if (deltap_bool) plotandsave_combined_hists(can_p_all, p_vec, leg, "deltap", ptname, norm_string, hist_all_addname, pt_max, logbins, true, -1, "p");
    if (deltapt_bool) plotandsave_combined_hists(can_pt_all, pt_vec, leg, "deltapt", ptname, norm_string, hist_all_addname, pt_max, logbins, true, -1, "pt");
    if (deltapl_bool) plotandsave_combined_hists(can_pl_all, pl_vec, leg, "deltapl", ptname, norm_string, hist_all_addname, pt_max, logbins, true, -1, "pl");
    if (deltajt_bool) plotandsave_combined_hists(can_jt_all, jt_vec, leg, "deltajt", ptname, norm_string, hist_all_addname, pt_max, logbins, true, -1, "jt");
    if (deltajl_bool) plotandsave_combined_hists(can_jl_all, jl_vec, leg, "deltajl", ptname, norm_string, hist_all_addname, pt_max, logbins, true, -1, "jl");
    
    // make graphs
    if (norm_string == "unnormalized" && rc_bool) {
        TCanvas *can_rc = new TCanvas();
        ProcessCanvas(can_rc);
        TGraphErrors *gr_rc = MakeFormatGraph(RLcenters_vec, rc_vec, kBlack, 1.0, markers[0], "R_{L}", "r_{c}", "rc", hist_all_addname);
        draw_save_del_hists(f_out, can_rc, gr_rc, "rc", ptname, norm_string, hist_all_addname, false, false);
        
        
        // save vectors here
        RL_vals.push_back(RLcenters_vec);
        rc_vals.push_back(rc_vec);
        rc_errors.push_back(rc_err_vec);
    }
}


//TODO: do something about norm_string!!
void analyze(TChain * JETINFO_tree, TChain * PAIRINFO_tree, 
             TFile *f_out, std::string weightstr, std::string jetRname, std::string thrname,
             std::string norm_string, const int pt_bins[], int n_bins, const double RL_bins[][7], int n_RLbins,
             bool include_RL0, bool include_RL1, bool debug, bool debug2 ) {

    
    // now look at observables and make histograms
    // don't separate by pt or RL bin
    if (norm_string == "unnormalized" && jetpt_bool) {
        TH1D * jetpt_hist = getObs1DHistFromTChain(JETINFO_tree, "jet_pt", 200, 0, 200, 0, 200, 0, 0);
        jetpt_hist->GetXaxis()->SetTitle("p_{T,jet}");
        TCanvas *can_jetpt = new TCanvas();
        draw_save_del_hists(f_out, can_jetpt, jetpt_hist, "jet_pt", "", "", weightstr + jetRname + thrname, false, true);
    
        TH1D * jet_const = getObs1DHistFromTChain(JETINFO_tree, "total_num_const", 20, 0, 20, 0, 200, 0, 0);
        jet_const->GetXaxis()->SetTitle("Number Constituents (total)");
        TCanvas *can_numconst = new TCanvas();
        draw_save_del_hists(f_out, can_numconst, jet_const, "total_num_const", "", "", weightstr + jetRname + thrname, false, true);
    
        TH1D * jet_const_aftercut = getObs1DHistFromTChain(JETINFO_tree, "num_const_aftercut", 20, 0, 20, 0, 200, 0, 0);
        jet_const_aftercut->GetXaxis()->SetTitle("Number Constituents (after threshold cut)");
        TCanvas *can_numconst_aftercut = new TCanvas();
        draw_save_del_hists(f_out, can_numconst_aftercut, jet_const_aftercut, "num_const_aftercut", "", "", weightstr + jetRname + thrname, false, true);
    
        
    
    }
    // return;

    //variables
    vector<vector<double>> RL_vals;
    vector<vector<double>> rc_vals;
    vector<double> ptcenter_bins;
    vector<vector<double>> rc_errors;

    const double pt_avgs[] = { 25.0009, 46.5139, 67.3 };

    // needs to be separated by pt and RL bin
    for ( int i = 0; i < n_bins; i++ ) {
        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];
        ptcenter_bins.push_back( (pt_min+pt_max)/2 );
           
        if (debug) cout << " in pt bin" << i << " with " << pt_min << " - " << pt_max << endl;
        analyze_ptbin(JETINFO_tree, PAIRINFO_tree, f_out, weightstr, jetRname, thrname, norm_string, pt_min, pt_max, pt_avgs[i], RL_bins[i], n_RLbins, RL_vals, rc_vals, rc_errors, include_RL0, include_RL1, debug, debug2);
        
    }


    // plot r_c as a function of RL
    if (norm_string == "unnormalized" && rc_bool) {
        cout << "checkpoint 4" << endl;
        cout << "size of RL_vals " << RL_vals.size() << endl;
        cout << "size of RL_vals[0] " << RL_vals[0].size() << endl;
        cout << "size of rc_vals " << rc_vals.size() << endl;
        cout << "size of rc_vals[0] " << rc_vals[0].size() << endl;
        cout << "size of ptcenter_bins " << ptcenter_bins.size() << endl;

        plot_rc(RL_vals, rc_vals, ptcenter_bins, rc_errors); //, leg_RLbins, leg_ptbins);
    }
    


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

    // update dir names
    // if (ptrl_bins == false) attempt_dir = Form("data_secondattempt/rebinx%d", rebin);

    attempt_dir = "data_fourthattempt/";
    if (ptrl_bins == true) attempt_dir = attempt_dir.substr(0,attempt_dir.size()-1) + "_ptrlbins/";
    if (logbins == true) {
        rebin = 0;
        attempt_dir += "logbins";
    } else attempt_dir += Form("rebinx%d", rebin);
    
    outdir = "/software/users/blianggi/mypyjetty/storage/dEEC/plots/" + attempt_dir;
    
    // ntuple/histogram names
    std::string JETINFO_name = "tn_JETINFO_R0.4_1.0";
    std::string PAIRINFO_name = "tn_pairlevel_R0.4_1.0";
    std::string jet1D_name = "h_1Djet_pt_JetPt_R0.4_1.0"; // this one is a histogram
        
    // filenames
    std::string filename = Form("~/Documents/research/othercorrelations/data_ntuples/AnalysisResults_0001.root");
    std::string base_filepath_perly = Form("/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/31843529");
    std::string base_filepath_hic = Form("/rstorage/alice/AnalysisResults/blianggi/dEEC/468247"); //442528");
    
    // Output file for binned results
    std::string root_outfile = "/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/" + attempt_dir + "/DataHists.root"; //plots/ntuples/DataHists.root"; //FinalDataHists.root
    TFile* f_out = new TFile(root_outfile.c_str(), "RECREATE");
    std::string add_name = ""; //"_othercorrel";


    // analysis variables
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
    std::ifstream filelist("/software/users/blianggi/mypyjetty/dEEC/filelist_datatuples_468247_shortname.txt");
    if (!filelist.is_open()) {
        std::cerr << "Error: Could not open /software/users/blianggi/mypyjetty/dEEC/filelist_datatuples_468247_shortname.txt" << std::endl;
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
    
    // analyze for plots
    // norm_string = "unnormalized";
    // if (unnormalized_bool) analyze(JETINFO_tree, PAIRINFO_tree, f_out, weightstr, jetRname, thrname, norm_string, pt_bins, n_bins, RL_bins, n_RLbins, include_RL0, include_RL1, debug, debug2);
    
    norm_string = "self_normalized";
    if (self_normalized_bool) analyze(JETINFO_tree, PAIRINFO_tree, f_out, weightstr, jetRname, thrname, norm_string, pt_bins, n_bins, RL_bins, n_RLbins, include_RL0, include_RL1, debug, debug2);
    
    // norm_string = "norm_by_jets";
    // if (norm_by_jets_bool) analyze(JETINFO_tree, PAIRINFO_tree, f_out, weightstr, jetRname, thrname, norm_string, pt_bins, n_bins, RL_bins, n_RLbins, include_RL0, include_RL1, debug, debug2);
    


    // delete objects after saving for new pt-hat bin
    delete JETINFO_tree;
    delete PAIRINFO_tree;

    f_out->Close();
    
  
}


