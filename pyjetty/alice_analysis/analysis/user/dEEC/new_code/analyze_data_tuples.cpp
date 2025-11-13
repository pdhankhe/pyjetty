// ROOT macro to analyze and plot data tuples 
// This version will only accomodate the FOUR PTRL bins
// Beatrice Liang-Gilman (beatrice_lg@berkeley.edu)

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <variant>

using namespace std;

// global variables
Double_t colors[16] = {kGray, kMagenta, kBlue, kOrange+1, kViolet+1, kGreen+2, kRed, kYellow+1, kCyan+1};
Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
Double_t marker_size = 1.5;

bool logbins = false;
std::string attempt_dir; // = Form("data_secondattempt/rebinx%d", rebin);
std::string outdir; // = "/software/users/blianggi/mypyjetty/storage/dEEC/plots/" + attempt_dir;

int filecounter_cutoff = 100; //-1; //total: 845 // IDK why I thought it was 7601?
bool write_to_root_file = true; 

// bool jetpt_bool = false;
bool deltap_bool = false;
bool deltapt_bool = false;
bool deltajt_bool = false;
bool ew_bool = false;
bool twoDhists_bool = false;
bool rc_bool = false; 
bool deltajt_vs_ptrl_bool = true;

bool p1_bool = false;
bool jt1_bool = false;

bool unnormalized_bool = false;
bool self_normalized_bool = true;
bool norm_by_jets_bool = false;

bool unweighted_bool = true;
bool weighted_bool = false;


class Observable {
public:
    std::string name;
    bool obs_bool;
    
    int num_bins;
    double min_bound;
    double max_bound;
    
    std::string axis_label;
    std::string cs_label; //cross section label
    std::string filepath_plots;

    std::vector<TH1D*> obs_vec;

    Observable(std::string name_val, bool obs_bool_val, int num_bins_val, double min_bound_val, double max_bound_val, 
               std::string axis_label_val, std::string cs_label_val) {
        name = name_val;
        obs_bool = obs_bool_val;
        
        num_bins = num_bins_val;
        min_bound = min_bound_val;
        max_bound = max_bound_val;

        axis_label = axis_label_val;
        cs_label = cs_label_val; //cross section label, in y axis

        filepath_plots = "/software/users/blianggi/mypyjetty/storage/dEEC/plots/" + attempt_dir + "%s/%s/" + name + "/%s"; // ptname, norm_string, filename
        if (name.find("jet_") != std::string::npos || name.find("const") != std::string::npos) filepath_plots = "/software/users/blianggi/mypyjetty/storage/dEEC/plots/" + attempt_dir + "%s"; // filename

    }

    void addHist(TH1D* hist) {
        obs_vec.push_back(hist);
    }

    void recreate_output_root_file() {
        std::string root_outfile = "/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/" + attempt_dir + "/DataHists_" + name + ".root";
        TFile * f_out = new TFile(root_outfile.c_str(), "RECREATE");
        f_out->Close();
    }

    TFile * get_output_root_file() {
        std::string root_outfile = "/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/" + attempt_dir + "/DataHists_" + name + ".root";
        TFile * f_out = new TFile(root_outfile.c_str(), "UPDATE");
        return f_out;
    }
};

class Observable2D {
public:
    Observable obsx;
    Observable obsy;
    std::string name;
    bool obs_bool;
    
    std::string filepath_plots;

    // Observable2D(const Observable& o1, const Observable& o2, bool obs_bool_val)
    //     : obs1(o1), obs2(o2), name(name_val), obs_bool(obs_bool_val) {}
    Observable2D(Observable ox, Observable oy, bool obs_bool_val)
        : obsx(ox), obsy(oy), obs_bool(obs_bool_val) // required for non-default-constructible members
    {
        name = obsy.name + "_vs_" + obsx.name;
        filepath_plots = "/software/users/blianggi/mypyjetty/storage/dEEC/plots/" + attempt_dir + "%s/%s/" + name + "/%s"; // ptname, norm_string, filename
    }

    void recreate_output_root_file() {
        std::string root_outfile = "/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/" + attempt_dir + "/DataHists_" + name + ".root";
        TFile * f_out = new TFile(root_outfile.c_str(), "RECREATE");
        f_out->Close();
    }

    TFile* get_output_root_file() {
        std::string root_outfile = "/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/" + attempt_dir + "/DataHists_" + name + ".root";
        TFile * f_out = new TFile(root_outfile.c_str(), "UPDATE");
        return f_out;
    }
};


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

bool legendContainsLabel(TLegend& leg, const TString& label) {
    TList* list = leg.GetListOfPrimitives();
    for (int i = 0; i < list->GetSize(); ++i) {
        TLegendEntry* entry = dynamic_cast<TLegendEntry*>(list->At(i));
        if (entry && entry->GetLabel() && label == entry->GetLabel()) {
            return true;
        }
    }
    return false;
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
TH1D * getObs1DHistFromTChain(TChain *chain, Observable obs, int pt_min, int pt_max, 
                              double ptRL_min, double ptRL_max, double ptavg=0.0, bool weighting=false,
                              bool logbinning=false, bool debug=false) {
    
    cout << "OBSERVABLE NAME1!!!!" << obs.name << " // " << endl;

    // initializing the histogram
    TH1D * hist1D;
    hist1D = new TH1D(Form("%s_hist", obs.name.c_str()), Form("%s_hist", obs.name.c_str()), obs.num_bins, obs.min_bound, obs.max_bound);
    
    // if (logbinning == true)
    //     double new_hist_xmin = (hist_xmin == 0.0) ? 1e-1 : hist_xmin;
    //     if (branch_name == "p1" || branch_name == "pt1" || branch_name == "jl1") new_hist_xmin = 0.8;
    //     std::vector<double> bins = makeLogBins(new_hist_xmin, hist_xmax, num_bins);
    //     hist1D = new TH1D(Form("%s_hist", branch_name.c_str()), Form("%s_hist", branch_name.c_str()), num_bins, bins.data());
    // }

    // filling the histogram
    if ( obs.name == "jet_pt" || obs.name == "total_num_const" || obs.name == "num_const_aftercut") {
        chain->Draw(Form("%s>>%s_hist", obs.name.c_str(), obs.name.c_str()), Form("jet_pt >= %d && jet_pt < %d", pt_min, pt_max), "e");
    } else {
        if (weighting == false) chain->Draw(Form("%s>>%s_hist", obs.name.c_str(), obs.name.c_str()), Form("jet_pt >= %d && jet_pt < %d && %f*RL >= %f && %f*RL < %f", pt_min, pt_max, ptavg, ptRL_min, ptavg, ptRL_max), "e");
        else {
            Float_t obs_name; Float_t weight;
            Float_t jetpt; Float_t rl;
            if (obs.name != "weights") chain->SetBranchAddress(Form("%s",obs.name.c_str()), &obs_name); // this won't work if obs.name == weights, bc it will overwrite in next row
            chain->SetBranchAddress("weights", &weight);
            chain->SetBranchAddress("jet_pt", &jetpt);
            chain->SetBranchAddress("RL", &rl);

            int entries = chain->GetEntries();
            for ( int i = 0; i < entries; i++ ) {
                chain->GetEntry(i);
                if (jetpt >= pt_min && jetpt < pt_max && (ptavg*rl) >= ptRL_min && (ptavg*rl) < ptRL_max) {
                    if (obs.name == "weights") obs_name = weight;
                    hist1D->Fill(obs_name, weight);
                }
            }

            chain->ResetBranchAddresses();
        }
        
    }

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
vector<double> getRcFromTChain(TChain *chain, Observable obs, int pt_min, int pt_max, 
                              double ptRL_min, double ptRL_max, double ptavg, bool weighting=false) 
{
    // draw regular charge histogram, where like sign = +1, and unlike sign = -1
    // TH1D *hist_charge = new TH1D("hist_charge", "hist_charge", obs.num_bins, obs.min_bound, obs.max_bound);
    // chain->Draw("q1q2>>hist_charge", Form("jet_pt >= %d && jet_pt < %d && %f*RL >= %f && %f*RL < %f", pt_min, pt_max, ptavg, ptRL_min, ptavg, ptRL_max), "e");
    TH1D *hist_charge = getObs1DHistFromTChain(chain, obs, pt_min, pt_max, ptRL_min, ptRL_max, ptavg, weighting);

    // get # of like sign and # of unlike sign
    double num_likesign = hist_charge->GetBinContent(hist_charge->FindBin(1));
    double num_unlikesign = hist_charge->GetBinContent(hist_charge->FindBin(-1));
    // if (debug) cout << "num like sign " << num_likesign << " num unlike sign " << num_unlikesign << endl;

    // calculate the rc value for this pt & RL bin
    double rc = (double)(num_likesign - num_unlikesign) / (double)(num_likesign + num_unlikesign);
    //if (debug) cout << "and that makes rc " << rc << endl;
    
    // calculate the rc error for this pt & RL bin
    double num_totalpairs = num_likesign + num_unlikesign;
    double rc_err = ( 2 * sqrt( num_totalpairs * num_likesign * num_unlikesign ) ) / (num_totalpairs * num_totalpairs);
    cout << "RC ERR IN FUNC IS " << rc_err << endl;

    delete hist_charge;
    return { rc, rc_err };
}



/* get a typical 2D histogram from the TChain */
TH2D * getObs2DHistFromTChain(TChain *chain, Observable obs_x, Observable obs_y, 
                              int pt_min, int pt_max, double ptRL_min, double ptRL_max, double ptavg) {

    TH2D * hist2D = new TH2D(Form("%s_vs_%s_hist", obs_y.name.c_str(), obs_x.name.c_str()), Form("%s_vs_%s_hist", obs_y.name.c_str(), obs_x.name.c_str()), obs_x.num_bins, obs_x.min_bound, obs_x.max_bound, obs_y.num_bins, obs_y.min_bound, obs_y.max_bound);
    if (obs_x.name == "z" && obs_y.name == "z") { // zi vs zj
        chain->Draw("pt2/jet_pt:pt1/jet_pt>>zj_vs_zi_hist", Form("jet_pt >= %d && jet_pt < %d && %f*RL >= %f && %f*RL < %f", pt_min, pt_max, ptavg, ptRL_min, ptavg, ptRL_max));
    } else if (obs_x.name == "zsmol") { //zbig vs zsmol
        cout << "this isn't implemented yet! (idk how to )" << endl;
        return hist2D;
    } else if (obs_x.name == "maxpt") { //weights vs maxpt
        cout << "this isn't implemented yet! (idk how to )" << endl;
        return hist2D;
    } else if ( obs_x.name == "ptrl" || obs_y.name == "deltajt") {
        cout << "this is delta jt vs ptrl" << endl;
        Float_t deltajt; Float_t weight; 
        Float_t jetpt; Float_t rl;
        chain->SetBranchAddress("deltajt", &deltajt); // this won't work if obs.name == weights, bc it will overwrite in next row
        chain->SetBranchAddress("weights", &weight);
        chain->SetBranchAddress("jet_pt", &jetpt);
        chain->SetBranchAddress("RL", &rl);

        int entries = chain->GetEntries();
        for ( int i = 0; i < entries; i++ ) {
            chain->GetEntry(i);
            if (jetpt >= pt_min && jetpt < pt_max) {
                hist2D->Fill(ptavg*rl, deltajt); //unweighted only!!!!
            }
        }

        chain->ResetBranchAddresses();
    } else {
        chain->Draw(Form("%s:%s>>%s_vs_%s_hist", obs_y.name.c_str(), obs_x.name.c_str(), obs_y.name.c_str(), obs_x.name.c_str()), Form("jet_pt >= %d && jet_pt < %d && %f*RL >= %f && %f*RL < %f", pt_min, pt_max, ptavg, ptRL_min, ptavg, ptRL_max));
    }
    return hist2D;
}


/* Format and adjust histograms */
void Format1DHist(Observable obs, TH1D *hist, TH1D *jetpt_hist, std::string norm_string, int markercolor, double markeralpha,
                  int markerstyle, std::string xtitle, std::string ytitle, TLegend& leg, TString leg_text, 
                  std::string hist_addname="") {

    hist->SetTitle(Form("h_%s%s", obs.name.c_str(), hist_addname.c_str()));
    hist->SetName(Form("h_%s%s", obs.name.c_str(), hist_addname.c_str()));

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
	if (obs.name == "weights") {
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
    if (!legendContainsLabel(leg, leg_text)) leg.AddEntry(hist, leg_text, "pl");

}

void Format2DHist(Observable obs_x, Observable obs_y, TH2D *hist2D, TH1D *jetpt_hist, std::string norm_string, std::string xtitle, std::string ytitle, 
                  std::string hist_addname="", bool restrictzaxis=false) {

    double ptvsew_normbounds[3][2] = { { 0, 500 }, { 1e-2, 5e2 }, { 1e-4, 2e2 }}; //TODO: make this less pt specific?
    // also TODO: potentially set all lower boundaries to 0 for data??

    hist2D->SetTitle(Form("h_%s_vs_%s%s", obs_y.name.c_str(), obs_x.name.c_str(), hist_addname.c_str()));
    hist2D->SetName(Form("h_%s_vs_%s%s", obs_y.name.c_str(), obs_x.name.c_str(), hist_addname.c_str()));

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

    // set z axis bounds
    // if (restrictzaxis) hist2D->GetZaxis()->SetRangeUser(ptvsew_normbounds[norm_index][0], ptvsew_normbounds[norm_index][1]);

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


void FormatGraphMarker(TGraphErrors * g, int markercolor, double markeralpha, int markerstyle, double markersize) {
    g->SetMarkerColorAlpha(markercolor, markeralpha);
    g->SetMarkerStyle(markerstyle);
    g->SetMarkerSize(markersize);
    g->SetLineColorAlpha(markercolor, markeralpha);
}

TGraphErrors * MakeFormatGraph(vector<double> xvals, vector<double> yvals, vector<double> yval_errors, int markercolor, double markeralpha,
                  int markerstyle, std::string xtitle, std::string ytitle, std::string obs_name, std::string hist_addname="") {

    vector<double> xval_errors;
    for (int i=0; i<xvals.size(); i++) xval_errors.push_back(0);
    
    TGraphErrors * graph = new TGraphErrors(xvals.size(), xvals.data(), yvals.data(), xval_errors.data(), yval_errors.data());
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

    graph->GetYaxis()->SetLabelFont(42);
	graph->GetYaxis()->SetTitleFont(42);
    graph->GetYaxis()->SetTitleOffset(1.05); 
	graph->GetYaxis()->SetTitleSize(0.06); //(0.042);
	graph->GetYaxis()->SetLabelSize(0.05); //(0.042);

    return graph;
}

TLine * drawVertLine(double x1, double y1, double y2, int color, int linestyle=2){
    auto fvertline = new TLine(x1, y1, x1, y2);
	fvertline->SetLineWidth(1);
    fvertline->SetLineColor(color);
    fvertline->SetLineStyle(linestyle);
    return fvertline;
}

TLine * drawHoriLine(double x1, double x2, double y1, int color, int linestyle=2){
    auto fhoriline = new TLine(x1, y1, x2, y1);
	fhoriline->SetLineWidth(1);
    fhoriline->SetLineColor(color);
    fhoriline->SetLineStyle(linestyle);
    return fhoriline;
}

/* Save and delete histograms */
// (TFile *fout, TCanvas *can, TH1D* hist, std::string obs_name, 
//                          std::string ptname, std::string norm_string, std::string hist_addname,
//                          bool logx, bool logy)
using ObservableVariant = std::variant<Observable, Observable2D>;
void draw_save_del_hists(Observable& obs, TObject* obj,  
                         std::string ptname, std::string norm_string, std::string hist_addname,
                         bool logx, bool logy, bool logz=false, //std::string obs_filename="", 
                         const double ptRL_bins[] = nullptr, int n_ptRLbins = 0,
                         bool include_RL0=true, bool include_RL1=true) {
    
    TCanvas *can = new TCanvas();
    can->cd();
    if (logx) gPad->SetLogx();
    if (logy) gPad->SetLogy();

    if (TH2* hist2D = dynamic_cast<TH2*>(obj)) { // put this first bc TH2 is a subclass of TH1!! (and it will go into the other loop :( )
        cout << "TH2D should not be running in this function for now!" << endl;
        return;
    } else if (TH1* hist = dynamic_cast<TH1*>(obj)) {
        hist->Draw();
    } else if (TGraphErrors* graph = dynamic_cast<TGraphErrors*>(obj)) {
        graph->Draw("ALP");
    } else if (TGraph* graph = dynamic_cast<TGraph*>(obj)) { // this is specifically for delta jt scatter
        gPad->DrawFrame(3E-2, 0., 35., 5.); //xmin, ymin, xmax, ymax
        graph->SetMinimum(0.);
        graph->SetMaximum(5.);
        graph->GetXaxis()->SetTitle(obs.axis_label.c_str()); //"#LTp_{T}#GTR_{L}");
        graph->GetXaxis()->SetTitle(obs.cs_label.c_str()); //"#Deltaj_{T}");
        graph->Draw("AP");
        cout << ptRL_bins << endl;
        if ( ptRL_bins != nullptr ) { // (n_ptRLbins != 0)
            cout << "in function" << endl;
            // int n_ptRLbins = sizeof(ptRL_bins) / sizeof(ptRL_bins[0]);
            int y1 = 0; int y2 = 5;
            for ( int j = 0; j < n_ptRLbins; j++ ) {
                int k = j;
                if (!include_RL0) {
                    k = j-1;
                    if (j == 0) continue; // can add something here to change the filename for ALL
                }
                if (!include_RL1 && j == n_ptRLbins-1) continue;
                
                cout << "RL BINS HERE!" << ptRL_bins[j] << endl;
                TLine * vertline = drawVertLine(ptRL_bins[j], y1, y2, kViolet, 1);
                vertline->Draw("SAME");
            }
        }
    } else {
        cout << "Error: Unsupported object type. Only TH1, TGraph, TGraphErrors, and TH2 are supported." << endl;
        return;
    }
    // hist->Draw();

    TFile * fout = obs.get_output_root_file();
    fout->cd();
    if (write_to_root_file) obj->Write(); //TODO: this might not be right! Might have to use the casted type
    fout->Close();


    // std::string outdir = "/software/users/blianggi/mypyjetty/storage/dEEC/plots/data_firstattempt"; // + ptbin_name + "/";//"plots/test/";
    /*std::string add_dir = "";
    if (obs_name != "jet_pt" && obs_name != "total_num_const" && obs_name != "num_const_aftercut") {
        if (obs_name == "rc") add_dir = "/" + ptname + "/" + norm_string + "/" + obs_name;
        else if (obs_name == "deltajt_scatter") add_dir = "/" + ptname;
        // else if (obs_name == "deltajt_scatter_ind") add_dir = "/" + ptname + "self_normalized/deltajt";
        else add_dir = "/" + ptname + "/" + norm_string + "/" + obs_name + "/individuals";
    }
    std::string usethisname = (obs_filename == "") ? obs_name : obs_filename;
    
    std::string fname_out = outdir + add_dir + "/corrhist_" + usethisname + hist_addname + ".pdf";*/

    std::string fname_out = "";
    if (obs.name.find("jet_") != std::string::npos || obs.name.find("const") != std::string::npos) fname_out = Form(obs.filepath_plots.c_str(), ("corrhist_" + obs.name + hist_addname + ".pdf").c_str());
    else {
        fname_out = Form(obs.filepath_plots.c_str(), ptname.c_str(), norm_string.c_str(), ("corrhist_" + obs.name + hist_addname + ".pdf").c_str());
    }
    can->SaveAs(fname_out.c_str());

    // delete hist;
    delete can;
}

void draw_save_del_hists2D(Observable2D& obs2D, TH2D* hist2D,  
                         std::string ptname, std::string norm_string, std::string hist_addname,
                         bool logx, bool logy, bool logz=false,
                         const double ptRL_bins[] = nullptr, int n_ptRLbins = 0,
                         bool include_RL0=true, bool include_RL1=true) {
    
    TCanvas *can_hist2D = new TCanvas("", "", 800, 500); //"can_weights_vs_deltap"
    can_hist2D->cd();
    if (logx) gPad->SetLogx();
    if (logy) gPad->SetLogy();

    gPad->SetRightMargin(0.12);
    if (logz) gPad->SetLogz();
    can_hist2D->SetFillColor(kWhite);
    hist2D->Draw("COLZ");
    

    TFile * fout = obs2D.get_output_root_file();
    fout->cd();
    if (write_to_root_file) hist2D->Write(); //TODO: this might not be right! Might have to use the casted type
    fout->Close();


    std::string fname_out = Form(obs2D.filepath_plots.c_str(), ptname.c_str(), norm_string.c_str(), ("corrhist_" + obs2D.name + hist_addname + ".pdf").c_str());
    if (obs2D.name.find("z") != std::string::npos) fname_out = Form(obs2D.filepath_plots.c_str(), ("corrhist_zj_vs_zi" + hist_addname + ".pdf").c_str());
    can_hist2D->SaveAs(fname_out.c_str());

    // delete hist;
    delete can_hist2D;
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
void plotandsave_combined_hists(Observable obs, TLegend *l, std::string ptname, std::string norm_string, 
                          std::string hist_addname, int pt_max, bool logx, bool logy, bool debug=false) {

    // go into canvas
    TCanvas *can_all = new TCanvas();
    can_all->cd();
    if (logx) gPad->SetLogx();
    if (logy) gPad->SetLogy();

    // set min and max
    size_t num_curves = obs.obs_vec.size();
    double max_val = 0.0; double min_val = 999.0;
    for (int j=0; j<num_curves; j++) {
        double max_temp = obs.obs_vec[j]->GetMaximum();
        double min_temp = obs.obs_vec[j]->GetMinimum();
        if (max_temp > max_val) max_val = max_temp;
        if (min_temp < min_val) min_val = min_temp;
    }
    obs.obs_vec[0]->SetMaximum(max_val*1.2);
    // obs.obs_vec[0]->SetMinimum(min_val*1.2);
    cout << "max value is: " << max_val << endl;
    cout << "min value is: " << min_val << endl;
    
    // draw!
    for (int j=0; j<num_curves; j++) {
        obs.obs_vec[j]->Draw("same");
    }

    l->Draw("same");

    //save as PDF
    std::string fname_out = Form(obs.filepath_plots.c_str(), ptname.c_str(), norm_string.c_str(), ("corrhist_" + obs.name + "_ALL" + hist_addname + ".pdf").c_str());
    can_all->SaveAs(fname_out.c_str());

    // deleteVecOfHists(h_vec);
    delete can_all;

}



// ======================================================= //
//                   SPECIFIC FUNCTIONS
// ======================================================= //

/* // TODO: work on this! commented out so i could run code
void func() {
    TCanvas *can = new TCanvas("can", "can", 750, 500);
    can->cd();
    for ( int i = 0; i < ptcenters_vec.size(); i++ ) {
        if ( i == 0 ) {
            rc_graphs_func_of_ptRL[i]->SetMinimum(-0.5);  // Lower y limit
            rc_graphs_func_of_ptRL[i]->SetMaximum(0.025); 
            rc_graphs_func_of_ptRL[i]->GetXaxis()->SetTitle("R_{L} bin center"); //change!!
            rc_graphs_func_of_ptRL[i]->GetYaxis()->SetTitle("r_{c}"); 
            rc_graphs_func_of_ptRL[i]->Draw("AP");
        } 

        for (int j = 0; j < ptRLcenters_vec.size(); j++) {
            rc_graphs_func_of_ptRL_ind[i][j]->Draw("P SAME");
        }
    }
    leg_RLbins.Draw("same");
    leg_ptbins.Draw("same");
    can->SaveAs(fname_func_of_RL_out.c_str());
    delete can;

}
*/

// void plot_rc(vector<double>& ptcenters_vec, vector<double>& ptRLcenters_vec, vector<vector<double>>& rc_vec,
//              vector<vector<double>>& rc_errors_vec) {
//             //  TLegend& leg_RLbins, TLegend& leg_ptbins) {
    
//     vector<TGraphErrors *> rc_graphs_func_of_ptRL;
//     vector<TGraphErrors *> rc_graphs_func_of_pT;
//     vector<vector<TGraphErrors*>> rc_graphs_func_of_ptRL_ind;
//     vector<vector<TGraphErrors*>> rc_graphs_func_of_pT_ind;

//     vector<double> ptRL_err;
//     vector<double> pt_err;
//     for (int i=0; i<ptcenters_vec.size(); i++) pt_err.push_back(0);
//     for (int j=0; j<ptRLcenters_vec.size(); j++) ptRL_err.push_back(0);

//     cout <<" ptRL_err size " << ptRL_err.size() << endl;
//     cout <<" pt_err size " << pt_err.size() << endl;

//     // std::string outdir = "/software/users/blianggi/mypyjetty/storage/dEEC/plots/data_firstattempt"; // + ptbin_name + "/";//"plots/test/";
//     std::string fname_func_of_RL_out = outdir + "/corrhist_rc_func_of_RL.pdf"; // could add jetR and threshold info later??, maybe not needed tho 
//     std::string fname_func_of_pT_out = outdir + "/corrhist_rc_func_of_pT.pdf";

//     TLegend leg_RLbins(0.2, 0.2, 0.4, 0.45); 
//     TLegend leg_ptbins(0.5, 0.3, 0.65, 0.45); 

//     leg_RLbins.SetTextSize(0.037);
//     leg_RLbins.SetBorderSize(0);
//     leg_ptbins.SetTextSize(0.037);
//     leg_ptbins.SetBorderSize(0);

//     // delete hist;

//     // get graphs of r_c as a function of RL
//     // loop over pt bins
//     for ( int i = 0; i < ptcenters_vec.size(); i++ ) { 

//         // for graphs as a function of ptRL
//         TGraphErrors *g = new TGraphErrors(ptRLcenters_vec.size(), ptRLcenters_vec.data(), rc_vec[i].data(), ptRL_err.data(), rc_errors_vec[i].data());
//         FormatGraphMarker(g, kBlack, 1.0, markers[i], marker_size);
//         leg_ptbins.AddEntry(g, Form("p_{T} = %d-%d", (int)ptcenters_vec[i]-10, (int)ptcenters_vec[i]+10), "P");
        
//         vector<TGraphErrors*> ind_temp_vec;
//         for (int j=0; j<ptRLcenters_vec.size(); j++){
//             int k=j+1;
//             TGraphErrors *g_ind = new TGraphErrors(1, &ptRLcenters_vec[j], &rc_vec[i][j], &ptRL_err[j], &rc_errors_vec[i][j]);
//             cout << "i: " << i << " j: " << j << " err: " << rc_errors_vec[i][j] << endl;
//             FormatGraphMarker(g_ind, colors[k], 1.0, markers[i], marker_size);
//             ind_temp_vec.push_back(g_ind);
            
//             if (i==0) leg_RLbins.AddEntry(g_ind, Form("R_{L} bin %d", j+1), "P");
            
//         }
//         rc_graphs_func_of_ptRL.push_back(g); 
//         rc_graphs_func_of_ptRL_ind.push_back(ind_temp_vec); 
                 

//     }

//     // plot r_c as a function of ptRL
//     func();
//     TCanvas *can_func_of_ptRL = new TCanvas("can_func_of_ptRL", "can_func_of_ptRL", 750, 500);
//     can_func_of_ptRL->cd();
//     for ( int i = 0; i < ptcenters_vec.size(); i++ ) {
//         if ( i == 0 ) {
//             rc_graphs_func_of_ptRL[i]->SetMinimum(-0.5);  // Lower y limit
//             rc_graphs_func_of_ptRL[i]->SetMaximum(0.025); 
//             rc_graphs_func_of_ptRL[i]->GetXaxis()->SetTitle("R_{L} bin center"); 
//             rc_graphs_func_of_ptRL[i]->GetYaxis()->SetTitle("r_{c}"); 
//             rc_graphs_func_of_ptRL[i]->Draw("AP");
//         } 

//         for (int j = 0; j < ptRLcenters_vec.size(); j++) {
//             rc_graphs_func_of_ptRL_ind[i][j]->Draw("P SAME");
//         }
//     }
//     leg_RLbins.Draw("same");
//     leg_ptbins.Draw("same");
//     can_func_of_ptRL->SaveAs(fname_func_of_RL_out.c_str());
//     delete can_func_of_ptRL;

//     //========================================================

//     // get graphs of r_c as a function of pT
//     vector<vector<double>> rc_vals_func_of_pT;
//     vector<vector<double>> rc_err_vals_func_of_pT;
//     // loop over RL bins
//     for ( int j = 0; j < ptRLcenters_vec.size(); j++ ) {
//         vector<double> temp_vec;
//         vector<TGraphErrors*> ind_temp_vec;
//         vector<double> err_temp_vec;

//         // save values into appropriate vectors
//         for ( int i = 0; i < ptcenters_vec.size(); i++ ) { 
//             temp_vec.push_back(rc_vec[i][j]);
//             err_temp_vec.push_back(rc_errors_vec[i][j]);
            
//             int k=j+1;
//             TGraphErrors *g_ind = new TGraphErrors(1, &ptcenters_vec[i], &rc_vec[i][j], &pt_err[i], &rc_errors_vec[i][j]);
//             cout << "i: " << i << " j: " << j << " err: " << rc_errors_vec[i][j] << endl;
//             FormatGraphMarker(g_ind, colors[k], 1.0, markers[i], marker_size);
//             ind_temp_vec.push_back(g_ind);
//         }
//         rc_vals_func_of_pT.push_back(temp_vec);
//         rc_err_vals_func_of_pT.push_back(err_temp_vec);
//         rc_graphs_func_of_pT_ind.push_back(ind_temp_vec); 

//         // for graphs as a function of RL
//         TGraphErrors *g = new TGraphErrors(ptcenters_vec.size(), ptcenters_vec.data(), rc_vals_func_of_pT[j].data(), pt_err.data(), rc_err_vals_func_of_pT[j].data());
//         for (int aa = 0; aa < ptcenters_vec.size(); aa++) {
//             // cout << "studying pt=" << ptcenters_vec[aa] << " // " << rc_vals_func_of_pT[j][aa] << endl;
//         }
//         rc_graphs_func_of_pT.push_back(g); 
//     }
    
//     // plot r_c as a function of pT
//     TCanvas *can_func_of_pT = new TCanvas("can_func_of_pT", "can_func_of_pT", 750, 500);
//     can_func_of_pT->cd();
//     for ( int j = 0; j < ptRLcenters_vec.size(); j++ ) {
//         int k = j+1;
//         // rc_graphs_func_of_pT[j]->SetMarkerSize(1.0);
//         // rc_graphs_func_of_pT[j]->SetMarkerStyle(markers[0]);
//         rc_graphs_func_of_pT[j]->SetMarkerColorAlpha(colors[k], 0.0);
//         if ( j == 0 ) {
//             rc_graphs_func_of_pT[j]->SetMinimum(-0.5);  // Lower y limit
//             rc_graphs_func_of_pT[j]->SetMaximum(0.025); 
//             rc_graphs_func_of_pT[j]->GetXaxis()->SetTitle("p_{T} bin center"); 
//             rc_graphs_func_of_pT[j]->GetYaxis()->SetTitle("r_{c}"); 
//             rc_graphs_func_of_pT[j]->Draw("AP");
//         } // else { 
//         //     rc_graphs_func_of_pT[j]->Draw("P SAME");
//         // }

//         for (int i = 0; i < ptcenters_vec.size(); i++) {
//             rc_graphs_func_of_pT_ind[j][i]->Draw("P SAME");
//         }
//     }
//     leg_RLbins.Draw("same");
//     leg_ptbins.Draw("same");
//     can_func_of_pT->SaveAs(fname_func_of_pT_out.c_str());
//     delete can_func_of_pT;



// }

// xaxis can be "pt" or "ptrl"
void idkyet(std::string xaxis, vector<double>& x_vec, vector<double>& x_err,
            vector<vector<double>>& rc_vec, vector<vector<double>>& rc_errors_vec) {

    int num_ptbins = rc_vec.size();
    int num_ptrlbins = rc_vec[0].size();

    for ( int i = 0; i < num_ptbins; i++ ) { // loop over number pt bins
        for ( int j = 0; j < num_ptrlbins; j++ ) { // loop over number of ptRL bins
            TGraphErrors * g_ind = new TGraphErrors(1, &x_vec[i], &rc_vec[i][j], &x_err[i], &rc_errors_vec[i][j]);
        }
    }


    

}

/* // TODO: work on this! commented out so i could run code
void plot_rc(vector<double>& ptcenters_vec, vector<double>& ptRLcenters_vec, vector<vector<double>>& rc_vec,
             vector<vector<double>>& rc_errors_vec) {
            //  TLegend& leg_RLbins, TLegend& leg_ptbins) {
    
    vector<TGraphErrors*> ptbin_graph_labels;
    vector<vector<TGraphErrors*>> rc_graphs_func_of_ptRL_ind;
    vector<vector<TGraphErrors*>> rc_graphs_func_of_pT_ind;

    vector<double> ptRL_err;
    vector<double> pt_err;
    for (int i=0; i<ptcenters_vec.size(); i++) pt_err.push_back(0);
    for (int j=0; j<ptRLcenters_vec.size(); j++) ptRL_err.push_back(0);

    cout <<" ptRL_err size " << ptRL_err.size() << endl;
    cout <<" pt_err size " << pt_err.size() << endl;

    // std::string outdir = "/software/users/blianggi/mypyjetty/storage/dEEC/plots/data_firstattempt"; // + ptbin_name + "/";//"plots/test/";
    std::string fname_func_of_RL_out = outdir + "/corrhist_rc_func_of_RL.pdf"; // could add jetR and threshold info later??, maybe not needed tho 
    std::string fname_func_of_pT_out = outdir + "/corrhist_rc_func_of_pT.pdf";

    TLegend leg_ptRLbins(0.2, 0.2, 0.4, 0.45); 
    TLegend leg_ptbins(0.5, 0.3, 0.65, 0.45); 

    leg_ptRLbins.SetTextSize(0.037);
    leg_ptRLbins.SetBorderSize(0);
    leg_ptbins.SetTextSize(0.037);
    leg_ptbins.SetBorderSize(0);



    // get graphs of r_c as a function of RL
    // idkyet("pt", ptcenters_vec, pt_err, rc_vec, rc_errors_vec); // TODO: work on this! commented out so i could run code
            // get rid of "pt"??

    // get graphs of r_c as a function of ptRL
    // idkyet("ptrl", ptRLcenters_vec, ptRL_err, rc_vec, rc_errors_vec); // TODO: work on this! commented out so i could run code



    //----====----====----====----====----====----====----====----====----====----====----====----====
    // get graphs of r_c as a function of RL
    // loop over pt bins
    for ( int i = 0; i < ptcenters_vec.size(); i++ ) { 

        // for graphs as a function of ptRL
        TGraphErrors *g = new TGraphErrors(ptRLcenters_vec.size(), ptRLcenters_vec.data(), rc_vec[i].data(), ptRL_err.data(), rc_errors_vec[i].data());
        FormatGraphMarker(g, kBlack, 1.0, markers[i], marker_size);
        leg_ptbins.AddEntry(g, Form("p_{T} = %d-%d", (int)ptcenters_vec[i]-10, (int)ptcenters_vec[i]+10), "P");
        
        vector<TGraphErrors*> ind_temp_vec;
        for (int j=0; j<ptRLcenters_vec.size(); j++){
            int k=j+1;
            TGraphErrors *g_ind = new TGraphErrors(1, &ptRLcenters_vec[j], &rc_vec[i][j], &ptRL_err[j], &rc_errors_vec[i][j]);
            cout << "i: " << i << " j: " << j << " err: " << rc_errors_vec[i][j] << endl;
            FormatGraphMarker(g_ind, colors[k], 1.0, markers[i], marker_size);
            ind_temp_vec.push_back(g_ind);
            
            if (i==0) leg_ptRLbins.AddEntry(g_ind, Form("R_{L} bin %d", j+1), "P");
            
        }
        rc_graphs_func_of_ptRL.push_back(g); 
        rc_graphs_func_of_ptRL_ind.push_back(ind_temp_vec); 
                 

    }

    // plot r_c as a function of ptRL
    func();
    TCanvas *can_func_of_ptRL = new TCanvas("can_func_of_ptRL", "can_func_of_ptRL", 750, 500);
    can_func_of_ptRL->cd();
    for ( int i = 0; i < ptcenters_vec.size(); i++ ) {
        if ( i == 0 ) {
            rc_graphs_func_of_ptRL[i]->SetMinimum(-0.5);  // Lower y limit
            rc_graphs_func_of_ptRL[i]->SetMaximum(0.025); 
            rc_graphs_func_of_ptRL[i]->GetXaxis()->SetTitle("R_{L} bin center"); 
            rc_graphs_func_of_ptRL[i]->GetYaxis()->SetTitle("r_{c}"); 
            rc_graphs_func_of_ptRL[i]->Draw("AP");
        } 

        for (int j = 0; j < ptRLcenters_vec.size(); j++) {
            rc_graphs_func_of_ptRL_ind[i][j]->Draw("P SAME");
        }
    }
    leg_ptRLbins.Draw("same");
    leg_ptbins.Draw("same");
    can_func_of_ptRL->SaveAs(fname_func_of_RL_out.c_str());
    delete can_func_of_ptRL;

    //========================================================

    // get graphs of r_c as a function of pT
    vector<vector<double>> rc_vals_func_of_pT;
    vector<vector<double>> rc_err_vals_func_of_pT;
    // loop over RL bins
    for ( int j = 0; j < ptRLcenters_vec.size(); j++ ) {
        vector<double> temp_vec;
        vector<TGraphErrors*> ind_temp_vec;
        vector<double> err_temp_vec;

        // save values into appropriate vectors
        for ( int i = 0; i < ptcenters_vec.size(); i++ ) { 
            temp_vec.push_back(rc_vec[i][j]);
            err_temp_vec.push_back(rc_errors_vec[i][j]);
            
            int k=j+1;
            TGraphErrors *g_ind = new TGraphErrors(1, &ptcenters_vec[i], &rc_vec[i][j], &pt_err[i], &rc_errors_vec[i][j]);
            cout << "i: " << i << " j: " << j << " err: " << rc_errors_vec[i][j] << endl;
            FormatGraphMarker(g_ind, colors[k], 1.0, markers[i], marker_size);
            ind_temp_vec.push_back(g_ind);
        }
        rc_vals_func_of_pT.push_back(temp_vec);
        rc_err_vals_func_of_pT.push_back(err_temp_vec);
        rc_graphs_func_of_pT_ind.push_back(ind_temp_vec); 

        // for graphs as a function of RL
        TGraphErrors *g = new TGraphErrors(ptcenters_vec.size(), ptcenters_vec.data(), rc_vals_func_of_pT[j].data(), pt_err.data(), rc_err_vals_func_of_pT[j].data());
        for (int aa = 0; aa < ptcenters_vec.size(); aa++) {
            // cout << "studying pt=" << ptcenters_vec[aa] << " // " << rc_vals_func_of_pT[j][aa] << endl;
        }
        rc_graphs_func_of_pT.push_back(g); 
    }
    
    // plot r_c as a function of pT
    TCanvas *can_func_of_pT = new TCanvas("can_func_of_pT", "can_func_of_pT", 750, 500);
    can_func_of_pT->cd();
    for ( int j = 0; j < ptRLcenters_vec.size(); j++ ) {
        int k = j+1;
        // rc_graphs_func_of_pT[j]->SetMarkerSize(1.0);
        // rc_graphs_func_of_pT[j]->SetMarkerStyle(markers[0]);
        rc_graphs_func_of_pT[j]->SetMarkerColorAlpha(colors[k], 0.0);
        if ( j == 0 ) {
            rc_graphs_func_of_pT[j]->SetMinimum(-0.5);  // Lower y limit
            rc_graphs_func_of_pT[j]->SetMaximum(0.025); 
            rc_graphs_func_of_pT[j]->GetXaxis()->SetTitle("p_{T} bin center"); 
            rc_graphs_func_of_pT[j]->GetYaxis()->SetTitle("r_{c}"); 
            rc_graphs_func_of_pT[j]->Draw("AP");
        } // else { 
        //     rc_graphs_func_of_pT[j]->Draw("P SAME");
        // }

        for (int i = 0; i < ptcenters_vec.size(); i++) {
            rc_graphs_func_of_pT_ind[j][i]->Draw("P SAME");
        }
    }
    leg_ptRLbins.Draw("same");
    leg_ptbins.Draw("same");
    can_func_of_pT->SaveAs(fname_func_of_pT_out.c_str());
    delete can_func_of_pT;

}

*/

void get_deltajt_scatter(TChain * PAIRINFO_tree, vector<double>& ptrl_vals, vector<double>& deltajt_vals,
                        int pt_min, int pt_max){
    
    Float_t jetpt; Float_t rl;
    Float_t delta_jt;
    PAIRINFO_tree->SetBranchAddress("jet_pt", &jetpt);
    PAIRINFO_tree->SetBranchAddress("RL", &rl); // is really ptrl but labeled wrong
    PAIRINFO_tree->SetBranchAddress("deltajt", &delta_jt);

    int entries = PAIRINFO_tree->GetEntries();
    for ( int i = 0; i < entries; i++ ) {
        PAIRINFO_tree->GetEntry(i);
        if (jetpt >= pt_min && jetpt < pt_max) {
            ptrl_vals.push_back(jetpt*rl);
            deltajt_vals.push_back(delta_jt);
        }
    } 
}

// ======================================================= //
//                     SOME FUNCTIONS
// ======================================================= //

// need an address on Observable so that the original is modified, not a copy
void analyze_1D_obs(TChain * PAIRINFO_tree, TH1D * jetpt_inptbin_hist, TLegend& leg, int j, int k,
                    Observable& obs, int pt_min, int pt_max, double ptRL_min, double ptRL_max, int pt_avg, 
                    bool weighting, bool logbins, std::string norm_string, std::string ytitle_norm,
                    std::string ytitle_weight_str, std::string ptRLname_leg, std::string hist_addname,
                    std::string ptname) {

    cout << " in analyze_1d_obs! with " << obs.name << endl;

    // get histograms
    TH1D * obs_hist = getObs1DHistFromTChain(PAIRINFO_tree, obs, pt_min, pt_max, ptRL_min, ptRL_max, pt_avg, weighting, logbins);
    
    // push to vectors
    obs.addHist((TH1D*) obs_hist->Clone(obs_hist->GetName()));

    // format histograms in vector
    Format1DHist(obs, obs.obs_vec[k], jetpt_inptbin_hist, norm_string, colors[j], 0.6, markers[0], obs.axis_label, ytitle_norm + obs.cs_label + ytitle_weight_str, leg, ptRLname_leg, hist_addname);
    
    // draw, save, and delete histograms
    draw_save_del_hists(obs, obs.obs_vec[k], ptname, norm_string, hist_addname, false, true);
    delete obs_hist;
    

}

void analyze_2D_obs(TChain * PAIRINFO_tree, TH1D * jetpt_inptbin_hist, Observable2D obs2D, int pt_min, int pt_max, 
                    double ptRL_min, double ptRL_max, int pt_avg, bool weighting, bool logbins, std::string norm_string, 
                    std::string ytitle_norm, std::string ytitle_weight_str, std::string hist_addname, std::string ptname) {

    Observable obs_x = obs2D.obsx;
    Observable obs_y = obs2D.obsy;
    cout << "in analyze_2d_obs! with " << obs_y.name << "_vs_" << obs_x.name << endl;

    // get histogram
    TH2D * hist2D = getObs2DHistFromTChain(PAIRINFO_tree, obs_x, obs_y, pt_min, pt_max, ptRL_min, ptRL_max, pt_avg);
    
    // format histograms
    Format2DHist(obs_x, obs_y, hist2D, jetpt_inptbin_hist, norm_string, ytitle_norm + obs_x.axis_label, ytitle_norm + obs_y.axis_label, hist_addname, true);
    
    // plot and save
    if (obs_x.name == "ptrl") draw_save_del_hists2D(obs2D, hist2D, ptname, norm_string, hist_addname, true, false, true);
    else draw_save_del_hists2D(obs2D, hist2D, ptname, norm_string, hist_addname, false, false, true);
    delete hist2D;

} 


// TODO: implement weight in the charge!
void analyze_ptbin(TChain * JETINFO_tree, TChain * PAIRINFO_tree, vector<Observable> obs_1D_list, vector<Observable2D> obs_2D_list,
             std::string weightstr, std::string jetRname, std::string thrname,
             std::string norm_string, int pt_min, int pt_max, double pt_avg, const double ptRL_bins[], int n_ptRLbins, 
             vector<double> deltajt_vals, bool include_RL0, bool include_RL1, bool debug, bool debug2) {
    
    std::string ptname = to_string(pt_min) + "-" + to_string(pt_max);

    bool weighting = (weightstr == "_Weighted");
    std::string ytitle_weight_str = (weighting) ? " #times #frac{p_{T,1}p_{T,2}}{p_{T,jet}^{2}}" : "";

    std::string ytitle_norm = "";
    if (norm_string == "self_normalized") ytitle_norm = "#frac{1}{N_{pair}} ";
    if (norm_string == "norm_by_jets") ytitle_norm = "#frac{1}{N_{jet}} ";
    

    TCanvas *cdumdum = new TCanvas(); // need a canvas so legend can be made
    TLegend* leg = new TLegend(0.5, 0.62, 0.85, 0.85);
    leg->SetTextSize(0.037);
    leg->SetBorderSize(0);
    leg->AddEntry("NULL",Form("%d #leq p_{T, jet} < %d", pt_min, pt_max),"h");
    
    for ( int j = 0; j < n_ptRLbins; j++ ) {
        int k = j;
        if (!include_RL0) {
            k = j-1;
            if (j == 0) continue; // can add something here to change the filename for ALL
        }
        if (!include_RL1 && j == n_ptRLbins-1) continue;
        
        double ptRL_min = ptRL_bins[j];
        double ptRL_max = ptRL_bins[j+1];
        std::string ptRLname = Form("_pTRL%.1f-%.1f", ptRL_min, ptRL_max);
        std::string ptRLname_leg = Form("#LT p_{T} #GT R_{L} = %.1f-%.1f", ptRL_min, ptRL_max);
        std::string hist_addname = weightstr + jetRname + thrname + "_pt" + ptname + ptRLname + "_" + norm_string;
        if (debug) cout << " in ptRL bin" << j << " with " << ptRL_min << " - " << ptRL_max << endl;
        

        // get histograms
        TH1D * jetpt_inptbin_hist;
        Observable obs_jetpt("jet_pt", true, 200, 0, 200, "p_{T,jet}", "#frac{dN}{dp_{T,jet}}");
        if (norm_by_jets_bool) jetpt_inptbin_hist = getObs1DHistFromTChain(JETINFO_tree, obs_jetpt, pt_min, pt_max, 0, 0, pt_avg);

        for (Observable& obs : obs_1D_list) { // passing reference to not make a copy!
            if (obs.obs_bool) analyze_1D_obs(PAIRINFO_tree, jetpt_inptbin_hist, *leg, j, k, obs, 
                              pt_min, pt_max, ptRL_min, ptRL_max, pt_avg, weighting, logbins, norm_string, 
                              ytitle_norm, ytitle_weight_str, ptRLname_leg, hist_addname, ptname);  
        }
                    
        
        for (Observable2D obs2D : obs_2D_list) {
            if (obs2D.obs_bool) analyze_2D_obs(PAIRINFO_tree, jetpt_inptbin_hist, obs2D, pt_min, pt_max, ptRL_min, ptRL_max, pt_avg, weighting, logbins,
                                norm_string, ytitle_norm, ytitle_weight_str, hist_addname, ptname);
        }
        // // TH2D * zbig_vs_zsmol_hist2D = getObs2DHistFromTChain(PAIRINFO_tree, "zsmol", "zbig", z_numbins, 0, 1, z_numbins, 0, 1, pt_min, pt_max, ptRL_min, ptRL_max, pt_avg);
        // // TH2D * weights_vs_maxpt_hist2D = getObs2DHistFromTChain(PAIRINFO_tree, "maxpt", "weights", deltap_numbins, 0, 84, weights_numbins, 0, 0.3, pt_min, pt_max, ptRL_min, ptRL_max, pt_avg);

        // getting z?
        /*if (deltajt_bool && just_done_once == false) {
            TCanvas *can_deltajt_ptrl_scatter_ind = new TCanvas();
            ProcessCanvas(can_deltajt_ptrl_scatter_ind);

            TH1D * h_deltajt_ptrl_scatter_ind = new TH1D("h_deltajt_ptrl_scatter_ind", "h_deltajt_ptrl_scatter_ind", 200, 0., 1.);
            h_deltajt_ptrl_scatter_ind->GetXaxis()->SetTitle("#Deltaj_{T} / #LTp_{T} #GTR_{L}");
            for (int l = 0; l < ptrl_vals.size(); l++) {
                if (ptrl_vals[l] >= ptRL_min && ptrl_vals[l] < ptRL_max) {
                    h_deltajt_ptrl_scatter_ind->Fill(deltajt_vals[l] / ptrl_vals[l]);
                }
            }
            draw_save_del_hists(f_out, can_deltajt_ptrl_scatter_ind, h_deltajt_ptrl_scatter_ind, "deltajt_scatter", ptname, norm_string, hist_addname, false, false, "scatter_deltajt_ptrl");

        }*/
    }
    
    /* do pt bin stuff here */
    std::string hist_all_addname = weightstr + jetRname + thrname + "_pt" + ptname;

	// combine RL plots to get 1 plot per pt bin
    for (Observable& obs : obs_1D_list) { // passing reference to not make a copy!
        if (obs.obs_bool) plotandsave_combined_hists(obs, leg, ptname, norm_string, hist_all_addname, pt_max, logbins, true);
    }
    
}

void analyze_rc(TChain * JETINFO_tree, TChain * PAIRINFO_tree, std::string weightstr, std::string jetRname, 
                std::string thrname, const int pt_bins[], int n_bins, const double pt_avgs[], 
                const double ptRL_bins[], int n_ptRLbins, bool include_RL0, bool include_RL1, bool debug = false) {

    Observable obs_q1q2("q1q2", rc_bool, 6, -3, 3, "", "");
    Observable obs_rc("rc", rc_bool, 0, 0, 0, "", "");
    vector<double> ptcenters_vec;
    vector<double> ptRLcenters_vec;

    vector<vector<double>> rc_vec;
    vector<vector<double>> rc_errors_vec;
    
    bool weighting = (weightstr == "_Weighted");
    std::string ytitle_weight_str = (weighting) ? " #times #frac{p_{T,1}p_{T,2}}{p_{T,jet}^{2}}" : "";
    
    for ( int i = 0; i < n_bins; i++ ) {
        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];
        ptcenters_vec.push_back( (pt_min+pt_max)/2 );
        std::string ptname = to_string(pt_min) + "-" + to_string(pt_max);

        vector<double> rc_temp_vec;
        vector<double> rc_err_temp_vec;
        
        for ( int j = 0; j < n_ptRLbins; j++ ) {
            int k = j;
            if (!include_RL0) {
                k = j-1;
                if (j == 0) continue; // can add something here to change the filename for ALL
            }
            if (!include_RL1 && j == n_ptRLbins-1) continue;
            
            double ptRL_min = ptRL_bins[j];
            double ptRL_max = ptRL_bins[j+1];
            if (i==0) ptRLcenters_vec.push_back( (ptRL_min+ptRL_max)/2 );
            // std::string ptRLname = Form("_pTRL%.1f-%.1f", ptRL_min, ptRL_max);
            // std::string ptRLname_leg = Form("#LT p_{T} #GT R_{L} = %.1f-%.1f", ptRL_min, ptRL_max);
            // std::string hist_addname = weightstr + jetRname + thrname + "_pt" + ptname + ptRLname + "_" + norm_string;
            if (debug) cout << " in ptRL bin" << j << " with " << ptRL_min << " - " << ptRL_max << endl;
            
            vector<double> rc_and_err = getRcFromTChain(PAIRINFO_tree, obs_q1q2, pt_min, pt_max, ptRL_min, ptRL_max, pt_avgs[i], weighting);
            double rc_value = rc_and_err[0];
            double rc_err = rc_and_err[1];
            cout << "RC IS " << rc_value << "and RC ERR IS " << rc_err << "(pt_min=" << pt_min << ", j=" << j << ")" <<endl;

            rc_temp_vec.push_back(rc_value);
            rc_err_temp_vec.push_back(rc_err);
            
        }

        /* do pt bin stuff here */
        std::string hist_all_addname = weightstr + jetRname + thrname + "_pt" + ptname;

        rc_vec.push_back(rc_temp_vec);
        rc_errors_vec.push_back(rc_err_temp_vec);

        // make graphs
        // TCanvas *can_rc = new TCanvas();
        // ProcessCanvas(can_rc);
        TGraphErrors *gr_rc = MakeFormatGraph(ptRLcenters_vec, rc_vec[i], rc_errors_vec[i], kBlack, 1.0, markers[0], "R_{L}", "r_{c}", "rc", hist_all_addname);
        draw_save_del_hists(obs_rc, gr_rc, ptname, "unnormalized", hist_all_addname, false, false);

    }


    // plot r_c as a function of RL
    // cout << "checkpoint 4" << endl;
    // cout << "size of RL_vals " << RL_vals.size() << endl;
    // cout << "size of RL_vals[0] " << RL_vals[0].size() << endl;
    // cout << "size of rc_vec " << rc_vec.size() << endl;
    // cout << "size of rc_vec[0] " << rc_vec[0].size() << endl;
    // cout << "size of ptcenters_vec " << ptcenters_vec.size() << endl;
    // plot_rc(ptcenters_vec, ptRLcenters_vec, rc_vec, rc_errors_vec); //, leg_RLbins, leg_ptbins); // TODO: work on this! commented out so i could run code

}


// don't separate by pt or RL bin
void analyze_jetlevel_observables(TChain * JETINFO_tree, std::string jetRname, std::string thrname,
                                  bool debug, bool debug2) {

    Observable obs_jet_pt("jet_pt", true, 200, 0, 200, "p_{T,jet}", "");
    Observable obs_jet_const("total_num_const", true, 20, 0, 20, "Number Constituents (total)", "");
    Observable obs_jet_const_aftercut("num_const_aftercut", true, 20, 0, 20, "Number Constituents (after threshold cut)", "");

    obs_jet_pt.recreate_output_root_file();
    obs_jet_const.recreate_output_root_file();
    obs_jet_const_aftercut.recreate_output_root_file();

    TH1D * jetpt_hist = getObs1DHistFromTChain(JETINFO_tree, obs_jet_pt, 0, 200, 0, 0);
    jetpt_hist->GetXaxis()->SetTitle(obs_jet_pt.axis_label.c_str());
    draw_save_del_hists(obs_jet_pt, jetpt_hist, "", "", jetRname + thrname, false, true);

    TH1D * jet_const = getObs1DHistFromTChain(JETINFO_tree, obs_jet_const, 0, 200, 0, 0);
    jet_const->GetXaxis()->SetTitle(obs_jet_const.axis_label.c_str());
    draw_save_del_hists(obs_jet_const, jet_const, "", "", jetRname + thrname, false, true);

    TH1D * jet_const_aftercut = getObs1DHistFromTChain(JETINFO_tree, obs_jet_const_aftercut, 0, 200, 0, 0);
    jet_const_aftercut->GetXaxis()->SetTitle(obs_jet_const_aftercut.axis_label.c_str());
    draw_save_del_hists(obs_jet_const_aftercut, jet_const_aftercut, "", "", jetRname + thrname, false, true);

    // return;
}


void analyze(TChain * JETINFO_tree, TChain * PAIRINFO_tree, vector<Observable> obs_1D_list, vector<Observable2D> obs_2D_list,
             std::string weightstr, std::string jetRname, std::string thrname, std::string norm_string, 
             const double pt_avgs[], const int pt_bins[], int n_bins, const double ptRL_bins[], int n_ptRLbins,
             bool include_RL0, bool include_RL1, bool debug, bool debug2 ) {
    
    // needs to be separated by pt and RL bin
    for ( int i = 0; i < n_bins; i++ ) {
        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];
           
        if (debug) cout << " in pt bin" << i << " with " << pt_min << " - " << pt_max << endl;

        // make delta jt vs ptrl scatter plot??
        vector<double> ptrl_vals;
        vector<double> deltajt_vals;
        /*if (deltajt_bool && just_done_once == false) {
            get_deltajt_scatter(PAIRINFO_tree, ptrl_vals, deltajt_vals, pt_min, pt_max);
            cout << "NUM IN PTRL VALS 2: " << ptrl_vals.size() << endl;
            cout << "NUM IN DELTA JT VALS 2: " << deltajt_vals.size() << endl;
            TGraph * scatter_deltajt_vs_ptrl = new TGraph(ptrl_vals.size(), ptrl_vals.data(), deltajt_vals.data());

            
            TCanvas * can_scatter = new TCanvas();
            std::string ptname = to_string(pt_min) + "-" + to_string(pt_max);
            std::string hist_addname = "_pt" + ptname;
            draw_save_del_hists(f_out, can_scatter, scatter_deltajt_vs_ptrl, "deltajt_scatter", ptname, 
                                norm_string, hist_addname, true, false, false, "scatter_deltajt_vs_ptrl", ptRL_bins, n_ptRLbins, include_RL0, include_RL1);
        
                
        }*/

        // do the rest of the observables
        analyze_ptbin(JETINFO_tree, PAIRINFO_tree, obs_1D_list, obs_2D_list, weightstr, jetRname, thrname, norm_string, pt_min, pt_max, pt_avgs[i], ptRL_bins, n_ptRLbins, deltajt_vals, include_RL0, include_RL1, debug, debug2);
        
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
    attempt_dir = "data_fifthattempt_ptrlbins/";
    if (logbins == true) {
        attempt_dir += "logbins";
    } 
    outdir = "/software/users/blianggi/mypyjetty/storage/dEEC/plots/" + attempt_dir;
    
    // ntuple/histogram names
    std::string JETINFO_name = "tn_JETINFO_R0.4_1.0";
    std::string PAIRINFO_name = "tn_pairlevel_R0.4_1.0";
    std::string jet1D_name = "h_1Djet_pt_JetPt_R0.4_1.0"; // this one is a histogram
        
    // filenames
    std::string filename = Form("~/Documents/research/othercorrelations/data_ntuples/AnalysisResults_0001.root");
    std::string base_filepath_perly = Form("/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/31843529");
    std::string base_filepath_hic = Form("/rstorage/alice/AnalysisResults/blianggi/dEEC/468247"); //442528");
    

    Observable obs_deltap("deltap", deltap_bool, 42, 0, 84, "#Deltap", "#frac{dN}{d#Deltap}"); // bin sizes of 2 GeV
    Observable obs_deltapt("deltapt", deltapt_bool, 42, 0, 84, "#Deltap_{T}", "#frac{dN}{d#Deltap_{T}}");
    Observable obs_deltajt("deltajt", deltajt_bool, 50, 0, 5, "#Deltaj_{T}", "#frac{dN}{d#Deltaj_{T}}");
    Observable obs_weights("weights", ew_bool, 60, 0, 0.3, "#frac{p_{T,1}p_{T,2}}{p_{T,jet}^{2}}", "#frac{dN}{d[EW]}"); // bin sizes of 0.005
    Observable obs_z("z", twoDhists_bool, 50, 0, 1, "z", "#frac{dN}{dz}");

    Observable obs_p1("p1", p1_bool, 42, 0, 84, "p", "#frac{dN}{dp}");
    Observable obs_jt1("jt1", jt1_bool, 50, 0, 5, "j{T}", "#frac{dN}{dj_{T}}");

    Observable2D obs_weights_vs_deltap(obs_deltap, obs_weights, twoDhists_bool);
    Observable2D obs_weights_vs_deltajt(obs_deltajt, obs_weights, twoDhists_bool);
    Observable2D obs_weights_vs_p1(obs_p1, obs_weights, twoDhists_bool);
    Observable2D obs_zj_vs_zi(obs_z, obs_z, false);

    Observable obs_ptrl("ptrl", deltajt_vs_ptrl_bool, 50, 0.2, 35, "#LTp_{T}#GTR_{L}", "");
    Observable2D obs_deltajt_vs_ptrl(obs_ptrl, obs_deltajt, deltajt_vs_ptrl_bool);
            

    vector<Observable> obs_1D_list = { obs_deltap, obs_deltapt, obs_deltajt, obs_weights, obs_p1, obs_jt1 };
    vector<Observable2D> obs_2D_list = { obs_weights_vs_deltap, obs_weights_vs_deltajt, obs_weights_vs_p1, obs_zj_vs_zi, obs_deltajt_vs_ptrl };
    for (Observable obs : obs_1D_list) {
        if (obs.obs_bool) obs.recreate_output_root_file();
    }
    for (Observable2D obs2D : obs_2D_list) {
        if (obs2D.obs_bool) obs2D.recreate_output_root_file();
    }
    

    // analysis variables
    const int pt_bins[] = { 20, 40, 60, 80 };
    const double pt_avgs[] = { 25.0009, 46.5139, 67.3 };
    const int n_bins = sizeof(pt_bins) / sizeof(pt_bins[0]) - 1; //3;
    
    double ptRL_bins[7] = { 0, 2e-1, 8e-1, 5.0, 10.0, 30.0, 100.0 };
    int n_ptRLbins = sizeof(ptRL_bins) / sizeof(ptRL_bins[0]) - 1; //gets the columns //6; //7; //5;
    
    for (int a=0; a<7; a++) {
        cout << ptRL_bins[a] << " ";
    }
    cout << endl;
    if (debug2) cout << "pt_bins " << n_bins << " n_ptRLbins " << n_ptRLbins << endl;
    

    std::string jetRname = "_R0.4"; // + jetR;
    std::string thrname = "_t1.0"; // + threshold;
    std::string weightstr = ""; //"_xx";
    std::string norm_string = "unnormalized";
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

    // analyze jet level for plots
    // analyze_jetlevel_observables(JETINFO_tree, jetRname, thrname, debug, debug2);

    // rc analysis
    // analyze_rc(JETINFO_tree, PAIRINFO_tree, weightstr, jetRname, thrname, pt_bins, n_bins, pt_avgs, ptRL_bins, n_ptRLbins, include_RL0, include_RL1);

    // analyze for plots
    if (unweighted_bool) {
        if (unnormalized_bool) analyze(JETINFO_tree, PAIRINFO_tree, obs_1D_list, obs_2D_list, weightstr, jetRname, thrname, norm_string, pt_avgs, pt_bins, n_bins, ptRL_bins, n_ptRLbins, include_RL0, include_RL1, debug, debug2);
        
        norm_string = "self_normalized";
        if (self_normalized_bool) analyze(JETINFO_tree, PAIRINFO_tree, obs_1D_list, obs_2D_list, weightstr, jetRname, thrname, norm_string, pt_avgs, pt_bins, n_bins, ptRL_bins, n_ptRLbins, include_RL0, include_RL1, debug, debug2);

        norm_string = "norm_by_jets";
        if (norm_by_jets_bool) analyze(JETINFO_tree, PAIRINFO_tree, obs_1D_list, obs_2D_list, weightstr, jetRname, thrname, norm_string, pt_avgs, pt_bins, n_bins, ptRL_bins, n_ptRLbins, include_RL0, include_RL1, debug, debug2);
    }

    if (weighted_bool) {
        weightstr = "_Weighted";
        norm_string = "unnormalized";
        if (unnormalized_bool) analyze(JETINFO_tree, PAIRINFO_tree, obs_1D_list, obs_2D_list, weightstr, jetRname, thrname, norm_string, pt_avgs, pt_bins, n_bins, ptRL_bins, n_ptRLbins, include_RL0, include_RL1, debug, debug2);
        
        norm_string = "self_normalized";
        if (self_normalized_bool) analyze(JETINFO_tree, PAIRINFO_tree, obs_1D_list, obs_2D_list, weightstr, jetRname, thrname, norm_string, pt_avgs, pt_bins, n_bins, ptRL_bins, n_ptRLbins, include_RL0, include_RL1, debug, debug2);

        norm_string = "norm_by_jets";
        if (norm_by_jets_bool) analyze(JETINFO_tree, PAIRINFO_tree, obs_1D_list, obs_2D_list, weightstr, jetRname, thrname, norm_string, pt_avgs, pt_bins, n_bins, ptRL_bins, n_ptRLbins, include_RL0, include_RL1, debug, debug2);
    }


    // delete objects after saving for new pt-hat bin
    delete JETINFO_tree;
    delete PAIRINFO_tree;    
  
}


