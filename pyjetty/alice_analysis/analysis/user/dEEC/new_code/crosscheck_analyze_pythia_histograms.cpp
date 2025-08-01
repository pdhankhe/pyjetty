// ROOT macro make as a crosscheck
// it is currently able to plot 5 TeV pythia histograms (old method)
// Beatrice Liang-Gilman (beatrice_lg@berkeley.edu)

#include <iostream>
#include <yaml-cpp/yaml.h>


// global variables
// Double_t colors[16] = {kGray, kMagenta, kGreen+2, kBlue, kOrange+1, kViolet+1, kRed, kYellow+1, kCyan+1};
Double_t colors[16] = {kGray, kMagenta, kBlue, kOrange+1, kViolet+1, kGreen+2, kRed, kYellow+1, kCyan+1};
Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
Double_t marker_size = 1.5;

std::string attempt_dir = Form("pythia5TeV_histograms_crosscheck"); //;
std::string outdir = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/plots/" + attempt_dir; //;

bool write_to_root_file = true; 

// bool jetpt_bool = false;
bool deltap_bool = false;
bool deltapt_bool = false;
bool deltajt_bool = false;
bool ew_bool = false;
bool twoDhists_bool = false;
bool rc_bool = true; 

bool p1_bool = false;
bool jt1_bool = false;

bool unnormalized_bool = false;
bool self_normalized_bool = true;
bool norm_by_jets_bool = false;

bool unweighted_bool = true;
bool weighted_bool = true;

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


void Format1DHist(TH1D *hist, TH1D *jetpt_hist, std::string norm_string, double x_left, double x_right,
                  int markercolor, double markeralpha, int markerstyle, std::string xtitle, std::string ytitle, 
                  TLegend& leg, TString leg_text, 
                  bool drawline=false, double linealpha=1., std::string obs_name="",
                  std::string hist_addname = "") {

    std::string new_name = Form("h_%s%s", obs_name.c_str(), hist_addname.c_str());        
    // ^ = Form("h_%s_R0.4_t1.0_pt%d-%d_RL%.3f-%.3f_%s", obsname.c_str(), pt_min, pt_max, RL_min, RL_max, norm.c_str());
    hist->SetNameTitle(new_name.c_str(), new_name.c_str());

    // set x range
    // hist->GetXaxis()->SetRangeUser(x_left, x_right); //comment out for now... //TODO: see if can be fully removed
    
    // normalization
    if ( norm_string == "self_normalized" ) {
        double selfnorm_value = hist->Integral();
        hist->Scale(1/selfnorm_value, "width");
    } else if ( norm_string == "norm_by_jets" ) { // this is pretty much unused
        double numjets = jetpt_hist->Integral();
        cout << "Number of jets in " << leg_text << ": " << numjets << endl;
        hist->Scale(1/numjets, "width");
    }

    // stylization
    hist->SetLineColorAlpha(markercolor, markeralpha);
    hist->SetMarkerColorAlpha(markercolor, markeralpha);
    hist->SetMarkerStyle(markerstyle);
    hist->SetMarkerSize(1.5);
    if (drawline) {

        for (int a=0; a < hist->GetNbinsX(); a++){
            hist->SetBinError(a+1, 0);
        }

        hist->SetMarkerStyle(20);
        hist->SetMarkerColorAlpha(markercolor, 0);

        hist->SetFillStyle(0);
        hist->SetLineColorAlpha(markercolor, linealpha);
        hist->SetFillColor(markercolor);
        hist->SetLineStyle(1);
        hist->SetLineWidth(3);
    }

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


void Format2DHist(TH2D *hist2D, TH1D *jetpt_hist, std::string norm_string,
                  std::string xtitle, std::string ytitle, std::string obs_name="", std::string hist_addname = "") {

    std::string new_name = Form("h_%s%s", obs_name.c_str(), hist_addname.c_str());        
    // ^ = Form("h_%s_R0.4_t1.0_pt%d-%d_RL%.3f-%.3f_%s", obsname.c_str(), pt_min, pt_max, RL_min, RL_max, norm.c_str());
    hist2D->SetNameTitle(new_name.c_str(), new_name.c_str());

    // set x and y ranges //comment out for now... //TODO: see if can be fully removed
    // hist2D->GetXaxis()->SetRangeUser(0, ptmax);
    // hist2D->GetZaxis()->SetRangeUser(bounds[0], bounds[1]);
    
    // normalization
    if ( norm_string == "self_normalized" ) {
        double selfnorm_value = hist2D->Integral();
        hist2D->Scale(1/selfnorm_value, "width");
    } else if ( norm_string == "norm_by_jets" ) { // this is pretty much unused
        double numjets = jetpt_hist->Integral();
        cout << "Number of jets in " << obs_name << ": " << numjets << endl;
        hist2D->Scale(1/numjets, "width");
    }

    // axes
    hist2D->GetXaxis()->SetLabelFont(42);
    hist2D->GetXaxis()->SetTitleFont(42);
	hist2D->GetXaxis()->SetTitleSize(0.06); //(0.042);
    hist2D->GetXaxis()->SetTitleOffset(1.0); 
	hist2D->GetXaxis()->SetLabelSize(0.05);
    hist2D->GetXaxis()->SetTitle(xtitle.c_str());

    hist2D->GetYaxis()->SetLabelFont(42);
	hist2D->GetYaxis()->SetTitleFont(42);
    if (obs_name.find("weights") != std::string::npos) { // if weights is found in the string
        hist2D->GetYaxis()->SetTitleSize(0.035);
        hist2D->GetYaxis()->SetTitleOffset(1.5);
    } else {
        hist2D->GetYaxis()->SetTitleSize(0.05); //(0.042);
        hist2D->GetYaxis()->SetTitleOffset(1.0);
    }	
	hist2D->GetYaxis()->SetLabelSize(0.05); //(0.042);
    hist2D->GetYaxis()->SetTitle(ytitle.c_str());

}

void FormatGraphMarker(TGraphErrors * g, int markercolor, double markeralpha, int markerstyle, double markersize) {
    g->SetMarkerColorAlpha(markercolor, markeralpha);
    g->SetMarkerStyle(markerstyle);
    g->SetMarkerSize(markersize);
    g->SetLineColorAlpha(markercolor, markeralpha);
}

// imported function from data
TGraphErrors * MakeFormatGraph(vector<double> xvals, vector<double> yvals, int markercolor, double markeralpha,
                  int markerstyle, std::string xtitle, std::string ytitle, std::string obs_name) {
    
                    TGraphErrors * graph = new TGraphErrors(xvals.size(), xvals.data(), yvals.data());
    graph->SetTitle(Form("Charge Ratio;%s;%s", xtitle.c_str(), ytitle.c_str())); // Set the title and axis labels

    // Set graph styles
    FormatGraphMarker(graph, markercolor, markeralpha, markerstyle, marker_size);

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


// find the full width at half max
// returns a vector that gives (halfmax_height, leftPos, rightPos, full_width)
vector<double> findWidthOfCurve(TH1* hist, double topofCurve_bin) {

    vector<double> FWHM_vec;
    double topofCurve = hist->GetBinContent(topofCurve_bin);
    double halfmax_height = topofCurve/2;
    double leftPos = 0;
    double rightPos = 0;
    double full_width = 0;

    //find the left point of the half max width
    for (int i=topofCurve_bin-1; i>0.; i--) {
        double binheight = hist->GetBinContent(i);
        if (binheight <= halfmax_height) {
            leftPos = hist->GetBinCenter(i);
            break;
        }
    }
    //find the right point of the half max width
    for (int i=topofCurve_bin+1; i<hist->GetNbinsX(); i++) {
        double binheight = hist->GetBinContent(i);
        if (binheight <= halfmax_height) {
            rightPos = hist->GetBinCenter(i);
            break;
        }
    }

    full_width = rightPos - leftPos;

    FWHM_vec.push_back(halfmax_height);
    FWHM_vec.push_back(leftPos);
    FWHM_vec.push_back(rightPos);
    FWHM_vec.push_back(full_width);
    
    return FWHM_vec;
    
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
 

// Function to get the bin edges of a histogram
std::vector<double> get_bin_edges(TH1D * hist) {

    // Get the X-axis
    TAxis *xAxis = hist->GetXaxis();

    // Get the number of bins, minimum, and maximum values
    int nBins = xAxis->GetNbins();
    double xmin = xAxis->GetXmin();
    double xmax = xAxis->GetXmax();

    // Calculate the bin width
    double binWidth = (xmax - xmin) / nBins;

    // Generate the array of bin edges
    std::vector<double> binEdges;
    for (int i = 0; i <= nBins; ++i) {
        binEdges.push_back(xmin + i * binWidth);
    }
    
    return binEdges;
}




// ======================================================= //
//                     SOME FUNCTIONS 
// ======================================================= //

//get histogram and clone it
TH3D * get3DHistAndClone(TFile *f, std::string histname) {
    TH3D *h3D = (TH3D*) f->Get(histname.c_str());
    cout << "HNAME: " << histname.c_str() << endl;
    std::string hn = h3D->GetName();
    hn += "_clone";
    TH3D *h3D_clone = (TH3D *) h3D->Clone( hn.c_str() ); //TODO: change this name!

    return h3D_clone;
}


//get histogram and clone it
THnSparse * getTHnSparseAndClone(TFile *f, std::string histname) {
    THnSparse *hnsparse = (THnSparse*) f->Get(histname.c_str());
    std::string hn = hnsparse->GetName();
    hn += "_clone";
    THnSparse *hnsparse_clone = (THnSparse *) hnsparse->Clone( hn.c_str() ); //TODO: change this name!

    return hnsparse_clone;
}



/* get just 1D histogram */
// takes in a 1D histogram
TH1D * getObs1DHist_direct(TFile *filename, std::string h_name, bool debug=false) {
    TH1D *hist = (TH1D*) filename->Get(h_name.c_str());
    return hist;
}

/* get a typical 1D histogram */
// takes in a 3D histogram w/ 3 axes: (obs, pTRL, jet pT)
TH1D * getObs1DHistFrom3D(TFile *filename, std::string h_name,
                    int pt_min, int pt_max, double pTRL_min, double pTRL_max, bool debug=false) {

    TH3D *h3d = get3DHistAndClone(filename, h_name);
    h3d->GetYaxis()->SetRangeUser(pTRL_min, pTRL_max);
    h3d->GetZaxis()->SetRangeUser(pt_min, pt_max);
    TH1D *hist1D = h3d->ProjectionX();

    return hist1D;
}

/* get a typical 2D histogram */
// takes in a THnSparse w/ 4 axes: (x_obs, y_obs, pTRL, jet pT)
TH2D * getObs2DHist(TFile *filename, std::string h_name,
                    int pt_min, int pt_max, double pTRL_min, double pTRL_max, bool debug=false) {

    THnSparse *hnsparse = getTHnSparseAndClone(filename, h_name);
    hnsparse->GetAxis(2)->SetRangeUser(pTRL_min, pTRL_max);
    hnsparse->GetAxis(3)->SetRangeUser(pt_min, pt_max);
    TH2D *hist2D = hnsparse->Projection(1, 0);

    return hist2D;
}





// make a clone for the fit
TH1D * makeFitForClone(TH1D * hist, double xbound_min, double xbound_max) {
    TH1D * hclone = (TH1D*) hist->Clone(hist->GetName());
    hclone->GetXaxis()->SetRangeUser(xbound_min, xbound_max);
    return hclone;
}


// make a fit function - linear function
double linFunc(double *x, double *par) {
    return par[0] * x[0] + par[1]; // Example of a linear function y = par[0] * x + par[1]
}
// make a fit function - exponential function
double expFunc(double *x, double *par) {
    return par[0] * TMath::Exp(par[1] * x[0] + par[2]);
}


// fit to histogram
void fit_hist(TCanvas *c, TH1D* hist, double fit_xmin, double fit_xmax) {

    cout << "CHECKPOINT 1A" << endl;

    // TF1 *expfit = new TF1("expfit", "[0] * exp([1] * x + [2])", 0, 40); // 2 parameters

    // TF1 * expo_func = new TF1("expo","expo", fit_xmin, fit_xmax);
    TF1 * pol3_func = new TF1("pol3_func","[0] + [1]*x + [2]*pow(x,2) + [3]*pow(x,3)", 0,5);


    cout << "CHECKPOINT 1B" << endl;

    // Define the fit range
    // double fit_min = 0.3;
    // double fit_max = 0.7;

    // expfit->SetParameters(0.001, -0.25, 13); // Initial guess for parameters
    pol3_func->SetParameters(1E-5, 1E-7, -2E-6, 2E-7);
    hist->Fit("pol3_func", "R", "WLM", fit_xmin, fit_xmax);

    cout << "CHECKPOINT 1C" << endl;

    double A_fit = pol3_func->GetParameter(0); // Get first parameter
    double B_fit = pol3_func->GetParameter(1); 
    double C_fit = pol3_func->GetParameter(2); 
    double D_fit = pol3_func->GetParameter(3); 
    // double p0_err = pol3_func->GetParError(0); // Get error of the first parameter
    cout << "CHECKPOINT 1D" << endl;
    // graph->Draw("AP"); // Draw the graph with axis and points
    c->cd();
    pol3_func->Draw("same"); // Draw the fit function on the same canvas
    cout << "A_fit " << A_fit << " and B_fit " << B_fit << " and C_fit " << C_fit << endl;
}



/* Save and delete histograms */
// (TFile *fout, TCanvas *can, TH1D* hist, std::string obs_name, 
//                          std::string ptname, std::string norm_string, std::string hist_addname,
//                          bool logx, bool logy)
void draw_save_del_hists(TFile *fout, TCanvas *can, TObject* obj, std::string obs_name, 
                         std::string ptname, std::string norm_string, std::string hist_addname,
                         bool logx, bool logy, bool logz=false) {
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
    } else if (TGraph* graph = dynamic_cast<TGraph*>(obj)) {
        graph->Draw("ALP");
    } else {
        cout << "Error: Unsupported object type. Only TH1, TGraph, TGraphErrors, and TH2 are supported." << endl;
    }

    fout->cd();
    // hist->Write();
    obj->Write(); //TODO: this might not be right! Might have to use the casted type

    

    // std::string outdir = "/software/users/blianggi/mypyjetty/storage/dEEC/plots/data_firstattempt"; // + ptbin_name + "/";//"plots/test/";
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
        // gPad->SetRightMargin(0.12);
        // if (logz) gPad->SetLogz();
        // can->SetFillColor(kWhite);
        // hist2D->Draw("COLZ");
    } else if (TH1* hist = dynamic_cast<TH1*>(obj)) {
        hist->Draw();
    } else if (TGraphErrors* graph = dynamic_cast<TGraphErrors*>(obj)) {
        graph->Draw("ALP");
    } else if (TGraph* graph = dynamic_cast<TGraph*>(obj)) { // this is specifically for delta jt scatter
        graph->Draw("ALP");
    } else {
        cout << "Error: Unsupported object type. Only TH1, TGraph, TGraphErrors, and TH2 are supported." << endl;
        return;
    }

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
                          int pt_max, bool logx, bool logy, double pl_axis_cut=-1, bool debug=false) {

    // go into canvas
    can_all->cd();
    if (logx) gPad->SetLogx();
    if (logy) gPad->SetLogy();

    // if momentum axis, adjust x bounds accordingly
    size_t length = h_vec.size();
    for (int j=0; j<length; j++) {
        // cout << j << ": " << pTRL_bin_width[j] << endl;
        // if (scalebyRLbinwidth) h_vec[j]->Scale(pTRL_bin_width[j]); // this needs to be done before normalization
        // if (mom_axis) {
        //     // h_vec[j]->Rebin(4);
        //     // h_vec[j]->GetXaxis()->SetRangeUser(0, pt_max+5);
        //     // cout << h_vec[j]->GetEntries() << endl;
        //     h_vec[j]->Scale(pTRL_bin_width[j]); // this needs to be done before normalization
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
//                   2ND MAIN FUNCTION 
// ======================================================= //

void analyze_ptbin(TFile * f_in, TFile * f_out, std::string weightstr, std::string jetR, std::string threshold, 
                   std::string norm_string, int pt_min, int pt_max, const double ptRL_bins[], int n_ptRLbins, 
                   vector<vector<double>>& pTRL_vals, vector<vector<double>>& rc_vals, vector<vector<double>>& rc_errors, 
                   bool include_RL0, bool include_RL1, bool debug, bool debug2) {

    std::string jetRname = Form("_R%s", jetR.c_str());
    std::string thrname = Form("_t%s", threshold.c_str());
    std::string ptname = to_string(pt_min) + "-" + to_string(pt_max);
    cout << "pt min " << pt_min << " and pt max " << pt_max << endl;
    
    for (int j = 0; j < n_ptRLbins; ++j) {
        cout << "RL j= " << j << " gives " << ptRL_bins[j] << endl;
    }
    
    double pTRL_bin_width[7] = {0}; 
    double pTRL_bin_centers[7] = {0};
    for (int j = 0; j < n_ptRLbins; ++j) {
        pTRL_bin_width[j] = ptRL_bins[j+1] - ptRL_bins[j];
        pTRL_bin_centers[j] = (ptRL_bins[j+1] + ptRL_bins[j])/2;
        // cout << "RL BIN WIDTH HERE" << pTRL_bin_width[i][j] << endl;
        // cout << " AND CENTERS " << pTRL_bin_centers[j] << endl;
    }

    std::string ytitle_norm = "";
    if (norm_string == "self_normalized") ytitle_norm = "#frac{1}{N_{pair}} ";
    else if (norm_string == "norm_by_jets") ytitle_norm = "#frac{1}{N_{jet}} ";
    
    vector<TH1D*> deltap_vec;
    vector<TH1D*> deltapt_vec;
    vector<TH1D*> deltajt_vec;
    vector<TH1D*> weights_vec;
    // vector<TH1D*> q1q2_vec;
    vector<double> rc_vec;
    vector<double> rc_err_vec;
    vector<double> pTRLcenters_vec;

    TLegend *leg = new TLegend(0.6, 0.6, 0.85, 0.87);
    leg->SetTextSize(0.037);
    leg->SetBorderSize(0);
    TLegend *leg_dummy = new TLegend();

    // Names of histograms in the file
    // std::string deltap_truth_name = Form("h_corr_deltap%s_JetPt_Truth_R%s_%sScaled", weightstr.c_str(), jetR.c_str(), threshold.c_str());
    std::string deltap_truth_name = Form("gen_deltap_unmatchedScaled");
    std::string deltapt_truth_name = Form("gen_deltapt_unmatchedScaled");
    std::string deltajt_truth_name = Form("gen_deltajt_unmatchedScaled");
    std::string weights_truth_name = Form("gen_weights_unmatchedScaled");
    std::string charge_truth_name = Form("gen_charge_unmatchedScaled");
        
    std::string jet_pt_truth_name = Form("h_1Djet_pt_JetPt_Truth_R%s_%sScaled", jetR.c_str(), threshold.c_str());
    // std::string deltaptvsEW_truth_name = Form("h_ptvsenergyweights_JetPt_Truth_R%s_%sScaled", jetR.c_str(), threshold.c_str()); //not implemented yet

    std::string weights_vs_deltap_truth_name = Form("h2D_weights_vs_deltap_Truth_R%s_%sScaled", jetR.c_str(), threshold.c_str());
    std::string weights_vs_p1_truth_name = Form("h2D_weights_vs_p1_Truth_R%s_%sScaled", jetR.c_str(), threshold.c_str());
    std::string weights_vs_deltajt_truth_name = Form("h2D_weights_vs_deltajt_Truth_R%s_%sScaled", jetR.c_str(), threshold.c_str());
    std::string zj_vs_zi_truth_name = Form("h2D_zj_vs_zi_Truth_R%s_%sScaled", jetR.c_str(), threshold.c_str());
    
    for ( int j = 0; j < n_ptRLbins; j++ ) {
        int k = j;
        if (!include_RL0) {
            k = j-1;
            if (j == 0) continue; // can add something here to change the filename for ALL
        }
        if (!include_RL1 && j == n_ptRLbins-1) continue;

        // int i = 0;
        // if (pt_min == 40) i = 1;
        // else if (pt_min == 60) i = 2;
        
        double pTRL_min = ptRL_bins[j];
        double pTRL_max = ptRL_bins[j+1];
        std::string pTRLname = Form("_pTRL%.3f-%.3f", pTRL_min, pTRL_max);
        std::string pTRLname_leg = Form("#LTpT#GTR_{L} = %.1f-%.1f", pTRL_min, pTRL_max);
        // if (debug) cout << " in pTRL bin" << j << " with " << pTRL_min << " - " << pTRL_max << endl;
        
        std::string hist_addname = weightstr + jetRname + thrname + "_pt" + ptname + pTRLname + "_" + norm_string;
        

        // get histograms
        // get jet pT range - no D0 reconstruction, so don't make D0 cuts
        TH1D *hcorr_jetpt_inptbin_hist = getObs1DHist_direct(f_in, jet_pt_truth_name);
        TH1D *hcorr_deltap_truth_hist = getObs1DHistFrom3D(f_in, deltap_truth_name, pt_min, pt_max, pTRL_min, pTRL_max);
        cout << "checkpoint 1 " << hcorr_deltap_truth_hist->GetEntries() << endl;
        TH1D *hcorr_deltapt_truth_hist = getObs1DHistFrom3D(f_in, deltapt_truth_name, pt_min, pt_max, pTRL_min, pTRL_max);
        TH1D *hcorr_deltajt_truth_hist = getObs1DHistFrom3D(f_in, deltajt_truth_name, pt_min, pt_max, pTRL_min, pTRL_max);
        TH1D *hcorr_weights_truth_hist = getObs1DHistFrom3D(f_in, weights_truth_name, pt_min, pt_max, pTRL_min, pTRL_max);

        
        double rc_value = 0.0;
        double rc_err = 0.0;
        
        // NOT IMPLEMENTED - Rc implementation
        // if (norm_string == "unnormalized") {
        //     TH1D *hcorr_charge_truth_hist = getObs1DHistFrom3D(f_in, charge_truth_name, i+1, 4, pt_min, pt_max, pTRL_min, pTRL_max);
        //     rc_value = getRcFromHists(hcorr_charge_truth_hist); //hcorr_oppcharge_truth, hcorr_samecharge_truth);
        //     rc_err = getRcErrFromHists(hcorr_charge_truth_hist); //PAIRINFO_tree, "rc", 6, -3, 3, pt_min, pt_max, pTRL_min, pTRL_max);
        //     cout << "RC ERR IS " << rc_err << "(pt_min=" << pt_min << ", j=" << j << ")" <<endl;

        //     rc_vec.push_back(rc_value);
        //     rc_err_vec.push_back(rc_err);
        //     pTRLcenters_vec.push_back( (pTRL_min+pTRL_max)/2 );
        // }

        // 2D histograms!
        TH2D * weights_vs_deltap_hist2D = getObs2DHist(f_in, weights_vs_deltap_truth_name, pt_min, pt_max, pTRL_min, pTRL_max);
        TH2D * weights_vs_p1_hist2D = getObs2DHist(f_in, weights_vs_p1_truth_name, pt_min, pt_max, pTRL_min, pTRL_max);
        TH2D * weights_vs_deltajt_hist2D = getObs2DHist(f_in, weights_vs_deltajt_truth_name, pt_min, pt_max, pTRL_min, pTRL_max);
        TH2D * zj_vs_zi_hist2D = getObs2DHist(f_in, zj_vs_zi_truth_name, pt_min, pt_max, pTRL_min, pTRL_max);


        // push to vectors
        deltap_vec.push_back((TH1D*) hcorr_deltap_truth_hist->Clone(hcorr_deltap_truth_hist->GetName()));
        deltapt_vec.push_back((TH1D*) hcorr_deltapt_truth_hist->Clone(hcorr_deltapt_truth_hist->GetName()));
        deltajt_vec.push_back((TH1D*) hcorr_deltajt_truth_hist->Clone(hcorr_deltajt_truth_hist->GetName()));
        weights_vec.push_back((TH1D*) hcorr_weights_truth_hist->Clone(hcorr_weights_truth_hist->GetName()));
        

        // format histograms in vector
        Format1DHist(deltap_vec[k], hcorr_jetpt_inptbin_hist, norm_string, 0, pt_max+5, colors[j], 0.6, markers[0], "#Deltap", ytitle_norm + "#frac{dN}{d#Deltap}", *leg, pTRLname_leg, true, 1.0, "deltap", hist_addname);
        Format1DHist(deltapt_vec[k], hcorr_jetpt_inptbin_hist, norm_string, 0, pt_max+5, colors[j], 0.6, markers[0], "#Deltap_{T}", ytitle_norm + "#frac{dN}{d#Deltap_{T}}", *leg_dummy, pTRLname_leg, true, 1.0, "deltapt", hist_addname);
        Format1DHist(deltajt_vec[k], hcorr_jetpt_inptbin_hist, norm_string, 0, pt_max/2, colors[j], 0.6, markers[0], "#Deltap_{L}", ytitle_norm + "#frac{dN}{d#Deltap_{L}}", *leg_dummy, pTRLname_leg, true, 1.0, "deltajt", hist_addname);
        Format1DHist(weights_vec[k], hcorr_jetpt_inptbin_hist, norm_string, 0, 0.3, colors[j], 0.6, markers[0], "#frac{p_{T,1}p_{T,2}}{p_{T,jet}^{2}}", ytitle_norm + "#frac{dN}{d[EW]}", *leg_dummy, pTRLname_leg, true, 1.0, "weights", hist_addname);
        
        Format2DHist(weights_vs_deltap_hist2D, hcorr_jetpt_inptbin_hist, norm_string, ytitle_norm + "#Deltap", ytitle_norm + "#frac{p_{T,1}p_{T,2}}{p_{T,jet}^{2}}", "weights_vs_deltap", hist_addname);
        Format2DHist(weights_vs_p1_hist2D, hcorr_jetpt_inptbin_hist, norm_string, ytitle_norm + "p_{1}", ytitle_norm + "#frac{p_{T,1}p_{T,2}}{p_{T,jet}^{2}}", "weights_vs_p1", hist_addname);
        Format2DHist(weights_vs_deltajt_hist2D, hcorr_jetpt_inptbin_hist, norm_string, ytitle_norm + "#Deltaj_{T}", ytitle_norm + "#frac{p_{T,1}p_{T,2}}{p_{T,jet}^{2}}", "weights_vs_deltajt", hist_addname);
        Format2DHist(zj_vs_zi_hist2D, hcorr_jetpt_inptbin_hist, norm_string, ytitle_norm + "z_{i}", ytitle_norm + "z_{j}", "zj_vs_zi", hist_addname);
        
        // Get the X-axis
        TAxis *xAxis = deltap_vec[k]->GetXaxis();

        // Get the number of bins, minimum, and maximum values
        int nBins = xAxis->GetNbins();
        double xmin = xAxis->GetXmin();
        double xmax = xAxis->GetXmax();

        // Calculate the bin width
        double binWidth = (xmax - xmin) / nBins;

        // Generate the array of bin edges
        std::vector<double> binEdges;
        std::cout << "NBINS: " << nBins << ", bin edges: " << endl;
        for (int i = 0; i <= nBins; ++i) {
            binEdges.push_back(xmin + i * binWidth);
            cout << binEdges[i] << " ";
        }
        cout << endl;

        f_out->cd();

        // draw, save, and delete histograms
        TCanvas *can_deltap = new TCanvas();
        TCanvas *can_deltapt = new TCanvas();
        TCanvas *can_deltajt = new TCanvas();
        TCanvas *can_weights = new TCanvas();

        draw_save_del_hists(f_out, can_deltap, deltap_vec[k], "deltap", ptname, norm_string, hist_addname, false, true);
        draw_save_del_hists(f_out, can_deltapt, deltapt_vec[k], "deltapt", ptname, norm_string, hist_addname, false, true);
        draw_save_del_hists(f_out, can_deltajt, deltajt_vec[k], "deltapl", ptname, norm_string, hist_addname, false, false);
        draw_save_del_hists(f_out, can_weights, weights_vec[k], "weights", ptname, norm_string, hist_addname, false, true);
        

        TCanvas *can_weights_vs_deltap = new TCanvas("can_weights_vs_deltap", "can_weights_vs_deltap", 800, 500);
        TCanvas *can_weights_vs_p1 = new TCanvas("can_weights_vs_p1", "can_weights_vs_p1", 800, 500);
        TCanvas *can_weights_vs_deltajt = new TCanvas("can_weights_vs_deltajt", "can_weights_vs_deltajt", 800, 500);
        TCanvas *can_zj_vs_zi = new TCanvas("can_zj_vs_zi", "can_zj_vs_zi", 800, 500);
 
        draw_save_del_hists(f_out, can_weights_vs_deltap, weights_vs_deltap_hist2D, "weights_vs_deltap", ptname, norm_string, hist_addname, false, false, true);
        draw_save_del_hists(f_out, can_weights_vs_p1, weights_vs_p1_hist2D, "weights_vs_p1", ptname, norm_string, hist_addname, false, false, true);
        draw_save_del_hists(f_out, can_weights_vs_deltajt, weights_vs_deltajt_hist2D, "weights_vs_deltajt", ptname, norm_string, hist_addname, false, false, true);
        draw_save_del_hists(f_out, can_zj_vs_zi, zj_vs_zi_hist2D, "zj_vs_zi", ptname, norm_string, hist_addname, false, false, true);
                
    }

    /* do pt bin stuff here */
    std::string hist_all_addname = weightstr + jetRname + thrname + "_pt" + ptname;

	// combine RL plots to get 1 plot per pt bin

    TCanvas *can_deltap_all = new TCanvas();
    TCanvas *can_deltapt_all = new TCanvas();
    TCanvas *can_deltajt_all = new TCanvas();
    TCanvas *can_weights_all = new TCanvas();

    // size_t length_deltap = deltap_vec.size();
    // cout << " LENGTH DELTA P " << length_deltap << endl;
	
    plotandsave_combined_hists(can_deltap_all, deltap_vec, leg, "deltap", ptname, norm_string, hist_all_addname, pt_max, false, true, -1);
    plotandsave_combined_hists(can_deltapt_all, deltapt_vec, leg, "deltapt", ptname, norm_string, hist_all_addname, pt_max, false, true, -1);
    plotandsave_combined_hists(can_deltajt_all, deltajt_vec, leg, "deltajt", ptname, norm_string, hist_all_addname, pt_max, false, true, -1);
    plotandsave_combined_hists(can_weights_all, weights_vec, leg, "weights", ptname, norm_string, hist_all_addname, pt_max, false, true, -1);

    // // make graphs
    // if (norm_string == "unnormalized") {
    //     TCanvas *can_rc = new TCanvas();
    //     ProcessCanvas(can_rc);
    //     TGraphErrors *gr_rc = MakeFormatGraph(pTRLcenters_vec, rc_vec, kBlack, 1.0, markers[0], "R_{L}", "r_{c}", "rc");
    //     draw_save_del_hists(f_out, can_rc, gr_rc, "rc", ptname, norm_string, hist_all_addname, false, false);
        
        
    //     // save vectors here
    //     pTRL_vals.push_back(pTRLcenters_vec);
    //     rc_vals.push_back(rc_vec);
    //     rc_errors.push_back(rc_err_vec);
    // }

}

// void plot_histograms(TFile* f, std::string add_name, int normed, bool weighted, bool include_RL0) {
void plot_histograms(TFile* f_in, TFile* f_out, bool weighted, std::string norm_string, 
                     const int pt_bins[], int n_bins, const double ptRL_bins[], int n_ptRLbins, 
                     bool include_RL0, bool include_RL1, bool debug, bool debug2) {
    
    // //    gROOT->SetBatch(); //prevents plots from showing up

 
    std::string weightstr = "";
    if (weighted == true) {
        weightstr = "_Weighted";
    }


    // std::string RL0_string = "";
    // if (include_RL0) RL0_string = "_withRL0";

    // add_name += norm_string + RL0_string;


    // lists
    std::string jetR_list[] = { "0.4" };
    std::string threshold_list[] = { "1.0" }; // "0.15", "0.5"

     // Jet r value
     for (std::string jetR : jetR_list) {
        std::string jetRname = Form("_R%s", jetR.c_str());
        
        for (std::string threshold : threshold_list) {
            std::string thrname = Form("_t%s", threshold.c_str());

            // Names of histograms in the file
            std::string jet_pt_truth_name = Form("h_1Djet_pt_JetPt_Truth_R%s_%sScaled", jetR.c_str(), threshold.c_str());
            std::string jet_num_const_truth_name = Form("h_Nconst_JetPt_Truth_R%s_%sScaled", jetR.c_str(), threshold.c_str());

            if (norm_string == "unnormalized") {
                TH1D * jetpt_hist = getObs1DHist_direct(f_in, jet_pt_truth_name);
                jetpt_hist->GetXaxis()->SetTitle("p_{T,jet}");
                TCanvas *can_jetpt = new TCanvas();
                draw_save_del_hists(f_out, can_jetpt, jetpt_hist, "jet_pt", "", "", weightstr + jetRname + thrname, false, true);
            
                TH1D * jet_const = getObs1DHist_direct(f_in, jet_num_const_truth_name);
                jet_const->GetXaxis()->SetTitle("Number Constituents"); // TODO: check if this is before or after the threshold cut!
                TCanvas *can_numconst = new TCanvas();
                draw_save_del_hists(f_out, can_numconst, jet_const, "total_num_const", "", "", weightstr + jetRname + thrname, false, true);
            
                // TH1D * jet_const_aftercut = getObs1DHist_direct(f_in, jet_pt_truth_name);
                // jet_const_aftercut->GetXaxis()->SetTitle("Number Constituents (after threshold cut)");
                // TCanvas *can_numconst_aftercut = new TCanvas();
                // draw_save_del_hists(f_out, can_numconst_aftercut, jet_const_aftercut, "num_const_aftercut", "", "", weightstr + jetRname + thrname, false, true);
            
            }

            //variables
            vector<vector<double>> pTRL_vals;
            vector<vector<double>> rc_vals;
            vector<double> ptcenter_bins;
            vector<vector<double>> rc_errors;
    
            for (int i = 0; i < n_bins; i++) {
                cout << "in pt bin" << i << endl;
                int pt_min = pt_bins[i];
                int pt_max = pt_bins[i+1];
                ptcenter_bins.push_back( (pt_min+pt_max)/2 );
                // std::string ptname = "_pt" + std::to_string(pt_min) + '-' + std::to_string(pt_max);

                if (debug) cout << " in pt bin" << i << " with " << pt_min << " - " << pt_max << endl;
                analyze_ptbin(f_in, f_out, weightstr, jetR, threshold, norm_string, pt_min, pt_max, ptRL_bins, n_ptRLbins, pTRL_vals, rc_vals, rc_errors, include_RL0, include_RL1, debug, debug2);
        

                
            } // pT bins loop

            // // plot r_c as a function of RL
            // if (norm_string == "unnormalized") {
            //     cout << "checkpoint 4" << endl;
            //     cout << "size of pTRL_vals " << pTRL_vals.size() << endl;
            //     cout << "size of pTRL_vals[0] " << pTRL_vals[0].size() << endl;
            //     cout << "size of rc_vals " << rc_vals.size() << endl;
            //     cout << "size of rc_vals[0] " << rc_vals[0].size() << endl;
            //     cout << "size of ptcenter_bins " << ptcenter_bins.size() << endl;

            //     plot_rc(pTRL_vals, rc_vals, ptcenter_bins, rc_errors); //, colors, markers); //, leg_RLbins, leg_ptbins);
            // }

        } // threshold loop
    } // jetR loop

    
}

// ======================================================= //
//                     MAIN FUNCTION 
// ======================================================= //

void crosscheck_analyze_pythia_histograms() {

    gStyle->SetOptStat(0);
    SetStyle();
    
    // CONTOL VARIABLES HERE
    // normed is 0 if unnormalized, 1 for self-normalization 
    // weighted is true if using "Weighted", false if using unweighted
    bool weighted = false;
    std::string norm_string = "";
    bool include_RL0 = false;
    bool include_RL1 = false;
    
    // setup variables
    bool debug = false;
    bool debug2 = false;

    // int filecounter = 0;
    // int filecounter_cutoff = 500; //total: 5000
    
    // Files
    // const char infile[] = "/rstorage/alice/AnalysisResults/blianggi/dEEC/445125/1132588/scaling/AnalysisResultsFinal.root"; //hiccup
    const char infile[] = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/40154566/1132588/scaling/AnalysisResultsFinal.root"; //perlmutter, after june 2024
    // const char infile[] = "/Volumes/NO NAME/AnalysisResultsFinal.root"; //local
    TFile* root_infile = new TFile(infile, "READ");

    // Output file for binned results
    // std::string outfile = "/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/" + attempt_dir + "/PYTHIAHists.root"; //plots/ntuples/DataHists.root"; //FinalDataHists.root
    std::string outfile = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/rootfiles/" + attempt_dir + "/PYTHIAHists.root"; // perlmutter
    TFile* root_outfile = new TFile(outfile.c_str(), "RECREATE");
    // std::string add_name = ""; // "_othercorrel";
    // cout << "output name will be " << add_name << endl;


    // analysis variables
    const int pt_bins[] = { 20, 40, 60, 80 };
    const int n_bins = sizeof(pt_bins) / sizeof(pt_bins[0]) - 1; //3;
    
    const double ptRL_bins[7] = { 0., 2e-1, 8e-1, 5.0, 10.0, 30.0, 100.0 };
    const int n_ptRLbins = sizeof(ptRL_bins) / sizeof(ptRL_bins[0]) - 1; //gets the columns //6; //7; //5;


    if (debug2) cout << "pt_bins " << n_bins << " n_ptRLbins " << n_ptRLbins << endl;

            
    // ====================================================================================


    // analyze for plots
    // norm_string = "unnormalized";
    // plot_histograms(root_infile, root_outfile, weighted, norm_string, pt_bins, n_bins, ptRL_bins, n_ptRLbins, include_RL0, include_RL1, debug, debug2);
    
    norm_string = "self_normalized";
    plot_histograms(root_infile, root_outfile, weighted, norm_string, pt_bins, n_bins, ptRL_bins, n_ptRLbins, include_RL0, include_RL1, debug, debug2);
        
    // norm_string = "norm_by_jets";
    // plot_histograms(root_infile, root_outfile, weighted, norm_string, pt_bins, n_bins, ptRL_bins, n_ptRLbins, include_RL0, include_RL1, debug, debug2);
    


    // // first do unnormalized
    // int normed = 0;
    // // plot unweighted hists
    // plot_histograms(f, add_name, normed, false, include_RL0);
    // // // plot weighted hists
    // // plot_histograms(f, add_name, normed, true, include_RL0);


    // // // now do self-normalized
    // // normed = 1;
    // // // plot unweighted hists
    // // plot_histograms(f, add_name, normed, false, include_RL0);
    // // // plot weighted hists
    // // plot_histograms(f, add_name, normed, true, include_RL0);


    // // // now do normalize by # of jets
    // // normed = 2;
    // // // plot unweighted hists
    // // plot_histograms(f, add_name, normed, false, include_RL0);
    // // // plot weighted hists
    // // plot_histograms(f, add_name, normed, true, include_RL0);



    root_infile->Close();
    delete root_infile;

    root_outfile->Close();


    return;
}
