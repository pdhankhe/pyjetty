// ROOT macro make as a crosscheck
// it is currently able to plot 5 TeV pythia histograms (old method), and also herwig
// To switch between pythia and herwig, change the generator/attempt_dir variables in Lines 18 + 19/20/21/22
// Beatrice Liang-Gilman (beatrice_lg@berkeley.edu)

// TODO: figure out how to make this possible with gen or reco, matched or unmatched

#include <iostream>
using namespace std;


// global variables
// Double_t colors[16] = {kGray, kMagenta, kGreen+2, kBlue, kOrange+1, kViolet+1, kRed, kYellow+1, kCyan+1};
Double_t colors[16] = {kGray, kMagenta, kBlue, kOrange+1, kViolet+1, kGreen+2, kRed, kYellow+1, kCyan+1};
Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
Double_t marker_size = 1.5;

std::string generator = "pythia_cteq"; // "pythia" or "herwig" or "anchMC" or "pythia_cteq"
// std::string attempt_dir = Form("pythia5TeV_histograms_crosscheck");
// std::string attempt_dir = Form("herwig_secondattempt");
// std::string attempt_dir = Form("anchMC_firstattempt");
std::string attempt_dir = Form("pythia_cteq_firstattempt");
std::string outdir = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/plots/" + attempt_dir; //; // PERLY_FIX: REVERT TO THIS
std::string outdir_rootfiles = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/rootfiles/" + attempt_dir; //; // PERLY_FIX: REVERT TO THIS
// std::string outdir = "~/Documents/research/code/dEEC/new_code/storage/plots/" + attempt_dir;
// std::string outdir_rootfiles = "~/Documents/research/code/dEEC/new_code/storage/rootfiles/" + attempt_dir;

bool write_to_root_file = true; 

// bool jetpt_bool = false;
bool deltap_bool = true;
bool deltajt_bool = true;
bool ew_bool = true;
bool twoDhists_bool = true;
bool rc_bool = false;

bool p1_bool = false;
bool jt1_bool = false; 
bool deltajt_vs_ptrl_bool = false; // this doesn't work yet bc histogram doesn't exist yet

bool unnormalized_bool = false;
bool self_normalized_bool = true;
bool norm_by_jets_bool = false;

bool unweighted_bool = true;
bool weighted_bool = true; 

class Observable {
public:
    std::string name;
    bool obs_bool;
    std::string hist_toextract_name; // this is the histogram taken from the AnalysisResults file, could be 1D, 3D, whatever

    int num_bins;
    double min_bound;
    double max_bound;
    
    std::string axis_label;
    std::string cs_label; //cross section label
    std::string filepath_plots;

    std::vector<TH1D*> obs_vec;

    Observable(std::string name_val, bool obs_bool_val, std::string hist_toextract_name_val, 
               int num_bins_val, double min_bound_val, double max_bound_val,
               std::string axis_label_val, std::string cs_label_val) {
        name = name_val;
        obs_bool = obs_bool_val;
        hist_toextract_name = hist_toextract_name_val;

        num_bins = num_bins_val;
        min_bound = min_bound_val;
        max_bound = max_bound_val;

        axis_label = axis_label_val;
        cs_label = cs_label_val; //cross section label, in y axis
        
        filepath_plots = outdir + "/%s/%s/" + name + "/%s"; // ptname, norm_string, filename
        if (name.find("jet_") != std::string::npos || name.find("const") != std::string::npos) filepath_plots = outdir + "/%s"; // filename

    }

    void addHist(TH1D* hist) {
        obs_vec.push_back(hist);
    }

    void recreate_output_root_file() {
        std::string root_outfile = outdir_rootfiles + "/PYTHIAHists_" + name + ".root";
        if (generator == "herwig") root_outfile = outdir_rootfiles + "/HERWIGHists_" + name + ".root";
        if (generator == "anchMC") root_outfile = outdir_rootfiles + "/ANCHMCHists_" + name + ".root";
        if (generator == "pythia_cteq") root_outfile = outdir_rootfiles + "/PYTHIACTEQHists_" + name + ".root";
        TFile * f_out = new TFile(root_outfile.c_str(), "RECREATE");
        f_out->Close();
    }

    TFile * get_output_root_file() {
        std::string root_outfile = outdir_rootfiles + "/PYTHIAHists_" + name + ".root";
        if (generator == "herwig") root_outfile = outdir_rootfiles + "/HERWIGHists_" + name + ".root";
        if (generator == "anchMC") root_outfile = outdir_rootfiles + "/ANCHMCHists_" + name + ".root";
        if (generator == "pythia_cteq") root_outfile = outdir_rootfiles + "/PYTHIACTEQHists_" + name + ".root";
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

    Observable2D(Observable ox, Observable oy, bool obs_bool_val)
        : obsx(ox), obsy(oy), obs_bool(obs_bool_val) // required for non-default-constructible members
    {
        name = obsy.name + "_vs_" + obsx.name;
        filepath_plots = outdir + "/%s/%s/" + name + "/%s"; // ptname, norm_string, filename
    }

    void recreate_output_root_file() {
        std::string root_outfile = outdir_rootfiles + "/PYTHIAHists_" + name + ".root";
        if (generator == "herwig") root_outfile = outdir_rootfiles + "/HERWIGHists_" + name + ".root";
        if (generator == "anchMC") root_outfile = outdir_rootfiles + "/ANCHMCHists_" + name + ".root";
        if (generator == "pythia_cteq") root_outfile = outdir_rootfiles + "/PYTHIACTEQHists_" + name + ".root";
        TFile * f_out = new TFile(root_outfile.c_str(), "RECREATE");
        f_out->Close();
    }

    TFile* get_output_root_file() {
        std::string root_outfile = outdir_rootfiles + "/PYTHIAHists_" + name + ".root";
        if (generator == "herwig") root_outfile = outdir_rootfiles + "/HERWIGHists_" + name + ".root";
        if (generator == "anchMC") root_outfile = outdir_rootfiles + "/ANCHMCHists_" + name + ".root";
        if (generator == "pythia_cteq") root_outfile = outdir_rootfiles + "/PYTHIACTEQHists_" + name + ".root";
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


/* Format and adjust histograms */
void Format1DHist(Observable obs, TH1D *hist, TH1D *jetpt_hist, std::string norm_string, int markercolor, double markeralpha,
                  int markerstyle, std::string xtitle, std::string ytitle, TLegend& leg, TString leg_text, bool drawline=false,
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
    if (drawline) {
        for (int a=0; a < hist->GetNbinsX(); a++){
            hist->SetBinError(a+1, 0);
        }

//        hist->SetMarkerStyle(20);
//        hist->SetMarkerColorAlpha(markercolor, 0);

        hist->SetFillStyle(0);
//        hist->SetFillColor(markercolor);
        hist->SetLineStyle(1);
        hist->SetLineWidth(3);
    }

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
TH3D * get3DHistAndClone(TFile *f, std::string histname, int pt_min, int pt_max, int pTRL_bin) {
    TH3D *h3D = (TH3D*) f->Get(histname.c_str());
    cout << "HNAME: " << histname.c_str() << endl;
    // std::string hn = Form("%s_clone_pt%d-%d_ptRLbin%d", h3D->GetName(), pt_min, pt_max, ptRL_bin);
    std::string hn = Form("%s_clone", h3D->GetName());
    TH3D *h3D_clone = (TH3D *) h3D->Clone( Form("%s_pt%d-%d_ptRLbin%d",hn.c_str(), pt_min, pt_max, pTRL_bin) );

    return h3D_clone;
}


//get histogram and clone it
THnSparse * getTHnSparseAndClone(TFile *f, std::string histname, int pt_min, int pt_max, int pTRL_bin) {
    THnSparse *hnsparse = (THnSparse*) f->Get(histname.c_str());
    std::string hn = hnsparse->GetName();
    hn += "_clone";
    THnSparse *hnsparse_clone = (THnSparse *) hnsparse->Clone( Form("%s_pt%d-%d_ptRLbin%d",hn.c_str(), pt_min, pt_max, pTRL_bin) );

    return hnsparse_clone;
}



/* get just 1D histogram */
// takes in a 1D histogram
TH1D * getObs1DHist_direct(TFile *filename, std::string h_name, bool apply_cuts = false, int x_min = 0, int x_max = 1000, bool debug=false) {
    TH1D *hist = (TH1D*) filename->Get(h_name.c_str());
    if (apply_cuts) hist->GetXaxis()->SetRangeUser(x_min, x_max);
    return hist;
}

/* get a typical 1D histogram */
// takes in a 3D histogram w/ 3 axes: (obs, pTRL, jet pT)
TH1D * getObs1DHistFrom3D(TFile *filename, Observable obs, std::string mc_level, 
                    int pt_min, int pt_max, double pTRL_min, double pTRL_max, int pTRL_bin = 0, bool weighting=false, bool debug=false) {

    // std::string h_name = obs.hist_toextract_name;
    // if (h_name.find("%s") != std::string::npos) h_name.replace(h_name.find("%s"), 2, weighting ? "Weighted" : "" );
    std::string weight_str = weighting ? "Weighted" : "";
    std::string mc_keyword = (mc_level == "truth") ? "gen" : "reco";
    std::string h_name = Form(obs.hist_toextract_name.c_str(), mc_keyword.c_str(), weight_str.c_str());// TODO: FIX WHAT IS BEING ACCESSED HERE, include jet pt and ptrl bin?

    TH3D *h3d = get3DHistAndClone(filename, h_name, pt_min, pt_max, pTRL_bin);
    cout << "making cuts " << pTRL_min << ", " << pTRL_max << " // " << pt_min << ", " << pt_max << endl;
    cout << " x axis name " << h3d->GetXaxis()->GetTitle() << endl;
    cout << " y axis name " << h3d->GetYaxis()->GetTitle() << endl;
    cout << " z axis name " << h3d->GetZaxis()->GetTitle() << endl;
    // h3d->GetYaxis()->SetRangeUser(pTRL_min, pTRL_max);
    // h3d->GetZaxis()->SetRangeUser(pt_min, pt_max);
    // TH1D *hist1D = h3d->ProjectionX(Form("%s_proj_pt%d-%d_ptRLbin%d", h3d->GetName(), pt_min, pt_max, pTRL_bin));

    // making cuts and projecting directly doesn't work... need to specify in ProjectionX() function
    int hist_ptrl_bin = h3d->GetYaxis()->FindBin(pTRL_min);
    int hist_pt_bin = h3d->GetZaxis()->FindBin(pt_min);
    cout << "using name: " << Form("h%s%s_proj_pt%d-%d_ptRLbin%d_%s", obs.name.c_str(), weight_str.c_str(), pt_min, pt_max, pTRL_bin, mc_keyword.c_str()) << " and hist_ptrl_bin: " << hist_ptrl_bin << " and hist_pt_bin: " << hist_pt_bin << endl;
    TH1D *hist1D = h3d->ProjectionX(Form("h%s%s_proj_pt%d-%d_ptRLbin%d_%s", obs.name.c_str(), weight_str.c_str(), pt_min, pt_max, pTRL_bin, mc_keyword.c_str()), hist_ptrl_bin, hist_ptrl_bin, hist_pt_bin, hist_pt_bin); // ybin1, ybin2, zbin1, zbin2
    cout << "integral: " << hist1D->Integral() << endl;

    return hist1D;
}

/* get a typical 2D histogram */
// takes in a THnSparse w/ 4 axes: (x_obs, y_obs, pTRL, jet pT)
TH2D * getObs2DHist(TFile *filename, Observable obs_x, Observable obs_y, std::string jetRname, std::string threshold,
                    int pt_min, int pt_max, int pTRL_bin, double pTRL_min, double pTRL_max, bool debug=false) {
    
    std::string h_name = Form("h2D_%s_vs_%s_Truth%s_%sScaled", obs_y.name.c_str(), obs_x.name.c_str(), jetRname.c_str(), threshold.c_str()); //i.e. "h2D_weights_vs_deltap_Truth%s_%sScaled"

    THnSparse *hnsparse = getTHnSparseAndClone(filename, h_name, pt_min, pt_max, pTRL_bin);
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
    if (write_to_root_file) obj->Write();
    fout->Close();

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
    if (write_to_root_file) hist2D->Write();
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
                          std::string hist_addname, bool weighting, bool logx, bool logy, bool debug=false) {

    if (weighting && obs.name == "weights") return;
    
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

// need an address on Observable so that the original is modified, not a copy
void analyze_1D_obs(TFile * f_in, TH1D * jetpt_inptbin_hist, TLegend& leg, int j, int k,
                    Observable& obs, std::string mc_level, int pt_min, int pt_max, double pTRL_min, double pTRL_max,
                    bool weighting, std::string norm_string, std::string ytitle_norm,
                    std::string ytitle_weight_str, std::string ptRLname_leg, std::string hist_addname,
                    std::string ptname) {

    if (weighting && obs.name == "weights") return;
    cout << " in analyze_1d_obs! with " << obs.name << endl;

    // get histograms
    TH1D *obs_hist = getObs1DHistFrom3D(f_in, obs, mc_level, pt_min, pt_max, pTRL_min, pTRL_max, j, weighting);
    
    // push to vectors
    obs.addHist((TH1D*) obs_hist->Clone(obs_hist->GetName()));

    // format histograms in vector
    Format1DHist(obs, obs.obs_vec[k], jetpt_inptbin_hist, norm_string, colors[j], 0.6, markers[0], obs.axis_label, ytitle_norm + obs.cs_label + ytitle_weight_str, leg, ptRLname_leg, true, hist_addname);
    
    // draw, save, and delete histograms
    draw_save_del_hists(obs, obs.obs_vec[k], ptname, norm_string, hist_addname, false, true);
    delete obs_hist;    

}


void analyze_2D_obs(TFile * f_in, TH1D * jetpt_inptbin_hist, Observable2D obs2D, int pt_min, int pt_max, int pTRL_bin,
                    double pTRL_min, double pTRL_max, bool weighting, std::string jetRname, std::string threshold, std::string norm_string,
                    std::string ytitle_norm, std::string ytitle_weight_str, std::string hist_addname, std::string ptname) {

    Observable obs_x = obs2D.obsx;
    Observable obs_y = obs2D.obsy;
    cout << "in analyze_2d_obs! with " << obs_y.name << "_vs_" << obs_x.name << endl;

    // get histogram
    TH2D * hist2D = getObs2DHist(f_in, obs_x, obs_y, jetRname, threshold, pt_min, pt_max, pTRL_bin, pTRL_min, pTRL_max);
    
    // format histograms
    Format2DHist(obs_x, obs_y, hist2D, jetpt_inptbin_hist, norm_string, ytitle_norm + obs_x.axis_label, ytitle_norm + obs_y.axis_label, hist_addname, true);
    
    // plot and save
    if (obs_x.name == "ptrl") draw_save_del_hists2D(obs2D, hist2D, ptname, norm_string, hist_addname, true, false, true);
    else draw_save_del_hists2D(obs2D, hist2D, ptname, norm_string, hist_addname, false, false, true);
    delete hist2D;

}


// ======================================================= //
//                   2ND MAIN FUNCTION 
// ======================================================= //
void analyze_ptbin(TFile * f_in, vector<Observable> obs_1D_list, vector<Observable2D> obs_2D_list, std::string mc_level, 
                   std::string weightstr, std::string jetRname, std::string threshold,
                   std::string norm_string, int pt_min, int pt_max, const double ptRL_bins[], int n_ptRLbins,
                   vector<vector<double>>& pTRL_vals, vector<vector<double>>& rc_vals, vector<vector<double>>& rc_errors, 
                   bool include_RL0, bool include_RL1, bool debug, bool debug2) {

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
    
    bool weighting = (weightstr == "_Weighted");
    std::string ytitle_weight_str = (weighting) ? " #times #frac{p_{T,1}p_{T,2}}{p_{T,jet}^{2}}" : "";

    std::string ytitle_norm = "";
    if (norm_string == "self_normalized") ytitle_norm = "#frac{1}{N_{pair}} ";
    else if (norm_string == "norm_by_jets") ytitle_norm = "#frac{1}{N_{jet}} ";
    
    vector<double> rc_vec;
    vector<double> rc_err_vec;
    vector<double> pTRLcenters_vec;

    TCanvas *cdumdum = new TCanvas(); // need a canvas so legend can be made
    TLegend *leg = new TLegend(0.6, 0.6, 0.85, 0.87);
    leg->SetTextSize(0.037);
    leg->SetBorderSize(0);
    leg->AddEntry("NULL",Form("%d #leq p_{T, jet} < %d", pt_min, pt_max),"h");
    TLegend *leg_dummy = new TLegend();

    // Names of histograms in the file
    std::string jet_pt_truth_name = Form("h_1Djet_pt_JetPt_Truth%s_%sScaled", jetRname.c_str(), threshold.c_str());
    // std::string deltaptvsEW_truth_name = Form("h_ptvsenergyweights_JetPt_Truth_R%s_%sScaled", jetR.c_str(), threshold.c_str()); //not implemented yet

//    std::string weights_vs_deltap_truth_name = Form("h2D_weights_vs_deltap_Truth%s_%sScaled", jetRname.c_str(), threshold.c_str());
//    std::string weights_vs_p1_truth_name = Form("h2D_weights_vs_p1_Truth%s_%sScaled", jetRname.c_str(), threshold.c_str());
//    std::string weights_vs_deltajt_truth_name = Form("h2D_weights_vs_deltajt_Truth%s_%sScaled", jetRname.c_str(), threshold.c_str());
//    std::string zj_vs_zi_truth_name = Form("h2D_zj_vs_zi_Truth%s_%sScaled", jetRname.c_str(), threshold.c_str());
    
    for ( int j = 0; j < n_ptRLbins; j++ ) {
        int k = j;
        if (!include_RL0) {
            k = j-1;
            if (j == 0) continue; // can add something here to change the filename for ALL
        }
        if (!include_RL1 && j == n_ptRLbins-1) continue;
        
        double pTRL_min = ptRL_bins[j];
        double pTRL_max = ptRL_bins[j+1];
        std::string pTRLname = Form("_pTRL%.1f-%.1f", pTRL_min, pTRL_max);
        std::string pTRLname_leg = Form("#LTpT#GTR_{L} = %.1f-%.1f", pTRL_min, pTRL_max);
        cout << " in pTRL bin" << j << " with " << pTRL_min << " - " << pTRL_max << endl;
        
        std::string hist_addname = weightstr + jetRname + thrname + "_pt" + ptname + pTRLname + "_" + norm_string + "_" + mc_level;
        
        std::string mc_keyword = (mc_level == "truth") ? "Truth" : "Det";
        std::string jet_pt_name = Form("h_1Djet_pt_JetPt_%s%s_%sScaled", mc_keyword.c_str(), jetRname.c_str(), threshold.c_str());
        
        // get histograms
        // get jet pT range - no D0 reconstruction, so don't make D0 cuts
        TH1D * hcorr_jetpt_inptbin_hist;
        Observable obs_jet_pt("jet_pt", true, jet_pt_name, 200, 0, 200, "p_{T,jet}", "#frac{dN}{dp_{T,jet}}");
        if (norm_by_jets_bool) hcorr_jetpt_inptbin_hist = getObs1DHist_direct(f_in, obs_jet_pt.name, true, pt_min, pt_max);
        for (Observable& obs : obs_1D_list) { // passing reference to not make a copy!
            if (obs.obs_bool) analyze_1D_obs(f_in, hcorr_jetpt_inptbin_hist, *leg, j, k, obs, mc_level,
                              pt_min, pt_max, pTRL_min, pTRL_max, weighting, norm_string,
                              ytitle_norm, ytitle_weight_str, pTRLname_leg, hist_addname, ptname);
        }
        
        for (Observable2D obs2D : obs_2D_list) {
            if (obs2D.obs_bool) analyze_2D_obs(f_in, hcorr_jetpt_inptbin_hist, obs2D, pt_min, pt_max, j, pTRL_min, pTRL_max, weighting, jetRname, threshold,
                                norm_string, ytitle_norm, ytitle_weight_str, hist_addname, ptname);
        }
        
        

        
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
//        TH2D * weights_vs_deltap_hist2D = getObs2DHist(f_in, weights_vs_deltap_truth_name, pt_min, pt_max, pTRL_min, pTRL_max);
//        TH2D * weights_vs_p1_hist2D = getObs2DHist(f_in, weights_vs_p1_truth_name, pt_min, pt_max, pTRL_min, pTRL_max);
//        TH2D * weights_vs_deltajt_hist2D = getObs2DHist(f_in, weights_vs_deltajt_truth_name, pt_min, pt_max, pTRL_min, pTRL_max);
//        TH2D * zj_vs_zi_hist2D = getObs2DHist(f_in, zj_vs_zi_truth_name, pt_min, pt_max, pTRL_min, pTRL_max);


//        // Get the X-axis - just printing here
//        TAxis *xAxis = deltap_vec[k]->GetXaxis();
//
//        // Get the number of bins, minimum, and maximum values
//        int nBins = xAxis->GetNbins();
//        double xmin = xAxis->GetXmin();
//        double xmax = xAxis->GetXmax();
//
//        // Calculate the bin width
//        double binWidth = (xmax - xmin) / nBins;
//
//        // Generate the array of bin edges
//        std::vector<double> binEdges;
//        std::cout << "NBINS: " << nBins << ", bin edges: " << endl;
//        for (int i = 0; i <= nBins; ++i) {
//            binEdges.push_back(xmin + i * binWidth);
//            cout << binEdges[i] << " ";
//        }
//        cout << endl;

//        f_out->cd();

    }

    /* do pt bin stuff here */
    std::string hist_all_addname = weightstr + jetRname + thrname + "_pt" + ptname;
    
    // combine RL plots to get 1 plot per pt bin
    for (Observable& obs : obs_1D_list) { // passing reference to not make a copy!
        if (obs.obs_bool) plotandsave_combined_hists(obs, leg, ptname, norm_string, hist_all_addname, weighting, false, true);
    }

}


// don't separate by pt or RL bin
void analyze_jetlevel_observables(TFile* f_in, std::string mc_level, std::string jetRname, std::string thrname,
                                  bool debug, bool debug2) {

    std::string threshold = "1.0";
    std::string mc_keyword = (mc_level == "truth") ? "Truth" : "Det";

    // Names of histograms in the file
    std::string jet_pt_name = Form("h_1Djet_pt_JetPt_%s%s_%sScaled", mc_keyword.c_str(), jetRname.c_str(), threshold.c_str());
    std::string jet_num_const_name = Form("h_Nconst_JetPt_%s%s_%sScaled", mc_keyword.c_str(), jetRname.c_str(), threshold.c_str());

    Observable obs_jet_pt("jet_pt", true, jet_pt_name, 200, 0, 200, "p_{T,jet}", "#frac{dN}{dp_{T,jet}}");
    // Observable obs_jet_const("total_num_const", true, "", 20, 0, 20, "Number Constituents (total)", "");
    Observable obs_jet_const_aftercut("num_const_aftercut", true, jet_num_const_name, 20, 0, 20, "Number Constituents (after threshold cut)", "");

    if (mc_level == "truth") { // could find a better way to do this, but for now... just recreate the file once (otherwise both versions won't save)
        obs_jet_pt.recreate_output_root_file();
        // obs_jet_const.recreate_output_root_file();
        obs_jet_const_aftercut.recreate_output_root_file();
    }


    TH1D * jetpt_hist = getObs1DHist_direct(f_in, obs_jet_pt.hist_toextract_name);
    jetpt_hist->GetXaxis()->SetTitle(obs_jet_pt.axis_label.c_str());
    jetpt_hist->SetNameTitle(Form("jet_pt_hist_%s", mc_level.c_str()), Form("jet_pt_hist_%s", mc_level.c_str()));
    draw_save_del_hists(obs_jet_pt, jetpt_hist, "", "", jetRname + thrname + "_" + mc_level, false, true);

    // TH1D * jet_const = getObs1DHist_direct(f_in, obs_jet_const.hist_toextract_name);
    // jet_const->GetXaxis()->SetTitle(obs_jet_const.axis_label.c_str()); 
    // jet_const->SetNameTitle(Form("jet_const_hist_%s", mc_keyword.c_str()), Form("jet_const_hist_%s", mc_keyword.c_str()));
    // draw_save_del_hists(obs_jet_const, jet_const, "", "", jetRname + thrname + "_" + mc_level, false, true);

    TH1D * jet_const_aftercut = getObs1DHist_direct(f_in, obs_jet_const_aftercut.hist_toextract_name);
    jet_const_aftercut->GetXaxis()->SetTitle(obs_jet_const_aftercut.axis_label.c_str());
    jet_const_aftercut->SetNameTitle(Form("jet_const_aftercut_hist_%s", mc_keyword.c_str()), Form("jet_const_aftercut_hist_%s", mc_keyword.c_str()));
    draw_save_del_hists(obs_jet_const_aftercut, jet_const_aftercut, "", "", jetRname + thrname + "_" + mc_level, false, true);

}


void plot_histograms(TFile* f_in, vector<Observable> obs_1D_list, vector<Observable2D> obs_2D_list,
            std::string weightstr, std::string jetRname, std::string thrname, std::string norm_string,
            const int pt_bins[], int n_bins, const double ptRL_bins[], int n_ptRLbins,
                 bool include_RL0, bool include_RL1, bool debug, bool debug2 ) {

    // // lists
    // std::string jetR_list[] = { "0.4" };
    // std::string threshold_list[] = { "1.0" }; // "0.15", "0.5"

    // // Jet r value
    // for (std::string jetR : jetR_list) {
    //     std::string jetRname = Form("_R%s", jetR.c_str());
        
    //     for (std::string threshold : threshold_list) {
    //         std::string thrname = Form("_t%s", threshold.c_str());
    
    std::string threshold = "1.0"; //TODO: figure out better way to do this?

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

        if (debug) cout << " in pt bin" << i << " with " << pt_min << " - " << pt_max << endl;

        std::string mc_level = "truth";
        analyze_ptbin(f_in, obs_1D_list, obs_2D_list, mc_level, weightstr, jetRname, threshold, norm_string, pt_min, pt_max, ptRL_bins, n_ptRLbins, pTRL_vals, rc_vals, rc_errors, include_RL0, include_RL1, debug, debug2);

        mc_level = "det";
        analyze_ptbin(f_in, obs_1D_list, obs_2D_list, mc_level, weightstr, jetRname, threshold, norm_string, pt_min, pt_max, ptRL_bins, n_ptRLbins, pTRL_vals, rc_vals, rc_errors, include_RL0, include_RL1, debug, debug2);
        
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

//        } // threshold loop
//    } // jetR loop

    
}

// ======================================================= //
//                     MAIN FUNCTION 
// ======================================================= //

void crosscheck_analyze_MC_histograms() {

    gStyle->SetOptStat(0);
    SetStyle();
    
    // CONTOL VARIABLES HERE
    // normed is 0 if unnormalized, 1 for self-normalization 
    // weighted is true if using "Weighted", false if using unweighted
//    bool weighted = false;
//    std::string norm_string = "";
    std::string jetRname = "_R0.4"; // + jetR;
    std::string thrname = "_t1.0"; // + threshold;
    std::string weightstr = ""; //"_xx";
    std::string norm_string = "unnormalized";
    bool include_RL0 = false;
    bool include_RL1 = false;
    
    // setup variables
    bool debug = false;
    bool debug2 = false;

    // int filecounter = 0;
    // int filecounter_cutoff = 500; //total: 5000
    
    // Files
    // const char infile[] = "/rstorage/alice/AnalysisResults/blianggi/dEEC/445125/1132588/scaling/AnalysisResultsFinal.root"; //hiccup
    TString infile;
    infile = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/46040507/1132588/scaling/AnalysisResultsFinal.root"; //perlmutter, after june 2024 // PERLY_FIX: REVERT TO THIS
    // const char infile[] = "~/Documents/research/code/dEEC/new_code/AnalysisResultsFinal.root"; //local
    // const char infile[] = "/Volumes/NO NAME/AnalysisResultsFinal.root"; //local
    if (generator == "herwig") infile = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/storage/herwig/519889-520889/260023/AnalysisResultsFinal_519889_520889.root"; //perlmutter, after dec 2025;
    if (generator == "anchMC") infile = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/47318545/scaling/AnalysisResultsFinal.root"; //perlmutter, after jan 2026;
    if (generator == "pythia_cteq") infile = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/47950666/47835863/scaling/AnalysisResultsFinal.root";

    TFile* root_infile = new TFile(infile, "READ");


    // analysis variables
    const int pt_bins[] = { 20, 40, 60, 80 };
    const int n_bins = sizeof(pt_bins) / sizeof(pt_bins[0]) - 1; //3;
    
    const double ptRL_bins[7] = { 0., 2e-1, 8e-1, 5.0, 10.0, 30.0, 100.0 };
    const int n_ptRLbins = sizeof(ptRL_bins) / sizeof(ptRL_bins[0]) - 1; //gets the columns //6; //7; //5;


    if (debug2) cout << "pt_bins " << n_bins << " n_ptRLbins " << n_ptRLbins << endl;

            
    // ====================================================================================
    
    Observable obs_deltap("deltap", deltap_bool, "%s_deltap%s_unmatchedScaled", 42, 0, 84, "#Deltap", "#frac{dN}{d#Deltap}"); // bin sizes of 2 GeV
    Observable obs_deltajt("deltajt", deltajt_bool, "%s_deltajt%s_unmatchedScaled", 50, 0, 5, "#Deltaj_{T}", "#frac{dN}{d#Deltaj_{T}}");
    Observable obs_weights("weights", ew_bool, "%s_weights_unmatchedScaled", 60, 0, 0.3, "#frac{p_{T,1}p_{T,2}}{p_{T,jet}^{2}}", "#frac{dN}{d[EW]}"); // bin sizes of 0.005
    Observable obs_z("z", twoDhists_bool, "", 50, 0, 1, "z", "#frac{dN}{dz}");

    Observable obs_p1("p1", p1_bool, "not implemented yet", 42, 0, 84, "p", "#frac{dN}{dp}");
    Observable obs_jt1("jt1", jt1_bool, "not implemented yet", 50, 0, 5, "j_{T}", "#frac{dN}{dj_{T}}");

    Observable2D obs_weights_vs_deltap(obs_deltap, obs_weights, twoDhists_bool);
    Observable2D obs_weights_vs_deltajt(obs_deltajt, obs_weights, twoDhists_bool);
    Observable2D obs_weights_vs_p1(obs_p1, obs_weights, twoDhists_bool);
    Observable2D obs_zj_vs_zi(obs_z, obs_z, false);
    
    Observable obs_ptrl("ptrl", deltajt_vs_ptrl_bool, "", 50, 0.2, 35, "#LTp_{T}#GTR_{L}", "");
    Observable2D obs_deltajt_vs_ptrl(obs_ptrl, obs_deltajt, deltajt_vs_ptrl_bool);
    
    vector<Observable> obs_1D_list = { obs_deltap, obs_deltajt, obs_weights };
    vector<Observable2D> obs_2D_list = { obs_weights_vs_deltap, obs_weights_vs_deltajt, obs_zj_vs_zi, obs_deltajt_vs_ptrl };
    for (Observable obs : obs_1D_list) {
        if (obs.obs_bool) obs.recreate_output_root_file();
    }
    for (Observable2D obs2D : obs_2D_list) {
        if (obs2D.obs_bool) obs2D.recreate_output_root_file();
    }
    
    // ====================================================================================

    // analyze jet level observables
    analyze_jetlevel_observables(root_infile, "truth", jetRname, thrname, debug, debug2);
    analyze_jetlevel_observables(root_infile, "det", jetRname, thrname, debug, debug2);

    // analyze for plots
    if (unweighted_bool) {
        if (unnormalized_bool) plot_histograms(root_infile, obs_1D_list, obs_2D_list, weightstr, jetRname, thrname, norm_string, pt_bins, n_bins, ptRL_bins, n_ptRLbins, include_RL0, include_RL1, debug, debug2);
        
        norm_string = "self_normalized";
        if (self_normalized_bool) plot_histograms(root_infile, obs_1D_list, obs_2D_list, weightstr, jetRname, thrname, norm_string, pt_bins, n_bins, ptRL_bins, n_ptRLbins, include_RL0, include_RL1, debug, debug2);

        norm_string = "norm_by_jets";
        if (norm_by_jets_bool) plot_histograms(root_infile, obs_1D_list, obs_2D_list, weightstr, jetRname, thrname, norm_string, pt_bins, n_bins, ptRL_bins, n_ptRLbins, include_RL0, include_RL1, debug, debug2);
    }

    if (weighted_bool) {
        weightstr = "_Weighted";
        norm_string = "unnormalized";
        if (unnormalized_bool) plot_histograms(root_infile, obs_1D_list, obs_2D_list, weightstr, jetRname, thrname, norm_string, pt_bins, n_bins, ptRL_bins, n_ptRLbins, include_RL0, include_RL1, debug, debug2);
        
        norm_string = "self_normalized";
        if (self_normalized_bool) plot_histograms(root_infile, obs_1D_list, obs_2D_list, weightstr, jetRname, thrname, norm_string, pt_bins, n_bins, ptRL_bins, n_ptRLbins, include_RL0, include_RL1, debug, debug2);

        norm_string = "norm_by_jets";
        if (norm_by_jets_bool) plot_histograms(root_infile, obs_1D_list, obs_2D_list, weightstr, jetRname, thrname, norm_string, pt_bins, n_bins, ptRL_bins, n_ptRLbins, include_RL0, include_RL1, debug, debug2);
    }

    root_infile->Close();
    delete root_infile;

    return;
}
