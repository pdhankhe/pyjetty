// ROOT macro to make plots of other correlations
// Beatrice Liang-Gilman (beatrice_lg@berkeley.edu)

void SetStyle(Bool_t graypalette=true) {
  	cout << "Setting style!" << endl;
  
  	gStyle->Reset("Plain");
  	gStyle->SetOptTitle(0);
  	gStyle->SetOptStat(0);
  	// if(graypalette) gStyle->SetPalette(8,0);
  	// else gStyle->SetPalette(1);
    gStyle->SetPalette(kRainbow);
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

void FormatHist(TLegend *l, TH1 *hist, TString text, int markercolor=1, int markerstyle=8, double linealpha=1., bool drawline=false) 
{
    if (drawline) {
        // for (int k=0; k < hist->GetNbinsX();k++){
        //     hist->SetBinError(k+1, 0);
        // }
        hist->SetMarkerStyle(20);
        hist->SetMarkerColorAlpha(markercolor, 0);

        hist->SetFillStyle(0);
        hist->SetLineColorAlpha(markercolor, linealpha);
        hist->SetFillColor(markercolor);
        hist->SetLineStyle(1);
        hist->SetLineWidth(3);
    } else {
        hist->SetLineColor(markercolor);
        hist->SetMarkerColor(markercolor);
        hist->SetMarkerStyle(markerstyle);
        hist->SetMarkerSize(1.5);
    }
    l->AddEntry(hist, text, "pl");

	//gPad->SetTickx(); 
	//gPad->SetTicky(); 
	// h->SetLineWidth(2);
	hist->GetYaxis()->SetTitleOffset(1.05); 
	hist->GetYaxis()->SetTitleSize(0.06); //(0.042);
	hist->GetYaxis()->SetLabelSize(0.05); //(0.042);
	hist->GetYaxis()->SetLabelFont(42);
	hist->GetXaxis()->SetLabelFont(42);
	hist->GetYaxis()->SetTitleFont(42);
	hist->GetXaxis()->SetTitleFont(42);
	hist->GetXaxis()->SetTitleOffset(1.0);
	hist->GetXaxis()->SetTitleSize(0.06); //(0.042);
	hist->GetXaxis()->SetLabelSize(0.05); //(0.042);


    return;
}

void FormatGraph(TGraph *gr, TLegend *l, TString l_text, TString xtitle, TString ytitle, 
                 double markersize=1.5, int markerstyle=20, int markercolor=kBlack, 
                 double markeralpha=1.0) // TString text, int markercolor=1, int markerstyle=8, double linealpha=1., bool drawline=false) 
{
    gr->SetMarkerSize(markersize);
    gr->SetMarkerStyle(markerstyle);
    gr->SetMarkerColorAlpha(markercolor, markeralpha);
    // hist->SetMarkerColorAlpha(markercolor, 0);

    // hist->SetFillStyle(0);
    // hist->SetLineColorAlpha(markercolor, linealpha);
    // hist->SetFillColor(markercolor);
    // hist->SetLineStyle(1);
    // hist->SetLineWidth(3);
    
    l->AddEntry(gr, l_text, "pl");

	//gPad->SetTickx(); 
	//gPad->SetTicky(); 
	// h->SetLineWidth(2);
	gr->GetYaxis()->SetTitleOffset(1.05); 
	gr->GetYaxis()->SetTitleSize(0.06); //(0.042);
	gr->GetYaxis()->SetLabelSize(0.05); //(0.042);
	gr->GetYaxis()->SetLabelFont(42);
	gr->GetXaxis()->SetLabelFont(42);
	gr->GetYaxis()->SetTitleFont(42);
	gr->GetXaxis()->SetTitleFont(42);
	gr->GetXaxis()->SetTitleOffset(1.0);
	gr->GetXaxis()->SetTitleSize(0.06); //(0.042);
	gr->GetXaxis()->SetLabelSize(0.05); //(0.042);

    gr->GetXaxis()->SetTitle(xtitle);
    gr->GetYaxis()->SetTitle(ytitle);


    return;
}

// ptrl is a boolean that says whether ptRL is being plotted (instead of RL)
// --> controls where the cutoff is to not look for the max point
//checkExtra is how many point-to-point slopes after finding a decreasing slope I want to check
int findTopOfCurve(TH1* hist, bool ptrl, int checkExtra=1) {
    
    // define the x axis start of bin search
    double startbinsearchat = ptrl ? 0.1 : 0.01;

    //for each point, find the slope from the 
    int numbins = hist->GetNbinsX();
    int binstart = hist->FindBin(startbinsearchat);
    bool falsealarm = false;

    for (int i=binstart; i<numbins; i++) {

        double y1 = hist->GetBinContent(i);
        double y2 = hist->GetBinContent(i+1);
        double slope_num = y2-y1;

        if (slope_num < 0) {
            for (int j=1; j<=checkExtra; j++) {
                y1 = hist->GetBinContent(i+j);
                y2 = hist->GetBinContent(i+1+j);
                slope_num = y2-y1;
                if (slope_num >= 0) {
                    falsealarm = true;
                    break;
                }
            }
            if (!falsealarm) return i; //return the bin number
        }
        falsealarm = false;
    }

    return 0; //return 0 if nothing found
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


// ======================================================= //
//                     SOME FUNCTIONS 
// ======================================================= //

//get histogram and clone it
THnSparse * getHistAndClone(TFile *f, std::string histname) {
    THnSparse *hsparsejet = (THnSparse*) f->Get(histname.c_str());
    std::string hn = hsparsejet->GetName();
    hn += "_clone";
    THnSparse *hsparsejet_clone = (THnSparse *) hsparsejet->Clone( hn.c_str() ); //TODO: change this name!

    return hsparsejet_clone;
    
}


void applyCuts(THnSparse *hsparse, int pt_min, int pt_max, int RLaxis, double RL_min, double RL_max, 
               bool EWaxis = false) { //todo: don't need paxis anymore??

    hsparse->GetAxis(0)->SetRangeUser(pt_min, pt_max);
    // cout << "RL AXIS is " << RLaxis << " w RL min: " << RL_min << " & RL max: " << RL_max << endl;
    if (RLaxis != -1) {
        cout << " here! RL AXIS is " << RLaxis << " w RL min: " << RL_min << " & RL max: " << RL_max << endl;
        hsparse->GetAxis(RLaxis)->SetRangeUser(RL_min, RL_max);
        // hsparse->GetAxis(RLaxis)->SetRangeUser(1e-5, 1.);
        
        //hsparse->GetAxis(RLaxis)->SetRangeUser(-2, -1); //9.9e-5, 1e-4);//RL_min, RL_max);
        // cout << "UNDERFLOW BIN HAS " << hsparse->Projection(RLaxis)->GetBinContent(0) << " entries" << endl;
        // TH1D* testhist = hsparse->Projection(RLaxis);
        // for (int aa=0; aa<testhist->GetNbinsX(); aa++) {
        //     cout << testhist->GetBinLowEdge(aa) << " ";
        // }
        // cout << endl;
        
    }
    if (EWaxis) { //energy weight axis
        hsparse->GetAxis(4)->SetRangeUser(0., 0.3);
    }

    // cout << "CHECK 1A LOOKING AT NUM EMNRTIRES HERE: " << hsparse->GetEntries() << endl;

    
}

// get the observable histogram
//usually obsaxis is 3, but in new histograms it is 4. For jet level histograms, use 0.
TH1D * getObsHist(TFile *filename, std::string h_name, std::string h_jet_name, int pt_min, int pt_max, 
                  int RLaxis, double RL_min, double RL_max, std::string newhistname, int obsaxis, std::string xtitle,
                  int normalized=0, bool scalebyRLbinwidth=false, double RL_bin_width=1.0, bool EWaxis=false, bool debug=false) { //int n_rebin_bins, double rebinbins[]
    
    cout << "------HNAME is " << h_name << endl;
    cout << "RL AXIS is " << RLaxis << " w RL min: " << RL_min << " & RL max: " << RL_max << endl;
    
    THnSparse *hsparse = getHistAndClone(filename, h_name);
    THnSparse *hsparse_jetlevel = getHistAndClone(filename, h_jet_name);

    // cout << "CHECK 1 LOOKING AT NUM EMNRTIRES HERE: " << hsparse->GetEntries() << endl;
    // cout << "CHECK 1 LOOKING AT NUM ENTRIES: " << hsparse_jetlevel->GetEntries() << endl;

    // if (debug)
    // TH1D *testhist = hsparse->Projection(obsaxis);
    // cout << "underflow entries: " << testhist->GetBinContent(0) << " & overflow entries:" << testhist->GetBinContent(testhist->GetNbinsX()+1) << endl;
    // int sum = 0;
    // for (int aa=0; aa <= testhist->GetNbinsX()+1; aa++) {
    //     sum += testhist->GetBinContent(aa);
    // }
    // cout << " and the total sum is " << sum << endl;
    // cout << " and the integral is " << testhist->Integral() << " & " << testhist->Integral("width") << endl;

    applyCuts(hsparse, pt_min, pt_max, RLaxis, RL_min, RL_max, EWaxis);
    applyCuts(hsparse_jetlevel, pt_min, pt_max, -1, RL_min, RL_max, false); //making no cuts on RL to keep it jet level??

    
    TH1D *h_proj = hsparse->Projection(obsaxis);
    TH1D *h_proj_jetlevel = hsparse_jetlevel->Projection(0); // jet pt axis

    // cout << "CHECK 3 LOOKING AT NUM EMNRTIRES HERE: " << h_proj->GetEntries() << endl;
    // cout << "CHECK 3 LOOKING AT NUM ENTRIES: " << h_proj_jetlevel->GetEntries() << endl;
    

    std::string hname = h_proj->GetName();
    hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_RL" + std::to_string(RL_min) + "-" + std::to_string(RL_max);
    h_proj->SetNameTitle(hname.c_str(), hname.c_str());
    // cout << "HISTOGRAM NAME IS " << h_proj->GetName() << " AND TITLE " << h_proj->GetTitle() << endl;

    // allow rebin or cloning here
    TH1D* hist;
    // int n_obs_bins = sizeof(rebinbins) / sizeof(rebinbins[0]) - 1;
    // if ( n_rebin_bins != 0) {
    //     hist = (TH1D*) h_proj->Rebin(n_rebin_bins, newhistname.c_str(), rebinbins); //h_proj->GetName()
    // } else {
    //     hist = (TH1D*) h_proj->Clone(newhistname.c_str());
    // }
    hist = (TH1D*) h_proj->Clone(newhistname.c_str());
    if (debug) {
        cout << "LOOKING AT NUM BINS: " << hist->GetNbinsX() << endl;
        cout << "LOOKING AT NUM ENTRIES: " << hist->GetEntries() << endl;
    }

    
    // scale by RL bin width if momentum/energy weight bin

    if (scalebyRLbinwidth) hist->Scale(RL_bin_width);

    // normalize
    cout << "IS THIS NORMALIZED? " << normalized << endl;
    if (normalized == 1) { // self normalization
        double numjets = hist->Integral();
        hist->Scale(1/numjets, "width");
        cout << "IN NORMALIZED1" << endl;
    } else if (normalized == 2) { //normalized by the number of jets
        double numjets = h_proj_jetlevel->Integral();
        hist->Scale(1/numjets, "width");
        cout << "IN NORMALIZED2" << endl;
    }
    // cout << "There are " << h_proj->GetEntries() << " pair entries in this pt bin" << endl;
    // cout << "There are " << h_proj_jetlevel->GetEntries() << " jet entries in this pt bin" << endl;

    // hist->GetXaxis()->SetTitle(xtitle.c_str()); //("#it{R}_{L}");
    // // THESE AXES LABELS ARE DEF WRONG!!
    // if (normalized != 0) hist->GetYaxis()->SetTitle( Form("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d%s}", xtitle.c_str()) );
    // else hist->GetYaxis()->SetTitle( Form("#frac{d#it{N}_{EEC}}{d%s}", xtitle.c_str()) ); 


    delete hsparse;
    delete hsparse_jetlevel;
    delete h_proj_jetlevel;
    delete h_proj;

    return hist;

}

// get 2D observable histogram
TH2D * get2DHist(TFile *filename, std::string h_name, std::string h_jet_name, int pt_min, int pt_max,
                 int RLaxis, double RL_min, double RL_max, std::string newhistname, int xaxis_axis, int yaxis_int,
                 std::string xtitle, std::string ytitle, int normalized=0, bool scalebyRLbinwidth=false, 
                 double RL_bin_width=1.0, bool debug=false) 
{
    cout << "------HNAME is " << h_name << endl;
    cout << "RL AXIS is " << RLaxis << " w RL min: " << RL_min << " & RL max: " << RL_max << endl;
    
    THnSparse *hsparse = getHistAndClone(filename, h_name);
    THnSparse *hsparse_jetlevel = getHistAndClone(filename, h_jet_name);

    applyCuts(hsparse, pt_min, pt_max, RLaxis, RL_min, RL_max, false);
    applyCuts(hsparse_jetlevel, pt_min, pt_max, -1, RL_min, RL_max, false); //making no cuts on RL to keep it jet level??
    hsparse->GetAxis(5)->SetRangeUser(0, 0.3); //for the energy weights axis

    TH2D *h_proj = hsparse->Projection(yaxis_int, xaxis_axis); //ydim, xdim
    TH1D *h_proj_jetlevel = hsparse_jetlevel->Projection(0); // jet pt axis
    
    std::string hname = h_proj->GetName();
    hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_RL" + std::to_string(RL_min) + "-" + std::to_string(RL_max);
    h_proj->SetNameTitle(hname.c_str(), hname.c_str());

    // allow rebin or cloning here
    TH2D* hist2D;
    hist2D = (TH2D*) h_proj->Clone(newhistname.c_str());

    
    // scale by RL bin width if momentum/energy weight bin
    if (scalebyRLbinwidth) {
        hist2D->Scale(RL_bin_width);
        // hist2D->GetYaxis()->Scale(RL_bin_width);
    }

    // normalize
    // cout << "IS THIS NORMALIZED? " << normalized << endl;
    if (normalized == 1) { // self normalization
        double numjets = hist2D->Integral(); // is this right for 2D?
        hist2D->Scale(1/numjets, "width");
    } else if (normalized == 2) { //normalized by the number of jets
        double numjets = h_proj_jetlevel->Integral();
        hist2D->Scale(1/numjets, "width");
    }

    // label axes
    hist2D->GetXaxis()->SetTitle(xtitle.c_str());
    hist2D->GetYaxis()->SetTitle(ytitle.c_str());

    delete hsparse;
    delete hsparse_jetlevel;
    delete h_proj;
    delete h_proj_jetlevel;

    return hist2D;
}

void fix_up_2D(TH2D * hist2D, double bounds[], double ptmax) {

    hist2D->GetXaxis()->SetRangeUser(0, ptmax);
    hist2D->GetZaxis()->SetRangeUser(bounds[0], bounds[1]);

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
    // pol3_func->SetParNames(“r”,“A”,“B”,“C”);
    // pol3_func->SetParLimits(0,0.002,1.5);
    // pol3_func->SetParLimits(1,0.0000013,-0.10);
    // pol3_func->SetParLimits(2,2.5,-0.20);
    // pol3_func->SetParLimits(3,0.001,0.50);


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


// plot all RL bins in one plot
void plot_combined_hists(TCanvas *can_all, vector<TH1D*> h_arr, TLegend *l, 
                          int pt_max, bool mom_axis, double RL_bin_width[],
                          double pl_axis_cut=-1, bool debug=false) {

    // go into canvas
    can_all->cd();

    // if momentum axis, adjust x bounds accordingly
    size_t length = h_arr.size();
    for (int j=0; j<length; j++) {
        if (mom_axis) {
            h_arr[j]->Rebin(4);
            h_arr[j]->GetXaxis()->SetRangeUser(0, pt_max+5);
            // h_arr[j]->Scale(RL_bin_width[j]); // this needs to be done before normalization
        }
        if (pl_axis_cut > 0) {
            h_arr[j]->GetXaxis()->SetRangeUser(0, pl_axis_cut);
        }
    }

    // set maximum based on maximum of all curves
    double max = 0;
    for (int j=0; j<length; j++) {
        double max_cand = h_arr[j]->GetMaximum();
        if (max_cand > max) max = max_cand;
    }
    if (debug) cout << "max is " << max << " which goes to " << max*1.5 << endl;
    h_arr[0]->SetMaximum( max * 1.5 );

    // draw!
    for (int j=0; j<length; j++) {
        h_arr[j]->Draw("same");
    }

    // can_all->Modified();
    // can_all->Update();
    l->Draw("same");
}

// save and delete canvas
void save_del_canvas(TCanvas *c, std::string outdir, std::string corr_name, std::string normstring, std::string weightstr, 
                     std::string fitrange, std::string ptname, std::string RLname, std::string jetR, std::string add_name, 
                     bool individual_obs) {
    
    std::string fname = "";
    if (individual_obs) fname = outdir + "individuals/" + corr_name + "/" + normstring + "/corrhist_" + corr_name + weightstr + "_all" + fitrange + ptname + RLname + "_R" + jetR + add_name + ".pdf";
    else fname = outdir + normstring + "/corrhist_" + corr_name + weightstr + "_all" + fitrange + ptname + "_R" + jetR + add_name + ".pdf";       
    const char* fname_c = fname.c_str();
    
    c->SaveAs(fname_c);
    delete c;
}

// ======================================================= //
//                   2ND MAIN FUNCTION 
// ======================================================= //

void plot_histograms(TFile* f, std::string add_name, int normed, bool weighted, bool include_RL0) {

    Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
    Double_t marker_size = 1.5;
    Double_t colors[16] = {kGray, kMagenta, kGreen+2, kBlue, kOrange+1, kViolet+1, kRed, kYellow+1, kCyan+1};
    // Double_t colors[16] = {kRed, kGreen+2, kBlue, kRed+1, kGreen+1, kBlue+1, kRed+2, kGreen+2, kBlue+2, kRed+3, kGreen+3, kBlue+3, kOrange+1, kViolet+1, kYellow+1, kCyan+1};

    //    gROOT->SetBatch(); //prevents plots from showing up
    const int pt_bins[] = { 20, 40, 60, 80 };
    const std::string pt_bin_names[] = {"20-40", "40-60", "60-80"};
    const double plaxis_cut[] = { 20, 30, 40 };
    const std::string normstring_dirs[] = {"unnormalized", "self_normalized", "norm_by_jets"};
    std::string normstring_dir = "";
    if (normed == 0) normstring_dir = normstring_dirs[0];
    else if (normed == 1) normstring_dir = normstring_dirs[1];
    else normstring_dir = normstring_dirs[2];
    // const double RL_bins[] = { 1e-2, 2e-2, 3e-2, 7.5e-2, 2e-1, 4e-1};
    // const double RL_bins[3][8] = { { 1e-5, 1e-2, 3e-2, 7e-2, 1.5e-1, 3e-1, 4e-1, 1 }, 
    //                              { 1e-5, 1e-2, 2.5e-2, 4e-2, 8e-2, 2.5e-1, 4e-1, 1 }, 
    //                              { 1e-5, 1e-2, 2.5e-2, 3e-2, 4.5e-2, 2e-1, 4e-1, 1 } };
    const int n_bins = 3;
    double RL_bins[3][8] = {{0}, {0}, {0}};
    double RL_bin_width[3][8] = {{0}, {0}, {0}}; //{0};
    double RL_bin_centers[3][8] = {{0}, {0}, {0}};
    int n_RLbins = 0;
    double temp[3][8] = { { 0, 1e-2, 3e-2, 7e-2, 1.5e-1, 3e-1, 4e-1, 1 }, 
                          { 0, 1e-2, 2.5e-2, 4e-2, 8e-2, 2.5e-1, 4e-1, 1 }, 
                          { 0, 1e-2, 2.5e-2, 3e-2, 4.5e-2, 2e-1, 4e-1, 1 } };   

    double ptvsew_normbounds[3][2] = { { 8e-8, 3e-1 }, { 1e-2, 5e2 }, { 1e-4, 2e2 }}; //normed: 0 = unnormalized, 1= self norm, 2 = norm by jets

    if (include_RL0) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 7; ++j) {
                RL_bins[i][j] = temp[i][j];
                RL_bin_width[i][j] = temp[i][j+1] - temp[i][j];
                RL_bin_centers[i][j] = (temp[i][j+1] + temp[i][j])/2;
            }
        }
        n_RLbins = 6; //7; //5;
    } else {
        for (int i = 0; i < 3; ++i) {
            // cout << "--- HERE LIES PT " << i << endl; 
            for (int j = 0; j < 6; ++j) {
                RL_bins[i][j] = temp[i][j+1];
                RL_bin_width[i][j] = temp[i][j+2] - temp[i][j+1];
                RL_bin_centers[i][j] = (temp[i][j+2] + temp[i][j+1])/2;
                // cout << "RL BIN WIDTH HERE" << RL_bins[i][j] << endl;
                // cout << " AND DIFF " << RL_bin_width[j] << endl;
            }
        }
        n_RLbins = 5;
        
        int temp_colors[16];
        std::copy(colors, colors + 16, temp_colors);
        for (int icol = 0; icol < 8; ++icol) {
            if ( icol < 7 ) colors[icol] = temp_colors[icol+1];
            else if ( icol == 7 ) colors[icol] = temp_colors[0]; 
        }
    }
    

    TString label1 = "";
    TString label2 = "";  
    // vector<TString> label;

    std::string weightstr = "";
    if (weighted == true) {
        weightstr = "_Weighted";
    }

    std::string norm_string = "";
    if (normed == 0) {
        norm_string = "_unnormalized";
    } else if (normed == 1) {
        norm_string = "_selfnorm";
    } else if (normed == 2) {
        norm_string = "_normbyjets";
    }

    std::string RL0_string = "";
    if (include_RL0) RL0_string = "_withRL0";

    add_name += norm_string + RL0_string;

    // Jet r value
    std::string jetR_list[] = { "0.4" };
    for (std::string jetR : jetR_list) {

        std::string threshold_list[] = { "1.0" }; // "0.15", "0.5"
        for (std::string threshold : threshold_list) {

            // Names of histograms in the file
            const std::string deltap_truth_name = Form("h_corr_deltap%s_JetPt_Truth_R0.4_%sScaled", weightstr.c_str(), threshold.c_str());
            const std::string deltapt_truth_name = Form("h_corr_deltapt%s_JetPt_Truth_R0.4_%sScaled", weightstr.c_str(), threshold.c_str());
            const std::string deltapl_truth_name = Form("h_corr_deltapl%s_JetPt_Truth_R0.4_%sScaled", weightstr.c_str(), threshold.c_str());
            const std::string oppcharge_truth_name = Form("h_corr_oppcharge%s_JetPt_Truth_R0.4_%sScaled", weightstr.c_str(), threshold.c_str());
            const std::string samecharge_truth_name = Form("h_corr_samecharge%s_JetPt_Truth_R0.4_%sScaled", weightstr.c_str(), threshold.c_str());
            const std::string unweightedRL_truth_name = Form("h_corr_unweightedRL%s_JetPt_Truth_R0.4_%sScaled", weightstr.c_str(), threshold.c_str());
            const std::string energyweights_truth_name = Form("h_corr_energyweights%s_JetPt_Truth_R0.4_%sScaled", weightstr.c_str(), threshold.c_str());
            const std::string baryonmeson_truth_name = Form("h_corr_baryonmeson%s_JetPt_Truth_R0.4_%sScaled", weightstr.c_str(), threshold.c_str());
            // std::cout << "THRESHOLD NAME TEST" << deltap_truth_name << std::endl;
            
            const std::string jet_pt_truth_name = Form("h_jet_pt_JetPt_Truth_R0.4_%sScaled", threshold.c_str());
            const std::string deltaptvsEW_truth_name = Form("h_ptvsenergyweights_JetPt_Truth_R0.4_%sScaled", threshold.c_str());

            for (int i = 0; i < n_bins; i++) {
                cout << "in pt bin" << i << endl;
                int pt_min = pt_bins[i];
                int pt_max = pt_bins[i+1];
                std::string ptname = "_pt" + std::to_string(pt_min) + '-' + std::to_string(pt_max);

                // Output directory
                std::string outdir = "plots/thirdattempt/"+pt_bin_names[i]+"/";//"plots/test/";
                // Output file for binned results
                std::string outfile = outdir + "AnalysisResultsFinal" + add_name + ".root";
                TFile* f_out = new TFile(outfile.c_str(), "RECREATE");

                // define pt related variables
                TString ptbin = TString::Format("%d #leq #it{p}_{T}^{ch. jet} < %d GeV/#it{c}, #font[122]{|}#it{#eta}_{jet}#font[122]{|} #leq 0.5", pt_min, pt_max);
                TString ptD = TString::Format("5 #leq #it{p}_{T}^{D^{0}} < %d GeV/#it{c}, #font[122]{|}#it{y}_{D^{0}}#font[122]{|} #leq 0.8", pt_max);
                

                TCanvas* c_deltap_all = new TCanvas();
                ProcessCanvas(c_deltap_all);
                gPad->SetLogy();

                TCanvas* c_deltap_all_fit05 = new TCanvas();
                ProcessCanvas(c_deltap_all_fit05);
                gPad->SetLogy();

                TCanvas* c_deltapt_all = new TCanvas();
                ProcessCanvas(c_deltapt_all);
                gPad->SetLogy();

                TCanvas* c_deltapl_all = new TCanvas();
                ProcessCanvas(c_deltapl_all);
                gPad->SetLogy();

                TCanvas* c_charge_all = new TCanvas();
                ProcessCanvas(c_charge_all);
                gPad->SetLogx();

                TCanvas* c_chargeratio_all = new TCanvas();
                ProcessCanvas(c_chargeratio_all);
                // gPad->SetLogx();

                TCanvas* c_energyweights_all = new TCanvas();
                ProcessCanvas(c_energyweights_all);
                // gPad->SetLogx();
                gPad->SetLogy();

                TCanvas* c_baryonmeson_all = new TCanvas();
                ProcessCanvas(c_baryonmeson_all);

                TCanvas* c_baryonmesonratios_all = new TCanvas();
                ProcessCanvas(c_baryonmesonratios_all);


                TLegend* l; // = new TLegend(0.17, 0.65, 0.5, 0.85);
                l = new TLegend(0.58, 0.6, 0.78, 0.87); //0.1797168,0.650741,0.4562155,0.8885185,""); //(0.17, 0.4, 0.5, 0.53);
                l->SetTextSize(0.037);
                l->SetBorderSize(0);

                TLegend* l_bottomleft; // = new TLegend(0.17, 0.65, 0.5, 0.85);
                l_bottomleft = new TLegend(0.2, 0.2, 0.4, 0.47); //0.1797168,0.650741,0.4562155,0.8885185,""); //(0.17, 0.4, 0.5, 0.53);
                l_bottomleft->SetTextSize(0.037);
                l_bottomleft->SetBorderSize(0);

                TLegend* l2; // = new TLegend(0.17, 0.65, 0.5, 0.85);
                l2 = new TLegend(0.1797168,0.700741,0.4562155,0.8885185,""); //(0.17, 0.4, 0.5, 0.53);
                // l2->SetTextSize(0.037);
                // l2->SetBorderSize(0);

                TLegend* l_bm = new TLegend(0.63, 0.6, 0.83, 0.87);
                TLegend* l_bmratios = new TLegend(0.58, 0.5, 0.88, 0.87);

                TLegend* l_right = new TLegend(0.75, 0.5, 0.85, 0.87);
                l_right->SetTextSize(0.037);
                l_right->SetBorderSize(0);

                // make histogram array bc otherwise they get overwritten and won't save S M H
                vector<TH1D*> hcorr_deltap_truth_arr;
                vector<TH1D*> hcorr_deltapt_truth_arr;
                vector<TH1D*> hcorr_deltapl_truth_arr;
                vector<TH1D*> hcorr_oppcharge_truth_arr;
                vector<TH1D*> hcorr_samecharge_truth_arr;
                vector<TH1D*> hcorr_chargeratio_ENC_truth_arr;
                vector<TH1D*> hcorr_oppcharge_ENC_truth_arr;
                vector<TH1D*> hcorr_samecharge_ENC_truth_arr;
                // vector<TH1D*> hcorr_unweightedRL_truth_arr;
                vector<TH1D*> hcorr_energyweights_truth_arr;
                vector<TH1D*> hcorr_baryonmeson_truth_arr;

                vector<TString> label;
                TH1D *hdummyRL;
                TH1D *hdummyRL2;
                double y_chargeratio_arr[n_RLbins];
                double y_baryonbaryon_arr[n_RLbins];
                double y_baryonmeson_arr[n_RLbins];
                double y_mesonmeson_arr[n_RLbins];
                    
                for (int j=0; j < n_RLbins; j++) {

                    double RL_min = RL_bins[i][j];
                    double RL_max = RL_bins[i][j+1];
                    cout << "in RL bin" << j << " with " << RL_min << " - " << RL_max << endl;

                    
                    double maxy = 0;                    

                    std::string jetRname = "_R" + jetR;
                    std::string thrname = "_t" + threshold;
                    std::string RLname = Form("_RL%.3f-%.3f", RL_min, RL_max); 
                    std::string hist_addname = jetRname + thrname + ptname + RLname;
                    // cout << "THE HIST ADDNAME IS " << hist_addname << endl;
                    // cout << "deltsap_truth_name " << deltap_truth_name << endl;
                    // double obs_rebins[] = {0.  , 0.05, 0.1 , 0.15, 0.2 , 0.25, 0.3 , 0.35, 0.4 , 0.45, 0.5 ,
                    //     0.55, 0.6 , 0.65, 0.7, 0.75, 0.8 , 0.85, 0.9 , 0.95, 1.  };
                    // double obs_rebins[] = { 0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10.,
                                        // 11.,  12.,  13.,  14.,  15.,  16.,  17.,  18.,  19.,  20.,  21.,
                                        // 22.,  23.,  24.,  25.,  26.,  27.,  28.,  29.,  30.,  31.,  32.,
                                        // 33.,  34.,  35.,  36.,  37.,  38.,  39.,  40.,  41.,  42.,  43.,
                                        // 44.,  45.,  46.,  47.,  48.,  49.,  50.,  51.,  52.,  53.,  54.,
                                        // 55.,  56.,  57.,  58.,  59.,  60.,  61.,  62.,  63.,  64.,  65.,
                                        // 66.,  67.,  68.,  69.,  70.,  71.,  72.,  73.,  74.,  75.,  76.,
                                        // 77.,  78.,  79.,  80.,  81.,  82.,  83.,  84.,  85.,  86.,  87.,
                                        // 88.,  89.,  90.,  91.,  92.,  93.,  94.,  95.,  96.,  97.,  98.,
                                        // 99., 100. };
                    // int n_obs_rebins = sizeof(obs_rebins) / sizeof(obs_rebins[0]) - 1;


                    // make a canvas for each pt range
                    TCanvas* c_deltap = new TCanvas();
                    ProcessCanvas(c_deltap);
                    gPad->SetLogy();

                    TCanvas* c_deltapt = new TCanvas();
                    ProcessCanvas(c_deltapt);
                    gPad->SetLogy();

                    TCanvas* c_deltapl = new TCanvas();
                    ProcessCanvas(c_deltapl);
                    gPad->SetLogy();

                    TCanvas* c_deltapt_vs_ew = new TCanvas();
                    ProcessCanvas(c_deltapt_vs_ew, true);
                    // gPad->SetLogy();


                    TCanvas* c_charge = new TCanvas();
                    ProcessCanvas(c_charge);

                    TCanvas* c_chargeratio = new TCanvas();
                    ProcessCanvas(c_chargeratio);
                    gPad->SetLogx(); // gPad->SetLogy();

                    TCanvas* c_unweightedRL = new TCanvas();
                    ProcessCanvas(c_unweightedRL);

                    TCanvas* c_energyweights = new TCanvas();
                    ProcessCanvas(c_energyweights);

                    
                    // Open histograms
                    TH1D *hcorr_deltap_truth = getObsHist(f, deltap_truth_name, jet_pt_truth_name, 
                                             pt_min, pt_max, i+1, RL_min, RL_max, "h_corr_deltap_Truth" + hist_addname, 
                                             4, "#Delta p", normed, true, RL_bin_width[i][j], false);
                    TH1D *hcorr_deltapt_truth = getObsHist(f, deltapt_truth_name, jet_pt_truth_name, 
                                             pt_min, pt_max, i+1, RL_min, RL_max, "h_corr_deltapt_Truth" + hist_addname, 
                                             4, "#Delta p_{T}", normed, true, RL_bin_width[i][j], false);
                    TH1D *hcorr_deltapl_truth = getObsHist(f, deltapl_truth_name, jet_pt_truth_name, 
                                             pt_min, pt_max, i+1, RL_min, RL_max, "h_corr_deltapl_Truth" + hist_addname, 
                                             4, "#Delta p_{L}", normed, true, RL_bin_width[i][j], false);
                    TH1D *hcorr_oppcharge_truth = getObsHist(f, oppcharge_truth_name, jet_pt_truth_name, 
                                             pt_min, pt_max, i+1, RL_min, RL_max, "h_corr_oppcharge_Truth" + hist_addname, 
                                             4, "q_{1}q_{2}", normed, true, RL_bin_width[i][j], false);
                    TH1D *hcorr_samecharge_truth = getObsHist(f, samecharge_truth_name, jet_pt_truth_name, 
                                             pt_min, pt_max, i+1, RL_min, RL_max, "h_corr_samecharge_Truth" + hist_addname, 
                                             4, "q_{1}q_{2}", normed, true, RL_bin_width[i][j], false);
                    TH1D *hcorr_energyweights_truth = getObsHist(f, energyweights_truth_name, jet_pt_truth_name, 
                                             pt_min, pt_max, i+1, RL_min, RL_max, "h_corr_energyweights_Truth" + hist_addname, 
                                             4, "p_{T,1}p_{T,2} / p_{T, jet}^{2}", normed, true, RL_bin_width[i][j], true);
                    TH1D *hcorr_baryonmeson_truth = getObsHist(f, baryonmeson_truth_name, jet_pt_truth_name, 
                                             pt_min, pt_max, i+1, RL_min, RL_max, "h_corr_baryonmeson_Truth" + hist_addname, 
                                             4, "baryon/meson pairs", normed, true, RL_bin_width[i][j], false);

                    TH1D *hcorr_oppcharg_ENC_truth_unnorm = getObsHist(f, oppcharge_truth_name, jet_pt_truth_name, 
                                             pt_min, pt_max, i+1, RL_min, RL_max, "h_corr_oppcharge_ENC_Truth" + hist_addname, 
                                             i+1, "R_{L}", normed, true, RL_bin_width[i][j], false); //make normalization = 2?
                    TH1D *hcorr_samecharge_ENC_truth_unnorm = getObsHist(f, samecharge_truth_name, jet_pt_truth_name, 
                                             pt_min, pt_max, i+1, RL_min, RL_max, "h_corr_samecharge_ENC_Truth" + hist_addname, 
                                             i+1, "R_{L}", normed, true, RL_bin_width[i][j], false); //make normalization = 2?
                    TH2D *hcorr_deltaptvsEW_truth = get2DHist(f, deltaptvsEW_truth_name, jet_pt_truth_name,
                                               pt_min, pt_max, i+1, RL_min, RL_max, "h_corr_deltaptvsEW" + hist_addname,
                                               4, 5, "#Delta p_{T}", "p_{T,1}p_{T,2} / p_{T, jet}^{2}", normed, true, 
                                               RL_bin_width[i][j]);

                    //------------//------------//------------//------------//------------//
                    // make fits
                    // TF1 *fitFunc = new TF1("fitFunc", myFunction, 0, 10, 2); // 2 parameters
                    // fitFunc->SetParameters(1, 0); // Initial guess for parameters
                    // hist->Fit("fitFunc");

                    // double p0 = fitFunc->GetParameter(0); // Get first parameter
                    // double p0_err = fitFunc->GetParError(0); // Get error of the first parameter

                    // graph->Draw("AP"); // Draw the graph with axis and points
                    // fitFunc->Draw("same"); // Draw the fit function on the same canvas

                    // TH2D *hcorr_deltapt_vs_ew = TH2D()
                    
                    //------------//------------//------------//------------//------------//

                    // get jet pT range - no D0 reconstruction, so don't make D0 cuts
                    // thnsparse axes: 0=jet pt, 1=RL (20 < pt < 40), 2=RL (40 < pt < 60), 3=RL (60 < pt < 80), 4 = observable 

                    

                    // dummy histograms!
                    THnSparse* hsparsejet_oppcharge_truth_dummy = (THnSparse*) f->Get(oppcharge_truth_name.c_str());
                    if (j==0) hdummyRL = hsparsejet_oppcharge_truth_dummy->Projection(1);
                    if (j==0) hdummyRL2 = hsparsejet_oppcharge_truth_dummy->Projection(1);




                    cout << "HOLA??  " << endl;

                    //print together
                    hcorr_deltap_truth_arr.push_back((TH1D*) hcorr_deltap_truth->Clone(hcorr_deltap_truth->GetName()));
                    hcorr_deltapt_truth_arr.push_back((TH1D*) hcorr_deltapt_truth->Clone(hcorr_deltapt_truth->GetName()));
                    hcorr_deltapl_truth_arr.push_back((TH1D*) hcorr_deltapl_truth->Clone(hcorr_deltapl_truth->GetName()));
                    hcorr_oppcharge_truth_arr.push_back((TH1D*) hcorr_oppcharge_truth->Clone(hcorr_oppcharge_truth->GetName()));
                    hcorr_samecharge_truth_arr.push_back((TH1D*) hcorr_samecharge_truth->Clone(hcorr_samecharge_truth->GetName()));
                    // hcorr_unweightedRL_truth_arr.push_back((TH1D*) hcorr_unweightedRL_truth_unnorm->Clone(hcorr_unweightedRL_truth_unnorm->GetName()));
                    hcorr_energyweights_truth_arr.push_back((TH1D*) hcorr_energyweights_truth->Clone(hcorr_energyweights_truth->GetName()));
                    hcorr_baryonmeson_truth_arr.push_back((TH1D*) hcorr_baryonmeson_truth->Clone(hcorr_baryonmeson_truth->GetName()));

                    label.push_back(Form("R_{L} = %.3f-%.3f", RL_min, RL_max));

                    for (int k=0; k < hcorr_deltap_truth_arr[j]->GetNbinsX();k++){
                        hcorr_deltap_truth_arr[j]->SetBinError(k+1, 0);
                        hcorr_deltapt_truth_arr[j]->SetBinError(k+1, 0);
                        hcorr_deltapl_truth_arr[j]->SetBinError(k+1, 0);
                    }
                    for (int k=0; k < hcorr_energyweights_truth_arr[j]->GetNbinsX();k++){
                        hcorr_energyweights_truth_arr[j]->SetBinError(k+1, 0);
                    }

                    // format histograms here
                    // regarding the legends - l is the default (top right)
                    // l_bottomleft is in the bottom left 
                    // only need to add one set of observables to legend - the rest get sent to l2 as the "rubbish legend" that doesn't get plotted
                    FormatHist(l, hcorr_deltap_truth_arr[j], label[j], colors[j], markers[0], true);
                    FormatHist(l_bottomleft, hcorr_deltapt_truth_arr[j], label[j], colors[j], markers[0], true);
                    FormatHist(l2, hcorr_deltapl_truth_arr[j], label[j], colors[j], markers[0], true);
                    FormatHist(l2, hcorr_oppcharge_truth_arr[j], label[j], colors[j], markers[0], false);
                    FormatHist(l2, hcorr_samecharge_truth_arr[j], label[j], colors[j], markers[2], false);
                    // FormatHist(l, hcorr_unweightedRL_truth, label[j], colors[j], markers[0]);
                    FormatHist(l2, hcorr_energyweights_truth_arr[j], label[j], colors[j], markers[0], true);
                    
                    // Format histograms for plotting (this order needed to keep legend in order and graphs lookin good)
                    cout << "hello??  " << endl;

                    // exponential fit function
                    // TF1 *expfit = new TF1("expfit", expFunc, 0, pt_max+5, 2); // 2 parameters

                    c_deltap->cd();
                    hcorr_deltap_truth_arr[j]->Draw();
                    cout << "CHECKPOINT 1" << endl;
                    // hcorr_deltap_truth_arr[j]->Fit("expfit", "", "", 0, 10);
                    cout << "CHECKPOINT 2" << endl;

                    c_deltap_all_fit05->cd();
                    TH1D *h_fitclone = makeFitForClone(hcorr_deltap_truth, 0, 5);
                    // if (j==0) 
                    hcorr_deltap_truth_arr[j]->Draw();
                    fit_hist(c_deltap_all_fit05, h_fitclone, 0, 5);

                    // double A_fit = expfit->GetParameter(0); // Get first parameter
                    // double B_fit = expfit->GetParameter(1); 
                    // // double p0_err = fitFunc->GetParError(0); // Get error of the first parameter

                    // // graph->Draw("AP"); // Draw the graph with axis and points
                    // c_deltap->cd();
                    // expfit->Draw("same"); // Draw the fit function on the same canvas
                    // cout << "A_fit " << A_fit << " and B_fit " << B_fit << endl;


                    c_deltapt->cd();
                    hcorr_deltapt_truth_arr[j]->Draw();
                    fit_hist(c_deltapt, hcorr_deltapt_truth_arr[j], 0, 10);

                    c_deltapl->cd();
                    hcorr_deltapl_truth_arr[j]->Draw();
                    // fit_hist(c_deltapt, hcorr_deltapt_truth_arr[j], 0, 10);

                    
                    c_deltapt_vs_ew->cd();
                    // gStyle->SetPadRightMargin(0.15);
                    gPad->SetLogz();
                    fix_up_2D(hcorr_deltaptvsEW_truth, ptvsew_normbounds[normed], pt_max);
                    hcorr_deltaptvsEW_truth->Draw("colz");



                    c_charge->cd();
                    hcorr_oppcharge_truth_arr[j]->Draw();
                    hcorr_samecharge_truth_arr[j]->Draw("same");

                    c_chargeratio->cd();
                    hcorr_oppcharge_ENC_truth_arr.push_back((TH1D*) hcorr_oppcharg_ENC_truth_unnorm->Clone(hcorr_oppcharg_ENC_truth_unnorm->GetName()));
                    hcorr_samecharge_ENC_truth_arr.push_back((TH1D*) hcorr_samecharge_ENC_truth_unnorm->Clone(hcorr_samecharge_ENC_truth_unnorm->GetName()));
                    FormatHist(l2, hcorr_oppcharge_ENC_truth_arr[j], label[j], colors[j], markers[0], false);
                    FormatHist(l2, hcorr_samecharge_ENC_truth_arr[j], label[j], colors[j], markers[2], false);

                    std::string charge_ratio_name = "hcorr_chargeratio_ENC" + ptname + RLname;
                    TH1D* hcorr_chargeratio_ENC_truth = (TH1D*) hcorr_oppcharge_ENC_truth_arr[j]->Clone(charge_ratio_name.c_str());
                    hcorr_chargeratio_ENC_truth->Divide(hcorr_samecharge_ENC_truth_arr[j]);
                    FormatHist(l2, hcorr_chargeratio_ENC_truth, label[j], colors[j], markers[2]);
                    hcorr_chargeratio_ENC_truth_arr.push_back((TH1D*) hcorr_chargeratio_ENC_truth->Clone(hcorr_chargeratio_ENC_truth->GetName()));
                    hcorr_chargeratio_ENC_truth_arr[j]->Draw();
                    // rp->GetLowYaxis()->SetNdivisions(505);
                    // c->Update();

                    c_energyweights->cd();
                    hcorr_energyweights_truth_arr[j]->Draw();
                    

                    // do combined graphs now!
                    
                    // this one has RL bin on the x axis, and y axis is #opp-sign-pairs/#same-sign-pairs
                    // later will save to c_chargeratio_all
                    double num_opp = hcorr_oppcharge_truth_arr[j]->GetBinContent(hcorr_oppcharge_truth_arr[j]->FindBin(-1));
                    double num_same = hcorr_samecharge_truth_arr[j]->GetBinContent(hcorr_samecharge_truth_arr[j]->FindBin(1));
                    y_chargeratio_arr[j] = num_opp/num_same;

                    //save baryon/meson pairs info
                    // 1 = p+p, 0 = p+pi, -1 = pi+pi, 2 = any other pairs
                    y_baryonbaryon_arr[j] = hcorr_baryonmeson_truth->GetBinContent(hcorr_baryonmeson_truth->FindBin(1));
                    y_baryonmeson_arr[j] = hcorr_baryonmeson_truth->GetBinContent(hcorr_baryonmeson_truth->FindBin(0));
                    y_mesonmeson_arr[j] = hcorr_baryonmeson_truth->GetBinContent(hcorr_baryonmeson_truth->FindBin(-1));
                    


                    // this one has charge EECs split into colors of RL bin
                    c_charge_all->cd();
                    if (j==0) {
                        c_charge_all->Divide(1,2);                    
                    }
                    c_charge_all->cd(1);
                    if (j==0) {
                        gPad->SetLogx();
                        // FormatHist(l2, hdummyRL, "", colors[j], markers[0]);
                        hdummyRL->SetMarkerColorAlpha(kBlue, 0);
                        hdummyRL->SetLineColorAlpha(kRed, 0);
                        hdummyRL->SetMinimum(0.);
                        hdummyRL->SetMaximum(3.5);
                        // hdummyRL->GetXaxis()->SetRangeUser(1e-4, 0.4); //(0.01, 0.4) //TODO: GOT RID OF FOR NOW
                        hdummyRL->Draw();
                    }
                    
                    hcorr_oppcharge_ENC_truth_arr[j]->Draw("same");
                    hcorr_samecharge_ENC_truth_arr[j]->Draw("same");

                    hcorr_oppcharge_ENC_truth_arr[j]->SetMarkerColorAlpha(colors[j], 0.6);
                    hcorr_oppcharge_ENC_truth_arr[j]->SetLineColorAlpha(colors[j], 0.6);
                    hcorr_samecharge_ENC_truth_arr[j]->SetMarkerColorAlpha(colors[j], 0.6);
                    hcorr_samecharge_ENC_truth_arr[j]->SetLineColorAlpha(colors[j], 0.6);
                    l_right->AddEntry(hcorr_oppcharge_ENC_truth_arr[j],label[j]);
                    // l_right->AddEntry(hcorr_samecharge_ENC_truth_arr[j],label[j]);
                    
                    // hcorr_oppcharge_truth_arr[j]->Draw("same");
                    // hcorr_samecharge_truth_arr[j]->Draw("same");

                    // hcorr_oppcharge_truth_arr[j]->SetMarkerColorAlpha(colors[j], 0.6);
                    // hcorr_oppcharge_truth_arr[j]->SetLineColorAlpha(colors[j], 0.6);
                    // hcorr_samecharge_truth_arr[j]->SetMarkerColorAlpha(colors[j], 0.6);
                    // hcorr_samecharge_truth_arr[j]->SetLineColorAlpha(colors[j], 0.6);
                    // l_right->AddEntry(hcorr_oppcharge_truth_arr[j],label[j]);
                    // // l_right->AddEntry(hcorr_samecharge_truth_arr[j],label[j]);
                    
                    l_right->Draw("same");
                    // gPad->BuildLegend();


                    c_charge_all->cd(2);
                    if (j==0) {
                        gPad->SetLogx();
                        hdummyRL2->SetMarkerColorAlpha(kBlue, 0);
                        hdummyRL2->SetLineColorAlpha(kRed, 0);
                        hdummyRL2->SetMinimum(0.);
                        hdummyRL2->SetMaximum(2.5);
                        // hdummyRL2->GetXaxis()->SetRangeUser(1e-4, 0.4); //(0.01, 0.4) //TODO: GOT RID OF FOR NOW
                        hdummyRL2->Draw();
                        drawHoriLine(0.01, 0.4, 1., kBlack);
                    }
                    hcorr_chargeratio_ENC_truth_arr[j]->Draw("same");

                    




                    
                    // cout << "about to format plots" << endl;
                    // FormatHist(l, hc, label1, markercolor1, markerstyle1);
                    // if (charmdecays) {
                    //     l->AddEntry("NULL","          D* decays on","h");
                    // } else {
                    //     l->AddEntry("NULL","          D* decays off","h");
                    // }
                    // FormatHist(l, hg, label2, markercolor2, markerstyle2);
                    // FormatHist(l, hl, label3, markercolor3, markerstyle3);
                    

                    // hc->Draw("L same");
                    // hg->Draw("L same");
                    // hl->Draw("L same");
                    
                


                    // make ratio plot
                    // auto rp = new TRatioPlot(hD0, hc);
                    // rp->Draw();
                    // rp->GetLowYaxis()->SetNdivisions(505);
                    // c->Update();


            


                    




                    // draw legend
                    // l->Draw("same");

                    save_del_canvas(c_deltap, outdir, "deltap", normstring_dir, weightstr, "", ptname, RLname, jetR, add_name, true);
                    save_del_canvas(c_deltapt, outdir, "deltapt", normstring_dir, weightstr, "", ptname, RLname, jetR, add_name, true);
                    save_del_canvas(c_deltapl, outdir, "deltapl", normstring_dir, weightstr, "", ptname, RLname, jetR, add_name, true);
                    if (normed != 1) save_del_canvas(c_charge, outdir, "charge", normstring_dir, weightstr, "", ptname, RLname, jetR, add_name, true);
                    if (normed != 1) save_del_canvas(c_chargeratio, outdir, "chargeratio", normstring_dir, weightstr, "", ptname, RLname, jetR, add_name, true);
                    save_del_canvas(c_energyweights, outdir, "energyweights", normstring_dir, weightstr, "", ptname, RLname, jetR, add_name, true);
                    
                    save_del_canvas(c_deltapt_vs_ew, outdir, "deltaptvsEW", normstring_dir, weightstr, "", ptname, RLname, jetR, add_name, true);


                    // delete expfit;

                    f_out->cd();
                    hcorr_deltap_truth->Write();
                    hcorr_deltapt_truth->Write();
                    hcorr_deltapl_truth->Write();
                    hcorr_oppcharge_truth->Write();
                    hcorr_samecharge_truth->Write();
                    // hcorr_unweightedRL_truth_unnorm->Write();
                    hcorr_chargeratio_ENC_truth->Write();
                    hcorr_energyweights_truth->Write();
                    hcorr_baryonmeson_truth->Write();

                    
                    delete hcorr_deltap_truth;
                    delete hcorr_deltapt_truth;
                    delete hcorr_deltapl_truth;
                    delete hcorr_oppcharge_truth;
                    delete hcorr_samecharge_truth;
                    // delete hcorr_unweightedRL_truth_unnorm;
                    delete hcorr_chargeratio_ENC_truth;
                    delete hcorr_energyweights_truth;
                    delete hcorr_baryonmeson_truth;

                    // delete h_pt_jetlevel_clone;
                    

                
                } // RL bins loop

                //----------------------------//
                /*c_charge_all->cd();
                for (int j=n_RLbins-1; j>=0; j--) {
                    
                    // c_charge_all->Divide(1,2);
                    // c_charge_all->cd(1);
                    // c_charge_all->cd(2);
                    // hcorr_chargeratio_truth_tt1_arr[j]->Draw("same");
                
                    if (j==4) {
                        hcorr_oppcharge_truth_tt1_arr[j]->Draw("same");
                        // hcorr_samecharge_truth_tt1_arr[j]->Draw("same");
                    }
                    if (j==4) {
                        hcorr_oppcharge_truth_tt1_arr[j]->SetMinimum(0.);
                        hcorr_oppcharge_truth_tt1_arr[j]->SetMaximum(2.); //65);
                        hcorr_oppcharge_truth_tt1_arr[j]->GetXaxis()->SetLimits(0.01, 0.4);
                        // hcorr_oppcharge_truth_tt1_arr[j]->SetAxisRange(0.01, 0.4,"X");
                        // hcorr_oppcharge_truth_tt1_arr[j]->GetXaxis()->SetRangeUser(0.01, 0.4);
                    }
                    l->Draw("same");
                }*/
                //----------------------------//

                // plot the combined graphs
                plot_combined_hists(c_deltap_all, hcorr_deltap_truth_arr, l_bottomleft, pt_max, true, RL_bin_width[i]);
                plot_combined_hists(c_deltapt_all, hcorr_deltapt_truth_arr, l_bottomleft, pt_max, true, RL_bin_width[i]);
                plot_combined_hists(c_deltapl_all, hcorr_deltapl_truth_arr, l_bottomleft, pt_max, true, RL_bin_width[i], plaxis_cut[i]);
                plot_combined_hists(c_energyweights_all, hcorr_energyweights_truth_arr, l, pt_max, false, RL_bin_width[i]);



                c_chargeratio_all->cd();
                TGraph *g = new TGraph(n_RLbins, RL_bin_centers[i], y_chargeratio_arr); 
                FormatGraph(g, l2, "", "R_{L}", "# opp sign pairs / # same sign pairs", 1.0, 20, kBlack, 0);
                g->Draw("AP");
                // std::vector<TGraph*> graph_inds;
                for (int j=0; j<n_RLbins; j++) {
                    double x_temp[1] = {RL_bin_centers[i][j]};
                    double y_temp[1] = {y_chargeratio_arr[j]};
                    // divide by the bin width??
                    // y_temp[0] /= RL_bin_width[j];
                    // TGraph *g_ind = new TGraph(1, &RL_bins[i][j], &y_chargeratio_arr[j]); //RL_bins[i] + j, y_chargeratio_arr + j); // gives what I thought would be TGraph(1, RL_bins[i][j], y_chargeratioarr[j]), but this syntax is wrong
                    TGraph *g_ind = new TGraph(1, x_temp, y_temp); //RL_bins[i] + j, y_chargeratio_arr + j); // gives what I thought would be TGraph(1, RL_bins[i][j], y_chargeratioarr[j]), but this syntax is wrong
                    FormatGraph(g_ind, l2, "", "R_{L}", "# opp sign pairs / # same sign pairs", 1.0, 20, colors[j], 1.0);
                    // graph_inds.push_back(g_ind);
                    g_ind->Draw("P same");
                    // cout << "the items here are " << RL_bins[i][j] << " vs " << RL_bins[i]+j << " and " << y_chargeratio_arr[j] << " vs " << y_chargeratio_arr + j << endl;
                }
                
                l->Draw("same");


                // do baryons + mesons
                c_baryonmeson_all->cd();
                TGraph *g_bb = new TGraph(n_RLbins, RL_bin_centers[i], y_baryonbaryon_arr); 
                TGraph *g_bm = new TGraph(n_RLbins, RL_bins[i], y_baryonmeson_arr); 
                TGraph *g_mm = new TGraph(n_RLbins, RL_bins[i], y_mesonmeson_arr); 
                // FormatGraph(g_bb, l2, "", "R_{L}", "# of pairs", 1.0, 20, kBlack, 0);
                FormatGraph(g_bb, l_bm, "baryon-baryon", "R_{L}", "# of pairs", 1.5, 20, kBlack, 1.0);
                FormatGraph(g_bm, l_bm, "baryon-meson", "R_{L}", "# of pairs", 1.5, kFullStar, kBlack, 1.0);
                FormatGraph(g_mm, l_bm, "meson-meson", "R_{L}", "# of pairs", 1.5, kFullSquare, kBlack, 1.0);
                    
                // FormatGraph(g_bm, "R_{L}", "# of pairs", 1.0, kFullStar, kBlack, 0); 
                // FormatGraph(g_mm, "R_{L}", "# of pairs", 1.0, kFullSquare, kBlack, 0);

                double max_bm = 0;
                for (int j=0; j<n_RLbins; j++) {
                    double max_bm_temp = hcorr_baryonmeson_truth_arr[j]->GetMaximum();
                    if (max_bm_temp > max_bm) max_bm = max_bm_temp;
                }
                g_bb->SetMaximum( max_bm * 1.2 );
                g_bb->Draw("AP");
                // g_bm->Draw("P same");
                // g_mm->Draw("P same");

                for (int j=0; j<n_RLbins; j++) {
                    double x_temp[1] = {RL_bin_centers[i][j]};
                    double y_bb_temp[1] = {y_baryonbaryon_arr[j]};
                    double y_bm_temp[1] = {y_baryonmeson_arr[j]};
                    double y_mm_temp[1] = {y_mesonmeson_arr[j]};
                    // TGraph *g_ind = new TGraph(1, &RL_bins[i][j], &y_chargeratio_arr[j]); //RL_bins[i] + j, y_chargeratio_arr + j); // gives what I thought would be TGraph(1, RL_bins[i][j], y_chargeratioarr[j]), but this syntax is wrong
                    TGraph *g_bb_ind = new TGraph(1, x_temp, y_bb_temp); //RL_bins[i] + j, y_chargeratio_arr + j); // gives what I thought would be TGraph(1, RL_bins[i][j], y_chargeratioarr[j]), but this syntax is wrong
                    TGraph *g_bm_ind = new TGraph(1, x_temp, y_bm_temp);
                    TGraph *g_mm_ind = new TGraph(1, x_temp, y_mm_temp);
                    FormatGraph(g_bb_ind, l2, "", "R_{L}", "# of pairs", 1.5, 20, colors[j], 1.0);
                    FormatGraph(g_bm_ind, l2, "", "R_{L}", "# of pairs", 1.5, kFullStar, colors[j], 1.0);
                    FormatGraph(g_mm_ind, l2, "", "R_{L}", "# of pairs", 1.5, kFullSquare, colors[j], 1.0);
                    g_bb_ind->Draw("P same");
                    g_bm_ind->Draw("P same");
                    g_mm_ind->Draw("P same");
                }
                l_bm->Draw("same");


                c_baryonmesonratios_all->cd();

                double max_bb_bm_mm = 0;
                for (int j=0; j<n_RLbins; j++) {
                    if (y_baryonbaryon_arr[j]/y_mesonmeson_arr[j] > max_bb_bm_mm) max_bb_bm_mm = y_baryonbaryon_arr[j]/y_mesonmeson_arr[j];
                    if (y_baryonbaryon_arr[j]/y_baryonmeson_arr[j] > max_bb_bm_mm) max_bb_bm_mm = y_baryonbaryon_arr[j]/y_baryonmeson_arr[j];
                    if (y_baryonmeson_arr[j]/y_mesonmeson_arr[j] > max_bb_bm_mm) max_bb_bm_mm = y_baryonmeson_arr[j]/y_mesonmeson_arr[j];
                }
                TGraph *g_bb_bm_mm_dummy = new TGraph(n_RLbins, RL_bin_centers[i], y_mesonmeson_arr); //doesn't matter what y values will be since will get scaled
                TGraph *g_bb_bm_mm = new TGraph(n_RLbins, RL_bin_centers[i], y_mesonmeson_arr); //doesn't matter what y values will be since will get scaled
                TGraph *g_bb_mm_mm = new TGraph(n_RLbins, RL_bin_centers[i], y_mesonmeson_arr); //doesn't matter what y values will be since will get scaled
                TGraph *g_bm_mm_mm = new TGraph(n_RLbins, RL_bin_centers[i], y_mesonmeson_arr); //doesn't matter what y values will be since will get scaled
                FormatGraph(g_bb_bm_mm_dummy, l2, "", "R_{L}", "# of pairs", 1.0, 20, kBlack, 0.0); //drawing transparently
                gStyle->SetLegendFont(30);
                FormatGraph(g_bb_bm_mm, l_bmratios, "#frac{baryon-baryon}{meson-meson}", "R_{L}", "# of pairs", 1.5, 20, kBlack, 1.0);
                FormatGraph(g_bb_mm_mm, l_bmratios, "#frac{baryon-baryon}{baryon-meson}", "R_{L}", "# of pairs", 1.5, kFullStar, kBlack, 1.0);
                FormatGraph(g_bm_mm_mm, l_bmratios, "#frac{baryon-meson}{meson-meson}", "R_{L}", "# of pairs", 1.5, kFullSquare, kBlack, 1.0);
                gStyle->SetLegendFont(42);
                g_bb_bm_mm_dummy->SetMaximum( max_bb_bm_mm * 1.2 );
                g_bb_bm_mm_dummy->Draw("AP");

                for (int j=0; j<n_RLbins; j++) {
                    double x_temp[1] = {RL_bin_centers[i][j]};
                    double y_bb_mm_temp[1] = {y_baryonbaryon_arr[j]/y_mesonmeson_arr[j]};
                    double y_bb_bm_temp[1] = {y_baryonbaryon_arr[j]/y_baryonmeson_arr[j]};
                    double y_bm_mm_temp[1] = {y_baryonmeson_arr[j]/y_mesonmeson_arr[j]};

                    TGraph *g_bb_mm_ind = new TGraph(1, x_temp, y_bb_mm_temp); //RL_bins[i] + j, y_chargeratio_arr + j); // gives what I thought would be TGraph(1, RL_bins[i][j], y_chargeratioarr[j]), but this syntax is wrong
                    TGraph *g_bb_bm_ind = new TGraph(1, x_temp, y_bb_bm_temp);
                    TGraph *g_bm_mm_ind = new TGraph(1, x_temp, y_bm_mm_temp);
                    FormatGraph(g_bb_mm_ind, l2, "", "R_{L}", "# of pairs", 1.5, 20, colors[j], 1.0);
                    FormatGraph(g_bb_bm_ind, l2, "", "R_{L}", "# of pairs", 1.5, kFullStar, colors[j], 1.0);
                    FormatGraph(g_bm_mm_ind, l2, "", "R_{L}", "# of pairs", 1.5, kFullSquare, colors[j], 1.0);
                    g_bb_mm_ind->Draw("P same");
                    g_bb_bm_ind->Draw("P same");
                    g_bm_mm_ind->Draw("P same");
                }
                // l_bmratios->Draw("same");

                


                // save combined plots as pdfs
                    
                save_del_canvas(c_deltap_all, outdir, "deltap", normstring_dir, weightstr, "", ptname, "", jetR, add_name, false);
                save_del_canvas(c_deltapt_all, outdir, "deltapt", normstring_dir, weightstr, "", ptname, "", jetR, add_name, false);
                save_del_canvas(c_deltapl_all, outdir, "deltapl", normstring_dir, weightstr, "", ptname, "", jetR, add_name, false);
                save_del_canvas(c_charge_all, outdir, "charge", normstring_dir, weightstr, "", ptname, "", jetR, add_name, false);
                if (normed != 1) save_del_canvas(c_chargeratio_all, outdir, "chargeratio", normstring_dir, weightstr, "", ptname, "", jetR, add_name, false);
                save_del_canvas(c_energyweights_all, outdir, "energyweights", normstring_dir, weightstr, "", ptname, "", jetR, add_name, false);
                if (normed != 1) save_del_canvas(c_baryonmeson_all, outdir, "baryonmeson", normstring_dir, weightstr, "", ptname, "", jetR, add_name, false);
                if (normed != 1) save_del_canvas(c_baryonmesonratios_all, outdir, "baryonmesonratios", normstring_dir, weightstr, "", ptname, "", jetR, add_name, false);
                save_del_canvas(c_deltap_all_fit05, outdir, "deltap", normstring_dir, weightstr, "_fit05", ptname, "", jetR, add_name, false);

               

            } // pT bins loop
        } // threshold loop
    } // jetR loop
}

// ======================================================= //
//                     MAIN FUNCTION 
// ======================================================= //

void old_make_corr_histograms() {

    gStyle->SetOptStat(0);
    SetStyle();
    

    // Files
    // const char infile[] = "/software/users/blianggi/mypyjetty/analysis/output/100k/AnalysisResultsFinal.root"; //hiccup
    const char infile[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/pythiagen/scaling/30174983/26652369/AnalysisResultsFinal.root"; //perlmutter, after june 2024
    // const char infile[] = "/Volumes/NO NAME/AnalysisResultsFinal.root"; //local

    // auto tree1 = new TChain("tree1");
    // for (int ifile=1; ifile<=10; file++) {
    //     std::string filename = indir + std::to_string(ifile) + "/AnalysisResultsFastSim_batch.root";
    //     tree1->Add(filename.c_str());
    // }


    TFile* f;
    std::string add_name;


    // CONTOL VARIABLES HERE
    // normed is 0 if unnormalized, 1 for self-normalization 
    // weighted is true if using "Weighted", false if using unweighted
    bool include_RL0 = false;


    f = new TFile(infile, "READ");
    add_name = "_othercorrel";
    cout << "output name will be " << add_name << endl;

    // first do unnormalized
    int normed = 0;
    // plot unweighted hists
    plot_histograms(f, add_name, normed, false, include_RL0);
    // plot weighted hists
    plot_histograms(f, add_name, normed, true, include_RL0);


    // now do self-normalized
    normed = 1;
    // plot unweighted hists
    plot_histograms(f, add_name, normed, false, include_RL0);
    // plot weighted hists
    plot_histograms(f, add_name, normed, true, include_RL0);


    // now do normalize by # of jets
    normed = 2;
    // plot unweighted hists
    plot_histograms(f, add_name, normed, false, include_RL0);
    // plot weighted hists
    plot_histograms(f, add_name, normed, true, include_RL0);



    f->Close();
    delete f;


    return;
}
