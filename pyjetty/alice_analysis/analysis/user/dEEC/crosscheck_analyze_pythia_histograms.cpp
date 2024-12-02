// ROOT macro make as a crosscheck
// it is currently able to plot 5 TeV pythia histograms (old method)
// Beatrice Liang-Gilman (beatrice_lg@berkeley.edu)

// global variables
Double_t colors[16] = {kGray, kMagenta, kGreen+2, kBlue, kOrange+1, kViolet+1, kRed, kYellow+1, kCyan+1};
Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
Double_t marker_size = 1.5;

std::string attempt_dir = "pythia5TeV_histograms_crosscheck";
std::string outdir = "/software/users/blianggi/mypyjetty/storage/dEEC/plots/" + attempt_dir;

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


void applyCuts(THnSparse *hsparse, int RLaxis, int pt_min, int pt_max, double RL_min, double RL_max, 
               bool EWaxis = false, bool BM_ENC = false, double bm_num = 0.) { //todo: don't need paxis anymore??

    hsparse->GetAxis(0)->SetRangeUser(pt_min, pt_max);
    // cout << "hello  " << BM_ENC << "//" << bm_num << endl;
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

    if (BM_ENC) {
        hsparse->GetAxis(4)->SetRangeUser(bm_num, bm_num+1); //i.e. if baryon-baryon, should get bin that is set at 1.
    }

    // cout << "CHECK 1A LOOKING AT NUM EMNRTIRES HERE: " << hsparse->GetEntries() << endl;

    
}

// get the observable histogram
//usually obsaxis is 3, but in new histograms it is 4. For jet level histograms, use 0.
TH1D * getObsHist(TFile *filename, std::string h_name, std::string h_jet_name, int pt_min, int pt_max, 
                  int RLaxis, double RL_min, double RL_max, std::string newhistname, int obsaxis, std::string xtitle,
                  int normalized=0, bool scalebyRLbinwidth=false, double RL_bin_width=1.0, bool EWaxis=false, 
                  bool BM_ENC=false, double bm_num=0., bool debug=false) { //int n_rebin_bins, double rebinbins[]
    
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

    applyCuts(hsparse, pt_min, pt_max, RLaxis, RL_min, RL_max, EWaxis, BM_ENC, bm_num);
    applyCuts(hsparse_jetlevel, pt_min, pt_max, -1, RL_min, RL_max, false, false, bm_num); //making no cuts on RL to keep it jet level??

    
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

/* get a typical 1D histogram */
// do we need jetlevelhist?, bool EWaxis=false, bool BM_ENC=false, double bm_num=0.
TH1D * getObs1DHist(TFile *filename, std::string h_name, int RLaxis, int obsaxis,
                    int pt_min, int pt_max, double RL_min, double RL_max, 
                    bool debug=false) {

    THnSparse *hsparse = getHistAndClone(filename, h_name);
    applyCuts(hsparse, RLaxis, pt_min, pt_max, RL_min, RL_max);
    TH1D *hist1D = hsparse->Projection(obsaxis);

    return hist1D;
}


void Format1DHist(TH1D *hist, TH1D *jetpt_hist, std::string norm_string, double x_left, double x_right,
                  int markercolor, double markeralpha, int markerstyle, std::string xtitle, std::string ytitle, 
                  TLegend& leg, TString leg_text, bool drawline=false, double linealpha=1., std::string obs_name="") {

    std::string new_name = std::string(hist->GetName()) + norm_string;
    hist->SetNameTitle(new_name.c_str(), new_name.c_str());
    // hist->SetName(Form("%s_%s", hist->GetName().c_str(), norm_string.c_str()));

    // set x range
    hist->GetXaxis()->SetRangeUser(x_left, x_right);
    
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

        for (int k=0; k < hist->GetNbinsX(); k++){
            hist->SetBinError(k+1, 0);
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




// get r_c
double getRcFromHists(TH1D *hist_charge) { //_oppsign, TH1D *hist_samesign) {
    // do i need to scale by the RL bin width here?? - I think this would be redundant.
    // if both like sign bin and unlike sign get scaled by RL bin width, then the ratio still stays the same

    // cout << "ENTRIES FOR SAME SIGH " << hist_samesign->GetEntries() << endl;
    // cout << "ENTRIES FOR OPP SIGH " << hist_oppsign->GetEntries() << endl;
    // cout << "BINS FOR SAME SIGH " << hist_samesign->GetNbinsX() << endl;
    // cout << "BINS FOR OPP SIGH " << hist_oppsign->GetNbinsX() << endl;

    // get # of like sign and # of unlike sign
    double num_likesign = hist_charge->GetBinContent(hist_charge->FindBin(1)); //GetBinContent(3);
    double num_unlikesign = hist_charge->GetBinContent(hist_charge->FindBin(-1)); //GetBinContent(1);
    // (if debug) cout << "num like sign " << num_likesign << " num unlike sign " << num_unlikesign << endl;

    // calculate the rc value for this pt & RL bin
    double rc = (double)(num_likesign - num_unlikesign) / (double)(num_likesign + num_unlikesign);
    //if (debug) cout << "and that makes rc " << rc << endl;
    if (num_likesign + num_unlikesign == 0) rc = 0;

    return rc;
}

// get r_c error
double getRcErrFromHists(TH1D *hist_charge) { //_oppsign, TH1D *hist_samesign) {
    // do i need to scale by the RL bin width here?? - I think this would be redundant.
    // if both like sign bin and unlike sign get scaled by RL bin width, then the ratio still stays the same

    // get # of like sign and # of unlike sign
    double num_likesign = hist_charge->GetBinContent(hist_charge->FindBin(1)); //GetBinContent(3);
    double num_unlikesign = hist_charge->GetBinContent(hist_charge->FindBin(-1)); //GetBinContent(1);
    // (if debug) cout << "num like sign " << num_likesign << " num unlike sign " << num_unlikesign << endl;

    // calculate the rc error for this pt & RL bin
    double num_totalpairs = num_likesign + num_unlikesign;
    double rc_err = ( 2 * sqrt( num_totalpairs * num_likesign * num_unlikesign ) ) / (num_totalpairs * num_totalpairs);

    // if (debug) cout << "RC ERR IN FUNC IS " << rc_err << endl;

    return rc_err;
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

void fix_up_2D(TH2D * hist2D, double bounds[], double ptmax, std::string xtitle, std::string ytitle) {

    hist2D->GetXaxis()->SetRangeUser(0, ptmax);
    hist2D->GetZaxis()->SetRangeUser(bounds[0], bounds[1]);

    // label axes
    hist2D->GetXaxis()->SetLabelFont(42);
    hist2D->GetXaxis()->SetTitleFont(42);
    hist2D->GetXaxis()->SetTitleSize(0.06); //(0.042);
    hist2D->GetXaxis()->SetTitleOffset(1.0);
	hist2D->GetXaxis()->SetLabelSize(0.05);
    hist2D->GetXaxis()->SetTitle(xtitle.c_str());

    hist2D->GetYaxis()->SetLabelFont(42);
	hist2D->GetYaxis()->SetTitleFont(42);
    // if (obs_name_y == "weights") {
    //     hist2D->GetYaxis()->SetTitleSize(0.035);
    //     hist2D->GetYaxis()->SetTitleOffset(1.5);
    // } else {
        hist2D->GetYaxis()->SetTitleSize(0.05); //0.06 //(0.042);
        hist2D->GetYaxis()->SetTitleOffset(1.0);
    // }
    hist2D->GetYaxis()->SetLabelSize(0.05); //(0.042);
    hist2D->GetYaxis()->SetTitle(ytitle.c_str());

}

// imported function from data
TGraphErrors * MakeFormatGraph(vector<double> xvals, vector<double> yvals, int markercolor, double markeralpha,
                  int markerstyle, std::string xtitle, std::string ytitle, std::string obs_name) {
    TGraphErrors * graph = new TGraphErrors(xvals.size(), xvals.data(), yvals.data());
    graph->SetTitle(Form("Charge Ratio;%s;%s", xtitle.c_str(), ytitle.c_str())); // Set the title and axis labels

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


// another function imported from data analysis code
// void plot_rc(vector<vector<double>>& RL_vals, vector<vector<double>>& rc_vals, vector<double>& ptcenter_bins,
//              Double_t *colors, Double_t *markers) {
//             //  TLegend& leg_RLbins, TLegend& leg_ptbins) {
    
//     vector<TGraph *> rc_graphs_func_of_RL;
//     vector<TGraph *> rc_graphs_func_of_pT;
//     vector<vector<TGraph*>> rc_graphs_func_of_RL_ind;
//     vector<vector<TGraph*>> rc_graphs_func_of_pT_ind;
//     // std::string outdir = "plots/thirdattempt"; // + ptbin_name + "/";//"plots/test/";
//     std::string fname_func_of_RL_out = outdir + "/corrhist_rc_func_of_RL.pdf"; // could add jetR and threshold info later??, maybe not needed tho 
//     std::string fname_func_of_pT_out = outdir + "/corrhist_rc_func_of_pT.pdf";

//     TLegend leg_RLbins(0.2, 0.6, 0.4, 0.85); 
//     TLegend leg_ptbins(0.5, 0.7, 0.65, 0.85); 

//     leg_RLbins.SetTextSize(0.037);
//     leg_RLbins.SetBorderSize(0);
//     leg_ptbins.SetTextSize(0.037);
//     leg_ptbins.SetBorderSize(0);

//     // delete hist;

//     // get graphs of r_c as a function of RL
//     // loop over pt bins
//     for ( int i = 0; i < ptcenter_bins.size(); i++ ) { 

//         // for graphs as a function of RL
//         TGraph *g = new TGraph(RL_vals[i].size(), RL_vals[i].data(), rc_vals[i].data());
//         g->SetMarkerStyle(markers[i]);
//         g->SetMarkerSize(1.5);
//         g->SetMarkerColorAlpha(kBlack, 1.0);
//         vector<TGraph*> ind_temp_vec;

//         for (int j=0; j<RL_vals[i].size(); j++){
//             int k=j+1;
//             TGraph *g_ind = new TGraph(1, &RL_vals[i][j], &rc_vals[i][j]);
//             g_ind->SetMarkerColorAlpha(colors[j], 1.0);
//             g_ind->SetMarkerSize(1.5);
//             g_ind->SetMarkerStyle(markers[i]);
//             ind_temp_vec.push_back(g_ind);
            
//             if (i==0) {
//                 leg_RLbins.AddEntry(g_ind, Form("R_{L} bin %d", j+1), "P");
//             }
//         }
//         rc_graphs_func_of_RL.push_back(g); 
//         rc_graphs_func_of_RL_ind.push_back(ind_temp_vec); 
        
//         // Should fix what gets subbed into %d so it is more flexible if the bins are not 20 GeV big?
//         leg_ptbins.AddEntry(g, Form("p_{T} = %d-%d", (int)ptcenter_bins[i]-10, (int)ptcenter_bins[i]+10), "P");
         

//     }

//     // plot r_c as a function of RL
//     TCanvas *can_func_of_RL = new TCanvas("can_func_of_RL", "can_func_of_RL", 750, 500);
//     can_func_of_RL->cd();
//     for ( int i = 0; i < ptcenter_bins.size(); i++ ) {
//         // rc_graphs_func_of_RL[i]->SetMarkerSize(1.0);
//         // rc_graphs_func_of_RL[i]->SetMarkerStyle(markers[i]);
//         if ( i == 0 ) {
//             rc_graphs_func_of_RL[i]->SetMinimum(-0.3);  // Lower y limit
//             rc_graphs_func_of_RL[i]->SetMaximum(0.1); 
//             rc_graphs_func_of_RL[i]->GetXaxis()->SetTitle("R_{L} bin center"); 
//             rc_graphs_func_of_RL[i]->GetYaxis()->SetTitle("r_{c}"); 
//             rc_graphs_func_of_RL[i]->Draw("AP");
//         } // else {
//         //     rc_graphs_func_of_RL[i]->Draw("P SAME");
//         // }

//         for (int j = 0; j < RL_vals[0].size(); j++) {
//             rc_graphs_func_of_RL_ind[i][j]->Draw("P SAME");
//         }
//     }
//     leg_RLbins.Draw("same");
//     leg_ptbins.Draw("same");
//     can_func_of_RL->SaveAs(fname_func_of_RL_out.c_str());
//     delete can_func_of_RL;

//     //========================================================

//     // get graphs of r_c as a function of pT
//     vector<vector<double>> rc_vals_func_of_pT;
//     // loop over RL bins
//     for ( int j = 0; j < RL_vals[0].size(); j++ ) {
//         vector<double> temp_vec;
//         vector<TGraph*> ind_temp_vec;

//         // save values into appropriate vectors
//         for ( int i = 0; i < ptcenter_bins.size(); i++ ) { 
//             temp_vec.push_back(rc_vals[i][j]);
            
//             int k=j+1;
//             TGraph *g_ind = new TGraph(1, &ptcenter_bins[i], &rc_vals[i][j]);
//             g_ind->SetMarkerColorAlpha(colors[j], 1.0);
//             g_ind->SetMarkerSize(1.5);
//             g_ind->SetMarkerStyle(markers[i]);
//             ind_temp_vec.push_back(g_ind);
//         }
//         rc_vals_func_of_pT.push_back(temp_vec);
//         rc_graphs_func_of_pT_ind.push_back(ind_temp_vec);  

//         // for graphs as a function of RL
//         TGraph *g = new TGraph(ptcenter_bins.size(), ptcenter_bins.data(), rc_vals_func_of_pT[j].data());
//         for (int aa = 0; aa < ptcenter_bins.size(); aa++) {
//             // cout << "studying pt=" << ptcenter_bins[aa] << " // " << rc_vals_func_of_pT[j][aa] << endl;
//         }
//         rc_graphs_func_of_pT.push_back(g); 
//     }
    
//     // plot r_c as a function of pT
//     TCanvas *can_func_of_pT = new TCanvas("can_func_of_pT", "can_func_of_pT", 750, 500);
//     can_func_of_pT->cd();
//     for ( int j = 0; j < RL_vals[0].size(); j++ ) {
//         int k = j+1;
//         // rc_graphs_func_of_pT[j]->SetMarkerSize(1.0);
//         // rc_graphs_func_of_pT[j]->SetMarkerStyle(markers[0]);
//         rc_graphs_func_of_pT[j]->SetMarkerColorAlpha(colors[j], 0.0);
//         if ( j == 0 ) {
//             rc_graphs_func_of_pT[j]->SetMinimum(-0.3);  // Lower y limit
//             rc_graphs_func_of_pT[j]->SetMaximum(0.1); 
//             rc_graphs_func_of_pT[j]->GetXaxis()->SetTitle("p_{T} bin center"); 
//             rc_graphs_func_of_pT[j]->GetYaxis()->SetTitle("r_{c}"); 
//             rc_graphs_func_of_pT[j]->Draw("AP");
//         } // else { 
//         //     rc_graphs_func_of_pT[j]->Draw("P SAME");
//         // }

//         for (int i = 0; i < ptcenter_bins.size(); i++) {
//             rc_graphs_func_of_pT_ind[j][i]->Draw("P SAME");
//         }
//     }
//     leg_RLbins.Draw("same");
//     leg_ptbins.Draw("same");
//     can_func_of_pT->SaveAs(fname_func_of_pT_out.c_str());
//     delete can_func_of_pT;



// }


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

// ====== TRASH FUNCTIONS SO FAR ====== 
// ====== KEEPING FOR POSTERITY =======
/*
void fitting() {
    //------------//------------//------------//------------//------------//
    // make fits
    // TF1 *fitFunc = new TF1("fitFunc", myFunction, 0, 10, 2); // 2 parameters
    // fitFunc->SetParameters(1, 0); // Initial guess for parameters
    // hist->Fit("fitFunc");

    // double p0 = fitFunc->GetParameter(0); // Get first parameter
    // double p0_err = fitFunc->GetParError(0); // Get error of the first parameter

    // graph->Draw("AP"); // Draw the graph with axis and points
    // fitFunc->Draw("same"); // Draw the fit function on the same canvas

    // exponential fit function
    // TF1 *expfit = new TF1("expfit", expFunc, 0, pt_max+5, 2); // 2 parameters
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

}


void plot_chargeotherthings() {
    double y_chargeratio_arr[n_RLbins];
    
    TH1D *hcorr_oppcharg_ENC_truth_unnorm = getObsHist(f, oppcharge_truth_name, jet_pt_truth_name, 
                                             pt_min, pt_max, i+1, RL_min, RL_max, "h_corr_oppcharge_ENC_Truth" + hist_addname, 
                                             i+1, "R_{L}", normed, false, RL_bin_width[i][j], false); //make normalization = 2?
                    
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

    // -----------------------------------------
    // do combined graphs now!
                    
    // this one has RL bin on the x axis, and y axis is #opp-sign-pairs/#same-sign-pairs
    // later will save to c_chargeratio_all
    double num_opp = hcorr_oppcharge_truth_arr[j]->GetBinContent(hcorr_oppcharge_truth_arr[j]->FindBin(-1));
    double num_same = hcorr_samecharge_truth_arr[j]->GetBinContent(hcorr_samecharge_truth_arr[j]->FindBin(1));
    y_chargeratio_arr[j] = num_opp/num_same;
    
    TLegend* l_right = new TLegend(0.75, 0.5, 0.85, 0.87);
    l_right->SetTextSize(0.037);
    l_right->SetBorderSize(0);
    
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
    l_right->Draw("same");  

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

    // make graphs
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

}

void plot_baryonsandmesons() {
    // do baryons + mesons

    TLegend* l_bm = new TLegend(0.63, 0.6, 0.83, 0.87);
    TLegend* l_bmratios = new TLegend(0.58, 0.5, 0.88, 0.87);
    
    TH1D *hcorr_baryonbaryon_ENC_truth_unnorm = getObsHist(f, baryonmeson_truth_name, jet_pt_truth_name, 
                                pt_min, pt_max, i+1, RL_min, RL_max, "h_corr_baryonbaryon_ENC_Truth" + hist_addname, 
                                i+1, "R_{L}", normed, false, RL_bin_width[i][j], false,
                                true, 1); //make normalization = 2?
    TH1D *hcorr_baryonmeson_ENC_truth_unnorm = getObsHist(f, baryonmeson_truth_name, jet_pt_truth_name, 
                                pt_min, pt_max, i+1, RL_min, RL_max, "h_corr_samecharge_baryonmeson_Truth" + hist_addname, 
                                i+1, "R_{L}", normed, false, RL_bin_width[i][j], false,
                                true, 0); //make normalization = 2?
    TH1D *hcorr_mesonmeson_ENC_truth_unnorm = getObsHist(f, baryonmeson_truth_name, jet_pt_truth_name, 
                                pt_min, pt_max, i+1, RL_min, RL_max, "h_corr_samecharge_mesonmeson_Truth" + hist_addname, 
                                i+1, "R_{L}", normed, false, RL_bin_width[i][j], false,
                                true, -1); //make normalization = 2?

    double y_baryonbaryon_arr[n_RLbins];
    double y_baryonmeson_arr[n_RLbins];
    double y_mesonmeson_arr[n_RLbins];

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

    // ============================================================
    c_bb_ENC->cd();
    hcorr_baryonbaryon_ENC_truth_unnorm->Draw();

    c_bm_ENC->cd();
    hcorr_baryonmeson_ENC_truth_unnorm->Draw();

    c_mm_ENC->cd();
    hcorr_mesonmeson_ENC_truth_unnorm->Draw();
    
    // -----------------------------------------
    
    //save baryon/meson pairs info
                    // 1 = p+p, 0 = p+pi, -1 = pi+pi, 2 = any other pairs
                    y_baryonbaryon_arr[j] = hcorr_baryonmeson_truth->GetBinContent(hcorr_baryonmeson_truth->FindBin(1));
                    y_baryonmeson_arr[j] = hcorr_baryonmeson_truth->GetBinContent(hcorr_baryonmeson_truth->FindBin(0));
                    y_mesonmeson_arr[j] = hcorr_baryonmeson_truth->GetBinContent(hcorr_baryonmeson_truth->FindBin(-1));
                    


}
*/
// ==================================== 


// ======================================================= //
//                   2ND MAIN FUNCTION 
// ======================================================= //

void analyze_ptbin(TFile * f_in, TFile * f_out, std::string weightstr, std::string jetR, std::string threshold, 
                   std::string norm_string, int pt_min, int pt_max, const double RL_bins[], int n_RLbins, 
                   vector<vector<double>>& RL_vals, vector<vector<double>>& rc_vals, vector<vector<double>>& rc_errors, 
                   bool include_RL0, bool include_RL1, bool debug, bool debug2) {

    std::string jetRname = Form("_R%s", jetR.c_str());
    std::string thrname = Form("_t%s", threshold.c_str());
    std::string ptname = to_string(pt_min) + "-" + to_string(pt_max);
    cout << "pt min " << pt_min << " and pt max " << pt_max << endl;
    
    for (int j = 0; j < n_RLbins; ++j) {
        cout << "RL j= " << j << " gives " << RL_bins[j] << endl;
    }
    
    double RL_bin_width[8] = {0}; 
    double RL_bin_centers[8] = {0};
    for (int j = 0; j < n_RLbins; ++j) {
        RL_bin_width[j] = RL_bins[j+1] - RL_bins[j];
        RL_bin_centers[j] = (RL_bins[j+1] + RL_bins[j])/2;
        // cout << "RL BIN WIDTH HERE" << RL_bin_width[i][j] << endl;
        // cout << " AND CENTERS " << RL_bin_centers[j] << endl;
    }

    std::string ytitle_norm = "#frac{1}{#DeltaR_{L}} ";
    if (norm_string == "self_normalized") ytitle_norm = "#frac{1}{N_{pair}#DeltaR_{L}} ";
    else if (norm_string == "norm_by_jets") ytitle_norm = "#frac{1}{N_{jet}#DeltaR_{L}} ";
    
    vector<TH1D*> deltap_vec;
    vector<TH1D*> deltapt_vec;
    vector<TH1D*> deltapl_vec;
    vector<TH1D*> weights_vec;
    // vector<TH1D*> q1q2_vec;
    vector<double> rc_vec;
    vector<double> rc_err_vec;
    vector<double> RLcenters_vec;

    TLegend *leg = new TLegend(0.6, 0.6, 0.85, 0.87);
    leg->SetTextSize(0.037);
    leg->SetBorderSize(0);
    TLegend *leg_dummy = new TLegend();

    // Names of histograms in the file
    std::string deltap_truth_name = Form("h_corr_deltap%s_JetPt_Truth_R%s_%sScaled", weightstr.c_str(), jetR.c_str(), threshold.c_str());
    std::string deltapt_truth_name = Form("h_corr_deltapt%s_JetPt_Truth_R%s_%sScaled", weightstr.c_str(), jetR.c_str(), threshold.c_str());
    std::string deltapl_truth_name = Form("h_corr_deltapl%s_JetPt_Truth_R%s_%sScaled", weightstr.c_str(), jetR.c_str(), threshold.c_str());
    std::string charge_truth_name = Form("h_corr_charge%s_JetPt_Truth_R%s_%sScaled", weightstr.c_str(), jetR.c_str(), threshold.c_str());
    // std::string unweightedRL_truth_name = Form("h_corr_unweightedRL%s_JetPt_Truth_R%s_%sScaled", weightstr.c_str(), jetR.c_str(), threshold.c_str()); //not implemented yet
    std::string energyweights_truth_name = Form("h_corr_energyweights%s_JetPt_Truth_R%s_%sScaled", weightstr.c_str(), jetR.c_str(), threshold.c_str());
    // std::string baryonmeson_truth_name = Form("h_corr_baryonmeson%s_JetPt_Truth_R%s_%sScaled", weightstr.c_str(), jetR.c_str(), threshold.c_str()); //not implemented yet
    
    // std::cout << "THRESHOLD NAME TEST" << deltap_truth_name << std::endl;
    
    std::string jet_pt_truth_name = Form("h_JETINFOjet_pt_Truth_R%s_%sScaled", jetR.c_str(), threshold.c_str());
    // std::string deltaptvsEW_truth_name = Form("h_ptvsenergyweights_JetPt_Truth_R%s_%sScaled", jetR.c_str(), threshold.c_str()); //not implemented yet

    
    for ( int j = 0; j < n_RLbins; j++ ) {
        int k = j;
        if (!include_RL0) {
            k = j-1;
            if (j == 0) continue; // can add something here to change the filename for ALL
        }
        if (!include_RL1 && j == n_RLbins-1) continue;

        int i = 0;
        if (pt_min == 40) i = 1;
        else if (pt_min == 60) i = 2;
        
        double RL_min = RL_bins[j];
        double RL_max = RL_bins[j+1];
        std::string RLname = Form("_RL%.3f-%.3f", RL_min, RL_max);
        std::string RLname_leg = Form("RL = %.3f-%.3f", RL_min, RL_max);
        // if (debug) cout << " in RL bin" << j << " with " << RL_min << " - " << RL_max << endl;
        
        std::string hist_addname = weightstr + jetRname + thrname + "_pt" + ptname + RLname;
        

        // get histograms
        // get jet pT range - no D0 reconstruction, so don't make D0 cuts
        // thnsparse axes: 0=jet pt, 1=RL (20 < pt < 40), 2=RL (40 < pt < 60), 3=RL (60 < pt < 80), 4 = observable 
        TH1D *hcorr_jetpt_inptbin_hist = getObs1DHist(f_in, jet_pt_truth_name, -1, 0, pt_min, pt_max, RL_min, RL_max);
        TH1D *hcorr_deltap_truth_hist = getObs1DHist(f_in, deltap_truth_name, i+1, 4, pt_min, pt_max, RL_min, RL_max);
        cout << "checkpoint 1 " << hcorr_deltap_truth_hist->GetEntries() << endl;
        TH1D *hcorr_deltapt_truth_hist = getObs1DHist(f_in, deltapt_truth_name, i+1, 4, pt_min, pt_max, RL_min, RL_max);
        TH1D *hcorr_deltapl_truth_hist = getObs1DHist(f_in, deltapl_truth_name, i+1, 4, pt_min, pt_max, RL_min, RL_max);
        TH1D *hcorr_energyweights_truth_hist = getObs1DHist(f_in, energyweights_truth_name, i+1, 4, pt_min, pt_max, RL_min, RL_max);
        // TH1D *hcorr_baryonmeson_truth_hist = getObs1DHist(f_in, baryonmeson_truth_name, i+1, 4, pt_min, pt_max, RL_min, RL_max);
        // TH2D *hcorr_deltaptvsEW_truth = get2DHist(f_in, deltaptvsEW_truth_name, jet_pt_truth_name,
                    //                            pt_min, pt_max, i+1, RL_min, RL_max, "h_corr_deltaptvsEW" + hist_addname,
                    //                            4, 5, "#Delta p_{T}", "p_{T,1}p_{T,2} / p_{T, jet}^{2}", normed, true, 
                    //                            RL_bin_width[i][j]);
                    

        
        double rc_value = 0.0;
        double rc_err = 0.0;
                    
        if (norm_string == "unnormalized") {
            TH1D *hcorr_charge_truth_hist = getObs1DHist(f_in, charge_truth_name, i+1, 4, pt_min, pt_max, RL_min, RL_max);
            rc_value = getRcFromHists(hcorr_charge_truth_hist); //hcorr_oppcharge_truth, hcorr_samecharge_truth);
            rc_err = getRcErrFromHists(hcorr_charge_truth_hist); //PAIRINFO_tree, "rc", 6, -3, 3, pt_min, pt_max, RL_min, RL_max);
            cout << "RC ERR IS " << rc_err << "(pt_min=" << pt_min << ", j=" << j << ")" <<endl;

            rc_vec.push_back(rc_value);
            rc_err_vec.push_back(rc_err);
            RLcenters_vec.push_back( (RL_min+RL_max)/2 );
        }

        // CURRENTLY NOT IMPLEMENTED
        // int nbins_2D = 50;
        // if (pt_min == 40 || pt_min == 60) nbins_2D = 30;
        // TH2D * weights_vs_deltapt_hist2D = getObs2DHistFromTChain(PAIRINFO_tree, "deltapt", "weights", nbins_2D, 0, pt_max+5, nbins_2D, 0, 0.3, pt_min, pt_max, RL_min, RL_max);


        // push to vectors
        deltap_vec.push_back((TH1D*) hcorr_deltap_truth_hist->Clone(hcorr_deltap_truth_hist->GetName()));
        deltapt_vec.push_back((TH1D*) hcorr_deltapt_truth_hist->Clone(hcorr_deltapt_truth_hist->GetName()));
        deltapl_vec.push_back((TH1D*) hcorr_deltapl_truth_hist->Clone(hcorr_deltapl_truth_hist->GetName()));
        weights_vec.push_back((TH1D*) hcorr_energyweights_truth_hist->Clone(hcorr_energyweights_truth_hist->GetName()));
        

        // format histograms in vector
        Format1DHist(deltap_vec[k], hcorr_jetpt_inptbin_hist, norm_string, 0, pt_max+5, colors[j], 0.6, markers[0], "#Deltap", ytitle_norm + "#frac{dN}{d#Deltap}", *leg, RLname_leg, true, 1.0);
        Format1DHist(deltapt_vec[k], hcorr_jetpt_inptbin_hist, norm_string, 0, pt_max+5, colors[j], 0.6, markers[0], "#Deltap_{T}", ytitle_norm + "#frac{dN}{d#Deltap_{T}}", *leg_dummy, RLname_leg, true, 1.0);
        Format1DHist(deltapl_vec[k], hcorr_jetpt_inptbin_hist, norm_string, 0, pt_max/2, colors[j], 0.6, markers[0], "#Deltap_{L}", ytitle_norm + "#frac{dN}{d#Deltap_{L}}", *leg_dummy, RLname_leg, true, 1.0);
        Format1DHist(weights_vec[k], hcorr_jetpt_inptbin_hist, norm_string, 0, 0.3, colors[j], 0.6, markers[0], "#frac{p_{T,1}p_{T,2}}{p_{T,jet}^{2}}", ytitle_norm + "#frac{dN}{d[EW]}", *leg_dummy, RLname_leg, true, 1.0, "weights");
        
        // Format2DHist(weights_vs_deltapt_hist2D, hcorr_jetpt_inptbin_hist, norm_string, ytitle_norm + "#Deltap_{T}", ytitle_norm + "p_{T,1}p_{T,2} / p_{T,jet}^{2}", true, RL_bin_width[j], "deltapt", "weights");

        // draw, save, and delete histograms
        TCanvas *can_deltap = new TCanvas();
        TCanvas *can_deltapt = new TCanvas();
        TCanvas *can_deltapl = new TCanvas();
        TCanvas *can_weights = new TCanvas();
        
        // TCanvas *can_weights_vs_deltapt = new TCanvas("can_weights_vs_deltapt", "can_weights_vs_deltapt", 800, 500);

        draw_save_del_hists(f_out, can_deltap, deltap_vec[k], "deltap", ptname, norm_string, hist_addname, false, true);
        draw_save_del_hists(f_out, can_deltapt, deltapt_vec[k], "deltapt", ptname, norm_string, hist_addname, false, true);
        draw_save_del_hists(f_out, can_deltapl, deltapl_vec[k], "deltapl", ptname, norm_string, hist_addname, false, false);
        draw_save_del_hists(f_out, can_weights, weights_vec[k], "weights", ptname, norm_string, hist_addname, false, true);
        
        // draw_save_del_hists(f_out, can_weights_vs_deltapt, weights_vs_deltapt_hist2D, "weights_vs_deltapt", ptname, norm_string, hist_addname, false, false, true);
        
    }

    /* do pt bin stuff here */
    std::string hist_all_addname = weightstr + jetRname + thrname + "_pt" + ptname;

	// combine RL plots to get 1 plot per pt bin

    TCanvas *can_deltap_all = new TCanvas();
    TCanvas *can_deltapt_all = new TCanvas();
    TCanvas *can_deltapl_all = new TCanvas();
    TCanvas *can_weights_all = new TCanvas();

    // size_t length_deltap = deltap_vec.size();
    // cout << " LENGTH DELTA P " << length_deltap << endl;
	
    plotandsave_combined_hists(can_deltap_all, deltap_vec, leg, "deltap", ptname, norm_string, hist_all_addname, pt_max, true, RL_bin_width, false, true, -1);
    plotandsave_combined_hists(can_deltapt_all, deltapt_vec, leg, "deltapt", ptname, norm_string, hist_all_addname, pt_max, true, RL_bin_width, false, true, -1);
    plotandsave_combined_hists(can_deltapl_all, deltapl_vec, leg, "deltapl", ptname, norm_string, hist_all_addname, pt_max, true, RL_bin_width, false, true, -1);
    plotandsave_combined_hists(can_weights_all, weights_vec, leg, "weights", ptname, norm_string, hist_all_addname, pt_max, true, RL_bin_width, false, true, -1);

    // make graphs
    if (norm_string == "unnormalized") {
        TCanvas *can_rc = new TCanvas();
        ProcessCanvas(can_rc);
        TGraphErrors *gr_rc = MakeFormatGraph(RLcenters_vec, rc_vec, kBlack, 1.0, markers[0], "R_{L}", "r_{c}", "rc");
        draw_save_del_hists(f_out, can_rc, gr_rc, "rc", ptname, norm_string, hist_all_addname, false, false);
        
        
        // save vectors here
        RL_vals.push_back(RLcenters_vec);
        rc_vals.push_back(rc_vec);
        rc_errors.push_back(rc_err_vec);
    }

}

// void plot_histograms(TFile* f, std::string add_name, int normed, bool weighted, bool include_RL0) {
void plot_histograms(TFile* f_in, TFile* f_out, bool weighted, std::string norm_string, 
                     const int pt_bins[], int n_bins, const double RL_bins[][8], int n_RLbins, 
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
            std::string jet_pt_truth_name = Form("h_JETINFOjet_pt_Truth_R%s_%sScaled", jetR.c_str(), threshold.c_str());

            if (norm_string == "unnormalized") {
                TH1D * jetpt_hist = getObs1DHist(f_in, jet_pt_truth_name, -1, 0, 0, 200, -1, -1, true);
                jetpt_hist->GetXaxis()->SetTitle("p_{T,jet}");
                TCanvas *can_jetpt = new TCanvas();
                draw_save_del_hists(f_out, can_jetpt, jetpt_hist, "jet_pt", "", "", weightstr + jetRname + thrname, false, true);
            
                TH1D * jet_const = getObs1DHist(f_in, jet_pt_truth_name, -1, 1, 0, 20, -1, -1, true);
                jet_const->GetXaxis()->SetTitle("Number Constituents (total)");
                TCanvas *can_numconst = new TCanvas();
                draw_save_del_hists(f_out, can_numconst, jet_const, "total_num_const", "", "", weightstr + jetRname + thrname, false, true);
            
                TH1D * jet_const_aftercut = getObs1DHist(f_in, jet_pt_truth_name, -1, 2, 0, 20, -1, -1, true);
                jet_const_aftercut->GetXaxis()->SetTitle("Number Constituents (after threshold cut)");
                TCanvas *can_numconst_aftercut = new TCanvas();
                draw_save_del_hists(f_out, can_numconst_aftercut, jet_const_aftercut, "num_const_aftercut", "", "", weightstr + jetRname + thrname, false, true);
            
            }

            //variables
            vector<vector<double>> RL_vals;
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
                analyze_ptbin(f_in, f_out, weightstr, jetR, threshold, norm_string, pt_min, pt_max, RL_bins[i], n_RLbins, RL_vals, rc_vals, rc_errors, include_RL0, include_RL1, debug, debug2);
        

                
            } // pT bins loop

            // plot r_c as a function of RL
            if (norm_string == "unnormalized") {
                cout << "checkpoint 4" << endl;
                cout << "size of RL_vals " << RL_vals.size() << endl;
                cout << "size of RL_vals[0] " << RL_vals[0].size() << endl;
                cout << "size of rc_vals " << rc_vals.size() << endl;
                cout << "size of rc_vals[0] " << rc_vals[0].size() << endl;
                cout << "size of ptcenter_bins " << ptcenter_bins.size() << endl;

                plot_rc(RL_vals, rc_vals, ptcenter_bins, rc_errors); //, colors, markers); //, leg_RLbins, leg_ptbins);
            }

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
    const char infile[] = "/rstorage/alice/AnalysisResults/blianggi/dEEC/445125/1132588/scaling/AnalysisResultsFinal.root"; //hiccup
    // const char infile[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/pythiagen/scaling/30174983/26652369/AnalysisResultsFinal.root"; //perlmutter, after june 2024
    // const char infile[] = "/Volumes/NO NAME/AnalysisResultsFinal.root"; //local
    TFile* root_infile = new TFile(infile, "READ");

    // Output file for binned results
    std::string outfile = "/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/" + attempt_dir + "/PYTHIAHists.root"; //plots/ntuples/DataHists.root"; //FinalDataHists.root
    TFile* root_outfile = new TFile(outfile.c_str(), "RECREATE");
    // std::string add_name = ""; // "_othercorrel";
    // cout << "output name will be " << add_name << endl;


    // analysis variables
    const int pt_bins[] = { 20, 40, 60, 80 };
    const int n_bins = sizeof(pt_bins) / sizeof(pt_bins[0]) - 1; //3;
    
    const double RL_bins[3][8] = { { 0, 1e-2, 3e-2, 7e-2, 1.5e-1, 3e-1, 4e-1, 1 },
                            { 0, 1e-2, 2.5e-2, 4e-2, 8e-2, 2.5e-1, 4e-1, 1 },
                            { 0, 1e-2, 2.5e-2, 3e-2, 4.5e-2, 2e-1, 4e-1, 1 } };
    const int n_RLbins = sizeof(RL_bins[0]) / sizeof(RL_bins[0][0]) - 1; //gets the columns //6; //7; //5;


    if (debug2) cout << "pt_bins " << n_bins << " n_RLbins " << n_RLbins << endl;

            
    // ====================================================================================


    // analyze for plots
    norm_string = "unnormalized";
    plot_histograms(root_infile, root_outfile, weighted, norm_string, pt_bins, n_bins, RL_bins, n_RLbins, include_RL0, include_RL1, debug, debug2);
    // analyze(root_infile, root_outfile, weightstr, jetRname, thrname, norm_string, pt_bins, n_bins, RL_bins, n_RLbins, include_RL0, include_RL1, debug, debug2);
    
    norm_string = "self_normalized";
    plot_histograms(root_infile, root_outfile, weighted, norm_string, pt_bins, n_bins, RL_bins, n_RLbins, include_RL0, include_RL1, debug, debug2);
        
    norm_string = "norm_by_jets";
    plot_histograms(root_infile, root_outfile, weighted, norm_string, pt_bins, n_bins, RL_bins, n_RLbins, include_RL0, include_RL1, debug, debug2);
    


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
