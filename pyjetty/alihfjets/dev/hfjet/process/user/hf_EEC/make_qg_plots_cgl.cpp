// ROOT macro to make quark-gluon jet plots
// this file is going to have EVERYTHING (g, l, i, c, D)
// the D will be no D* 
// all the latest files
// this one specifically for plotting&saving PT*RL charm-/gluon-/light-init jets from THnSparse
// Beatrice Liang-Gilman (beatrice_lg@berkeley.edu)

void SetStyle(Bool_t graypalette=true) {
  	cout << "Setting style!" << endl;
  
  	gStyle->Reset("Plain");
  	gStyle->SetOptTitle(0);
  	gStyle->SetOptStat(0);
  	if(graypalette) gStyle->SetPalette(8,0);
  	else gStyle->SetPalette(1);
  	gStyle->SetCanvasColor(10);
  	gStyle->SetCanvasBorderMode(0);
  	gStyle->SetFrameLineWidth(1);
  	gStyle->SetFrameFillColor(kWhite);
  	gStyle->SetPadColor(10);
  	gStyle->SetPadTickX(1);
  	gStyle->SetPadTickY(1);
  	gStyle->SetPadBottomMargin(0.15);
 	gStyle->SetPadLeftMargin(0.175);
    // gStyle->SetPadRightMargin(0.1);
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

void ProcessCanvas(TCanvas *Canvas) { 
	gStyle->SetOptStat(0);
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


void FormatHist(TLegend *l, TH1 *hist, TString text, int markercolor=1, int markerstyle=8, double markeralpha=1.,
                double xtitlesize=0.06, double xlabelsize=0.05, double xoffset=1.0,
                double ytitlesize=0.06, double ylabelsize=0.05, double yoffset=1.05, double markersize=1.5) 
                // double xtitlesize=0.04, double xlabelsize=0.04, double xoffset=1.2,
                // double ytitlesize=0.04, double ylabelsize=0.04, double yoffset=1.0) 
{
    hist->SetLineColor(markercolor);
    hist->SetMarkerColorAlpha(markercolor, markeralpha);
    hist->SetMarkerStyle(markerstyle);
    hist->SetMarkerSize(markersize);
    l->AddEntry(hist, text, "pl");

	//gPad->SetTickx(); 
	//gPad->SetTicky(); 
	// h->SetLineWidth(2);
	hist->GetYaxis()->SetTitleOffset(yoffset); //(1.05); 
	hist->GetYaxis()->SetTitleSize(ytitlesize); //(0.042); //the axis number labels
	hist->GetYaxis()->SetLabelSize(ylabelsize); //(0.042);
	hist->GetYaxis()->SetLabelFont(42);
	hist->GetXaxis()->SetLabelFont(42);
	hist->GetYaxis()->SetTitleFont(42);
	hist->GetXaxis()->SetTitleFont(42);
	hist->GetXaxis()->SetTitleOffset(xoffset);
	hist->GetXaxis()->SetTitleSize(xtitlesize); //(0.042);
	hist->GetXaxis()->SetLabelSize(xlabelsize); //(0.042);


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



//get histogram and clone it
THnSparse * getHistAndClone(TFile *f, std::string histname) {
    THnSparse *hsparsejet = (THnSparse*) f->Get(histname.c_str());
    THnSparse *hsparsejet_clone = (THnSparse *) hsparsejet->Clone("hsparsejet_c_clone"); //TODO: change this name!

    return hsparsejet_clone;
    
}


void applyCuts(THnSparse *hsparse, int pt_min, int pt_max, int d0_pt_cut, bool d0cuts=false) {
    hsparse->GetAxis(0)->SetRangeUser(pt_min, pt_max);
    if (d0cuts) {
        hsparse->GetAxis(1)->SetRangeUser(d0_pt_cut, pt_max); //d0_pt_cuts[i], pt_max); // apply cut on Dmeson pt
        hsparse->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
    }
}

// get the observable histogram
//usually obsaxis is 3, but in new histograms it is 4. For jet level histograms, use 0.
TH1D * getObsHist(TFile *filename, std::string h_name, std::string h_jet_name, int pt_min, int pt_max, int d0_pt_cut, std::string newhistname, bool d0cuts=false, int obsaxis=3, bool ptrl=false) {
    
    THnSparse *hsparse = getHistAndClone(filename, h_name);
    THnSparse *hsparse_jetlevel = getHistAndClone(filename, h_jet_name);

    applyCuts(hsparse, pt_min, pt_max, d0_pt_cut, d0cuts);
    applyCuts(hsparse_jetlevel, pt_min, pt_max, d0_pt_cut, d0cuts);
    
    TH1D *h_proj = hsparse->Projection(obsaxis);
    TH1D *h_proj_jetlevel = hsparse_jetlevel->Projection(0); // jet pt axis

    std::string hname = h_proj->GetName();
    hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
    h_proj->SetNameTitle(hname.c_str(), hname.c_str());

    // allow rebin or cloning here
    TH1D* hist = (TH1D*) h_proj->Clone(newhistname.c_str());

    //normalize
    double numjets = h_proj_jetlevel->Integral();
    hist->Scale(1/numjets, "width");
    // cout << "There are " << h_proj->GetEntries() << " pair entries in this pt bin" << endl;
    // cout << "There are " << h_proj_jetlevel->GetEntries() << " jet entries in this pt bin" << endl;

    if (ptrl) {
        hist->GetXaxis()->SetTitle("#it{p}_{T}#it{R}_{L}");
    } else {
        hist->GetXaxis()->SetTitle("#it{R}_{L}");
    }
    hist->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");

    return hist;
}

double drawWithLineAtMax(TH1D *hist, TLegend *l, std::string labeltext, int markercolor, int markerstyle, double lowx,
                    bool ptrl=false, bool removexaxis=false, double markeralpha=1, double xoffset=1.0, bool verbosity=false) {
    // FormatHist(l, hist, "gluon-init jets", markercolor_g, markerstyle_g);
    FormatHist(l, hist, labeltext, markercolor, markerstyle, markeralpha,
               0.06, 0.05, xoffset, 0.06, 0.05, 1.05);
    hist->GetXaxis()->SetRangeUser(lowx, 1.);
    hist->Draw("L same");
    if (removexaxis) hist->GetXaxis()->SetLabelSize(0);
    
    double hist_top_binpos = findTopOfCurve(hist, ptrl);
    drawVertLine(hist->GetBinCenter(hist_top_binpos), 0, hist->GetBinContent(hist_top_binpos), markercolor, 1)->Draw();

    if (verbosity) {
        cout << hist->GetName() << " top bincenter " << hist->GetBinCenter(hist_top_binpos) << ", " << hist->GetName() << " top binpos" << hist->GetBinContent(hist_top_binpos) << endl;
    }   

    return hist->GetBinCenter(hist_top_binpos);
}

void drawNoLine(TH1D *hist, TLegend *l, std::string labeltext, int markercolor, int markerstyle, double lowx,
                  float markersize, bool removexaxis=false, double markeralpha=1, double xoffset=1.0, bool verbosity=false) {
    // FormatHist(l, hist, "gluon-init jets", markercolor_g, markerstyle_g);
    FormatHist(l, hist, labeltext, markercolor, markerstyle, markeralpha,
               0.06, 0.05, xoffset, 0.06, 0.05, 1.05, markersize);
    hist->GetXaxis()->SetRangeUser(lowx, 1.);
    hist->Draw("L same");
    if (removexaxis) hist->GetXaxis()->SetLabelSize(0);
    
    if (verbosity) {
        cout << "no prints here!" << endl;
    }   

    return;
}

double getPeak(TH1D *hist, bool ptrl=false, bool verbosity=false) {
    
    double hist_top_binpos = findTopOfCurve(hist, ptrl);
    if (verbosity) {
        cout << hist->GetName() << " top bincenter " << hist->GetBinCenter(hist_top_binpos) << ", " << hist->GetName() << " top binpos" << hist->GetBinContent(hist_top_binpos) << endl;
    }   

    return hist->GetBinCenter(hist_top_binpos);
}

double getPeakErr(TH1D *hist, bool ptrl=false, bool verbosity=false) {
    
    double hist_top_binpos = findTopOfCurve(hist, ptrl);
    if (verbosity) {
        cout << hist->GetName() << " top bincenter " << hist->GetBinCenter(hist_top_binpos) << ", " << hist->GetName() << " top binpos" << hist->GetBinContent(hist_top_binpos) << endl;
    }   

    return hist->GetBinWidth(hist_top_binpos)/2;
}

void plotGraph(TMultiGraph *mg, TLegend *l, int n_bins, double *pt_x, double *peaks_y, double *err_y, TString text, int markercolor, int markerstyle, double markeralpha=1.0) {
    // TGraph *g = new TGraph(n_bins,pt_x,peaks_y);
    double err_x[n_bins];
    for (int i=0; i<n_bins; i++) {
        err_x[i] = 0;
    }
    TGraphErrors *g = new TGraphErrors(n_bins,pt_x,peaks_y, err_x, err_y);
    g->SetMarkerColorAlpha(markercolor, markeralpha);
    g->SetMarkerStyle(markerstyle);
    g->SetMarkerSize(1.5);
    g->SetLineColorAlpha(markercolor, markeralpha);
    // if (markerstyle == kFullDiamond) g->SetMarkerSize(1.75);
    if (markercolor == kRed && markerstyle == kOpenCircle) g->SetMarkerSize(1);

    l->AddEntry(g, text, "pl");

    g->GetXaxis()->SetTitle("#it{p}_{T}");
    g->GetYaxis()->SetTitle("#it{R}_{L} peak position");
    g->SetMaximum(5.);

    if (text == "#frac{b-init full jets}{l-init full jets}") {
        cout << "here!!!!!" << endl;
        g->SetMaximum(5.);
    }
    
   
    // g->Draw("ap SAME");
    // return g;
    mg->Add(g);

}

// TPad * makeTopPad(TH1D *hdummy, double ymax, double botmarg=0.31) {
TPad * makeTopPad(double botmarg=0.31) {
    TPad *pad1 = new TPad("pad1","pad1",0.,0.,1.,1.);
    pad1->SetLogx();
    pad1->SetFillColor(0);
    pad1->SetFillStyle(0);
    pad1->SetTopMargin(0.025);
    pad1->SetBottomMargin(botmarg);
    pad1->Draw();
    pad1->cd();

    // hdummy->SetMaximum(ymax);
    // hdummy->Draw();

    return pad1;
}

TPad * makeBottomPad(double topmarg=0.71, double botmarg=0.975) {
    TPad *pad2 = new TPad("pad1","",0.,0.,1.,1.);
    pad2->SetTopMargin(topmarg); //0.71);
    pad2->SetBottomMargin(botmarg);
    pad2->SetFillColor(0);
    pad2->SetFillStyle(0);
    pad2->Draw();
    pad2->SetLogx();
    // pad2->SetGridy();
    pad2->cd();

    return pad2;
}

TPad * plotRatio(TH1D *h1, TH1D *h2, std::string ratio_name, TLegend *l, int markercolor, int markerstyle, std::string yaxislabel, 
               double ymin=0.5, double ymax=1.5, double topmarg=0.71, double botmarg=0.975, double lowx=1e-3, std::string legendlabel="ratio",
               bool removexaxis=false, bool drawlineatone=false, bool justratio=false, int ndiv=5) {   

    TPad *pad2;
    if (!justratio) pad2 = makeBottomPad(topmarg, botmarg);

    TH1D* hratio = (TH1D*) h1->Clone(ratio_name.c_str());
    hratio->Divide(h2);
    hratio->SetMinimum(ymin);
    hratio->SetMaximum(ymax);
    
    //xtitlesize, xlabelsize, xoffset, ytitlesize, ylabelsize, yoffset
    FormatHist(l, hratio, legendlabel, markercolor, markerstyle, 1.0, 0.05, 0.04, 1.2, 0.035, 0.03, 1.5);
    hratio->GetYaxis()->SetNdivisions(5);
    if (removexaxis) hratio->GetXaxis()->SetLabelSize(0);
    else hratio->GetXaxis()->SetTitleOffset(0.95);
    
    hratio->GetXaxis()->SetRangeUser(lowx,1);

    hratio->Draw();

    hratio->GetYaxis()->SetTitleSize(0);

    if (!justratio) {
        double ypos = ((1-topmarg) + botmarg) / 2;
        TLatex *t = new TLatex(0.075,ypos,yaxislabel.c_str());
        t->SetTextAlign(22);
        t->SetTextColor(kBlack);
        t->SetTextFont(43);
        t->SetTextSize(14);
        // t->SetTextAngle(45);
        t->SetNDC(kTRUE);
        t->Draw();

        // draw line at 1
        
        if (drawlineatone) drawHoriLine(lowx, 1., 1., kGray+2, 1)->Draw();
        // drawHoriLine(1e-3, 1., 0.9, kGray+2, 6)->Draw();
        // drawHoriLine(1e-3, 1., 1.1, kGray+2, 6)->Draw();
    }

    return pad2;

}

// this function specifically plotting 2 ratios - one that is h1/h3, and one that is h2/h3
TPad * plotRatio2(TH1D *h1, TH1D *h2, TH1D *h3, std::string ratio_name1, std::string ratio_name2, 
               int markercolor1, int markerstyle1, double markersize1, double markeralpha1,
               int markercolor2, int markerstyle2, double markersize2, double markeralpha2, 
               std::string yaxislabel, double ymin=0.5, double ymax=1.5, 
               double topmarg=0.71, double botmarg=0.975, double lowx=1e-3, 
               std::string legendlabel1="ratio", std::string legendlabel2="ratio", 
               bool removexaxis=false, bool justratio=false, int ndiv=5) {   

    TPad *pad2;
    if (!justratio) pad2 = makeBottomPad(topmarg, botmarg);

    TLegend* l = new TLegend(0.50,0.300741,0.80,0.3585185,"");

    TH1D* hratio = (TH1D*) h1->Clone(ratio_name1.c_str());
    hratio->Divide(h3);
    hratio->SetMinimum(ymin);
    hratio->SetMaximum(ymax);

    TH1D* hratio2 = (TH1D*) h2->Clone(ratio_name2.c_str());
    hratio2->Divide(h3);
    
    //xtitlesize, xlabelsize, xoffset, ytitlesize, ylabelsize, yoffset                
    FormatHist(l, hratio, legendlabel1, markercolor1, markerstyle1, markeralpha1, 0.05, 0.04, 1.2, 0.035, 0.03, 1.05, markersize1);
    FormatHist(l, hratio2, legendlabel2, markercolor2, markerstyle2, markeralpha2, 0.05, 0.04, 1.2, 0.035, 0.03, 1.05, markersize2);
    hratio->GetYaxis()->SetNdivisions(5);
    hratio2->GetYaxis()->SetNdivisions(5);
    if (removexaxis) hratio->GetXaxis()->SetLabelSize(0);
    else hratio->GetXaxis()->SetTitleOffset(0.95);
    
    hratio->GetXaxis()->SetRangeUser(lowx,1);
    hratio2->GetXaxis()->SetRangeUser(lowx,1);

    hratio->Draw("same");
    hratio2->Draw("same");

    // hratio->GetYaxis()->SetTitleSize(0);
    hratio2->GetYaxis()->SetTitleSize(0);
    hratio->GetYaxis()->SetTitle(yaxislabel.c_str());
    hratio->GetYaxis()->SetTitleSize(0.06);
    hratio->GetYaxis()->SetTitleOffset(0.9);

    /*if (!justratio) {
        double ypos = ((1-topmarg) + botmarg) / 2;
        TLatex *t = new TLatex(0.075,ypos,yaxislabel.c_str());
        t->SetTextAlign(22);
        t->SetTextColor(kBlack);
        t->SetTextFont(43);
        t->SetTextSize(14);
        // t->SetTextAngle(45);
        t->SetNDC(kTRUE);
        t->Draw();

        // draw line at 1
        
        // drawHoriLine(1e-3, 1., 0.9, kGray+2, 6)->Draw();
        // drawHoriLine(1e-3, 1., 1.1, kGray+2, 6)->Draw();
    }
    */

    drawHoriLine(lowx, 1., 1., kGray+2, 1)->Draw();
    drawVertLine(0.4, ymin, ymax, kBlack, 2)->Draw();
    l->Draw("same");

    return pad2;

}


//===========================================================================//

void make_qg_plots_cgl() {

//    gROOT->SetBatch(); //prevents plots from showing up
    gStyle->SetOptStat(0);
    SetStyle();
    Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
    Double_t marker_size = 1.5;
    Double_t colors[16] = {kRed, kGreen+2, kBlue, kRed+1, kGreen+1, kBlue+1, kRed+2, kGreen+2, kBlue+2, kRed+3, kGreen+3, kBlue+3, kOrange+1, kViolet+1, kYellow+1, kCyan+1};

    //=====================================================================//
    // File containing quark vs gluon histograms
    // the c/b enhanced means those files will only have c/b jets
    // most are charged jets unless otherwise specified

    const TString basedir = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/";
    const TString root_filename = "/AnalysisResultsFinal.root";

    const TString infile_ptrl_charmdecaysON = basedir + "24691538" + root_filename; //--nocharmdecay 0 --giveptRL 1 **
    const TString infile_ptrl_charmdecaysOFF = basedir + "24692158" + root_filename; //--nocharmdecay 1 --giveptRL 1 **
    const TString infile_ptrl_c_enhanced_D0 = basedir + "24690700" + root_filename; //--replaceKP 1 --chinitscat 3 --giveptRL 1**
        
    //---------------------------------------------------------------------//
    // these histograms have D0 pt and z information:
    const TString infile_charmdecaysON = basedir + "25735401" + root_filename; //--nocharmdecay 0
    const TString infile_charmdecaysOFF = basedir + "25735384" + root_filename; //--nocharmdecay 1
    const TString infile_charmdecaysON_higherpthat = basedir + "25963614" + root_filename; //--nocharmdecay 0, higher pt-hat
    const TString infile_charmdecaysOFF_higherpthat = basedir + "25963640" + root_filename; //--nocharmdecay 1, higher pt-hat
    
    const TString infile_c_enhanced_D0 = basedir + "25536401" + root_filename; //--replaceKP 1 --chinitscat 3
    const TString infile_c_enhanced_D0_higherpthat = basedir + "25963105" + root_filename; //--replaceKP 1 --chinitscat 3, higher pt-hat

    const TString infile_c_enhanced_D0_fulljets = basedir + "26404382" + root_filename; //--replaceKP 1 --chinitscat 3 --fulljets 1

    const TString infile_c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons = basedir + "25536163" + root_filename; //--nocharmdecay 1 --chinitscat 1 //now --fulljets 2?
    const TString infile_b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons = basedir + "25536177" + root_filename; //--nobeautydecay 1 --chinitscat 5 //now --fulljets 2?
    
    const TString infile_c_enhanced_charmdecaysOFF_fulljets = basedir + "25620729" + root_filename; //--nocharmdecay 1 --chinitscat 1 --fulljets 1
    const TString infile_b_enhanced_beautydecaysOFF_fulljets = basedir + "25620738" + root_filename; //--nobeautydecay 1 --chinitscat 5 --fulljets 1
    const TString infile_c_enhanced_charmdecaysOFF_fulljets_higherpthat = basedir + "25963688" + root_filename; //--nocharmdecay 1 --chinitscat 1 --fulljets 1, higher pt-hat
    const TString infile_b_enhanced_beautydecaysOFF_fulljets_higherpthat = basedir + "25963694" + root_filename; //--nobeautydecay 1 --chinitscat 5 --fulljets 1, higher pt-hat
    
    const TString infile_charmdecaysON_fulljets = basedir + "25735421" + root_filename; //--nocharmdecay 0 --fulljets 1
    const TString infile_charmdecaysOFF_fulljets = basedir + "25735417" + root_filename; //--nocharmdecay 1 --fulljets 1
    const TString infile_charmdecaysON_fulljets_higherpthat = basedir + "25963982" + root_filename; //--nocharmdecay 0 --fulljets 1, higher pt-hat
    const TString infile_charmdecaysOFF_fulljets_higherpthat = basedir + "25963983" + root_filename; //--nocharmdecay 1 --fulljets 1, higher pt-hat
    
    //---------------------------------------------------------------------//
    // leading pt cut = 5 GeV
    const TString infile_charmdecaysON_leadpt5 = basedir + "26735768" + root_filename; //--nocharmdecay 0 --leadingptcut 5.0 
    const TString infile_c_enhanced_D0_leadpt5 = basedir + "26735755" + root_filename; //--replaceKP 1 --chinitscat 3 --leadingptcut 5.0 
    
    const TString infile_charmdecaysON_fulljets_leadpt5 = basedir + "26735877" + root_filename; //--nocharmdecay 0 --fulljets 1 --leadingptcut 5.0 
    const TString infile_c_enhanced_charmdecaysOFF_fulljets_leadpt5 = basedir + "26735879" + root_filename; //--nocharmdecay 1 --chinitscat 1 --fulljets 1 --leadingptcut 5.0 
    const TString infile_b_enhanced_beautydecaysOFF_fulljets_leadpt5 = basedir + "26735885" + root_filename; //--nobeautydecay 1 --chinitscat 5 --fulljets 1 --leadingptcut 5.0 

    //---------------------------------------------------------------------//
    // other
    const TString infile_c_enhanced_charmdecaysOFF = basedir + "24735947" + root_filename; //--nocharmdecay 1 --chinitscat 1
    const TString infile_c_enhanced_charmdecaysON = basedir + "24693346" + root_filename; //--nocharmdecay 0 --chinitscat 1 //currently not used

    const TString infile_c_enhanced_D0_rigidcone04 = basedir + "26834828" + root_filename; //--nocharmdecay 1 --chinitscat 1
    const TString infile_c_enhanced_D0_rigidcone1 = basedir + "26834904" + root_filename; //--nocharmdecay 0 --chinitscat 1 //currently not used

    //=====================================================================//

    // open the files
    TFile *f_charmdecaysON = new TFile(infile_charmdecaysON, "READ");
    TFile *f_charmdecaysOFF = new TFile(infile_charmdecaysOFF, "READ");
    TFile *f_charmdecaysON_fulljets = new TFile(infile_charmdecaysON_fulljets, "READ");
    TFile *f_charmdecaysOFF_fulljets = new TFile(infile_charmdecaysOFF_fulljets, "READ");
    TFile *f_ptrl_charmdecaysON = new TFile(infile_ptrl_charmdecaysON, "READ");
    TFile *f_ptrl_charmdecaysOFF = new TFile(infile_ptrl_charmdecaysOFF, "READ");
    TFile *f_ptrl_c_enhanced_D0 = new TFile(infile_ptrl_c_enhanced_D0, "READ");

    TFile *f_c_enhanced_D0 = new TFile(infile_c_enhanced_D0, "READ");
    TFile *f_c_enhanced_D0_fulljets = new TFile(infile_c_enhanced_D0_fulljets, "READ");
    TFile *f_c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons = new TFile(infile_c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons, "READ");
    TFile *f_b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons = new TFile(infile_b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons, "READ");
    TFile *f_c_enhanced_charmdecaysOFF_fulljets = new TFile(infile_c_enhanced_charmdecaysOFF_fulljets, "READ");
    TFile *f_b_enhanced_beautydecaysOFF_fulljets = new TFile(infile_b_enhanced_beautydecaysOFF_fulljets, "READ");
    TFile *f_c_enhanced_charmdecaysOFF = new TFile(infile_c_enhanced_charmdecaysOFF, "READ");
    TFile *f_c_enhanced_charmdecaysON = new TFile(infile_c_enhanced_charmdecaysON, "READ");

    TFile *f_charmdecaysON_leadpt5 = new TFile(infile_charmdecaysON_leadpt5, "READ");
    TFile *f_charmdecaysON_fulljets_leadpt5 = new TFile(infile_charmdecaysON_fulljets_leadpt5, "READ");
    TFile *f_c_enhanced_D0_leadpt5 = new TFile(infile_c_enhanced_D0_leadpt5, "READ");
    TFile *f_c_enhanced_charmdecaysOFF_fulljets_leadpt5 = new TFile(infile_c_enhanced_charmdecaysOFF_fulljets_leadpt5, "READ");
    TFile *f_b_enhanced_beautydecaysOFF_fulljets_leadpt5 = new TFile(infile_b_enhanced_beautydecaysOFF_fulljets_leadpt5, "READ");
    
    
    TFile *f_c_enhanced_D0_rigidcone04 = new TFile(infile_c_enhanced_D0_rigidcone04, "READ");
    TFile *f_c_enhanced_D0_rigidcone1 = new TFile(infile_c_enhanced_D0_rigidcone1, "READ");
    

    //CONTOL VARIABLES HERE
    // plot cases:
    //final:
    // 0: plot c, g, l, i --NOT IMPLEMENTED
        // g, l, i from charmdecaysON; c from c_enhanced_charmdecaysOFF
    // 1: plot D0, g, l, i, include ratio of D0/inclusive **
        // g, l, i from charmdecaysON; D0 from c_enhanced_D0 - thought about this, want ON bc we don't want to control any decays for the other samples
    // 2: plot D0, g, l, i but pt*RL -- TODO: already done in other file, port over here later
        // g, l, i from ptrl_charmdecaysON; D0 from ptrl_c_enhanced_D0
    // 3: plot l, c, b full jets **
        // l from charmdecaysON_fulljets; c from c_enhanced_charmdecaysOFF_fulljets, b from b_enhanced_beautydecaysOFF_fulljets
    // 4: RL edge effect
    // 100: plot the peak positions as a function of pt

    //exploring:
    // 10: plot c, g, l, i, b from full jets -- not neccesary...
        // g, l, i from charmdecaysON_fulljets; c from c_enhanced_charmdecaysOFF_fulljets; b from b_enhanced_beautydecaysOFF_fulljets
    // 11: plot l, c, b full jets and charged jets + neutrals hadrons, and charged jets
        // g, l, i from charmdecaysON_fulljets; c from c_enhanced_charmdecaysOFF_fulljets; b from b_enhanced_beautydecaysOFF_fulljets
        // don't have light for charged jets+neutral, but c from c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons; b from b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons
        // also add light charged jets from infile_charmdecaysON, c charged jets from infile_c_enhanced_charmdecaysOFF, don't have b charged jets
    // 12: plot_case 3 and also with charm decays OFF for full light jets
    // 13: plot_case 3 but with all having leading pt cut of 5 gev
    // 20: plot_case 1 and also with charm decays OFF from g, l, i
    // 21: plot_case 1 but with all having leading pt cut of 5 gev
    // 14: plot_case 1 and 21 combined
    int plot_case = 14;

    // std::vector<TFile*> files;
    std::string add_name;

    add_name = "_plotcase" + std::to_string(plot_case) + ".pdf";
    cout << "output name will be " << add_name << endl;

    // Output directory
    std::string outdir= "plots/final/cgl/"; //testing/";//"plots/test/";

    // Output file for binned results
    std::string outfile = outdir + "AnalysisResultsFinal_fromcgl.root";
    TFile* f_out = new TFile(outfile.c_str(), "RECREATE");

    // save pdf plots in more specific directories
    if (plot_case == 2) {
        outdir += "ptrl/";
    } else if (plot_case == 3 or plot_case == 11 or plot_case == 12 or plot_case == 13) {
        outdir += "beauty/";
    } else if (plot_case == 4) {
        outdir += "RLedge/";
    }

    // Jet r value
    std::string jetR_list[] = { "0.4" };
    for (std::string jetR : jetR_list) {

        // Names of histograms in the file (quark, charm, gluon)
        const std::string hc_name = "h_EEC_JetPt_charm_R" + jetR;
        const std::string hl_name = "h_EEC_JetPt_light_R" + jetR;
        const std::string hg_name = "h_EEC_JetPt_gluon_R" + jetR;
        const std::string hi_name = "h_EEC_JetPt_inclusive_R" + jetR;
        const std::string hb_name = "h_EEC_JetPt_beauty_R" + jetR;
        const std::string hc_jet_name = "h_JetPt_charm_R" + jetR + "_jetlevel";
        const std::string hl_jet_name = "h_JetPt_light_R" + jetR + "_jetlevel";
        const std::string hg_jet_name = "h_JetPt_gluon_R" + jetR + "_jetlevel";
        const std::string hi_jet_name = "h_JetPt_inclusive_R" + jetR + "_jetlevel";
        const std::string hb_jet_name = "h_JetPt_beauty_R" + jetR + "_jetlevel";

        //-------------------------------------------------//
        //TODO: IF THNSPARSE CLONES ARE BEFORE THE PT BIN LOOP, MOVE THEM INTO IT!!
        //-------------------------------------------------//
        // Make canvases that save across pt bins
        TCanvas* c_D0 = new TCanvas();
        ProcessCanvas(c_D0);
        c_D0->cd();
        gPad->SetLogx();
        // gPad->SetLogy();

        TCanvas* c_rat1 = new TCanvas();
        ProcessCanvas(c_rat1);
        c_rat1->cd();
        gPad->SetLogx();

        TCanvas *c_peaks = new TCanvas();
        ProcessCanvas(c_peaks);

        TCanvas *c_peaks_frac = new TCanvas();
        ProcessCanvas(c_peaks_frac);
        


        //-------------------------------------------------//

        //const int pt_bins[] = { 10, 20, 40, 60, 80, 100, 150 };
        const int pt_bins[] = { 7, 10, 15, 30, 50, 70, 100, 150, 200 }; //, 100, 150 }; //{ 10, 20, 40 };
        const int d0_pt_cuts[] = { 3, 5, 5, 5, 5, 5, 5, 5 }; //, 5, 5 };
        const int n_bins = 8; //7;
        double pt_x[n_bins] = { 10., 15., 30., 50., 70., 100., 150., 200.0 };
        
        double peaks_l_full_y[n_bins];
        double peaks_c_full_y[n_bins];
        double peaks_b_full_y[n_bins];
        double peaks_D0_full_y[n_bins];
        
        double peaks_l_ch_y[n_bins];
        double peaks_g_ch_y[n_bins];
        double peaks_incl_ch_y[n_bins];
        double peaks_D0_ch_y[n_bins];

        double peakserr_l_full_y[n_bins];
        double peakserr_c_full_y[n_bins];
        double peakserr_b_full_y[n_bins];
        double peakserr_D0_full_y[n_bins];
        
        double peakserr_l_ch_y[n_bins];
        double peakserr_g_ch_y[n_bins];
        double peakserr_incl_ch_y[n_bins];
        double peakserr_D0_ch_y[n_bins];

        //Format color and style
        int markercolor_c = kRed; //charm
        int markerstyle_c = kFullCircle;
        int markercolor_g = kViolet+2; //gluon
        int markerstyle_g = 33; //diamond
        int markercolor_l = kGreen+2; //light
        int markerstyle_l = 21; //square
        int markercolor_i = kMagenta-6; //inclusive
        int markerstyle_i = 29; //star
        int markercolor_D0 = kOrange+7; //kGreen-5; //D0
        int markerstyle_D0 = kFullCircle;
        int markercolor_b = kAzure; //beauty
        int markerstyle_b = 34;

        int markerstyle_c_ch = kOpenCircle; //charged jets when being compared to full
        int markerfill_c_chn = 3944; //charged jets with neutral hadrons

        int markerstyle_l_ch = kOpenSquare; //charged jets when being compared to full
        int markerfill_b_chn = 3944; //charged jets with neutral hadrons

        
        for (int i = 0; i < n_bins; i++) {
            cout << "in pt bin" << i << endl;
            int pt_min = pt_bins[i];
            int pt_max = pt_bins[i+1];

            // define pt related variables
            TString ptbin = TString::Format("%d #leq #it{p}_{T}^{ch. jet} < %d GeV/#it{c}, #font[122]{|}#it{#eta}_{jet}#font[122]{|} #leq 0.5", pt_min, pt_max);
            TString ptD = TString::Format("%d #leq #it{p}_{T}^{D^{0}} < %d GeV/#it{c}, #font[122]{|}#it{y}_{D^{0}}#font[122]{|} #leq 0.8", d0_pt_cuts[i], pt_max);
            if (plot_case == 3 or plot_case == 11 or plot_case == 12) {
                ptbin = TString::Format("%d #leq #it{p}_{T}^{full jet} < %d GeV/#it{c}, #font[122]{|}#it{#eta}_{jet}#font[122]{|} #leq 0.5", pt_min, pt_max);
            } else if (plot_case == 13 or plot_case == 21) {
                ptD = TString::Format("%d #leq #it{p}_{T}^{D^{0}} < %d GeV/#it{c}, #font[122]{|}#it{y}_{D^{0}}#font[122]{|} #leq 0.8", 5, pt_max);
            }

            // make a canvas for each pt range
            TCanvas* c = new TCanvas();
            ProcessCanvas(c);
            c->cd();
            gPad->SetLogx();
            // if (plot_case == 1 or plot_case == 20) {
            //     c->SetCanvasSize(900, 900);
            // }
            // gPad->SetLogy();

            // if (plot_case == 4) {
            //     TCanvas* c_ratio_RLedge = new TCanvas();
            //     ProcessCanvas(c_ratio_RLedge);
            //     c_ratio_RLedge->cd();
            //     gPad->SetLogx();
            // }

            TLegend* l; // = new TLegend(0.17, 0.65, 0.5, 0.85);
            TLegend* l2 = new TLegend(0.1797168,0.400741,0.4562155,0.8885185,""); //dummy legend

            double maxy = 0;


            // Open histograms


            // TODO: figure out plot_case 13 and 21!!
            l = new TLegend(0.1957168,0.600741,0.462155,0.9485185,"");
            // else l = new TLegend(0.1797168,0.400741,0.4562155,0.8885185,""); //(0.17, 0.4, 0.5, 0.53);
            l->SetTextSize(0.045);
            l->AddEntry("NULL","PYTHIA 8 Monash 2013","h");
            l->AddEntry("NULL","pp, #sqrt{#it{s}} = 13 TeV","h");
            if (plot_case == 3 or plot_case == 11 or plot_case == 12) l->AddEntry("NULL","anti-#it{k}_{T}, #it{R} = 0.4","h");
            else {
                l->AddEntry("NULL","D^{0} #rightarrow K^{#minus} #pi^{+} and charge conj.","h");
                l->AddEntry("NULL","in charged jets, anti-#it{k}_{T}, #it{R} = 0.4","h");
            }
            l->AddEntry("NULL",ptbin,"h");
            if (!(plot_case == 3 or plot_case == 11 or plot_case == 12)) l->AddEntry("NULL",ptD,"h");
            // if (plot_case == 1 or plot_case == 20) 
            l->SetTextSize(0.028);
            // else l->SetTextSize(0.037);
            l->SetBorderSize(0);
            l->SetFillStyle(0);
            // l->Draw("same");


            std::string pt_name = "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);


            //format: getObsHist(TFile *filename, std::string h_name, std::string h_jet_name, int pt_min, int pt_max, int d0_pt_cut, std::string newhistname, bool d0cuts=false, int obsaxis=3)
            TH1D *h_charmdecaysOFF_g = getObsHist(f_charmdecaysOFF, hg_name, hg_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_charmdecaysOFF_g" + pt_name, false, 4);
            TH1D *h_charmdecaysOFF_l = getObsHist(f_charmdecaysOFF, hl_name, hl_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_charmdecaysOFF_l" + pt_name, false, 4);
            TH1D *h_charmdecaysOFF_i = getObsHist(f_charmdecaysOFF, hi_name, hi_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_charmdecaysOFF_i" + pt_name, false, 4);
            TH1D *h_c_enhanced_charmdecaysOFF_c = getObsHist(f_c_enhanced_charmdecaysOFF, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_c_enhanced_charmdecaysOFF_c" + pt_name, false, 3);
            
            TH1D *h_charmdecaysON_g = getObsHist(f_charmdecaysON, hg_name, hg_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_charmdecaysON_g" + pt_name, false, 4);
            TH1D *h_charmdecaysON_l = getObsHist(f_charmdecaysON, hl_name, hl_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_charmdecaysON_l" + pt_name, false, 4);
            TH1D *h_charmdecaysON_i = getObsHist(f_charmdecaysON, hi_name, hi_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_charmdecaysON_i" + pt_name, false, 4);
            
            TH1D *h_ptrl_charmdecaysOFF_g = getObsHist(f_ptrl_charmdecaysOFF, hg_name, hg_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_ptrl_charmdecaysOFF_g" + pt_name, false, 3, true);
            TH1D *h_ptrl_charmdecaysOFF_l = getObsHist(f_ptrl_charmdecaysOFF, hl_name, hl_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_ptrl_charmdecaysOFF_l" + pt_name, false, 3, true);
            TH1D *h_ptrl_charmdecaysOFF_i = getObsHist(f_ptrl_charmdecaysOFF, hi_name, hi_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_ptrl_charmdecaysOFF_i" + pt_name, false, 3, true);
            TH1D *h_ptrl_charmdecaysOFF_c = getObsHist(f_ptrl_charmdecaysOFF, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_ptrl_charmdecaysOFF_c" + pt_name, false, 3, true);
            
            TH1D *h_ptrl_charmdecaysON_g = getObsHist(f_ptrl_charmdecaysON, hg_name, hg_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_ptrl_charmdecaysON_g" + pt_name, false, 3, true);
            TH1D *h_ptrl_charmdecaysON_l = getObsHist(f_ptrl_charmdecaysON, hl_name, hl_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_ptrl_charmdecaysON_l" + pt_name, false, 3, true);
            TH1D *h_ptrl_charmdecaysON_i = getObsHist(f_ptrl_charmdecaysON, hi_name, hi_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_ptrl_charmdecaysON_i" + pt_name, false, 3, true);
            
            TH1D *h_c_enhanced_D0 = getObsHist(f_c_enhanced_D0, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_c_enhanced_D0" + pt_name, true, 4);

            TH1D *h_c_enhanced_D0_fulljets = getObsHist(f_c_enhanced_D0_fulljets, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_c_enhanced_D0_fulljets" + pt_name, true, 4);

            TH1D *h_ptrl_c_enhanced_D0 = getObsHist(f_ptrl_c_enhanced_D0, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_ptrl_c_enhanced_D0" + pt_name, true, 3, true);

            TH1D *h_charmdecaysOFF_fulljets_g = getObsHist(f_charmdecaysOFF_fulljets, hg_name, hg_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_charmdecaysOFF_fulljets_g" + pt_name, false, 4);
            TH1D *h_charmdecaysOFF_fulljets_l = getObsHist(f_charmdecaysOFF_fulljets, hl_name, hl_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_charmdecaysOFF_fulljets_l" + pt_name, false, 4);
            TH1D *h_c_enhanced_charmdecaysOFF_fulljets_c = getObsHist(f_c_enhanced_charmdecaysOFF_fulljets, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_c_enhanced_charmdecaysOFF_fulljets_c" + pt_name, false, 4);
            TH1D *h_b_enhanced_beautydecaysOFF_fulljets_b = getObsHist(f_b_enhanced_beautydecaysOFF_fulljets, hb_name, hb_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_b_enhanced_beautydecaysOFF_fulljets_b" + pt_name, false, 4);

            TH1D *h_charmdecaysON_fulljets_g = getObsHist(f_charmdecaysON_fulljets, hg_name, hg_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_charmdecaysON_fulljets_g" + pt_name, false, 4);
            TH1D *h_charmdecaysON_fulljets_l = getObsHist(f_charmdecaysON_fulljets, hl_name, hl_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_charmdecaysON_fulljets_l" + pt_name, false, 4);
            
            TH1D *h_c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons_c = getObsHist(f_c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons_c" + pt_name, false, 4);
            TH1D *h_b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons_b = getObsHist(f_b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons, hb_name, hb_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons_b" + pt_name, false, 4);

            TH1D *h_charmdecaysON_leadpt5_g = getObsHist(f_charmdecaysON_leadpt5, hg_name, hg_jet_name, pt_min, pt_max, 5, "h_charmdecaysON_leadpt5_g" + pt_name, false, 4);
            TH1D *h_charmdecaysON_leadpt5_l = getObsHist(f_charmdecaysON_leadpt5, hl_name, hl_jet_name, pt_min, pt_max, 5, "h_charmdecaysON_leadpt5_l" + pt_name, false, 4);
            TH1D *h_charmdecaysON_leadpt5_i = getObsHist(f_charmdecaysON_leadpt5, hi_name, hi_jet_name, pt_min, pt_max, 5, "h_charmdecaysON_leadpt5_i" + pt_name, false, 4);
            TH1D *h_charmdecaysON_fulljets_leadpt5_l = getObsHist(f_charmdecaysON_fulljets_leadpt5, hl_name, hl_jet_name, pt_min, pt_max, 5, "h_charmdecaysON_fulljets_leadpt5_l" + pt_name, false, 4);
            TH1D *h_c_enhanced_D0_leadpt5_c = getObsHist(f_c_enhanced_D0_leadpt5, hc_name, hc_jet_name, pt_min, pt_max, 5, "h_c_enhanced_D0_leadpt5_c" + pt_name, false, 4);
            TH1D *h_c_enhanced_charmdecaysOFF_fulljets_leadpt5_c = getObsHist(f_c_enhanced_charmdecaysOFF_fulljets_leadpt5, hc_name, hc_jet_name, pt_min, pt_max, 5, "h_c_enhanced_charmdecaysOFF_fulljets_leadpt5_c" + pt_name, false, 4);
            TH1D *h_b_enhanced_beautydecaysOFF_fulljets_leadpt5_b = getObsHist(f_b_enhanced_beautydecaysOFF_fulljets_leadpt5, hb_name, hb_jet_name, pt_min, pt_max, 5, "h_b_enhanced_beautydecaysOFF_fulljets_leadpt5_b" + pt_name, false, 4);
    
            TH1D *h_c_enhanced_D0_rigidcone04 = getObsHist(f_c_enhanced_D0_rigidcone04, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_c_enhanced_D0_rigidcone04" + pt_name, false, 4);
            TH1D *h_c_enhanced_D0_rigidcone1 = getObsHist(f_c_enhanced_D0_rigidcone1, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_c_enhanced_D0_rigidcone1" + pt_name, false, 4);

            //make dummy histogram
            TH1D *h_dummy = getObsHist(f_charmdecaysOFF, hg_name, hg_jet_name, pt_min, pt_max, d0_pt_cuts[i], "hdummy", false, 4);

            


            // // Format histograms for plotting (this order needed to keep legend in order and graphs lookin good)

            // if (charmdecays) {
            //     l->AddEntry("NULL","          D* decays on","h");
            // } else {
            //     l->AddEntry("NULL","          D* decays off","h");
            // }




            // FormatHist(l, h_ptrl_charmdecaysOFF_c, "charm-init jets", markercolor_c, markerstyle_c);

            // FormatHist(l, h_charmdecaysOFF_fulljets_g, "gluon-init full jets", markercolor_g, markerstyle_g);

            TPad *pad1;
            TPad *pad2;

            // h_dummy->GetXaxis()->SetRangeUser(3e-4, 1.);
            // h_dummy->SetMarkerColor(0);

            if (plot_case == 0) {

            } else if (plot_case == 1 or plot_case == 20) { //D0, g, l, i in charged jets

                // pad1 = makeTopPad(h_dummy, h_charmdecaysON_l->GetMaximum()*1.2, 0.41);
                pad1 = makeTopPad(0.41);
               
                double in_yaxis_min=0;
                double in_yaxis_max=0;
                double in_xaxis_min=0;
                double in_xaxis_max=0;
                
                double l_yaxis_min=0;
                double l_yaxis_max=0;
                double l_xaxis_min=0;
                double l_xaxis_max=0;
                double lowx=0;
                
                if (pt_min==7){in_yaxis_min=0.;in_yaxis_max=0.6; l_yaxis_min=0.;l_yaxis_max=0.7;lowx=1e-2;}
                if (pt_min==10){in_yaxis_min=0.;in_yaxis_max=0.6; l_yaxis_min=0.;l_yaxis_max=0.8;lowx=0.6*1e-2;}
                if (pt_min==15){in_yaxis_min=0.0;in_yaxis_max=0.7; l_yaxis_min=0.;l_yaxis_max=1;lowx=0.6*1e-2;}
                if (pt_min==30){in_yaxis_min=0.3;in_yaxis_max=1.1; l_yaxis_min=0.;l_yaxis_max=1;lowx=0.5*1e-2;}
                if (pt_min==50){in_yaxis_min=0.0;in_yaxis_max=1.1; l_yaxis_min=0.3;l_yaxis_max=1;lowx=1e-3;}
                if (pt_min==70){in_yaxis_min=0.0;in_yaxis_max=1.1; l_yaxis_min=0.4;l_yaxis_max=1.1;lowx=1e-3;}
                if (pt_min==100){in_yaxis_min=0.0;in_yaxis_max=1.1; l_yaxis_min=0.4;l_yaxis_max=1.1;lowx=0.5*1e-3;}
                if (pt_min==150){in_yaxis_min=0.0;in_yaxis_max=0; l_yaxis_min=0.;l_yaxis_max=1.1;}
   
                
                h_charmdecaysON_l->SetMaximum(h_charmdecaysON_l->GetMaximum()*1.2);
                
                drawWithLineAtMax(h_charmdecaysON_l, l, "light-init jets", markercolor_l, markerstyle_l, lowx, false, true, 0.8);
                drawWithLineAtMax(h_charmdecaysON_g, l, "gluon-init jets", markercolor_g, markerstyle_g, lowx, false, true, 0.8);
                drawWithLineAtMax(h_charmdecaysON_i, l, "inclusive jets", markercolor_i, markerstyle_i, lowx, false, true, 0.8);
                drawWithLineAtMax(h_c_enhanced_D0, l, "D^{0}-tagged, c-init jets", markercolor_D0, markerstyle_D0, lowx, false, true, 0.8);
                
                //make ratio plots!
                std::string yaxislabel = "#frac{D^{0}-tagged jets}{inclusive jets}";
                plotRatio(h_c_enhanced_D0, h_charmdecaysON_i, "hratio_Djet_incl" + pt_name + "_R" + jetR + add_name, l2, kBlack, markers[2], yaxislabel, in_yaxis_min, in_yaxis_max, 0.61, 0.25, lowx, "", true);
                drawHoriLine(lowx, 1., 0.5, kGray+2, 6)->Draw();
                
                yaxislabel = "#frac{D^{0}-tagged jets}{light jets}";
                plotRatio(h_c_enhanced_D0, h_charmdecaysON_l, "hratio_Djet_light" + pt_name + "_R" + jetR + add_name, l2, kBlack, markers[2], yaxislabel,l_yaxis_min, l_yaxis_max, 0.75, 0.11, lowx);
                drawHoriLine(lowx, 1., 0.5, kGray+2, 6)->Draw();


            }
            if (plot_case == 21) { //D0, g, l, i in charged jets, 5 GeV leading pt cut

                // pad1 = makeTopPad(h_dummy, h_charmdecaysON_l->GetMaximum()*1.2, 0.41);
                if (plot_case == 21) pad1 = makeTopPad(0.41);
               
                double in_yaxis_min=0;
                double in_yaxis_max=0;
                double in_xaxis_min=0;
                double in_xaxis_max=0;
                
                double l_yaxis_min=0;
                double l_yaxis_max=0;
                double l_xaxis_min=0;
                double l_xaxis_max=0;
                double lowx=0;
            
                if (pt_min==7){in_yaxis_min=0.;in_yaxis_max=0.6; l_yaxis_min=0.;l_yaxis_max=0.7;lowx=1e-2;}
                if (pt_min==10){in_yaxis_min=0.;in_yaxis_max=0.6; l_yaxis_min=0.;l_yaxis_max=0.8;lowx=0.6*1e-2;}
                if (pt_min==15){in_yaxis_min=0.0;in_yaxis_max=0.7; l_yaxis_min=0.;l_yaxis_max=1;lowx=0.6*1e-2;}
                if (pt_min==30){in_yaxis_min=0.3;in_yaxis_max=1.1; l_yaxis_min=0.;l_yaxis_max=1;lowx=0.5*1e-2;}
                if (pt_min==50){in_yaxis_min=0.0;in_yaxis_max=1.1; l_yaxis_min=0.3;l_yaxis_max=1;lowx=1e-3;}
                if (pt_min==70){in_yaxis_min=0.0;in_yaxis_max=1.1; l_yaxis_min=0.4;l_yaxis_max=1.1;lowx=1e-3;}
                if (pt_min==100){in_yaxis_min=0.0;in_yaxis_max=1.1; l_yaxis_min=0.4;l_yaxis_max=1.1;lowx=0.5*1e-3;}
                if (pt_min==150){in_yaxis_min=0.0;in_yaxis_max=0; l_yaxis_min=0.;l_yaxis_max=1.1;}
   
                
                h_charmdecaysON_leadpt5_l->SetMaximum(h_charmdecaysON_leadpt5_l->GetMaximum()*1.2);
                
                // drawWithLineAtMax(h_charmdecaysON_l, l, "light-init jets", markercolor_l, markerstyle_l, lowx, false, true, 0.8);
                // drawWithLineAtMax(h_charmdecaysON_g, l, "gluon-init jets", markercolor_g, markerstyle_g, lowx, false, true, 0.8);
                // drawWithLineAtMax(h_charmdecaysON_i, l, "inclusive jets", markercolor_i, markerstyle_i, lowx, false, true, 0.8);
                // drawWithLineAtMax(h_c_enhanced_D0, l, "D^{0}-tagged, c-init jets", markercolor_D0, markerstyle_D0, lowx, false, true, 0.8);
                drawWithLineAtMax(h_charmdecaysON_leadpt5_l, l, "light-init jets", markercolor_l, kOpenSquare, lowx, false, true, 0.8);
                drawWithLineAtMax(h_charmdecaysON_leadpt5_g, l, "gluon-init jets", markercolor_g, kOpenDiamond, lowx, false, true, 0.8);
                drawWithLineAtMax(h_charmdecaysON_leadpt5_i, l, "inclusive jets", markercolor_i, kOpenStar, lowx, false, true, 0.8);
                drawWithLineAtMax(h_c_enhanced_D0_leadpt5_c, l, "D^{0}-tagged, c-init jets", markercolor_D0, kOpenCircle, lowx, false, true, 0.8);
                
                //make ratio plots!
                // std::string yaxislabel = "#frac{D^{0}-tagged jets}{inclusive jets}";
                // plotRatio(h_c_enhanced_D0_leadpt5_c, h_charmdecaysON_leadpt5_i, "hratio_Djet_incl" + pt_name + "_R" + jetR + add_name, l2, kBlack, markers[2], yaxislabel, in_yaxis_min, in_yaxis_max, 0.61, 0.25, lowx, "", true);
                // drawHoriLine(lowx, 1., 0.5, kGray+2, 6)->Draw();
                
                // yaxislabel = "#frac{D^{0}-tagged jets}{light jets}";
                // plotRatio(h_c_enhanced_D0_leadpt5_c, h_charmdecaysON_leadpt5_l, "hratio_Djet_light" + pt_name + "_R" + jetR + add_name, l2, kBlack, markers[2], yaxislabel,l_yaxis_min, l_yaxis_max, 0.75, 0.11, lowx);
                // drawHoriLine(lowx, 1., 0.5, kGray+2, 6)->Draw();


            } if (plot_case == 14) { //D0, g, l, i in charged jets

                // pad1 = makeTopPad(h_dummy, h_charmdecaysON_l->GetMaximum()*1.2, 0.41);
                pad1 = makeTopPad(0.41);
               
                double in_yaxis_min=0;
                double in_yaxis_max=0;
                double in_xaxis_min=0;
                double in_xaxis_max=0;
                
                double l_yaxis_min=0;
                double l_yaxis_max=0;
                double l_xaxis_min=0;
                double l_xaxis_max=0;
                double lowx=0;
                
                if (pt_min==7){in_yaxis_min=0.;in_yaxis_max=0.6; l_yaxis_min=0.;l_yaxis_max=0.7;lowx=1e-2;}
                if (pt_min==10){in_yaxis_min=0.;in_yaxis_max=0.7; l_yaxis_min=0.;l_yaxis_max=0.9;lowx=0.6*1e-2;}
                if (pt_min==15){in_yaxis_min=0.0;in_yaxis_max=0.79; l_yaxis_min=0.;l_yaxis_max=1.1;lowx=0.6*1e-2;}
                if (pt_min==30){in_yaxis_min=0.3;in_yaxis_max=1.1; l_yaxis_min=0.;l_yaxis_max=1.1;lowx=0.5*1e-2;}
                if (pt_min==50){in_yaxis_min=0.0;in_yaxis_max=1.1; l_yaxis_min=0.3;l_yaxis_max=1.1;lowx=1e-3;}
                if (pt_min==70){in_yaxis_min=0.0;in_yaxis_max=1.1; l_yaxis_min=0.4;l_yaxis_max=1.1;lowx=1e-3;}
                if (pt_min==100){in_yaxis_min=0.0;in_yaxis_max=1.1; l_yaxis_min=0.4;l_yaxis_max=1.1;lowx=0.5*1e-3;}
                if (pt_min==150){in_yaxis_min=0.0;in_yaxis_max=0; l_yaxis_min=0.;l_yaxis_max=1.1;}
   
                
                h_charmdecaysON_l->SetMaximum(h_charmdecaysON_l->GetMaximum()*1.2);
                
                drawWithLineAtMax(h_charmdecaysON_l, l, "light-init jets", markercolor_l, markerstyle_l, lowx, false, true, 0.6);
                drawWithLineAtMax(h_charmdecaysON_leadpt5_l, l2, "light-init jets", markercolor_l, kOpenSquare, lowx, false, true, 1.0);
                
                drawWithLineAtMax(h_charmdecaysON_g, l, "gluon-init jets", markercolor_g, markerstyle_g, lowx, false, true, 0.6);
                drawWithLineAtMax(h_charmdecaysON_leadpt5_g, l2, "gluon-init jets", markercolor_g, kOpenDiamond, lowx, false, true, 1.0);
                
                drawWithLineAtMax(h_charmdecaysON_i, l, "inclusive jets", markercolor_i, markerstyle_i, lowx, false, true, 0.6);
                drawWithLineAtMax(h_charmdecaysON_leadpt5_i, l2, "inclusive jets", markercolor_i, kOpenStar, lowx, false, true, 1.0);
                
                drawWithLineAtMax(h_c_enhanced_D0, l, "D^{0}-tagged, c-init jets", markercolor_D0, markerstyle_D0, lowx, false, true, 0.6);
                drawWithLineAtMax(h_c_enhanced_D0_leadpt5_c, l2, "D^{0}-tagged, c-init jets", markercolor_D0, kOpenCircle, lowx, false, true, 1.0);
                
                //purely for the legend
                // l->AddEntry("NULL", "", "h");
                TLegend* l3 = new TLegend(0.5797168,0.890741,0.875,0.9485185,""); //dummy legend
                l3->SetTextSize(0.028);
                l3->SetBorderSize(0);
                l3->SetFillStyle(0);
                Double_t xedges[2] = {lowx, 1.0};
                TH1D* h_dummy2 = new TH1D("h_dummy2", "h_dummy2", 1, xedges);
                TH1D* h_dummy2_clone = (TH1D*) h_dummy2->Clone("");
                drawNoLine(h_dummy2, l3, "no leading track p_{T} cut", kBlack, kFullDiamond, lowx, 1.5, true, 1.0);
                drawNoLine(h_dummy2_clone, l3, "leading track p_{T} cut = 5 GeV/c", kBlack, kOpenDiamond, lowx, 1.5, true, 1.0);
                l3->Draw("same");

                //make ratio plots!
                std::string yaxislabel = "#frac{D^{0}-tagged jets}{inclusive jets}";
                pad2 = plotRatio(h_c_enhanced_D0, h_charmdecaysON_i, "hratio_Djet_incl" + pt_name + "_R" + jetR + add_name, l2, kBlack, markers[2], yaxislabel, in_yaxis_min, in_yaxis_max, 0.61, 0.25, lowx, "", true);
                
                pad2->cd();
                TH1D* hratio2 = (TH1D*) h_c_enhanced_D0_leadpt5_c->Clone("");
                hratio2->Divide(h_charmdecaysON_leadpt5_i);
                                FormatHist(l2, hratio2, "", kBlack, kOpenDiamond); //, 1.0, 0.05, 0.04, 1.2, 0.035, 0.03, 1.5);
                hratio2->Draw("same");
                drawHoriLine(lowx, 1., 0.5, kGray+2, 6)->Draw();
                
                yaxislabel = "#frac{D^{0}-tagged jets}{light jets}";
                TPad * pad3 = plotRatio(h_c_enhanced_D0, h_charmdecaysON_l, "hratio_Djet_light" + pt_name + "_R" + jetR + add_name, l2, kBlack, markers[2], yaxislabel,l_yaxis_min, l_yaxis_max, 0.75, 0.11, lowx);
                pad3->cd();
                TH1D* hratio3 = (TH1D*) h_c_enhanced_D0_leadpt5_c->Clone("");
                hratio3->Divide(h_charmdecaysON_leadpt5_l);
                                FormatHist(l2, hratio3, "", kBlack, kOpenDiamond); //, 1.0, 0.05, 0.04, 1.2, 0.035, 0.03, 1.5);
                hratio3->Draw("same");
                drawHoriLine(lowx, 1., 0.5, kGray+2, 6)->Draw();

    

            } else if (plot_case == 2) { //D0, g, l, i pT*RL
                double lowx = 1e-4;
                h_ptrl_charmdecaysOFF_l->SetMaximum(h_ptrl_charmdecaysOFF_l->GetMaximum()*1.2);                
                
                drawWithLineAtMax(h_ptrl_charmdecaysON_l, l, "light-init jets", markercolor_l, markerstyle_l, lowx, true, false, 0.8);
                drawWithLineAtMax(h_ptrl_charmdecaysON_g, l, "gluon-init jets", markercolor_g, markerstyle_g, lowx, true, false, 0.8);
                drawWithLineAtMax(h_ptrl_charmdecaysON_i, l, "inclusive jets", markercolor_i, markerstyle_i, lowx, true, false, 0.8);
                // drawWithLineAtMax(h_ptrl_c_enhanced_charmdecaysOFF_chargedjets_c, l, "charm-init jets", markercolor_c, markerstyle_c, lowx, true, false, 0.8);
                drawWithLineAtMax(h_ptrl_c_enhanced_D0, l, "D^{0}-tagged, c-init jets", markercolor_D0, markerstyle_D0, lowx, true, false, 0.8);

                // TH1D* h_ptrl_c_enhanced_D0
            } else if (plot_case == 3 or plot_case == 11 or plot_case == 12) { //l, c, b
                pad1 = makeTopPad(0.41); //0.31);
                
                h_charmdecaysON_fulljets_l->SetMaximum(h_charmdecaysON_fulljets_l->GetMaximum()*1.2);
                

               
                // pad2 = plotRatio(h_c_enhanced_charmdecaysOFF_fulljets_c, h_charmdecaysON_fulljets_l, "hratio_charm_light" + pt_name + "_R" + jetR + add_name, l2, markercolor_c, markers[2], yaxislabel, 0., 3., 0.71, 0.11, true);
                // pad2->cd();
                // plotRatio2(h_b_enhanced_beautydecaysOFF_fulljets_b, h_charmdecaysON_fulljets_l, "hratio_beauty_light" + pt_name + "_R" + jetR + add_name, l2, markercolor_b, markers[2], "", 0., 1.48, 0.71, 0.11, true, true)->Draw();
                // drawHoriLine(1e-3, 1., 0.5, kGray+2, 6)->Draw();
                
                // yaxislabel = "#frac{D^{0}-tagged jets}{light jets}";
                // plotRatio(h_c_enhanced_D0, h_charmdecaysON_l, "hratio_Djet_light" + pt_name + "_R" + jetR + add_name, l2, markers[2], yaxislabel, 0., 1.48, 0.75, 0.11);
                // drawHoriLine(1e-3, 1., 0.5, kGray+2, 6)->Draw();
                
                double c_yaxis_min=0;
                double c_yaxis_max=0;
                double c_xaxis_min=0;
                double c_xaxis_max=0;
                
                double b_yaxis_min=0;
                double b_yaxis_max=0;
                double b_xaxis_min=0;
                double b_xaxis_max=0;
                double lowx=0;
                
                if (pt_min==7){c_yaxis_min=0.3;c_yaxis_max=0.82; b_yaxis_min=0.;b_yaxis_max=0.5;lowx=1e-2;}
                if (pt_min==10){c_yaxis_min=0.35;c_yaxis_max=0.82; b_yaxis_min=0.;b_yaxis_max=0.5;lowx=0.6*1e-2;}
                if (pt_min==15){c_yaxis_min=0.35;c_yaxis_max=1; b_yaxis_min=0.;b_yaxis_max=0.8;lowx=0.6*1e-2;}
                if (pt_min==30){c_yaxis_min=0.35;c_yaxis_max=1.1; b_yaxis_min=0.;b_yaxis_max=0.9;lowx=0.3*1e-2;}
                if (pt_min==50){c_yaxis_min=0.35;c_yaxis_max=1.1; b_yaxis_min=0.;b_yaxis_max=1;lowx=1e-3;}
                if (pt_min==70){c_yaxis_min=0.35;c_yaxis_max=1.1; b_yaxis_min=0.;b_yaxis_max=1.1;lowx=1e-3;}
                if (pt_min==100){c_yaxis_min=0.35;c_yaxis_max=1.1; b_yaxis_min=0.;b_yaxis_max=1.1;lowx=0.5*1e-3;}
                if (pt_min==150){c_yaxis_min=0.35;c_yaxis_max=1.1; b_yaxis_min=0.;b_yaxis_max=1.1;}
   
                
                drawWithLineAtMax(h_charmdecaysON_fulljets_l, l, "light-init full jets", markercolor_l, markerstyle_l, lowx, false, true, 0.8);
                drawWithLineAtMax(h_c_enhanced_charmdecaysOFF_fulljets_c, l, "charm-init full jets", markercolor_c, markerstyle_c, lowx, false, true, 0.8);
                drawWithLineAtMax(h_b_enhanced_beautydecaysOFF_fulljets_b, l, "beauty-init full jets", markercolor_b, markerstyle_b, lowx, false, true, 0.8);
                
                //make ratio plots!
                std::string yaxislabel = "#frac{charm jets}{light jets}";//"#frac{heavy quark jets}{light quark jets}";
                
               
                plotRatio(h_c_enhanced_charmdecaysOFF_fulljets_c, h_charmdecaysON_fulljets_l, "hratio_charm_light" + pt_name + "_R" + jetR + add_name, l2, markercolor_c, markers[2], yaxislabel, c_yaxis_min, c_yaxis_max, 0.61, 0.25, lowx, "", true);
                drawHoriLine(lowx, 1., 0.5, kGray+2, 6)->Draw();
                
                yaxislabel = "#frac{beauty jets}{light jets}";
                plotRatio(h_b_enhanced_beautydecaysOFF_fulljets_b, h_charmdecaysON_fulljets_l, "hratio_beauty_light" + pt_name + "_R" + jetR + add_name, l2, markercolor_b, markers[2], yaxislabel, b_yaxis_min, b_yaxis_max, 0.75, 0.11, lowx);
                drawHoriLine(lowx, 1., 0.5, kGray+2, 6)->Draw();


            } else if (plot_case == 13) { //l, c, b
                pad1 = makeTopPad(0.41); //0.31);
                
                h_charmdecaysON_fulljets_leadpt5_l->SetMaximum(h_charmdecaysON_fulljets_leadpt5_l->GetMaximum()*1.2);
                
                double c_yaxis_min=0;
                double c_yaxis_max=0;
                double c_xaxis_min=0;
                double c_xaxis_max=0;
                
                double b_yaxis_min=0;
                double b_yaxis_max=0;
                double b_xaxis_min=0;
                double b_xaxis_max=0;
                double lowx=0;
                
                if (pt_min==7){c_yaxis_min=0.3;c_yaxis_max=0.82; b_yaxis_min=0.;b_yaxis_max=0.5;lowx=1e-2;}
                if (pt_min==10){c_yaxis_min=0.35;c_yaxis_max=0.82; b_yaxis_min=0.;b_yaxis_max=0.5;lowx=0.6*1e-2;}
                if (pt_min==15){c_yaxis_min=0.35;c_yaxis_max=1; b_yaxis_min=0.;b_yaxis_max=0.8;lowx=0.6*1e-2;}
                if (pt_min==30){c_yaxis_min=0.35;c_yaxis_max=1.1; b_yaxis_min=0.;b_yaxis_max=0.9;lowx=0.3*1e-2;}
                if (pt_min==50){c_yaxis_min=0.35;c_yaxis_max=1.1; b_yaxis_min=0.;b_yaxis_max=1;lowx=1e-3;}
                if (pt_min==70){c_yaxis_min=0.35;c_yaxis_max=1.1; b_yaxis_min=0.;b_yaxis_max=1.1;lowx=1e-3;}
                if (pt_min==100){c_yaxis_min=0.35;c_yaxis_max=1.1; b_yaxis_min=0.;b_yaxis_max=1.1;lowx=0.5*1e-3;}
                if (pt_min==150){c_yaxis_min=0.35;c_yaxis_max=1.1; b_yaxis_min=0.;b_yaxis_max=1.1;}
   
                
                // drawWithLineAtMax(h_charmdecaysON_fulljets_l, l, "light-init full jets", markercolor_l, markerstyle_l, lowx, false, true, 0.8);
                // drawWithLineAtMax(h_c_enhanced_charmdecaysOFF_fulljets_c, l, "charm-init full jets", markercolor_c, markerstyle_c, lowx, false, true, 0.8);
                // drawWithLineAtMax(h_b_enhanced_beautydecaysOFF_fulljets_b, l, "beauty-init full jets", markercolor_b, markerstyle_b, lowx, false, true, 0.8);
                drawWithLineAtMax(h_charmdecaysON_fulljets_leadpt5_l, l, "light-init full jets", markercolor_l, kOpenSquare, lowx, false, true, 0.8);
                drawWithLineAtMax(h_c_enhanced_charmdecaysOFF_fulljets_leadpt5_c, l, "charm-init full jets", markercolor_c, kOpenCircle, lowx, false, true, 0.8);
                drawWithLineAtMax(h_b_enhanced_beautydecaysOFF_fulljets_leadpt5_b, l, "beauty-init full jets", markercolor_b, kOpenCross, lowx, false, true, 0.8);
                
                //make ratio plots!
                std::string yaxislabel = "#frac{charm jets}{light jets}";//"#frac{heavy quark jets}{light quark jets}";
                
               
                plotRatio(h_c_enhanced_charmdecaysOFF_fulljets_leadpt5_c, h_charmdecaysON_fulljets_leadpt5_l, "hratio_charm_light" + pt_name + "_R" + jetR + add_name, l2, markercolor_c, markers[2], yaxislabel, c_yaxis_min, c_yaxis_max, 0.61, 0.25, lowx, "", true);
                drawHoriLine(lowx, 1., 0.5, kGray+2, 6)->Draw();
                
                yaxislabel = "#frac{beauty jets}{light jets}";
                plotRatio(h_b_enhanced_beautydecaysOFF_fulljets_leadpt5_b, h_charmdecaysON_fulljets_leadpt5_l, "hratio_beauty_light" + pt_name + "_R" + jetR + add_name, l2, markercolor_b, markers[2], yaxislabel, b_yaxis_min, b_yaxis_max, 0.75, 0.11, lowx);
                drawHoriLine(lowx, 1., 0.5, kGray+2, 6)->Draw();

            } else if (plot_case == 4) {
                pad1 = makeTopPad(0.41);

                double lowx = 1e-4;
                h_c_enhanced_D0->SetMaximum(h_c_enhanced_D0->GetMaximum()*1.2);                
                
                drawNoLine(h_c_enhanced_D0, l, "D^{0}-tagged, c-init jets, jet constituents", kBlack, kOpenSquare, lowx, 1.5, true, 0.8);
                drawNoLine(h_c_enhanced_D0_rigidcone04, l, "#Delta R_{part-axis} < R", kBlue, kFullCircle, lowx, 1.5, true, 0.6);
                drawNoLine(h_c_enhanced_D0_rigidcone1, l, "#Delta R_{part-axis} < 1", kRed, kFullCircle, lowx, 1.2, true, 0.6);

                if (i == 1) drawVertLine(0.4, 0, 0.98, kBlack, 2)->Draw("same");
                else drawVertLine(0.4, 0, 1.4, kBlack, 2)->Draw("same");
                l->Draw("same");

                // now make the ratio plots
                // c_ratio_RLedge->cd();
                std::string yaxislabel = "Ratio   "; //"#frac{charm jets}{light jets}";
                double rc04_yaxis_min = 0.5; double rc04_yaxis_max = 1.5;
                // plotRatio(h_c_enhanced_D0_rigidcone04, h_c_enhanced_D0, "hratio_rigidcone04_jetconst" + pt_name + "_R" + jetR + add_name, l2, kBlue, kFullCircle, 
                //     yaxislabel, rc04_yaxis_min, rc04_yaxis_max, 0.75, 0.11, lowx, "#Delta R_{part-axis} < R / jet constituents");

                // plotRatio(h_c_enhanced_D0_rigidcone1, h_c_enhanced_D0, "hratio_rigidcone1_jetconst" + pt_name + "_R" + jetR + add_name, l2, kRed, kFullCircle, 
                //     yaxislabel, rc04_yaxis_min, rc04_yaxis_max, 0.75, 0.11, lowx, "#Delta R_{part-axis} < 1 / jet constituents");
                
                std::string ratio_name1 = "hratio_rigidcone04_jetconst" + pt_name + "_R" + jetR + add_name;
                std::string ratio_name2 = "hratio_rigidcone1_jetconst" + pt_name + "_R" + jetR + add_name;
                
                plotRatio2(h_c_enhanced_D0_rigidcone04, h_c_enhanced_D0_rigidcone1, h_c_enhanced_D0, ratio_name1, ratio_name2, 
                    kBlue, kFullCircle, 1.5, 0.6, kRed, kFullCircle, 1.2, 0.6, 
                    yaxislabel, rc04_yaxis_min, rc04_yaxis_max, 0.61, 0.11, lowx, 
                    "#Delta R_{part-axis} < R / jet constituents", "#Delta R_{part-axis} < 1 / jet constituents");

            }

            if (plot_case == 11) {
                double lowx = 1e-4;

                drawWithLineAtMax(h_c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons_c, l, "charm-init ch jets + neutral hadrons", markercolor_c, markerstyle_c, lowx, false, false, 0.4);
                drawWithLineAtMax(h_b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons_b, l, "beauty-init ch jets + neutral hadrons", markercolor_b, markerstyle_b, lowx, false, false, 0.4);
                
                drawWithLineAtMax(h_charmdecaysON_l, l, "light-init charged jets", markercolor_l, markerstyle_l_ch, lowx, false, false, 0.8);
                drawWithLineAtMax(h_c_enhanced_charmdecaysOFF_c, l, "charm-init charged jets, charm decays off", markercolor_c, markerstyle_c_ch, lowx, false, false, 0.8);

            } else if (plot_case == 12) {
                double lowx = 1e-4;
                drawWithLineAtMax(h_charmdecaysOFF_fulljets_l, l, "light-init full jets, charm decays off", markercolor_l, markerstyle_l, lowx, false, false, 0.4);
            }

            if (plot_case == 20) {
                double lowx = 1e-4;
                pad1->cd();
                drawWithLineAtMax(h_charmdecaysOFF_l, l, "light-init jets, charm decays off", markercolor_l, markerstyle_l, lowx, false, true, 0.4);
                drawWithLineAtMax(h_charmdecaysOFF_g, l, "gluon-init jets, charm decays off", markercolor_g, markerstyle_g, lowx, false, true, 0.4);
                drawWithLineAtMax(h_charmdecaysOFF_i, l, "inclusive jets, charm decays off", markercolor_i, markerstyle_i, lowx, false, true, 0.4);
            }


            // save peak positions as a function of pT
            if ( plot_case == 100) {
                // pt_x[i] = pt_max;
                // getPeak(hist, ptrl)
                // FULL JETS L/C/B
                peaks_l_full_y[i] = getPeak(h_charmdecaysON_fulljets_l, false);
                peaks_c_full_y[i] = getPeak(h_c_enhanced_charmdecaysOFF_fulljets_c, false);
                peaks_b_full_y[i] = getPeak(h_b_enhanced_beautydecaysOFF_fulljets_b, false);
                peaks_D0_full_y[i] = getPeak(h_c_enhanced_D0_fulljets, false);

                peaks_l_ch_y[i] = getPeak(h_charmdecaysON_l, false);
                peaks_g_ch_y[i] = getPeak(h_charmdecaysON_g, false);
                peaks_incl_ch_y[i] = getPeak(h_charmdecaysON_i, false);
                peaks_D0_ch_y[i] = getPeak(h_c_enhanced_D0, false);

                peakserr_l_full_y[i] = getPeakErr(h_charmdecaysON_fulljets_l, false);
                peakserr_c_full_y[i] = getPeakErr(h_c_enhanced_charmdecaysOFF_fulljets_c, false);
                peakserr_b_full_y[i] = getPeakErr(h_b_enhanced_beautydecaysOFF_fulljets_b, false);
                peakserr_D0_full_y[i] = getPeakErr(h_c_enhanced_D0_fulljets, false);

                peakserr_l_ch_y[i] = getPeakErr(h_charmdecaysON_l, false);
                peakserr_g_ch_y[i] = getPeakErr(h_charmdecaysON_g, false);
                peakserr_incl_ch_y[i] = getPeakErr(h_charmdecaysON_i, false);
                peakserr_D0_ch_y[i] = getPeakErr(h_c_enhanced_D0, false);

            }


            // vector<double> fullwidth_vec = findWidthOfCurve(hD0,  hD0_top_binpos);
            // drawHoriLine(fullwidth_vec[1], fullwidth_vec[2], fullwidth_vec[0], kMagenta+3, 1)->Draw();
            



            // } //end of file loop?

            // draw legend
            l->Draw("same");


            std::string fname = outdir + "QG_comp" + pt_name + "_R" + jetR + add_name; //"_charmdecaysONcomparison.pdf"; // + "_normbytype.pdf"; //"_nonorm.pdf";
            const char* fnamec = fname.c_str();
            if (plot_case != 100 ) c->SaveAs(fnamec);
            delete c;

            f_out->cd();
            h_charmdecaysOFF_g->Write();
            h_charmdecaysOFF_l->Write();
            h_charmdecaysOFF_i->Write();
            h_c_enhanced_charmdecaysOFF_c->Write();
            
            h_ptrl_charmdecaysOFF_g->Write();
            h_ptrl_charmdecaysOFF_l->Write();
            h_ptrl_charmdecaysOFF_i->Write();
            h_ptrl_charmdecaysOFF_c->Write();
            
            h_c_enhanced_D0->Write();
            h_ptrl_c_enhanced_D0->Write();

            h_charmdecaysOFF_fulljets_g->Write();
            h_charmdecaysOFF_fulljets_l->Write();
            h_c_enhanced_charmdecaysOFF_fulljets_c->Write();
            h_b_enhanced_beautydecaysOFF_fulljets_b->Write();

            h_c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons_c->Write();
            h_b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons_b->Write();

            h_charmdecaysON_leadpt5_l->Write();
            h_charmdecaysON_leadpt5_g->Write();
            h_charmdecaysON_leadpt5_i->Write();
            h_c_enhanced_D0_leadpt5_c->Write();


            delete h_charmdecaysOFF_g;
            delete h_charmdecaysOFF_l;
            delete h_charmdecaysOFF_i;
            delete h_c_enhanced_charmdecaysOFF_c;
            
            delete h_ptrl_charmdecaysOFF_g;
            delete h_ptrl_charmdecaysOFF_l;
            delete h_ptrl_charmdecaysOFF_i;
            delete h_ptrl_charmdecaysOFF_c;
            
            delete h_c_enhanced_D0;
            delete h_ptrl_c_enhanced_D0;

            delete h_charmdecaysOFF_fulljets_g;
            delete h_charmdecaysOFF_fulljets_l;
            delete h_c_enhanced_charmdecaysOFF_fulljets_c;
            delete h_b_enhanced_beautydecaysOFF_fulljets_b;

            delete h_c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons_c;
            delete h_b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons_b;

            delete h_charmdecaysON_leadpt5_l;
            delete h_charmdecaysON_leadpt5_g;
            delete h_charmdecaysON_leadpt5_i;
            delete h_c_enhanced_D0_leadpt5_c;
            
             
        } // pT bins loop

        // plot peak positions as a function of pt
        if (plot_case == 100) {
            c_peaks->cd();
            // c_peaks->SetLogy();
            auto mg = new TMultiGraph();
            TLegend* l_peaks = new TLegend(0.5497168,0.550741,0.8562155,0.8785185,""); //dummy legend
            plotGraph(mg, l_peaks, n_bins, pt_x, peaks_l_full_y, peakserr_l_full_y, "light-init full jets", markercolor_l, kOpenSquare);
            plotGraph(mg, l_peaks, n_bins, pt_x, peaks_c_full_y, peakserr_c_full_y, "charm-init full jets", markercolor_c, kOpenCircle);
            plotGraph(mg, l_peaks, n_bins, pt_x, peaks_D0_full_y, peakserr_D0_full_y, "D0-tagged full jets", markercolor_D0, kOpenCircle);
            plotGraph(mg, l_peaks, n_bins, pt_x, peaks_b_full_y, peakserr_b_full_y, "beauty-init full jets", markercolor_b, kOpenCross);
            
            plotGraph(mg, l_peaks, n_bins, pt_x, peaks_l_ch_y, peakserr_l_ch_y, "light-init ch jets", markercolor_l, markerstyle_l, 0.65);
            plotGraph(mg, l_peaks, n_bins, pt_x, peaks_D0_ch_y, peakserr_D0_ch_y, "D0-tagged ch jets", markercolor_D0, markerstyle_D0, 0.65);
            plotGraph(mg, l_peaks, n_bins, pt_x, peaks_g_ch_y, peakserr_g_ch_y, "gluon-init ch jets", markercolor_g, markerstyle_g, 0.65);
            plotGraph(mg, l_peaks, n_bins, pt_x, peaks_incl_ch_y, peakserr_incl_ch_y, "inclusive ch jets", markercolor_i, markerstyle_i, 0.65);
            
            mg->GetXaxis()->SetTitle("#it{p}_{T, jet}");
            mg->GetYaxis()->SetTitle("Most probable #it{R}_{L} value");
        
            //mg->SetMaximum(5.);
            mg->Draw("ap same");
            l_peaks->Draw("same");

            std::string fname_peaks = outdir + "PEAKPOS" + "_R" + jetR + add_name; //"_charmdecaysONcomparison.pdf"; // + "_normbytype.pdf"; //"_nonorm.pdf";
            const char* fname_peaksc = fname_peaks.c_str();
            c_peaks->SaveAs(fname_peaksc);
            delete c_peaks;

            //--------------------------------------
            c_peaks_frac->cd();

            auto mg2 = new TMultiGraph();
            TLegend* l_peaks_frac = new TLegend(0.5497168,0.550741,0.8562155,0.8785185,""); //dummy legend
            l_peaks_frac->SetTextSize(0.028);
            double frac_b_c_y[n_bins];
            double frac_c_l_y[n_bins];
            double frac_b_l_y[n_bins];
            double err_y_dummy[n_bins];
            for (int j=0; j<n_bins; j++) {
                frac_b_c_y[j] = peaks_b_full_y[j]/peaks_c_full_y[j];
                frac_c_l_y[j] = peaks_c_full_y[j]/peaks_l_full_y[j];
                frac_b_l_y[j] = peaks_b_full_y[j]/peaks_l_full_y[j];
                err_y_dummy[j] = 0;

            }
            TFile::Open("MyFile.root", "RECREATE");
            plotGraph(mg2, l_peaks_frac, n_bins, pt_x, frac_b_l_y, err_y_dummy, "#frac{b-init full jets}{l-init full jets}", markercolor_l, markerstyle_l);
            plotGraph(mg2, l_peaks_frac, n_bins, pt_x, frac_b_c_y, err_y_dummy, "#frac{b-init full jets}{c-init full jets}", markercolor_b, markerstyle_b);
            plotGraph(mg2, l_peaks_frac, n_bins, pt_x, frac_c_l_y, err_y_dummy, "#frac{c-init full jets}{l-init full jets}", markercolor_c, markerstyle_c);
            
            mg2->SetMaximum(5.);
            mg2->GetXaxis()->SetTitle("#it{p}_{T, jet}");
            mg2->GetYaxis()->SetTitle("Fraction of most probable #it{R}_{L} value");
        
            mg2->Draw("ap");
            l_peaks_frac->Draw("same");

            std::string fname_peaksfrac = outdir + "PEAKPOS_FRAC" + "_R" + jetR + add_name; //"_charmdecaysONcomparison.pdf"; // + "_normbytype.pdf"; //"_nonorm.pdf";
            const char* fname_peaksfracc = fname_peaksfrac.c_str();
            c_peaks_frac->SaveAs(fname_peaksfracc);
            delete c_peaks_frac;
            
            
            
            
        }


    } // jetR loop

    f_charmdecaysON->Close();
    f_charmdecaysOFF->Close();
    f_charmdecaysON_fulljets->Close();
    f_charmdecaysOFF_fulljets->Close();
    f_ptrl_charmdecaysON->Close();
    f_ptrl_charmdecaysOFF->Close();
    f_ptrl_c_enhanced_D0->Close();

    f_c_enhanced_D0->Close();
    f_c_enhanced_D0_fulljets->Close();
    f_c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons->Close();
    f_b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons->Close();
    f_c_enhanced_charmdecaysOFF_fulljets->Close();
    f_b_enhanced_beautydecaysOFF_fulljets->Close();
    f_c_enhanced_charmdecaysOFF->Close();
    f_c_enhanced_charmdecaysON->Close();

    delete f_charmdecaysON;
    delete f_charmdecaysOFF;
    delete f_charmdecaysON_fulljets;
    delete f_charmdecaysOFF_fulljets;
    delete f_ptrl_charmdecaysON;
    delete f_ptrl_charmdecaysOFF;
    delete f_ptrl_c_enhanced_D0;

    delete f_c_enhanced_D0;
    delete f_c_enhanced_D0_fulljets;
    delete f_c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons;
    delete f_b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons;
    delete f_c_enhanced_charmdecaysOFF_fulljets;
    delete f_b_enhanced_beautydecaysOFF_fulljets;
    delete f_c_enhanced_charmdecaysOFF;
    delete f_c_enhanced_charmdecaysON;



    return;
}
