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

// void FormatHist(TLegend *l, TH1 *hist, TString text, int markercolor=1, int markerstyle=8) 
// {
//     hist->SetLineColor(markercolor);
//     hist->SetMarkerColor(markercolor);
//     hist->SetMarkerStyle(markerstyle);
//     hist->SetMarkerSize(1.5);
//     l->AddEntry(hist, text, "pl");

// 	//gPad->SetTickx(); 
// 	//gPad->SetTicky(); 
// 	// h->SetLineWidth(2);
// 	hist->GetYaxis()->SetTitleOffset(1.05); 
// 	hist->GetYaxis()->SetTitleSize(0.06); //(0.042);
// 	hist->GetYaxis()->SetLabelSize(0.05); //(0.042);
// 	hist->GetYaxis()->SetLabelFont(42);
// 	hist->GetXaxis()->SetLabelFont(42);
// 	hist->GetYaxis()->SetTitleFont(42);
// 	hist->GetXaxis()->SetTitleFont(42);
// 	hist->GetXaxis()->SetTitleOffset(1.0);
// 	hist->GetXaxis()->SetTitleSize(0.06); //(0.042);
// 	hist->GetXaxis()->SetLabelSize(0.05); //(0.042);


//     return;
// }

void FormatHist(TLegend *l, TH1 *hist, TString text, int markercolor=1, int markerstyle=8, double markeralpha=1., double markerfill=-1.,
                double xtitlesize=0.06, double xlabelsize=0.05, double xoffset=1.0,
                double ytitlesize=0.06, double ylabelsize=0.05, double yoffset=1.05) 
                // double xtitlesize=0.04, double xlabelsize=0.04, double xoffset=1.2,
                // double ytitlesize=0.04, double ylabelsize=0.04, double yoffset=1.0) 
{
    hist->SetLineColor(markercolor);
    hist->SetMarkerColorAlpha(markercolor, markeralpha);
    hist->SetMarkerStyle(markerstyle);
    hist->SetMarkerSize(1.5);
    if (markerfill != -1) {
        hist->SetFillStyle(markerfill);
    }
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

    if (ptrl) {
        hist->GetXaxis()->SetTitle("#it{p}_{T}#it{R}_{L}");
    } else {
        hist->GetXaxis()->SetTitle("#it{R}_{L}");
    }
    hist->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");

    return hist;
}

void drawWithLineAtMax(TH1D *hist, TLegend *l, std::string labeltext, int markercolor, int markerstyle, 
                    bool ptrl=false, bool removexaxis=false, double markeralpha=1, double markerfill=-1, bool verbosity=false) {
    // FormatHist(l, hist, "gluon-init jets", markercolor_g, markerstyle_g);
    FormatHist(l, hist, labeltext, markercolor, markerstyle, markeralpha, markerfill);
    hist->Draw("L same");
    if (removexaxis) hist->GetXaxis()->SetLabelSize(0);
    
    double hist_top_binpos = findTopOfCurve(hist, ptrl);
    drawVertLine(hist->GetBinCenter(hist_top_binpos), 0, hist->GetBinContent(hist_top_binpos), markercolor, 1)->Draw();

    if (verbosity) {
        cout << hist->GetName() << " top bincenter " << hist->GetBinCenter(hist_top_binpos) << ", " << hist->GetName() << " top binpos" << hist->GetBinContent(hist_top_binpos) << endl;
    }   
}

void makeTopPad() {
    TPad *pad1 = new TPad("pad1","pad1",0.,0.,1.,1.);
    pad1->SetLogx();
    pad1->SetFillColor(0);
    pad1->SetFillStyle(0);
    pad1->SetBottomMargin(0.31);
    pad1->Draw();
    pad1->cd();
}

void makeBottomPad() {
    TPad *pad2 = new TPad("pad1","",0.,0.,1.,1.);
    pad2->SetTopMargin(0.71);
    pad2->SetFillColor(0);
    pad2->SetFillStyle(0);
    pad2->Draw();
    pad2->SetLogx();
    pad2->cd();
}

void plotRatio(TH1D *h1, TH1D *h2, std::string ratio_name, TLegend *l, int markerstyle, std::string yaxislabel, int ndiv=5) {    
    makeBottomPad();

    TH1D* hratio = (TH1D*) h1->Clone(ratio_name.c_str());
    hratio->Divide(h2);
    // hratio->SetMinimum(0.5);
    // hratio->SetMaximum(1.5);
    
    //xtitlesize, xlabelsize, xoffset, ytitlesize, ylabelsize, yoffset
    FormatHist(l, hratio, "ratio", kBlack, markerstyle, 1.0, 0.05, 0.04, 1.2, 0.035, 0.03, 1.5);
    hratio->GetYaxis()->SetTitle(yaxislabel.c_str());
    hratio->GetYaxis()->SetNdivisions(5);

    hratio->Draw();

    // draw line at 1
    drawHoriLine(1e-4, 1., 1., kGray+2, 1)->Draw();
    // drawHoriLine(1e-4, 1., 0.9, kGray+2, 6)->Draw();
    // drawHoriLine(1e-4, 1., 1.1, kGray+2, 6)->Draw();

}



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

    const TString infile_ptrl_charmdecaysON = basedir + "24691538" + root_filename; //--nocharmdecay 0 --giveptRL 1 //currently not used**
    const TString infile_ptrl_charmdecaysOFF = basedir + "24692158" + root_filename; //--nocharmdecay 1 --giveptRL 1 **
    const TString infile_ptrl_c_enhanced_D0 = basedir + "24690700" + root_filename; //--replaceKP 1 --chinitscat 3 --giveptRL 1**
        
    //---------------------------------------------------------------------//
    // these histograms have D0 pt and z information:
    const TString infile_charmdecaysON = basedir + "25735401" + root_filename; //--nocharmdecay 0 //currently not used**
    const TString infile_charmdecaysOFF = basedir + "25735384" + root_filename; //--nocharmdecay 1**
    
    const TString infile_c_enhanced_D0 = basedir + "25536401" + root_filename; //--replaceKP 1 --chinitscat 3

    const TString infile_c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons = basedir + "25536163" + root_filename; //--nocharmdecay 1 --chinitscat 1 //now --fulljets 2?
    const TString infile_b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons = basedir + "25536177" + root_filename; //--nobeautydecay 1 --chinitscat 5 //now --fulljets 2?
    
    const TString infile_c_enhanced_charmdecaysOFF_fulljets = basedir + "25620729" + root_filename; //--nocharmdecay 1 --chinitscat 1 --fulljets 1
    const TString infile_b_enhanced_beautydecaysOFF_fulljets = basedir + "25620738" + root_filename; //--nobeautydecay 1 --chinitscat 5 --fulljets 1
    
    const TString infile_charmdecaysON_fulljets = basedir + "25735421" + root_filename; //--nocharmdecay 0 --fulljets 1 //currently not used
    const TString infile_charmdecaysOFF_fulljets = basedir + "25735417" + root_filename; //--nocharmdecay 1 --fulljets 1
    
    //---------------------------------------------------------------------//
    // other
    const TString infile_c_enhanced_charmdecaysOFF = basedir + "24735947" + root_filename; //--nocharmdecay 1 --chinitscat 1
    const TString infile_c_enhanced_charmdecaysON = basedir + "24693346" + root_filename; //--nocharmdecay 0 --chinitscat 1 //currently not used

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
    TFile *f_c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons = new TFile(infile_c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons, "READ");
    TFile *f_b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons = new TFile(infile_b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons, "READ");
    TFile *f_c_enhanced_charmdecaysOFF_fulljets = new TFile(infile_c_enhanced_charmdecaysOFF_fulljets, "READ");
    TFile *f_b_enhanced_beautydecaysOFF_fulljets = new TFile(infile_b_enhanced_beautydecaysOFF_fulljets, "READ");
    TFile *f_c_enhanced_charmdecaysOFF = new TFile(infile_c_enhanced_charmdecaysOFF, "READ");
    TFile *f_c_enhanced_charmdecaysON = new TFile(infile_c_enhanced_charmdecaysON, "READ");



    //CONTOL VARIABLES HERE
    // plot cases:
    //final:
    // 0: plot c, g, l, i --NOT IMPLEMENTED
        // g, l, i from charmdecaysOFF; c from c_enhanced_charmdecaysOFF
    // 1: plot D0, g, l, i, include ratio of D0/inclusive **
        // g, l, i from charmdecaysOFF; D0 from c_enhanced_D0
    // 2: plot D0, g, l, i but pt*RL -- TODO: already done in other file, port over here later
        // g, l, i from ptrl_charmdecaysOFF; D0 from ptrl_c_enhanced_D0
    // 3: plot l, c, b full jets **
        // l from charmdecaysOFF_fulljets; c from c_enhanced_charmdecaysOFF_fulljets, b from b_enhanced_beautydecaysOFF_fulljets
    
    //exploring:
    // 10: plot c, g, l, i, b from full jets -- not neccesary...
        // g, l, i from charmdecaysOFF_fulljets; c from c_enhanced_charmdecaysOFF_fulljets; b from b_enhanced_beautydecaysOFF_fulljets
    // 11: plot l, c, b full jets and charged jets + neutrals hadrons, and charged jets
        // g, l, i from charmdecaysOFF_fulljets; c from c_enhanced_charmdecaysOFF_fulljets; b from b_enhanced_beautydecaysOFF_fulljets
        // don't have light for charged jets+neutral, but c from c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons; b from b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons
        // also add light charged jets from infile_charmdecaysOFF, c charged jets from infile_c_enhanced_charmdecaysOFF, don't have b charged jets
    int plot_case = 3;

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
    } else if (plot_case == 3 or plot_case == 11) {
        outdir += "beauty/";
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


        //const int pt_bins[] = { 10, 20, 40, 60, 80, 100, 150 };
        const int pt_bins[] = { 7, 10, 15, 30, 50, 70 }; //, 100, 150 }; //{ 10, 20, 40 };
        const int d0_pt_cuts[] = { 3, 5, 5, 5, 5 }; //, 5, 5 };
        const int n_bins = 5; //7;
        for (int i = 0; i < n_bins; i++) {
            cout << "in pt bin" << i << endl;
            int pt_min = pt_bins[i];
            int pt_max = pt_bins[i+1];

            // define pt related variables
            TString ptbin = TString::Format("%d #leq #it{p}_{T}^{ch. jet} < %d GeV/#it{c}, #font[122]{|}#it{#eta}_{jet}#font[122]{|} #leq 0.5", pt_min, pt_max);
            TString ptD = TString::Format("%d #leq #it{p}_{T}^{D^{0}} < %d GeV/#it{c}, #font[122]{|}#it{y}_{D^{0}}#font[122]{|} #leq 0.8", d0_pt_cuts[i], pt_max);
            

            // make a canvas for each pt range
            TCanvas* c = new TCanvas();
            ProcessCanvas(c);
            c->cd();
            gPad->SetLogx();
            // gPad->SetLogy();

            TLegend* l; // = new TLegend(0.17, 0.65, 0.5, 0.85);
            TLegend* l2 = new TLegend(0.1797168,0.400741,0.4562155,0.8885185,""); //dummy legend

            double maxy = 0;


            // Open histograms


            l = new TLegend(0.1797168,0.400741,0.4562155,0.8885185,""); //(0.17, 0.4, 0.5, 0.53);
            l->SetTextSize(0.045);
            l->AddEntry("NULL","PYTHIA 8 Monash 2013","h");
            l->AddEntry("NULL","pp, #sqrt{#it{s}} = 13 TeV","h");
            if (plot_case == 3 or plot_case == 11) l->AddEntry("NULL","anti-#it{k}_{T}, #it{R} = 0.4","h");
            else {
                l->AddEntry("NULL","D^{0} #rightarrow K^{#minus} #pi^{+} and charge conj.","h");
                l->AddEntry("NULL","in charged jets, anti-#it{k}_{T}, #it{R} = 0.4","h");
            }
            l->AddEntry("NULL",ptbin,"h");
            l->AddEntry("NULL",ptD,"h");
            l->SetTextSize(0.037);
            l->SetBorderSize(0);
            // l->Draw("same");


            std::string pt_name = "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);


            //format: getObsHist(TFile *filename, std::string h_name, std::string h_jet_name, int pt_min, int pt_max, int d0_pt_cut, std::string newhistname, bool d0cuts=false, int obsaxis=3)
            TH1D *h_charmdecaysOFF_g = getObsHist(f_charmdecaysOFF, hg_name, hg_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_charmdecaysOFF_g" + pt_name, false, 4);
            TH1D *h_charmdecaysOFF_l = getObsHist(f_charmdecaysOFF, hl_name, hl_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_charmdecaysOFF_l" + pt_name, false, 4);
            TH1D *h_charmdecaysOFF_i = getObsHist(f_charmdecaysOFF, hi_name, hi_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_charmdecaysOFF_i" + pt_name, false, 4);
            TH1D *h_c_enhanced_charmdecaysOFF_c = getObsHist(f_c_enhanced_charmdecaysOFF, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_c_enhanced_charmdecaysOFF_c" + pt_name, false, 3);
            
            TH1D *h_ptrl_charmdecaysOFF_g = getObsHist(f_ptrl_charmdecaysOFF, hg_name, hg_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_ptrl_charmdecaysOFF_g" + pt_name, false, 3, true);
            TH1D *h_ptrl_charmdecaysOFF_l = getObsHist(f_ptrl_charmdecaysOFF, hl_name, hl_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_ptrl_charmdecaysOFF_l" + pt_name, false, 3, true);
            TH1D *h_ptrl_charmdecaysOFF_i = getObsHist(f_ptrl_charmdecaysOFF, hi_name, hi_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_ptrl_charmdecaysOFF_i" + pt_name, false, 3, true);
            TH1D *h_ptrl_charmdecaysOFF_c = getObsHist(f_ptrl_charmdecaysOFF, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_ptrl_charmdecaysOFF_c" + pt_name, false, 3, true);
            
            TH1D *h_c_enhanced_D0 = getObsHist(f_c_enhanced_D0, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_c_enhanced_D0" + pt_name, true, 4);

            TH1D *h_ptrl_c_enhanced_D0 = getObsHist(f_ptrl_c_enhanced_D0, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_ptrl_c_enhanced_D0" + pt_name, true, 3, true);

            TH1D *h_charmdecaysOFF_fulljets_g = getObsHist(f_charmdecaysOFF_fulljets, hg_name, hg_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_charmdecaysOFF_fulljets_g" + pt_name, false, 4);
            TH1D *h_charmdecaysOFF_fulljets_l = getObsHist(f_charmdecaysOFF_fulljets, hl_name, hl_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_charmdecaysOFF_fulljets_l" + pt_name, false, 4);
            TH1D *h_c_enhanced_charmdecaysOFF_fulljets_c = getObsHist(f_c_enhanced_charmdecaysOFF_fulljets, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_c_enhanced_charmdecaysOFF_fulljets_c" + pt_name, false, 4);
            TH1D *h_b_enhanced_beautydecaysOFF_fulljets_b = getObsHist(f_b_enhanced_beautydecaysOFF_fulljets, hb_name, hb_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_b_enhanced_beautydecaysOFF_fulljets_b" + pt_name, false, 4);

            TH1D *h_c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons_c = getObsHist(f_c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons_c" + pt_name, false, 4);
            TH1D *h_b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons_b = getObsHist(f_b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons, hb_name, hb_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons_b" + pt_name, false, 4);



            // Rebin
            int n_obs_bins = 50; //-1;
            double obs_bins[] = {1.00000000e-04, 1.25892541e-04, 1.58489319e-04, 1.99526231e-04,
                2.51188643e-04, 3.16227766e-04, 3.98107171e-04, 5.01187234e-04,
                6.30957344e-04, 7.94328235e-04, 1.00000000e-03, 1.25892541e-03,
                1.58489319e-03, 1.99526231e-03, 2.51188643e-03, 3.16227766e-03,
                3.98107171e-03, 5.01187234e-03, 6.30957344e-03, 7.94328235e-03,
                1.00000000e-02, 1.25892541e-02, 1.58489319e-02, 1.99526231e-02,
                2.51188643e-02, 3.16227766e-02, 3.98107171e-02, 5.01187234e-02,
                6.30957344e-02, 7.94328235e-02, 1.00000000e-01, 1.25892541e-01,
                1.58489319e-01, 1.99526231e-01, 2.51188643e-01, 3.16227766e-01,
                3.98107171e-01, 5.01187234e-01, 6.30957344e-01, 7.94328235e-01,
                1.00000000e+00};


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



            // // Format histograms for plotting (this order needed to keep legend in order and graphs lookin good)

            // if (charmdecays) {
            //     l->AddEntry("NULL","          D* decays on","h");
            // } else {
            //     l->AddEntry("NULL","          D* decays off","h");
            // }




            // FormatHist(l, h_ptrl_charmdecaysOFF_c, "charm-init jets", markercolor_c, markerstyle_c);

            // FormatHist(l, h_charmdecaysOFF_fulljets_g, "gluon-init full jets", markercolor_g, markerstyle_g);



            if (plot_case == 0) {

            } else if (plot_case == 1) { //D0, g, l, i in charged jets

                makeTopPad();
                h_charmdecaysOFF_l->SetMaximum(h_charmdecaysOFF_l->GetMaximum()*1.2);
                
                drawWithLineAtMax(h_charmdecaysOFF_l, l, "light-init jets", markercolor_l, markerstyle_l, false, true, 0.8);
                drawWithLineAtMax(h_charmdecaysOFF_g, l, "gluon-init jets", markercolor_g, markerstyle_g, false, true, 0.8);
                drawWithLineAtMax(h_charmdecaysOFF_i, l, "inclusive jets", markercolor_i, markerstyle_i, false, true, 0.8);
                drawWithLineAtMax(h_c_enhanced_D0, l, "D^{0}-tagged, c-init jets", markercolor_D0, markerstyle_D0, false, true, 0.8);

                //make ratio plot!
                std::string yaxislabel = "#frac{D^{0}-tagged jets}{inclusive jets}";
                plotRatio(h_c_enhanced_D0, h_charmdecaysOFF_i, "hratio" + pt_name + "_R" + jetR + add_name, l2, markers[2], yaxislabel);



            } else if (plot_case == 2) { //D0, g, l, i pT*RL
                h_ptrl_charmdecaysOFF_l->SetMaximum(h_ptrl_charmdecaysOFF_l->GetMaximum()*1.2);

                drawWithLineAtMax(h_ptrl_charmdecaysOFF_l, l, "light-init jets", markercolor_l, markerstyle_l, true, false, 0.8);
                drawWithLineAtMax(h_ptrl_charmdecaysOFF_g, l, "gluon-init jets", markercolor_g, markerstyle_g, true, false, 0.8);
                drawWithLineAtMax(h_ptrl_charmdecaysOFF_i, l, "inclusive jets", markercolor_i, markerstyle_i, true, false, 0.8);
                // drawWithLineAtMax(h_ptrl_c_enhanced_charmdecaysOFF_chargedjets_c, l, "charm-init jets", markercolor_c, markerstyle_c, true, false, 0.8);
                drawWithLineAtMax(h_ptrl_c_enhanced_D0, l, "D^{0}-tagged, c-init jets", markercolor_D0, markerstyle_D0, true, false, 0.8);

                // TH1D* h_ptrl_c_enhanced_D0
            } else if (plot_case == 3 or plot_case == 11) { //l, c, b
                h_charmdecaysOFF_fulljets_l->SetMaximum(h_charmdecaysOFF_fulljets_l->GetMaximum()*1.2);

                drawWithLineAtMax(h_charmdecaysOFF_fulljets_l, l, "light-init full jets", markercolor_l, markerstyle_l, false, false, 0.8);
                drawWithLineAtMax(h_c_enhanced_charmdecaysOFF_fulljets_c, l, "charm-init full jets", markercolor_c, markerstyle_c, false, false, 0.8);
                drawWithLineAtMax(h_b_enhanced_beautydecaysOFF_fulljets_b, l, "beauty-init full jets", markercolor_b, markerstyle_b, false, false, 0.8);
                
            }

            if (plot_case == 11) {

                drawWithLineAtMax(h_c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons_c, l, "charm-init ch jets + neutral hadrons", markercolor_c, markerstyle_c, false, false, 0.4);
                drawWithLineAtMax(h_b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons_b, l, "beauty-init ch jets + neutral hadrons", markercolor_b, markerstyle_b, false, false, 0.4);
                
                drawWithLineAtMax(h_charmdecaysOFF_l, l, "light-init charged jets", markercolor_l, markerstyle_l_ch, false, false, 0.8);
                drawWithLineAtMax(h_c_enhanced_charmdecaysOFF_c, l, "charm-init charged jets", markercolor_c, markerstyle_c_ch, false, false, 0.8);

            }



            
            // vector<double> fullwidth_vec = findWidthOfCurve(hD0,  hD0_top_binpos);
            // drawHoriLine(fullwidth_vec[1], fullwidth_vec[2], fullwidth_vec[0], kMagenta+3, 1)->Draw();
            
            
            // Add legend about D0 info
            


            // TODO: ADD RATIO PLOT
            // auto rp = new TRatioPlot(hD0, hc);
            // rp->Draw();
            // rp->GetLowYaxis()->SetNdivisions(505);
            // c->Update();


    


            



            // } //end of file loop?

            // draw legend
            l->Draw("same");


            std::string fname = outdir + "QG_comp" + pt_name + "_R" + jetR + add_name; //"_charmdecaysONcomparison.pdf"; // + "_normbytype.pdf"; //"_nonorm.pdf";
            const char* fnamec = fname.c_str();
            c->SaveAs(fnamec);
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
            
             
        } // pT bins loop
    } // jetR loop

    f_charmdecaysON->Close();
    f_charmdecaysOFF->Close();
    f_charmdecaysON_fulljets->Close();
    f_charmdecaysOFF_fulljets->Close();
    f_ptrl_charmdecaysON->Close();
    f_ptrl_charmdecaysOFF->Close();
    f_ptrl_c_enhanced_D0->Close();

    f_c_enhanced_D0->Close();
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
    delete f_c_enhanced_charmdecaysOFF_chargedjets_and_neutralhadrons;
    delete f_b_enhanced_beautydecaysOFF_chargedjets_and_neutralhadrons;
    delete f_c_enhanced_charmdecaysOFF_fulljets;
    delete f_b_enhanced_beautydecaysOFF_fulljets;
    delete f_c_enhanced_charmdecaysOFF;
    delete f_c_enhanced_charmdecaysON;



    return;
}
