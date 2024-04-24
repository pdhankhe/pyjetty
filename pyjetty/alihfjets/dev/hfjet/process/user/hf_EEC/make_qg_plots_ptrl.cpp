// ROOT macro to make quark-gluon jet plots
// this one specifically for plotting&saving PT*RL plots from THnSparse
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

void FormatHist(TLegend *l, TH1 *hist, TString text, int markercolor=1, int markerstyle=8, double markeralpha=1) 
{
    hist->SetLineColor(markercolor);
    hist->SetMarkerColorAlpha(markercolor, markeralpha);
    hist->SetMarkerStyle(markerstyle);
    hist->SetMarkerSize(1.5);
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

void FormatHistwithLine(TLegend *l, TH1 *hist, TString text, int linecolor=1, int linestyle=1, double linealpha=1) 
{
    hist->SetMarkerStyle(20);
    hist->SetMarkerColorAlpha(linecolor, 0);

    hist->SetFillStyle(0);
    hist->SetLineColorAlpha(linecolor, linealpha);
    hist->SetFillColor(linecolor);
    // hist->SetLineStyle(linestyle);
    hist->SetLineWidth(3);
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


void make_qg_plots_ptrl() {

//    gROOT->SetBatch(); //prevents plots from showing up
    gStyle->SetOptStat(0);
    SetStyle();
    Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
    Double_t marker_size = 1.5;
    Double_t colors[16] = {kRed, kGreen+2, kBlue, kRed+1, kGreen+1, kBlue+1, kRed+2, kGreen+2, kBlue+2, kRed+3, kGreen+3, kBlue+3, kOrange+1, kViolet+1, kYellow+1, kCyan+1};


    // we want a couple of plots
    // 1. D-jet vs g- vs l- vs i- (decays on) - one plot for each pt ranges
    // 2. D-jet comparison of pt ranges
    // so we need two files: the one with ptrl for dif parton-initiated jets, and one with D-jet ptrl info

     //FOR WHEN WEIGHTED/UNWEIGHTED IN SAME FILE
    const char infile_ptrl_charmdecaysON[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/24691538/AnalysisResultsFinal.root";
    const char infile_ptrl_D0[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/24690700/AnalysisResultsFinal.root"; 
    

    TFile* f_parton = new TFile(infile_ptrl_charmdecaysON, "READ");
    TFile* f_D0 = new TFile(infile_ptrl_D0, "READ");
    // std::vector<TFile*> files;

    TString label1 = "";
    TString label2 = "";  


    // Output directory
    std::string outdir = "plots/final/ptrl/";//"plots/test/";
    // Output file for binned results
    std::string outfile = outdir + "AnalysisResultsFinal_ptrl.root";
    TFile* f_out = new TFile(outfile.c_str(), "RECREATE");

    // Jet r value
    std::string jetR_list[] = { "0.4" };
    for (std::string jetR : jetR_list) {

        // Names of histograms in the file (quark, charm, gluon)
        const std::string hc_name = "h_EEC_JetPt_charm_R" + jetR;
        const std::string hl_name = "h_EEC_JetPt_light_R" + jetR;
        const std::string hg_name = "h_EEC_JetPt_gluon_R" + jetR;
        const std::string hi_name = "h_EEC_JetPt_inclusive_R" + jetR;
        const std::string hc_jet_name = "h_JetPt_charm_R" + jetR + "_jetlevel";
        const std::string hl_jet_name = "h_JetPt_light_R" + jetR + "_jetlevel";
        const std::string hg_jet_name = "h_JetPt_gluon_R" + jetR + "_jetlevel";
        const std::string hi_jet_name = "h_JetPt_inclusive_R" + jetR + "_jetlevel";

        const std::string hD0KpiNjets_name = "hD0KpiNjets"; //TODO later: = "hD0KpiNehD0KpiNjetsvents" for run 16729583

        // set up variables across pts
        vector<TH1D*> hD0_vec;
        Double_t pt_D0_colors[] = { kRed+2, kRed, kRed-9 };
        // const TString pt_labels[] = {}
        TCanvas* c_D0_allpt = new TCanvas();
        ProcessCanvas(c_D0_allpt);
        gPad->SetLogx();

        TLegend* l_2 = new TLegend(0.1797168,0.400741,0.4562155,0.8885185,"");
        l_2->SetTextSize(0.045);
        l_2->AddEntry("NULL","PYTHIA 8 Monash 2013","h");
        l_2->AddEntry("NULL","pp, #sqrt{#it{s}} = 13 TeV","h");
        l_2->AddEntry("NULL","D^{0} #rightarrow K^{#minus} #pi^{+} and charge conj.","h");
        l_2->AddEntry("NULL","in charged jets, anti-#it{k}_{T}, #it{R} = 0.4","h");
        // l_2->AddEntry("NULL",ptbin,"h");
        // l_2->AddEntry("NULL",ptD,"h");
        l_2->SetTextSize(0.037);
        l_2->SetBorderSize(0);


        //const int pt_bins[] = { 10, 20, 40, 60, 80, 100, 150 };
        const int pt_bins[] = { 7, 10, 15, 30 }; //{ 10, 20, 40 };
        const int d0_pt_cuts[] = { 3, 5, 5 };
        const int n_bins = 3;
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
            // c->cd();
            gPad->SetLogx();
            // gPad->SetLogy();

            TLegend* l; // = new TLegend(0.17, 0.65, 0.5, 0.85);
            

            double maxy = 0;


            // Open histograms


            l = new TLegend(0.1797168,0.400741,0.4562155,0.8885185,""); //(0.17, 0.4, 0.5, 0.53);
            l->SetTextSize(0.045);
            // TLegend *leg = new TLegend(0.1797168,0.5390741,0.4562155,0.8885185,"");
            l->AddEntry("NULL","PYTHIA 8 Monash 2013","h");
            l->AddEntry("NULL","pp, #sqrt{#it{s}} = 13 TeV","h");
            l->AddEntry("NULL","D^{0} #rightarrow K^{#minus} #pi^{+} and charge conj.","h");
            l->AddEntry("NULL","in charged jets, anti-#it{k}_{T}, #it{R} = 0.4","h");
            l->AddEntry("NULL",ptbin,"h");
            l->AddEntry("NULL",ptD,"h");
            l->SetTextSize(0.037);
            l->SetBorderSize(0);
            // l->Draw("same");
            

            //-------------------------------------------------//
            // find D0 reconstruction through charm
            THnSparse* hsparsejet_c = (THnSparse*) f_D0->Get(hc_name.c_str());
            THnSparse* hsparsejet_c_jetlevel = (THnSparse*) f_D0->Get(hc_jet_name.c_str());
            THnSparse* hsparsejet_g = (THnSparse*) f_parton->Get(hg_name.c_str());
            THnSparse* hsparsejet_g_jetlevel = (THnSparse*) f_parton->Get(hg_jet_name.c_str());
            THnSparse* hsparsejet_l = (THnSparse*) f_parton->Get(hl_name.c_str());
            THnSparse* hsparsejet_l_jetlevel = (THnSparse*) f_parton->Get(hl_jet_name.c_str());
            THnSparse* hsparsejet_i = (THnSparse*) f_parton->Get(hi_name.c_str());
            THnSparse* hsparsejet_i_jetlevel = (THnSparse*) f_parton->Get(hi_jet_name.c_str());
            


            // testing - look at # jets before cuts
            // cout << "numDtaggedjets from hist before cuts " << hsparsejet_c_jetlevel->Projection(0)->GetEntries() << endl;

            // for THnSparse: make clone to work with, make cuts, get projection
            cout << "checkpoint 0" << endl;
            THnSparse *hsparsejet_c_clone = (THnSparse *) hsparsejet_c->Clone("hsparsejet_c_clone");
            THnSparse *hsparsejet_c_jetlevel_clone = (THnSparse *) hsparsejet_c_jetlevel->Clone("hsparsejet_c_jetlevel_clone");
            THnSparse *hsparsejet_g_clone = (THnSparse *) hsparsejet_g->Clone("hsparsejet_g_clone");
            THnSparse *hsparsejet_g_jetlevel_clone = (THnSparse *) hsparsejet_g_jetlevel->Clone("hsparsejet_g_jetlevel_clone");
            THnSparse *hsparsejet_l_clone = (THnSparse *) hsparsejet_l->Clone("hsparsejet_l_clone");
            THnSparse *hsparsejet_l_jetlevel_clone = (THnSparse *) hsparsejet_l_jetlevel->Clone("hsparsejet_l_jetlevel_clone");
            THnSparse *hsparsejet_i_clone = (THnSparse *) hsparsejet_i->Clone("hsparsejet_i_clone");
            THnSparse *hsparsejet_i_jetlevel_clone = (THnSparse *) hsparsejet_i_jetlevel->Clone("hsparsejet_i_jetlevel_clone");

            cout << "checkpoint 1" << endl;

            // get jet pT range - no D0 reconstruction, so don't make D0 cuts
            hsparsejet_c_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
            hsparsejet_c_jetlevel_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
            hsparsejet_g_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max);
            hsparsejet_g_jetlevel_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max);
            hsparsejet_l_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max);
            hsparsejet_l_jetlevel_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max);
            hsparsejet_i_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max);
            hsparsejet_i_jetlevel_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max);

            //make additional D0 cuts
            hsparsejet_c_clone->GetAxis(1)->SetRangeUser(d0_pt_cuts[i], pt_max); // apply cut on Dmeson pt
            hsparsejet_c_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
            hsparsejet_c_jetlevel_clone->GetAxis(1)->SetRangeUser(d0_pt_cuts[i], pt_max);
            hsparsejet_c_jetlevel_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8);
            
            cout << "checkpoint 2" << endl;


            // Project onto observable axis
            TH1D *hc_proj = hsparsejet_c_clone->Projection(3); //CALL THESE TH1*????
            TH1D *hc1D_jet = hsparsejet_c_jetlevel_clone->Projection(0);
            TH1D *hg_proj = hsparsejet_g_clone->Projection(3);
            TH1D *hg1D_jet = hsparsejet_g_jetlevel_clone->Projection(0);
            TH1D *hl_proj = hsparsejet_l_clone->Projection(3);
            TH1D *hl1D_jet = hsparsejet_l_jetlevel_clone->Projection(0);
            TH1D *hi_proj = hsparsejet_i_clone->Projection(3);
            TH1D *hi1D_jet = hsparsejet_i_jetlevel_clone->Projection(0);
            cout << "checkpoint 3" << endl;

            // Set to appropriate name
            std::string hname = hc_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hc_proj->SetNameTitle(hname.c_str(), hname.c_str());

            hname = hg_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hg_proj->SetNameTitle(hname.c_str(), hname.c_str());

            hname = hl_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hl_proj->SetNameTitle(hname.c_str(), hname.c_str());

            hname = hi_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hi_proj->SetNameTitle(hname.c_str(), hname.c_str());

            cout << "checkpoint 4" << endl;


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

            TH1D* hD0 = (TH1D*) hc_proj->Clone((hc_name + "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max)).c_str());
            TH1D* hg = (TH1D*) hg_proj->Clone((hg_name + "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max)).c_str());
            TH1D* hl = (TH1D*) hl_proj->Clone((hl_name + "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max)).c_str());
            TH1D* hi = (TH1D*) hi_proj->Clone((hi_name + "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max)).c_str());
            

            // Find normalization factor
            double numjets_charm = hc1D_jet->Integral();
            double numjets_gluon = hg1D_jet->Integral();
            double numjets_light = hl1D_jet->Integral();
            double numjets_inclusive = hi1D_jet->Integral();

            // Set normalization
            hD0->Scale(1/numjets_charm, "width");
            hg->Scale(1/numjets_gluon, "width");
            hl->Scale(1/numjets_light, "width");
            hi->Scale(1/numjets_inclusive, "width");


            // // Find maximum
            // maxy = hi->GetMaximum() * 1.1;
            // hi->SetMaximum(maxy); 



            //Format color and style
            int markercolor1 = kRed; //charm
            int markerstyle1 = kFullCircle;
            int markercolor2 = kViolet+2; //gluon
            int markerstyle2 = 33;
            int markercolor3 = kGreen+2; //light
            int markerstyle3 = 21;
            int markercolor4 = kMagenta-6; //inclusive
            int markerstyle4 = 29; //star
            label1 = "D0-tagged jets";
            label2 = "gluon-init jets";
            TString label3 = "light-init jets";
            TString label4 = "inclusive jets";

            // Format histograms for plotting (this order needed to keep legend in order and graphs lookin good)
            hD0->GetXaxis()->SetTitle("#it{p}_{T}#it{R}_{L}");
            hD0->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
            hg->GetXaxis()->SetTitle("#it{p}_{T}#it{R}_{L}");
            hg->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
            hl->GetXaxis()->SetTitle("#it{p}_{T}#it{R}_{L}");
            hl->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
            hi->GetXaxis()->SetTitle("#it{p}_{T}#it{R}_{L}");
            hi->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
            
            cout << "about to format plots" << endl;
            FormatHist(l, hD0, label1, markercolor1, markerstyle1);
            FormatHist(l, hg, label2, markercolor2, markerstyle2);
            FormatHist(l, hl, label3, markercolor3, markerstyle3);
            FormatHist(l, hi, label4, markercolor4, markerstyle4);
            
            hg->SetMaximum(hg->GetMaximum()*1.2);

            hg->Draw("L same");
            hl->Draw("L same");
            hi->Draw("L same");
            hD0->Draw("L same");

            // draw legend
            l->Draw("same");
                        

            
            double hD0_top_binpos = findTopOfCurve(hD0, true);
            cout << "hD0 top bincenter" << hD0->GetBinCenter(hD0_top_binpos) << "hD0 top binpos" << hD0->GetBinContent(hD0_top_binpos) << endl;
            drawVertLine(hD0->GetBinCenter(hD0_top_binpos), 0, hD0->GetBinContent(hD0_top_binpos), markercolor1, 1)->Draw();
            double hg_top_binpos = findTopOfCurve(hg, true);
            cout << "hg top bincenter" << hg->GetBinCenter(hg_top_binpos) << "hg top binpos" << hg->GetBinContent(hg_top_binpos) << endl;
            drawVertLine(hg->GetBinCenter(hg_top_binpos), 0, hg->GetBinContent(hg_top_binpos), markercolor2, 1)->Draw();
            double hl_top_binpos = findTopOfCurve(hl, true);
            cout << "hl top bincenter" << hl->GetBinCenter(hl_top_binpos) << "hl top binpos" << hl->GetBinContent(hl_top_binpos) << endl;
            drawVertLine(hl->GetBinCenter(hl_top_binpos), 0, hl->GetBinContent(hl_top_binpos), markercolor3, 1)->Draw();
            double hi_top_binpos = findTopOfCurve(hi, true);
            cout << "hi top bincenter" << hi->GetBinCenter(hi_top_binpos) << "hi top binpos" << hi->GetBinContent(hi_top_binpos) << endl;
            drawVertLine(hi->GetBinCenter(hi_top_binpos), 0, hi->GetBinContent(hi_top_binpos), markercolor4, 1)->Draw();


            
            // vector<double> fullwidth_vec = findWidthOfCurve(hD0,  hD0_top_binpos);
            // drawHoriLine(fullwidth_vec[1], fullwidth_vec[2], fullwidth_vec[0], kMagenta+3, 1)->Draw();
            
            
            

            // make ratio plot
            // auto rp = new TRatioPlot(hD0, hc);
            // rp->Draw();
            // rp->GetLowYaxis()->SetNdivisions(505);
            // c->Update();




            std::string fname = outdir + "QG_comp_pt" + std::to_string(pt_min) + '-' + std::to_string(pt_max) + "_R" + jetR + "Dgli_ptrl_comparison.pdf";
            const char* fnamec = fname.c_str();
            c->SaveAs(fnamec);
            delete c;


            // now save to D0 pt comparison plot
            TString D0_label = std::to_string(pt_min) + " #leq #it{p}_{T}^{ch. jet} < " + std::to_string(pt_max);
            hD0_vec.push_back( (TH1D*) hD0->Clone( hD0->GetName() ) );
            std::cout << "adding to legend" << D0_label << endl;
            for (int j=0; j < hD0_vec[i]->GetNbinsX();j++){
                hD0_vec[i]->SetBinError(j+1, 0);
            }
            FormatHistwithLine(l_2, hD0_vec[i], D0_label, pt_D0_colors[i], markerstyle1, 0.60);
            c_D0_allpt->cd();
            hD0_vec[i]->Draw("L, same");
            if (i==2) {
                l_2->Draw();
            }
            


            f_out->cd();
            hD0->Write();
            hg->Write();
            hl->Write();
            hi->Write();

            delete hD0;
            delete hl;
            delete hg;
            delete hi;
            delete hc1D_jet;
            delete hg1D_jet;
            delete hl1D_jet;
            delete hi1D_jet;
            
             
        } // pT bins loop

        std::string fname_D = outdir + "QG_comp_R" + jetR + "D-tagged_ptrl_comparison.pdf";
        const char* fname_Dc = fname_D.c_str();
        c_D0_allpt->SaveAs(fname_Dc);
        delete c_D0_allpt;

    } // jetR loop

    f_parton->Close();
    delete f_parton;
    f_D0->Close();
    delete f_D0;


    return;
}
