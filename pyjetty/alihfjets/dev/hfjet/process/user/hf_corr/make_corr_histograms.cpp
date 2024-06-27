// ROOT macro to make plots of other correlations
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

void FormatHist(TLegend *l, TH1 *hist, TString text, int markercolor=1, int markerstyle=8) 
{
    hist->SetLineColor(markercolor);
    hist->SetMarkerColor(markercolor);
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


void make_corr_histograms() {

//    gROOT->SetBatch(); //prevents plots from showing up
    gStyle->SetOptStat(0);
    SetStyle();
    Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
    Double_t marker_size = 1.5;
    Double_t colors[16] = {kRed, kGreen+2, kBlue, kOrange+1, kViolet+1, kYellow+1, kCyan+1};
    // Double_t colors[16] = {kRed, kGreen+2, kBlue, kRed+1, kGreen+1, kBlue+1, kRed+2, kGreen+2, kBlue+2, kRed+3, kGreen+3, kBlue+3, kOrange+1, kViolet+1, kYellow+1, kCyan+1};

    const int pt_bins[] = { 20, 40, 60, 80 };
    const std::string pt_bin_names[] = {"20-40", "40-60", "60-80"};
    const double RL_bins[] = { 1e-2, 2e-2, 3e-2, 7.5e-2, 2e-1, 4e-1};
    // const double RL_bins_end[] = { 2*10e-2, 3*10e-2, 7.5*10e-2, 2*10e-1, 4*10e-1 };
    const int n_bins = 3;
    const int n_RLbins = 5;



    // Files
    // const char infile[] = "/software/users/blianggi/mypyjetty/analysis/output/100k/AnalysisResultsFinal.root"; //hiccup
    const char infile[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/pythiagen/byhand/23204141/AnalysisResultsFinal_1percent.root"; //perlmutter, before june 2024
    // const char infile[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/pythiagen/scaling/26652369??/AnalysisResultsFinal.root"; //perlmutter, after june 2024
    // const char infile[] = "/software/users/blianggi/mypyjetty/analysis/output/AnalysisResults.root";
    // std::string indir = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/pythiagen/23121850/";

    // auto tree1 = new TChain("tree1");
    // for (int ifile=1; ifile<=10; file++) {
    //     std::string filename = indir + std::to_string(ifile) + "/AnalysisResultsFastSim_batch.root";
    //     tree1->Add(filename.c_str());
    // }

    //CONTOL VARIABLES HERE
    // bool charmdecays = false; // true = charm decays on, false = charm decays off
    // bool ptrl = true; // true = Plot the pT*RL

    TFile* f;
    std::string add_name;

    TString label1 = "";
    TString label2 = "";  
    vector<TString> label;
    
    f = new TFile(infile, "READ");
    add_name = "_othercorrel";
    // if (ptrl) {
    //     if (charmdecays) {
    //         f = new TFile(infile_ptrl_charmdecaysON, "READ");
    //         add_name = "_cgl_ptrl_charmdecaysON.pdf";
    //     } else {
    //         f = new TFile(infile_ptrl_charmdecaysOFF, "READ");
    //         add_name = "_cgl_ptrl_charmdecaysOFF.pdf";
    //     }
    // } else {
    //     cout << "This didn't work because no files are defined" << endl;
    // }

    cout << "output name will be " << add_name << endl;

    // Jet r value
    std::string jetR_list[] = { "0.4" };
    for (std::string jetR : jetR_list) {

        // Names of histograms in the file (quark, charm, gluon)
        const std::string deltap_truth_tt015_name = "h_corr_deltap2_JetPt_Truth_R0.4_0.15";
        const std::string deltap_truth_tt005_name = "h_corr_deltap2_JetPt_Truth_R0.4_0.05";
        const std::string deltap_truth_tt1_name = "h_corr_deltap2_JetPt_Truth_R0.4_1.0";

        const std::string deltapt_truth_tt015_name = "h_corr_deltapt2_JetPt_Truth_R0.4_0.15";
        const std::string deltapt_truth_tt005_name = "h_corr_deltapt2_JetPt_Truth_R0.4_0.05";
        const std::string deltapt_truth_tt1_name = "h_corr_deltapt2_JetPt_Truth_R0.4_1.0";

        const std::string oppcharge_truth_tt015_name = "h_corr_oppcharge2_JetPt_Truth_R0.4_0.15";
        const std::string oppcharge_truth_tt005_name = "h_corr_oppcharge2_JetPt_Truth_R0.4_0.05";
        const std::string oppcharge_truth_tt1_name = "h_corr_oppcharge2_JetPt_Truth_R0.4_1.0";

        const std::string samecharge_truth_tt015_name = "h_corr_samecharge2_JetPt_Truth_R0.4_0.15";
        const std::string samecharge_truth_tt005_name = "h_corr_samecharge2_JetPt_Truth_R0.4_0.05";
        const std::string samecharge_truth_tt1_name = "h_corr_samecharge2_JetPt_Truth_R0.4_1.0";

        const std::string unweightedRL_truth_tt015_name = "h_corr_unweightedRL2_JetPt_Truth_R0.4_0.15";
        const std::string unweightedRL_truth_tt005_name = "h_corr_unweightedRL2_JetPt_Truth_R0.4_0.05";
        const std::string unweightedRL_truth_tt1_name = "h_corr_unweightedRL2_JetPt_Truth_R0.4_1.0";

        const std::string jet_pt_truth_tt1_name = "h_jet_pt_JetPt_Truth_R0.4_1.0";


        for (int i = 0; i < n_bins; i++) {
            cout << "in pt bin" << i << endl;
            int pt_min = pt_bins[i];
            int pt_max = pt_bins[i+1];
            std::string ptname = std::to_string(pt_min) + '-' + std::to_string(pt_max);

            // Output directory
            std::string outdir = "plots/secondattempt/"+pt_bin_names[i]+"/";//"plots/test/";
            // Output file for binned results
            std::string outfile = outdir + "AnalysisResultsFinal" + add_name + ".root";
            TFile* f_out = new TFile(outfile.c_str(), "RECREATE");

            // define pt related variables
            TString ptbin = TString::Format("%d #leq #it{p}_{T}^{ch. jet} < %d GeV/#it{c}, #font[122]{|}#it{#eta}_{jet}#font[122]{|} #leq 0.5", pt_min, pt_max);
            TString ptD = TString::Format("5 #leq #it{p}_{T}^{D^{0}} < %d GeV/#it{c}, #font[122]{|}#it{y}_{D^{0}}#font[122]{|} #leq 0.8", pt_max);
            

            TCanvas* c_deltap_all = new TCanvas();
            ProcessCanvas(c_deltap_all);
            // c_deltap_all->SetMaximum(1.5);

            TCanvas* c_deltapt_all = new TCanvas();
            ProcessCanvas(c_deltapt_all);
            // c_deltapt_all->SetMaximum(1.5);

            TCanvas* c_charge_all = new TCanvas();
            ProcessCanvas(c_charge_all);
            gPad->SetLogx();

            TLegend* l; // = new TLegend(0.17, 0.65, 0.5, 0.85);
            l = new TLegend(0.1797168,0.650741,0.4562155,0.8885185,""); //(0.17, 0.4, 0.5, 0.53);
            l->SetTextSize(0.037);
            l->SetBorderSize(0);

            TLegend* l2; // = new TLegend(0.17, 0.65, 0.5, 0.85);
            l2 = new TLegend(0.1797168,0.700741,0.4562155,0.8885185,""); //(0.17, 0.4, 0.5, 0.53);
            // l2->SetTextSize(0.037);
            // l2->SetBorderSize(0);

            TLegend* l_right = new TLegend(0.75, 0.5, 0.85, 0.87);
            l_right->SetTextSize(0.037);
            l_right->SetBorderSize(0);

            // make histogram array bc otherwise they get overwritten and won't save S M H
            vector<TH1D*> hcorr_deltap_truth_tt1_arr;
            vector<TH1D*> hcorr_deltapt_truth_tt1_arr;
            vector<TH1D*> hcorr_oppcharge_truth_tt1_arr;
            vector<TH1D*> hcorr_samecharge_truth_tt1_arr;
            vector<TH1D*> hcorr_chargeratio_truth_tt1_arr;
            TH1D *hdummyRL;
            TH1D *hdummyRL2;
                
            for (int j=0; j < n_RLbins; j++) {

                double RL_min = RL_bins[j];
                double RL_max = RL_bins[j+1];
                cout << "in RL bin" << j << " with " << RL_min << " - " << RL_max << endl;

                
                double maxy = 0;


                // Open histograms
                

                //-------------------------------------------------//
                // find D0 reconstruction through charm
                THnSparse* hsparsejet_deltap_truth_tt1 = (THnSparse*) f->Get(deltap_truth_tt1_name.c_str());
                THnSparse* hsparsejet_deltapt_truth_tt1 = (THnSparse*) f->Get(deltapt_truth_tt1_name.c_str());
                THnSparse* hsparsejet_oppcharge_truth_tt1 = (THnSparse*) f->Get(oppcharge_truth_tt1_name.c_str());
                THnSparse* hsparsejet_samecharge_truth_tt1 = (THnSparse*) f->Get(samecharge_truth_tt1_name.c_str());
                // THnSparse* hsparsejet_unweightedRL_truth_tt1 = (THnSparse*) f->Get(unweightedRL_truth_tt1_name.c_str());

                TH1D* h_pt_tt1_jetlevel = (TH1D*) f->Get(jet_pt_truth_tt1_name.c_str());
                // cout << "checkpoint 0" << endl;


                // for THnSparse: make clone to work with, make cuts, get projection
                THnSparse *hsparsejet_deltap_truth_tt1_clone = (THnSparse *) hsparsejet_deltap_truth_tt1->Clone("hsparsejet_deltap_truth_tt1_clone");
                THnSparse *hsparsejet_deltapt_truth_tt1_clone = (THnSparse *) hsparsejet_deltapt_truth_tt1->Clone("hsparsejet_deltapt_truth_tt1_clone");
                THnSparse *hsparsejet_oppcharge_truth_tt1_clone = (THnSparse *) hsparsejet_oppcharge_truth_tt1->Clone("hsparsejet_oppcharge_truth_tt1_clone");
                THnSparse *hsparsejet_samecharge_truth_tt1_clone = (THnSparse *) hsparsejet_samecharge_truth_tt1->Clone("hsparsejet_samecharge_truth_tt1_clone");
                // THnSparse *hsparsejet_unweightedRL_truth_tt1_clone = (THnSparse *) hsparsejet_unweightedRL_truth_tt1->Clone("hsparsejet_unweightedRL_truth_tt1_clone");

                TH1D *h_pt_tt1_jetlevel_clone = (TH1D *) h_pt_tt1_jetlevel->Clone("h_pt_tt1_jetlevel_clone");
                // cout << "checkpoint 1" << endl;
                

                // get jet pT range - no D0 reconstruction, so don't make D0 cuts
                hsparsejet_deltap_truth_tt1_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
                hsparsejet_deltap_truth_tt1_clone->GetAxis(1)->SetRangeUser(RL_min, RL_max); // apply cut on RL
                hsparsejet_deltapt_truth_tt1_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); 
                hsparsejet_deltapt_truth_tt1_clone->GetAxis(1)->SetRangeUser(RL_min, RL_max);
                hsparsejet_oppcharge_truth_tt1_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); 
                hsparsejet_oppcharge_truth_tt1_clone->GetAxis(1)->SetRangeUser(RL_min, RL_max);
                hsparsejet_samecharge_truth_tt1_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); 
                hsparsejet_samecharge_truth_tt1_clone->GetAxis(1)->SetRangeUser(RL_min, RL_max);
                // hsparsejet_unweightedRL_truth_tt1_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); 
                // hsparsejet_unweightedRL_truth_tt1_clone->GetAxis(1)->SetRangeUser(RL_min, RL_max);

                h_pt_tt1_jetlevel_clone->GetXaxis()->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
                // cout << "checkpoint 2" << endl;
                

                // Project onto observable axis
                TH1D *hcorr_deltap_truth_tt1_proj = hsparsejet_deltap_truth_tt1_clone->Projection(2); //CALL THESE TH1*????
                TH1D *hcorr_deltapt_truth_tt1_proj = hsparsejet_deltapt_truth_tt1_clone->Projection(2);
                TH1D *hcorr_oppcharge_truth_tt1_proj = hsparsejet_oppcharge_truth_tt1_clone->Projection(1);
                TH1D *hcorr_samecharge_truth_tt1_proj = hsparsejet_samecharge_truth_tt1_clone->Projection(1);
                // TH1D *hcorr_unweightedRL_truth_tt1_proj = hsparsejet_unweightedRL_truth_tt1_clone->Projection(2);
                if (j==0) hdummyRL = hsparsejet_oppcharge_truth_tt1->Projection(1);
                if (j==0) hdummyRL2 = hsparsejet_oppcharge_truth_tt1->Projection(1);

                // cout << "checkpoint 3" << endl;

                // Set to appropriate name
                std::string hname = hcorr_deltap_truth_tt1_proj->GetName();
                hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_RL" + std::to_string(RL_min) + "-" + std::to_string(RL_max);
                hcorr_deltap_truth_tt1_proj->SetNameTitle(hname.c_str(), hname.c_str());

                hname = hcorr_deltapt_truth_tt1_proj->GetName();
                hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_RL" + std::to_string(RL_min) + "-" + std::to_string(RL_max);
                hcorr_deltapt_truth_tt1_proj->SetNameTitle(hname.c_str(), hname.c_str());

                hname = hcorr_oppcharge_truth_tt1_proj->GetName();
                hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_RL" + std::to_string(RL_min) + "-" + std::to_string(RL_max);
                hcorr_oppcharge_truth_tt1_proj->SetNameTitle(hname.c_str(), hname.c_str());

                hname = hcorr_samecharge_truth_tt1_proj->GetName();
                hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_RL" + std::to_string(RL_min) + "-" + std::to_string(RL_max);
                hcorr_samecharge_truth_tt1_proj->SetNameTitle(hname.c_str(), hname.c_str());

                // hname = hcorr_unweightedRL_truth_tt1_proj->GetName();
                // hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_RL" + std::to_string(RL_min) + "-" + std::to_string(RL_max);
                // hcorr_unweightedRL_truth_tt1_proj->SetNameTitle(hname.c_str(), hname.c_str());

                // cout << "checkpoint 4" << endl;

                // Rebin
                int n_obs_bins = 20; //50; //-1;
                //         double obs_bins[] = {1.00000000e-04, 1.25892541e-04, 1.58489319e-04, 1.99526231e-04,
                // 2.51188643e-04, 3.16227766e-04, 3.98107171e-04, 5.01187234e-04,
                // 6.30957344e-04, 7.94328235e-04, 1.00000000e-03, 1.25892541e-03,
                // 1.58489319e-03, 1.99526231e-03, 2.51188643e-03, 3.16227766e-03,
                // 3.98107171e-03, 5.01187234e-03, 6.30957344e-03, 7.94328235e-03,
                // 1.00000000e-02, 1.25892541e-02, 1.58489319e-02, 1.99526231e-02,
                // 2.51188643e-02, 3.16227766e-02, 3.98107171e-02, 5.01187234e-02,
                // 6.30957344e-02, 7.94328235e-02, 1.00000000e-01, 1.25892541e-01,
                // 1.58489319e-01, 1.99526231e-01, 2.51188643e-01, 3.16227766e-01,
                // 3.98107171e-01, 5.01187234e-01, 6.30957344e-01, 7.94328235e-01,
                // 1.00000000e+00};
                double obs_bins[] = {0.  , 0.05, 0.1 , 0.15, 0.2 , 0.25, 0.3 , 0.35, 0.4 , 0.45, 0.5 ,
                    0.55, 0.6 , 0.65, 0.7, 0.75, 0.8 , 0.85, 0.9 , 0.95, 1.  };


                //if rebinning

                TH1D* hcorr_deltap_truth_tt1 = (TH1D*) hcorr_deltap_truth_tt1_proj->Rebin(n_obs_bins, hcorr_deltap_truth_tt1_proj->GetName(), obs_bins);
                TH1D* hcorr_deltapt_truth_tt1 = (TH1D*) hcorr_deltapt_truth_tt1_proj->Rebin(n_obs_bins, hcorr_deltapt_truth_tt1_proj->GetName(), obs_bins); 
                
                // TH1D* hcorr_oppcharge_truth_tt1 = (TH1D*) hcorr_oppcharge_truth_tt1_proj->Rebin(n_obs_bins, hcorr_oppcharge_truth_tt1_proj->GetName(), obs_bins);
                // TH1D* hcorr_samecharge_truth_tt1 = (TH1D*) hcorr_samecharge_truth_tt1_proj->Rebin(n_obs_bins, hcorr_samecharge_truth_tt1_proj->GetName(), obs_bins);
                // TH1D* hcorr_unweightedRL_truth_tt1 = (TH1D*) hcorr_unweightedRL_truth_tt1_proj->Rebin(n_obs_bins, hcorr_unweightedRL_truth_tt1_proj->GetName(), obs_bins);                

                
                // TH1D* hcorr_deltap_truth_tt1 = (TH1D*) hcorr_deltap_truth_tt1_proj->Clone(hcorr_deltap_truth_tt1_proj->GetName()); //(hc_name + "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max)).c_str());
                // TH1D* hcorr_deltapt_truth_tt1 = (TH1D*) hcorr_deltapt_truth_tt1_proj->Clone(hcorr_deltapt_truth_tt1_proj->GetName());
                
                TH1D* hcorr_oppcharge_truth_tt1 = (TH1D*) hcorr_oppcharge_truth_tt1_proj->Clone(hcorr_oppcharge_truth_tt1_proj->GetName());
                TH1D* hcorr_samecharge_truth_tt1 = (TH1D*) hcorr_samecharge_truth_tt1_proj->Clone(hcorr_samecharge_truth_tt1_proj->GetName());
                // TH1D* hcorr_unweightedRL_truth_tt1 = (TH1D*) hcorr_unweightedRL_truth_tt1_proj->Clone(hcorr_unweightedRL_truth_tt1_proj->GetName());                
                

                // Find normalization factor
                double numjets_tt1 = h_pt_tt1_jetlevel_clone->Integral();

                // cout << "checkpoint 5" << endl;

                // Set normalization
                hcorr_deltap_truth_tt1->Scale(1/numjets_tt1, "width");
                hcorr_deltapt_truth_tt1->Scale(1/numjets_tt1, "width");
                hcorr_oppcharge_truth_tt1->Scale(1/numjets_tt1, "width");
                hcorr_samecharge_truth_tt1->Scale(1/numjets_tt1, "width");
                // hcorr_unweightedRL_truth_tt1->Scale(1/numjets_tt1, "width");

                // cout << "checkpoint 6" << endl;
                


                // // Find maximum
                // maxy = hi->GetMaximum() * 1.1;
                // hi->SetMaximum(maxy); 




                // make a canvas for each pt range
                TCanvas* c_deltap = new TCanvas();
                ProcessCanvas(c_deltap);
                // gPad->SetLogx(); // gPad->SetLogy();

                TCanvas* c_deltapt = new TCanvas();
                ProcessCanvas(c_deltapt);
                // c_deltapt->cd();
                // gPad->SetLogx(); // gPad->SetLogy();

                TCanvas* c_oppcharge = new TCanvas();
                ProcessCanvas(c_oppcharge);
                // c_oppcharge->cd();
                gPad->SetLogx(); // gPad->SetLogy();

                TCanvas* c_samecharge = new TCanvas();
                ProcessCanvas(c_samecharge);
                // c_samecharge->cd();
                gPad->SetLogx(); // gPad->SetLogy();

                TCanvas* c_chargeratio = new TCanvas();
                ProcessCanvas(c_chargeratio);
                // c_chargeratio->cd();
                gPad->SetLogx(); // gPad->SetLogy();

                TCanvas* c_unweightedRL = new TCanvas();
                ProcessCanvas(c_unweightedRL);
                // c_unweightedRL->cd();
                // gPad->SetLogx(); // gPad->SetLogy();



                
                // TLegend* l; // = new TLegend(0.17, 0.65, 0.5, 0.85);
                // l = new TLegend(0.1797168,0.400741,0.4562155,0.8885185,""); //(0.17, 0.4, 0.5, 0.53);
                // /*
                // l->SetTextSize(0.045);
                // // TLegend *leg = new TLegend(0.1797168,0.5390741,0.4562155,0.8885185,"");
                // l->AddEntry("NULL","PYTHIA 8 Monash 2013","h");
                // l->AddEntry("NULL","pp, #sqrt{#it{s}} = 13 TeV","h");
                // l->AddEntry("NULL","D^{0} #rightarrow K^{#minus} #pi^{+} and charge conj.","h");
                // l->AddEntry("NULL","in charged jets, anti-#it{k}_{T}, #it{R} = 0.4","h");
                // l->AddEntry("NULL",ptbin,"h");
                // l->AddEntry("NULL",ptD,"h");
                // */
                // l->SetTextSize(0.037);
                // l->SetBorderSize(0);
                // l->Draw("same");
                



                //Format color and style
                // int markercolor1 = kRed; //charm
                // int markerstyle1 = kFullCircle;
                // int markercolor2 = kViolet+2; //gluon
                // int markerstyle2 = 33;
                // int markercolor3 = kGreen+2; //light
                // int markerstyle3 = 21;
                // label1 = "R_{L} = " + std::to_string(RL_min) + "-" + std::to_string(RL_max);
                // label2 = "gluon-init jets";
                // TString label3 = "light-init jets";

                cout << "HOLA??  " << endl;

                // FormatHist(l, hcorr_deltap_truth_tt1, label[j], colors[j], markers[0]);
                // FormatHist(l2, hcorr_deltapt_truth_tt1, label[j], colors[j], markers[0]);
                // FormatHist(l2, hcorr_oppcharge_truth_tt1, label[j], colors[j], markers[0]);
                // FormatHist(l2, hcorr_samecharge_truth_tt1, label[j], colors[j], markers[2]);
                // FormatHist(l, hcorr_unweightedRL_truth_tt1, label[j], colors[j], markers[0]);

                cout << "hiiiiii??  " << endl;
                //print together
                hcorr_deltap_truth_tt1_arr.push_back((TH1D*) hcorr_deltap_truth_tt1->Clone(hcorr_deltap_truth_tt1->GetName()));
                hcorr_deltapt_truth_tt1_arr.push_back((TH1D*) hcorr_deltapt_truth_tt1->Clone(hcorr_deltapt_truth_tt1->GetName()));
                hcorr_oppcharge_truth_tt1_arr.push_back((TH1D*) hcorr_oppcharge_truth_tt1->Clone(hcorr_oppcharge_truth_tt1->GetName()));
                hcorr_samecharge_truth_tt1_arr.push_back((TH1D*) hcorr_samecharge_truth_tt1->Clone(hcorr_samecharge_truth_tt1->GetName()));
                
                label.push_back(Form("R_{L} = %.2f-%.2f", RL_min, RL_max));

                FormatHist(l, hcorr_deltap_truth_tt1_arr[j], label[j], colors[j], markers[0]);
                FormatHist(l2, hcorr_deltapt_truth_tt1_arr[j], label[j], colors[j], markers[0]);
                FormatHist(l2, hcorr_oppcharge_truth_tt1_arr[j], label[j], colors[j], markers[0]);
                FormatHist(l2, hcorr_samecharge_truth_tt1_arr[j], label[j], colors[j], markers[2]);
                // FormatHist(l, hcorr_unweightedRL_truth_tt1, label[j], colors[j], markers[0]);

                // Format histograms for plotting (this order needed to keep legend in order and graphs lookin good)
                // hc->GetXaxis()->SetTitle("#it{p}_{T}#it{R}_{L}");
                // hc->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
                // hg->GetXaxis()->SetTitle("#it{p}_{T}#it{R}_{L}");
                // hg->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
                // hl->GetXaxis()->SetTitle("#it{p}_{T}#it{R}_{L}");
                // hl->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
                cout << "hello??  " << endl;

                c_deltap->cd();
                hcorr_deltap_truth_tt1_arr[j]->Draw();

                c_deltapt->cd();
                hcorr_deltapt_truth_tt1_arr[j]->Draw();

                c_oppcharge->cd();
                hcorr_oppcharge_truth_tt1_arr[j]->Draw();

                c_samecharge->cd();
                hcorr_samecharge_truth_tt1_arr[j]->Draw();

                c_chargeratio->cd();
                std::string charge_ratio_name = "hcorr_chargeratio_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_RL" + std::to_string(RL_min) + "-" + std::to_string(RL_max);
                TH1D* hcorr_chargeratio_truth_tt1 = (TH1D*) hcorr_oppcharge_truth_tt1_arr[j]->Clone(charge_ratio_name.c_str());
                hcorr_chargeratio_truth_tt1->Divide(hcorr_samecharge_truth_tt1_arr[j]);
                FormatHist(l2, hcorr_chargeratio_truth_tt1, label[j], colors[j], markers[2]);
                hcorr_chargeratio_truth_tt1_arr.push_back((TH1D*) hcorr_chargeratio_truth_tt1->Clone(hcorr_chargeratio_truth_tt1->GetName()));
                hcorr_chargeratio_truth_tt1_arr[j]->Draw();
                // rp->GetLowYaxis()->SetNdivisions(505);
                // c->Update();

                // plot the combined graphs
                c_deltap_all->cd();
                // std::cout << "HERE!!!" << std::endl;
                hcorr_deltap_truth_tt1_arr[j]->Draw("same");
                if (j==0) {
                    hcorr_deltap_truth_tt1_arr[j]->SetMaximum(0.1); //65);
                }
                // c_deltap_all->Modified();
                // c_deltap_all->Update();
                l->Draw("same");

                c_deltapt_all->cd();
                hcorr_deltapt_truth_tt1_arr[j]->Draw("same");
                if (j==0) {
                    hcorr_deltapt_truth_tt1_arr[j]->SetMaximum(0.1); //65);
                }
                l->Draw("same");

                c_charge_all->cd();
                
                    // c_charge_all->cd(2);
                    // hcorr_chargeratio_truth_tt1_arr[j]->Draw("same");
                if (j==0) {
                    c_charge_all->Divide(1,2);                    
                }
                c_charge_all->cd(1);
                if (j==0) {
                    // gPad->SetLogx();
                    // FormatHist(l2, hdummyRL, "", colors[j], markers[0]);
                    hdummyRL->SetMarkerColorAlpha(kBlue, 0);
                    hdummyRL->SetLineColorAlpha(kRed, 0);
                    hdummyRL->SetMinimum(0.);
                    hdummyRL->SetMaximum(3.5);
                    hdummyRL->GetXaxis()->SetRangeUser(1e-4, 0.4); //(0.01, 0.4)
                    hdummyRL->Draw();
                }
                
                hcorr_oppcharge_truth_tt1_arr[j]->Draw("same");
                hcorr_samecharge_truth_tt1_arr[j]->Draw("same");

                hcorr_oppcharge_truth_tt1_arr[j]->SetMarkerColorAlpha(colors[j], 0.6);
                hcorr_oppcharge_truth_tt1_arr[j]->SetLineColorAlpha(colors[j], 0.6);
                hcorr_samecharge_truth_tt1_arr[j]->SetMarkerColorAlpha(colors[j], 0.6);
                hcorr_samecharge_truth_tt1_arr[j]->SetLineColorAlpha(colors[j], 0.6);
                l_right->AddEntry(hcorr_oppcharge_truth_tt1_arr[j],label[j]);
                // l_right->AddEntry(hcorr_samecharge_truth_tt1_arr[j],label[j]);
                
                l_right->Draw("same");
                // gPad->BuildLegend();


                c_charge_all->cd(2);
                if (j==0) {
                    // gPad->SetLogx();
                    hdummyRL2->SetMarkerColorAlpha(kBlue, 0);
                    hdummyRL2->SetLineColorAlpha(kRed, 0);
                    hdummyRL2->SetMinimum(0.);
                    hdummyRL2->SetMaximum(2.5);
                    hdummyRL2->GetXaxis()->SetRangeUser(1e-4, 0.4); //(0.01, 0.4)
                    hdummyRL2->Draw();
                    drawHoriLine(0.01, 0.4, 1., kBlack);
                }
                hcorr_chargeratio_truth_tt1_arr[j]->Draw("same");

                




                
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

                std::string RLname = "_RL" + std::to_string(RL_min) + "-" + std::to_string(RL_max);

                std::string fname_deltap = outdir + "corrhist_deltap_pt" + ptname + RLname + "_R" + jetR + add_name + ".pdf";
                std::string fname_deltapt = outdir + "corrhist_deltapt_pt" + ptname + RLname + "_R" + jetR + add_name + ".pdf";
                std::string fname_oppcharge = outdir + "corrhist_oppcharge_pt" + ptname + RLname + "_R" + jetR + add_name + ".pdf";
                std::string fname_samecharge = outdir + "corrhist_samecharge_pt" + ptname + RLname + "_R" + jetR + add_name + ".pdf";
                std::string fname_unweightedRL = outdir + "corrhist_unweightedRL_pt" + ptname + RLname + "_R" + jetR + add_name + ".pdf";
                std::string fname_chargeratio = outdir + "corrhist_chargeratio_pt" + ptname + RLname + "_R" + jetR + add_name + ".pdf";

                const char* fname_deltapc = fname_deltap.c_str();
                const char* fname_deltaptc = fname_deltapt.c_str();
                const char* fname_oppchargec = fname_oppcharge.c_str();
                const char* fname_samechargec = fname_samecharge.c_str();
                const char* fname_unweightedRLc = fname_unweightedRL.c_str();
                const char* fname_chargeratioc = fname_chargeratio.c_str();

                c_deltap->SaveAs(fname_deltapc);
                c_deltapt->SaveAs(fname_deltaptc);
                c_oppcharge->SaveAs(fname_oppchargec);
                c_samecharge->SaveAs(fname_samechargec);
                // c_unweightedRL->SaveAs(fname_unweightedRLc);
                c_chargeratio->SaveAs(fname_chargeratioc);

                delete c_deltap;
                delete c_deltapt;
                delete c_oppcharge;
                delete c_samecharge;
                // delete c_unweightedRL;
                delete c_chargeratio;

                f_out->cd();
                hcorr_deltap_truth_tt1->Write();
                hcorr_deltapt_truth_tt1->Write();
                hcorr_oppcharge_truth_tt1->Write();
                hcorr_samecharge_truth_tt1->Write();
                // hcorr_unweightedRL_truth_tt1->Write();
                hcorr_chargeratio_truth_tt1->Write();

                
                delete hcorr_deltap_truth_tt1;
                delete hcorr_deltapt_truth_tt1;
                delete hcorr_oppcharge_truth_tt1;
                delete hcorr_samecharge_truth_tt1;
                // delete hcorr_unweightedRL_truth_tt1;
                delete hcorr_chargeratio_truth_tt1;

                delete h_pt_tt1_jetlevel_clone;
                

            
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

            std::string fname_deltap_all = outdir + "corrhist_deltap_all_pt" + ptname + "_R" + jetR + add_name + ".pdf";
            std::string fname_deltapt_all = outdir + "corrhist_deltapt_all_pt" + ptname + "_R" + jetR + add_name + ".pdf";
            std::string fname_charge_all = outdir + "corrhist_charge_all_pt" + ptname + "_R" + jetR + add_name + ".pdf";
            
            const char* fname_deltapc_all = fname_deltap_all.c_str();
            const char* fname_deltaptc_all = fname_deltapt_all.c_str();
            const char* fname_chargec_all = fname_charge_all.c_str();
            
            c_deltap_all->SaveAs(fname_deltapc_all);
            c_deltapt_all->SaveAs(fname_deltaptc_all);
            c_charge_all->SaveAs(fname_chargec_all);

            delete c_deltap_all;
            delete c_deltapt_all;
            delete c_charge_all;

        } // pT bins loop
    } // jetR loop

    f->Close();
    delete f;


    return;
}
