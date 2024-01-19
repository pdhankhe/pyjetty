// ROOT macro to make quark-gluon jet plots
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

//checkExtra is how many point-to-point slopes after finding a decreasing slope I want to check
int findTopOfCurve(TH1* hist, int checkExtra=1) {
    
    //for each point, find the slope from the 
    int numbins = hist->GetNbinsX();
    int binstart = hist->FindBin(0.01);
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


void make_qg_plots_write() {

//    gROOT->SetBatch(); //prevents plots from showing up
    gStyle->SetOptStat(0);
    SetStyle();
    Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
    Double_t marker_size = 1.5;
    Double_t colors[16] = {kRed, kGreen+2, kBlue, kRed+1, kGreen+1, kBlue+1, kRed+2, kGreen+2, kBlue+2, kRed+3, kGreen+3, kBlue+3, kOrange+1, kViolet+1, kYellow+1, kCyan+1};


    // File containing quark vs gluon histograms

     //FOR WHEN WEIGHTED/UNWEIGHTED IN SAME FILE
    const char infile_D0_Preeti[] = "/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alihfjets/dev/hfjet/process/user/hf_EEC/D0jet_EEC_15_30_ForBeatrice.root";
    const char infile_D0[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/17651853/AnalysisResultsFinal.root"; //this is using thnsparse
    const char infile_Dstar[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/18008705/AnalysisResultsFinal.root"; //this is using thnsparse
    
    const char infile_D0_difNorm[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/17652140/AnalysisResultsFinal.root";
    const char infile_Dstar_difNorm[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/18008709/AnalysisResultsFinal.root";
    const char infile_D0wDstar_difNorm[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/18044499/AnalysisResultsFinal.root";
    const char infile_Dstar_difNorm_softpionremoved[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/18044503/AnalysisResultsFinal.root";
    const char infile_Dstar_difNorm_justsoftpion_noD0[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/18864762/AnalysisResultsFinal.root";
    const char infile_Dstar_difNorm_justsoftpionwD0[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/18864953/AnalysisResultsFinal.root";
    const char infile_Dstar_difNorm_justsoftpion[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/18865141/AnalysisResultsFinal.root";

    const char infile_noD0replacement_charmdecaysOFF[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/14680819/AnalysisResultsFinal.root"
    const char infile_noD0replacement_charmdecaysON[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/14680822/AnalysisResultsFinal.root"
    
    // int plot_case:
    // 0 = D0, D* decays turned off

    //CONTOL VARIABLES HERE
    int plot_case = 0;
    bool unnormalized = false; //default = false

    TFile* f;
    std::string logstring = "_log"; //not currently using this or have a bool

    TString label1 = "";

    f = new TFile(infile_D0, "READ");
    // cout << "output name will be " << add_name << endl;

    // Output directory
    std::string outdir = "plots/final/";//"plots/test/"; 
    // Output file for binned results
    std::string outfile = outdir + "AnalysisResultsFinal_D0_noDstardecays.root"; 
    TFile* f_out = new TFile(outfile.c_str(), "RECREATE");

    // Jet R value
    std::string jetR_list[] = { "0.4" };
    for (std::string jetR : jetR_list) {

        // Names of histograms in the file (quark, charm, gluon)
        const std::string hc_name = "h_EEC_JetPt_charm_R" + jetR;
        const std::string hc_unweighted_name = "h_EEC_JetPt_charm_R" + jetR + "_unweighted";
        const std::string hc_jet_name = "h_JetPt_charm_R" + jetR + "_jetlevel";

        const std::string hD0KpiNjets_name = "hD0KpiNjets"; //TODO later: = "hD0KpiNehD0KpiNjetsvents" for run 16729583



        const int pt_bins[] = { 10, 15, 30 }; //{ 10, 20, 40 }; // CHANGE HERE!!
        const int n_bins = 2; 
        for (int i = 0; i < n_bins; i++) {
            cout << "in pt bin" << i << endl;
            int pt_min = pt_bins[i];
            int pt_max = pt_bins[i+1];

            // define pt related variables
            TString ptbin = TString::Format("%d #leq #it{p}_{T}^{ch. jet} < %d GeV/#it{c}, #font[122]{|}#it{#eta}_{jet}#font[122]{|} #leq 0.5", pt_min, pt_max);
            TString ptD = TString::Format("5 #leq #it{p}_{T}^{D^{0}} < %d GeV/#it{c}, #font[122]{|}#it{y}_{D^{0}}#font[122]{|} #leq 0.8", pt_max);


            // make a canvas for each pt range
            TCanvas* c = new TCanvas();
            ProcessCanvas(c);
            c->cd();
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

            TH1* hD0KpiNjets = (TH1*) f->Get(hD0KpiNjets_name.c_str());
            

            //-------------------------------------------------//
            // find D0 reconstruction through charm
            THnSparse* hsparsejet_c = (THnSparse*) f->Get(hc_name.c_str());
            THnSparse* hsparsejet_c_jetlevel = (THnSparse*) f->Get(hc_jet_name.c_str());

            // testing - look at # jets before cuts
            cout << "numDtaggedjets from hist before cuts " << hsparsejet_c_jetlevel->Projection(0)->GetEntries() << endl;

            // for THnSparse: make clone to work with, make cuts, get projection
            THnSparse *hsparsejet_c_clone = (THnSparse *) hsparsejet_c->Clone("hsparsejet_c_clone");
            THnSparse *hsparsejet_c_jetlevel_clone = (THnSparse *) hsparsejet_c_jetlevel->Clone("hsparsejet_c_jetlevel_clone");

            // get jet pT range
            hsparsejet_c_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
            hsparsejet_c_clone->GetAxis(1)->SetRangeUser(5., pt_max); // apply cut on Dmeson pt
            hsparsejet_c_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
            hsparsejet_c_jetlevel_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
            hsparsejet_c_jetlevel_clone->GetAxis(1)->SetRangeUser(5., pt_max); // apply cut on Dmeson pt
            hsparsejet_c_jetlevel_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity

            // Project onto observable axis
            TH1D *hD0_proj = hsparsejet_c_clone->Projection(3); //CALL THESE TH1*????
            TH1D *hc1D_jet = hsparsejet_c_jetlevel_clone->Projection(0);

            // Set to appropriate name
            std::string hname = hD0_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hD0_proj->SetNameTitle(hname.c_str(), hname.c_str());

            // Find normalization factor
            double numjets_charm = hc1D_jet->Integral();
            double numDtaggedjets = hD0KpiNjets->GetEntries(); // this should match the # from projection before cuts
            double numDtaggedjets_fromhist = hc1D_jet->GetEntries();
            double numDtaggedjets_fromhist_eff = hc1D_jet->GetEffectiveEntries();

            cout << "numjets_charm " << numjets_charm << endl;
            cout << "numDtaggedjets " << numDtaggedjets << endl;
            cout << "numDtaggedjets_fromhist " << numDtaggedjets_fromhist << endl;
            cout << "numDtaggedjets_fromhist_eff " << numDtaggedjets_fromhist_eff << endl;

            // Set normalization
            // hD0_proj->Scale(1/numjets_charm, "width");


            // Rebin
            int n_obs_bins = 13; //-1;
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
            double obs_bins[14] = {0.0001, 0.00020892961308540387, 0.0004365158322401661, 0.0009120108393559096, 
                    0.0019054607179632482, 0.005754399373371567, 0.017378008287493762, 0.03630780547701014, 
                    0.07585775750291836, 0.15848931924611143, 0.3311311214825911, 0.47863009232263853, 0.6918309709189363, 1.0};

            cout << "About to rebin" << endl;
            TH1* hD0 = (TH1*) hD0_proj->Rebin(n_obs_bins, (hname + "rebin").c_str(), obs_bins);
            cout << "Rebin done" << endl;

            hD0->Scale(1/numjets_charm, "width");




            //Format color and style
            int markercolor1 = kGreen-5; //D0, me
            int markerstyle1 = kFullCircle;
            label1 = "D^{0}-tagged, c-init jets";

            

            // Format histograms for plotting (this order needed to keep legend in order and graphs lookin good)
            hD0->GetXaxis()->SetTitle("#it{R}_{L}");
            hD0->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
            cout << "about to format D0" << endl;
            FormatHist(l, hD0, label1, markercolor1, markerstyle1); //FormatHist(l, hD0, "D^{0}-tagged, c-init jets", kMagenta+3, kOpenSquare);
            l->AddEntry("NULL","          D* decays off","h");
            hD0->Draw("L same");
            
            
            double hD0_top_binpos = findTopOfCurve(hD0);
            drawVertLine(hD0->GetBinCenter(hD0_top_binpos), 0, hD0->GetBinContent(hD0_top_binpos), markercolor1, 1)->Draw();
            

            
            // vector<double> fullwidth_vec = findWidthOfCurve(hD0,  hD0_top_binpos);
            // drawHoriLine(fullwidth_vec[1], fullwidth_vec[2], fullwidth_vec[0], kMagenta+3, 1)->Draw();
            
            

            // draw legend
            l->Draw("same");



            // std::string fname = outdir + "QG_comp_pt" + std::to_string(pt_min) + '-' + std::to_string(pt_max) + "_R" + jetR + add_name; //"_charmdecaysONcomparison.pdf"; // + "_normbytype.pdf"; //"_nonorm.pdf";
            // const char* fnamec = fname.c_str();
            // c->SaveAs(fnamec);
            delete c;

            // if (pt_min == 10) { //} && grooming == "") {
                // Write rebinned histograms to root file
            f_out->cd();
            hD0->Write();
            // }
             
        } // pT bins loop
    } // jetR loop

    f->Close();
    delete f;

    return;
}
