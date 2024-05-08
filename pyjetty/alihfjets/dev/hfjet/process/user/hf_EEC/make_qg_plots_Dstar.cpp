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


void make_qg_plots_Dstar() {

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

    const char infile_D0_feeddown[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/21046102/AnalysisResultsFinal.root";
    
    // //test
    // const char infile_D0[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/17460087/AnalysisResultsFinal.root"; //this is using thnsparse
    // const char infile_Dstar[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/17460100/AnalysisResultsFinal.root"; //this is using thnsparse
    
    // const char infile_D0_difNorm[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/17460092/AnalysisResultsFinal.root";
    // const char infile_Dstar_difNorm[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/17460102/AnalysisResultsFinal.root";
    // const char infile_D0wDstar_difNorm[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/17460118/AnalysisResultsFinal.root";
    // const char infile_Dstar_difNorm_softpionremoved[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/17460125/AnalysisResultsFinal.root";

    bool include_gluon = true; //true = draw gluon, false = do not draw gluon
    bool unweighted = false; //true = draw unweighted, false = do not draw unweighted
    bool compareD0 = true; //true = compare the D0 reconstruction between me and Preeti - SET THIS TO TRUE IT IS THE PURPOSE OF THIS FILE
    
    // int plot_case:
    // 0 = compare Preeti's D0 w my D0 - DONE
    // 1 = compare D0 (norm with D0) with D* (norm with D*) - DONE
    // 2 = compare D0 (norm with D0) with D* (norm with (D0+D*)) - DONE
    // 3 = compare D0 (norm with D0) with D0+D* (norm with (D0+D*))
    // 4 = compare D* (norm with (D0+D*)) with D*, soft pion removed (norm with (D0+D*))
    // 5 = compare D* (norm with (D0+D*)) with just soft pion correlations (without corelating to D0) - DONE
    // 6 = compare D* (norm with (D0+D*)) with soft pion + D0 correlations - DONE
    // 7 = compare D* (norm with (D0+D*)) with just soft pion correlations
    // 8 = compare D0 (norm with D0) with the D0 feeddown

    //CONTOL VARIABLES HERE
    int plot_case = 3;
    bool unnormalized = false; //default = false

    TFile* f;
    TFile* f2;
    // std::vector<TFile*> files;
    std::string add_name;
    std::string quarkstring = include_gluon ? "" : "_justquarks";
    std::string unweightedstring = unweighted ? "_samefilecomp" : "";
    std::string logstring = "_log"; //not currently using this or have a bool
    // TString ptbin = compareD0 ? "15 #leq #it{p}_{T}^{ch. jet} < 30 GeV/#it{c}, #font[122]{|}#it{#eta}_{jet}#font[122]{|} #leq 0.5" : "15 #leq #it{p}_{T}^{ch. jet} < 30 GeV/#it{c}";
    // TString ptD = "5 #leq #it{p}_{T}^{D^{0}} < 30 GeV/#it{c}, #font[122]{|}#it{y}_{D^{0}}#font[122]{|} #leq 0.8";

    TString label1 = "";
    TString label2 = "";  
    // if (compareD0) { 
    // f = new TFile(infile_D0, "READ");
    // f2 = new TFile(infile_Dstar, "READ");
    // add_name = "_Dstar_17412030.pdf";
    if (plot_case == 0) {
        f = new TFile(infile_D0, "READ");
        f2 = new TFile(infile_D0_Preeti, "READ");
        add_name = "_Dstar_plot_case0.pdf";
    } else if (plot_case == 1) {
        f = new TFile(infile_D0, "READ");
        f2 = new TFile(infile_Dstar, "READ");
        if (unnormalized == true) {
            add_name = "_Dstar_plot_case1_unnormalized.pdf";
        } else {
            add_name = "_Dstar_plot_case1.pdf";
        }
    } else if (plot_case == 2) {
        f = new TFile(infile_D0, "READ");
        f2 = new TFile(infile_Dstar_difNorm, "READ");
        if (unnormalized == true) {
            add_name = "_Dstar_plot_case2_unnormalized.pdf";
        } else {
            add_name = "_Dstar_plot_case2.pdf";
        }
    } else if (plot_case == 3) {
        f = new TFile(infile_D0, "READ");
        f2 = new TFile(infile_D0wDstar_difNorm, "READ");
        add_name = "_Dstar_plot_case3.pdf";
    } else if (plot_case == 4) {
        f = new TFile(infile_Dstar_difNorm, "READ");
        f2 = new TFile(infile_Dstar_difNorm_softpionremoved, "READ");
        add_name = "_Dstar_plot_case4.pdf";
    } else if (plot_case == 5) {
        f = new TFile(infile_Dstar_difNorm, "READ");
        f2 = new TFile(infile_Dstar_difNorm_justsoftpion_noD0, "READ");
        add_name = "_Dstar_plot_case5.pdf";
    } else if (plot_case == 6) {
        f = new TFile(infile_Dstar_difNorm, "READ");
        f2 = new TFile(infile_Dstar_difNorm_justsoftpionwD0, "READ");
        add_name = "_Dstar_plot_case6.pdf";
    } else if (plot_case == 7) {
        f = new TFile(infile_Dstar_difNorm, "READ");
        f2 = new TFile(infile_Dstar_difNorm_justsoftpion, "READ");
        add_name = "_Dstar_plot_case7.pdf";
    } else if (plot_case == 8) {
        f = new TFile(infile_D0, "READ");
        f2 = new TFile(infile_D0_feeddown, "READ");
        add_name = "_Dstar_plot_case8.pdf";
    }
    cout << "output name will be " << add_name << endl;

    // Output directory
    //std::string outdir = "/rstorage/alice/AnalysisResults/ang/1224559/plots/";
    std::string outdir = "plots/final/15-30/";//"plots/test/";
    std::string outdir_root = "plots/final/";
    // Output file for binned results
    std::string outfile = outdir_root + "AnalysisResultsFinal_afteranalysis_Dstar_plotcase" + std::to_string(plot_case) + ".root";
    //std::string outfile = outdir + "AnalysisResultsFinal_test.root";
    TFile* f_out = new TFile(outfile.c_str(), "RECREATE");

    // Angularity alpha value
//    std::string alpha_list[] = { "1", "1.5", "2", "3" };
//    std::string grooming_list[] = { "", "_SD_zcut02_B0" };
    std::string jetR_list[] = { "0.4" };
    for (std::string jetR : jetR_list) {
        //        for (std::string alpha : alpha_list) {
        //            for (std::string grooming : grooming_list) {

        // Names of histograms in the file (quark, charm, gluon)
        const std::string hi_name = "h_EEC_JetPt_inclusive_R" + jetR;
        const std::string hi_unweighted_name = "h_EEC_JetPt_inclusive_R" + jetR + "_unweighted";
        const std::string hi_jet_name = "h_JetPt_inclusive_R" + jetR + "_jetlevel";
        const std::string hc_name = "h_EEC_JetPt_charm_R" + jetR;
        const std::string hl_name = "h_EEC_JetPt_light_R" + jetR;
        const std::string hg_name = "h_EEC_JetPt_gluon_R" + jetR;
        const std::string hc_unweighted_name = "h_EEC_JetPt_charm_R" + jetR + "_unweighted";
        const std::string hl_unweighted_name = "h_EEC_JetPt_light_R" + jetR + "_unweighted";
        const std::string hg_unweighted_name = "h_EEC_JetPt_gluon_R" + jetR + "_unweighted";
        const std::string hc_jet_name = "h_JetPt_charm_R" + jetR + "_jetlevel";
        const std::string hl_jet_name = "h_JetPt_light_R" + jetR + "_jetlevel";
        const std::string hg_jet_name = "h_JetPt_gluon_R" + jetR + "_jetlevel";
        const std::string hb_name = "h_EEC_JetPt_beauty_R" + jetR;
        const std::string hb_jet_name = "h_JetPt_beauty_R" + jetR + "_jetlevel";
        

        const std::string hD0KpiNjets_name = "hD0KpiNjets"; //TODO later: = "hD0KpiNehD0KpiNjetsvents" for run 16729583

        const std::string D0_jet_name = "D0jet_EEC_15_30";

        // find D0 reconstruction through charm
        THnSparse* hsparsejet_c = (THnSparse*) f->Get(hc_name.c_str());
        THnSparse* hsparsejet_c_jetlevel = (THnSparse*) f->Get(hc_jet_name.c_str());
        // TH1* hc1D_jet = (TH1*) f->Get(hc_jet_name.c_str()); 

        // testing - look at # jets before cuts
        cout << "numDtaggedjets from hist before cuts " << hsparsejet_c_jetlevel->Projection(0)->GetEntries() << endl;

        
        // save D0 pt spectrum
        TH1D *hD0_pT = hsparsejet_c_jetlevel->Projection(1); //use jet level
        hD0_pT->SetNameTitle("hD0_pt", "hD0_pt");




        //const int pt_bins[] = { 10, 20, 40, 60, 80, 100, 150 };
        const int pt_bins[] = { 7, 10, 15, 30 }; //{ 10, 20, 40 };
        const int d0_pt_cuts[] = { 3, 5, 5 };
        const int n_bins = 3; //2;
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

            double maxy = 0;

            // start files loop
            // for (int ifile=0; ifile<files.size(); ifile++) {


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

            // THnSparse* hsparsejet_i = (THnSparse*) f->Get(hi_name.c_str());
            // THnSparse* hsparsejet_i_jetlevel = (THnSparse*) f->Get(hi_jet_name.c_str());
            // // TH1* hi1D_jet = (TH1*) f->Get(hi_jet_name.c_str()); 
            TH1* hD0KpiNjets = (TH1*) f->Get(hD0KpiNjets_name.c_str());
            // TH1* hD0 = (TH1*) f2->Get(D0_jet_name.c_str());
            
    
            // // for THnSparse: make clone to work with, make cuts, get projection
            // THnSparse *hsparsejet_i_clone = (THnSparse *) hsparsejet_i->Clone("hsparsejet_i_clone");
            // THnSparse *hsparsejet_i_jetlevel_clone = (THnSparse *) hsparsejet_i_jetlevel->Clone("hsparsejet_i_jetlevel_clone");
            
            

            // get jet pT range
            // hsparsejet_i_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
            // hsparsejet_i_clone->GetAxis(1)->SetRangeUser(5., pt_max); // apply cut on Dmeson pt
            // hsparsejet_i_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
            // hsparsejet_i_jetlevel_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
            // hsparsejet_i_jetlevel_clone->GetAxis(1)->SetRangeUser(5., pt_max); // apply cut on Dmeson pt
            // hsparsejet_i_jetlevel_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
            // hi1D_jet->GetXaxis()->SetRangeUser(pt_min, pt_max);

            // // Project onto observable axis
            // // TH1* hi = (TH1*) hi2D->ProjectionY();
            // TH1D *hi = hsparsejet_i_clone->Projection(3); //CALL THESE TH1*????
            // TH1D *hi1D_jet = hsparsejet_i_jetlevel_clone->Projection(0); //project onto jet pt axis??

            // // Set to appropriate name
            // std::string hname = hi->GetName();
            // hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            // hi->SetNameTitle(hname.c_str(), hname.c_str());

            //-------------------------------------------------//

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
            // hc2D->GetXaxis()->SetRangeUser(pt_min, pt_max);
            // hc1D_jet->GetXaxis()->SetRangeUser(pt_min, pt_max);

            // Project onto observable axis
            TH1D *hD0 = hsparsejet_c_clone->Projection(3); //CALL THESE TH1*????
            TH1D *hc1D_jet = hsparsejet_c_jetlevel_clone->Projection(0);
            // TH1* hc = (TH1*) hc2D->ProjectionY();

            // Set to appropriate name
            std::string hname = hD0->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hD0->SetNameTitle(hname.c_str(), hname.c_str());

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
            if (unnormalized == false) {
                hD0->Scale(1/numjets_charm, "width");
            }
            // hc->Scale(numDtaggedjets, "width"); // this is wrong - maybe bc EEC is scaled when saved to root file

            
            TH1* hDstar;
            if (plot_case == 0) {
                hDstar = (TH1*) f2->Get(D0_jet_name.c_str());
            } else {
                // Find D* histogram
                THnSparse* hsparsejet_c_Dstar = (THnSparse*) f2->Get(hc_name.c_str());
                THnSparse* hsparsejet_c_Dstar_jetlevel = (THnSparse*) f2->Get(hc_jet_name.c_str());
                // TH1* hc1D_jet = (TH1*) f->Get(hc_jet_name.c_str()); 

                // for THnSparse: make clone to work with, make cuts, get projection
                THnSparse *hsparsejet_c_Dstar_clone = (THnSparse *) hsparsejet_c_Dstar->Clone("hsparsejet_c_Dstar_clone");
                THnSparse *hsparsejet_c_Dstar_jetlevel_clone = (THnSparse *) hsparsejet_c_Dstar_jetlevel->Clone("hsparsejet_c_Dstar_jetlevel_clone");

                // get jet pT range
                hsparsejet_c_Dstar_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
                hsparsejet_c_Dstar_clone->GetAxis(1)->SetRangeUser(d0_pt_cuts[i], pt_max); // apply cut on Dmeson pt
                hsparsejet_c_Dstar_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
                hsparsejet_c_Dstar_jetlevel_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
                hsparsejet_c_Dstar_jetlevel_clone->GetAxis(1)->SetRangeUser(d0_pt_cuts[i], pt_max); // apply cut on Dmeson pt
                hsparsejet_c_Dstar_jetlevel_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
                // hc2D->GetXaxis()->SetRangeUser(pt_min, pt_max);
                // hc1D_jet->GetXaxis()->SetRangeUser(pt_min, pt_max);

                // Project onto observable axis
                hDstar = hsparsejet_c_Dstar_clone->Projection(3); //CALL THESE TH1*????
                TH1D *hc1D_Dstar_jet = hsparsejet_c_Dstar_jetlevel_clone->Projection(0);
                // TH1* hc = (TH1*) hc2D->ProjectionY();

                // Set to appropriate name
                hname = hDstar->GetName();
                hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
                hDstar->SetNameTitle(hname.c_str(), hname.c_str());

                // Find normalization factor
                double numjets_charm_Dstar = hc1D_Dstar_jet->Integral();
                // double numDtaggedjets = hD0KpiNjets->GetEntries();
                double numDtaggedjets_fromhist_charm = hc1D_Dstar_jet->GetEntries();

                cout << "numDtaggedjets_fromhist_Dstar " << numDtaggedjets_fromhist_charm << endl;

                // Set normalization
                if (unnormalized == false) {
                    hDstar->Scale(1/numjets_charm_Dstar, "width");
                }

            }



            // Rebin
            int n_obs_bins = 50; //-1;
//            double obs_bins[51]; //[60];
            // double obs_bins[] = {1.00000000e-04, 1.20226443e-04, 1.44543977e-04, 1.73780083e-04, 2.08929613e-04, 2.51188643e-04, 3.01995172e-04, 3.63078055e-04, 4.36515832e-04, 5.24807460e-04, 6.30957344e-04, 7.58577575e-04, 9.12010839e-04, 1.09647820e-03, 1.31825674e-03, 1.58489319e-03, 1.90546072e-03, 2.29086765e-03, 2.75422870e-03, 3.31131121e-03, 3.98107171e-03, 4.78630092e-03, 5.75439937e-03, 6.91830971e-03, 8.31763771e-03, 1.00000000e-02, 1.20226443e-02, 1.44543977e-02, 1.73780083e-02, 2.08929613e-02, 2.51188643e-02, 3.01995172e-02, 3.63078055e-02, 4.36515832e-02, 5.24807460e-02, 6.30957344e-02, 7.58577575e-02, 9.12010839e-02, 1.09647820e-01, 1.31825674e-01, 1.58489319e-01, 1.90546072e-01, 2.29086765e-01, 2.75422870e-01, 3.31131121e-01, 3.98107171e-01, 4.78630092e-01, 5.75439937e-01, 6.91830971e-01, 8.31763771e-01, 1.00000000e+00};
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

            // // Find normalization factor
            // double numjets_inclusive = hi1D_jet->Integral();

            // // Set normalization
            // hi->Scale(1/numjets_inclusive, "width");

            // // Find maximum
            // maxy = hi->GetMaximum() * 1.1;
            // hi->SetMaximum(maxy); 


            // print # of events and jets
            TH1I *hD0Nevents1 = (TH1I*) f->Get("hD0Nevents");
            TH1I *hD0KpiNevents1 = (TH1I*) f->Get("hD0KpiNevents");
            TH1I *hD0KpiNjets1 = (TH1I*) f->Get("hD0KpiNjets");
            TH1I *hDstarNjets1 = (TH1I*) f->Get("hDstarNjets");

            cout << "------------------------------- " << endl;
            cout << "D0 Nevents 1 " << hD0Nevents1->GetEntries() << endl;
            cout << "D0->Kpi Nevents 1 " << hD0KpiNevents1->GetEntries() << endl;
            cout << "D0->Kpi (primordial) Njets 1 " << hD0KpiNjets1->GetEntries() << endl;
            cout << "D* Njets 1 " << hDstarNjets1->GetEntries() << endl;
            cout << "------------------------------- " << endl;
            if (plot_case != 0 ) {
                TH1I *hD0Nevents2 = (TH1I*) f2->Get("hD0Nevents");
                TH1I *hD0KpiNevents2 = (TH1I*) f2->Get("hD0KpiNevents");
                TH1I *hD0KpiNjets2 = (TH1I*) f2->Get("hD0KpiNjets");
                TH1I *hDstarNjets2 = (TH1I*) f2->Get("hDstarNjets");

                cout << "D0 Nevents 2 " << hD0Nevents2->GetEntries() << endl;
                cout << "D0->Kpi Nevents 2 " << hD0KpiNevents2->GetEntries() << endl;
                cout << "D0->Kpi (primordial) Njets 2 " << hD0KpiNjets2->GetEntries() << endl;
                cout << "D* Njets 2 " << hDstarNjets2->GetEntries() << endl;
                cout << "------------------------------- " << endl;
            }



            //Format color and style
            int markercolor1 = kGreen-5; //D0, me
            int markerstyle1 = kFullCircle;
            int markercolor2 = 0;
            int markerstyle2 = 0;
            if (plot_case == 0) { 
                markercolor2 = kMagenta+3; //D0, Preeti
                markerstyle2 = kOpenSquare;
                label1 = "D^{0}-tagged, c-init jets, Beatrice";
                label2 = "D^{0}-tagged, c-init jets, Preeti";
            } else if (plot_case == 1 or plot_case == 2) {
                markercolor2 = kRed-3; //D*
                markerstyle2 = 29;
                label1 = "D^{0}-tagged, c-init jets";
                label2 = "D*-tagged, c-init jets";
            } else if (plot_case == 3){
                markercolor2 = kRed-3; //D*+D0
                markerstyle2 = 29;
                label1 = "D^{0}-tagged, c-init jets";
                label2 = "(D^{0}+D*)-tagged, c-init jets";
            } else if (plot_case >= 4){
                markercolor1 = kRed-3; //D*+D0
                markerstyle1 = 29;
                markercolor2 = kViolet+2;
                markerstyle2 = 33;
                label1 = "D*-tagged, c-init jets";
                if (plot_case == 4) {
                    label2 = "D*-tagged, soft #pi removed, c-init jets";
                } else if (plot_case == 5) {
                    label2 = "D*-tagged, only soft #pi";
                } else if (plot_case == 6) {
                    label2 = "D*-tagged, soft #pi w/ D0, c-init jets";
                } else if (plot_case == 7) {
                    label2 = "D*-tagged, only soft #pi, c-init jets";
                }
            } 
            

            // Format histograms for plotting (this order needed to keep legend in order and graphs lookin good)
            hDstar->GetXaxis()->SetTitle("#it{R}_{L}");
            hDstar->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
            hD0->GetXaxis()->SetTitle("#it{R}_{L}");
            hD0->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
            cout << "about to format D0" << endl;
            FormatHist(l, hD0, label1, markercolor1, markerstyle1); //FormatHist(l, hD0, "D^{0}-tagged, c-init jets", kMagenta+3, kOpenSquare);
            l->AddEntry("NULL","          D* decays off","h");
            cout << "about to format Dstar" << endl;
            FormatHist(l, hDstar, label2, markercolor2, markerstyle2); //"D*-tagged, c-init jets", kRed-7, 29);
            if (plot_case == 0) l->AddEntry("NULL","          D* decays off","h");
            if (plot_case == 5) l->AddEntry("NULL","           (w/out D^{0} correl), c-init jets","h");
            if ((plot_case == 1 && unnormalized == false && pt_min == 15) or plot_case == 7 or plot_case == 0 or plot_case == 4) { //} or plot_case == 6) {
                hDstar->Draw("L same");
                hD0->Draw("L same");
            } else {
                hD0->Draw("L same");
                hDstar->Draw("L same");
            }
            
            
            // hDstar->Draw("L same");
            

            // double hc_top_binpos = findTopOfCurve(hc);
            // drawVertLine(hc->GetBinCenter(hc_top_binpos), 0, hc->GetBinContent(hc_top_binpos), kRed-7, 1)->Draw();
            double hD0_top_binpos = findTopOfCurve(hD0);
            drawVertLine(hD0->GetBinCenter(hD0_top_binpos), 0, hD0->GetBinContent(hD0_top_binpos), markercolor1, 1)->Draw();
            double hDstar_top_binpos = findTopOfCurve(hDstar);
            drawVertLine(hDstar->GetBinCenter(hDstar_top_binpos), 0, hDstar->GetBinContent(hDstar_top_binpos), markercolor2, 1)->Draw();


            
            // vector<double> fullwidth_vec = findWidthOfCurve(hD0,  hD0_top_binpos);
            // drawHoriLine(fullwidth_vec[1], fullwidth_vec[2], fullwidth_vec[0], kMagenta+3, 1)->Draw();
            
            
            // Add legend about D0 info
            
            // TLegend *leg = new TLegend(0.1797168,0.5390741,0.4562155,0.8885185,"");
            // leg->AddEntry("NULL","PYTHIA 8","h");
            // leg->AddEntry("NULL","pp, #sqrt{#it{s}} = 13 TeV","h");
            // leg->AddEntry("NULL","D^{0} #rightarrow K^{#minus} #pi^{+} and charge conj.","h");
            // leg->AddEntry("NULL","in charged jets, anti-#it{k}_{T}, #it{R} = 0.4","h");
            // leg->AddEntry("NULL",ptbin,"h");
            // leg->AddEntry("NULL",ptD,"h");
            // leg->SetTextSize(0.037);
            // leg->SetBorderSize(0);
            // leg->Draw("same");
            // l->AddEntry("NULL","          Charm decay off","h");
            




            // make ratio plot
            // auto rp = new TRatioPlot(hD0, hc);
            // rp->Draw();
            // rp->GetLowYaxis()->SetNdivisions(505);
            // c->Update();


    


            



            // } //end of file loop?

            // draw legend
            l->Draw("same");


            // Label plots
            TText* t = new TText(0.0003, maxy/1.05, "PYTHIA 8 Monash 2013");
            t->SetTextAlign(22);
            t->SetTextFont(43);
            t->SetTextSize(20);
            // t->Draw("same");
            //std::string text = "#it{#alpha} = " + alpha;
            //t->SetText(0.3, 10, text.c_str());

            std::string fname = outdir + "QG_comp_pt" + std::to_string(pt_min) + '-' + std::to_string(pt_max) + "_R" + jetR + add_name; //"_charmdecaysONcomparison.pdf"; // + "_normbytype.pdf"; //"_nonorm.pdf";
            const char* fnamec = fname.c_str();
            c->SaveAs(fnamec);
            delete c;
            delete t;

            // if (pt_min == 10) { //} && grooming == "") {
                // Write rebinned histograms to root file
                f_out->cd();
                // if (!inclusive) {
                hD0->Write();
                hDstar->Write();
                // }
            // }
             
        } // pT bins loop

        hD0_pT->Write();
//            } // grooming loop
//        } // alpha loop
    } // jetR loop

    f->Close();
    delete f;
    if (compareD0) {
        f2->Close();
        delete f2;
    }


    return;
}
