// ROOT macro to make quark-gluon jet plots
// Ezra Lesser (elesser@berkeley.edu)

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

TLine * drawVertLine(double x1, double y1, double y2, int color, int linestyle=2){
    auto fvertline = new TLine(x1, y1, x1, y2);
	fvertline->SetLineWidth(1);
    fvertline->SetLineColor(color);
    fvertline->SetLineStyle(linestyle);
    return fvertline;

}


void make_qg_plots_samefilecomp() {

//    gROOT->SetBatch(); //prevents plots from showing up
    gStyle->SetOptStat(0);
    SetStyle();
    Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
    Double_t marker_size = 1.5;
    Double_t colors[16] = {kRed, kGreen+2, kBlue, kRed+1, kGreen+1, kBlue+1, kRed+2, kGreen+2, kBlue+2, kRed+3, kGreen+3, kBlue+3, kOrange+1, kViolet+1, kYellow+1, kCyan+1};


    // File containing quark vs gluon histograms

     //FOR WHEN WEIGHTED/UNWEIGHTED IN SAME FILE
    const char infile_charmOFF[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/14680819/AnalysisResultsFinal.root"; //perlmutter
    const char infile_charmON[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/14680822/AnalysisResultsFinal.root"; //perlmutter 
    const char infile_replaceKPON[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/13777236/AnalysisResultsFinal.root";
    const char infile_D0_fromPreeti[] = "/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alihfjets/dev/hfjet/process/user/hf_EEC/D0jet_EEC_15_30_ForBeatrice.root";
    const char infile_D0[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/15442613/AnalysisResultsFinal.root";
    
    bool include_gluon = true; //true = draw gluon, false = do not draw gluon
    bool unweighted = false; //true = draw unweighted, false = do not draw unweighted
    bool inclusive = true; //true = access and draw only inclusive, false = do not draw
    bool selffoundD0 = true; // true = use the D0 reconstruction I did, false = use Preeti's
    bool switchON = false; //true = charm ON, false = charm OFF
    bool replaceKPON = false; //true = D0 reconstruction implemented
    TFile* f;
    TFile* f2;
    // std::vector<TFile*> files;
    std::string add_name;
    std::string quarkstring = include_gluon ? "" : "_justquarks";
    std::string unweightedstring = unweighted ? "_samefilecomp" : "";
    std::string logstring = "_log"; //not currently using this or have a bool
    TString ptbin = inclusive ? "15 #leq #it{p}_{T}^{ch. jet} < 30 GeV/#it{c}, #font[122]{|}#it{#eta}_{jet}#font[122]{|} #leq 0.5" : "15 #leq #it{p}_{T}^{ch. jet} < 30 GeV/#it{c}";
    TString ptD = "5 #leq #it{p}_{T}^{D^{0}} < 30 GeV/#it{c}, #font[122]{|}#it{y}_{D^{0}}#font[122]{|} #leq 0.8";
                
    if (inclusive) { 
        f = new TFile(infile_charmOFF, "READ");
        f2 = new TFile(infile_D0, "READ");
        add_name = "_inclusive.pdf";
    } else if (replaceKPON) {
        f = new TFile(infile_replaceKPON, "READ");
        // f2 = new TFile(infile_charmON, "READ");
        // f3 = new Tfile(infile_charmOFF, "READ");
        // files.push_back(f); files.push_back(f2); files.push_back(f3);

        // bool compall = true;
        // std::string compallstring = "_compall";
        add_name = "_replaceKP" + quarkstring + unweightedstring + ".pdf"; //+ compallstring + ".pdf"; 

    } else {
        if (switchON) {
            f = new TFile(infile_charmON, "READ");
            add_name = "_charmdecaysON" + quarkstring + unweightedstring + ".pdf"; //_samefilecomp_log.pdf"; 
        } else {
            f = new TFile(infile_charmOFF, "READ");
            add_name = "_charmdecaysOFF" + quarkstring + unweightedstring + ".pdf"; //_samefilecomp_log.pdf";
        }
        // files.push_back(f);
    }
    cout << "output name will be " << add_name << endl;

    // Output directory
    //std::string outdir = "/rstorage/alice/AnalysisResults/ang/1224559/plots/";
    std::string outdir = "plots/";
    // Output file for binned results
    std::string outfile = outdir + "AnalysisResultsFinal_afteranalysis.root";
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

        const std::string D0_jet_name = "D0jet_EEC_15_30";


        //const int pt_bins[] = { 10, 20, 40, 60, 80, 100, 150 };
        const int pt_bins[] = { 15, 30 }; //{ 10, 20, 40 };
        const int n_bins = 1; //2;
        for (int i = 0; i < n_bins; i++) {
            cout << "in pt bin" << i << endl;
            int pt_min = pt_bins[i];
            int pt_max = pt_bins[i+1];

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
            if (inclusive) {

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

                TH2* hi2D = (TH2*) f->Get(hi_name.c_str());
                TH1* hi1D_jet = (TH1*) f->Get(hi_jet_name.c_str()); 
                TH1* hD0;
                if (!selffoundD0) {
                    hD0 = (TH1*) f2->Get(D0_jet_name.c_str());
                }

                // get jet pT range
                hi2D->GetXaxis()->SetRangeUser(pt_min, pt_max);
                hi1D_jet->GetXaxis()->SetRangeUser(pt_min, pt_max);

                // Project onto observable axis
                TH1* hi = (TH1*) hi2D->ProjectionY();

                // Set to appropriate name
                std::string hname = hi->GetName();
                hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
                hi->SetNameTitle(hname.c_str(), hname.c_str());

                // find D0 reconstruction through charm
                if (selffoundD0) {
                    TH2* hc2D = (TH2*) f2->Get(hc_name.c_str());
                    TH1* hc1D_jet = (TH1*) f2->Get(hc_jet_name.c_str()); 

                    // get jet pT range
                    hc2D->GetXaxis()->SetRangeUser(pt_min, pt_max);
                    hc1D_jet->GetXaxis()->SetRangeUser(pt_min, pt_max);

                    // Project onto observable axis
                    hD0 = (TH1*) hc2D->ProjectionY();

                    // Set to appropriate name
                    std::string hname = hD0->GetName();
                    hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
                    hD0->SetNameTitle(hname.c_str(), hname.c_str());

                    // Find normalization factor
                    double numjets_charm = hc1D_jet->Integral();

                    // Set normalization
                    hD0->Scale(1/numjets_charm, "width");

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

                // Find normalization factor
                double numjets_inclusive = hi1D_jet->Integral();

                // Set normalization
                hi->Scale(1/numjets_inclusive, "width");

                // Find maximum
                maxy = hi->GetMaximum() * 1.1;
                hi->SetMaximum(maxy); 

                // Format histograms for plotting (this order needed to keep legend in order and graphs lookin good)
                hi->GetXaxis()->SetTitle("#it{R}_{L}");
                hi->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
                FormatHist(l, hD0, "D^{0}-tagged jets", kRed+1, kFullSquare);
                FormatHist(l, hi, "Inclusive jets", kViolet+1, 29);
                hi->Draw("L same");
                hD0->Draw("L same");

                

                double hi_top_binpos = findTopOfCurve(hi);
                drawVertLine(hi->GetBinCenter(hi_top_binpos), 0, hi->GetBinContent(hi_top_binpos), kViolet+1, 1)->Draw();
                double hD0_top_binpos = findTopOfCurve(hD0);
                drawVertLine(hD0->GetBinCenter(hD0_top_binpos), 0, hD0->GetBinContent(hD0_top_binpos), kRed+1, 1)->Draw();
                
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
                l->AddEntry("NULL","          Charm decay off","h");
                

            } else {

                l = new TLegend(0.1797168,0.4040741,0.4562155,0.8885185,""); //(0.17, 0.35, 0.5, 0.65);
                // TLegend *leg = new TLegend(0.1797168,0.6640741,0.4562155,0.8885185,"");
                l->AddEntry("NULL","PYTHIA 8 Monash 2013","h");
                l->AddEntry("NULL","pp, #sqrt{#it{s}} = 13 TeV","h");
                l->AddEntry("NULL","in charged jets, anti-#it{k}_{T}, #it{R} = 0.4","h");
                l->AddEntry("NULL",ptbin,"h");
                l->SetTextSize(0.037);
                l->SetBorderSize(0);

                TH2* hc2D = (TH2*) f->Get(hc_name.c_str());
                TH2* hl2D = (TH2*) f->Get(hl_name.c_str());
                TH2* hg2D = (TH2*) f->Get(hg_name.c_str());
                TH2* hc2D_unweighted = (TH2*) f->Get(hc_unweighted_name.c_str());//->Get(hc_unweighted_name.c_str());
                TH2* hl2D_unweighted = (TH2*) f->Get(hl_unweighted_name.c_str());//->Get(hl_unweighted_name.c_str());
                TH2* hg2D_unweighted = (TH2*) f->Get(hg_unweighted_name.c_str());//->Get(hg_unweighted_name.c_str());
                TH1* hc1D_jet = (TH1*) f->Get(hc_jet_name.c_str());
                TH1* hl1D_jet = (TH1*) f->Get(hl_jet_name.c_str());
                TH1* hg1D_jet = (TH1*) f->Get(hg_jet_name.c_str());

                // Set jet pT range
                hc2D->GetXaxis()->SetRangeUser(pt_min, pt_max);
                hl2D->GetXaxis()->SetRangeUser(pt_min, pt_max);
                hg2D->GetXaxis()->SetRangeUser(pt_min, pt_max);
                hc2D_unweighted->GetXaxis()->SetRangeUser(pt_min, pt_max);
                hl2D_unweighted->GetXaxis()->SetRangeUser(pt_min, pt_max);
                hg2D_unweighted->GetXaxis()->SetRangeUser(pt_min, pt_max);
                hc1D_jet->GetXaxis()->SetRangeUser(pt_min, pt_max);
                hl1D_jet->GetXaxis()->SetRangeUser(pt_min, pt_max);
                hg1D_jet->GetXaxis()->SetRangeUser(pt_min, pt_max);
                            
    //            TCanvas *c2 = new TCanvas();
    //            hq2D->Draw();
    //            c2->SaveAs("plots/test.pdf");

                // Project onto observable axis
                TH1* hc_proj = (TH1*) hc2D->ProjectionY();
                TH1* hl_proj = (TH1*) hl2D->ProjectionY();
                TH1* hg_proj = (TH1*) hg2D->ProjectionY();
                TH1* hc_unweighted = (TH1*) hc2D_unweighted->ProjectionY();
                TH1* hl_unweighted = (TH1*) hl2D_unweighted->ProjectionY();
                TH1* hg_unweighted = (TH1*) hg2D_unweighted->ProjectionY();

    //            TCanvas *c2 = new TCanvas();
    //            hq->Draw();
    //            c2->SaveAs("plots/test.pdf");

                // Set to appropriate name
                std::string hname = hc_proj->GetName();
                hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
                hc_proj->SetNameTitle(hname.c_str(), hname.c_str());
                hname = hl_proj->GetName();
                hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
                hl_proj->SetNameTitle(hname.c_str(), hname.c_str());
                hname = hg_proj->GetName();
                hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
                hg_proj->SetNameTitle(hname.c_str(), hname.c_str());

                hname = hc_unweighted->GetName();
                hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
                hc_unweighted->SetNameTitle(hname.c_str(), hname.c_str());
                hname = hl_unweighted->GetName();
                hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
                hl_unweighted->SetNameTitle(hname.c_str(), hname.c_str());
                hname = hg_unweighted->GetName();
                hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
                hg_unweighted->SetNameTitle(hname.c_str(), hname.c_str());

                // Rebin
                int n_obs_bins = 46; //-1;
    //            double obs_bins[51]; //[60];
                // double obs_bins[] = {1.00000000e-04, 1.20226443e-04, 1.44543977e-04, 1.73780083e-04, 2.08929613e-04, 2.51188643e-04, 3.01995172e-04, 3.63078055e-04, 4.36515832e-04, 5.24807460e-04, 6.30957344e-04, 7.58577575e-04, 9.12010839e-04, 1.09647820e-03, 1.31825674e-03, 1.58489319e-03, 1.90546072e-03, 2.29086765e-03, 2.75422870e-03, 3.31131121e-03, 3.98107171e-03, 4.78630092e-03, 5.75439937e-03, 6.91830971e-03, 8.31763771e-03, 1.00000000e-02, 1.20226443e-02, 1.44543977e-02, 1.73780083e-02, 2.08929613e-02, 2.51188643e-02, 3.01995172e-02, 3.63078055e-02, 4.36515832e-02, 5.24807460e-02, 6.30957344e-02, 7.58577575e-02, 9.12010839e-02, 1.09647820e-01, 1.31825674e-01, 1.58489319e-01, 1.90546072e-01, 2.29086765e-01, 2.75422870e-01, 3.31131121e-01, 3.98107171e-01, 4.78630092e-01, 5.75439937e-01, 6.91830971e-01, 8.31763771e-01, 1.00000000e+00};
                double obs_bins[] = {1.00000000e-04, 1.22167735e-04, 1.49249555e-04, 1.82334800e-04,
       2.22754295e-04, 2.72133877e-04, 3.32459793e-04, 4.06158599e-04,
       4.96194760e-04, 6.06189899e-04, 7.40568469e-04, 9.04735724e-04,
       1.10529514e-03, 1.35031404e-03, 1.64964807e-03, 2.01533769e-03,
       2.46209240e-03, 3.00788252e-03, 3.67466194e-03, 4.48925126e-03,
       5.48441658e-03, 6.70018750e-03, 8.18546731e-03, 1.00000000e-02,
       1.22167735e-02, 1.49249555e-02, 1.82334800e-02, 2.22754295e-02,
       2.72133877e-02, 3.32459793e-02, 4.06158599e-02, 4.96194760e-02,
       6.06189899e-02, 7.40568469e-02, 9.04735724e-02, 1.10529514e-01,
       1.35031404e-01, 1.64964807e-01, 2.01533769e-01, 2.46209240e-01,
       3.00788252e-01, 3.67466194e-01, 4.48925126e-01, 5.48441658e-01,
       6.70018750e-01, 8.18546731e-01, 1.00000000e+00};
                cout << "hc " << hc_proj->GetNbinsX() << endl;
                cout << "hl " << hl_proj->GetNbinsX() << endl;
                cout << "hg " << hg_proj->GetNbinsX() << endl;
                
                // TH1* hc = (TH1*) hc_proj->Rebin(n_obs_bins, (hc_name + "rebin").c_str(), obs_bins);
                // TH1* hl = (TH1*) hl_proj->Rebin(n_obs_bins, (hl_name + "rebin").c_str(), obs_bins);
                // TH1* hg = (TH1*) hg_proj->Rebin(n_obs_bins, (hg_name + "rebin").c_str(), obs_bins);
                // TH1* hc_unweighted = (TH1*) hc_unweighted_proj->Rebin(n_obs_bins, (hg_name + "rebin").c_str(), obs_bins);
                //unnecessary renaming for EEC...
                TH1* hc = (TH1*) hc_proj->Clone((hc_name + "rebin").c_str());
                TH1* hl = (TH1*) hl_proj->Clone((hl_name + "rebin").c_str());
                TH1* hg = (TH1*) hg_proj->Clone((hg_name + "rebin").c_str());

                // Find normalization factor
                double numjets_charm = hc1D_jet->Integral(); //GetEntries();
                double numjets_light = hl1D_jet->Integral(); //GetEntries();
                double numjets_gluon = hg1D_jet->Integral(); //GetEntries();

                // Set normalization
                //            hq->Scale(1/hq->Integral(), "width");
    //            hc->Scale(1/hc->Integral(), "width");
    //            hg->Scale(1/hg->Integral(), "width");

                hc->Scale(1/numjets_charm, "width");
                hl->Scale(1/numjets_light, "width");
                hg->Scale(1/numjets_gluon, "width");
                hc_unweighted->Scale(1/numjets_charm, "width");
                hl_unweighted->Scale(1/numjets_light, "width");
                hg_unweighted->Scale(1/numjets_gluon, "width");
                // cout << "num jets quark " << numjets_quark << endl;
                cout << "num jets charm " << numjets_charm << endl;
                cout << "num jets light " << numjets_light << endl;
                cout << "num jets gluon " << numjets_gluon << endl;


                // Calculate where the slope changes sign to find the max
                cout << "top of hc" << hc->GetBinCenter(findTopOfCurve(hc)) << endl;
                cout << "top of hl" << hl->GetBinCenter(findTopOfCurve(hl)) << endl;
                cout << "top of hg" << hg->GetBinCenter(findTopOfCurve(hg)) << endl;
                cout << "top of hc_unweighted" << hc_unweighted->GetBinCenter(findTopOfCurve(hc_unweighted)) << endl;
                cout << "top of hl_unweighted" << hl_unweighted->GetBinCenter(findTopOfCurve(hl_unweighted)) << endl;
                cout << "top of hg_unweighted" << hg_unweighted->GetBinCenter(findTopOfCurve(hg_unweighted)) << endl;

            


                /*
                // Find plotting range -- use gluon (broadest) distribution
                int maxbin = 0;
                double somevalue = 0.000001;
                int nxbins = hg->GetXaxis()->GetNbins();
                for (int j = 0; j < nxbins; j++) {
                    if (hg->GetBinContent(nxbins-j) > somevalue) {
                        maxbin = nxbins-j;
                        break;
                    }
                }
                hq->GetXaxis()->SetRange(1, maxbin);
                hc->GetXaxis()->SetRange(1, maxbin);
                hg->GetXaxis()->SetRange(1, maxbin);
                */
                // Set minimum / maximum y-axis value
    //            hq->SetMinimum(1e-3);

                maxy = include_gluon ? hg->GetMaximum() : 0;
                if (hc->GetMaximum() > maxy) maxy = hc->GetMaximum();
                if (hl->GetMaximum() > maxy) maxy = hl->GetMaximum();
                if (unweighted){
                    if (hc_unweighted->GetMaximum() > maxy) maxy = hc_unweighted->GetMaximum();
                    if (hl_unweighted->GetMaximum() > maxy) maxy = hl_unweighted->GetMaximum();
                    if (hg_unweighted->GetMaximum() > maxy) maxy = hg_unweighted->GetMaximum();
                }
                maxy *= 1.2; //1.5;
                hc->SetMaximum(maxy); //TODO: this needs to be fixed for when we have multiple files...
                cout << "maximums, hl: " << hl->GetMaximum() << " hc: " << hc->GetMaximum() << " hg: " << hg->GetMaximum() << endl;
                cout << "maximums, hl_unweighted: " << hl_unweighted->GetMaximum() << " hc_unweighted: " << hc_unweighted->GetMaximum() << " hg_unweighted: " << hg_unweighted->GetMaximum() << endl;
                cout << "chosen max: " << maxy << " " << maxy/1.5 << endl;
            

                // // Fix the x-axis label
                hc->GetXaxis()->SetTitle("#it{R}_{L}");
                hc->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");

                // Format histograms for plotting
                FormatHist(l, hc, "Charm-initiated jets", colors[0], markers[0]);
                if (switchON) {
                    l->AddEntry("NULL","          Charm decay on","h");
                } else {
                    l->AddEntry("NULL","          Charm decay off","h");
                }
                FormatHist(l, hl, "Light-initiated jets", colors[1], markers[1]);

                // Create plot
                hc->Draw("L same");
                hl->Draw("L same");

                // draw center
                double hc_top_binpos = findTopOfCurve(hc);
                double hl_top_binpos = findTopOfCurve(hl);
                drawVertLine(hc->GetBinCenter(hc_top_binpos), 0, hc->GetBinContent(hc_top_binpos), colors[0], 1)->Draw();
                drawVertLine(hl->GetBinCenter(hl_top_binpos), 0, hl->GetBinContent(hl_top_binpos), colors[1], 1)->Draw();


                // if we want to plot the gluon
                if (include_gluon) {
                    FormatHist(l, hg, "Gluon-initiated jets", kViolet+2, markers[2]);
                    hg->Draw("L same");
                    double hg_top_binpos = findTopOfCurve(hg);
                    drawVertLine(hg->GetBinCenter(hg_top_binpos), 0, hg->GetBinContent(hg_top_binpos), kViolet+2, 1)->Draw();

                }  

                // if we want unweighted plotted
                if (unweighted) {
                    // cout << "shouldn't be in here" << endl;
                    FormatHist(l, hc_unweighted, "Charm-initiated jets, unweighted", colors[6], 24); //fix the colors and markers on these
                    FormatHist(l, hl_unweighted, "Light-initiated jets, unweighted", colors[7], 26);
                    FormatHist(l, hg_unweighted, "Gluon-initiated jets, unweighted", kViolet+2, 27); //colors[8]

                    hc_unweighted->Draw("L same");
                    hl_unweighted->Draw("L same");
                    hg_unweighted->Draw("L same");

                    double hc_unweighted_top_binpos = findTopOfCurve(hc_unweighted);
                    double hl_unweighted_top_binpos = findTopOfCurve(hl_unweighted);
                    double hg_unweighted_top_binpos = findTopOfCurve(hg_unweighted);
                    drawVertLine(hc_unweighted->GetBinCenter(hc_unweighted_top_binpos), 0, hc_unweighted->GetBinContent(hc_unweighted_top_binpos), colors[6])->Draw();
                    drawVertLine(hl_unweighted->GetBinCenter(hl_unweighted_top_binpos), 0, hl_unweighted->GetBinContent(hl_unweighted_top_binpos), colors[7])->Draw();
                    drawVertLine(hg_unweighted->GetBinCenter(hg_unweighted_top_binpos), 0, hg_unweighted->GetBinContent(hg_unweighted_top_binpos), colors[8])->Draw();
                } 


                // TLegend *leg = new TLegend(0.1797168,0.6640741,0.4562155,0.8885185,"");
                // leg->AddEntry("NULL","PYTHIA 8","h");
                // leg->AddEntry("NULL","pp, #sqrt{#it{s}} = 13 TeV","h");
                // leg->AddEntry("NULL","in charged jets, anti-#it{k}_{T}, #it{R} = 0.4","h");
                // leg->AddEntry("NULL",ptbin,"h");
                // leg->SetTextSize(0.037);
                // leg->SetBorderSize(0);
                // leg->Draw("same");




            }


    


            



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

            if (pt_min == 10) { //} && grooming == "") {
                // Write rebinned histograms to root file
                f_out->cd();
                // if (!inclusive) {
                //     hl->Write();
                //     hg->Write();
                // }
            }
             
        } // pT bins loop
//            } // grooming loop
//        } // alpha loop
    } // jetR loop

    f->Close();
    delete f;
    if (inclusive) {
        f2->Close();
        delete f2;
    }


    return;
}
