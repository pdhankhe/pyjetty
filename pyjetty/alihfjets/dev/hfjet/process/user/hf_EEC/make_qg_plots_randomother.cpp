// ROOT macro to make quark-gluon jet plots
// this code is to write random other stuff
// right now:
//      - make comparison between D-tagged jets and charm-initiated jets with no decays
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


void make_qg_plots_randomother() {

//    gROOT->SetBatch(); //prevents plots from showing up
    gStyle->SetOptStat(0);
    SetStyle();
    Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
    Double_t marker_size = 1.5;
    Double_t colors[16] = {kRed, kGreen+2, kBlue, kRed+1, kGreen+1, kBlue+1, kRed+2, kGreen+2, kBlue+2, kRed+3, kGreen+3, kBlue+3, kOrange+1, kViolet+1, kYellow+1, kCyan+1};


    // File containing quark vs gluon histograms

     //FOR WHEN WEIGHTED/UNWEIGHTED IN SAME FILE
    const char infile_D0_Preeti[] = "/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alihfjets/dev/hfjet/process/user/hf_EEC/D0jet_EEC_15_30_ForBeatrice.root";
    const char infile_D0_1[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/17651853/AnalysisResultsFinal.root"; //this is using thnsparse - the original version
    const char infile_D0_2[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/18063568/AnalysisResultsFinal.root"; //this is using thnsparse - the second iteration
    const char infile_D0_3[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/24722441/AnalysisResultsFinal.root"; //this is using thnsparse - the third iteration
    const char infile_D0[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/24737979/AnalysisResultsFinal.root"; //this is using thnsparse - latest generation - USE THIS (it matches with #3)
    const char infile_Dstar[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/18008705/AnalysisResultsFinal.root"; //this is using thnsparse
    
    const char infile_D0_difNorm[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/17652140/AnalysisResultsFinal.root";
    const char infile_Dstar_difNorm[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/18008709/AnalysisResultsFinal.root";
    const char infile_D0wDstar_difNorm[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/18044499/AnalysisResultsFinal.root";
    const char infile_Dstar_difNorm_softpionremoved[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/18044503/AnalysisResultsFinal.root";
    const char infile_Dstar_difNorm_justsoftpion_noD0[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/18864762/AnalysisResultsFinal.root";
    const char infile_Dstar_difNorm_justsoftpionwD0[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/18864953/AnalysisResultsFinal.root";
    const char infile_Dstar_difNorm_justsoftpion[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/18865141/AnalysisResultsFinal.root";

    const char infile_noD0replacement_charmdecaysOFF[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/14680819/AnalysisResultsFinal.root"; //not thnsparse
    const char infile_noD0replacement_charmdecaysON[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/14680822/AnalysisResultsFinal.root";
    const char infile_noD0replacement_charmdecaysOFF_thnsparse_higherstats[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/24735947/AnalysisResultsFinal.root";
    const char infile_D0_feeddown[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/21046102/AnalysisResultsFinal.root";
    const char infile_phi[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/21077103/AnalysisResultsFinal.root";
    const char infile_D0_feeddown_ptrl[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/21078054/AnalysisResultsFinal.root";


    // int plot_case:
    // 0 = D0 vs c-init, charm decays off
    // 1 = D0 vs c-init, charm decays off, with multiple iterations of the jobs
        // also tryna figure out the difference bw the two D0 files

    //CONTOL VARIABLES HERE
    int plot_case = 1;

    TString label_D0 = "";
    TString label_D01 = "";
    TString label_D02 = "";
    TString label_D03 = "";
    TString label_c_nothnsparse = "";
    TString label_c = "";

    // const char file_list[4][] = {infile_D0, infile_D0_feeddown, infile_phi, infile_D0_feeddown_ptrl};
    TFile* f_D0 = new TFile(infile_D0, "READ");
    TFile* f_D01 = new TFile(infile_D0_1, "READ");
    TFile* f_D02 = new TFile(infile_D0_2, "READ");
    TFile* f_D03 = new TFile(infile_D0_3, "READ");
    TFile* f_charmdecaysOFF_nothnsparse = new TFile(infile_noD0replacement_charmdecaysOFF, "READ");
    TFile* f_charmdecaysOFF = new TFile(infile_noD0replacement_charmdecaysOFF_thnsparse_higherstats, "READ");
    // TFile* f_charmdecaysOFF = new TFile(infile_noD0replacement_charmdecaysOFF, "READ");
    
    // std::string add_name = std::string(file_list[plot_case]).substr(6);
    std::string add_name = "_D0_charmNOdecays_comparison";

    cout << "output name will be " << add_name << endl;

    // Output directory
    std::string outdir = "plots/final/other/";//"plots/test/"; 
    // Output file for binned results
    std::string outfile = outdir + "AnalysisResultsFinal" + add_name + ".root"; 
    TFile* f_out = new TFile(outfile.c_str(), "RECREATE");

    // Jet R value
    std::string jetR_list[] = { "0.4" };
    for (std::string jetR : jetR_list) {

        // Names of histograms in the file (quark, charm, gluon)
        const std::string hc_name = "h_EEC_JetPt_charm_R" + jetR;
        const std::string hc_unweighted_name = "h_EEC_JetPt_charm_R" + jetR + "_unweighted";
        const std::string hc_jet_name = "h_JetPt_charm_R" + jetR + "_jetlevel";
        
        const std::string hb_name = "h_EEC_JetPt_beauty_R" + jetR;
        const std::string hb_jet_name = "h_JetPt_beauty_R" + jetR + "_jetlevel";
        const std::string hi_name = "h_EEC_JetPt_inclusive_R" + jetR;
        const std::string hi_jet_name = "h_JetPt_inclusive_R" + jetR + "_jetlevel";

        const std::string hD0KpiNjets_name = "hD0KpiNjets"; //TODO later: = "hD0KpiNehD0KpiNjetsvents" for run 16729583



        const int pt_bins[] = { 7, 10, 15, 30 }; //{ 10, 20, 40 }; // CHANGE HERE!!
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
            

            
            //-------------------------------------------------//
            // find D0 reconstruction through charm
            
            THnSparse* hsparsejet_D0 = (THnSparse*) f_D0->Get(hc_name.c_str());
            THnSparse* hsparsejet_D0_jetlevel = (THnSparse*) f_D0->Get(hc_jet_name.c_str());
            THnSparse* hsparsejet_D01 = (THnSparse*) f_D01->Get(hc_name.c_str());
            THnSparse* hsparsejet_D01_jetlevel = (THnSparse*) f_D01->Get(hc_jet_name.c_str());
            THnSparse* hsparsejet_D02 = (THnSparse*) f_D02->Get(hc_name.c_str());
            THnSparse* hsparsejet_D02_jetlevel = (THnSparse*) f_D02->Get(hc_jet_name.c_str());
            THnSparse* hsparsejet_D03 = (THnSparse*) f_D03->Get(hc_name.c_str());
            THnSparse* hsparsejet_D03_jetlevel = (THnSparse*) f_D03->Get(hc_jet_name.c_str());

            TH2F* hcharm_nothnsparse = (TH2F*) f_charmdecaysOFF_nothnsparse->Get(hc_name.c_str());
            TH1F* hcharm_nothnsparse_jetlevel = (TH1F*) f_charmdecaysOFF_nothnsparse->Get(hc_jet_name.c_str());

            THnSparse* hsparsejet_charm = (THnSparse*) f_charmdecaysOFF->Get(hc_name.c_str());
            THnSparse* hsparsejet_charm_jetlevel = (THnSparse*) f_charmdecaysOFF->Get(hc_jet_name.c_str());

            // for THnSparse: make clone to work with, make cuts, get projection
            THnSparse *hsparsejet_D0_clone = (THnSparse *) hsparsejet_D0->Clone("hsparsejet_D0_clone");
            THnSparse *hsparsejet_D0_jetlevel_clone = (THnSparse *) hsparsejet_D0_jetlevel->Clone("hsparsejet_D0_jetlevel_clone");
            THnSparse *hsparsejet_D01_clone = (THnSparse *) hsparsejet_D01->Clone("hsparsejet_D01_clone");
            THnSparse *hsparsejet_D01_jetlevel_clone = (THnSparse *) hsparsejet_D01_jetlevel->Clone("hsparsejet_D01_jetlevel_clone");
            THnSparse *hsparsejet_D02_clone = (THnSparse *) hsparsejet_D02->Clone("hsparsejet_D02_clone");
            THnSparse *hsparsejet_D02_jetlevel_clone = (THnSparse *) hsparsejet_D02_jetlevel->Clone("hsparsejet_D02_jetlevel_clone");
            THnSparse *hsparsejet_D03_clone = (THnSparse *) hsparsejet_D03->Clone("hsparsejet_D03_clone");
            THnSparse *hsparsejet_D03_jetlevel_clone = (THnSparse *) hsparsejet_D03_jetlevel->Clone("hsparsejet_D03_jetlevel_clone");
            THnSparse *hsparsejet_charm_clone = (THnSparse *) hsparsejet_charm->Clone("hsparsejet_charm_clone");
            THnSparse *hsparsejet_charm_jetlevel_clone = (THnSparse *) hsparsejet_charm_jetlevel->Clone("hsparsejet_charm_jetlevel_clone");


            // get jet pT range
            hsparsejet_D0_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
            hsparsejet_D0_jetlevel_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
            hsparsejet_D01_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); 
            hsparsejet_D01_jetlevel_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max);
            hsparsejet_D02_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); 
            hsparsejet_D02_jetlevel_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max);
            hsparsejet_D03_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); 
            hsparsejet_D03_jetlevel_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max);

            hsparsejet_charm_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max);
            hsparsejet_charm_jetlevel_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max);
            hcharm_nothnsparse->GetXaxis()->SetRangeUser(pt_min, pt_max);
            hcharm_nothnsparse_jetlevel->GetXaxis()->SetRangeUser(pt_min, pt_max);
                
            // D meson cuts
            hsparsejet_D0_clone->GetAxis(1)->SetRangeUser(d0_pt_cuts[i], pt_max); // apply cut on Dmeson pt
            hsparsejet_D0_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
            hsparsejet_D0_jetlevel_clone->GetAxis(1)->SetRangeUser(d0_pt_cuts[i], pt_max);
            hsparsejet_D0_jetlevel_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8);
            hsparsejet_D01_clone->GetAxis(1)->SetRangeUser(d0_pt_cuts[i], pt_max);
            hsparsejet_D01_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8);
            hsparsejet_D01_jetlevel_clone->GetAxis(1)->SetRangeUser(d0_pt_cuts[i], pt_max);
            hsparsejet_D01_jetlevel_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8);
            hsparsejet_D02_clone->GetAxis(1)->SetRangeUser(d0_pt_cuts[i], pt_max);
            hsparsejet_D02_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8);
            hsparsejet_D02_jetlevel_clone->GetAxis(1)->SetRangeUser(d0_pt_cuts[i], pt_max);
            hsparsejet_D02_jetlevel_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8);
            hsparsejet_D03_clone->GetAxis(1)->SetRangeUser(d0_pt_cuts[i], pt_max);
            hsparsejet_D03_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8);
            hsparsejet_D03_jetlevel_clone->GetAxis(1)->SetRangeUser(d0_pt_cuts[i], pt_max);
            hsparsejet_D03_jetlevel_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8);


            // Project onto observable axis
            TH1D *hD0_proj = hsparsejet_D0_clone->Projection(3); //CALL THESE TH1*????
            TH1D *hc1D_jet = hsparsejet_D0_jetlevel_clone->Projection(0);
            TH1D *hD01_proj = hsparsejet_D01_clone->Projection(3); 
            TH1D *hc1D1_jet = hsparsejet_D01_jetlevel_clone->Projection(0);
            TH1D *hD02_proj = hsparsejet_D02_clone->Projection(3); 
            TH1D *hc1D2_jet = hsparsejet_D02_jetlevel_clone->Projection(0);
            TH1D *hD03_proj = hsparsejet_D03_clone->Projection(3); 
            TH1D *hc1D3_jet = hsparsejet_D03_jetlevel_clone->Projection(0);

            TH1D *hcharm_proj = hsparsejet_charm_clone->Projection(3); 
            TH1D *hcharm_jet = hsparsejet_charm_jetlevel_clone->Projection(0);
            TH1* hc_nothnsparse_proj = (TH1*) hcharm_nothnsparse->ProjectionY();

            // Set to appropriate name
            std::string hname = hD0_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hD0_proj->SetNameTitle(hname.c_str(), hname.c_str());

            hname = hD01_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hD01_proj->SetNameTitle(hname.c_str(), hname.c_str());

            hname = hD02_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hD02_proj->SetNameTitle(hname.c_str(), hname.c_str());

            hname = hD03_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hD03_proj->SetNameTitle(hname.c_str(), hname.c_str());

            hname = hcharm_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hcharm_proj->SetNameTitle(hname.c_str(), hname.c_str());

            hname = hc_nothnsparse_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hc_nothnsparse_proj->SetNameTitle(hname.c_str(), hname.c_str());


            // Rebin
            int n_obs_bins = 13; //-1;
            double obs_bins[14] = {0.0001, 0.00020892961308540387, 0.0004365158322401661, 0.0009120108393559096, 
                    0.0019054607179632482, 0.005754399373371567, 0.017378008287493762, 0.03630780547701014, 
                    0.07585775750291836, 0.15848931924611143, 0.3311311214825911, 0.47863009232263853, 0.6918309709189363, 1.0};

            cout << "About to rebin" << endl;
            // TH1* hD0 = (TH1*) hD0_proj->Rebin(n_obs_bins, (hname + "rebin").c_str(), obs_bins);
            TH1* hD0 = (TH1*) hD0_proj->Clone( hD0_proj->GetName() );
            TH1* hD01 = (TH1*) hD01_proj->Clone( hD01_proj->GetName() );
            TH1* hD02 = (TH1*) hD02_proj->Clone( hD02_proj->GetName() );
            TH1* hD03 = (TH1*) hD03_proj->Clone( hD03_proj->GetName() );
            TH1* hc = (TH1*) hcharm_proj->Clone( hcharm_proj->GetName() );
            TH1* hc_nothnsparse = (TH1*) hc_nothnsparse_proj->Clone( hc_nothnsparse_proj->GetName() );
            cout << "Rebin done" << endl;

            // Find normalization factor
            double numjets_D0 = hc1D_jet->Integral();
            double numjets_D01 = hc1D1_jet->Integral();
            double numjets_D02 = hc1D2_jet->Integral();
            double numjets_D03 = hc1D3_jet->Integral();
            double numjets_charm = hcharm_jet->Integral();
            double numjets_charm_nothnsparse = hcharm_nothnsparse_jetlevel->Integral();
            hD0->Scale(1/numjets_D0, "width");
            hD01->Scale(1/numjets_D01, "width");
            hD02->Scale(1/numjets_D02, "width");
            hD03->Scale(1/numjets_D03, "width");
            hc->Scale(1/numjets_charm, "width");
            hc_nothnsparse->Scale(1/numjets_charm_nothnsparse, "width");




            //Format color and style
            int markercolor0 = kGreen-5; //D0
            int markerstyle0 = kFullCircle;
            int markercolor1 = kGreen-7; //D01
            int markerstyle1 = kFullCircle;
            int markercolor2 = kGreen+2; //D02
            int markerstyle2 = kOpenCircle;
            int markercolor3 = kTeal; //D03
            int markerstyle3 = kOpenCircle;
            int markercolor4 = kRed; //charm
            int markerstyle4 = kFullCircle;
            int markercolor5 = kRed-9; //charm - no thnsparse
            int markerstyle5 = kOpenCircle;
            
            label_D0 = "D^{0}-tagged, c-init jets";
            label_D01 = "D^{0}-tagged [1], c-init jets";
            label_D02 = "D^{0}-tagged [2], c-init jets";
            label_D03 = "D^{0}-tagged [3], c-init jets";
            label_c = "c-init jets, charm decays off";
            label_c_nothnsparse = "c-init jets [lower stats], charm decays off";

            

            // Format histograms for plotting (this order needed to keep legend in order and graphs lookin good)
            hD0->GetXaxis()->SetTitle("#it{R}_{L}");
            hD01->GetXaxis()->SetTitle("#it{R}_{L}");
            hD02->GetXaxis()->SetTitle("#it{R}_{L}");
            hD03->GetXaxis()->SetTitle("#it{R}_{L}");
            hc->GetXaxis()->SetTitle("#it{R}_{L}");
            hc_nothnsparse->GetXaxis()->SetTitle("#it{R}_{L}");
            hD0->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
            hD01->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
            hD02->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
            hD03->GetXaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
            hc->GetXaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
            hc_nothnsparse->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");

            cout << "about to format D0" << endl;
            FormatHist(l, hD0, label_D0, markercolor0, markerstyle0); //FormatHist(l, hD0, "D^{0}-tagged, c-init jets", kMagenta+3, kOpenSquare);
            if ( plot_case == 1 ) {
                FormatHist(l, hD01, label_D01, markercolor1, markerstyle1);
                FormatHist(l, hD02, label_D02, markercolor2, markerstyle2);
                FormatHist(l, hD03, label_D03, markercolor3, markerstyle3);
                FormatHist(l, hc_nothnsparse, label_c_nothnsparse, markercolor5, markerstyle5);
            }
            FormatHist(l, hc, label_c, markercolor4, markerstyle4);
            if (plot_case == 0){
                l->AddEntry("NULL","          D* decays off","h");
            } 
            hD0->SetMaximum(hD0->GetMaximum()*1.5);
            hD0->Draw("L same");
            if ( plot_case == 1 ) {
                hD01->Draw("L same");
                hD02->Draw("L same");
                hD03->Draw("L same");
                hc_nothnsparse->Draw("L same");
            }
            hc->Draw("L same");

            
            
            double hD0_top_binpos = findTopOfCurve(hD0);
            drawVertLine(hD0->GetBinCenter(hD0_top_binpos), 0, hD0->GetBinContent(hD0_top_binpos), markercolor1, 1)->Draw();
            if ( plot_case == 1 ) {
                double hD01_top_binpos = findTopOfCurve(hD01);
                drawVertLine(hD01->GetBinCenter(hD01_top_binpos), 0, hD01->GetBinContent(hD01_top_binpos), markercolor2, 1)->Draw();
                double hD02_top_binpos = findTopOfCurve(hD02);
                drawVertLine(hD02->GetBinCenter(hD02_top_binpos), 0, hD02->GetBinContent(hD02_top_binpos), markercolor2, 1)->Draw();
                double hD03_top_binpos = findTopOfCurve(hD03);
                drawVertLine(hD03->GetBinCenter(hD03_top_binpos), 0, hD03->GetBinContent(hD03_top_binpos), markercolor3, 1)->Draw();
                double hc_nothnsparse_top_binpos = findTopOfCurve(hc_nothnsparse);
                drawVertLine(hc_nothnsparse->GetBinCenter(hc_nothnsparse_top_binpos), 0, hc_nothnsparse->GetBinContent(hc_nothnsparse_top_binpos), markercolor5, 1)->Draw();
            }
            double hc_top_binpos = findTopOfCurve(hc);
            drawVertLine(hc->GetBinCenter(hc_top_binpos), 0, hc->GetBinContent(hc_top_binpos), markercolor4, 1)->Draw();
            

            
            // vector<double> fullwidth_vec = findWidthOfCurve(hD0,  hD0_top_binpos);
            // drawHoriLine(fullwidth_vec[1], fullwidth_vec[2], fullwidth_vec[0], kMagenta+3, 1)->Draw();
            
            

            // draw legend
            l->Draw("same");



            std::string fname = outdir + "QG_comp_pt" + std::to_string(pt_min) + '-' + std::to_string(pt_max) + "_R" + jetR + add_name + "_plotcase" + std::to_string(plot_case) + ".pdf";
            const char* fnamec = fname.c_str();
            c->SaveAs(fnamec);
            delete c;

            // if (pt_min == 10) { //} && grooming == "") {
                // Write rebinned histograms to root file
            f_out->cd();
            hD0->Write();
            hD01->Write();
            hD02->Write();
            hD03->Write();
            hc_nothnsparse->Write();
            hc->Write();
            // }
             
        } // pT bins loop
    } // jetR loop

    f_D0->Close();
    delete f_D0;
    f_D01->Close();
    delete f_D01;
    f_D02->Close();
    delete f_D02;
    f_D03->Close();
    delete f_D03;
    f_charmdecaysOFF->Close();
    delete f_charmdecaysOFF;
    f_charmdecaysOFF_nothnsparse->Close();
    delete f_charmdecaysOFF_nothnsparse;

    return;
}
