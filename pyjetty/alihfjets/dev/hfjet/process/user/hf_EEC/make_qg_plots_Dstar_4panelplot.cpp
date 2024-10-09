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

void FormatHist(TLegend *l, TH1 *hist, TString text, int markercolor=1, int markerstyle=8, bool yaxisvisible=true) 
{
    hist->SetLineColor(markercolor);
    hist->SetMarkerColor(markercolor);
    hist->SetMarkerStyle(markerstyle);
    hist->SetMarkerSize(1.0);
    l->AddEntry(hist, text, "pl");

	//gPad->SetTickx(); 
	//gPad->SetTicky(); 
	// h->SetLineWidth(2);
    if (yaxisvisible) {
        hist->GetYaxis()->SetTitleOffset(1.1); 
        hist->GetYaxis()->SetTitleSize(0.06); //(0.042);
        hist->GetYaxis()->SetLabelSize(0.05); //(0.042);
        hist->GetYaxis()->SetLabelFont(42);
        hist->GetYaxis()->SetTitleFont(42);
    } else {
        hist->GetYaxis()->SetTitleOffset(1.1);
        hist->GetYaxis()->SetTitleSize(0); //(0.042);
        hist->GetYaxis()->SetLabelSize(0); //(0.042);
    }
	hist->GetXaxis()->SetLabelFont(42);
	hist->GetXaxis()->SetTitleFont(42);
	hist->GetXaxis()->SetTitleOffset(1.05);
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

TLine * drawVertLine(double x1, double y1, double y2, int color, int linestyle=2, double linewidth=1){
    auto fvertline = new TLine(x1, y1, x1, y2);
	fvertline->SetLineWidth(linewidth);
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
TH1D * getObsHist(TFile *filename, std::string h_name, std::string h_jet_name, int pt_min, int pt_max, int d0_pt_cut, std::string newhistname, 
                  bool d0cuts=false, int obsaxis=3, bool ptrl=false, double cutonxlow=1e-4, double cutonxup=1.) {
    
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
    // hist->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
    hist->GetYaxis()->SetTitle("#Sigma_{EEC}(#it{R}_{L})");
    hist->GetXaxis()->SetRangeUser(cutonxlow, cutonxup);

    return hist;
}


void FormatPad(TPad * pad, double leftcut, double bottomcut, double rightcut, double topcut, bool logx=false, bool logy=false) {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // gStyle->SetOptLabel(0);
    
    pad->SetLeftMargin(leftcut);
    pad->SetBottomMargin(bottomcut);
    pad->SetRightMargin(rightcut);
    pad->SetTopMargin(topcut);

    if (logx) pad->SetLogx();
    if (logy) pad->SetLogy();
}

void make_qg_plots_Dstar_4panelplot() {

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
    // 50 = put cases 2,3,4,7 into one 4-panel plot --> go to new file make_qg_plots_Dstar_4panelplot.cpp

    //CONTOL VARIABLES HERE
    int plot_case = 50;
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

    // if (plot_case == 2) {
    TFile *f_D0 = new TFile(infile_D0, "READ"); //2, 3
    TFile *f_Dstar_difNorm = new TFile(infile_Dstar_difNorm, "READ"); //2, 4, 7
    TFile *f_D0wDstar_difNorm = new TFile(infile_D0wDstar_difNorm, "READ"); //3
    TFile *f_Dstar_difNorm_softpionremoved = new TFile(infile_Dstar_difNorm_softpionremoved, "READ"); //4
    TFile *f_Dstar_difNorm_justsoftpion = new TFile(infile_Dstar_difNorm_justsoftpion, "READ"); //7

    add_name = "_Dstar_plot_case50.pdf";
    cout << "output name will be " << add_name << endl;

    // Output directory
    //std::string outdir = "/rstorage/alice/AnalysisResults/ang/1224559/plots/";
    std::string outdir = "plots/final/poster_plots/";//"plots/test/";
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


        // testing - look at # jets before cuts
        // cout << "numDtaggedjets from hist before cuts " << hsparsejet_c_jetlevel->Projection(0)->GetEntries() << endl;



        //const int pt_bins[] = { 10, 20, 40, 60, 80, 100, 150 };
        const int pt_bins[] = { 7, 10, 15, 30 }; //{ 10, 20, 40 };
        const int d0_pt_cuts[] = { 3, 5, 5 };
        const int n_bins = 3; //2;
        for (int i = 0; i < n_bins; i++) {
            cout << "in pt bin" << i << endl;
            int pt_min = pt_bins[i];
            int pt_max = pt_bins[i+1];
            std::string pt_name = "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);

            // define pt related variables
            TString ptbin = TString::Format("%d #leq #it{p}_{T}^{ch. jet} < %d GeV/#it{c}, #font[122]{|}#it{#eta}_{jet}#font[122]{|} #leq 0.5", pt_min, pt_max);
            TString ptD = TString::Format("%d #leq #it{p}_{T}^{D^{0}} < %d GeV/#it{c}, #font[122]{|}#it{y}_{D^{0}}#font[122]{|} #leq 0.8", d0_pt_cuts[i], pt_max);


            // make a canvas for each pt range
            TCanvas* c_4panel = new TCanvas();
            ProcessCanvas(c_4panel);
            // c_4panel->Divide(2,2);
            

            TLegend* l; // = new TLegend(0.17, 0.65, 0.5, 0.85);

            double maxy = 0;

            // start files loop
            // for (int ifile=0; ifile<files.size(); ifile++) {


            // Open histograms


            l = new TLegend(0.2,0.4,0.4,0.75,""); //(0.17, 0.4, 0.5, 0.53);
            l->SetTextSize(0.04);
            l->SetBorderSize(0);
            l->AddEntry("NULL","PYTHIA 8 Monash 2013","h");
            l->AddEntry("NULL","pp, #sqrt{#it{s}} = 13 TeV","h");
            l->AddEntry("NULL","  ","h");
            l->AddEntry("NULL","D^{0} #rightarrow K^{#minus} #pi^{+} and charge conj.","h");
            l->AddEntry("NULL","anti-#it{k}_{T} ch. jets, #it{R} = 0.4","h");
            l->AddEntry("NULL",ptbin,"h");

            TLegend *l_main2 = new TLegend(0.2,0.27,0.44,0.38,""); //(0.17, 0.4, 0.5, 0.53);
            l_main2->SetTextSize(0.04);
            l_main2->SetBorderSize(0);
            l_main2->AddEntry("NULL",ptD,"h");
            l_main2->AddEntry("NULL","  ","h");
            
            // l->Draw("same");


            //-------------------------------------------------//



            // Find normalization factor
            // double numjets_charm = hc1D_jet->Integral();
            // double numDtaggedjets = hD0KpiNjets->GetEntries(); // this should match the # from projection before cuts
            // double numDtaggedjets_fromhist = hc1D_jet->GetEntries();
            // double numDtaggedjets_fromhist_eff = hc1D_jet->GetEffectiveEntries();

            // cout << "numjets_charm " << numjets_charm << endl;
            // cout << "numDtaggedjets " << numDtaggedjets << endl;
            // cout << "numDtaggedjets_fromhist " << numDtaggedjets_fromhist << endl;
            // cout << "numDtaggedjets_fromhist_eff " << numDtaggedjets_fromhist_eff << endl;

            
            // ------------------------ DOING WORK ----------------------------------
            TH1D *h_D0 = getObsHist(f_D0, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_D0" + pt_name, true, 3, false, 1.3e-4, 1.1); //note projection axis
            TH1D *h_Dstar_difNorm = getObsHist(f_Dstar_difNorm, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_Dstar_difNorm" + pt_name, true, 3, false, 1.3e-4, 1.1); //note projection axis
            TH1D *h_D0wDstar_difNorm = getObsHist(f_D0wDstar_difNorm, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_D0wDstar_difNorm" + pt_name, true, 3, false, 1.3e-4, 1.1); //note projection axis
            TH1D *h_Dstar_difNorm_softpionremoved = getObsHist(f_Dstar_difNorm_softpionremoved, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_Dstar_difNorm_softpionremoved" + pt_name, true, 3, false, 1.3e-4, 1.1); //note projection axis
            TH1D *h_Dstar_difNorm_justsoftpion = getObsHist(f_Dstar_difNorm_justsoftpion, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_Dstar_difNorm_justsoftpion" + pt_name, true, 3, false, 1.3e-4, 1.1); //note projection axis

            // ----------------------------------------------------------------------


            // cout << "numDtaggedjets_fromhist_Dstar " << numDtaggedjets_fromhist_charm << endl;




            //Format color and style
            int markercolor_D0 = kGreen-5; //D0, me
            int markerstyle_D0 = kFullCircle;
            std::string label_D0 = "D^{0}-tagged"; //, c-init jets";
            int markercolor_Dstar_difNorm = kRed-3; //D*
            int markerstyle_Dstar_difNorm = 29;
            std::string label_Dstar_difNorm = "D*-tagged"; //, c-init jets";
            int markercolor_D0wDstar_difNorm = kRed-3; //D*+D0
            int markerstyle_D0wDstar_difNorm = 29;
            std::string label_D0wDstar_difNorm = "(D^{0}+D*)-tagged"; //, c-init jets";

            int markercolor_Dstar_difNorm_softpionremoved = kViolet+2; //D*, soft pion removed
            int markerstyle_Dstar_difNorm_softpionremoved = 33;
            std::string label_Dstar_difNorm_softpionremoved = "D*-tagged, soft #pi removed"; //, c-init jets";
            int markercolor_Dstar_difNorm_justsoftpion = kViolet+2; //D*, only soft pion
            int markerstyle_Dstar_difNorm_justsoftpion = 33;
            std::string label_Dstar_difNorm_justsoftpion = "D*-tagged w/ only soft #pi"; //, c-init jets";

            TLegend *l_topleft = new TLegend(0.18,0.15,0.38,0.3,"");
            TLegend *l_topright = new TLegend(0.03,0.6,0.18,0.75,"");
            TLegend *l_bottomleft = new TLegend(0.18,0.85,0.38,0.95,"");
            TLegend *l_bottomright = new TLegend(0.03,0.85,0.18,0.95,"");
            l_topleft->SetTextSize(0.04);
            l_topright->SetTextSize(0.04);
            l_bottomleft->SetTextSize(0.04);
            l_bottomright->SetTextSize(0.04);

        

            // calculate max peak position lines
            double D0_top_binpos = findTopOfCurve(h_D0);
            double Dstar_difNorm_top_binpos = findTopOfCurve(h_Dstar_difNorm);
            double D0wDstar_difNorm_top_binpos = findTopOfCurve(h_D0wDstar_difNorm);
            double Dstar_difNorm_softpionremoved_top_binpos = findTopOfCurve(h_Dstar_difNorm_softpionremoved);
            double Dstar_difNorm_justsoftpion_top_binpos = findTopOfCurve(h_Dstar_difNorm_justsoftpion);
            

            // draw in panels
            c_4panel->cd();
            // c_4panel->cd(3);
            // TPad *p3 = new TPad("p3","p3", 0., 0., 1., 1.);
            TPad *p3 = new TPad("p3","p3", 0., 0., 0.5, 0.5);
            p3->Draw();
            cout << "canv 3" << endl;
            // FormatPad(p3, 0, 0, 0.5, 0.5, true, false);
            FormatPad(p3, 0.16, 0.2, 0, 0, true, false);
            p3->cd();

            FormatHist(l_bottomleft, h_Dstar_difNorm, label_Dstar_difNorm, markercolor_Dstar_difNorm, markerstyle_Dstar_difNorm); //"D*-tagged, c-init jets", kRed-7, 29);
            FormatHist(l_bottomleft, h_Dstar_difNorm_softpionremoved, label_Dstar_difNorm_softpionremoved, markercolor_Dstar_difNorm_softpionremoved, markerstyle_Dstar_difNorm_softpionremoved); //"D*-tagged, c-init jets", kRed-7, 29);
            h_Dstar_difNorm->Draw("L same");
            h_Dstar_difNorm_softpionremoved->Draw("L same");
            drawVertLine(h_Dstar_difNorm->GetBinCenter(Dstar_difNorm_top_binpos), 0, h_Dstar_difNorm->GetBinContent(Dstar_difNorm_top_binpos), markercolor_Dstar_difNorm, 3, 1)->Draw();
            drawVertLine(h_Dstar_difNorm_softpionremoved->GetBinCenter(Dstar_difNorm_softpionremoved_top_binpos), 0, h_Dstar_difNorm_softpionremoved->GetBinContent(Dstar_difNorm_softpionremoved_top_binpos), markercolor_Dstar_difNorm_softpionremoved, 3, 1)->Draw();
            l_bottomleft->Draw("same");


            c_4panel->cd();
            // c_4panel->cd(4);
            // TPad *p4 = new TPad("p4","p4", 0., 0., 1., 1.);
            TPad *p4 = new TPad("p4","p4", 0.5, 0., 1., 0.5);
            p4->Draw();
            cout << "canv 4" << endl;cout << "canv 4" << endl;
            // FormatPad(p4, 0.5, 0, 0, 0.5, true, false);
            FormatPad(p4, 0, 0.2, 0.16, 0, true, false);
            p4->cd();

            FormatHist(l_bottomright, h_Dstar_difNorm, label_Dstar_difNorm, markercolor_Dstar_difNorm, markerstyle_Dstar_difNorm); //"D*-tagged, c-init jets", kRed-7, 29);
            FormatHist(l_bottomright, h_Dstar_difNorm_justsoftpion, label_Dstar_difNorm_justsoftpion, markercolor_Dstar_difNorm_justsoftpion, markerstyle_Dstar_difNorm_justsoftpion); //"D*-tagged, c-init jets", kRed-7, 29);
            h_Dstar_difNorm->Draw("L same");
            h_Dstar_difNorm_justsoftpion->Draw("L same");
            drawVertLine(h_Dstar_difNorm->GetBinCenter(Dstar_difNorm_top_binpos), 0, h_Dstar_difNorm->GetBinContent(Dstar_difNorm_top_binpos), markercolor_Dstar_difNorm, 3, 1)->Draw();
            drawVertLine(h_Dstar_difNorm_justsoftpion->GetBinCenter(Dstar_difNorm_justsoftpion_top_binpos), 0, h_Dstar_difNorm_justsoftpion->GetBinContent(Dstar_difNorm_justsoftpion_top_binpos), markercolor_Dstar_difNorm_justsoftpion, 3, 1)->Draw();
            l_bottomright->Draw("same");

            h_D0->SetMaximum(h_D0->GetMaximum()*1.1);
            h_Dstar_difNorm->SetMaximum(h_Dstar_difNorm->GetMaximum()*1.1);


            c_4panel->cd();
            // TPad *p1 = new TPad("p1","p1", 0., 0., 1., 1.);
            TPad *p1 = new TPad("p1","p1", 0., 0.5, 0.5, 1.);
            p1->Draw();
            cout << "canv 1" << endl;
            // FormatPad(p1, 0, 0.5, 0.5, 0, true, false);
            FormatPad(p1, 0.16, 0, 0, 0.2, true, false);
            p1->cd();

            FormatHist(l_topleft, h_D0, label_D0, markercolor_D0, markerstyle_D0); //FormatHist(l, hD0, "D^{0}-tagged, c-init jets", kMagenta+3, kOpenSquare);
            l_topleft->AddEntry("NULL","          D* decays off","h");
            // cout << "about to format Dstar" << endl;
            FormatHist(l_topleft, h_Dstar_difNorm, label_Dstar_difNorm, markercolor_Dstar_difNorm, markerstyle_Dstar_difNorm); //"D*-tagged, c-init jets", kRed-7, 29);
            h_D0->Draw("L same");
            h_Dstar_difNorm->Draw("L same");
            drawVertLine(h_D0->GetBinCenter(D0_top_binpos), 0, h_D0->GetBinContent(D0_top_binpos), markercolor_D0, 3, 1)->Draw();
            drawVertLine(h_Dstar_difNorm->GetBinCenter(Dstar_difNorm_top_binpos), 0, h_Dstar_difNorm->GetBinContent(Dstar_difNorm_top_binpos), markercolor_Dstar_difNorm, 3, 1)->Draw();
            
            // draw legend
            l->Draw("same");
            l_main2->Draw("same");
            l_topleft->Draw("same");


            c_4panel->cd();
            // c_4panel->cd(2);
            // TPad *p2 = new TPad("p2","p2", 0., 0., 1., 1.);
            TPad *p2 = new TPad("p2","p2", 0.5, 0.5, 1., 1.);
            p2->Draw();
            cout << "canv 2" << endl;
            // FormatPad(p2, 0.5, 0.5, 0, 0, true, false);
            FormatPad(p2, 0, 0, 0.16, 0.2, true, false);
            p2->cd();

            FormatHist(l_topright, h_D0, label_D0, markercolor_D0, markerstyle_D0); //FormatHist(l, hD0, "D^{0}-tagged, c-init jets", kMagenta+3, kOpenSquare);
            l_topright->AddEntry("NULL","          D* decays off","h");
            FormatHist(l_topright, h_D0wDstar_difNorm, label_D0wDstar_difNorm, markercolor_D0wDstar_difNorm, markerstyle_D0wDstar_difNorm); //"D*-tagged, c-init jets", kRed-7, 29);
            h_D0->Draw("L same");
            h_D0wDstar_difNorm->Draw("L same");
            drawVertLine(h_D0->GetBinCenter(D0_top_binpos), 0, h_D0->GetBinContent(D0_top_binpos), markercolor_D0, 3, 1)->Draw();
            drawVertLine(h_D0wDstar_difNorm->GetBinCenter(D0wDstar_difNorm_top_binpos), 0, h_D0wDstar_difNorm->GetBinContent(D0wDstar_difNorm_top_binpos), markercolor_D0wDstar_difNorm, 3, 1)->Draw();
            l_topright->Draw("same");


            c_4panel->cd();
            TPaveText *text_stupid = new TPaveText(0.07, 0.45, 0.075, 0.55, "NDC");
            text_stupid->SetFillColorAlpha(kWhite, 0);
            text_stupid->SetTextFont(42);
            text_stupid->SetTextSize(0.025); //35);
            text_stupid->SetBorderSize(0.);
            text_stupid->AddText("0");
            text_stupid->Draw();


            


        
            



            // } //end of file loop?

            


            std::string fname = outdir + "QG_comp_pt" + std::to_string(pt_min) + '-' + std::to_string(pt_max) + "_R" + jetR + add_name; //"_charmdecaysONcomparison.pdf"; // + "_normbytype.pdf"; //"_nonorm.pdf";
            const char* fnamec = fname.c_str();
            c_4panel->SaveAs(fnamec);
            delete c_4panel;

            // if (pt_min == 10) { //} && grooming == "") {
                // Write rebinned histograms to root file
                f_out->cd();

                h_D0->Write();
                h_Dstar_difNorm->Write();
                h_D0wDstar_difNorm->Write();
                h_Dstar_difNorm_softpionremoved->Write();
                h_Dstar_difNorm_justsoftpion->Write();
                            // }
            // }
             
        } // pT bins loop

//            } // grooming loop
//        } // alpha loop
    } // jetR loop

    f_D0->Close();
    f_Dstar_difNorm->Close();
    f_D0wDstar_difNorm->Close();
    f_Dstar_difNorm_softpionremoved->Close();
    f_Dstar_difNorm_justsoftpion->Close();
            
    delete f_D0;
    delete f_Dstar_difNorm;
    delete f_D0wDstar_difNorm;
    delete f_Dstar_difNorm_softpionremoved;
    delete f_Dstar_difNorm_justsoftpion;



    return;
}
