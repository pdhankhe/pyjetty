// ROOT macro to make quark-gluon jet plots
// Beatrice Liang-Gilman (beatrice_lg@berkeley.edu)

#include <typeinfo>

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
  	gStyle->SetLabelOffset(0.01,"y"); //(0.005,"y"); 
  	gStyle->SetLabelOffset(0.01,"x"); //(0.005,"x");
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

void addLegendInfo(TLegend *l, TString ptbin, TString ptD, bool D0recon) {
    l->SetTextSize(0.045);
    // TLegend *leg = new TLegend(0.1797168,0.5390741,0.4562155,0.8885185,"");
    l->AddEntry("NULL","PYTHIA 8 Monash 2013","h");
    l->AddEntry("NULL","pp, #sqrt{#it{s}} = 13 TeV","h");
    if (D0recon) l->AddEntry("NULL","D^{0} #rightarrow K^{#minus} #pi^{+} and charge conj.","h");
    l->AddEntry("NULL","in charged jets, anti-#it{k}_{T}, #it{R} = 0.4","h");
    l->AddEntry("NULL",ptbin,"h");
    if (D0recon) l->AddEntry("NULL",ptD,"h");
    l->SetTextSize(0.037);
    l->SetBorderSize(0);
    l->SetFillStyle(0); // turn legend transparent
    // l->Draw("same");
}


void make_qg_plots_write() {

//    gROOT->SetBatch(); //prevents plots from showing up
    gStyle->SetOptStat(0);
    SetStyle();
    Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
    Double_t marker_size = 1.5;
    Double_t colors[16] = {kRed, kGreen+2, kBlue, kRed+1, kGreen+1, kBlue+1, kRed+2, kGreen+2, kBlue+2, kRed+3, kGreen+3, kBlue+3, kOrange+1, kViolet+1, kYellow+1, kCyan+1};


    // int plot_case:
    // 0 = D0 replacement (primordial + from D* D0's)
    // 1 = charm decays DISABLED
    // 2 = charm decays ENABLED
    // int zcut_case:
    // 1 = zcut=0.1
    // 2 = zcut=0.2

    //CONTOL VARIABLES HERE
    int plot_case = 0;
    int zcut_case = 1;

    // File containing quark vs gluon histograms

     //FOR WHEN WEIGHTED/UNWEIGHTED IN SAME FILE
    // const char infile_D0[128];
    // const char infile_charmOFF[128];
    // const char infile_charmON[128];
    // static const char * const infile_D0[];
    // static const char * const infile_charmOFF[];
    // static const char * const infile_charmON[];
    //zcut=0.1
    // if (zcut_case == 1) {
        const char infile_D0[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jetaxis/22078282/AnalysisResultsFinal.root";
        const char infile_charmOFF[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jetaxis/22078315/AnalysisResultsFinal.root"; //perlmutter
        const char infile_charmON[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jetaxis/22078314/AnalysisResultsFinal.root"; //perlmutter 
    // } else if (zcut_case == 2) { //zcut=0.2
        // const char infile_D0[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jetaxis/22060892/AnalysisResultsFinal.root";
        // const char infile_charmOFF[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jetaxis/22060858/AnalysisResultsFinal.root"; //perlmutter
        // const char infile_charmON[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jetaxis/22060888/AnalysisResultsFinal.root"; //perlmutter 
    // }

    TFile* f;
    TString label1 = "";

    TString file_list[3] = {string(infile_D0), string(infile_charmOFF), 
                                string(infile_charmON)};
    std::string output_dirs[] = { "D0", "charmOFF", "charmON" };
    f = new TFile(file_list[plot_case], "READ");

    std::string add_name;
    if (plot_case == 0) {
        add_name = "_D0";
    } else if (plot_case == 1) {
        add_name = "_charmOFF";
    } else if (plot_case == 2) {
        add_name = "_charmON";
    } 

    cout << "output name will be " << add_name << endl;

    // Output directory
    std::string outdir = "plots/" + output_dirs[plot_case] + "/"; //final/";//"plots/test/"; 
    // Output file for binned results
    std::string outfile = outdir + "AnalysisResultsFinal" + add_name + ".root"; 
    TFile* f_out = new TFile(outfile.c_str(), "RECREATE");

    const int numjetaxes = 6;

    // Jet R value
    std::string jetR_list[] = { "0.4" };
    for (std::string jetR : jetR_list) {

        // Names of histograms in the file (quark, charm, gluon)
        // jetaxis_names[6];
        // if (zcut_case == 1) {
        // const std::string jetaxis_names[] = { "SD-D_SD_zcut02_B0", "STD-D",
        //         "STD-SD_SD_zcut02_B0", "STD-WTA", "WTA-D", "WTA-SD_SD_zcut02_B0" };
        // } else if (zcut_case == 2) {
        const std::string jetaxis_names[] = { "SD-D_SD_zcut01_B0", "STD-D",
                "STD-SD_SD_zcut01_B0", "STD-WTA", "WTA-D", "WTA-SD_SD_zcut01_B0" };
        // }

        std::vector<std::string> hc_names;
        std::vector<std::string> hg_names;
        std::vector<std::string> hl_names;
        std::vector<std::string> hi_names;
        for (int iobs=0; iobs<numjetaxes; iobs++) {
            const std::string hc_name = "h_jet_axis_JetPt_charm_R" + jetR + "_" + jetaxis_names[iobs];            
            const std::string hg_name = "h_jet_axis_JetPt_gluon_R" + jetR + "_" + jetaxis_names[iobs];            
            const std::string hl_name = "h_jet_axis_JetPt_light_R" + jetR + "_" + jetaxis_names[iobs];            
            const std::string hi_name = "h_jet_axis_JetPt_inclusive_R" + jetR + "_" + jetaxis_names[iobs];            
            
            hc_names.push_back(hc_name);
            hg_names.push_back(hg_name);
            hl_names.push_back(hl_name);
            hi_names.push_back(hi_name);
        }

        // const std::string htest_name = "h_jet_axis_JetPt_charm_R" + jetR + "_" + jetaxis_names[0];            
        // cout << "test name" << htest_name << endl;

        const std::string hD0KpiNjets_name = "hD0KpiNjets";



        const int pt_bins[] = { 5, 7, 10, 20, 50 }; //{ 10, 15, 30 }; //{ 10, 20, 40 }; // CHANGE HERE!!
        const int n_bins = 4; 
        const std::string pt_bin_folders[] = { "5-7", "7-10", "10-20", "20-50" };
        const int pt_D0_mincut[] = { 2, 3, 5, 12 };
        for (int i = 0; i < n_bins; i++) {
            cout << "in pt bin" << i << endl;
            int pt_min = pt_bins[i];
            int pt_max = pt_bins[i+1];
            int pt_D0_min = pt_D0_mincut[i];

            // define pt related variables
            TString ptbin = TString::Format("%d #leq #it{p}_{T}^{ch. jet} < %d GeV/#it{c}, #font[122]{|}#it{#eta}_{jet}#font[122]{|} #leq 0.5", pt_min, pt_max);
            TString ptD = TString::Format("%d #leq #it{p}_{T}^{D^{0}} < %d GeV/#it{c}, #font[122]{|}#it{y}_{D^{0}}#font[122]{|} #leq 0.8", pt_D0_min, pt_max);
            
            std::string pdf_outdir = "plots/" + output_dirs[plot_case] + "/" + pt_bin_folders[i] + "/";

            // // make a canvas for each pt range
            // TCanvas* c = new TCanvas();
            // ProcessCanvas(c);
            // c->cd();
            // gPad->SetLogx();
            // // gPad->SetLogy();

            // TLegend* l; // = new TLegend(0.17, 0.65, 0.5, 0.85);

            // double maxy = 0;


            // Open histograms


            // l = new TLegend(0.1797168,0.400741,0.4562155,0.8885185,""); //(0.17, 0.4, 0.5, 0.53);
            // l->SetTextSize(0.045);
            // // TLegend *leg = new TLegend(0.1797168,0.5390741,0.4562155,0.8885185,"");
            // l->AddEntry("NULL","PYTHIA 8 Monash 2013","h");
            // l->AddEntry("NULL","pp, #sqrt{#it{s}} = 13 TeV","h");
            // l->AddEntry("NULL","D^{0} #rightarrow K^{#minus} #pi^{+} and charge conj.","h");
            // l->AddEntry("NULL","in charged jets, anti-#it{k}_{T}, #it{R} = 0.4","h");
            // l->AddEntry("NULL",ptbin,"h");
            // l->AddEntry("NULL",ptD,"h");
            // l->SetTextSize(0.037);
            // l->SetBorderSize(0);
            // // l->Draw("same");

            TH1* hD0KpiNjets = (TH1*) f->Get(hD0KpiNjets_name.c_str());
            

            //-------------------------------------------------//
            // find D0 reconstruction through charm
            std::vector<THnSparse*> hsparsejet_c;
            std::vector<THnSparse*> hsparsejet_g;
            std::vector<THnSparse*> hsparsejet_l;
            std::vector<THnSparse*> hsparsejet_i;
            
            for (int iobs=0; iobs<numjetaxes; iobs++) {
                hsparsejet_c.push_back( (THnSparse*) f->Get(hc_names[iobs].c_str()) );
                hsparsejet_g.push_back( (THnSparse*) f->Get(hg_names[iobs].c_str()) );
                hsparsejet_l.push_back( (THnSparse*) f->Get(hl_names[iobs].c_str()) );
                hsparsejet_i.push_back( (THnSparse*) f->Get(hi_names[iobs].c_str()) );
            }
            // THnSparse* hsparsejet_test = (THnSparse*) f->Get(htest_name.c_str());

            // testing - look at # jets before cuts
            // cout << "numDtaggedjets from hist before cuts " << hsparsejet_c_jetlevel->Projection(0)->GetEntries() << endl;

            // for THnSparse: make clone to work with, make cuts, get projection
            std::vector<THnSparse*> hsparsejet_c_clones;
            std::vector<THnSparse*> hsparsejet_g_clones;
            std::vector<THnSparse*> hsparsejet_l_clones;
            std::vector<THnSparse*> hsparsejet_i_clones;

            for (int iobs=0; iobs<numjetaxes; iobs++) {
                hsparsejet_c_clones.push_back( (THnSparse *) hsparsejet_c[iobs]->Clone(Form("hsparsejet_c_%s_clone", jetaxis_names[iobs].c_str())) );
                hsparsejet_g_clones.push_back( (THnSparse *) hsparsejet_g[iobs]->Clone(Form("hsparsejet_g_%s_clone", jetaxis_names[iobs].c_str())) );
                hsparsejet_l_clones.push_back( (THnSparse *) hsparsejet_l[iobs]->Clone(Form("hsparsejet_l_%s_clone", jetaxis_names[iobs].c_str())) );
                hsparsejet_i_clones.push_back( (THnSparse *) hsparsejet_i[iobs]->Clone(Form("hsparsejet_i_%s_clone", jetaxis_names[iobs].c_str())) );
            }

            // THnSparse* hsparsejet_test_clone = (THnSparse *) hsparsejet_test->Clone(Form("hsparsejet_c_%s_clone", jetaxis_names[0].c_str()));
            // hsparsejet_test_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max);
            // hsparsejet_test_clone->GetAxis(1)->SetRangeUser(5., pt_max);
            // hsparsejet_test_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8);

            // get jet pT range
            for (int iobs=0; iobs<numjetaxes; iobs++) {
                hsparsejet_c_clones[iobs]->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
                hsparsejet_g_clones[iobs]->GetAxis(0)->SetRangeUser(pt_min, pt_max);
                hsparsejet_l_clones[iobs]->GetAxis(0)->SetRangeUser(pt_min, pt_max);
                hsparsejet_i_clones[iobs]->GetAxis(0)->SetRangeUser(pt_min, pt_max);
            }
           
            if (plot_case == 0) { //D0 reconstruction
                for (int iobs=0; iobs<numjetaxes; iobs++) {
                    hsparsejet_c_clones[iobs]->GetAxis(1)->SetRangeUser(pt_D0_min, pt_max); // apply cut on Dmeson pt
                    hsparsejet_g_clones[iobs]->GetAxis(1)->SetRangeUser(pt_D0_min, pt_max);
                    hsparsejet_l_clones[iobs]->GetAxis(1)->SetRangeUser(pt_D0_min, pt_max);
                    hsparsejet_i_clones[iobs]->GetAxis(1)->SetRangeUser(pt_D0_min, pt_max);
                    
                    hsparsejet_c_clones[iobs]->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
                    hsparsejet_g_clones[iobs]->GetAxis(2)->SetRangeUser(-0.8, 0.8);
                    hsparsejet_l_clones[iobs]->GetAxis(2)->SetRangeUser(-0.8, 0.8);
                    hsparsejet_i_clones[iobs]->GetAxis(2)->SetRangeUser(-0.8, 0.8);
                }
            }

            // Project onto observable axis
            std::vector<TH1D*> hc_projs;
            std::vector<TH1D*> hg_projs;
            std::vector<TH1D*> hl_projs;
            std::vector<TH1D*> hi_projs;
            for (int iobs=0; iobs<numjetaxes; iobs++) {
                // TH1D *hD0_projs = hsparsejet_c_clones[iobs]->Projection(3);
                hc_projs.push_back( hsparsejet_c_clones[iobs]->Projection(3) );
                hg_projs.push_back( hsparsejet_g_clones[iobs]->Projection(3) );
                hl_projs.push_back( hsparsejet_l_clones[iobs]->Projection(3) );
                hi_projs.push_back( hsparsejet_i_clones[iobs]->Projection(3) );
            }

            // TH1D *htest_proj = hsparsejet_test_clone->Projection(3);
            // std::string hname = htest_proj->GetName();
            // hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            // htest_proj->SetNameTitle(hname.c_str(), hname.c_str());
            // cout << "hello " << ", " << hname.c_str() << endl;
              
            // Set to appropriate name
            std::string hname;
            for (int iobs=0; iobs<numjetaxes; iobs++) {
                hname = hc_projs[iobs]->GetName();
                hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
                hc_projs[iobs]->SetNameTitle(hname.c_str(), hname.c_str());
                hname = hg_projs[iobs]->GetName();
                hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
                hg_projs[iobs]->SetNameTitle(hname.c_str(), hname.c_str());
                hname = hl_projs[iobs]->GetName();
                hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
                hl_projs[iobs]->SetNameTitle(hname.c_str(), hname.c_str());
                hname = hi_projs[iobs]->GetName();
                hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
                hi_projs[iobs]->SetNameTitle(hname.c_str(), hname.c_str());
                cout << "iter " << iobs << ", " << hname.c_str() << endl;
            }


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
            std::vector<TH1*> hc;
            std::vector<TH1*> hg;
            std::vector<TH1*> hl;
            std::vector<TH1*> hi;

            for (int iobs=0; iobs<numjetaxes; iobs++) {
                std::string hc_newname = hc_projs[iobs]->GetName();
                std::string hg_newname = hg_projs[iobs]->GetName();
                std::string hl_newname = hl_projs[iobs]->GetName();
                std::string hi_newname = hi_projs[iobs]->GetName();

                hc.push_back( (TH1*) hc_projs[iobs]->Clone( hc_newname.c_str() )); //hc_newname + "rebin").c_str() );
                hg.push_back( (TH1*) hg_projs[iobs]->Clone( hg_newname.c_str() )); //hg_newname + "rebin").c_str() );
                hl.push_back( (TH1*) hl_projs[iobs]->Clone( hl_newname.c_str() )); //hl_newname + "rebin").c_str() );
                hi.push_back( (TH1*) hi_projs[iobs]->Clone( hi_newname.c_str() )); //hi_newname + "rebin").c_str() );
            }
            // TH1* hD0 = (TH1*) hD0_proj->Rebin(n_obs_bins, (hname + "rebin").c_str(), obs_bins);
            cout << "Rebin done" << endl;


            // Find normalization factor
            std::vector<double> numjets_charm;
            std::vector<double> numjets_gluon;
            std::vector<double> numjets_light;
            std::vector<double> numjets_inclusive;

            for (int iobs=0; iobs<numjetaxes; iobs++) {
                numjets_charm.push_back( hc[iobs]->Integral() );
                numjets_gluon.push_back( hg[iobs]->Integral() );
                numjets_light.push_back( hl[iobs]->Integral() );
                numjets_inclusive.push_back( hi[iobs]->Integral() );
            }


            // Set normalization
            /*
            for (int iobs=0; iobs<numjetaxes; iobs++) {
                cout << "num charm jets iter " << iobs << "; " << numjets_charm[iobs] << endl;
                hc[iobs]->Scale(1/numjets_charm[iobs], "width");
                hg[iobs]->Scale(1/numjets_gluon[iobs], "width");
                hl[iobs]->Scale(1/numjets_light[iobs], "width");
                hi[iobs]->Scale(1/numjets_inclusive[iobs], "width");
            }
            */




            //Format color and style
            int markercolor1 = kRed; //charm
            int markerstyle1 = kFullCircle;
            int markercolor2 = kViolet+2; //gluon
            int markerstyle2 = 33;
            int markercolor3 = kGreen+2; //light
            int markerstyle3 = 21;
            int markercolor4 = kMagenta-7; //inclusive
            int markerstyle4 = 29;
            // label1 = "charm-init jets";
            // label2 = "gluon-init jets";
            // TString label3 = "light-init jets";


            // if (plot_case == 0) {
            //     label1 = "D^{0}-tagged, c-init jets";
            // } else if (plot_case == 1 or plot_case == 3) {
            //     label1 = "D^{0}-tagged feeddown, b-init jets";
            // } else if (plot_case == 2) {
            //     label1 = "#phi-tagged, inclusive jets";
            // }

            // htest_proj->Draw();
            // std::string fnametest = outdir + "test.pdf";
            // const char* fnametestc = fnametest.c_str();
            // c->SaveAs(fnametestc);

            // Format histograms for plotting (this order needed to keep legend in order and graphs lookin good)
            for (int iobs=0; iobs<numjetaxes; iobs++) {

                //make a canvas
                TCanvas* c = new TCanvas();
                ProcessCanvas(c);
                c->cd();

                hc[iobs]->GetXaxis()->SetTitle("#DeltaR");
                hc[iobs]->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #frac{d#it{N}}{d#DeltaR}");
                cout << "about to format D0" << endl;


                TLegend* l = new TLegend(0.5597168,0.400741,0.8562155,0.8885185,""); //(0.17, 0.4, 0.5, 0.53);
                
                if (plot_case == 0) {
                    addLegendInfo(l, ptbin, ptD, true);
                    FormatHist(l, hc[iobs], "D0-tagged jets", markercolor1, markerstyle1);    
                
                } else if (plot_case > 0) { //no D0 reconstruction
                    addLegendInfo(l, ptbin, ptD, false);
                    FormatHist(l, hc[iobs], "Charm-initiated jets", markercolor1, markerstyle1);

                    if (plot_case == 1) l->AddEntry("NULL","          charm decays off","h");
                    else if (plot_case == 2) l->AddEntry("NULL","          charm decays on","h");

                    FormatHist(l, hg[iobs], "Gluon-initiated jets", markercolor2, markerstyle2);
                    FormatHist(l, hl[iobs], "Light-initiated jets", markercolor3, markerstyle3);
                    FormatHist(l, hi[iobs], "Inclusive jets", markercolor4, markerstyle4);

                    // find max to plot nicely
                    // cout << "max " << hc[iobs]->GetMaximum() << endl;
                    // cout << "max " << hg[iobs]->GetMaximum() << endl;
                    // cout << "max " << hl[iobs]->GetMaximum() << endl;
                    // cout << "max " << hi[iobs]->GetMaximum() << endl;
                    double arr_of_maxes[] = { hc[iobs]->GetMaximum(), hg[iobs]->GetMaximum(), hl[iobs]->GetMaximum(), hi[iobs]->GetMaximum() };
                    double &maxy = *std::max_element(arr_of_maxes, arr_of_maxes+4); //bc there are 4 elements in arr_of_maxes
                    // maxy = arr_of_maxes[maxindex]
                    cout << "the max is " << maxy << endl;
                    maxy *= 1.2;
                    hc[iobs]->SetMaximum(maxy);
                }


                hc[iobs]->Draw("L same");
                if (plot_case > 0) { //no D0 reconstruction
                    hg[iobs]->Draw("L same");
                    hl[iobs]->Draw("L same");
                    hi[iobs]->Draw("L same");
                }

                // draw legend
                l->Draw("same");

                std::string fname = pdf_outdir + "QG_comp_" + jetaxis_names[iobs] + "_pt" + std::to_string(pt_min) + '-' + std::to_string(pt_max) + "_R" + jetR + add_name + ".pdf";
                const char* fnamec = fname.c_str();
                c->SaveAs(fnamec);
                delete c;
                delete l;

                // Write rebinned histograms to root file
                f_out->cd();
                hc[iobs]->Write();
                hg[iobs]->Write();
                hl[iobs]->Write();
                hi[iobs]->Write();

            } //obs bins loop
            
             
        } // pT bins loop
    } // jetR loop

    f->Close();
    delete f;

    return;
}
