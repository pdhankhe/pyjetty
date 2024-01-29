// ROOT macro to make quark-gluon jet plots
// this one specifically for plotting&saving charm-/gluon-/light-init jets from THnSparse
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


void make_qg_plots_cgl() {

//    gROOT->SetBatch(); //prevents plots from showing up
    gStyle->SetOptStat(0);
    SetStyle();
    Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
    Double_t marker_size = 1.5;
    Double_t colors[16] = {kRed, kGreen+2, kBlue, kRed+1, kGreen+1, kBlue+1, kRed+2, kGreen+2, kBlue+2, kRed+3, kGreen+3, kBlue+3, kOrange+1, kViolet+1, kYellow+1, kCyan+1};


    // File containing quark vs gluon histograms

     //FOR WHEN WEIGHTED/UNWEIGHTED IN SAME FILE
    const char infile_ptrl_charmdecaysON[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/20780159/AnalysisResultsFinal.root";
    const char infile_ptrl_charmdecaysOFF[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/20780184/AnalysisResultsFinal.root"; 
    

    //CONTOL VARIABLES HERE
    bool charmdecays = true; // true = charm decays on, false = charm decays off
    bool ptrl = true; // true = Plot the pT*RL

    TFile* f;
    // std::vector<TFile*> files;
    std::string add_name;

    TString label1 = "";
    TString label2 = "";  

    if (ptrl) {
        if (charmdecays) {
            f = new TFile(infile_ptrl_charmdecaysON, "READ");
            add_name = "_cgl_ptrl_charmdecaysON.pdf";
        } else {
            f = new TFile(infile_ptrl_charmdecaysOFF, "READ");
            add_name = "_cgl_ptrl_charmdecaysOFF.pdf";
        }
    } else {
        cout << "This didn't work because no files are defined" << endl;
    }

    cout << "output name will be " << add_name << endl;

    // Output directory
    std::string outdir = "plots/final/cgl/ptrl/";//"plots/test/";
    // Output file for binned results
    std::string outfile = outdir + "AnalysisResultsFinal" + add_name + ".root";
    TFile* f_out = new TFile(outfile.c_str(), "RECREATE");

    // Jet r value
    std::string jetR_list[] = { "0.4" };
    for (std::string jetR : jetR_list) {

        // Names of histograms in the file (quark, charm, gluon)
        const std::string hc_name = "h_EEC_JetPt_charm_R" + jetR;
        const std::string hl_name = "h_EEC_JetPt_light_R" + jetR;
        const std::string hg_name = "h_EEC_JetPt_gluon_R" + jetR;
        const std::string hc_jet_name = "h_JetPt_charm_R" + jetR + "_jetlevel";
        const std::string hl_jet_name = "h_JetPt_light_R" + jetR + "_jetlevel";
        const std::string hg_jet_name = "h_JetPt_gluon_R" + jetR + "_jetlevel";

        const std::string hD0KpiNjets_name = "hD0KpiNjets"; //TODO later: = "hD0KpiNehD0KpiNjetsvents" for run 16729583


        //const int pt_bins[] = { 10, 20, 40, 60, 80, 100, 150 };
        const int pt_bins[] = { 10, 15, 30 }; //{ 10, 20, 40 };
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
            

            //-------------------------------------------------//
            // find D0 reconstruction through charm
            THnSparse* hsparsejet_c = (THnSparse*) f->Get(hc_name.c_str());
            THnSparse* hsparsejet_c_jetlevel = (THnSparse*) f->Get(hc_jet_name.c_str());
            THnSparse* hsparsejet_g = (THnSparse*) f->Get(hg_name.c_str());
            THnSparse* hsparsejet_g_jetlevel = (THnSparse*) f->Get(hg_jet_name.c_str());
            THnSparse* hsparsejet_l = (THnSparse*) f->Get(hl_name.c_str());
            THnSparse* hsparsejet_l_jetlevel = (THnSparse*) f->Get(hl_jet_name.c_str());


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

            cout << "checkpoint 1" << endl;

            // get jet pT range - no D0 reconstruction, so don't make D0 cuts
            hsparsejet_c_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
            hsparsejet_c_jetlevel_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
            hsparsejet_g_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
            hsparsejet_g_jetlevel_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
            hsparsejet_l_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
            hsparsejet_l_jetlevel_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
            
            cout << "checkpoint 2" << endl;


            // Project onto observable axis
            TH1D *hc_proj = hsparsejet_c_clone->Projection(3); //CALL THESE TH1*????
            TH1D *hc1D_jet = hsparsejet_c_jetlevel_clone->Projection(0);
            TH1D *hg_proj = hsparsejet_g_clone->Projection(3);
            TH1D *hg1D_jet = hsparsejet_g_jetlevel_clone->Projection(0);
            TH1D *hl_proj = hsparsejet_l_clone->Projection(3);
            TH1D *hl1D_jet = hsparsejet_l_jetlevel_clone->Projection(0);
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

            TH1D* hc = (TH1D*) hc_proj->Clone((hc_name + "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max)).c_str());
            TH1D* hl = (TH1D*) hl_proj->Clone((hl_name + "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max)).c_str());
            TH1D* hg = (TH1D*) hg_proj->Clone((hg_name + "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max)).c_str());
            

            // Find normalization factor
            double numjets_charm = hc1D_jet->Integral();
            double numjets_gluon = hg1D_jet->Integral();
            double numjets_light = hl1D_jet->Integral();

            // Set normalization
            hc->Scale(1/numjets_charm, "width");
            hg->Scale(1/numjets_gluon, "width");
            hl->Scale(1/numjets_light, "width");


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
            label1 = "charm-init jets";
            label2 = "gluon-init jets";
            TString label3 = "light-init jets";

            // Format histograms for plotting (this order needed to keep legend in order and graphs lookin good)
            hc->GetXaxis()->SetTitle("#it{p}_{T}#it{R}_{L}");
            hc->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
            hg->GetXaxis()->SetTitle("#it{p}_{T}#it{R}_{L}");
            hg->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
            hl->GetXaxis()->SetTitle("#it{p}_{T}#it{R}_{L}");
            hl->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
            
            cout << "about to format plots" << endl;
            FormatHist(l, hc, label1, markercolor1, markerstyle1);
            if (charmdecays) {
                l->AddEntry("NULL","          D* decays on","h");
            } else {
                l->AddEntry("NULL","          D* decays off","h");
            }
            FormatHist(l, hg, label2, markercolor2, markerstyle2);
            FormatHist(l, hl, label3, markercolor3, markerstyle3);
            
            if (charmdecays == false) {
                hl->Draw("L same");
            }
            hc->Draw("L same");
            hg->Draw("L same");
            if (charmdecays) {
                hl->Draw("L same");
            }
            
            
            
            // hDstar->Draw("L same");
            

            double hc_top_binpos = findTopOfCurve(hc, ptrl);
            cout << "hc top bincenter" << hc->GetBinCenter(hc_top_binpos) << "hc top binpos" << hc->GetBinContent(hc_top_binpos) << endl;
            drawVertLine(hc->GetBinCenter(hc_top_binpos), 0, hc->GetBinContent(hc_top_binpos), markercolor1, 1)->Draw();
            double hg_top_binpos = findTopOfCurve(hg, ptrl);
            cout << "hg top bincenter" << hg->GetBinCenter(hg_top_binpos) << "hg top binpos" << hg->GetBinContent(hg_top_binpos) << endl;
            drawVertLine(hg->GetBinCenter(hg_top_binpos), 0, hg->GetBinContent(hg_top_binpos), markercolor2, 1)->Draw();
            double hl_top_binpos = findTopOfCurve(hl, ptrl);
            cout << "hl top bincenter" << hl->GetBinCenter(hl_top_binpos) << "hc top binpos" << hc->GetBinContent(hc_top_binpos) << endl;
            drawVertLine(hl->GetBinCenter(hl_top_binpos), 0, hl->GetBinContent(hl_top_binpos), markercolor3, 1)->Draw();


            
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



            std::string fname = outdir + "QG_comp_pt" + std::to_string(pt_min) + '-' + std::to_string(pt_max) + "_R" + jetR + add_name; //"_charmdecaysONcomparison.pdf"; // + "_normbytype.pdf"; //"_nonorm.pdf";
            const char* fnamec = fname.c_str();
            c->SaveAs(fnamec);
            delete c;

            f_out->cd();
            hc->Write();
            hg->Write();
            hl->Write();

            delete hc;
            delete hl;
            delete hg;
            delete hc1D_jet;
            delete hg1D_jet;
            delete hl1D_jet;
            
             
        } // pT bins loop
    } // jetR loop

    f->Close();
    delete f;


    return;
}
