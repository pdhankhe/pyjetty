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


void make_qg_plots_inclusive() {

//    gROOT->SetBatch(); //prevents plots from showing up
    gStyle->SetOptStat(0);
    SetStyle();
    Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
    Double_t marker_size = 1.5;
    Double_t colors[16] = {kRed, kGreen+2, kBlue, kRed+1, kGreen+1, kBlue+1, kRed+2, kGreen+2, kBlue+2, kRed+3, kGreen+3, kBlue+3, kOrange+1, kViolet+1, kYellow+1, kCyan+1};


    // File names
    const char infile_inclusive_original[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/23209925/AnalysisResultsFinal.root";
    const char infile_inclusive_only_13TeV[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/23135855/AnalysisResultsFinal.root"; //this is using thnsparse
    const char infile_inclusive_only_5TeV[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/23209945/AnalysisResultsFinal.root"; //this is using thnsparse
    const char infile_inclusive_original_idktrackcuts[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/14680822/AnalysisResultsFinal.root"; //this is using thnsparse
    


    // int plot_case:
    // 0 = D0, D* decays turned off
    // 1 = D0 feeddown
    // 2 = phi meson
    // 3 = D0 feeddown, pT*RL

    //CONTOL VARIABLES HERE
    int plot_case = 3;

    TFile* f1;
    TFile* f2;
    TFile* f3;
    TFile* f4;
    TString label1 = "";
    TString label2 = "";
    TString label3 = "";
    TString label4 = "";

    f1 = new TFile(infile_inclusive_original, "READ");
    f2 = new TFile(infile_inclusive_only_13TeV, "READ");
    f3 = new TFile(infile_inclusive_only_5TeV, "READ");
    f4 = new TFile(infile_inclusive_original_idktrackcuts, "READ");

    std::string add_name;
    add_name = "_inclusivecomparison";
    cout << "output name will be " << add_name << endl;

    // Output directory
    std::string outdir = "plots/";//"plots/test/"; 
    // Output file for binned results
    std::string outfile = outdir + "AnalysisResultsFinal" + add_name + ".root"; 
    TFile* f_out = new TFile(outfile.c_str(), "RECREATE");

    //Format color and style
    int markercolor1 = kBlue-4;//inclusive, old
    int markerstyle1 = kFullCircle;
    int markercolor2 = kRed-3; //inclusive, new, 13 TeV
    int markerstyle2 = kFullCircle;
    int markercolor3 = kGreen-3;  //kViolet-1; //inclusive, new, 5 TeV
    int markerstyle3 = kFullCircle;
    int markercolor4 = kTeal; //kViolet-1; //inclusive, old, unknown track cuts
    int markerstyle4 = kFullCircle;
    label1 = "Inclusive jets, old, 13 TeV";
    label2 = "Inclusive jets, new, 13 TeV";
    label3 = "Inclusive jets, new, 5 TeV";
    label4 = "Inclusive jets, old, 13 TeV, unknown track cuts";


    // Jet R value
    std::string jetR_list[] = { "0.4" };
    for (std::string jetR : jetR_list) {

        // Names of histograms in the file 
        const std::string hi_name = "h_EEC_JetPt_inclusive_R" + jetR;
        const std::string hi_jet_name = "h_JetPt_inclusive_R" + jetR + "_jetlevel";

        const std::string hD0KpiNjets_name = "hD0KpiNjets"; //TODO later: = "hD0KpiNehD0KpiNjetsvents" for run 16729583

        // make a canvas for jet pt spectrum
        TCanvas* c_jetpt = new TCanvas();
        ProcessCanvas(c_jetpt);
        c_jetpt->cd();
        gPad->SetLogy();

        THnSparse* hsparsejet_i1_jet = (THnSparse*) f1->Get(hi_jet_name.c_str());
        THnSparse* hsparsejet_i2_jet = (THnSparse*) f2->Get(hi_jet_name.c_str());
        THnSparse* hsparsejet_i3_jet = (THnSparse*) f3->Get(hi_jet_name.c_str());
        TH1* hi4_jet = (TH1*) f4->Get(hi_jet_name.c_str());

        TH1D *hi1_jet = hsparsejet_i1_jet->Projection(0);
        TH1D *hi2_jet = hsparsejet_i2_jet->Projection(0);
        TH1D *hi3_jet = hsparsejet_i3_jet->Projection(0);

        // FormatHist(l, hi1_jet, label1, markercolor1, markerstyle1, 0.5);
        // FormatHist(l, hi2_jet, label2, markercolor2, markerstyle2, 0.5);
        // FormatHist(l, hi3_jet, label3, markercolor3, markerstyle3, 0.6);
        // FormatHist(l, hi4_jet, label4, markercolor4, markerstyle4, 1.);
        hi1_jet->SetMarkerColor(markercolor1); hi1_jet->SetLineColor(markercolor1); hi1_jet->SetLineStyle(1); hi1_jet->SetLineWidth(3);
        hi2_jet->SetMarkerColor(markercolor2); hi2_jet->SetLineColor(markercolor2); hi2_jet->SetLineStyle(1); hi2_jet->SetLineWidth(3);
        hi3_jet->SetMarkerColor(markercolor3); hi3_jet->SetLineColor(markercolor3); hi3_jet->SetLineStyle(1); hi3_jet->SetLineWidth(3);
        hi4_jet->SetMarkerColor(markercolor4); hi4_jet->SetLineColor(markercolor4); hi4_jet->SetLineStyle(1); hi4_jet->SetLineWidth(3);

        hi1_jet->Draw("C same");
        hi2_jet->Draw("C same");
        hi3_jet->Draw("C same");
        hi4_jet->Draw("C same");

        std::string fnamejetpt = outdir + "QG_comp_jetptspectrum" + "_R" + jetR + add_name + ".pdf";
        const char* fnamejetptc = fnamejetpt.c_str();
        c_jetpt->SaveAs(fnamejetptc);
        delete c_jetpt;




        const int pt_bins[] = { 20, 40, 60, 80 }; //{ 10, 15, 30 }; //{ 10, 20, 40 }; // CHANGE HERE!!
        const int n_bins = 3; //2; 
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


            l = new TLegend(0.1797168,0.500741,0.4562155,0.8885185,""); //(0.17, 0.4, 0.5, 0.53);
            l->SetTextSize(0.045);
            // TLegend *leg = new TLegend(0.1797168,0.5390741,0.4562155,0.8885185,"");
            l->AddEntry("NULL","PYTHIA 8 Monash 2013","h");
            l->AddEntry("NULL","pp, #sqrt{#it{s}} = 13 TeV or 5.02 TeV","h");
            // l->AddEntry("NULL","D^{0} #rightarrow K^{#minus} #pi^{+} and charge conj.","h");
            l->AddEntry("NULL","in charged jets, anti-#it{k}_{T}, #it{R} = 0.4","h");
            l->AddEntry("NULL",ptbin,"h");
            // l->AddEntry("NULL",ptD,"h");
            l->SetTextSize(0.037);
            l->SetBorderSize(0);
            // l->Draw("same");

            // TH1* hD0KpiNjets = (TH1*) f->Get(hD0KpiNjets_name.c_str());
            

            //-------------------------------------------------//
            // find D0 reconstruction through charm
            THnSparse* hsparsejet_i1;
            THnSparse* hsparsejet_i1_jetlevel;
            THnSparse* hsparsejet_i2;
            THnSparse* hsparsejet_i2_jetlevel;
            THnSparse* hsparsejet_i3;
            THnSparse* hsparsejet_i3_jetlevel;

            hsparsejet_i1 = (THnSparse*) f1->Get(hi_name.c_str());
            hsparsejet_i1_jetlevel = (THnSparse*) f1->Get(hi_jet_name.c_str());
            hsparsejet_i2 = (THnSparse*) f2->Get(hi_name.c_str());
            hsparsejet_i2_jetlevel = (THnSparse*) f2->Get(hi_jet_name.c_str());
            hsparsejet_i3 = (THnSparse*) f3->Get(hi_name.c_str());
            hsparsejet_i3_jetlevel = (THnSparse*) f3->Get(hi_jet_name.c_str());
            TH2* hi4_2D = (TH2*) f4->Get(hi_name.c_str());
            TH1* hi4_1D_jet = (TH1*) f4->Get(hi_jet_name.c_str());

            // testing - look at # jets before cuts
            cout << "numDtaggedjets from hist before cuts " << hsparsejet_i1_jetlevel->Projection(0)->GetEntries() << endl;

            // for THnSparse: make clone to work with, make cuts, get projection
            THnSparse *hsparsejet_i1_clone = (THnSparse *) hsparsejet_i1->Clone("hsparsejet_i1_clone");
            THnSparse *hsparsejet_i1_jetlevel_clone = (THnSparse *) hsparsejet_i1_jetlevel->Clone("hsparsejet_i1_jetlevel_clone");
            THnSparse *hsparsejet_i2_clone = (THnSparse *) hsparsejet_i2->Clone("hsparsejet_i2_clone");
            THnSparse *hsparsejet_i2_jetlevel_clone = (THnSparse *) hsparsejet_i2_jetlevel->Clone("hsparsejet_i2_jetlevel_clone");
            THnSparse *hsparsejet_i3_clone = (THnSparse *) hsparsejet_i3->Clone("hsparsejet_i3_clone");
            THnSparse *hsparsejet_i3_jetlevel_clone = (THnSparse *) hsparsejet_i3_jetlevel->Clone("hsparsejet_i3_jetlevel_clone");

            // get jet pT range
            hsparsejet_i1_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
            hsparsejet_i1_jetlevel_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
            hsparsejet_i2_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); 
            hsparsejet_i2_jetlevel_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max);
            hsparsejet_i3_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); 
            hsparsejet_i3_jetlevel_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max);
            hi4_2D->GetXaxis()->SetRangeUser(pt_min, pt_max);
            hi4_1D_jet->GetXaxis()->SetRangeUser(pt_min, pt_max);


            // Project onto observable axis
            TH1D *hi1_proj = hsparsejet_i1_clone->Projection(3); //CALL THESE TH1*????
            TH1D *hi1_jet_proj = hsparsejet_i1_jetlevel_clone->Projection(0);
            TH1D *hi2_proj = hsparsejet_i2_clone->Projection(3);
            TH1D *hi2_jet_proj = hsparsejet_i2_jetlevel_clone->Projection(0);
            TH1D *hi3_proj = hsparsejet_i3_clone->Projection(3);
            TH1D *hi3_jet_proj = hsparsejet_i3_jetlevel_clone->Projection(0);
            TH1D *hi4_proj = hi4_2D->ProjectionY(); //(TH1*) hi2D->ProjectionY();

            // Set to appropriate name
            std::string hname = hi1_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hi1_proj->SetNameTitle(hname.c_str(), hname.c_str());
            hname = hi2_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hi2_proj->SetNameTitle(hname.c_str(), hname.c_str());
            hname = hi3_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hi3_proj->SetNameTitle(hname.c_str(), hname.c_str());
            hname = hi4_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hi4_proj->SetNameTitle(hname.c_str(), hname.c_str());

            // Find normalization factor
            double numjets_incl1 = hi1_jet_proj->Integral();
            double numjets_incl2 = hi2_jet_proj->Integral();
            double numjets_incl3 = hi3_jet_proj->Integral();
            double numjets_incl4 = hi4_1D_jet->Integral();

            // Set normalization
            // hD0_proj->Scale(1/numjets_charm, "width");


            // Rebin
            int n_obs_bins = 13; //-1;
            double obs_bins[14] = {0.0001, 0.00020892961308540387, 0.0004365158322401661, 0.0009120108393559096, 
                    0.0019054607179632482, 0.005754399373371567, 0.017378008287493762, 0.03630780547701014, 
                    0.07585775750291836, 0.15848931924611143, 0.3311311214825911, 0.47863009232263853, 0.6918309709189363, 1.0};

            cout << "About to rebin" << endl;

            // TH1* hi1 = (TH1*) hi1_proj->Rebin(n_obs_bins, (hname + "rebin").c_str(), obs_bins);
            TH1* hi1 = (TH1*) hi1_proj->Clone( hi1_proj->GetName() );
            TH1* hi2 = (TH1*) hi2_proj->Clone( hi2_proj->GetName() );
            TH1* hi3 = (TH1*) hi3_proj->Clone( hi3_proj->GetName() );
            TH1* hi4 = (TH1*) hi4_proj->Clone( hi4_proj->GetName() );

            cout << "Rebin done" << endl;


            hi1->Scale(1/numjets_incl1, "width");
            hi2->Scale(1/numjets_incl2, "width");
            hi3->Scale(1/numjets_incl3, "width");
            hi4->Scale(1/numjets_incl4, "width");


            

            // Format histograms for plotting (this order needed to keep legend in order and graphs lookin good)
            hi4->GetXaxis()->SetTitle("#it{R}_{L}");
            hi4->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
            cout << "about to format D0" << endl;
            FormatHist(l, hi4, label4, markercolor4, markerstyle4, 1.);
            FormatHist(l, hi1, label1, markercolor1, markerstyle1, 0.5);
            FormatHist(l, hi2, label2, markercolor2, markerstyle2, 0.5);
            FormatHist(l, hi3, label3, markercolor3, markerstyle3, 0.6);
            hi4->Draw("L same");
            hi1->Draw("L same");
            hi2->Draw("L same");
            hi3->Draw("L same");
            
            
            double hi1_top_binpos = findTopOfCurve(hi1);
            drawVertLine(hi1->GetBinCenter(hi1_top_binpos), 0, hi1->GetBinContent(hi1_top_binpos), markercolor1, 1)->Draw();
            double hi2_top_binpos = findTopOfCurve(hi2);
            drawVertLine(hi2->GetBinCenter(hi2_top_binpos), 0, hi2->GetBinContent(hi2_top_binpos), markercolor2, 1)->Draw();
            double hi3_top_binpos = findTopOfCurve(hi3);
            drawVertLine(hi3->GetBinCenter(hi3_top_binpos), 0, hi3->GetBinContent(hi3_top_binpos), markercolor3, 1)->Draw();
            double hi4_top_binpos = findTopOfCurve(hi4);
            drawVertLine(hi4->GetBinCenter(hi4_top_binpos), 0, hi4->GetBinContent(hi4_top_binpos), markercolor3, 1)->Draw();
            

            
            // vector<double> fullwidth_vec = findWidthOfCurve(hD0,  hD0_top_binpos);
            // drawHoriLine(fullwidth_vec[1], fullwidth_vec[2], fullwidth_vec[0], kMagenta+3, 1)->Draw();
            
            

            // draw legend
            l->Draw("same");



            std::string fname = outdir + "QG_comp_pt" + std::to_string(pt_min) + '-' + std::to_string(pt_max) + "_R" + jetR + add_name + ".pdf";
            const char* fnamec = fname.c_str();
            c->SaveAs(fnamec);
            delete c;

            // if (pt_min == 10) { //} && grooming == "") {
                // Write rebinned histograms to root file
            f_out->cd();
            hi1->Write();
            hi2->Write();
            hi3->Write();
            hi4->Write();
            // }
             
        } // pT bins loop
    } // jetR loop

    f1->Close();
    delete f1;
    f2->Close();
    delete f2;
    f3->Close();
    delete f3;
    f4->Close();
    delete f4;

    return;
}
