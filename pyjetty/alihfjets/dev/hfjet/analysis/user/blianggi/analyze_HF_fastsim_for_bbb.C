// This code will take in HF Herwig fastsim (both truth and detector level)
// and find the bin-by-bin corrections


const int pt_bins[] = { 10, 15, 30 }; //, 100, 150 }; //{ 10, 20, 40 };
const int d0_pt_cuts[] = { 5, 5 }; //, 5, 5 };
const int n_bins = 2;

std::string output_add_name = "";

// Set bool for what to run
bool include_dstar = true;


// format canvas!!
void SetStyle(Bool_t graypalette=true) {
    cout << "Setting style!" << endl;
  
    gStyle->Reset("Plain");
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    // if(graypalette) gStyle->SetPalette(8,0);
    // else gStyle->SetPalette(1);
    gStyle->SetPalette(kRainbow); //kBird
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


// THnSparse axes are: jet pt, D0 pt, D0 y, D0 z (, RL)
void applyCuts(THnSparse *hsparse, int pt_min, int pt_max, int d0_pt_cut, bool d0cuts=false) {
    hsparse->GetAxis(0)->SetRangeUser(pt_min, pt_max);
    if (d0cuts) {
        hsparse->GetAxis(1)->SetRangeUser(d0_pt_cut, pt_max); //d0_pt_cuts[i], pt_max); // apply cut on Dmeson pt
        hsparse->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
    }
}

void FormatHist(TLegend *l, TH1 *hist, TString text, int markercolor=1, int markerstyle=8, double markeralpha=1.,
                std::string xtitle = "", std::string ytitle = "",
                double xtitlesize=0.06, double xlabelsize=0.05, double xoffset=1.0,
                double ytitlesize=0.06, double ylabelsize=0.05, double yoffset=1.05, double markersize=1.0) 
                // double xtitlesize=0.04, double xlabelsize=0.04, double xoffset=1.2,
                // double ytitlesize=0.04, double ylabelsize=0.04, double yoffset=1.0) 
{
    hist->SetLineColor(markercolor);
    hist->SetMarkerColorAlpha(markercolor, markeralpha);
    hist->SetMarkerStyle(markerstyle);
    hist->SetMarkerSize(markersize);
    l->AddEntry(hist, text, "pl");

    hist->GetXaxis()->SetTitle(xtitle.c_str());
    hist->GetYaxis()->SetTitle(ytitle.c_str());

	//gPad->SetTickx(); 
	//gPad->SetTicky(); 
	// h->SetLineWidth(2);
	hist->GetYaxis()->SetTitleOffset(yoffset); //(1.05); 
	hist->GetYaxis()->SetTitleSize(ytitlesize); //(0.042); //the axis number labels
	hist->GetYaxis()->SetLabelSize(ylabelsize); //(0.042);
	hist->GetYaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetTitleFont(42);

	hist->GetXaxis()->SetTitleOffset(xoffset);
	hist->GetXaxis()->SetTitleSize(xtitlesize); //(0.042);
	hist->GetXaxis()->SetLabelSize(xlabelsize); //(0.042);
    hist->GetXaxis()->SetLabelFont(42);
	hist->GetXaxis()->SetTitleFont(42);


    return;
}

TLine * drawHoriLine(double x1, double x2, double y1, int color, int linestyle=2){
    auto fhoriline = new TLine(x1, y1, x2, y1);
	fhoriline->SetLineWidth(1);
    fhoriline->SetLineColor(color);
    fhoriline->SetLineStyle(linestyle);
    return fhoriline;

}


void plot_pt_comparisons(THnSparse * h_truth_jet_pt_thnsparse, THnSparse * h_det_jet_pt_thnsparse, 
                         THnSparse * h_matched_truth_jet_pt_thnsparse, THnSparse * h_matched_det_jet_pt_thnsparse, 
                         std::string type_pt_name) {

    
    TCanvas * can_pt = new TCanvas();
    // Define pad heights (top bigger, bottom smaller)
    float padSplit = 0.3;  // fraction for bottom pad

    // Create top pad
    TPad *pad1 = new TPad("pad1", "pad1", 0, padSplit, 1, 1);
    pad1->SetBottomMargin(0); // remove bottom margin for top pad
    pad1->SetTicks(1,1);
    pad1->Draw();

    TLegend * leg_pt = new TLegend(0.6, 0.6, 0.8, 0.8);
    TLegend * leg_pt_ratio = new TLegend(0.6, 0.7, 0.8, 0.9);

    int proj_axis = 0; // if (type_pt_name == "jet_pt") 
    std::string xtitle = "p_{T, jet}";
    // if (type_pt_name == "D0_pt") {
    if (type_pt_name.find("D0_pt") != std::string::npos){
        proj_axis = 1;
        xtitle = "p_{T, D^{0}}";
    }
    
    // Get pt distributions
    TH1D * h_pt_truth_all = h_truth_jet_pt_thnsparse->Projection(proj_axis);
    TH1D * h_pt_det_all = h_det_jet_pt_thnsparse->Projection(proj_axis);
    TH1D * h_pt_matched_truth_all = h_matched_truth_jet_pt_thnsparse->Projection(proj_axis);
    TH1D * h_pt_matched_det_all = h_matched_det_jet_pt_thnsparse->Projection(proj_axis);

    // Format histograms
    FormatHist(leg_pt, h_pt_truth_all, "All truth", kBlue, kFullStar, 0.8, xtitle.c_str(), "Counts");
    FormatHist(leg_pt, h_pt_det_all, "All det", kBlue, kFullCrossX, 0.8, xtitle.c_str(), "Counts");
    FormatHist(leg_pt, h_pt_matched_truth_all, "Matched truth", kRed, kFullStar, 0.8, xtitle.c_str(), "Counts");
    FormatHist(leg_pt, h_pt_matched_det_all, "Matched det", kRed, kFullCrossX, 0.8, xtitle.c_str(), "Counts");
    h_pt_matched_det_all->SetMaximum(h_pt_truth_all->GetMaximum() * 1.2);

    // Get and format ratios
    TH1D * h_pt_all_ratio = (TH1D *) h_pt_det_all->Clone("h_pt_all_ratio");
    h_pt_all_ratio->Divide(h_pt_truth_all);
    TH1D * h_pt_matched_all_ratio = (TH1D *) h_pt_matched_det_all->Clone("h_pt_matched_all_ratio");
    h_pt_matched_all_ratio->Divide(h_pt_matched_truth_all);
    FormatHist(leg_pt_ratio, h_pt_all_ratio, "All", kBlue, kFullCross, 0.8, xtitle.c_str(), "Det / Truth");
    FormatHist(leg_pt_ratio, h_pt_matched_all_ratio, "Matched truth", kRed, kFullCross, 0.8, xtitle.c_str(), "Det / Truth");
    h_pt_all_ratio->SetMaximum(2);
    
    // Plot comparison of pt distributions
    pad1->cd();
    gPad->SetLogy();
    h_pt_matched_det_all->Draw(); // drawing this first to get the minimum right
    h_pt_det_all->Draw("SAME");
    h_pt_matched_truth_all->Draw("SAME");
    h_pt_truth_all->Draw("SAME");
    leg_pt->Draw();

    // Plot ratio of pt distributions
    can_pt->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, padSplit);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2); // leave room for x-axis labels
    pad2->SetTicks(1,1);
    pad2->Draw();

    pad2->cd();
    h_pt_all_ratio->Draw();
    h_pt_matched_all_ratio->Draw("SAME");
    leg_pt_ratio->Draw();


    // Save
    can_pt->SaveAs(Form("/software/users/blianggi/mypyjetty/storage/HF_EEC/plots/%s_distribution%s.pdf", type_pt_name.c_str(), output_add_name.c_str()));
}

void get_and_plot_bbb(TFile *file) {
    
    // Get jet pt histograms
    THnSparse * h_truth_jet_pt_thnsparse = (THnSparse *) file->Get("h_jet_pt_JetPt_Truth_R0.4_1.0Scaled");
    THnSparse * h_det_jet_pt_thnsparse = (THnSparse *) file->Get("h_jet_pt_JetPt_Det_R0.4_1.0Scaled");
    THnSparse * h_matched_truth_jet_pt_thnsparse = (THnSparse *) file->Get("h_matched_jet_pt_JetPt_Truth_R0.4_1.0Scaled");
    THnSparse * h_matched_det_jet_pt_thnsparse = (THnSparse *) file->Get("h_matched_jet_pt_JetPt_Det_R0.4_1.0Scaled");

    // Get EEC histograms
    THnSparse * h_truth_ENC_RL2_thnsparse = (THnSparse *) file->Get("h_jet_ENC_RL2_JetPt_Truth_R0.4_1.0Scaled");
    THnSparse * h_det_ENC_RL2_thnsparse = (THnSparse *) file->Get("h_jet_ENC_RL2_JetPt_Det_R0.4_1.0Scaled");
    THnSparse * h_matched_truth_ENC_RL2_thnsparse = (THnSparse *) file->Get("h_matched_jet_ENC_RL2_JetPt_Truth_R0.4_1.0Scaled");
    THnSparse * h_matched_det_ENC_RL2_thnsparse = (THnSparse *) file->Get("h_matched_jet_ENC_RL2_JetPt_Det_R0.4_1.0Scaled");

    // Print statements on number of entries
    cout << endl;
    cout << "  # of jets in truth     all: " << h_truth_jet_pt_thnsparse->GetEntries() << endl;
    cout << "  # of jets in det       all: " << h_det_jet_pt_thnsparse->GetEntries() << endl;
    cout << "  # of jets in truth matched: " << h_matched_truth_jet_pt_thnsparse->GetEntries() << endl;
    cout << "  # of jets in det   matched: " << h_matched_det_jet_pt_thnsparse->GetEntries() << endl;
    cout << endl;
    cout << "  # of pairs in truth     all: " << h_truth_ENC_RL2_thnsparse->GetEntries() << endl;
    cout << "  # of pairs in det       all: " << h_det_ENC_RL2_thnsparse->GetEntries() << endl;
    cout << "  # of pairs in truth matched: " << h_matched_truth_ENC_RL2_thnsparse->GetEntries() << endl;
    cout << "  # of pairs in det   matched: " << h_matched_det_ENC_RL2_thnsparse->GetEntries() << endl;
    cout << endl;

    // Plot pt comparisons
    plot_pt_comparisons(h_truth_jet_pt_thnsparse, h_det_jet_pt_thnsparse, h_matched_truth_jet_pt_thnsparse, h_matched_det_jet_pt_thnsparse, "jet_pt");
    plot_pt_comparisons(h_truth_jet_pt_thnsparse, h_det_jet_pt_thnsparse, h_matched_truth_jet_pt_thnsparse, h_matched_det_jet_pt_thnsparse, "D0_pt");
    

    // Make canvas for text box
    TCanvas * can_text = new TCanvas();
    TPaveText * pt_text = new TPaveText(0.05, 0.05, 0.95, 0.95, "NDC"); // (x1, y1, x2, y2)
    pt_text->SetTextSize(0.04);

    for ( int i = 0; i < n_bins; i ++ ) {

        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];
        int d0_pt_cut = d0_pt_cuts[i];

        cout << " in pt bin " << i << ": " << pt_min <<"-" << pt_max << endl;
        
        // Initialize canvas and legend
        TCanvas * can_eecs = new TCanvas();
        TCanvas * can_ratios = new TCanvas();
        TLegend * leg_eecs = new TLegend(0.3, 0.6, 0.5, 0.8);
        TLegend * leg_ratios = new TLegend(0.6, 0.6, 0.8, 0.8);

        // Format legend
        leg_eecs->SetBorderSize(0);
        leg_ratios->SetBorderSize(0);

        // Make clones.
        THnSparse * h_truth_jet_pt_clone = (THnSparse *) h_truth_jet_pt_thnsparse->Clone((std::string(h_truth_jet_pt_thnsparse->GetName()) + "_clone").c_str());
        THnSparse * h_det_jet_pt_clone = (THnSparse *) h_det_jet_pt_thnsparse->Clone((std::string(h_det_jet_pt_thnsparse->GetName()) + "_clone").c_str());
        THnSparse * h_matched_truth_jet_pt_clone = (THnSparse *) h_matched_truth_jet_pt_thnsparse->Clone((std::string(h_matched_truth_jet_pt_thnsparse->GetName()) + "_clone").c_str());
        THnSparse * h_matched_det_jet_pt_clone = (THnSparse *) h_matched_det_jet_pt_thnsparse->Clone((std::string(h_matched_det_jet_pt_thnsparse->GetName()) + "_clone").c_str());

        THnSparse * h_truth_ENC_RL2_clone = (THnSparse *) h_truth_ENC_RL2_thnsparse->Clone((std::string(h_truth_ENC_RL2_thnsparse->GetName()) + "_clone").c_str());
        THnSparse * h_det_ENC_RL2_clone = (THnSparse *) h_det_ENC_RL2_thnsparse->Clone((std::string(h_det_ENC_RL2_thnsparse->GetName()) + "_clone").c_str());
        THnSparse * h_matched_truth_ENC_RL2_clone = (THnSparse *) h_matched_truth_ENC_RL2_thnsparse->Clone((std::string(h_matched_truth_ENC_RL2_thnsparse->GetName()) + "_clone").c_str());
        THnSparse * h_matched_det_ENC_RL2_clone = (THnSparse *) h_matched_det_ENC_RL2_thnsparse->Clone((std::string(h_matched_det_ENC_RL2_thnsparse->GetName()) + "_clone").c_str());

        // Make cuts
        applyCuts(h_truth_jet_pt_clone, pt_min, pt_max, d0_pt_cut, true);
        applyCuts(h_det_jet_pt_clone, pt_min, pt_max, d0_pt_cut, true);
        applyCuts(h_matched_truth_jet_pt_clone, pt_min, pt_max, d0_pt_cut, true);
        applyCuts(h_matched_det_jet_pt_clone, pt_min, pt_max, d0_pt_cut, true);

        applyCuts(h_truth_ENC_RL2_clone, pt_min, pt_max, d0_pt_cut, true);
        applyCuts(h_det_ENC_RL2_clone, pt_min, pt_max, d0_pt_cut, true);
        applyCuts(h_matched_truth_ENC_RL2_clone, pt_min, pt_max, d0_pt_cut, true);
        applyCuts(h_matched_det_ENC_RL2_clone, pt_min, pt_max, d0_pt_cut, true);

        // Plot pt comparisons -- DELETE LATER
        plot_pt_comparisons(h_truth_ENC_RL2_clone, h_det_ENC_RL2_clone, h_matched_truth_ENC_RL2_clone, h_matched_det_ENC_RL2_clone, Form("jet_pt_pairs_pt%d-%d",pt_min,pt_max));
        plot_pt_comparisons(h_truth_ENC_RL2_clone, h_det_ENC_RL2_clone, h_matched_truth_ENC_RL2_clone, h_matched_det_ENC_RL2_clone, Form("D0_pt_pairs_pt%d-%d",pt_min,pt_max));

        // Project
        TH1D * h_truth_jet_pt = h_truth_jet_pt_clone->Projection(0);
        TH1D * h_det_jet_pt = h_det_jet_pt_clone->Projection(0);
        TH1D * h_matched_truth_jet_pt = h_matched_truth_jet_pt_clone->Projection(0);
        TH1D * h_matched_det_jet_pt = h_matched_det_jet_pt_clone->Projection(0);

        TH1D * h_truth_EEC = h_truth_ENC_RL2_clone->Projection(4);
        TH1D * h_det_EEC = h_det_ENC_RL2_clone->Projection(4);

        TH1D * h_matched_truth_EEC = h_matched_truth_ENC_RL2_clone->Projection(4);
        TH1D * h_matched_det_EEC = h_matched_det_ENC_RL2_clone->Projection(4);
        
        // Scale
        double numjets_truth = h_truth_jet_pt->Integral();
        h_truth_EEC->Scale(1/numjets_truth, "width");

        double numjets_det = h_det_jet_pt->Integral();
        h_det_EEC->Scale(1/numjets_det, "width");

        double numjets_truth_matched = h_matched_truth_jet_pt->Integral();
        h_matched_truth_EEC->Scale(1/numjets_truth_matched, "width");

        double numjets_det_matched = h_matched_det_jet_pt->Integral();
        h_matched_det_EEC->Scale(1/numjets_det_matched, "width");\

        // Record number of jets
        pt_text->AddText(Form("in jet pt=%d-%d:", pt_min, pt_max));
        pt_text->AddText(Form("  # of jets in truth     all: %.3f", h_truth_jet_pt->Integral()));
        pt_text->AddText(Form("  # of jets in det       all: %.3f", h_det_jet_pt->Integral()));
        pt_text->AddText(Form("  # of jets in truth matched: %.3f", h_matched_truth_jet_pt->Integral()));
        pt_text->AddText(Form("  # of jets in det   matched: %.3f", h_matched_det_jet_pt->Integral()));
        pt_text->AddText("  ");
        pt_text->AddText(Form("  # of pairs in truth     all: %.3f", h_truth_EEC->Integral()));
        pt_text->AddText(Form("  # of pairs in det       all: %.3f", h_det_EEC->Integral()));
        pt_text->AddText(Form("  # of pairs in truth matched: %.3f", h_matched_truth_EEC->Integral()));
        pt_text->AddText(Form("  # of pairs in det   matched: %.3f", h_matched_det_EEC->Integral()));

        // Take ratios
        TH1D * h_ratio_unmatched = (TH1D *) h_det_EEC->Clone("h_ratio_unmatched");
        h_ratio_unmatched->Divide(h_truth_EEC);

        TH1D * h_ratio_matched = (TH1D *) h_matched_det_EEC->Clone("h_ratio_matched");
        h_ratio_matched->Divide(h_matched_truth_EEC);

        // Format
        FormatHist(leg_eecs, h_truth_EEC, "All truth", kBlue, kFullStar, 1.0, "R_{L}", "#SigmaEEC");
        FormatHist(leg_eecs, h_det_EEC, "All det", kBlue, kFullCrossX, 1.0, "R_{L}", "#SigmaEEC");
        FormatHist(leg_eecs, h_matched_truth_EEC, "Matched truth", kRed, kFullStar, 1.0, "R_{L}", "#SigmaEEC");
        FormatHist(leg_eecs, h_matched_det_EEC, "Matched det", kRed, kFullCrossX, 1.0, "R_{L}", "#SigmaEEC");
        
        FormatHist(leg_ratios, h_ratio_unmatched, "All bbb", kBlue, kFullCircle, 1.0, "R_{L}", "det/truth");
        FormatHist(leg_ratios, h_ratio_matched, "Matched bbb", kRed, kFullCircle, 1.0, "R_{L}", "det/truth");

        // Reset maximum for plotting purposes
        double max_height = std::max({h_truth_EEC->GetMaximum(), h_det_EEC->GetMaximum(), h_matched_truth_EEC->GetMaximum(), h_matched_det_EEC->GetMaximum()});
        h_truth_EEC->SetMaximum(max_height * 1.1);

        h_ratio_unmatched->SetMaximum(2.5);

        // Plot EECs
        can_eecs->cd();
        gPad->SetLogx();

        h_truth_EEC->Draw();
        h_det_EEC->Draw("SAME");
        h_matched_truth_EEC->Draw("SAME");
        h_matched_det_EEC->Draw("SAME");
        leg_eecs->Draw();

        // Plot ratios
        can_ratios->cd();
        gPad->SetLogx();
        h_ratio_unmatched->Draw();
        h_ratio_matched->Draw("SAME");

        leg_ratios->Draw();
        drawHoriLine(1e-4, 1, 1, kBlack, 3)->Draw();

        // Save
        can_eecs->SaveAs(Form("/software/users/blianggi/mypyjetty/storage/HF_EEC/plots/eecs_pt%d_%d%s.pdf", pt_min, pt_max, output_add_name.c_str()));
        can_ratios->SaveAs(Form("/software/users/blianggi/mypyjetty/storage/HF_EEC/plots/binbybin_ratios_pt%d_%d%s.pdf", pt_min, pt_max, output_add_name.c_str()));

        // Delete hists
        delete h_truth_EEC;
        delete h_det_EEC;
        delete h_matched_truth_EEC;
        delete h_matched_det_EEC;
        delete h_ratio_unmatched;
        delete h_ratio_matched;


    }

    // Plot text
    can_text->cd();
    pt_text->Draw();

    // Save 
    can_text->SaveAs(Form("/software/users/blianggi/mypyjetty/storage/HF_EEC/plots/num_jets_pairs%s.pdf", output_add_name.c_str()));

}



void analyze_HF_fastsim_for_bbb() {

    // Set canvas style
    SetStyle();

    // Read in the input file
    TFile * herwig_fastim_file;
    if (include_dstar == false) {
        herwig_fastim_file= new TFile("/rstorage/generators/herwig_alice/scaling/500341/299990/AnalysisResultsFinal.root", "READ");
    } else {
        herwig_fastim_file= new TFile("/rstorage/generators/herwig_alice/scaling/500841/299990/AnalysisResultsFinal.root", "READ");
        output_add_name = "_withDstar";
    }

    // get and plot herwig bbb corrections
    get_and_plot_bbb(herwig_fastim_file);

}