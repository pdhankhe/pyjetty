// This code will take in HF Herwig fastsim (both truth and detector level)
// and find the bin-by-bin corrections


const int pt_bins[] = { 10, 15, 30 }; //, 100, 150 }; //{ 10, 20, 40 };
const int d0_pt_cuts[] = { 5, 5 }; //, 5, 5 };
const int n_bins = 2;

std::string output_add_name = "";
std::string output_dir = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/HF_EEC/plots"; //hiccup
// std::string output_dir = "/software/users/blianggi/mypyjetty/storage/HF_EEC/plots"; //hiccup

// Set bool for what to run
bool include_dstar = true; // right now true is running on perlmutter. if switching to hiccup, need to find some additional files

class MCHistCollection {
public:
    std::string name;
    
    TH1D * h_jet_pt_gen_all;
    TH1D * h_jet_pt_det_all;
    TH1D * h_jet_pt_gen_matched;
    TH1D * h_jet_pt_det_matched;

    TH1D * h_D0_pt_gen_all;
    TH1D * h_D0_pt_det_all;
    TH1D * h_D0_pt_gen_matched;
    TH1D * h_D0_pt_det_matched;

    TH1D * h_D0_z_gen_all;
    TH1D * h_D0_z_det_all;
    TH1D * h_D0_z_gen_matched;
    TH1D * h_D0_z_det_matched;

    std::vector<TH1D *> h_EEC_gen_all;
    std::vector<TH1D *> h_EEC_det_all;
    std::vector<TH1D *> h_EEC_gen_matched;
    std::vector<TH1D *> h_EEC_det_matched;

    std::vector<TH1D *> h_bbb_ratio_all;
    std::vector<TH1D *> h_bbb_ratio_matched;


    MCHistCollection(std::string name_val) {
        name = name_val;
        
        // filepath_plots = outdir + "/%s/%s/" + name + "/%s"; // ptname, norm_string, filename
        // if (name.find("jet_") != std::string::npos || name.find("const") != std::string::npos) filepath_plots = outdir + "/%s"; // filename

    }

    void addHist(std::vector<TH1D *>& hist_vec, TH1D* hist) {
        hist_vec.push_back(hist);
    }
};


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

void get_jet_quantities(MCHistCollection& MCHists, THnSparse * h_truth_jet_pt_thnsparse, THnSparse * h_det_jet_pt_thnsparse, 
                         THnSparse * h_matched_truth_jet_pt_thnsparse, THnSparse * h_matched_det_jet_pt_thnsparse, 
                         std::string type_pt_name, TLegend * leg_pt, int proj_axis, std::string xtitle, int markerstyle=kFullCircle) {
    
    TH1D * h_truth_all = h_truth_jet_pt_thnsparse->Projection(proj_axis);
    TH1D * h_det_all = h_det_jet_pt_thnsparse->Projection(proj_axis);
    TH1D * h_truth_matched = h_matched_truth_jet_pt_thnsparse->Projection(proj_axis);
    TH1D * h_det_matched = h_matched_det_jet_pt_thnsparse->Projection(proj_axis);

    // Normlize by number of jets (self-normalize)
    double num_truth_all_jets = h_truth_all->Integral();
    double num_det_all_jets = h_det_all->Integral();
    double num_truth_matched_jets = h_truth_matched->Integral();
    double num_det_matched_jets = h_det_matched->Integral();

    h_truth_all->Scale(1/num_truth_all_jets, "width");
    h_det_all->Scale(1/num_det_all_jets, "width");
    h_truth_matched->Scale(1/num_truth_matched_jets, "width");
    h_det_matched->Scale(1/num_det_matched_jets, "width");

    h_truth_all->GetXaxis()->SetRangeUser(0,40);
    h_det_all->GetXaxis()->SetRangeUser(0,40);
    h_truth_matched->GetXaxis()->SetRangeUser(0,40);
    h_det_matched->GetXaxis()->SetRangeUser(0,40);

    // Format histograms
    FormatHist(leg_pt, h_truth_all, "All truth", kRed, markerstyle, 0.8, xtitle.c_str(), "Counts");
    FormatHist(leg_pt, h_det_all, "All det", kBlue, markerstyle, 0.8, xtitle.c_str(), "Counts");
    FormatHist(leg_pt, h_truth_matched, "Matched truth", kRed-7, markerstyle, 0.8, xtitle.c_str(), "Counts");
    FormatHist(leg_pt, h_det_matched, "Matched det", kBlue-9, markerstyle, 0.8, xtitle.c_str(), "Counts");
    // h_det_matched->SetMaximum(h_truth_all->GetMaximum() * 1.2);

    if (type_pt_name.find("jet_pt") != std::string::npos){
        MCHists.h_jet_pt_gen_all = (TH1D*) h_truth_all->Clone(h_truth_all->GetName());
        MCHists.h_jet_pt_det_all = (TH1D*) h_det_all->Clone(h_det_all->GetName());
        MCHists.h_jet_pt_gen_matched = (TH1D*) h_truth_matched->Clone(h_truth_matched->GetName());
        MCHists.h_jet_pt_det_matched = (TH1D*) h_det_matched->Clone(h_det_matched->GetName());
    } else if (type_pt_name.find("D0_pt") != std::string::npos) {
        MCHists.h_D0_pt_gen_all = (TH1D*) h_truth_all->Clone(h_truth_all->GetName());
        MCHists.h_D0_pt_det_all = (TH1D*) h_det_all->Clone(h_det_all->GetName());
        MCHists.h_D0_pt_gen_matched = (TH1D*) h_truth_matched->Clone(h_truth_matched->GetName());
        MCHists.h_D0_pt_det_matched = (TH1D*) h_det_matched->Clone(h_det_matched->GetName());
    } else if (type_pt_name.find("D0_z") != std::string::npos) {
        MCHists.h_D0_z_gen_all = (TH1D*) h_truth_all->Clone(h_truth_all->GetName());
        MCHists.h_D0_z_det_all = (TH1D*) h_det_all->Clone(h_det_all->GetName());
        MCHists.h_D0_z_gen_matched = (TH1D*) h_truth_matched->Clone(h_truth_matched->GetName());
        MCHists.h_D0_z_det_matched = (TH1D*) h_det_matched->Clone(h_det_matched->GetName());
    }

}

void fill_MC_Hists_jet_quantities(MCHistCollection& MCHists, THnSparse * h_truth_jet_pt_thnsparse, THnSparse * h_det_jet_pt_thnsparse, 
                         THnSparse * h_matched_truth_jet_pt_thnsparse, THnSparse * h_matched_det_jet_pt_thnsparse, int markerstyle) {
    // make dummy canvas and tlegend ugh
    TCanvas * c_dum = new TCanvas();
    TLegend * l_dum = new TLegend();

    // jet pt
    std::string type_pt_name = "jet_pt";
    int proj_axis = 0; // if (type_pt_name == "jet_pt") 
    std::string xtitle = "p_{T, jet}";
    get_jet_quantities(MCHists, h_truth_jet_pt_thnsparse, h_det_jet_pt_thnsparse, h_matched_truth_jet_pt_thnsparse, h_matched_det_jet_pt_thnsparse, type_pt_name, l_dum, proj_axis, xtitle, markerstyle);

    // D0 pt
    type_pt_name = "D0_pt";
    proj_axis = 1;
    xtitle = "p_{T, D^{0}}";
    get_jet_quantities(MCHists, h_truth_jet_pt_thnsparse, h_det_jet_pt_thnsparse, h_matched_truth_jet_pt_thnsparse, h_matched_det_jet_pt_thnsparse, type_pt_name, l_dum, proj_axis, xtitle, markerstyle);

    // D0 z
    type_pt_name = "D0_z";
    proj_axis = 3;
    xtitle = "z_{D^{0}}";
    get_jet_quantities(MCHists, h_truth_jet_pt_thnsparse, h_det_jet_pt_thnsparse, h_matched_truth_jet_pt_thnsparse, h_matched_det_jet_pt_thnsparse, type_pt_name, l_dum, proj_axis, xtitle, markerstyle);
    
}


// top_panel_hists_vec is all the histograms that will go in the top panel
// the bottom panel will include ratios of the histograms in the top plot
// the ratios are calculated in this function, so to specify which histograms are divided by which - 
// ratio_num_index_vec will hold the indices of the histograms that will be in the numerator of the ratio, 
// and correspondingly ratio_den_index_vec will hold the indices of the histograms to be put in the denominator
void plot_two_panels(std::vector<TH1D *> top_panel_hists_vec, std::vector<int> ratio_num_index_vec, std::vector<int> ratio_den_index_vec, 
                    std::string xtitle, std::string ratio_ylabel, 
                    std::vector<std::string> top_panel_leg_labels_vec, std::vector<std::string> bottom_panel_leg_labels_vec,
                    std::vector<int> ratio_markercolor_vec, std::vector<int> ratio_markerstyle_vec, std::string output_filepath,
                    int ratio_max=2.) { //}, std::string type_pt_name, std::string gen_name) {

    
    TCanvas * can_twopanel = new TCanvas();
    // Define pad heights (top bigger, bottom smaller)
    float padSplit = 0.3;  // fraction for bottom pad

    // Create top pad
    TPad *pad1 = new TPad("pad1", "pad1", 0, padSplit, 1, 1);
    pad1->SetBottomMargin(0); // remove bottom margin for top pad
    pad1->SetTicks(1,1);
    pad1->Draw();

    TLegend * leg_top = new TLegend(0.6, 0.6, 0.8, 0.8);
    TLegend * leg_bottom = new TLegend(0.6, 0.7, 0.8, 0.9);
    
    // // Get and format ratios
    // TH1D * h_pt_all_ratio = (TH1D *) MCHists.h_jet_pt_det_all->Clone("h_pt_all_ratio");
    // h_pt_all_ratio->Divide(MCHists.h_jet_pt_gen_all);
    // TH1D * h_pt_matched_all_ratio = (TH1D *) MCHists.h_jet_pt_det_matched->Clone("h_pt_matched_all_ratio");
    // h_pt_matched_all_ratio->Divide(MCHists.h_jet_pt_gen_matched);
    // FormatHist(leg_pt_ratio, h_pt_all_ratio, "All", kBlue, kFullCross, 0.8, xtitle.c_str(), "Det / Truth");
    // FormatHist(leg_pt_ratio, h_pt_matched_all_ratio, "Matched truth", kRed, kFullCross, 0.8, xtitle.c_str(), "Det / Truth");
    // h_pt_all_ratio->SetMaximum(2);

    cout << "Size of array! " << ratio_num_index_vec.size() << endl;
    
    // Plot comparison of pt distributions
    pad1->cd();
    gPad->SetLogy();
    for ( int i = 0; i < top_panel_hists_vec.size(); i++ ) {
        leg_top->AddEntry(top_panel_hists_vec[i], top_panel_leg_labels_vec[i].c_str(), "pl");
        top_panel_hists_vec[i]->Draw("SAME");
    }
    leg_top->Draw();

    // Plot ratio of pt distributions
    can_twopanel->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, padSplit);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2); // leave room for x-axis labels
    pad2->SetTicks(1,1);
    pad2->Draw();

    pad2->cd();
    for ( int j = 0; j < ratio_num_index_vec.size(); j++ ) {
        int index_num = ratio_num_index_vec[j];
        int index_den = ratio_den_index_vec[j];
        TH1D * h_ratio = (TH1D *) top_panel_hists_vec[index_num]->Clone(Form("h_ratio_%d", j));
        h_ratio->Divide(top_panel_hists_vec[index_den]);
        
        FormatHist(leg_bottom, h_ratio, bottom_panel_leg_labels_vec[j], ratio_markercolor_vec[j], ratio_markerstyle_vec[j], 0.8, xtitle.c_str(), ratio_ylabel);
        // leg_bottom->AddEntry(h_ratio, bottom_panel_leg_labels_vec[j].c_str(), "pl");
        h_ratio->SetMinimum(0.);
        h_ratio->SetMaximum(ratio_max);
        h_ratio->Draw("SAME");

    }
    leg_bottom->Draw();

    // Save
    // can_twopanel->SaveAs(Form("%s/%s/%s_distribution%s.pdf", output_dir.c_str(), gen_name.c_str(), type_pt_name.c_str(), output_add_name.c_str()));
    can_twopanel->SaveAs(output_filepath.c_str());

}


void plot_pt_comparisons(MCHistCollection& MCHists, THnSparse * h_truth_jet_pt_thnsparse, THnSparse * h_det_jet_pt_thnsparse, 
                         THnSparse * h_matched_truth_jet_pt_thnsparse, THnSparse * h_matched_det_jet_pt_thnsparse, 
                         std::string type_pt_name, std::string gen_name, std::string pr_or_nonpr) {

    
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
    get_jet_quantities(MCHists, h_truth_jet_pt_thnsparse, h_det_jet_pt_thnsparse, h_matched_truth_jet_pt_thnsparse, h_matched_det_jet_pt_thnsparse, type_pt_name, leg_pt, proj_axis, xtitle);

    // Get and format ratios
    TH1D * h_pt_all_ratio = (TH1D *) MCHists.h_jet_pt_det_all->Clone("h_pt_all_ratio");
    h_pt_all_ratio->Divide(MCHists.h_jet_pt_gen_all);
    TH1D * h_pt_matched_all_ratio = (TH1D *) MCHists.h_jet_pt_det_matched->Clone("h_pt_matched_all_ratio");
    h_pt_matched_all_ratio->Divide(MCHists.h_jet_pt_gen_matched);
    FormatHist(leg_pt_ratio, h_pt_all_ratio, "All", kBlue, kFullCross, 0.8, xtitle.c_str(), "Det / Truth");
    FormatHist(leg_pt_ratio, h_pt_matched_all_ratio, "Matched truth", kRed, kFullCross, 0.8, xtitle.c_str(), "Det / Truth");
    h_pt_all_ratio->SetMaximum(2);
    
    // Plot comparison of pt distributions
    pad1->cd();
    gPad->SetLogy();
    MCHists.h_jet_pt_det_matched->Draw(); // drawing this first to get the minimum right
    MCHists.h_jet_pt_det_all->Draw("SAME");
    MCHists.h_jet_pt_gen_matched->Draw("SAME");
    MCHists.h_jet_pt_gen_all->Draw("SAME");
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
    can_pt->SaveAs(Form("%s/%s/%s_distribution%s_%s.pdf", output_dir.c_str(), gen_name.c_str(), type_pt_name.c_str(), output_add_name.c_str(), pr_or_nonpr.c_str()));

    // delete h_pt_matched_det_all;
    // delete h_pt_det_all;
    // delete h_pt_matched_truth_all;
    // delete h_pt_truth_all;
    // delete h_pt_all_ratio;
    // delete h_pt_matched_all_ratio;
}



// returns a vector of histograms of bin by bin factors, 
// where every two indices is one pt bin
// and they are in the order all jets, matched jets
std::vector<TH1D *> get_and_plot_bbb(TFile *file, std::string gen_name, MCHistCollection& MCHists, std::string pr_or_nonpr) {

    // initialize return vector
    std::vector<TH1D *> vec_bbb_factors;
    
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

    // fill jet quantities of mc collection hists
    int markerstyle = kFullStar; // if herwig
    if ( gen_name == "pythia" ) markerstyle = kFullCross;
    fill_MC_Hists_jet_quantities(MCHists, h_truth_jet_pt_thnsparse, h_det_jet_pt_thnsparse, h_matched_truth_jet_pt_thnsparse, h_matched_det_jet_pt_thnsparse, markerstyle);

    // // Plot pt comparisons
    // plot_pt_comparisons(h_truth_jet_pt_thnsparse, h_det_jet_pt_thnsparse, h_matched_truth_jet_pt_thnsparse, h_matched_det_jet_pt_thnsparse, "jet_pt", gen_name, pr_or_nonpr);
    // plot_pt_comparisons(h_truth_jet_pt_thnsparse, h_det_jet_pt_thnsparse, h_matched_truth_jet_pt_thnsparse, h_matched_det_jet_pt_thnsparse, "D0_pt", gen_name, pr_or_nonpr);
    

    // // Make canvas for text box
    // TCanvas * can_text = new TCanvas();
    // TPaveText * pt_text = new TPaveText(0.05, 0.05, 0.95, 0.95, "NDC"); // (x1, y1, x2, y2)
    // pt_text->SetTextSize(0.04);

    for ( int i = 0; i < n_bins; i ++ ) {

        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];
        int d0_pt_cut = d0_pt_cuts[i];

        cout << "individual gen // in pt bin " << i << ": " << pt_min <<"-" << pt_max << endl;
        
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

        // // Plot pt comparisons -- DELETE LATER
        // plot_pt_comparisons(h_truth_ENC_RL2_clone, h_det_ENC_RL2_clone, h_matched_truth_ENC_RL2_clone, h_matched_det_ENC_RL2_clone, Form("jet_pt_pairs_pt%d-%d",pt_min,pt_max), gen_name, pr_or_nonpr);
        // plot_pt_comparisons(h_truth_ENC_RL2_clone, h_det_ENC_RL2_clone, h_matched_truth_ENC_RL2_clone, h_matched_det_ENC_RL2_clone, Form("D0_pt_pairs_pt%d-%d",pt_min,pt_max), gen_name, pr_or_nonpr);

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

        // // Record number of jets
        // pt_text->AddText(Form("in jet pt=%d-%d:", pt_min, pt_max));
        // pt_text->AddText(Form("  # of jets in truth     all: %.3f", h_truth_jet_pt->Integral()));
        // pt_text->AddText(Form("  # of jets in det       all: %.3f", h_det_jet_pt->Integral()));
        // pt_text->AddText(Form("  # of jets in truth matched: %.3f", h_matched_truth_jet_pt->Integral()));
        // pt_text->AddText(Form("  # of jets in det   matched: %.3f", h_matched_det_jet_pt->Integral()));
        // pt_text->AddText("  ");
        // pt_text->AddText(Form("  # of pairs in truth     all: %.3f", h_truth_EEC->Integral()));
        // pt_text->AddText(Form("  # of pairs in det       all: %.3f", h_det_EEC->Integral()));
        // pt_text->AddText(Form("  # of pairs in truth matched: %.3f", h_matched_truth_EEC->Integral()));
        // pt_text->AddText(Form("  # of pairs in det   matched: %.3f", h_matched_det_EEC->Integral()));

        // Take ratios
        TH1D * h_ratio_unmatched = (TH1D *) h_det_EEC->Clone("h_ratio_unmatched_clone");
        h_ratio_unmatched->SetDirectory(nullptr); 
        h_ratio_unmatched->Divide(h_truth_EEC);

        TH1D * h_ratio_matched = (TH1D *) h_matched_det_EEC->Clone("h_ratio_matched_clone");
        h_ratio_matched->SetDirectory(nullptr); 
        h_ratio_matched->Divide(h_matched_truth_EEC);

        // Format
        FormatHist(leg_eecs, h_truth_EEC, "All truth", kBlue, kFullStar, 1.0, "R_{L}", "#SigmaEEC");
        FormatHist(leg_eecs, h_det_EEC, "All det", kBlue, kFullCrossX, 1.0, "R_{L}", "#SigmaEEC");
        FormatHist(leg_eecs, h_matched_truth_EEC, "Matched truth", kRed, kFullStar, 1.0, "R_{L}", "#SigmaEEC");
        FormatHist(leg_eecs, h_matched_det_EEC, "Matched det", kRed, kFullCrossX, 1.0, "R_{L}", "#SigmaEEC");
        
        FormatHist(leg_ratios, h_ratio_unmatched, "All bbb", kBlue, kFullCircle, 1.0, "R_{L}", "det/truth");
        FormatHist(leg_ratios, h_ratio_matched, "Matched bbb", kRed, kFullCircle, 1.0, "R_{L}", "det/truth");

        // prevent root from auto-deleting histograms
        // h_truth_EEC->SetDirectory(nullptr); 
        // h_det_EEC->SetDirectory(nullptr); 
        // h_matched_truth_EEC->SetDirectory(nullptr); 
        // h_matched_det_EEC->SetDirectory(nullptr); 
        // h_ratio_unmatched->SetDirectory(nullptr); 
        // h_ratio_matched->SetDirectory(nullptr);

        // Add to collection
        MCHists.addHist(MCHists.h_EEC_gen_all, (TH1D*) h_truth_EEC->Clone(h_truth_EEC->GetName()));
        MCHists.addHist(MCHists.h_EEC_det_all, (TH1D*) h_det_EEC->Clone(h_det_EEC->GetName()));
        MCHists.addHist(MCHists.h_EEC_gen_matched, (TH1D*) h_matched_truth_EEC->Clone(h_matched_truth_EEC->GetName()));
        MCHists.addHist(MCHists.h_EEC_det_matched, (TH1D*) h_matched_det_EEC->Clone(h_matched_det_EEC->GetName()));
        MCHists.addHist(MCHists.h_bbb_ratio_all, (TH1D*) h_ratio_unmatched->Clone(h_ratio_unmatched->GetName()));
        MCHists.addHist(MCHists.h_bbb_ratio_matched, (TH1D*) h_ratio_matched->Clone(h_ratio_matched->GetName()));

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
        vec_bbb_factors.push_back((TH1D*)h_ratio_unmatched->Clone(Form("h_ratio_unmatched_pt%d-%d", pt_min, pt_max)));
        vec_bbb_factors.push_back((TH1D*)h_ratio_matched->Clone(Form("h_ratio_matched%d-%d", pt_min, pt_max)));

        leg_ratios->Draw();
        drawHoriLine(1e-4, 1, 1, kBlack, 3)->Draw();

        // Save
        can_eecs->SaveAs(Form("%s/%s/eecs_pt%d_%d%s_%s.pdf", output_dir.c_str(), gen_name.c_str(), pt_min, pt_max, output_add_name.c_str(), pr_or_nonpr.c_str()));
        can_ratios->SaveAs(Form("%s/%s/binbybin_ratios_pt%d_%d%s_%s.pdf", output_dir.c_str(), gen_name.c_str(), pt_min, pt_max, output_add_name.c_str(), pr_or_nonpr.c_str()));

        // Delete hists
        delete h_truth_jet_pt_clone;
        delete h_det_jet_pt_clone;
        delete h_matched_truth_jet_pt_clone;
        delete h_matched_det_jet_pt_clone;
        delete h_truth_ENC_RL2_clone;
        delete h_det_ENC_RL2_clone;
        delete h_matched_truth_ENC_RL2_clone;
        delete h_matched_det_ENC_RL2_clone;

        delete h_truth_jet_pt;
        delete h_det_jet_pt;
        delete h_matched_truth_jet_pt;
        delete h_matched_det_jet_pt;

        delete h_truth_EEC;
        delete h_det_EEC;
        delete h_matched_truth_EEC;
        delete h_matched_det_EEC;
        delete h_ratio_unmatched;
        delete h_ratio_matched;


    }

    // // Plot text
    // can_text->cd();
    // pt_text->Draw();

    // // Save 
    // can_text->SaveAs(Form("%s/%s/num_jets_pairs%s.pdf", output_dir.c_str(), gen_name.c_str(), output_add_name.c_str()));

    return vec_bbb_factors;
}

void plot_different_generators_together(std::vector<TH1D *> vec_bbb_herwig, std::vector<TH1D *> vec_bbb_pythia, MCHistCollection herwigHists, MCHistCollection pythiaHists, std::string pr_or_nonpr) {
    
    for ( int i = 0; i < n_bins; i ++ ) {

        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];
        int d0_pt_cut = d0_pt_cuts[i];
        
        // Initialize canvas and legend
        TCanvas * can_all_ratios = new TCanvas();
        TLegend * leg_all_ratios = new TLegend(0.6, 0.6, 0.8, 0.8);
        can_all_ratios->cd();
        gPad->SetLogx();

        // Format legend
        leg_all_ratios->SetBorderSize(0);

        // Format histograms
        FormatHist(leg_all_ratios, vec_bbb_herwig[i*2], "Herwig all jets", kBlue, kFullCrossX, 1.0, "R_{L}", "det/truth");
        FormatHist(leg_all_ratios, vec_bbb_herwig[i*2+1], "Herwig matched jets", kBlue, kFullCross, 1.0, "R_{L}", "det/truth");
        FormatHist(leg_all_ratios, vec_bbb_pythia[i*2], "Pythia all jets", kRed, kFullCrossX, 1.0, "R_{L}", "det/truth");
        FormatHist(leg_all_ratios, vec_bbb_pythia[i*2+1], "Pythia matched jets", kRed, kFullCross, 1.0, "R_{L}", "det/truth");

        // Draw histograms
        vec_bbb_herwig[i*2]->Draw("SAME");
        vec_bbb_herwig[i*2+1]->Draw("SAME");
        vec_bbb_pythia[i*2]->Draw("SAME");
        vec_bbb_pythia[i*2+1]->Draw("SAME");
        drawHoriLine(1e-4, 1, 1, kBlack, 3)->Draw();

        // Add text box with information
        TLegend* text_imp = new TLegend(0.1, 0.6, 0.5, 0.85);
        text_imp->SetBorderSize(0);
        text_imp->SetFillStyle(0);
        text_imp->SetTextSize(0.03);
        text_imp->AddEntry((TObject*)0, "pp, #sqrt{s} = 13 TeV", "");
        text_imp->AddEntry((TObject*)0, "prompt D^{0}-tagged ch. jets", "");
        text_imp->AddEntry((TObject*)0, "anti-k_{T}, R = 0.4", "");
        text_imp->AddEntry((TObject*)0, Form("%d #leq p_{T}^{ch. jet} < %d GeV/c, |#eta_{jet}| #leq 0.5", pt_min, pt_max), "");
        text_imp->AddEntry((TObject*)0, Form("%d #leq p_{T}^{D^{0}} < %d GeV/c, |y_{D^{0}}| #leq 0.8", d0_pt_cut, pt_max), "");
        text_imp->Draw();
        leg_all_ratios->Draw();

        can_all_ratios->SaveAs(Form("%s/pythiaherwig_binbybin_ratios_pt%d_%d%s_%s.pdf", output_dir.c_str(), pt_min, pt_max, output_add_name.c_str(), pr_or_nonpr.c_str()));

    }

    // do jet quantities here
    std::vector<int> ratio_num_index_vec = { 0, 1, 2, 3 };
    std::vector<int> ratio_den_index_vec = { 4, 5, 6, 7 }; 
    std::string ratio_ylabel = "Herwig / Pythia";
    std::vector<std::string> top_panel_leg_labels_vec = { "Herwig gen all", "Herwig det all", "Herwig gen matched", "Herwig det matched", "Pythia gen all", "Pythia det all", "Pythia gen matched", "Pythia det matched" };
    std::vector<std::string> bottom_panel_leg_labels_vec = { "gen all jets", "det all jets", "gen matched jets", "det matched jets" };
    std::vector<int> ratio_markercolor_vec = { kRed, kBlue, kRed-7, kBlue-9 };
    std::vector<int> ratio_markerstyle_vec = { kFullCrossX, kFullCrossX, kFullCrossX, kFullCrossX };
    
    // jet pt
    std::vector<TH1D *> jet_pt_hists_vec = { herwigHists.h_jet_pt_gen_all, herwigHists.h_jet_pt_det_all, herwigHists.h_jet_pt_gen_matched, herwigHists.h_jet_pt_det_matched, 
                                              pythiaHists.h_jet_pt_gen_all, pythiaHists.h_jet_pt_det_all, pythiaHists.h_jet_pt_gen_matched, pythiaHists.h_jet_pt_det_matched };
    std::string jet_pt_xtitle = "p_{T, jet}";
    std::string jet_pt_output_filepath = Form("%s/pythiaherwig_jet_pt%s_%s.pdf", output_dir.c_str(), output_add_name.c_str(), pr_or_nonpr.c_str());
    plot_two_panels(jet_pt_hists_vec, ratio_num_index_vec, ratio_den_index_vec, jet_pt_xtitle, ratio_ylabel, top_panel_leg_labels_vec, bottom_panel_leg_labels_vec, ratio_markercolor_vec, ratio_markerstyle_vec, jet_pt_output_filepath, 5);

    // D0 pt
    std::vector<TH1D *> D0_pt_hists_vec = { herwigHists.h_D0_pt_gen_all, herwigHists.h_D0_pt_det_all, herwigHists.h_D0_pt_gen_matched, herwigHists.h_D0_pt_det_matched, 
                                              pythiaHists.h_D0_pt_gen_all, pythiaHists.h_D0_pt_det_all, pythiaHists.h_D0_pt_gen_matched, pythiaHists.h_D0_pt_det_matched };
    std::string D0_pt_xtitle = "p_{T, D^{0}}";
    std::string D0_pt_output_filepath = Form("%s/pythiaherwig_D0_pt%s_%s.pdf", output_dir.c_str(), output_add_name.c_str(), pr_or_nonpr.c_str());
    plot_two_panels(D0_pt_hists_vec, ratio_num_index_vec, ratio_den_index_vec, D0_pt_xtitle, ratio_ylabel, top_panel_leg_labels_vec, bottom_panel_leg_labels_vec, ratio_markercolor_vec, ratio_markerstyle_vec, D0_pt_output_filepath, 5);

    // D0 z
    std::vector<TH1D *> D0_z_hists_vec = { herwigHists.h_D0_z_gen_all, herwigHists.h_D0_z_det_all, herwigHists.h_D0_z_gen_matched, herwigHists.h_D0_z_det_matched, 
                                              pythiaHists.h_D0_z_gen_all, pythiaHists.h_D0_z_det_all, pythiaHists.h_D0_z_gen_matched, pythiaHists.h_D0_z_det_matched };
    std::string D0_z_xtitle = "z_{D^{0}}";
    std::string D0_z_output_filepath = Form("%s/pythiaherwig_D0_z%s_%s.pdf", output_dir.c_str(), output_add_name.c_str(), pr_or_nonpr.c_str());
    plot_two_panels(D0_z_hists_vec, ratio_num_index_vec, ratio_den_index_vec, D0_z_xtitle, ratio_ylabel, top_panel_leg_labels_vec, bottom_panel_leg_labels_vec, ratio_markercolor_vec, ratio_markerstyle_vec, D0_z_output_filepath);
    
    
}

void analyze_files(TFile * herwig_fastsim_file, TFile * pythia_fastsim_file, MCHistCollection& herwigHists, MCHistCollection& pythiaHists, std::string pr_or_nonpr) {
    // get and plot herwig bbb corrections
    std::vector<TH1D *> vec_bbb_herwig = get_and_plot_bbb(herwig_fastsim_file, "herwig", herwigHists, pr_or_nonpr);

    // get and plot pythia bbb corrections
    std::vector<TH1D *> vec_bbb_pythia = get_and_plot_bbb(pythia_fastsim_file, "pythia", pythiaHists, pr_or_nonpr);

    // plot together
    plot_different_generators_together(vec_bbb_herwig, vec_bbb_pythia, herwigHists, pythiaHists, pr_or_nonpr);

    cout << "size check! " << herwigHists.h_bbb_ratio_all.size() << " vs " << pythiaHists.h_bbb_ratio_all.size() << endl;
}


void format_hist_for_allcombinedplot(TH1D * hist, TLegend * leg, int markercolor, int markerstyle, std::string label) {
    hist->SetMarkerColor(markercolor);
    hist->SetLineColor(markercolor);
    hist->SetMarkerStyle(markerstyle);
    // hist->SetMarkerSize(1.25);

    leg->AddEntry(hist, label.c_str(), "lp");

}

void plot_prompt_and_nonprompt(MCHistCollection herwigHists, MCHistCollection pythiaHists, MCHistCollection herwigNPHists, MCHistCollection pythiaNPHists) {
    for ( int i = 0; i < n_bins; i ++ ) {
        cout << "in pt bin " << i << ": " << pt_bins[i] <<"-" << pt_bins[i+1] << endl;

        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];

        TCanvas * can = new TCanvas();
        can->cd();
        gPad->SetLogx();
        TLegend * leg_ev = new TLegend(0.15, 0.6, 0.4, 0.85);
        TLegend * leg = new TLegend(0.6, 0.55, 0.8, 0.85);
        leg_ev->SetBorderSize(0);
        leg->SetBorderSize(0);

        cout << "herwigHists.h_bbb_ratio_all.size: " << herwigHists.h_bbb_ratio_all.size() << endl;

        leg_ev->AddEntry((TObject*)0, "pp, #sqrt{s} = 13 TeV", "");
        leg_ev->AddEntry((TObject*)0, "D^{0}-tagged ch. jets", "");
        leg_ev->AddEntry((TObject*)0, "anti-k_{T}, R = 0.4", "");
        leg_ev->AddEntry((TObject*)0, Form("%d #leq p_{T}^{ch. jet} < %d GeV/c, |#eta_{jet}| #leq 0.5", pt_min, pt_max), "");
        leg_ev->AddEntry((TObject*)0, Form("%d #leq p_{T}^{D^{0}} < %d GeV/c, |y_{D^{0}}| #leq 0.8", d0_pt_cuts[i], pt_max), "");

        leg->AddEntry((TObject*)0, "Herwig", "");
        format_hist_for_allcombinedplot(herwigHists.h_bbb_ratio_all[i], leg, kRed-4, 20, "Prompt");
        format_hist_for_allcombinedplot(herwigHists.h_bbb_ratio_matched[i], leg, kRed-4, 21, "Prompt matched");
        format_hist_for_allcombinedplot(herwigNPHists.h_bbb_ratio_all[i], leg, kOrange-3, 20, "Non-prompt");
        format_hist_for_allcombinedplot(herwigNPHists.h_bbb_ratio_matched[i], leg, kOrange-3, 21, "Non-prompt matched");
        leg->AddEntry((TObject*)0, "PYTHIA 8", ""); // empty line
        format_hist_for_allcombinedplot(pythiaHists.h_bbb_ratio_all[i], leg, kAzure-2, 20, "Prompt");
        format_hist_for_allcombinedplot(pythiaHists.h_bbb_ratio_matched[i], leg, kAzure-2, 21, "Prompt matched");
        format_hist_for_allcombinedplot(pythiaNPHists.h_bbb_ratio_all[i], leg, kCyan-3, 20, "Non-prompt");
        format_hist_for_allcombinedplot(pythiaNPHists.h_bbb_ratio_matched[i], leg, kCyan-3, 21, "Non-prompt matched");

        herwigHists.h_bbb_ratio_all[i]->SetMaximum(2.5);
        
        herwigHists.h_bbb_ratio_all[i]->Draw();
        herwigHists.h_bbb_ratio_matched[i]->Draw("SAME");
        herwigNPHists.h_bbb_ratio_all[i]->Draw("SAME");
        herwigNPHists.h_bbb_ratio_matched[i]->Draw("SAME");
        pythiaHists.h_bbb_ratio_all[i]->Draw("SAME");
        pythiaHists.h_bbb_ratio_matched[i]->Draw("SAME");
        pythiaNPHists.h_bbb_ratio_all[i]->Draw("SAME");
        pythiaNPHists.h_bbb_ratio_matched[i]->Draw("SAME");
        leg_ev->Draw();
        leg->Draw();
        drawHoriLine(1e-4, 1, 1, kBlack, 3)->Draw();

        can->SaveAs(Form("%s/ALL_pythiaherwig_promptnonprompt_binbybin_ratios_pt%d_%d%s.pdf", output_dir.c_str(), pt_min, pt_max, output_add_name.c_str()));
    }
}



void analyze_HF_fastsim_for_bbb() {

    // Set canvas style
    SetStyle();

    // Read in the input file
    TFile * herwig_fastsim_file;
    TFile * pythia_fastsim_file;
    if (include_dstar == false) {
        herwig_fastsim_file = new TFile("/rstorage/generators/herwig_alice/scaling/500341/299990/AnalysisResultsFinal.root", "READ"); // only available on hiccup rn
        cout << "no pythia fastsim file specified for this option! Exiting..." << endl;
        return;
    } else {
        herwig_fastsim_file = new TFile("/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/storage/herwig/500841/299990/AnalysisResultsFinal.root", "READ"); // perlmutter link
        pythia_fastsim_file = new TFile("/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/pythiagen/scaling/45190202/45154942/AnalysisResultsFinal.root", "READ"); // perlmutter link
        output_add_name = "_withDstar";
    }

    MCHistCollection pythiaHists("pythia");
    MCHistCollection herwigHists("herwig");

    // get and plot herwig + pythia bbb corrections -- prompt
    analyze_files(herwig_fastsim_file, pythia_fastsim_file, herwigHists, pythiaHists, "prompt");

    // get and plot herwig + pythia bbb corrections -- nonprompt
    TFile * herwig_fastsim_nonprompt_file;
    TFile * pythia_fastsim_nonprompt_file;
    MCHistCollection pythiaNPHists("pythia_nonprompt");
    MCHistCollection herwigNPHists("herwig_nonprompt");
    if (include_dstar == true) {
        herwig_fastsim_nonprompt_file = new TFile("/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/storage/herwig/517789/515788/AnalysisResultsFinal.root", "READ"); // perlmutter link
        pythia_fastsim_nonprompt_file = new TFile("/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/pythiagen/scaling/46372190/46293548/AnalysisResultsFinal.root", "READ"); // perlmutter link
        output_add_name = "_withDstar";
    }
    analyze_files(herwig_fastsim_nonprompt_file, pythia_fastsim_nonprompt_file, herwigNPHists, pythiaNPHists, "nonprompt");

    cout << "done individuals. now plotting prompt vs non-prompt together..." << endl;

    plot_prompt_and_nonprompt(herwigHists, pythiaHists, herwigNPHists, pythiaNPHists);

}
