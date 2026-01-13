// This script is to read in the (already made) data and generator histograms
// and plot them together!
// with ratios!


Double_t colors[16] = {kMagenta, kBlue, kOrange+1, kViolet+1, kGreen+2, kRed, kYellow+1, kCyan+1}; //kGray, 
Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
Double_t marker_size = 1.5;

const double ptRL_bins[5] = { 2e-1, 8e-1, 5.0, 10.0, 30.0 };
const int n_ptRLbins = 4;
std::string attempt_dir = "data_with_generators";

// Define class Observable
class Observable {
public:
    std::string name;
    
    std::string axis_label;
    std::string cs_label; //cross section label
    std::string histname_base;
    std::string filepath_plots;

    TFile * input_data_file;
    TFile * input_pythia_file;
    TFile * input_herwig_file;
    TFile * input_anchmc_file;

    std::vector<TH1D*> obs_vec;

    Observable(std::string name_val, //int num_bins_val, double min_bound_val, double max_bound_val, 
               std::string axis_label_val, std::string cs_label_val) {
        name = name_val;
        
        axis_label = axis_label_val;
        cs_label = cs_label_val; //cross section label, in y axis

        histname_base = "h_" + name + "%s_R%s_t%s_pt%d-%d_pTRL%.1f-%.1f_%s%s"; // will be filled with: weight_str, jetR, threshold, pt_min, pt_max, pTRL_min, pTRL_max, norm_string

        filepath_plots = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/plots/" + attempt_dir;
        // filepath_plots = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/plots/" + attempt_dir + "%s/%s/" + name + "/%s"; // ptname, norm_string, filename
        // if (name.find("jet_") != std::string::npos || name.find("const") != std::string::npos) filepath_plots = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/plots/" + attempt_dir + "%s"; // filename

        input_data_file = new TFile(Form("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/rootfiles/data_fifthattempt_ptrlbins/DataHists_%s.root", name.c_str()), "READ");
        input_pythia_file = new TFile(Form("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/rootfiles/pythia5TeV_histograms_crosscheck/PYTHIAHists_%s.root", name.c_str()), "READ");
        input_herwig_file = new TFile(Form("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/rootfiles/herwig_firstattempt/HERWIGHists_%s.root", name.c_str()), "READ");
        input_anchmc_file = new TFile(Form("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/rootfiles/anchMC_firstattempt/ANCHMCHists_%s.root", name.c_str()), "READ");
    }

    void addHist(TH1D* hist) {
        obs_vec.push_back(hist);
    }
};


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

// Function to format histogram
void Format1DHist(Observable obs, TH1 *hist, TLegend* leg, TString leg_text, bool drawline=false, int markercolor=1, int markerstyle=8, double markeralpha=1.,
                double xtitlesize=0.06, double xlabelsize=0.05, double xoffset=1.0,
                double ytitlesize=0.06, double ylabelsize=0.05, double yoffset=1.05, double markersize=1.0) 
                // double xtitlesize=0.04, double xlabelsize=0.04, double xoffset=1.2,
                // double ytitlesize=0.04, double ylabelsize=0.04, double yoffset=1.0) 
{
    hist->SetLineColorAlpha(markercolor, markeralpha);
    hist->SetMarkerColorAlpha(markercolor, markeralpha);
    hist->SetMarkerStyle(markerstyle);
    hist->SetMarkerSize(markersize);
    leg->AddEntry(hist, leg_text, "pl");

    if (drawline) {

        for (int a=0; a < hist->GetNbinsX(); a++) hist->SetBinError(a+1, 0);
        
        hist->SetMarkerStyle(0); //20);
        // hist->SetMarkerColorAlpha(markercolor, 0);
        hist->SetMarkerSize(0);

        // hist->SetFillStyle(0);
        hist->SetLineColorAlpha(markercolor, markeralpha);
        // hist->SetFillColor(markercolor);
        hist->SetLineStyle(markerstyle); //1);
        hist->SetLineWidth(3);
    }

    hist->GetXaxis()->SetTitle(obs.axis_label.c_str());
    hist->GetYaxis()->SetTitle(obs.cs_label.c_str());

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


TH1D * get_1D_histogram(TFile * f, std::string histname) {
    TH1D * hist = (TH1D *) f->Get(histname.c_str());
    return hist;
}

TH2D * get_2D_histogram(TFile * f, std::string histname) {
    TH2D * hist = (TH2D *) f->Get(histname.c_str()); 
    return hist;
}

// For each observable, plot data and pythia together.
// For now, this is the uncorrected data with the rec-level MC
// Can also do corrected data with the truth-level MC
void plot_1D_obs(Observable obs, int pt_min, int pt_max, //int ptrl_bin, 
                 std::string weight_str, std::string jetR, std::string threshold, std::string norm_string, 
                 int option) {

    cout << "in obs " << obs.name << ", pt: " << pt_min << "-" << pt_max << ", option " << option << endl;

    TCanvas * can_obs_all = new TCanvas();
    can_obs_all->cd();
    gPad->SetLogy();

    TLegend *leg_ev = new TLegend(0.18, 0.2, 0.33, 0.35);
    leg_ev->SetTextSize(0.03);
    leg_ev->AddEntry((TObject*)0, "pp, #sqrt{s} = 5.02 TeV", "");
    leg_ev->AddEntry((TObject*)0, "all ch. jets", "");
    leg_ev->AddEntry((TObject*)0, "anti-k_{T}, R = 0.4", "");
    leg_ev->AddEntry((TObject*)0, Form("%d #leq p_{T}^{ch. jet} < %d GeV/c, |#eta_{jet}| #leq 0.5", pt_min, pt_max), "");

    TLegend *leg_obs_all = new TLegend(0.65, 0.5, 0.88, 0.88);
    TLegend *leg_ratio_all = new TLegend(0.7, 0.7, 0.88, 0.88);
    leg_obs_all->SetTextSize(0.03);
    

    for ( int j = 0; j < n_ptRLbins; j++ ) { 

        double pTRL_min = ptRL_bins[j];
        double pTRL_max = ptRL_bins[j+1];

        // get histogram for data and pythia
        std::string histname_1;
        std::string histname_2;
        std::string histname_3;
        std::string histname_4;
        TFile * infile_1;
        TFile * infile_2;
        TFile * infile_3;
        TFile * infile_4;
        std::string label_1;
        std::string label_2;
        std::string label_3;
        std::string label_4;

        if ( option == 1 ) { // data (raw) vs pythia (det) vs herwig (det)
            histname_1 = Form(obs.histname_base.c_str(), weight_str.c_str(), jetR.c_str(), threshold.c_str(), pt_min, pt_max, pTRL_min, pTRL_max, norm_string.c_str(), "");// TODO: FIX WHAT IS BEING ACCESSED HERE, include jet pt and ptrl bin?
            histname_2 = Form(obs.histname_base.c_str(), weight_str.c_str(), jetR.c_str(), threshold.c_str(), pt_min, pt_max, pTRL_min, pTRL_max, norm_string.c_str(), "_det");// TODO: FIX WHAT IS BEING ACCESSED HERE, include jet pt and ptrl bin?
            histname_3 = Form(obs.histname_base.c_str(), weight_str.c_str(), jetR.c_str(), threshold.c_str(), pt_min, pt_max, pTRL_min, pTRL_max, norm_string.c_str(), "_det");// TODO: FIX WHAT IS BEING ACCESSED HERE, include jet pt and ptrl bin?
            histname_4 = Form(obs.histname_base.c_str(), weight_str.c_str(), jetR.c_str(), threshold.c_str(), pt_min, pt_max, pTRL_min, pTRL_max, norm_string.c_str(), "_det");// TODO: FIX WHAT IS BEING ACCESSED HERE, include jet pt and ptrl bin?
            infile_1 = obs.input_data_file;
            infile_2 = obs.input_pythia_file;
            infile_3 = obs.input_herwig_file;
            infile_4 = obs.input_anchmc_file;
            label_1 = "Raw data";
            label_2 = "det-level PYTHIA 8";
            label_3 = "det-level HERWIG 7";
            label_4 = "det-level anchored MC LHC23a3";
        } else if ( option == 2 ) { // pythia (gen) vs pythia (det)
            histname_1 = Form(obs.histname_base.c_str(), weight_str.c_str(), jetR.c_str(), threshold.c_str(), pt_min, pt_max, pTRL_min, pTRL_max, norm_string.c_str(), "_truth");// TODO: FIX WHAT IS BEING ACCESSED HERE, include jet pt and ptrl bin?
            histname_2 = Form(obs.histname_base.c_str(), weight_str.c_str(), jetR.c_str(), threshold.c_str(), pt_min, pt_max, pTRL_min, pTRL_max, norm_string.c_str(), "_det");// TODO: FIX WHAT IS BEING ACCESSED HERE, include jet pt and ptrl bin?
            infile_1 = obs.input_pythia_file;
            infile_2 = obs.input_pythia_file;
            label_1 = "truth-level PYTHIA 8";
            label_2 = "det-level PYTHIA 8";
        } else if ( option == 3 ) { // herwig (gen) vs herwig (det)
            histname_1 = Form(obs.histname_base.c_str(), weight_str.c_str(), jetR.c_str(), threshold.c_str(), pt_min, pt_max, pTRL_min, pTRL_max, norm_string.c_str(), "_truth");// TODO: FIX WHAT IS BEING ACCESSED HERE, include jet pt and ptrl bin?
            histname_2 = Form(obs.histname_base.c_str(), weight_str.c_str(), jetR.c_str(), threshold.c_str(), pt_min, pt_max, pTRL_min, pTRL_max, norm_string.c_str(), "_det");// TODO: FIX WHAT IS BEING ACCESSED HERE, include jet pt and ptrl bin?
            infile_1 = obs.input_herwig_file;
            infile_2 = obs.input_herwig_file;
            label_1 = "truth-level HERWIG 7";
            label_2 = "det-level HERWIG 7";
        } else if ( option == 4 ) { // pythia (det) vs anch mc (det)
            histname_1 = Form(obs.histname_base.c_str(), weight_str.c_str(), jetR.c_str(), threshold.c_str(), pt_min, pt_max, pTRL_min, pTRL_max, norm_string.c_str(), "_det");// TODO: FIX WHAT IS BEING ACCESSED HERE, include jet pt and ptrl bin?
            histname_2 = Form(obs.histname_base.c_str(), weight_str.c_str(), jetR.c_str(), threshold.c_str(), pt_min, pt_max, pTRL_min, pTRL_max, norm_string.c_str(), "_det");// TODO: FIX WHAT IS BEING ACCESSED HERE, include jet pt and ptrl bin?
            infile_1 = obs.input_pythia_file;
            infile_2 = obs.input_anchmc_file;
            label_1 = "det-level PYTHIA 8";
            label_2 = "det-level anchored MC LHC23a3";
        } else if ( option == 5 ) { // pythia (gen) vs anch mc (gen)
            histname_1 = Form(obs.histname_base.c_str(), weight_str.c_str(), jetR.c_str(), threshold.c_str(), pt_min, pt_max, pTRL_min, pTRL_max, norm_string.c_str(), "_truth");// TODO: FIX WHAT IS BEING ACCESSED HERE, include jet pt and ptrl bin?
            histname_2 = Form(obs.histname_base.c_str(), weight_str.c_str(), jetR.c_str(), threshold.c_str(), pt_min, pt_max, pTRL_min, pTRL_max, norm_string.c_str(), "_truth");// TODO: FIX WHAT IS BEING ACCESSED HERE, include jet pt and ptrl bin?
            infile_1 = obs.input_pythia_file;
            infile_2 = obs.input_anchmc_file;
            label_1 = "gen-level PYTHIA 8";
            label_2 = "gen-level anchored MC LHC23a3";
        }

        TH1D * obs_hist_1 = get_1D_histogram(infile_1, histname_1);
        TH1D * obs_hist_2 = get_1D_histogram(infile_2, histname_2);
        TH1D * obs_hist_3;
        TH1D * obs_hist_4;
        if ( option == 1 ) {
            obs_hist_3 = get_1D_histogram(infile_3, histname_3);
            obs_hist_4 = get_1D_histogram(infile_4, histname_4);
        }

        // Calculate ratios
        TH1D * hist_obs_ratio = (TH1D *) obs_hist_2->Clone(Form("hratio_%s_pt%d-%d_ptrlbin%d", obs.name.c_str(), pt_min, pt_max, j));
        hist_obs_ratio->Divide(obs_hist_1);

        TH1D * hist_obs_ratio_31;
        TH1D * hist_obs_ratio_41;
        if ( option == 1 ) {
            hist_obs_ratio_31 = (TH1D *) obs_hist_3->Clone(Form("hratio31_%s_pt%d-%d_ptrlbin%d", obs.name.c_str(), pt_min, pt_max, j)); //31 means hist 3 / hist 1
            hist_obs_ratio_31->Divide(obs_hist_1);
            hist_obs_ratio_41 = (TH1D *) obs_hist_4->Clone(Form("hratio41_%s_pt%d-%d_ptrlbin%d", obs.name.c_str(), pt_min, pt_max, j)); //41 means hist 4 / hist 1
            hist_obs_ratio_41->Divide(obs_hist_1);
        }

        // Make canvas and legend
        TCanvas * can_obs = new TCanvas();
        can_obs->cd();
        gPad->SetLogy();
        TLegend *leg_obs = new TLegend(0.45, 0.75, 0.88, 0.85);
        TLegend *leg_ratio = new TLegend(0.7, 0.7, 0.88, 0.88);
        leg_obs->SetTextSize(0.027);

        // Format hists, and add to legend
        leg_obs_all->AddEntry("NULL", Form("p_{T}R_{L} = %.1f - %.1f", ptRL_bins[j], ptRL_bins[j+1]), "h");
        Format1DHist(obs, obs_hist_1, leg_obs_all, label_1, false, colors[j], kFullCircle);
        Format1DHist(obs, obs_hist_2, leg_obs_all, label_2, true, colors[j], 1, 0.6);
        if ( option == 1 ) {
            Format1DHist(obs, obs_hist_3, leg_obs_all, label_3, true, colors[j], 2, 0.6);
            Format1DHist(obs, obs_hist_4, leg_obs_all, label_4, true, colors[j], 3, 0.6);
        }
        Format1DHist(obs, hist_obs_ratio, leg_ratio_all, Form("PYTHIA / data p_{T}R_{L} = %.1f - %.1f", ptRL_bins[j], ptRL_bins[j+1]), true, colors[j], 1);
        if ( option == 1 ) {
            Format1DHist(obs, hist_obs_ratio_31, leg_ratio_all, Form("HERWIG / data p_{T}R_{L} = %.1f - %.1f", ptRL_bins[j], ptRL_bins[j+1]), true, colors[j], 2);
            Format1DHist(obs, hist_obs_ratio_41, leg_ratio_all, Form("LHC23a3 / data p_{T}R_{L} = %.1f - %.1f", ptRL_bins[j], ptRL_bins[j+1]), true, colors[j], 3);
        }
        leg_obs->AddEntry(obs_hist_1, Form("%s p_{T}R_{L} = %.1f - %.1f", label_1.c_str(), ptRL_bins[j], ptRL_bins[j+1]), "pl");
        leg_obs->AddEntry(obs_hist_2, Form("%s p_{T}R_{L} = %.1f - %.1f", label_2.c_str(), ptRL_bins[j], ptRL_bins[j+1]), "pl");
        if ( option == 1 ) {
            leg_obs->AddEntry(obs_hist_3, Form("%s p_{T}R_{L} = %.1f - %.1f", label_3.c_str(), ptRL_bins[j], ptRL_bins[j+1]), "pl");
            leg_obs->AddEntry(obs_hist_4, Form("%s p_{T}R_{L} = %.1f - %.1f", label_4.c_str(), ptRL_bins[j], ptRL_bins[j+1]), "pl");
        }
        leg_ratio->AddEntry(hist_obs_ratio, "ratio", "pl");
                

        // Draw individual pTRL bins
        obs_hist_2->Draw("HIST");
        if ( option == 1 ) {
            obs_hist_3->Draw("HIST SAME");
            obs_hist_4->Draw("HIST SAME");
        }
        obs_hist_1->Draw("P SAME");
        leg_obs->Draw();
        leg_ev->Draw();

        // Save individual pTRL bins
        can_obs->SaveAs(Form("%s/individuals/%s_pt%d-%d_ptrlbin%d_option%d.pdf", obs.filepath_plots.c_str(), obs.name.c_str(), pt_min, pt_max, j, option));
    
        // Now draw all pTRL bins together
        can_obs_all->cd();
        obs_hist_2->Draw("HIST SAME");
        if ( option == 1 ) {
            obs_hist_3->Draw("HIST SAME");
            obs_hist_4->Draw("HIST SAME");
        }
        obs_hist_1->Draw("P SAME");
        
    }

    // Save all pTRL bins together
    can_obs_all->cd();
    leg_obs_all->Draw();
    leg_ev->Draw();
    can_obs_all->SaveAs(Form("%s/%s_pt%d-%d_option%d.pdf", obs.filepath_plots.c_str(), obs.name.c_str(), pt_min, pt_max, option));
    

}

void plot_pt_bins(std::string weight_str, std::string jetR, std::string threshold, std::string norm_string) {

    const int pt_bins[] = { 20, 40, 60, 80 };
    const int n_bins = sizeof(pt_bins) / sizeof(pt_bins[0]) - 1;

    Observable obs_deltap("deltap", "#Deltap", "#frac{dN}{d#Deltap}"); 
    Observable obs_deltajt("deltajt", "#Deltaj_{T}", "#frac{dN}{d#Deltaj_{T}}");

    for ( int i = 0; i < n_bins; i++ ) {
        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];

        // option 1 = plot data (raw) vs pythia (det) vs herwig (det) vs anchMC (det)
        plot_1D_obs(obs_deltap, pt_min, pt_max, weight_str, jetR, threshold, norm_string, 1);
        plot_1D_obs(obs_deltajt, pt_min, pt_max, weight_str, jetR, threshold, norm_string, 1);

        // option 2 = plot pythia (gen) vs pythia (det)
        plot_1D_obs(obs_deltap, pt_min, pt_max, weight_str, jetR, threshold, norm_string, 2);
        plot_1D_obs(obs_deltajt, pt_min, pt_max, weight_str, jetR, threshold, norm_string, 2);

        // option 3 = plot herwig (gen) vs herwig (det)
        plot_1D_obs(obs_deltap, pt_min, pt_max, weight_str, jetR, threshold, norm_string, 3);
        plot_1D_obs(obs_deltajt, pt_min, pt_max, weight_str, jetR, threshold, norm_string, 3);

        // option 4 = plot pythia (det) vs anchMC (det)
        plot_1D_obs(obs_deltap, pt_min, pt_max, weight_str, jetR, threshold, norm_string, 4);
        plot_1D_obs(obs_deltajt, pt_min, pt_max, weight_str, jetR, threshold, norm_string, 4);

        // option 5 = plot pythia (gen) vs anchMC (gen)
        plot_1D_obs(obs_deltap, pt_min, pt_max, weight_str, jetR, threshold, norm_string, 5);
        plot_1D_obs(obs_deltajt, pt_min, pt_max, weight_str, jetR, threshold, norm_string, 5);
    }
}

void plot_pt_distribution() {
    Observable obs_jet_pt("jet_pt", "p_{T, jet}", "Counts");

    TCanvas * can_pt_all = new TCanvas();
    // Define pad heights (top bigger, bottom smaller)
    float padSplit = 0.3;  // fraction for bottom pad

    // Create top pad
    TPad *pad1 = new TPad("pad1", "pad1", 0, padSplit, 1, 1);
    pad1->SetBottomMargin(0); // remove bottom margin for top pad
    pad1->SetTicks(1,1);
    pad1->Draw();

    TLegend *leg_ev = new TLegend(0.5, 0.6, 0.6, 0.88);
    leg_ev->SetTextSize(0.03);
    leg_ev->AddEntry((TObject*)0, "pp, #sqrt{s} = 5.02 TeV", "");
    leg_ev->AddEntry((TObject*)0, "all ch. jets", "");
    leg_ev->AddEntry((TObject*)0, "anti-k_{T}, R = 0.4", "");


    TLegend *leg_pt_all = new TLegend(0.65, 0.5, 0.88, 0.88);
    TLegend *leg_ratio_all = new TLegend(0.7, 0.7, 0.88, 0.88);
    leg_pt_all->SetTextSize(0.03);

    // get all jet pt histograms
    TH1D * hist_jet_pt_data = get_1D_histogram(obs_jet_pt.input_data_file, "jet_pt_hist");
    TH1D * hist_jet_pt_pythia_truth = get_1D_histogram(obs_jet_pt.input_pythia_file, "jet_pt_hist_truth");
    TH1D * hist_jet_pt_pythia_det = get_1D_histogram(obs_jet_pt.input_pythia_file, "jet_pt_hist_det");
    TH1D * hist_jet_pt_herwig_truth = get_1D_histogram(obs_jet_pt.input_herwig_file, "jet_pt_hist_truth");
    TH1D * hist_jet_pt_herwig_det = get_1D_histogram(obs_jet_pt.input_herwig_file, "jet_pt_hist_det");
    TH1D * hist_jet_pt_anchmc_truth = get_1D_histogram(obs_jet_pt.input_anchmc_file, "jet_pt_hist_truth");
    TH1D * hist_jet_pt_anchmc_det = get_1D_histogram(obs_jet_pt.input_anchmc_file, "jet_pt_hist_det");

    // Make ratios
    TH1D * hist_ratio_pythia_truth_over_data = (TH1D *) hist_jet_pt_pythia_truth->Clone("hist_ratio_pythia_truth_over_data");
    TH1D * hist_ratio_pythia_det_over_data = (TH1D *) hist_jet_pt_pythia_det->Clone("hist_ratio_pythia_det_over_data");
    TH1D * hist_ratio_herwig_truth_over_data = (TH1D *) hist_jet_pt_herwig_truth->Clone("hist_ratio_herwig_truth_over_data");
    TH1D * hist_ratio_herwig_det_over_data = (TH1D *) hist_jet_pt_herwig_det->Clone("hist_ratio_herwig_det_over_data");
    TH1D * hist_ratio_anchmc_truth_over_data = (TH1D *) hist_jet_pt_anchmc_truth->Clone("hist_ratio_anchmc_truth_over_data");
    TH1D * hist_ratio_anchmc_det_over_data = (TH1D *) hist_jet_pt_anchmc_det->Clone("hist_ratio_anchmc_det_over_data");
    hist_ratio_pythia_truth_over_data->Divide(hist_jet_pt_data);
    hist_ratio_pythia_det_over_data->Divide(hist_jet_pt_data);
    hist_ratio_herwig_truth_over_data->Divide(hist_jet_pt_data);
    hist_ratio_herwig_det_over_data->Divide(hist_jet_pt_data);
    hist_ratio_anchmc_truth_over_data->Divide(hist_jet_pt_data);
    hist_ratio_anchmc_det_over_data->Divide(hist_jet_pt_data);

    // Format hists, and add to legend
    // leg_pt_all->AddEntry("NULL", Form("p_{T}R_{L} = %.1f - %.1f", ptRL_bins[j], ptRL_bins[j+1]), "h");
    Format1DHist(obs_jet_pt, hist_jet_pt_data, leg_pt_all, "Raw data", false, kBlack, kFullCircle);
    Format1DHist(obs_jet_pt, hist_jet_pt_pythia_truth, leg_pt_all, "truth-level PYTHIA 8", true, kBlue, 1); //, 1, 0.6);
    Format1DHist(obs_jet_pt, hist_jet_pt_pythia_det, leg_pt_all, "det-level PYTHIA 8", true, kRed, 1); //, 1, 0.6);
    Format1DHist(obs_jet_pt, hist_jet_pt_herwig_truth, leg_pt_all, "truth-level HERWIG 7", true, kGreen+2, 1); //, 2, 0.6);
    Format1DHist(obs_jet_pt, hist_jet_pt_herwig_det, leg_pt_all, "det-level HERWIG 7", true, kOrange+1, 1); //, 2, 0.6);
    Format1DHist(obs_jet_pt, hist_jet_pt_anchmc_truth, leg_pt_all, "truth-level anchored MC LHC23a3", true, kViolet-3, 1); //, 2, 0.6);
    Format1DHist(obs_jet_pt, hist_jet_pt_anchmc_det, leg_pt_all, "det-level anchored MC LHC23a3", true, kMagenta-4, 1); //, 2, 0.6);

    hist_jet_pt_data->GetXaxis()->SetTitle(obs_jet_pt.axis_label.c_str());
    hist_jet_pt_pythia_truth->GetXaxis()->SetTitle(obs_jet_pt.axis_label.c_str());
    hist_jet_pt_pythia_det->GetXaxis()->SetTitle(obs_jet_pt.axis_label.c_str());
    hist_jet_pt_herwig_truth->GetXaxis()->SetTitle(obs_jet_pt.axis_label.c_str());
    hist_jet_pt_herwig_det->GetXaxis()->SetTitle(obs_jet_pt.axis_label.c_str());
    hist_jet_pt_anchmc_truth->GetXaxis()->SetTitle(obs_jet_pt.axis_label.c_str());
    hist_jet_pt_anchmc_det->GetXaxis()->SetTitle(obs_jet_pt.axis_label.c_str());
    
    Format1DHist(obs_jet_pt, hist_ratio_pythia_truth_over_data, leg_ratio_all, "truth-level PYTHIA 8", true, kBlue, 1);
    Format1DHist(obs_jet_pt, hist_ratio_pythia_det_over_data, leg_ratio_all, "det-level PYTHIA 8", true, kRed, 1);
    Format1DHist(obs_jet_pt, hist_ratio_herwig_truth_over_data, leg_ratio_all, "truth-level HERWIG 7", true, kGreen+2, 1);
    Format1DHist(obs_jet_pt, hist_ratio_herwig_det_over_data, leg_ratio_all, "det-level HERWIG 7", true, kOrange+1, 1);
    Format1DHist(obs_jet_pt, hist_ratio_anchmc_truth_over_data, leg_ratio_all, "truth-level anchored MC LHC23a3", true, kViolet-3, 1);
    Format1DHist(obs_jet_pt, hist_ratio_anchmc_det_over_data, leg_ratio_all, "det-level anchored MC LHC23a3", true, kMagenta-4, 1);

    hist_ratio_pythia_truth_over_data->GetXaxis()->SetTitle(obs_jet_pt.axis_label.c_str());
    hist_ratio_pythia_det_over_data->GetXaxis()->SetTitle(obs_jet_pt.axis_label.c_str());
    hist_ratio_herwig_truth_over_data->GetXaxis()->SetTitle(obs_jet_pt.axis_label.c_str());
    hist_ratio_herwig_det_over_data->GetXaxis()->SetTitle(obs_jet_pt.axis_label.c_str());
    hist_ratio_anchmc_truth_over_data->GetXaxis()->SetTitle(obs_jet_pt.axis_label.c_str());
    hist_ratio_anchmc_det_over_data->GetXaxis()->SetTitle(obs_jet_pt.axis_label.c_str());
    // leg_ratio->AddEntry(hist_obs_ratio, "ratio", "pl");

    hist_ratio_pythia_truth_over_data->GetYaxis()->SetTitle("MC/Data");
    hist_ratio_pythia_det_over_data->GetYaxis()->SetTitle("MC/Data");
    hist_ratio_herwig_truth_over_data->GetYaxis()->SetTitle("MC/Data");
    hist_ratio_herwig_det_over_data->GetYaxis()->SetTitle("MC/Data");
    hist_ratio_anchmc_truth_over_data->GetYaxis()->SetTitle("MC/Data");
    hist_ratio_anchmc_det_over_data->GetYaxis()->SetTitle("MC/Data");


    // Draw individual pTRL bins
    pad1->cd();
    gPad->SetLogy();
    hist_jet_pt_data->Draw("PE");
    hist_jet_pt_pythia_truth->Draw("L SAME");
    hist_jet_pt_pythia_det->Draw("L SAME");
    hist_jet_pt_herwig_truth->Draw("L SAME");
    hist_jet_pt_herwig_det->Draw("L SAME");
    hist_jet_pt_anchmc_truth->Draw("L SAME");
    hist_jet_pt_anchmc_det->Draw("L SAME");
    leg_ev->Draw();
    leg_pt_all->Draw();
    

    // Make bottom ratio panel
    can_pt_all->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, padSplit);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.28); // leave room for x-axis labels
    pad2->SetTicks(1,1);
    pad2->Draw();

    pad2->cd();
    gPad->SetLogy();
    hist_ratio_pythia_truth_over_data->SetMinimum(0.01);
    hist_ratio_pythia_truth_over_data->Draw("L");
    hist_ratio_pythia_det_over_data->Draw("L SAME");
    hist_ratio_herwig_truth_over_data->Draw("L SAME");
    hist_ratio_herwig_det_over_data->Draw("L SAME");
    hist_ratio_anchmc_truth_over_data->Draw("L SAME");
    hist_ratio_anchmc_det_over_data->Draw("L SAME");
    leg_ratio_all->Draw();

    // Save individual pTRL bins
    can_pt_all->SaveAs(Form("%s/%s_all.pdf", obs_jet_pt.filepath_plots.c_str(), obs_jet_pt.name.c_str()));
    
}


void plot_data_with_generators() {
    // TFile * file_data_hists = new TFile("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/rootfiles/data_fifthattempt_ptrlbins/DataHists_%s.root", "READ");
    // TFile * pythia_data_hists = new TFile("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/rootfiles/pythia5TeV_histograms_crosscheck/PYTHIAHists_%s.root", "READ");

    SetStyle();

    plot_pt_distribution();

    std::string weight_str = "";
    std::string jetR = "0.4";
    std::string threshold = "1.0";
    std::string norm_string = "self_normalized";
    plot_pt_bins(weight_str, jetR, threshold, norm_string);
}