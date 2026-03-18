// i want to compare generator and detector level HF pythia and herwig
// pp 13 TeV
// there are four configurations: pythia_prompt, pythia_nonprompt, herwig_prompt, herwig_nonprompt
// are there differences in herwig and pythia, even at particle level?
// branches: ParticlePt, ParticleEta, ParticlePhi, ParticlePID
// this takes the already made histograms and plots them

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <typeinfo>

std::string base_filepath_header = "/global/cfs/projectdirs/alice/alicepro/hiccup";

TLine * drawHoriLine(double x1, double x2, double y1, int color, int linestyle=2){
    auto fhoriline = new TLine(x1, y1, x2, y1);
	fhoriline->SetLineWidth(1);
    fhoriline->SetLineColor(color);
    fhoriline->SetLineStyle(linestyle);
    return fhoriline;
}


class Observable {
public:
    std::string name;
    std::string numtype;
    bool logy;
    std::string xtitle;
    bool cs; // true if the cross section is included in the y axis
    bool cs_method2;
    bool onlycharged;

    std::string ytitle;

    double min_xval = 0;
    double max_xval = -1;
    double min_yval = -1;
    double max_yval_ratio = -1;

    TH1 * hist_pythia_prompt;
    TH1 * hist_pythia_nonprompt;
    TH1 * hist_herwig_prompt;
    TH1 * hist_herwig_nonprompt;

    TH1 * hist_pythia_notforcedD0toKPi_prompt;
    TH1 * hist_pythia_notforcedD0toKPi_nonprompt;

    Observable(std::string name_val, std::string numtype_val, bool logy_val, std::string xtitle_val, bool cs_val, bool cs_method2_val=false, bool onlycharged_val=false) {
        name = name_val;
        numtype = numtype_val;
        logy = logy_val;
        xtitle = xtitle_val;
        cs = cs_val;
        cs_method2 = cs_method2_val;
        onlycharged = onlycharged_val;

        assign_specific_bounds();
    }

    void assign_specific_bounds() {
        if (name == "Pt") {
            max_xval = 100;
            max_yval_ratio = 2.25;
            if (cs) max_yval_ratio = 10; //1e-5;
        } else if (name == "D0_Pt") {
            max_xval = 60;
            max_yval_ratio = 50;
            if (cs) max_yval_ratio = 150;
        } else if (name == "Phi") {
            max_yval_ratio = 3.5;
        } else if (name == "D0_Phi") {
            max_yval_ratio = 30;
        } else if (name == "Eta") {
            max_yval_ratio = 4;
        } else if (name == "D0_Eta" || name == "D0_Rap") {
            min_yval = 2e3;
            max_yval_ratio = 25; //idk why this height is showing as double
        } else if (name == "PID") {
            max_yval_ratio = 50; 
        } else if (name == "D0_MPID") {
            min_xval = -600;
            max_xval = 600;
            max_yval_ratio = 50; 
        } else if (name == "numD0s") {
            min_yval = 1e2;
            max_yval_ratio = 40;
        } 
    }

    void style_hists() {
        hist_pythia_prompt->SetMarkerColorAlpha(kBlue, 1.0);
        hist_pythia_prompt->SetLineColorAlpha(kBlue, 1.0);
        hist_pythia_prompt->SetMarkerStyle(kFullCircle);

        hist_pythia_nonprompt->SetMarkerColorAlpha(kBlue-7, 1.0);
        hist_pythia_nonprompt->SetLineColorAlpha(kBlue-7, 1.0);
        hist_pythia_nonprompt->SetMarkerStyle(kOpenCircle);
        hist_pythia_nonprompt->SetLineStyle(9);

        hist_herwig_prompt->SetMarkerColorAlpha(kRed, 1.0);
        hist_herwig_prompt->SetLineColorAlpha(kRed, 1.0);
        hist_herwig_prompt->SetMarkerStyle(kFullCircle);

        hist_herwig_nonprompt->SetMarkerColorAlpha(kRed-7, 1.0);
        hist_herwig_nonprompt->SetLineColorAlpha(kRed-7, 1.0);
        hist_herwig_nonprompt->SetMarkerStyle(kOpenCircle);
        hist_herwig_nonprompt->SetLineStyle(9);

        hist_pythia_notforcedD0toKPi_prompt->SetMarkerColorAlpha(kOrange-3, 1.0);
        hist_pythia_notforcedD0toKPi_prompt->SetLineColorAlpha(kOrange-3, 1.0);
        hist_pythia_notforcedD0toKPi_prompt->SetMarkerStyle(kFullCircle);

        hist_pythia_notforcedD0toKPi_nonprompt->SetMarkerColorAlpha(kGreen-3, 1.0);
        hist_pythia_notforcedD0toKPi_nonprompt->SetLineColorAlpha(kGreen-3, 1.0);
        hist_pythia_notforcedD0toKPi_nonprompt->SetMarkerStyle(kOpenCircle);
        hist_pythia_notforcedD0toKPi_nonprompt->SetLineStyle(9);

        if ( cs or cs_method2 or onlycharged ) ytitle = "#frac{d#sigma}{d" + xtitle + "}";  
        else ytitle = "#frac{dN}{d" + xtitle + "}";
              
        hist_pythia_prompt->GetYaxis()->SetTitle(ytitle.c_str());
        hist_pythia_nonprompt->GetYaxis()->SetTitle(ytitle.c_str());
        hist_herwig_prompt->GetYaxis()->SetTitle(ytitle.c_str());
        hist_herwig_nonprompt->GetYaxis()->SetTitle(ytitle.c_str());
    }

    // -------- DRAW --------
    void drawPairForObs(bool normalized=false, bool add_third_curve=false) {
    
        std::string c_name = Form("c%s", name.c_str());
        TCanvas *c = new TCanvas(c_name.c_str(), c_name.c_str(), 800, 600);

        // --- pads ---
        TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
        TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.3);

        pad1->SetBottomMargin(0.0); //0.02);
        pad2->SetTopMargin(0.0); //0.02);
        pad2->SetBottomMargin(0.3);

        pad1->Draw();
        pad2->Draw();

        // ======================
        // Top pad: histograms
        // ======================
        pad1->cd();
        if ( logy ) gPad->SetLogy();

        if (normalized) {
            cout << "sigmas! " << endl;
            cout << "hist pythia prompt " << hist_pythia_prompt->Integral() << endl;
            cout << "hist pythia nonprompt " << hist_pythia_nonprompt->Integral() << endl;
            cout << "hist herwig prompt " << hist_herwig_prompt->Integral() << endl;
            cout << "hist herwig nonprompt " << hist_herwig_nonprompt->Integral() << endl;
            
            hist_pythia_prompt->Scale(1.0 / hist_pythia_prompt->Integral());
            hist_pythia_nonprompt->Scale(1.0 / hist_pythia_nonprompt->Integral());
            hist_herwig_prompt->Scale(1.0 / hist_herwig_prompt->Integral());
            hist_herwig_nonprompt->Scale(1.0 / hist_herwig_nonprompt->Integral());

            if (add_third_curve) {
                hist_pythia_notforcedD0toKPi_prompt->Scale(1.0 / hist_pythia_notforcedD0toKPi_prompt->Integral());
                hist_pythia_notforcedD0toKPi_nonprompt->Scale(1.0 / hist_pythia_notforcedD0toKPi_nonprompt->Integral()); 
            }

            hist_pythia_prompt->GetYaxis()->SetTitle(("#frac{1}{#sigma} " + ytitle).c_str());
        }

        // determine max
        double max_height = std::max( {hist_pythia_prompt->GetMaximum(), hist_pythia_nonprompt->GetMaximum(), hist_herwig_prompt->GetMaximum(), hist_herwig_nonprompt->GetMaximum()} );
        hist_pythia_prompt->SetMaximum(max_height*1.2);
        if (min_yval > 0) hist_pythia_prompt->SetMinimum(min_yval);

        hist_pythia_prompt->Draw("hist");
        hist_pythia_nonprompt->Draw("hist same");
        hist_herwig_prompt->Draw("hist same");
        hist_herwig_nonprompt->Draw("hist same");
        if (add_third_curve) {
            hist_pythia_notforcedD0toKPi_prompt->Draw("hist same");
            hist_pythia_notforcedD0toKPi_nonprompt->Draw("hist same"); 
        }

        TLegend *leg = new TLegend(0.55, 0.75, 0.8, 0.88); //0.15, 0.75, 0.4, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(hist_pythia_prompt, "Pythia prompt", "l");
        leg->AddEntry(hist_pythia_nonprompt, "Pythia non-prompt", "l");
        leg->AddEntry(hist_herwig_prompt, "Herwig prompt", "l");
        leg->AddEntry(hist_herwig_nonprompt, "Herwig non-prompt", "l");
        if (add_third_curve) {
            leg->AddEntry(hist_pythia_notforcedD0toKPi_prompt, "Pythia prompt no D0->Kpi forcing", "l");
            leg->AddEntry(hist_pythia_notforcedD0toKPi_nonprompt, "Pythia non-prompt no D0->Kpi forcing", "l");
        }
        leg->Draw();

        // ======================
        // Bottom pad: ratio
        // ======================
        pad2->cd();

        if ( cs ) gPad->SetLogy();

        TH1D *h_ratio_prompt = (TH1D *)hist_pythia_prompt->Clone(Form("h_ratio_%s_prompt", name.c_str()));
        h_ratio_prompt->Divide(hist_herwig_prompt);
        TH1D *h_ratio_nonprompt = (TH1D *)hist_pythia_nonprompt->Clone(Form("h_ratio_%s_nonprompt", name.c_str()));
        h_ratio_nonprompt->Divide(hist_herwig_nonprompt);

        TH1D *h_ratio_notforcedD0toKPi_prompt;
        TH1D *h_ratio_notforcedD0toKPi_nonprompt;
        if (add_third_curve) {
            h_ratio_notforcedD0toKPi_prompt = (TH1D *)hist_pythia_notforcedD0toKPi_prompt->Clone(Form("h_ratio_%s_notforcedD0toKPi_prompt", name.c_str()));
            h_ratio_notforcedD0toKPi_prompt->Divide(hist_herwig_prompt);
            h_ratio_notforcedD0toKPi_nonprompt = (TH1D *)hist_pythia_notforcedD0toKPi_nonprompt->Clone(Form("h_ratio_%s_notforcedD0toKPi_nonprompt", name.c_str()));
            h_ratio_notforcedD0toKPi_nonprompt->Divide(hist_herwig_nonprompt);
        }

        // h_ratio->SetTitle(Form("Ratio %s", name.c_str()));
        h_ratio_prompt->GetYaxis()->SetTitle("PYTHIA / HERWIG");
        h_ratio_prompt->GetYaxis()->SetNdivisions(505);
        h_ratio_prompt->GetYaxis()->SetTitleSize(0.10);
        h_ratio_prompt->GetYaxis()->SetLabelSize(0.08);
        h_ratio_prompt->GetYaxis()->SetTitleOffset(0.5);

        h_ratio_prompt->GetXaxis()->SetTitle(h_ratio_prompt->GetXaxis()->GetTitle());
        h_ratio_prompt->GetXaxis()->SetTitleSize(0.12);
        h_ratio_prompt->GetXaxis()->SetLabelSize(0.10);

        h_ratio_prompt->SetStats(0);
        h_ratio_nonprompt->SetStats(0);
        if (add_third_curve) {
            h_ratio_notforcedD0toKPi_prompt->SetStats(0);
            h_ratio_notforcedD0toKPi_nonprompt->SetStats(0);
        }

        if (max_yval_ratio > 0 ) h_ratio_prompt->SetMaximum(max_yval_ratio);
        if (name == "D0_Pt") h_ratio_prompt->SetMinimum(1e-1);
        if (normalized && name == "Pt") h_ratio_prompt->SetMaximum(1.2);
        if (normalized && name == "D0_Pt") h_ratio_prompt->SetMaximum(2.0);

        h_ratio_prompt->Draw("hist");
        h_ratio_nonprompt->Draw("hist same");
        if (add_third_curve) {
            h_ratio_notforcedD0toKPi_prompt->Draw("hist same");
            h_ratio_notforcedD0toKPi_nonprompt->Draw("hist same");
        }
        if (min_xval != 0 || max_xval != -1) drawHoriLine(min_xval, max_xval, 1, kGray+2)->Draw("same");

        TLegend *leg_ratio = new TLegend(0.78, 0.65, 0.93, 0.9);
        leg_ratio->SetBorderSize(0);
        leg_ratio->SetFillStyle(0);
        leg_ratio->AddEntry(h_ratio_prompt, "prompt", "l");
        leg_ratio->AddEntry(h_ratio_nonprompt, "non-prompt", "l");
        if (add_third_curve) {
            leg_ratio->AddEntry(h_ratio_notforcedD0toKPi_prompt, "prompt no forced D0->Kpi", "l");
            leg_ratio->AddEntry(h_ratio_notforcedD0toKPi_nonprompt, "non-prompt no forced D0->Kpi", "l");
        }
        leg_ratio->Draw();

        std::string file_plot_output = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/HF_EEC/plots/HF_particle_comparisons/" + name + "_comparison.pdf";
        if ( cs ) {
            file_plot_output = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/HF_EEC/plots/HF_particle_comparisons/" + name + "_crosssection_comparison.pdf";
            if ( normalized ) file_plot_output = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/HF_EEC/plots/HF_particle_comparisons/" + name + "_crosssection_normalized_comparison.pdf";
        }
        if ( cs_method2 ) file_plot_output = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/HF_EEC/plots/HF_particle_comparisons/" + name + "_crosssection_method2_comparison.pdf";
        if ( onlycharged ) file_plot_output = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/HF_EEC/plots/HF_particle_comparisons/" + name + "_onlycharged_crosssection_comparison.pdf";
        c->SaveAs(file_plot_output.c_str());
    }
    
};



// template <typename TH1> 
TH1 * read_histogram( TFile * fin, Observable obs, std::string gen_choice, std::string p_or_np, std::string gen_or_det ) {
    // static_assert(std::is_base_of<TH1, TH>::value, "TH must inherit from TH1");

    std::string histname = Form( "h%s_%s_%s_%s", obs.name.c_str(), gen_choice.c_str(), p_or_np.c_str(), gen_or_det.c_str() );
    if (obs.cs) histname = Form( "h%s_%s_%s_%s_crosssection", obs.name.c_str(), gen_choice.c_str(), p_or_np.c_str(), gen_or_det.c_str() );
    if (obs.cs_method2) histname = Form( "h%s_%s_%s_%s_crosssection_method2", obs.name.c_str(), gen_choice.c_str(), p_or_np.c_str(), gen_or_det.c_str() );
    if (obs.onlycharged && gen_choice == "pythia") Form( "h%s_onlycharged_%s_%s_%s_crosssection", obs.name.c_str(), gen_choice.c_str(), p_or_np.c_str(), gen_or_det.c_str() );

    TH1 * hist = dynamic_cast<TH1 *>(fin->Get(histname.c_str()));
    // TH1 * hist = (TH1 *) fin->Get(histname.c_str());
    cout << "got histogram!" << endl;
    if (!hist) {
        Error("read_histogram", "Histogram %s not found or wrong type", histname.c_str());
        return nullptr;
    }

    if ( obs.max_xval > 0 ) hist->GetXaxis()->SetRangeUser(0, obs.max_xval); // i guess this is only applicable for pt right now

    return hist;
}

void draw_multiple_testing(std::vector<Observable> obs_vec) {

    std::string c_name = Form("c%s_multiple", obs_vec[0].name.c_str());
    TCanvas *c = new TCanvas(c_name.c_str(), c_name.c_str(), 800, 600);
    gPad->SetLogy();

    // determine max
    double max_height = 0;
    for ( int i = 0; i < obs_vec.size(); i++ ) {
        double temp_max = std::max( {obs_vec[i].hist_pythia_prompt->GetMaximum(), obs_vec[i].hist_pythia_nonprompt->GetMaximum(), obs_vec[i].hist_herwig_prompt->GetMaximum(), obs_vec[i].hist_herwig_nonprompt->GetMaximum()} );
        max_height = std::max(max_height, temp_max);
    }
    obs_vec[0].hist_pythia_prompt->SetMaximum(max_height*1.2);
    obs_vec[0].hist_pythia_prompt->SetMinimum(0.);
    // obs_vec[1].hist_pythia_prompt->SetMaximum(max_height*1.2);
    // obs_vec[1].hist_pythia_prompt->SetMinimum(0.);

    TLegend *leg = new TLegend(0.15, 0.75, 0.4, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    Double_t colors[12] = {kBlue, kBlue-7, kRed, kRed-7, kOrange, kOrange-4, kGreen, kGreen-7, kMagenta, kMagenta-9, kCyan, kCyan-9};
    
    for ( int i = 0; i < obs_vec.size(); i++ ) {
        if ( i == 0 ) continue;
        obs_vec[i].hist_pythia_prompt->SetLineColorAlpha(colors[4*i], 0.5);
        obs_vec[i].hist_pythia_nonprompt->SetLineColorAlpha(colors[4*i+1], 0.5);
        obs_vec[i].hist_herwig_prompt->SetLineColorAlpha(colors[4*i+2], 0.5);
        obs_vec[i].hist_herwig_nonprompt->SetLineColorAlpha(colors[4*i+3], 0.5);
        
        obs_vec[i].hist_pythia_prompt->Draw("hist same");
        obs_vec[i].hist_pythia_nonprompt->Draw("hist same");
        obs_vec[i].hist_herwig_prompt->Draw("hist same");
        obs_vec[i].hist_herwig_nonprompt->Draw("hist same");

        leg->AddEntry(obs_vec[i].hist_pythia_prompt, "Pythia prompt", "l");
        leg->AddEntry(obs_vec[i].hist_pythia_nonprompt, "Pythia non-prompt", "l");
        leg->AddEntry(obs_vec[i].hist_herwig_prompt, "Herwig prompt", "l");
        leg->AddEntry(obs_vec[i].hist_herwig_nonprompt, "Herwig non-prompt", "l");
    }
    leg->Draw();
    std::string file_plot_output = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/HF_EEC/plots/HF_particle_comparisons/TESTING_" + obs_vec[0].name + "_multiple.pdf";
    c->SaveAs(file_plot_output.c_str());
}

void comparing_pythia_pt_spectra() {

}


void analyze(TFile * fin_pythia, TFile * fin_herwig, TFile * fin_pythia_notforcedD0toKPi = nullptr) {

    bool add_third_curve = (fin_pythia_notforcedD0toKPi != nullptr);

    // obs = Pt, Eta, Phi, PID
    // obs = D0_Pt, D0_Eta, D0_Phi, D0_Rap, D0_MPID, numD0s
    Observable obs_pt("Pt", "double", true, "p_{T}", false);
    Observable obs_eta("Eta", "double", true, "#eta", false);
    Observable obs_phi("Phi", "double", true, "#phi", false);
    Observable obs_pid("PID", "int", true, "PID", false);

    Observable obs_D0_pt("D0_Pt", "double", true, "p_{T}", false);
    Observable obs_D0_eta("D0_Eta", "double", true, "#eta", false);
    Observable obs_D0_phi("D0_Phi", "double", false, "#phi", false);
    Observable obs_D0_rap("D0_Rap", "double", true, "y", false);
    Observable obs_D0_mpid("D0_MPID", "int", true, "Mother PID", false);
    Observable obs_numD0s("numD0s", "int", true, "# D0s per event", false);

    Observable obs_pt_cs("Pt", "double", true, "p_{T}", true);
    Observable obs_D0_pt_cs("D0_Pt", "double", true, "p_{T}", true);
    Observable obs_D0_rap_cs("D0_Rap", "double", true, "y", true);

    // Observable obs_pt_cs_method2("Pt", "double", true, "p_{T}", false, true);
    // Observable obs_D0_pt_cs_method2("D0_Pt", "double", true, "p_{T}", false, true);
    // Observable obs_D0_rap_cs_method2("D0_Rap", "double", true, "y", false, true);

    // Observable obs_pt_onlycharged_cs("Pt", "double", true, "p_{T}", false, false, true);

    Observable* obs_list[13] = { &obs_pt, &obs_eta, &obs_phi, &obs_pid, &obs_D0_pt, &obs_D0_eta, &obs_D0_phi, &obs_D0_rap, &obs_D0_mpid, &obs_numD0s, &obs_pt_cs, &obs_D0_pt_cs, &obs_D0_rap_cs }; //, &obs_pt_onlycharged_cs };

    for ( Observable* obs : obs_list ) {
        cout << "Running observable " << obs->name << endl;

        obs->hist_pythia_prompt = read_histogram(fin_pythia, *obs, "pythia", "prompt", "gen");
        obs->hist_pythia_nonprompt = read_histogram(fin_pythia, *obs, "pythia", "nonprompt", "gen");
        obs->hist_herwig_prompt = read_histogram(fin_herwig, *obs, "herwig", "prompt", "gen");
        obs->hist_herwig_nonprompt = read_histogram(fin_herwig, *obs, "herwig", "nonprompt", "gen");

        if (add_third_curve) {
            obs->hist_pythia_notforcedD0toKPi_prompt = read_histogram(fin_pythia_notforcedD0toKPi, *obs, "pythia", "prompt", "gen");
            obs->hist_pythia_notforcedD0toKPi_nonprompt = read_histogram(fin_pythia_notforcedD0toKPi, *obs, "pythia", "nonprompt", "gen");
        }

        cout << "styling hists now" << endl;
        obs->style_hists();
        cout << "drawing hists now" << endl;
        if (add_third_curve) obs->drawPairForObs(false, true);
        else obs->drawPairForObs();
        if (obs->cs) {
            if (add_third_curve) obs->drawPairForObs(true, true); // normalize by 1/sigma
            else obs->drawPairForObs(true); // normalize by 1/sigma
        }
    }

    /* // this is debugging method!
    std::vector<Observable> obs_pt_vec = { obs_pt, obs_pt_cs, obs_pt_cs_method2 }; // could use pointers + addresses here, but I don't want to update anything in Observable
    draw_multiple_testing(obs_pt_vec);
    std::vector<Observable> obs_D0_pt_vec = { obs_D0_pt, obs_D0_pt_cs, obs_D0_pt_cs_method2 };
    draw_multiple_testing(obs_D0_pt_vec);
    */
}

// ============ CHECKING NUM OF ENTRIES ============
class GenEntries {
public:
    std::string gen_name;
    std::string txtfile;

    long long num_particles_prompt;
    long long num_particles_nonprompt;
    long long num_D0s_prompt;
    long long num_D0s_nonprompt;

    long long num_particles_prompt_arr[10];
    long long num_particles_nonprompt_arr[10];
    long long num_D0s_prompt_arr[10];
    long long num_D0s_nonprompt_arr[10];


    GenEntries(std::string gen_name_val, std::string txtfile_val) {
        gen_name = gen_name_val;
        txtfile = txtfile_val;
    }
    
};

// Helper function to get the stream of numbers after the ":"
std::stringstream getNumbersStream(std::ifstream& file) {
    std::string line;
    if (std::getline(file, line)) {
        size_t colonPos = line.find(':');
        if (colonPos != std::string::npos) {
            // Return a stream starting after the colon
            return std::stringstream(line.substr(colonPos + 1));
        }
    }
    return std::stringstream(""); // Return empty stream if failed
}

void read_num_entries_txtfile(GenEntries& genentries) {
    std::ifstream infile(genentries.txtfile);

    // first 4 files
    getNumbersStream(infile) >> genentries.num_particles_prompt;
    getNumbersStream(infile) >> genentries.num_particles_nonprompt;
    getNumbersStream(infile) >> genentries.num_D0s_prompt;
    getNumbersStream(infile) >> genentries.num_D0s_nonprompt;

    // next 4 lines
    std::stringstream ss = getNumbersStream(infile);
    for ( int i = 0; i < 10; i++ ) ss >> genentries.num_particles_prompt_arr[i];
    ss = getNumbersStream(infile);
    for ( int i = 0; i < 10; i++ ) ss >> genentries.num_particles_nonprompt_arr[i];
    ss = getNumbersStream(infile);
    for ( int i = 0; i < 10; i++ ) ss >> genentries.num_D0s_prompt_arr[i];
    ss = getNumbersStream(infile);
    for ( int i = 0; i < 10; i++ ) ss >> genentries.num_D0s_nonprompt_arr[i];

    infile.close();

}

void check_numbers(GenEntries genentries) {
    long long sum_of_ptbins = 0;
    for ( int i = 0; i < 10; i++ ) sum_of_ptbins += genentries.num_particles_prompt_arr[i];
    if ( sum_of_ptbins == genentries.num_particles_prompt ) cout << genentries.gen_name <<  " prompt particles confirmed at " << sum_of_ptbins << endl;

    sum_of_ptbins = 0;
    for ( int i = 0; i < 10; i++ ) sum_of_ptbins += genentries.num_particles_nonprompt_arr[i];
    if ( sum_of_ptbins == genentries.num_particles_nonprompt ) cout << genentries.gen_name <<  " nonprompt particles confirmed at " << sum_of_ptbins << endl;

    sum_of_ptbins = 0;
    for ( int i = 0; i < 10; i++ ) sum_of_ptbins += genentries.num_D0s_prompt_arr[i];
    if ( sum_of_ptbins == genentries.num_D0s_prompt ) cout << genentries.gen_name <<  " prompt D0s confirmed at " << sum_of_ptbins << endl;

    sum_of_ptbins = 0;
    for ( int i = 0; i < 10; i++ ) sum_of_ptbins += genentries.num_D0s_nonprompt_arr[i];
    if ( sum_of_ptbins == genentries.num_D0s_nonprompt ) cout << genentries.gen_name <<  " nonprompt D0s confirmed at " << sum_of_ptbins << endl;

}

void check_num_entries() {
    std::string pythia_numentries_txt = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/HF_EEC/plots/HF_particle_comparisons/number_of_entries_pythia.txt";
    std::string herwig_numentries_txt = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/HF_EEC/plots/HF_particle_comparisons/number_of_entries_herwig.txt";

    GenEntries genentries_pythia("pythia", pythia_numentries_txt);
    GenEntries genentries_herwig("herwig", herwig_numentries_txt);

    read_num_entries_txtfile(genentries_pythia);
    read_num_entries_txtfile(genentries_herwig);
    
    check_numbers(genentries_pythia);
    check_numbers(genentries_herwig);

}

// ================================================

void compare_mc_HF_particle_level() {

    // -------- INPUT HISTOGRAMS -------- -- no det level at this point in time
    TString pythia_hists = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/HF_EEC/rootfiles/HF_particle_comparisons/HF_particle_comparisons_pythia.root"; 
    TString herwig_hists = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/HF_EEC/rootfiles/HF_particle_comparisons/HF_particle_comparisons_herwig.root";
    TString pythia_notforcedD0toKPi_hists = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/HF_EEC/rootfiles/HF_particle_comparisons/HF_particle_comparisons_pythia_notforcedD0toKPi.root"; 

    // -------- OPEN FILES --------
    TFile * file_pythia_hists = new TFile(pythia_hists, "READ");
    TFile * file_herwig_hists = new TFile(herwig_hists, "READ");
    TFile * file_pythia_notforcedD0toKPi_hists = new TFile(pythia_notforcedD0toKPi_hists, "READ");

    // -------- PLOT --------
    // analyze(file_pythia_hists, file_herwig_hists);
    analyze(file_pythia_hists, file_herwig_hists, file_pythia_notforcedD0toKPi_hists);

    file_pythia_hists->Close();
    file_herwig_hists->Close();
    file_pythia_notforcedD0toKPi_hists->Close();

    // -------- DOUBLE CHECK # OF ENTRIES --------
    check_num_entries();
}


