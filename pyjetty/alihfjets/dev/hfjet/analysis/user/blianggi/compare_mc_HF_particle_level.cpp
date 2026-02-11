// i want to compare generator and detector level HF pythia and herwig
// pp 13 TeV
// there are four configurations: pythia_prompt, pythia_nonprompt, herwig_prompt, herwig_nonprompt
// are there differences in herwig and pythia, even at particle level?
// branches: ParticlePt, ParticleEta, ParticlePhi, ParticlePID
// this takes the already made histograms and plots them


std::string base_filepath_header = "/global/cfs/projectdirs/alice/alicepro/hiccup";


class Observable {
public:
    std::string name;
    std::string numtype;
    bool logy;
    std::string xtitle;

    std::string ytitle_entries;
    std::string ytitle_crosssection;

    double max_xval = -1;
    double min_yval = -1;
    double max_yval_ratio = -1;

    TH1 * hist_pythia_prompt;
    TH1 * hist_pythia_nonprompt;
    TH1 * hist_herwig_prompt;
    TH1 * hist_herwig_nonprompt;

    Observable(std::string name_val, std::string numtype_val, bool logy_val, std::string xtitle_val) {
        name = name_val;
        numtype = numtype_val;
        logy = logy_val;
        xtitle = xtitle_val;

        ytitle_entries = "#frac{dN}{d" + xtitle + "}";
        ytitle_crosssection = "#frac{d#sigma}{d" + xtitle + "}";

        assign_specific_bounds();
    }

    void assign_specific_bounds() {
        if (name == "Pt") {
            max_xval = 100;
            max_yval_ratio = 2.25;
        } else if (name == "D0_Pt") {
            max_xval = 60;
            max_yval_ratio = 50;
        } else if (name == "Phi") {
            max_yval_ratio = 3.5;
        } else if (name == "D0_Phi") {
            max_yval_ratio = 30;
        } else if (name == "Eta") {
            max_yval_ratio = 4;
        } else if (name == "D0_Eta" || name == "D0_Rap") {
            min_yval = 2e3;
            max_yval_ratio = 25; //idk why this height is showing as double
        } else if (name == "PID" || name == "D0_MPID") {
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

        hist_pythia_prompt->GetYaxis()->SetTitle(ytitle_entries.c_str());
        hist_pythia_nonprompt->GetYaxis()->SetTitle(ytitle_entries.c_str());
        hist_herwig_prompt->GetYaxis()->SetTitle(ytitle_entries.c_str());
        hist_herwig_nonprompt->GetYaxis()->SetTitle(ytitle_entries.c_str());
    }

    // -------- DRAW --------
    void drawPairForObs() {
    
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

        // determine max
        double max_height = std::max( {hist_pythia_prompt->GetMaximum(), hist_pythia_nonprompt->GetMaximum(), hist_herwig_prompt->GetMaximum(), hist_herwig_nonprompt->GetMaximum()} );
        hist_pythia_prompt->SetMaximum(max_height*1.2);
        if (min_yval > 0) hist_pythia_prompt->SetMinimum(min_yval);

        hist_pythia_prompt->Draw("hist");
        hist_pythia_nonprompt->Draw("hist same");
        hist_herwig_prompt->Draw("hist same");
        hist_herwig_nonprompt->Draw("hist same");

        TLegend *leg = new TLegend(0.15, 0.75, 0.4, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(hist_pythia_prompt, "Pythia prompt", "l");
        leg->AddEntry(hist_pythia_nonprompt, "Pythia non-prompt", "l");
        leg->AddEntry(hist_herwig_prompt, "Herwig prompt", "l");
        leg->AddEntry(hist_herwig_nonprompt, "Herwig non-prompt", "l");
        leg->Draw();

        // ======================
        // Bottom pad: ratio
        // ======================
        pad2->cd();

        TH1D *h_ratio_prompt = (TH1D *)hist_pythia_prompt->Clone(Form("h_ratio_%s_prompt", name.c_str()));
        h_ratio_prompt->Divide(hist_herwig_prompt);
        TH1D *h_ratio_nonprompt = (TH1D *)hist_pythia_nonprompt->Clone(Form("h_ratio_%s_nonprompt", name.c_str()));
        h_ratio_nonprompt->Divide(hist_herwig_nonprompt);

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

        if (max_yval_ratio > 0 ) h_ratio_prompt->SetMaximum(max_yval_ratio);

        h_ratio_prompt->Draw("hist");
        h_ratio_nonprompt->Draw("hist same");

        TLegend *leg_ratio = new TLegend(0.8, 0.42, 0.9, 0.5);
        leg_ratio->SetBorderSize(0);
        leg_ratio->SetFillStyle(0);
        leg_ratio->AddEntry(h_ratio_prompt, "prompt", "l");
        leg_ratio->AddEntry(h_ratio_nonprompt, "non-prompt", "l");
        leg_ratio->Draw();

        c->SaveAs(Form("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/HF_EEC/plots/HF_particle_comparisons/%s_comparison.pdf", name.c_str()));
    }
    
};



// template <typename TH1> 
TH1 * read_histogram( TFile * fin, Observable obs, std::string gen_choice, std::string p_or_np, std::string gen_or_det ) {
    // static_assert(std::is_base_of<TH1, TH>::value, "TH must inherit from TH1");

    std::string histname = Form( "h%s_%s_%s_%s", obs.name.c_str(), gen_choice.c_str(), p_or_np.c_str(), gen_or_det.c_str() );

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







    // drawPair(fout_root, gen1, gen2, hPt_1,  hPt_2, "Pt"  max_yval_ratio = 2.25;+ gen1.gen_or_det);
    // drawPair(fout_root, gen1, gen2, hEta_1, hEta_2, "Eta" + gen1.gen_or_det);
    // drawPair(fout_root, gen1, gen2, hPhi_1, hPhi_2, "Phi" + gen1.gen_or_det);


void analyze(TFile * fin_pythia, TFile * fin_herwig) {
    // obs = Pt, Eta, Phi, PID
    // obs = D0_Pt, D0_Eta, D0_Phi, D0_Rap, D0_MPID, numD0s
    Observable obs_pt("Pt", "double", true, "p_{T}");
    Observable obs_eta("Eta", "double", true, "#eta");
    Observable obs_phi("Phi", "double", true, "#phi");
    Observable obs_pid("PID", "int", true, "PID");

    Observable obs_D0_pt("D0_Pt", "double", true, "p_{T}");
    Observable obs_D0_eta("D0_Eta", "double", true, "#eta");
    Observable obs_D0_phi("D0_Phi", "double", false, "#phi");
    Observable obs_D0_rap("D0_Rap", "double", true, "y");
    Observable obs_D0_mpid("D0_MPID", "int", true, "Mother PID");
    Observable obs_numD0s("numD0s", "int", true, "# D0s per event");
    Observable obs_list[10] = { obs_pt, obs_eta, obs_phi, obs_pid, obs_D0_pt, obs_D0_eta, obs_D0_phi, obs_D0_rap, obs_D0_mpid, obs_numD0s };

    for ( Observable obs : obs_list ) {
        cout << "Running observable " << obs.name << endl;
        // if ( obs.numtype == "int" ) {
        obs.hist_pythia_prompt = read_histogram(fin_pythia, obs, "pythia", "prompt", "gen");
        obs.hist_pythia_nonprompt = read_histogram(fin_pythia, obs, "pythia", "nonprompt", "gen");
        obs.hist_herwig_prompt = read_histogram(fin_herwig, obs, "herwig", "prompt", "gen");
        obs.hist_herwig_nonprompt = read_histogram(fin_herwig, obs, "herwig", "nonprompt", "gen");
        // } else {
        //     obs.hist_pythia_prompt = read_histogram(fin_pythia, obs, "pythia", "prompt", "gen");
        //     obs.hist_pythia_nonprompt = read_histogram(fin_pythia, obs, "pythia", "nonprompt", "gen");
        //     obs.hist_herwig_prompt = read_histogram(fin_herwig, obs, "herwig", "prompt", "gen");
        //     obs.hist_herwig_nonprompt = read_histogram(fin_herwig, obs, "herwig", "nonprompt", "gen");
        // }

        cout << "styling hists now" << endl;
        obs.style_hists();
        cout << "drawing hists now" << endl;
        obs.drawPairForObs();
    }
}


void compare_mc_HF_particle_level() {

    // -------- INPUT HISTOGRAMS -------- -- no det level at this point in time
    TString pythia_hists = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/HF_EEC/rootfiles/HF_particle_comparisons/HF_particle_comparisons_pythia.root"; 
    TString herwig_hists = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/HF_EEC/rootfiles/HF_particle_comparisons/HF_particle_comparisons_herwig.root";

    // -------- DEFINE GENERATOR --------
    // Generator gen_anchmc("anchmc", "", anchmc_filepaths, "PWGHF_TreeCreator/tree_Particle_gen", "gen", "GEN ANCHORED MC LHC23a3");
    // Generator gen_pythiafastsim("pythiafastsim", base_filepath_header, pythiafastsim_filepaths, "tree_Particle_gen", "gen", "GEN PYTHIA FASTSIM 1143757");
    // Generator det_anchmc("anchmc", "", anchmc_filepaths, "PWGHF_TreeCreator/tree_Particle", "det", "DET ANCHORED MC LHC23a3");
    // Generator det_pythiafastsim("pythiafastsim", base_filepath_header, pythiafastsim_filepaths, "tree_Particle", "det", "DET PYTHIA FASTSIM 1143757");

    // -------- OPEN FILES --------
    TFile * file_pythia_hists = new TFile(pythia_hists, "READ");
    TFile * file_herwig_hists = new TFile(herwig_hists, "READ");

    // -------- COMPARE GENERATORS --------
    // compareParticleBranches_TChain(outfile, fout_root, gen_anchmc, gen_pythiafastsim);
    // compareParticleBranches_TChain(outfile, fout_root, det_anchmc, det_pythiafastsim);

    analyze(file_pythia_hists, file_herwig_hists);

    file_pythia_hists->Close();
    file_herwig_hists->Close();
}


