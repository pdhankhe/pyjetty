// i want to compare generator and detector level HF pythia and herwig
// pp 13 TeV
// there are four configurations: pythia_prompt, pythia_nonprompt, herwig_prompt, herwig_nonprompt
// are there differences in herwig and pythia, even at particle level?
// branches: ParticlePt, ParticleEta, ParticlePhi, ParticlePID

std::string base_filepath_header = "/global/cfs/projectdirs/alice/alicepro/hiccup";


class Generator {
public:
    std::string gen_type;
    std::string basepath;
    std::string pathtofiles;
    std::string treename;

    std::string gen_or_det;
    std::string label;

    Generator(std::string gen_type_val, std::string basepath_val, std::string pathtofiles_val, std::string treename_val, std::string gen_or_det_val, std::string label_val) {
        gen_type = gen_type_val;
        basepath = basepath_val;
        pathtofiles = pathtofiles_val;
        treename = treename_val;

        gen_or_det = gen_or_det_val;
        label = label_val;
    }
    
};


TChain *makeChain(Generator gen_mc, int filecounter_cutoff=-1) {

    cout << "IN MAKE CHAIN!!!" << endl;
    // TChain *chain = new TChain("PWGHF_TreeCreator/tree_Particle_gen");

    // TSystemDirectory dir("inputDir", dirPath.c_str());
    // TList *files = dir.GetListOfFiles();
    // if (!files) return chain;

    // for (auto obj : *files) {
    //     TSystemFile *file = (TSystemFile *)obj;
    //     std::string fname = file->GetName();

    //     if (!file->IsDirectory() && fname.find(".root") != std::string::npos) {
    //         std::string fullPath = dirPath + "/" + fname;
    //         chain->Add(fullPath.c_str());
    //     }
    // }
    // return chain;

    // make TChains
    TChain *chain = new TChain(gen_mc.treename.c_str());

    std::ifstream filelist(gen_mc.pathtofiles.c_str());
    if (!filelist.is_open()) {
        std::cerr << "Error: Could not open " << gen_mc.pathtofiles << std::endl;
        return chain;
    }

    std::string ntuple_filename;
    int filecounter = 0;
    
    // Loop through each line in filelist
    while (std::getline(filelist, ntuple_filename)) {

        if (filecounter == filecounter_cutoff) break;
        // cout << "file " << filecounter << endl;
        // cout << "basepath " << gen_mc.basepath << endl;
        // cout << "ntuple " << ntuple_filename << endl;
        

        // cout << Form("%s/%sPWGHF_TreeCreator/tree_Particle_gen", base_filepath.c_str(), (ntuple_filename + "/").c_str()) << endl;
        std::string fulltreename = Form("%s%s", gen_mc.basepath.c_str(), (ntuple_filename).c_str());
        chain->Add(fulltreename.c_str());
        

        if (filecounter%100 == 0) {
            cout << "num tree entries " << chain->GetEntries() << endl;
        }
        

        filecounter++;
    }

    // Close the filelist.txt file
    filelist.close();

    return chain;
}

void fillHistsFromChain_float( TChain *chain, TH1D *hPt, TH1D *hEta, TH1D *hPhi ) {
    if (!chain || chain->GetEntries() == 0) {
        std::cerr << "fillHistsFromChain: empty or null chain!" << std::endl;
        return;
    }

    // ---- Branch variables (MATCH TREE TYPES EXACTLY) ----
    float pt, eta, phi;

    // ---- Branch setup ----
    chain->SetBranchStatus("*", 0);

    chain->SetBranchStatus("ParticlePt",  1);
    chain->SetBranchStatus("ParticleEta", 1);
    chain->SetBranchStatus("ParticlePhi", 1);

    chain->SetBranchAddress("ParticlePt",  &pt);
    chain->SetBranchAddress("ParticleEta", &eta);
    chain->SetBranchAddress("ParticlePhi", &phi);

    // ---- Loop ----
    const Long64_t nEntries = chain->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        chain->GetEntry(i);

        hPt->Fill(pt);
        hEta->Fill(eta);
        hPhi->Fill(phi);
    }
}

void fillHistsFromChain_double( TChain *chain, TH1D *hPt, TH1D *hEta, TH1D *hPhi ) {
    // chain->Draw("ParticlePt >> hPt");
    // chain->Draw("ParticleEta >> hEta");
    // chain->Draw("ParticlePhi >> hPhi");
    
    if (!chain || chain->GetEntries() == 0) {
        std::cerr << "fillHistsFromChain: empty or null chain!" << std::endl;
        return;
    }

    // ---- Branch variables (MATCH TREE TYPES EXACTLY) ----
    double pt, eta, phi;

    // ---- Branch setup ----
    chain->SetBranchStatus("*", 0);

    chain->SetBranchStatus("ParticlePt",  1);
    chain->SetBranchStatus("ParticleEta", 1);
    chain->SetBranchStatus("ParticlePhi", 1);

    chain->SetBranchAddress("ParticlePt",  &pt);
    chain->SetBranchAddress("ParticleEta", &eta);
    chain->SetBranchAddress("ParticlePhi", &phi);

    // ---- Loop ----
    const Long64_t nEntries = chain->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        chain->GetEntry(i);

        hPt->Fill(pt);
        hEta->Fill(eta);
        hPhi->Fill(phi);
    }
}

void compareParticleBranches_TChain(std::ofstream &outfile, TFile * fout_root, Generator gen1, Generator gen2) {

    // -------- CHAINS --------
    TChain *chain1 = makeChain(gen1); //gen_anchmc);
    TChain *chain2 = makeChain(gen2); //gen_pythiafastsim);

    std::cout << "Chain1 entries: " << chain1->GetEntries() << std::endl;
    std::cout << "Chain2 entries: " << chain2->GetEntries() << std::endl;

    // save number of entries to a file
    outfile << "Number of entries in " << gen1.label << ": " << chain1->GetEntries() << std::endl;
    outfile << "Number of entries in " << gen2.label << ": " << chain2->GetEntries() << std::endl;

    // -------- HISTOGRAMS --------
    TH1D *hPt_1  = new TH1D(Form("hPt_%s_%s", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()),  Form("Particle p_{T} %s;p_{T};Entries", gen1.gen_or_det.c_str()), 200, 0, 200);
    TH1D *hPt_2  = new TH1D(Form("hPt_%s_%s", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()),  Form("Particle p_{T} %s;p_{T};Entries", gen2.gen_or_det.c_str()), 200, 0, 200);

    TH1D *hEta_1 = new TH1D(Form("hEta_%s_%s", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()), Form("Particle #eta %s;#eta;Entries", gen1.gen_or_det.c_str()), 100, -5, 5);
    TH1D *hEta_2 = new TH1D(Form("hEta_%s_%s", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()), Form("Particle #eta %s;#eta;Entries", gen2.gen_or_det.c_str()), 100, -5, 5);

    TH1D *hPhi_1 = new TH1D(Form("hPhi_%s_%s", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()), Form("Particle #phi %s;#phi;Entries", gen1.gen_or_det.c_str()), 64, -3.3, 6.5); //-TMath::Pi(), 2*TMath::Pi());
    TH1D *hPhi_2 = new TH1D(Form("hPhi_%s_%s", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()), Form("Particle #phi %s;#phi;Entries", gen2.gen_or_det.c_str()), 64, -3.3, 6.5); //-TMath::Pi(), 2*TMath::Pi());


    // -------- FILL --------
    cout << "filling first file hists " << endl;
    fillHistsFromChain_float(chain1, hPt_1, hEta_1, hPhi_1);
    cout << "filling second file hists " << endl;
    fillHistsFromChain_double(chain2, hPt_2, hEta_2, hPhi_2);

    // -------- STYLE --------
    hPt_1->SetLineColor(kRed);
    hPt_2->SetLineColor(kBlue);

    hEta_1->SetLineColor(kRed);
    hEta_2->SetLineColor(kBlue);

    hPhi_1->SetLineColor(kRed);
    hPhi_2->SetLineColor(kBlue);

    // -------- DRAW --------
    auto drawPair = [](TFile * fout_root, Generator gen1, Generator gen2, TH1D *h1, TH1D *h2, std::string name) {
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
        if ( name.find("Pt") != std::string::npos ) gPad->SetLogy();

        // determine max
        double max_height = max(h1->GetMaximum(), h2->GetMaximum());
        h1->SetMaximum(max_height*1.2);

        h1->Draw("hist");
        h2->Draw("hist same");

        TLegend *leg = new TLegend(0.15, 0.75, 0.4, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(h1, gen1.label.c_str(), "l"); //"ANCH MC LHC23a3", "l");
        leg->AddEntry(h2, gen2.label.c_str(), "l"); //"PYTHIA FASTSIM 1143757", "l");
        leg->Draw();

        // ======================
        // Bottom pad: ratio
        // ======================
        pad2->cd();

        TH1D *h_ratio = (TH1D *)h2->Clone(Form("h_ratio_%s", name.c_str()));
        h_ratio->Divide(h1);

        h_ratio->SetTitle(Form("Ratio %s", name.c_str()));
        h_ratio->GetYaxis()->SetTitle("PYTHIA / ANCH MC");
        h_ratio->GetYaxis()->SetNdivisions(505);
        h_ratio->GetYaxis()->SetTitleSize(0.10);
        h_ratio->GetYaxis()->SetLabelSize(0.08);
        h_ratio->GetYaxis()->SetTitleOffset(0.5);

        h_ratio->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
        h_ratio->GetXaxis()->SetTitleSize(0.12);
        h_ratio->GetXaxis()->SetLabelSize(0.10);

        // h_ratio->SetMinimum(0.5);
        // h_ratio->SetMaximum(1.5);

        h_ratio->Draw("hist");

        c->SaveAs(Form("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/plots/compare_pythia_anchmc/%s_comparison.pdf", name.c_str()));

        // Save to root file
        fout_root->cd();
        h1->Write();
        h2->Write();
        h_ratio->Write();
    };

    drawPair(fout_root, gen1, gen2, hPt_1,  hPt_2, "Pt" + gen1.gen_or_det);
    drawPair(fout_root, gen1, gen2, hEta_1, hEta_2, "Eta" + gen1.gen_or_det);
    drawPair(fout_root, gen1, gen2, hPhi_1, hPhi_2, "Phi" + gen1.gen_or_det);
}




void compare_mc_particle_level() {

    // -------- INPUT DIRECTORIES --------
    // post eff smearing -- generator + detector level
    std::string anchmc_filepaths = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/data/LHC23a3/806/files.txt"; // << doesnt need header
    std::string pythiafastsim_filepaths = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/generators/pythia_alice/tree_fastsim/1143757/files.txt"; // << needs /global/cfs/projectdirs/alice/alicepro/hiccup 

    // -------- DEFINE GENERATOR --------
    Generator gen_anchmc("anchmc", "", anchmc_filepaths, "PWGHF_TreeCreator/tree_Particle_gen", "gen", "GEN ANCHORED MC LHC23a3");
    Generator gen_pythiafastsim("pythiafastsim", base_filepath_header, pythiafastsim_filepaths, "tree_Particle_gen", "gen", "GEN PYTHIA FASTSIM 1143757");
    Generator det_anchmc("anchmc", "", anchmc_filepaths, "PWGHF_TreeCreator/tree_Particle", "det", "DET ANCHORED MC LHC23a3");
    Generator det_pythiafastsim("pythiafastsim", base_filepath_header, pythiafastsim_filepaths, "tree_Particle", "det", "DET PYTHIA FASTSIM 1143757");

    // -------- OPEN OUTPUT FILEs --------
    std::ofstream outfile("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/plots/compare_pythia_anchmc/number_of_entries.txt");
    TFile * fout_root = new TFile("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/rootfiles/compare_pythia_anchmc/compare_pythia_anchmc.root", "RECREATE");

    // -------- COMPARE GENERATORS --------
    compareParticleBranches_TChain(outfile, fout_root, gen_anchmc, gen_pythiafastsim);
    compareParticleBranches_TChain(outfile, fout_root, det_anchmc, det_pythiafastsim);

    // CLOSE OUTPUT FILE
    outfile.close();
    fout_root->Close();
}


