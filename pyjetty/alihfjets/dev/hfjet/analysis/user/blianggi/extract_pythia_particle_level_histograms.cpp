// i want to compare generator and detector level HF pythia and herwig
// pp 13 TeV
// this will extract pythia curves into histograms
// two configurations: pythia_prompt, pythia_nonprompt
// particle gen branches: ParticlePt, ParticleEta, ParticlePhi, ParticlePID, MotherPID
// D0 gen branches: ParticlePt, ParticleEta, ParticlePhi, ParticleRapidity, ParticleMCIndex, ParticlePID, MotherPID
// histogram outputs: 
// - particle gen spectra: pT, eta, phi, PID
// - particle det spectra: pT, eta, phi, PID
// - D0 gen spectra: pT, eta, phi, rapidity, motherPID, # D0's per event
// - D0 det spectra: pT, eta, phi, rapidity, motherPID, # D0's per event (but this needs to be done from reconstructed k+pi)

// NOTE: DETECTOR LEVEL IS MORE COMPLICATED -- CAN IMPLEMENT THAT LATER

// std::string base_filepath_header = "/global/cfs/projectdirs/alice/alicepro/hiccup";


class Generator {
public:
    std::string gen_type;
    std::string pathtofiles;

    std::string particle_treename;
    std::string D0_treename;

    std::string gen_or_det;
    std::string label;

    Generator(std::string gen_type_val, std::string pathtofiles_val, std::string particle_treename_val, std::string D0_treename_val, std::string gen_or_det_val, std::string label_val) {
        gen_type = gen_type_val;
        pathtofiles = pathtofiles_val;

        particle_treename = particle_treename_val;
        D0_treename = D0_treename_val;

        gen_or_det = gen_or_det_val;
        label = label_val;
    }
    
};


TChain *makeChain(Generator gen_mc, std::string whichtree, int filecounter_cutoff=-1) {

    cout << "IN MAKE CHAIN!!!" << endl;

    // make TChains
    TChain *chain;
    if (whichtree == "particle") chain = new TChain(gen_mc.particle_treename.c_str());
    else if (whichtree == "D0") chain = new TChain(gen_mc.D0_treename.c_str());

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
        std::string fulltreename = Form("%s", (ntuple_filename).c_str());
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

void fillParticleHistsFromChain( TChain *chain, TH1D *hPt, TH1D *hEta, TH1D *hPhi, TH1I *hPID ) {
    // chain->Draw("ParticlePt >> hPt");
    // chain->Draw("ParticleEta >> hEta");
    // chain->Draw("ParticlePhi >> hPhi");
    
    if (!chain || chain->GetEntries() == 0) {
        std::cerr << "fillParticleHistsFromChain: empty or null chain!" << std::endl;
        return;
    }

    // ---- Branch variables (MATCH TREE TYPES EXACTLY) ----
    double pt, eta, phi;
    int pid;

    // ---- Branch setup ----
    chain->SetBranchStatus("*", 0);

    chain->SetBranchStatus("ParticlePt",  1);
    chain->SetBranchStatus("ParticleEta", 1);
    chain->SetBranchStatus("ParticlePhi", 1);
    chain->SetBranchStatus("ParticlePID", 1);

    chain->SetBranchAddress("ParticlePt",  &pt);
    chain->SetBranchAddress("ParticleEta", &eta);
    chain->SetBranchAddress("ParticlePhi", &phi);
    chain->SetBranchAddress("ParticlePID", &pid);

    // ---- Loop ----
    const Long64_t nEntries = chain->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        chain->GetEntry(i);

        hPt->Fill(pt);
        hEta->Fill(eta);
        hPhi->Fill(phi);
        hPID->Fill(pid);
    }
}

void fillgenD0HistsFromChain( TChain *chain, TH1D *hPt, TH1D *hEta, TH1D *hPhi, TH1D *hRapidity, TH1I *hMotherPID, TH1I *hNumD0s ) {
    
    if (!chain || chain->GetEntries() == 0) {
        std::cerr << "fillD0HistsFromChain: empty or null chain!" << std::endl;
        return;
    }

    // ---- Branch variables (MATCH TREE TYPES EXACTLY) ----
    int evid;
    double pt, eta, phi, rap;
    int mpid;

    // ---- Branch setup ----
    chain->SetBranchStatus("*", 0);

    chain->SetBranchStatus("ev_id",  1);
    chain->SetBranchStatus("ParticlePt",  1);
    chain->SetBranchStatus("ParticleEta", 1);
    chain->SetBranchStatus("ParticlePhi", 1);
    chain->SetBranchStatus("ParticleRapidity", 1);
    chain->SetBranchStatus("MotherPID", 1);

    chain->SetBranchAddress("ev_id",  &evid);
    chain->SetBranchAddress("ParticlePt",  &pt);
    chain->SetBranchAddress("ParticleEta", &eta);
    chain->SetBranchAddress("ParticlePhi", &phi);
    chain->SetBranchAddress("ParticleRapidity", &rap);
    chain->SetBranchAddress("MotherPID", &mpid);

    // ---- Loop ----
    const Long64_t nEntries = chain->GetEntries();
    std::vector<int> ev_id_list;
    for (Long64_t i = 0; i < nEntries; ++i) {
        chain->GetEntry(i);

        hPt->Fill(pt);
        hEta->Fill(eta);
        hPhi->Fill(phi);
        hRapidity->Fill(rap);
        hMotherPID->Fill(mpid);

        ev_id_list.push_back(evid);
    }

    // find how many D0s exist per event
    int current = ev_id_list[0]; 
    int count = 0; 

    for (int x : ev_id_list) {
        if (x == current) {
            count++;
        } else {
            hNumD0s->Fill(count);
            for (int missing = current + 1; missing < x; ++missing) {
                hNumD0s->Fill(0);
            }
            current = x;
            count = 1;
        }
    }
    hNumD0s->Fill(count); // Fill in the last event
}

void filldetD0HistsFromChain( TChain *particlechain, TChain *D0chain, TH1D *hPt, TH1D *hEta, TH1D *hPhi, TH1D *hRapidity, TH1I *hMotherPID ) {
    
    if (!D0chain || D0chain->GetEntries() == 0) {
        std::cerr << "filldetD0HistsFromChain: empty or null D0chain!" << std::endl;
        return;
    }

    // ---- Branch variables (MATCH TREE TYPES EXACTLY) ----
    int evid, D0_evid;
    double pt, D0_pt, D0_eta, D0_phi, D0_rap;
    int pid, mpid, D0_mpid;

    // ---- Branch setup ----
    particlechain->SetBranchStatus("*", 0);
    D0chain->SetBranchStatus("*", 0);

    particlechain->SetBranchStatus("ev_id",  1);
    particlechain->SetBranchStatus("ParticlePt",  1);
    particlechain->SetBranchStatus("ParticlePID", 1);
    particlechain->SetBranchStatus("MotherPID", 1);

    D0chain->SetBranchStatus("ParticlePt",  1);
    D0chain->SetBranchStatus("ParticleEta", 1);
    D0chain->SetBranchStatus("ParticlePhi", 1);
    D0chain->SetBranchStatus("ParticleRapidity", 1);
    D0chain->SetBranchStatus("MotherPID", 1);

    particlechain->SetBranchAddress("ev_id", &evid);
    particlechain->SetBranchAddress("ParticlePt", &pt);
    particlechain->SetBranchAddress("ParticlePID", &pid);
    particlechain->SetBranchAddress("MotherPID", &mpid);

    D0chain->SetBranchAddress("ev_id",  &D0_evid);
    D0chain->SetBranchAddress("ParticlePt",  &D0_pt);
    D0chain->SetBranchAddress("ParticleEta", &D0_eta);
    D0chain->SetBranchAddress("ParticlePhi", &D0_phi);
    D0chain->SetBranchAddress("ParticleRapidity", &D0_rap);
    D0chain->SetBranchAddress("MotherPID", &D0_mpid);

    // ---- Loop ----
    // const Long64_t nEntries = chain->GetEntries();
    // for (Long64_t i = 0; i < nEntries; ++i) {
    //     chain->GetEntry(i);

    //     // find daughter k+pi
    //     // fill det

    //     // hPt->Fill(pt);
    //     // hEta->Fill(eta);
    //     // hPhi->Fill(phi);
    //     // hRapidity->Fill(rap);
    //     // hMotherPID->Fill(mpid);
    // }
}

void compareParticleBranches_TChain(std::ofstream &outfile, TFile * fout_root, Generator gen1, Generator gen2) {

    // -------- CHAINS --------
    TChain *chain1_particle = makeChain(gen1, "particle"); //pythia prompt
    TChain *chain2_particle = makeChain(gen2, "particle"); //pythia nonprompt
    TChain *chain1_D0 = makeChain(gen1, "D0"); //pythia prompt
    TChain *chain2_D0 = makeChain(gen2, "D0"); //pythia nonprompt

    std::cout << "particle chain1 entries: " << chain1_particle->GetEntries() << std::endl;
    std::cout << "particle chain2 entries: " << chain2_particle->GetEntries() << std::endl;
    std::cout << "D0 chain1 entries: " << chain1_D0->GetEntries() << std::endl;
    std::cout << "D0 chain2 entries: " << chain2_D0->GetEntries() << std::endl;

    // save number of entries to a file
    outfile << "Number of entries in particle " << gen1.label << ": " << chain1_particle->GetEntries() << std::endl;
    outfile << "Number of entries in particle " << gen2.label << ": " << chain2_particle->GetEntries() << std::endl;
    outfile << "Number of entries in D0 " << gen1.label << ": " << chain1_D0->GetEntries() << std::endl;
    outfile << "Number of entries in D0 " << gen2.label << ": " << chain2_D0->GetEntries() << std::endl;

    // -------- HISTOGRAMS --------
    TH1D *hPt_1  = new TH1D(Form("hPt_%s_%s", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()),  Form("Particle p_{T} %s;p_{T};Entries", gen1.gen_or_det.c_str()), 200, 0, 200);
    TH1D *hPt_2  = new TH1D(Form("hPt_%s_%s", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()),  Form("Particle p_{T} %s;p_{T};Entries", gen2.gen_or_det.c_str()), 200, 0, 200);

    TH1D *hEta_1 = new TH1D(Form("hEta_%s_%s", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()), Form("Particle #eta %s;#eta;Entries", gen1.gen_or_det.c_str()), 100, -5, 5);
    TH1D *hEta_2 = new TH1D(Form("hEta_%s_%s", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()), Form("Particle #eta %s;#eta;Entries", gen2.gen_or_det.c_str()), 100, -5, 5);

    TH1D *hPhi_1 = new TH1D(Form("hPhi_%s_%s", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()), Form("Particle #phi %s;#phi;Entries", gen1.gen_or_det.c_str()), 64, -3.3, 6.5); //-TMath::Pi(), 2*TMath::Pi());
    TH1D *hPhi_2 = new TH1D(Form("hPhi_%s_%s", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()), Form("Particle #phi %s;#phi;Entries", gen2.gen_or_det.c_str()), 64, -3.3, 6.5); //-TMath::Pi(), 2*TMath::Pi());

    TH1I *hPID_1 = new TH1I(Form("hPID_%s_%s", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()), Form("Particle PID %s;PID;Entries", gen1.gen_or_det.c_str()), 5000, -2500, 2500); //-TMath::Pi(), 2*TMath::Pi());
    TH1I *hPID_2 = new TH1I(Form("hPID_%s_%s", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()), Form("Particle PID %s;PID;Entries", gen2.gen_or_det.c_str()), 5000, -2500, 2500); //-TMath::Pi(), 2*TMath::Pi());

    // D0 hists
    TH1D *hD0_Pt_1  = new TH1D(Form("hD0_Pt_%s_%s", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()),  Form("D0 p_{T} %s;p_{T};Entries", gen1.gen_or_det.c_str()), 200, 0, 200);
    TH1D *hD0_Pt_2  = new TH1D(Form("hD0_Pt_%s_%s", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()),  Form("D0 p_{T} %s;p_{T};Entries", gen2.gen_or_det.c_str()), 200, 0, 200);

    TH1D *hD0_Eta_1 = new TH1D(Form("hD0_Eta_%s_%s", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()), Form("D0 #eta %s;#eta;Entries", gen1.gen_or_det.c_str()), 100, -5, 5);
    TH1D *hD0_Eta_2 = new TH1D(Form("hD0_Eta_%s_%s", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()), Form("D0 #eta %s;#eta;Entries", gen2.gen_or_det.c_str()), 100, -5, 5);

    TH1D *hD0_Phi_1 = new TH1D(Form("hD0_Phi_%s_%s", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()), Form("D0 #phi %s;#phi;Entries", gen1.gen_or_det.c_str()), 64, -3.3, 6.5); //-TMath::Pi(), 2*TMath::Pi());
    TH1D *hD0_Phi_2 = new TH1D(Form("hD0_Phi_%s_%s", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()), Form("D0 #phi %s;#phi;Entries", gen2.gen_or_det.c_str()), 64, -3.3, 6.5); //-TMath::Pi(), 2*TMath::Pi());

    TH1D *hD0_Rap_1 = new TH1D(Form("hD0_Rap_%s_%s", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()), Form("D0 PID %s;PID;Entries", gen1.gen_or_det.c_str()), 100, -5, 5); //-TMath::Pi(), 2*TMath::Pi());
    TH1D *hD0_Rap_2 = new TH1D(Form("hD0_Rap_%s_%s", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()), Form("D0 PID %s;PID;Entries", gen2.gen_or_det.c_str()), 100, -5, 5); //-TMath::Pi(), 2*TMath::Pi());

    TH1I *hD0_MPID_1 = new TH1I(Form("hD0_MPID_%s_%s", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()), Form("D0 Mother PID %s;MPID;Entries", gen1.gen_or_det.c_str()), 5000, -2500, 2500); //-TMath::Pi(), 2*TMath::Pi());
    TH1I *hD0_MPID_2 = new TH1I(Form("hD0_MPID_%s_%s", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()), Form("D0 Mother PID %s;MPID;Entries", gen2.gen_or_det.c_str()), 5000, -2500, 2500); //-TMath::Pi(), 2*TMath::Pi());

    TH1I *hnumD0s_1 = new TH1I(Form("hnumD0s_%s_%s", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()), Form("D0 PID %s;PID;Entries", gen1.gen_or_det.c_str()), 10, 0, 10); //-TMath::Pi(), 2*TMath::Pi());
    TH1I *hnumD0s_2 = new TH1I(Form("hnumD0s_%s_%s", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()), Form("D0 PID %s;PID;Entries", gen2.gen_or_det.c_str()), 10, 0, 10); //-TMath::Pi(), 2*TMath::Pi());


    // -------- FILL --------
    cout << "filling first file particle hists " << endl;
    fillParticleHistsFromChain(chain1_particle, hPt_1, hEta_1, hPhi_1, hPID_1);
    cout << "filling second file particle hists " << endl;
    fillParticleHistsFromChain(chain2_particle, hPt_2, hEta_2, hPhi_2, hPID_2);
    cout << "filling first file D0 hists " << endl;
    if ( gen1.gen_or_det == "gen" ) fillgenD0HistsFromChain(chain1_D0, hD0_Pt_1, hD0_Eta_1, hD0_Phi_1, hD0_Rap_1, hD0_MPID_1, hnumD0s_1);
    // else if ( gen1.gen_or_det == "det" ) fillDetD0HistsFromChain();
    cout << "filling second file D0 hists " << endl;
    if ( gen2.gen_or_det == "gen" ) fillgenD0HistsFromChain(chain2_D0, hD0_Pt_1, hD0_Eta_1, hD0_Phi_1, hD0_Rap_1, hD0_MPID_1, hnumD0s_2);
    // else if ( gen2.gen_or_det == "det" ) fillDetD0HistsFromChain();

    // -------- STYLE --------
    hPt_1->SetLineColor(kRed);
    hPt_2->SetLineColor(kBlue);

    hEta_1->SetLineColor(kRed);
    hEta_2->SetLineColor(kBlue);

    hPhi_1->SetLineColor(kRed);
    hPhi_2->SetLineColor(kBlue);

    hPID_1->SetLineColor(kRed);
    hPID_2->SetLineColor(kBlue);

    // D0 hists
    hD0_Pt_1 ->SetLineColor(kRed);
    hD0_Pt_2 ->SetLineColor(kBlue);

    hD0_Eta_1->SetLineColor(kRed);
    hD0_Eta_2->SetLineColor(kBlue);

    hD0_Phi_1->SetLineColor(kRed);
    hD0_Phi_2->SetLineColor(kBlue);

    hD0_Rap_1->SetLineColor(kRed);
    hD0_Rap_2->SetLineColor(kBlue);

    hD0_MPID_1->SetLineColor(kRed);
    hD0_MPID_2->SetLineColor(kBlue);

    hnumD0s_1->SetLineColor(kRed);
    hnumD0s_2->SetLineColor(kBlue);

    // -------- DRAW --------
    auto savePair = [](TFile * fout_root, Generator gen1, Generator gen2, TH1D *h1, TH1D *h2, std::string name) {
        
        
        // leg->AddEntry(h1, gen1.label.c_str(), "l"); //"ANCH MC LHC23a3", "l");
        // leg->AddEntry(h2, gen2.label.c_str(), "l"); //"PYTHIA FASTSIM 1143757", "l");

        
        TH1D *h_ratio = (TH1D *)h2->Clone(Form("h_ratio_%s", name.c_str()));
        h_ratio->Divide(h1);

        h_ratio->SetTitle(Form("Ratio %s", name.c_str()));
        h_ratio->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
        h_ratio->GetYaxis()->SetTitle("PYTHIA / ANCH MC");

        // Save to root file
        fout_root->cd();
        h1->Write();
        h2->Write();
        h_ratio->Write();
    };

    auto savePairI = [](TFile * fout_root, Generator gen1, Generator gen2, TH1I *h1, TH1I *h2, std::string name) {
        
        
        // leg->AddEntry(h1, gen1.label.c_str(), "l"); //"ANCH MC LHC23a3", "l");
        // leg->AddEntry(h2, gen2.label.c_str(), "l"); //"PYTHIA FASTSIM 1143757", "l");

        
        TH1D *h_ratio = (TH1D *)h2->Clone(Form("h_ratio_%s", name.c_str()));
        h_ratio->Divide(h1);

        h_ratio->SetTitle(Form("Ratio %s", name.c_str()));
        h_ratio->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
        h_ratio->GetYaxis()->SetTitle("PYTHIA / ANCH MC");

        // Save to root file
        fout_root->cd();
        h1->Write();
        h2->Write();
        h_ratio->Write();
    };

    savePair(fout_root, gen1, gen2, hPt_1,  hPt_2, "Pt" + gen1.gen_or_det);
    savePair(fout_root, gen1, gen2, hEta_1, hEta_2, "Eta" + gen1.gen_or_det);
    savePair(fout_root, gen1, gen2, hPhi_1, hPhi_2, "Phi" + gen1.gen_or_det);
    savePairI(fout_root, gen1, gen2, hPID_1, hPID_2, "PID" + gen1.gen_or_det);

    savePair(fout_root, gen1, gen2, hD0_Pt_1,  hD0_Pt_2, "D0_Pt" + gen1.gen_or_det);
    savePair(fout_root, gen1, gen2, hD0_Eta_1, hD0_Eta_2, "D0_Eta" + gen1.gen_or_det);
    savePair(fout_root, gen1, gen2, hD0_Phi_1, hD0_Phi_2, "D0_Phi" + gen1.gen_or_det);
    savePair(fout_root, gen1, gen2, hD0_Rap_1, hD0_Rap_2, "D0_Rap" + gen1.gen_or_det);
    savePairI(fout_root, gen1, gen2, hD0_MPID_1, hD0_MPID_2, "D0_MPID" + gen1.gen_or_det);
    savePairI(fout_root, gen1, gen2, hnumD0s_1, hnumD0s_2, "numD0s" + gen1.gen_or_det);
}




void extract_pythia_particle_level_histograms() {

    // -------- INPUT DIRECTORIES --------
    // post eff smearing -- generator + detector level
    std::string pythia_prompt_filepaths = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/pythiagen/tree_fastsim/45178629/45154942/files.txt"; 
    std::string pythia_nonprompt_filepaths = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/pythiagen/tree_fastsim/46306341/46293548/files.txt"; 

    // -------- DEFINE GENERATOR --------
    Generator gen_pythia_prompt("pythia_prompt", pythia_prompt_filepaths, "tree_Particle_gen", "tree_D0_gen", "gen", "Pythia prompt, gen");
    Generator gen_pythia_nonprompt("pythia_nonprompt", pythia_nonprompt_filepaths, "tree_Particle_gen", "tree_D0_gen", "gen", "Pythia non-prompt, gen");
    // Generator det_pythia_prompt("pythia_prompt", pythia_prompt_filepaths, "tree_Particle", "tree_D0", "det", "Pythia prompt, det");
    // Generator det_pythia_nonprompt("pythia_nonprompt", pythia_nonprompt_filepaths, "tree_Particle", "tree_D0", "det", "Pythia non-prompt,  det");
    
    // -------- OPEN OUTPUT FILEs --------
    std::ofstream outfile("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/HF_EEC/plots/HF_particle_comparisons/number_of_entries.txt");
    TFile * fout_root = new TFile("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/HF_EEC/rootfiles/HF_particle_comparisons/HF_particle_comparisons.root", "RECREATE");

    // -------- COMPARE GENERATORS --------
    compareParticleBranches_TChain(outfile, fout_root, gen_pythia_prompt, gen_pythia_nonprompt);
    // compareParticleBranches_TChain(outfile, fout_root, det_anchmc, det_pythiafastsim);

    // CLOSE OUTPUT FILE
    outfile.close();
    fout_root->Close();
}


