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

// RUN LIKE SO: root -l 'extract_pythia_particle_level_histograms.cpp("h")' for herwig or ("p") for pythia
// ALERT! I KNOW THIS IS REALLY ANNOYING but if "x" selected with "p", then need to change branch types in fillParticleHistsFromChain(), fillgenD0HistsFromChain() to be all floats!

// std::string base_filepath_header = "/global/cfs/projectdirs/alice/alicepro/hiccup";

#include <iostream>
#include <ctime>

class Generator {
public:
    std::string gen_type;
    std::string pathtofiles;

    std::string particle_treename;
    std::string D0_treename;

    std::string gen_or_det;
    std::string label;

    std::string sf_filepath; //scale factor file
    std::string sf_individuals_filepath; // scale factor list for every individual file

    Generator(std::string gen_type_val, std::string pathtofiles_val, std::string particle_treename_val, std::string D0_treename_val, std::string gen_or_det_val, std::string label_val, std::string sf_filepath_val, std::string sf_individuals_filepath_val) {
        gen_type = gen_type_val;
        pathtofiles = pathtofiles_val;

        particle_treename = particle_treename_val;
        D0_treename = D0_treename_val;

        gen_or_det = gen_or_det_val;
        label = label_val;

        sf_filepath = sf_filepath_val;
        sf_individuals_filepath = sf_individuals_filepath_val;
    }

    // this is assuming only 10 bins (for HF)
    void scale_hists(std::vector<TH1D *>& vec_hists) {

        std::ifstream sf_file(sf_filepath);
        int id;
        char colon;
        double scale;
        for ( int i = 0; i < 10; i++ ) {

            // The stream reads: [Integer] -> [Char] -> [Double]
            sf_file >> id >> colon >> scale;
            if (sf_file.fail()) break;
            
            cout << i << ": scaling by! " << scale << endl;
            vec_hists[i]->Scale(scale, "width");
        }

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

// Returns the bin number from the path. 
// Assuming it is always 3rd slash from the end (i.e. <run_number>/5/40/AnalysisResults.root)
int extractBin(std::string path) {
    size_t p1 = path.rfind('/');
    size_t p2 = path.rfind('/', p1 - 1);
    size_t p3 = path.rfind('/', p2 - 1);
    int bin_num = std::stoi(path.substr(p3 + 1, p2 - p3 - 1));
    // cout << "file! " << path << endl;
    // cout << "bin num! " << bin_num << endl;
    return bin_num;
}

std::vector<TChain *> makeChain_WithCS(Generator gen_mc, std::string whichtree, int filecounter_cutoff=-1) {

    cout << "IN MAKE CHAIN WITH CS!!!" << endl;

    // make TChains
    const int NCHAINS = 10;
    std::vector<TChain *> chains;
    chains.reserve(NCHAINS);
    for (int i = 0; i < NCHAINS; ++i) {
        if (whichtree == "particle") chains.push_back( new TChain(gen_mc.particle_treename.c_str()) );
        else if (whichtree == "D0") chains.push_back( new TChain(gen_mc.D0_treename.c_str()) );
    }


    std::string ntuple_filename;
    int filecounter = 0;

    // Restart file reading
    std::ifstream filelist(gen_mc.pathtofiles.c_str());
    filelist.clear();                 // clear EOF + error flags
    filelist.seekg(0, std::ios::beg); // go back to start of file
    
    // Loop through each line in filelist
    while (std::getline(filelist, ntuple_filename)) {

        if (filecounter == filecounter_cutoff) break;

        int bin_num = extractBin(ntuple_filename);
        std::string fulltreename = Form("%s", (ntuple_filename).c_str());
        chains[bin_num-1]->Add(fulltreename.c_str());

        if (filecounter%100 == 0) {
            cout << "num tree entries " << chains[0]->GetEntries();
            for (int i=1; i<10; i++) cout << " " << chains[i]->GetEntries();
            cout << endl;
        }
        filecounter++;
    }

    // Close the filelist.txt file
    filelist.close();

    return chains;
}

TH1D * addHists(std::vector<TH1D*> histVector, std::string histname) {

    TH1D* hcomb = (TH1D*)histVector[0]->Clone(histname.c_str());
    for (int i = 1; i < 10; i++) {
        hcomb->Add(histVector[i]);
    }
    return hcomb;
}

void fillParticleHistsFromChain( TChain *chain, TH1D *hPt, TH1D *hEta, TH1D *hPhi, TH1I *hPID, bool fillOnlyPt = false, bool include_neutrals = true ) {
    
    if (!chain || chain->GetEntries() == 0) {
        std::cerr << "fillParticleHistsFromChain: empty or null chain!" << std::endl;
        return;
    }

    chain->ResetBranchAddresses();
    // ---- Branch variables (MATCH TREE TYPES EXACTLY) ----
    double pt, eta, phi;
    Long64_t pid;
    // float pt, eta, phi, pid;

    // ---- Branch setup ----
    chain->SetBranchStatus("*", 0);

    chain->SetBranchStatus("ParticlePt",  1);
    if (!fillOnlyPt) {
        chain->SetBranchStatus("ParticleEta", 1);
        chain->SetBranchStatus("ParticlePhi", 1);
        chain->SetBranchStatus("ParticlePID", 1);
    }

    chain->SetBranchAddress("ParticlePt",  &pt);
    if (!fillOnlyPt) {
        chain->SetBranchAddress("ParticleEta", &eta);
        chain->SetBranchAddress("ParticlePhi", &phi);
        chain->SetBranchAddress("ParticlePID", &pid);
    }

    // ---- Loop ----
    const Long64_t nEntries = chain->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        chain->GetEntry(i);

        if (!include_neutrals) {
            // if (i < 50 ) cout << "checking pid " << pid << " pt: " << pt << "eta" << eta << "phi" << phi << endl;
            if (abs(pid)==22 or abs(pid)==12 or abs(pid)==14 or abs(pid)==16 or abs(pid)==130 or abs(pid)==2112) {
                // cout << pid << " neutral here!" << endl;
                continue;
            }
        }

        hPt->Fill(pt);
        if (!fillOnlyPt) {
            hEta->Fill(eta);
            hPhi->Fill(phi);
            hPID->Fill(pid);
        }
    }
}

void fillgenD0HistsFromChain( TChain *chain, TH1D *hPt, TH1D *hEta, TH1D *hPhi, TH1D *hRapidity, TH1I *hMotherPID, TH1I *hNumD0s, bool fillOnlyPt = false ) {
    
    if (!chain || chain->GetEntries() == 0) {
        std::cerr << "fillD0HistsFromChain: empty or null chain!" << std::endl;
        return;
    }

    chain->ResetBranchAddresses();

    // ---- Branch variables (MATCH TREE TYPES EXACTLY) - these are all D0 branches ----
    Long64_t evid;
    double pt, eta, phi, rap, mpid;
    // float evid, pt, eta, phi, rap, mpid;

    // ---- Branch setup ----
    chain->SetBranchStatus("*", 0);

    chain->SetBranchStatus("ParticlePt",  1);
    chain->SetBranchStatus("ParticleRapidity", 1);
    if (!fillOnlyPt) {
        chain->SetBranchStatus("ev_id",  1);
        chain->SetBranchStatus("ParticleEta", 1);
        chain->SetBranchStatus("ParticlePhi", 1);
        chain->SetBranchStatus("MotherPID", 1);
    }

    
    chain->SetBranchAddress("ParticlePt",  &pt);
    chain->SetBranchAddress("ParticleRapidity", &rap);
    if (!fillOnlyPt) {
        chain->SetBranchAddress("ev_id",  &evid);
        chain->SetBranchAddress("ParticleEta", &eta);
        chain->SetBranchAddress("ParticlePhi", &phi);
        chain->SetBranchAddress("MotherPID", &mpid);
    }

    // ---- Loop ----
    const Long64_t nEntries = chain->GetEntries();
    std::vector<int> ev_id_list;
    for (Long64_t i = 0; i < nEntries; ++i) {
        chain->GetEntry(i);

        hPt->Fill(pt);
        hRapidity->Fill(rap);
        if (!fillOnlyPt) {
            hEta->Fill(eta);
            hPhi->Fill(phi);
            hMotherPID->Fill(mpid);

            ev_id_list.push_back(evid);
        }
    }

    // find how many D0s exist per event
    if (!fillOnlyPt) {
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
}

/*
void filldetD0HistsFromChain( TChain *particlechain, TChain *D0chain, TH1D *hPt, TH1D *hEta, TH1D *hPhi, TH1D *hRapidity, TH1I *hMotherPID ) {
    
    if (!D0chain || D0chain->GetEntries() == 0) {
        std::cerr << "filldetD0HistsFromChain: empty or null D0chain!" << std::endl;
        return;
    }
    chain->ResetBranchAddresses();

    // ---- Branch variables (MATCH TREE TYPES EXACTLY) ----
    Long64_t evid, D0_evid;
    double pt, D0_pt, D0_eta, D0_phi, D0_rap, D0_mpid; // D0 MPID was saved as double for some reason
    Long64_t pid, mpid;

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
*/

void fillHistsPerFile(Generator gen_mc, std::string whichtree, TH1D * hPt, TH1D * hRapidity, bool fillOnlyPt = false) {

    std::ifstream fileList(gen_mc.pathtofiles.c_str());
    std::ifstream scaleList(gen_mc.sf_individuals_filepath.c_str());
    
    std::string fileName;
    double weight;

    while (fileList >> fileName && scaleList >> weight) {
        TFile* f = TFile::Open(fileName.c_str());
        if (!f || f->IsZombie()) continue;

        // Get the tree - change "events" to your actual tree name
        TTree* tree;
        if (whichtree == "particle") tree = (TTree*)f->Get(gen_mc.particle_treename.c_str());
        else if (whichtree == "D0") tree = (TTree*)f->Get(gen_mc.D0_treename.c_str());
        if (!tree) {
            std::cout << "Tree not found in " << fileName << std::endl;
            f->Close();
            continue;
        }

        // Setup the branch address
        // Using double or float depending on how you saved your tree
        tree->ResetBranchAddresses();
        
        double pt; 
        double rap;
        tree->SetBranchAddress("ParticlePt", &pt);
        if (!fillOnlyPt) tree->SetBranchAddress("ParticleRapidity", &rap);

        // 4. Loop over entries in this specific pT-hat bin
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; i++) {
            tree->GetEntry(i);
            
            // Fill with the weight: (sigma / nAccepted)
            hPt->Fill(pt, weight);
            if (!fillOnlyPt) hRapidity->Fill(rap, weight);
        }

        // std::cout << "Finished " << fileName << " (" << nEntries << " entries) scaled by " << weight << std::endl;
        f->Close();
    }
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
    TH1D *hD0_Pt_1  = new TH1D(Form("hD0_Pt_%s_%s", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()),  Form("D0 p_{T} %s;p_{T};Entries", gen1.gen_or_det.c_str()), 100, 0, 100);
    TH1D *hD0_Pt_2  = new TH1D(Form("hD0_Pt_%s_%s", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()),  Form("D0 p_{T} %s;p_{T};Entries", gen2.gen_or_det.c_str()), 100, 0, 100);

    TH1D *hD0_Eta_1 = new TH1D(Form("hD0_Eta_%s_%s", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()), Form("D0 #eta %s;#eta;Entries", gen1.gen_or_det.c_str()), 100, -5, 5);
    TH1D *hD0_Eta_2 = new TH1D(Form("hD0_Eta_%s_%s", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()), Form("D0 #eta %s;#eta;Entries", gen2.gen_or_det.c_str()), 100, -5, 5);

    TH1D *hD0_Phi_1 = new TH1D(Form("hD0_Phi_%s_%s", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()), Form("D0 #phi %s;#phi;Entries", gen1.gen_or_det.c_str()), 64, -3.3, 6.5); //-TMath::Pi(), 2*TMath::Pi());
    TH1D *hD0_Phi_2 = new TH1D(Form("hD0_Phi_%s_%s", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()), Form("D0 #phi %s;#phi;Entries", gen2.gen_or_det.c_str()), 64, -3.3, 6.5); //-TMath::Pi(), 2*TMath::Pi());

    TH1D *hD0_Rap_1 = new TH1D(Form("hD0_Rap_%s_%s", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()), Form("D0 y %s;y;Entries", gen1.gen_or_det.c_str()), 100, -5, 5); //-TMath::Pi(), 2*TMath::Pi());
    TH1D *hD0_Rap_2 = new TH1D(Form("hD0_Rap_%s_%s", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()), Form("D0 y %s;y;Entries", gen2.gen_or_det.c_str()), 100, -5, 5); //-TMath::Pi(), 2*TMath::Pi());

    TH1I *hD0_MPID_1 = new TH1I(Form("hD0_MPID_%s_%s", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()), Form("D0 Mother PID %s;MPID;Entries", gen1.gen_or_det.c_str()), 5000, -2500, 2500); //-TMath::Pi(), 2*TMath::Pi());
    TH1I *hD0_MPID_2 = new TH1I(Form("hD0_MPID_%s_%s", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()), Form("D0 Mother PID %s;MPID;Entries", gen2.gen_or_det.c_str()), 5000, -2500, 2500); //-TMath::Pi(), 2*TMath::Pi());

    TH1I *hnumD0s_1 = new TH1I(Form("hnumD0s_%s_%s", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()), Form("D0 PID %s;PID;Entries", gen1.gen_or_det.c_str()), 10, 0, 10); //-TMath::Pi(), 2*TMath::Pi());
    TH1I *hnumD0s_2 = new TH1I(Form("hnumD0s_%s_%s", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()), Form("D0 PID %s;PID;Entries", gen2.gen_or_det.c_str()), 10, 0, 10); //-TMath::Pi(), 2*TMath::Pi());


    // -------- FILL --------
    cout << "filling first file particle hists " << endl;
    fillParticleHistsFromChain(chain1_particle, hPt_1, hEta_1, hPhi_1, hPID_1, false, false);

    cout << "filling second file particle hists " << endl;
    fillParticleHistsFromChain(chain2_particle, hPt_2, hEta_2, hPhi_2, hPID_2, false, false);

    cout << "filling first file D0 hists " << endl;
    if ( gen1.gen_or_det == "gen" ) fillgenD0HistsFromChain(chain1_D0, hD0_Pt_1, hD0_Eta_1, hD0_Phi_1, hD0_Rap_1, hD0_MPID_1, hnumD0s_1);
    // else if ( gen1.gen_or_det == "det" ) fillDetD0HistsFromChain();

    cout << "filling second file D0 hists " << endl;
    if ( gen2.gen_or_det == "gen" ) fillgenD0HistsFromChain(chain2_D0, hD0_Pt_2, hD0_Eta_2, hD0_Phi_2, hD0_Rap_2, hD0_MPID_2, hnumD0s_2);
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
    auto savePair = [](TFile * fout_root, Generator gen1, Generator gen2, TH1 *h1, TH1 *h2, std::string name) {
        
        TH1D *h_ratio = (TH1D *)h2->Clone(Form("h_ratio_%s", name.c_str()));
        h_ratio->Divide(h1);

        h_ratio->SetTitle(Form("Ratio %s", name.c_str()));
        h_ratio->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
        h_ratio->GetYaxis()->SetTitle("NON-PROMPT / PROMPT");

        // Save to root file
        fout_root->cd();
        h1->Write();
        h2->Write();
        h_ratio->Write();
    };

    savePair(fout_root, gen1, gen2, hPt_1,  hPt_2, "Pt" + gen1.gen_or_det);
    savePair(fout_root, gen1, gen2, hEta_1, hEta_2, "Eta" + gen1.gen_or_det);
    savePair(fout_root, gen1, gen2, hPhi_1, hPhi_2, "Phi" + gen1.gen_or_det);
    savePair(fout_root, gen1, gen2, hPID_1, hPID_2, "PID" + gen1.gen_or_det); //savePairI

    savePair(fout_root, gen1, gen2, hD0_Pt_1,  hD0_Pt_2, "D0_Pt" + gen1.gen_or_det);
    savePair(fout_root, gen1, gen2, hD0_Eta_1, hD0_Eta_2, "D0_Eta" + gen1.gen_or_det);
    savePair(fout_root, gen1, gen2, hD0_Phi_1, hD0_Phi_2, "D0_Phi" + gen1.gen_or_det);
    savePair(fout_root, gen1, gen2, hD0_Rap_1, hD0_Rap_2, "D0_Rap" + gen1.gen_or_det);
    savePair(fout_root, gen1, gen2, hD0_MPID_1, hD0_MPID_2, "D0_MPID" + gen1.gen_or_det); //savePairI
    savePair(fout_root, gen1, gen2, hnumD0s_1, hnumD0s_2, "numD0s" + gen1.gen_or_det); //savePairI
}


void compareParticleBranches_WithCS_TChain(std::ofstream &outfile, TFile * fout_root, Generator gen1, Generator gen2) {

    // -------- CHAINS --------
    std::vector<TChain *> chains1_particle = makeChain_WithCS(gen1, "particle"); //pythia prompt
    std::vector<TChain *> chains2_particle = makeChain_WithCS(gen2, "particle"); //pythia nonprompt
    std::vector<TChain *> chains1_D0 = makeChain_WithCS(gen1, "D0"); //pythia prompt
    std::vector<TChain *> chains2_D0 = makeChain_WithCS(gen2, "D0"); //pythia nonprompt

    // -------- save number of entries to a file --------
    outfile << "Number of entries in particle " << gen1.label << ": " << chains1_particle[0]->GetEntries();
    for (int i=1; i<10; i++) outfile << " " << chains1_particle[i]->GetEntries();
    outfile << endl;

    outfile << "Number of entries in particle " << gen2.label << ": " << chains2_particle[0]->GetEntries();
    for (int i=1; i<10; i++) outfile << " " << chains2_particle[i]->GetEntries();
    outfile << endl;

    outfile << "Number of entries in D0 " << gen1.label << ": " << chains1_D0[0]->GetEntries();
    for (int i=1; i<10; i++) outfile << " " << chains1_D0[i]->GetEntries();
    outfile << endl;

    outfile << "Number of entries in D0 " << gen2.label << ": " << chains2_D0[0]->GetEntries();
    for (int i=1; i<10; i++) outfile << " " << chains2_D0[i]->GetEntries();
    outfile << endl;

    // -------- HISTOGRAMS --------
    std::vector<TH1D *> vec_hPt_1;
    std::vector<TH1D *> vec_hPt_2;
    // std::vector<TH1D *> vec_hPt_onlycharged_1;
    // std::vector<TH1D *> vec_hPt_onlycharged_2;
    for ( int i = 0; i < 10; i++ ) {
        TH1D *hPt_temp_1 = new TH1D(Form("hPt_%s_%s_bin%d", gen1.gen_type.c_str(), gen1.gen_or_det.c_str(), i),  Form("Particle p_{T} %s;p_{T};#frac{d#sigma}{dp_{T}}", gen1.gen_or_det.c_str()), 200, 0, 200);
        TH1D *hPt_temp_2 = new TH1D(Form("hPt_%s_%s_bin%d", gen2.gen_type.c_str(), gen2.gen_or_det.c_str(), i),  Form("Particle p_{T} %s;p_{T};#frac{d#sigma}{dp_{T}}", gen2.gen_or_det.c_str()), 200, 0, 200);
        vec_hPt_1.push_back(hPt_temp_1);
        vec_hPt_2.push_back(hPt_temp_2);

        // TH1D *hPt_onlychargedtemp_1 = new TH1D(Form("hPt_onlycharged%s_%s_bin%d", gen1.gen_type.c_str(), gen1.gen_or_det.c_str(), i),  Form("Particle p_{T} %s;p_{T};#frac{d#sigma}{dp_{T}}", gen1.gen_or_det.c_str()), 200, 0, 200);
        // TH1D *hPt_onlychargedtemp_2 = new TH1D(Form("hPt_onlycharged%s_%s_bin%d", gen2.gen_type.c_str(), gen2.gen_or_det.c_str(), i),  Form("Particle p_{T} %s;p_{T};#frac{d#sigma}{dp_{T}}", gen2.gen_or_det.c_str()), 200, 0, 200);
        // vec_hPt_onlycharged_1.push_back(hPt_onlychargedtemp_1);
        // vec_hPt_onlycharged_2.push_back(hPt_onlychargedtemp_2);
    }
    
    // D0 hists
    std::vector<TH1D *> vec_hD0_Pt_1;
    std::vector<TH1D *> vec_hD0_Pt_2;
    for ( int i = 0; i < 10; i++ ) {
        TH1D *hD0_Pt_temp_1 = new TH1D(Form("hD0_Pt_%s_%s_bin%d", gen1.gen_type.c_str(), gen1.gen_or_det.c_str(), i),  Form("D0 p_{T} %s;p_{T};#frac{d#sigma}{dp_{T}}", gen1.gen_or_det.c_str()), 100, 0, 100);
        TH1D *hD0_Pt_temp_2 = new TH1D(Form("hD0_Pt_%s_%s_bin%d", gen2.gen_type.c_str(), gen2.gen_or_det.c_str(), i),  Form("D0 p_{T} %s;p_{T};#frac{d#sigma}{dp_{T}}", gen2.gen_or_det.c_str()), 100, 0, 100);
        vec_hD0_Pt_1.push_back(hD0_Pt_temp_1);
        vec_hD0_Pt_2.push_back(hD0_Pt_temp_2);
    }

    std::vector<TH1D *> vec_hD0_Rap_1;
    std::vector<TH1D *> vec_hD0_Rap_2;
    for ( int i = 0; i < 10; i++ ) {
        TH1D *hD0_Rap_temp_1 = new TH1D(Form("hD0_Rap_%s_%s_bin%d", gen1.gen_type.c_str(), gen1.gen_or_det.c_str(), i),  Form("D0 y %s;y;#frac{d#sigma}{dy}", gen1.gen_or_det.c_str()), 100, -5, 5);
        TH1D *hD0_Rap_temp_2 = new TH1D(Form("hD0_Rap_%s_%s_bin%d", gen2.gen_type.c_str(), gen2.gen_or_det.c_str(), i),  Form("D0 y %s;y;#frac{d#sigma}{dy}", gen2.gen_or_det.c_str()), 100, -5, 5);
        vec_hD0_Rap_1.push_back(hD0_Rap_temp_1);
        vec_hD0_Rap_2.push_back(hD0_Rap_temp_2);
    }

    // -------- FILL --------
    cout << "filling first and second file particle hists " << endl;
    TH1D * hdummy_d;
    TH1I * hdummy_i;
    for ( int i = 0; i < 10; i++ ) {
        fillParticleHistsFromChain(chains1_particle[i], vec_hPt_1[i], hdummy_d, hdummy_d, hdummy_i, true, false); // looking at only charged!
        fillParticleHistsFromChain(chains2_particle[i], vec_hPt_2[i], hdummy_d, hdummy_d, hdummy_i, true, false); // looking at only charged!

        // // charged particles only
        // fillParticleHistsFromChain(chains1_particle[i], vec_hPt_onlycharged_1[i], hdummy_d, hdummy_d, hdummy_i, true, false);
        // fillParticleHistsFromChain(chains2_particle[i], vec_hPt_onlycharged_2[i], hdummy_d, hdummy_d, hdummy_i, true, false);
    }

    cout << "filling first and second file D0 hists " << endl;
    for ( int i = 0; i < 10; i++ ) {
        if ( gen1.gen_or_det == "gen" ) {
            fillgenD0HistsFromChain(chains1_D0[i], vec_hD0_Pt_1[i], hdummy_d, hdummy_d, vec_hD0_Rap_1[i], hdummy_i, hdummy_i, true);
            fillgenD0HistsFromChain(chains2_D0[i], vec_hD0_Pt_2[i], hdummy_d, hdummy_d, vec_hD0_Rap_2[i], hdummy_i, hdummy_i, true);
        }
    }

    // -------- SCALE --------
    gen1.scale_hists(vec_hPt_1);
    gen2.scale_hists(vec_hPt_2);
    gen1.scale_hists(vec_hD0_Pt_1);
    gen2.scale_hists(vec_hD0_Pt_2);
    gen1.scale_hists(vec_hD0_Rap_1);
    gen2.scale_hists(vec_hD0_Rap_2);

    // gen1.scale_hists(vec_hPt_onlycharged_1);
    // gen2.scale_hists(vec_hPt_onlycharged_2);

    // -------- ADD HISTOGRAMS --------
    TH1D * hPt_comb_1 = addHists(vec_hPt_1, Form("hPt_%s_%s_crosssection", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()) );
    TH1D * hPt_comb_2 = addHists(vec_hPt_2, Form("hPt_%s_%s_crosssection", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()) );
    TH1D * hD0_Pt_comb_1 = addHists(vec_hD0_Pt_1, Form("hD0_Pt_%s_%s_crosssection", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()) );
    TH1D * hD0_Pt_comb_2 = addHists(vec_hD0_Pt_2, Form("hD0_Pt_%s_%s_crosssection", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()) );
    TH1D * hD0_Rap_comb_1 = addHists(vec_hD0_Rap_1, Form("hD0_Rap_%s_%s_crosssection", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()) );
    TH1D * hD0_Rap_comb_2 = addHists(vec_hD0_Rap_2, Form("hD0_Rap_%s_%s_crosssection", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()) );

    // TH1D * hPt_onlycharged_comb_1 = addHists(vec_hPt_onlycharged_1, Form("hPt_onlycharged_%s_%s_crosssection", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()) );
    // TH1D * hPt_onlycharged_comb_2 = addHists(vec_hPt_onlycharged_2, Form("hPt_onlycharged_%s_%s_crosssection", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()) );

    // -------- STYLE --------
    hPt_comb_1->SetLineColor(kRed);
    hPt_comb_2->SetLineColor(kBlue);
    hD0_Pt_comb_1->SetLineColor(kRed);
    hD0_Pt_comb_2->SetLineColor(kBlue);
    hD0_Rap_comb_1->SetLineColor(kRed);
    hD0_Rap_comb_2->SetLineColor(kBlue);

    // hPt_onlycharged_comb_1->SetLineColor(kRed);
    // hPt_onlycharged_comb_2->SetLineColor(kBlue);

    // -------- DRAW --------
    auto savePair = [](TFile * fout_root, Generator gen1, Generator gen2, TH1 *h1, TH1 *h2, std::string name) {
        
        TH1D *h_ratio = (TH1D *)h2->Clone(Form("h_ratio_%s_crosssection", name.c_str()));
        h_ratio->Divide(h1);

        h_ratio->SetTitle(Form("Ratio %s", name.c_str()));
        h_ratio->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
        h_ratio->GetYaxis()->SetTitle("NON-PROMPT / PROMPT");

        // Save to root file
        fout_root->cd();
        h1->Write();
        h2->Write();
        h_ratio->Write();
    };

    savePair(fout_root, gen1, gen2, hPt_comb_1,  hPt_comb_2, "Pt" + gen1.gen_or_det);
    savePair(fout_root, gen1, gen2, hD0_Pt_comb_1,  hD0_Pt_comb_2, "D0_Pt" + gen1.gen_or_det);
    savePair(fout_root, gen1, gen2, hD0_Rap_comb_1,  hD0_Rap_comb_2, "D0_Rap" + gen1.gen_or_det);

    // savePair(fout_root, gen1, gen2, hPt_onlycharged_comb_1,  hPt_onlycharged_comb_2, "Pt_onlycharged" + gen1.gen_or_det);
}



void compareParticleBranches_WithCS_TChain_Method2(TFile * fout_root, Generator gen1, Generator gen2) {

    // -------- HISTOGRAMS --------
    TH1D *hPt_temp_1_method2 = new TH1D(Form("hPt_%s_%s_crosssection_method2", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()),  Form("Particle p_{T} %s;p_{T};#frac{d#sigma}{dp_{T}}", gen1.gen_or_det.c_str()), 200, 0, 200);
    TH1D *hPt_temp_2_method2 = new TH1D(Form("hPt_%s_%s_crosssection_method2", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()),  Form("Particle p_{T} %s;p_{T};#frac{d#sigma}{dp_{T}}", gen2.gen_or_det.c_str()), 200, 0, 200);
    
    // D0 hists
    TH1D *hD0_Pt_temp_1_method2 = new TH1D(Form("hD0_Pt_%s_%s_crosssection_method2", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()),  Form("D0 p_{T} %s;p_{T};#frac{d#sigma}{dp_{T}}", gen1.gen_or_det.c_str()), 100, 0, 100);
    TH1D *hD0_Pt_temp_2_method2 = new TH1D(Form("hD0_Pt_%s_%s_crosssection_method2", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()),  Form("D0 p_{T} %s;p_{T};#frac{d#sigma}{dp_{T}}", gen2.gen_or_det.c_str()), 100, 0, 100);
    
    TH1D *hD0_Rap_temp_1_method2 = new TH1D(Form("hD0_Rap_%s_%s_crosssection_method2", gen1.gen_type.c_str(), gen1.gen_or_det.c_str()),  Form("D0 y %s;y;#frac{d#sigma}{dy}", gen1.gen_or_det.c_str()), 100, -5, 5);
    TH1D *hD0_Rap_temp_2_method2 = new TH1D(Form("hD0_Rap_%s_%s_crosssection_method2", gen2.gen_type.c_str(), gen2.gen_or_det.c_str()),  Form("D0 y %s;y;#frac{d#sigma}{dy}", gen2.gen_or_det.c_str()), 100, -5, 5);
    

    // -------- FILL HISTOGRAMS --------
    TH1D * hdummy;
    fillHistsPerFile(gen1, "particle", hPt_temp_1_method2, hdummy, true);
    fillHistsPerFile(gen2, "particle", hPt_temp_2_method2, hdummy, true);
    fillHistsPerFile(gen1, "D0", hD0_Pt_temp_1_method2, hD0_Rap_temp_1_method2, false);
    fillHistsPerFile(gen2, "D0", hD0_Pt_temp_2_method2, hD0_Rap_temp_2_method2, false);

    // -------- NORMALIZE HISTOGRAMS --------
    hPt_temp_1_method2->Scale(1.0, "width");
    hPt_temp_2_method2->Scale(1.0, "width");
    hD0_Pt_temp_1_method2->Scale(1.0, "width");
    hD0_Pt_temp_2_method2->Scale(1.0, "width");
    hD0_Rap_temp_1_method2->Scale(1.0, "width");
    hD0_Rap_temp_2_method2->Scale(1.0, "width");


    // -------- STYLE --------
    hPt_temp_1_method2->SetLineColor(kRed);
    hPt_temp_2_method2->SetLineColor(kBlue);
    hD0_Pt_temp_1_method2->SetLineColor(kRed);
    hD0_Pt_temp_2_method2->SetLineColor(kBlue);
    hD0_Rap_temp_1_method2->SetLineColor(kRed);
    hD0_Rap_temp_2_method2->SetLineColor(kBlue);

    // -------- SAVE --------
    auto savePair = [](TFile * fout_root, Generator gen1, Generator gen2, TH1 *h1, TH1 *h2, std::string name) {
        
        TH1D *h_ratio = (TH1D *)h2->Clone(Form("h_ratio_%s_crosssection_method2", name.c_str()));
        h_ratio->Divide(h1);

        h_ratio->SetTitle(Form("Ratio %s", name.c_str()));
        h_ratio->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
        h_ratio->GetYaxis()->SetTitle("NON-PROMPT / PROMPT");

        // Save to root file
        fout_root->cd();
        h1->Write();
        h2->Write();
        h_ratio->Write();
    };

    savePair(fout_root, gen1, gen2, hPt_temp_1_method2,  hPt_temp_2_method2, "Pt" + gen1.gen_or_det);
    savePair(fout_root, gen1, gen2, hD0_Pt_temp_1_method2,  hD0_Pt_temp_2_method2, "D0_Pt" + gen1.gen_or_det);
    savePair(fout_root, gen1, gen2, hD0_Rap_temp_1_method2,  hD0_Rap_temp_2_method2, "D0_Rap" + gen1.gen_or_det);

}


void extract_pythia_particle_level_histograms(const char *opts = "") {

    TH1::SetDefaultSumw2(kTRUE);

    // ------- CHOOSE PYTHIA OR HERWIG -------
    TString options(opts);
    std::string generator_choice;
    if (options.Contains("h")) {
        generator_choice = "herwig"; // "herwig"
    } else if (options.Contains("p")) {
        generator_choice = "pythia"; // "pythia"
    } else {
        std::cerr << "Invalid option! Use 'h' for herwig or 'p' for pythia as per instructions at top of this file." << std::endl;
        return;
    }

    bool option_notforcedD0toKPi;
    std::string str_notforcedD0toKPi;
    if (options.Contains("x")) {
        option_notforcedD0toKPi = true;
        str_notforcedD0toKPi = "_notforcedD0toKPi"; 
    } else {
        option_notforcedD0toKPi = false; 
        str_notforcedD0toKPi = ""; 
    }

    std::string basepath;
    if ( generator_choice == "pythia" ) basepath = "/global/cfs/cdirs/alice/blianggi";
    else if ( generator_choice == "herwig" ) basepath = "/software/users/blianggi";

    // -------- INPUT DIRECTORIES --------
    // post eff smearing -- generator + detector level
    std::string pythia_prompt_filepaths = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/pythiagen/tree_fastsim/45178629/45154942/files.txt"; 
    std::string pythia_nonprompt_filepaths = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/pythiagen/tree_fastsim/46306341/46293548/files.txt"; 
    std::string herwig_prompt_filepaths = "/rstorage/generators/herwig_alice/tree_fastsim/492678/299990/files.txt"; 
    std::string herwig_nonprompt_filepaths = "/rstorage/generators/herwig_alice/tree_fastsim/516788/515788/files.txt"; 
    std::string pythia_prompt_notforcedD0toKPi_filepaths = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/pythiagen/tree_gen/49538712/49538712/files.txt"; 
    std::string pythia_nonprompt_notforcedD0toKPi_filepaths = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/pythiagen/tree_gen/49538735/49538735/files.txt"; 

    std::string pythia_prompt_scalefactor_filepaths = basepath + "/mypyjetty/analysis/scalefactors/PYTHIA_fastsim_HF_scaleFactors.yaml"; // only on perlmutter
    std::string pythia_nonprompt_scalefactor_filepaths = basepath + "/mypyjetty/analysis/scalefactors/PYTHIA_fastsim_nonprompt_D0_scaleFactors.yaml"; // only on perlmutter
    std::string herwig_prompt_scalefactor_filepaths = basepath + "/mypyjetty/analysis/scalefactors/herwig_HF_299990_scaleFactors.yaml"; // only on hiccup
    std::string herwig_nonprompt_scalefactor_filepaths = basepath + "/mypyjetty/analysis/scalefactors/herwig_bbbar_515788_scaleFactors.yaml"; // only on hiccup
    std::string pythia_prompt_notforcedD0toKPi_scalefactor_filepaths = basepath + "/mypyjetty/analysis/scalefactors/PYTHIA_fastsim_HF_notforcedD0toKPi_49538712_scaleFactors.yaml";
    std::string pythia_nonprompt_notforcedD0toKPi_scalefactor_filepaths = basepath + "/mypyjetty/analysis/scalefactors/PYTHIA_fastsim_nonprompt_notforcedD0toKPi_49538735_scaleFactors.yaml";

    std::string pythia_prompt_sf_ind_filepaths = "/global/cfs/cdirs/alice/blianggi/mypyjetty/analysis/scalefactors/PYTHIA_fastsim_HF_45154942_individualScaleFactors.txt";
    std::string pythia_nonprompt_sf_ind_filepaths = "/global/cfs/cdirs/alice/blianggi/mypyjetty/analysis/scalefactors/PYTHIA_fastsim_nonprompt_D0_46293548_individualScaleFactors.txt";
    std::string herwig_prompt_sf_ind_filepaths = "/software/users/blianggi/mypyjetty/analysis/scalefactors/herwig_HF_299990_individualScaleFactors.txt";
    std::string herwig_nonprompt_sf_ind_filepaths = "/software/users/blianggi/mypyjetty/analysis/scalefactors/herwig_bbbar_515788_individualScaleFactors.txt";
    
    // -------- DEFINE GENERATOR --------
    Generator gen_pythia_prompt("pythia_prompt", pythia_prompt_filepaths, "tree_Particle_gen", "tree_D0_gen", "gen", "Pythia prompt, gen", pythia_prompt_scalefactor_filepaths, pythia_prompt_sf_ind_filepaths);
    Generator gen_pythia_nonprompt("pythia_nonprompt", pythia_nonprompt_filepaths, "tree_Particle_gen", "tree_D0_gen", "gen", "Pythia non-prompt, gen", pythia_nonprompt_scalefactor_filepaths, pythia_nonprompt_sf_ind_filepaths);
    Generator gen_herwig_prompt("herwig_prompt", herwig_prompt_filepaths, "tree_Particle_gen", "tree_D0_gen", "gen", "Herwig prompt, gen", herwig_prompt_scalefactor_filepaths, herwig_prompt_sf_ind_filepaths);
    Generator gen_herwig_nonprompt("herwig_nonprompt", herwig_nonprompt_filepaths, "tree_Particle_gen", "tree_D0_gen", "gen", "Herwig non-prompt, gen", herwig_nonprompt_scalefactor_filepaths, herwig_nonprompt_sf_ind_filepaths);

    Generator gen_pythia_prompt_notforcedD0toKPi("pythia_prompt", pythia_prompt_notforcedD0toKPi_filepaths, "PWGHF_TreeCreator/tree_Particle_gen", "PWGHF_TreeCreator/tree_D0_gen", "gen", "Pythia prompt no forced D0->KPi, gen", pythia_prompt_notforcedD0toKPi_scalefactor_filepaths, "");
    Generator gen_pythia_nonprompt_notforcedD0toKPi("pythia_nonprompt", pythia_nonprompt_notforcedD0toKPi_filepaths, "PWGHF_TreeCreator/tree_Particle_gen", "PWGHF_TreeCreator/tree_D0_gen", "gen", "Pythia non-prompt not forced D0->KPi, gen", pythia_nonprompt_notforcedD0toKPi_scalefactor_filepaths, "");
    
    // Generator det_pythia_prompt("pythia_prompt", pythia_prompt_filepaths, "tree_Particle", "tree_D0", "det", "Pythia prompt, det");
    // Generator det_pythia_nonprompt("pythia_nonprompt", pythia_nonprompt_filepaths, "tree_Particle", "tree_D0", "det", "Pythia non-prompt,  det");
    
    // -------- OPEN OUTPUT FILEs --------
    std::ofstream outfile(Form("%s/mypyjetty/storage/HF_EEC/plots/HF_particle_comparisons/number_of_entries_%s%s.txt",basepath.c_str(), generator_choice.c_str(), str_notforcedD0toKPi.c_str()));
    TFile * fout_root = new TFile(Form("%s/mypyjetty/storage/HF_EEC/rootfiles/HF_particle_comparisons/HF_particle_comparisons_%s%s.root", basepath.c_str(), generator_choice.c_str(), str_notforcedD0toKPi.c_str()), "RECREATE");

    // -------- COMPARE GENERATORS --------
    if (generator_choice == "pythia") {
        if (option_notforcedD0toKPi) compareParticleBranches_TChain(outfile, fout_root, gen_pythia_prompt_notforcedD0toKPi, gen_pythia_nonprompt_notforcedD0toKPi);
        else compareParticleBranches_TChain(outfile, fout_root, gen_pythia_prompt, gen_pythia_nonprompt);
        // compareParticleBranches_TChain(outfile, fout_root, det_anchmc, det_pythiafastsim);
    }
    else if (generator_choice == "herwig") {
        compareParticleBranches_TChain(outfile, fout_root, gen_herwig_prompt, gen_herwig_nonprompt);
    }

    // -------- GET CROSS SECTIONS PER FILE --------
    if (generator_choice == "pythia") {
        if (option_notforcedD0toKPi) compareParticleBranches_WithCS_TChain(outfile, fout_root, gen_pythia_prompt_notforcedD0toKPi, gen_pythia_nonprompt_notforcedD0toKPi);
        else compareParticleBranches_WithCS_TChain(outfile, fout_root, gen_pythia_prompt, gen_pythia_nonprompt);
        // compareParticleBranches_WithCS_TChain_Method2(fout_root, gen_pythia_prompt, gen_pythia_nonprompt); // this was done as a check - method 2 is longer but more robust
    } else if (generator_choice == "herwig") {
        compareParticleBranches_WithCS_TChain(outfile, fout_root, gen_herwig_prompt, gen_herwig_nonprompt);
        // compareParticleBranches_WithCS_TChain_Method2(fout_root, gen_herwig_prompt, gen_herwig_nonprompt); // this was done as a check - method 2 is longer but more robust
    }

    // CLOSE OUTPUT FILE
    outfile.close();
    fout_root->Close();

    // PRINT OUT CLOSING TIME
    std::time_t t = std::time(0);   // Get current time_t
    std::cout << std::ctime(&t);    // Convert to string and print
    
}


