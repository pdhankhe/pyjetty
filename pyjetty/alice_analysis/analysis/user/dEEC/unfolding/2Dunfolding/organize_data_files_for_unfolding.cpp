


//===========================================================================
// This function extracts the observable histogram, changes the name, 
// and saves it to its own file in the correct directory.
//===========================================================================
void extract_2DRM(TFile *fin, TFile *fout, std::string observable, int ptbin, int RLbin) {

    // TODO: still need to make these 2d histograms!!
    std::string original_histname = Form("h_?_R0.4_1.0");

    fin->cd();
    TH2D * hist = (TH2D *)gDirectory->Get(original_histname.c_str());

    // change the name - { h_[obs]_JetPt_R[R]_[subobs]_[grooming setting] }
    std::string new_histname = Form("h_%s_JetPt_R0.4_1.0", observable);
    hist->SetNameTitle();

    fout->cd();
    hist->Write();

}


//===========================================================================
// This file takes in the following arguments:
// observable, PTBIN, RLBIN
// THIS IS GOOD FOR HICCUP
//===========================================================================
void organize_data_files_for_unfolding() {

    // Reading arguments
    std::string observable = argv[1];  // Keep argument as a string
    int ptbin = atoi(argv[2]);     // Convert argument to an int
    int rlbin = atoi(argv[3]);      // Convert argument to an int
   
    // TODO: change the input files and name!
    std::string base_slurmoutput_path = "/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/data_firstattempt/";
    std::string infile = base_slurmoutput_path + "DataHists.root";

    std::string new_filedir = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/unfolding/";
    std::string new_filedir_addition = Form("%s/PTBIN%d/RLBIN%d/", observable, ptbin, rlbin);
    std::string new_filename = new_filedir + new_filedir_addition + "AnalysisResults_data.root";

    TFile* root_infile = new TFile(infile.c_str(), "READ");
    TFile* root_outfile = new TFile(new_filename.c_str(), "RECREATE");

    extract_obshist(root_infile, root_outfile, observable, ptbin, rlbin);

}