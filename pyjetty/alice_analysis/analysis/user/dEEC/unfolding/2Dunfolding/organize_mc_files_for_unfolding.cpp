

//===========================================================================
// This function extracts the thnsparse response matrix, changes the name, 
// and saves it to its own file in the correct directory.
//===========================================================================
void extract_2DRM(TFile *fin, TFile *fout, std::string observable, int ptbin, int RLbin) {

    std::string original_histname = Form("hResponse_JetPt_%s_PTBIN%d_RLBIN%d_R0.4_1.0", observable, ptbin, RLbin);

    fin->cd();
    THnSparse * hist = (THnSparse *)gDirectory->Get(original_histname.c_str());

    // change the name - { hResponse_JetPt_[obs]_R[R]_[subobs]_[grooming setting] }
    std::string new_histname = Form("hResponse_JetPt_%s_R0.4_1.0", observable);
    hist->SetNameTitle();

    fout->cd();
    hist->Write();

}

//===========================================================================
// This file takes in the following arguments:
// observable, PTBIN, RLBIN
// THIS IS GOOD FOR PERLMUTTER
//===========================================================================
void organize_mc_files_for_unfolding(int argc, char **argv) {

    // Reading arguments
    std::string observable = argv[1];  // Keep argument as a string
    int ptbin = atoi(argv[2]);     // Convert argument to an int
    int rlbin = atoi(argv[3]);      // Convert argument to an int
   
    // std::string base_slurmoutput_path = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/33860840/";
    std::string base_slurmoutput_path = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/33818618/scaling"; //wrong for now   std::string infile = base_slurmoutput_path + "some_merged_filename.root";
    std::string infile = base_slurmoutput_path + "AnalysisResultsFinal.root";

    std::string new_filedir = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/unfolding/";
    std::string new_filedir_addition = Form("%s/PTBIN%d/RLBIN%d/", observable, ptbin, rlbin);
    std::string new_filename = new_filedir + new_filedir_addition + "AnalysisResults_Response.root";

    TFile* root_infile = new TFile(infile.c_str(), "READ");
    TFile* root_outfile = new TFile(new_filename.c_str(), "RECREATE");

    extract_2DRM(root_infile, root_outfile, observable, ptbin, rlbin);


}