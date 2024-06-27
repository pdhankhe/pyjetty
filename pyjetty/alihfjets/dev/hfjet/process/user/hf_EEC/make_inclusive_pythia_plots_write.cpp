


void make_inclusive_pythia_plots_write() {

    // files
    std::string infile_inclusive_2MC = "/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alihfjets/dev/hfjet/process/user/hf_EEC/mc_incl/AnalysisResultsFinal_inclusive_2MC.root";
    std::string infile_pythia_stuffs = "/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alihfjets/dev/hfjet/process/user/hf_EEC/plots/final/cgl/AnalysisResultsFinal_fromcgl.root";

    TFile* f_inclusive_2MC = new TFile(infile_inclusive_2MC.c_str(), "READ");
    TFile* f_pythia_stuffs = new TFile(infile_pythia_stuffs.c_str(), "READ");

    // Output directory
    std::string outdir = "plots/mc_incl/";
    // Output file for binned results
    std::string outfile = outdir + "AnalysisResultsFinal_inclusive_allMC" + ".root"; 
    TFile* f_out = new TFile(outfile.c_str(), "RECREATE");

    const int pt_bins[] = { 10, 15, 30 }; //{ 7, 10, 15, 30, 50, 70 }; //{ 10, 20, 40 }; // CHANGE HERE!!
    const int d0_pt_cuts[] = { 3, 5, 5 }; //{ 3, 5, 5, 5, 5 };
    const int n_bins = 2; //5; 
    for (int i = 0; i < n_bins; i++) {
        cout << "in pt bin" << i << endl;
        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];
        
        //hist names
        std::string i_herwig_leadcut5_name = "i_herwig_leadcut5_EEC_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_R0.4_trkthrd1.0";
        std::string i_sherpa_ahadic_leadcut5_name = "i_sherpa_jetpt" + std::to_string(pt_min) + "_leadcut5_EEC_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_R0.4_trkthrd1.0";
        std::string i_sherpa_lund_leadcut5_name = "i_sherpa_jetpt" + std::to_string(pt_min) + "_lund_leadcut5_EEC_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_R0.4_trkthrd1.0";
        std::string i_pythia_leadcut5_name = "h_charmdecaysON_leadpt5_i_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);

        TH1D* hi_herwig_leadcut5 = (TH1D*) f_inclusive_2MC->Get(i_herwig_leadcut5_name.c_str());
        TH1D* hi_sherpa_ahadic_leadcut5 = (TH1D*) f_inclusive_2MC->Get(i_sherpa_ahadic_leadcut5_name.c_str());
        TH1D* hi_sherpa_lund_leadcut5 = (TH1D*) f_inclusive_2MC->Get(i_sherpa_lund_leadcut5_name.c_str());
        TH1D* hi_pythia_leadcut5 = (TH1D*) f_pythia_stuffs->Get(i_pythia_leadcut5_name.c_str());


        hi_herwig_leadcut5->SetTitle(i_herwig_leadcut5_name.c_str());
        hi_sherpa_ahadic_leadcut5->SetTitle(i_sherpa_ahadic_leadcut5_name.c_str());
        hi_sherpa_lund_leadcut5->SetTitle(i_sherpa_lund_leadcut5_name.c_str());
        hi_pythia_leadcut5->SetTitle(("i_pythia_leadcut5_EEC_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_R0.4_trkthrd1.0").c_str());

        hi_pythia_leadcut5->SetName(("i_pythia_leadcut5_EEC_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_R0.4_trkthrd1.0").c_str());

        // write to output
        f_out->cd();
        hi_herwig_leadcut5->Write();
        hi_sherpa_ahadic_leadcut5->Write();
        hi_sherpa_lund_leadcut5->Write();
        hi_pythia_leadcut5->Write();


    }

    f_inclusive_2MC->Close();
    delete f_inclusive_2MC;
    f_pythia_stuffs->Close();
    delete f_pythia_stuffs;

    return;

}