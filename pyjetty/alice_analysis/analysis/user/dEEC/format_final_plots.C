

bool jetpt_bool = true;

bool deltap_bool = true;
bool p_bool = true;
bool deltajt_bool = true;
bool jt_bool = true;
bool ew_bool = true;
bool twoDhists_bool = false;
bool rc_bool = false;

std::string output_dir = "/software/users/blianggi/mypyjetty/storage/dEEC/plots/data_fourthattempt_ptrlbins/rebinx4/";


void plot_observable(std::string observable) {
    // cycle through pt and rl bins
    // cut on the ∆p axis to limit the range to appropriate amount
}

void format_final_plots() {

    // take input file
    TString input_file_name = "/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/data_fourthattempt_ptrlbins/rebinx4/DataHists.root";
    TFile * input_file = new TFile(input_file_name, "READ");

    // for each observable, plot 
    if (deltap_bool) plot_observable("deltap");

}