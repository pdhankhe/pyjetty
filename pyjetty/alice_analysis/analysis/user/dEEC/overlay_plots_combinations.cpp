// ROOT macro that can take in different files and overlay histograms of the same name
// on top of each other 
// Beatrice Liang-Gilman (beatrice_lg@berkeley.edu)


// global variables
Double_t colors[16] = {kGray, kMagenta, kGreen+2, kBlue, kOrange+1, kViolet+1, kRed, kYellow+1, kCyan+1};
Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
Double_t marker_size = 1.5;

std::string attempt_dir = "raw_data__data_correctionfactors_comparison";
std::string outdir = "/software/users/blianggi/mypyjetty/storage/dEEC/plots/" + attempt_dir;

// comp_cases
// 1: hic, correcting with EEC fcorr; a bit out of date now
// 2: perly, bin by bin corrections -- matched, rebinx4 -- not tested
// 3: hic, pythia vs data
int comp_case = 2;

// setting style
void SetStyle(Bool_t graypalette=true) {
  	cout << "Setting style!" << endl;
  
  	gStyle->Reset("Plain");
  	gStyle->SetOptTitle(0);
  	gStyle->SetOptStat(0);
  	// if(graypalette) gStyle->SetPalette(8,0);
  	// else gStyle->SetPalette(1);
    gStyle->SetPalette(kRainbow);
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


// prepare observable canvas
TCanvas * prepare_canvas(bool logx=false, bool logy=false) {
    TCanvas *can = new TCanvas();
    can->cd();

    if (logx) gPad->SetLogx();
    if (logy) gPad->SetLogy();

    return can;
}

// function to plot a histogram on a given canvas
void plot_histogram_on_given_canvas(TCanvas *can, TH1D *hist, int markerstyle, int markercolor) {

    // format hist with appropriate marker style
    hist->SetMarkerStyle(markerstyle);
    int color = hist->GetMarkerColor();
    cout << "color " << color << endl;
    hist->SetMarkerColorAlpha(markercolor, 0.75);
    // int color2 = hist->GetMarkerColor();
    // cout << "color2 " << color2 << endl;
    // hist->SetMarkerSize(1.5);

    // draw
    can->cd();
    hist->Draw("SAME");
    // can->Update();

    // delete hist;

    return;
}

void addleg_and_save_combined_hists(TCanvas *can_all, TLegend *l, 
                          std::string obs_name, std::string ptname, 
                          std::string norm_string, std::string hist_addname) {

    can_all->cd();
    l->Draw("same");

    //save as PDF
    std::string add_dir = "/" + ptname + "/" + norm_string + "/" + obs_name;
    std::string fname_out = outdir + add_dir + "/corrhist_" + obs_name + "_ALL" + hist_addname + ".pdf";
    can_all->SaveAs(fname_out.c_str());

    delete can_all;
    cout << "canvas deleted! " << endl;

}


// analyze one observable
void plot_one_observable(TFile *file1, TFile *file2, std::string obsname, int markerstyle1, int markerstyle2, bool logx, bool logy,
                         const int pt_bins[], int n_bins, const double RL_bins[][8], int n_RLbins, bool include_RL0, bool include_RL1) {

    std::string norms[] = { "unnormalized", "self_normalized", "norm_by_jets" };
    std::string weightstr = "";
    std::string jetRname = Form("_R%s", "0.4"); //jetR.c_str());
    std::string thrname = Form("_t%s", "1.0"); //threshold.c_str());
    
    for (std::string norm : norms) { // normalization loop

        for ( int i = 0; i < n_bins; i++) { // jet pt bins loop

            cout << "in pt bin" << i << endl;
            int pt_min = pt_bins[i];
            int pt_max = pt_bins[i+1];
            std::string ptname = to_string(pt_min) + "-" + to_string(pt_max);

            TCanvas *can = prepare_canvas(logx, logy);
            TLegend *leg = new TLegend(0.6, 0.6, 0.85, 0.87);
            leg->SetTextSize(0.037);
            leg->SetBorderSize(0);
            vector<TH1D*> hist_vec; //save all hists to this vector, where even indices are file 1 and odds are file 2

            for ( int j = 0; j < n_RLbins; j++) { // RL bins loop
            // for ( int j = n_RLbins - 1; j >= 0; j--) { // RL bins loop

                cout << "in RL bin" << j << endl;
                if (!include_RL0) {
                    if (j == 0) continue; // can add something here to change the filename for ALL
                }
                if (!include_RL1 && j == n_RLbins-1) continue;

                double RL_min = RL_bins[i][j];
                double RL_max = RL_bins[i][j+1];
                std::string RLname = Form("_RL%.3f-%.3f", RL_min, RL_max);
                std::string RLname_leg = Form("RL = %.3f-%.3f", RL_min, RL_max);
                // std::string hist_addname = weightstr + jetRname + thrname + "_pt" + ptname + RLname;
                cout << " RL min: " << RL_min << " // RL max: " << RL_max << endl;


                std::string histname = Form("h_%s_R0.4_t1.0_pt%d-%d_RL%.3f-%.3f_%s", obsname.c_str(), pt_min, pt_max, RL_min, RL_max, norm.c_str());

                TH1D * hist_file1 = (TH1D*) file1->Get(histname.c_str());
                TH1D * hist_file2 = (TH1D*) file2->Get(histname.c_str());

                hist_vec.push_back((TH1D*) hist_file1->Clone());
                hist_vec.push_back((TH1D*) hist_file2->Clone());

                leg->AddEntry(hist_file1, RLname_leg.c_str(), "pl");
                
                // plot_histogram_on_given_canvas(can, hist_file1, markerstyle1, colors[j]);
                // plot_histogram_on_given_canvas(can, hist_file2, markerstyle2, colors[j]);

            } // end RL bins loop

            // draw all plots
            for ( int j = 0; j < n_RLbins; j++) { 
                
                int k = j;
                if (!include_RL0) {
                    k = j-1;
                    if (j == 0) continue; // can add something here to change the filename for ALL
                }
                if (!include_RL1 && j == n_RLbins-1) continue;

                if ( k == 0) {
                    // set maximum based on maximum of all curves
                    double max = 0;
                    for (int a=0; a<hist_vec.size(); a++) {
                        double max_cand = hist_vec[a]->GetMaximum();
                        if (max_cand > max) max = max_cand;
                    }
                    // cout << "max is " << max << " which goes to " << max*1.5 << endl;
                    hist_vec[0]->SetMaximum( max * 1.5 );
                }

                plot_histogram_on_given_canvas(can, hist_vec[2*k], markerstyle1, colors[j]);
                plot_histogram_on_given_canvas(can, hist_vec[2*k+1], markerstyle2, colors[j]);
                cout << "plotted 2k = " << 2*k << " and color: " << colors[j] << " and markerstyle1: " << markerstyle1 << endl;
                cout << "plotted 2k+1 = " << 2*k + 1 << " and markerstyle2: " << markerstyle2 << endl;
            }
            
            // save canvas
            std::string hist_all_addname = weightstr + jetRname + thrname + "_pt" + ptname;
            addleg_and_save_combined_hists(can, leg, obsname, ptname, norm, hist_all_addname);

        } // end jet pt  bins loop

    } // end normalization loop
    
    


}

// function to plot raw data and corrected data together
void plot_comparisons(TFile *raw_data_infile, TFile *data_correctionfactors_infile, const int pt_bins[], int n_bins, 
                      const double RL_bins[][8], int n_RLbins, bool include_RL0, bool include_RL1) {
    // so we want to plot all the delta p curves, then all the delta pt and delta pl and EW curves. 2D hists don't need to be overlayed though
    // I should make projections of them at some point. I can also plot rc together.

    int markerstyle1 = kFullCircle;
    int markerstyle2 = kOpenCircle;
    // if (comp_case == 3) {
    //     markerstyle1 = 
    //
    // }
    
    // starting with deltap
    plot_one_observable(raw_data_infile, data_correctionfactors_infile, "deltap", markerstyle1, markerstyle2, false, true,
                         pt_bins, n_bins, RL_bins, n_RLbins, include_RL0, include_RL1);
    
    plot_one_observable(raw_data_infile, data_correctionfactors_infile, "deltapt", markerstyle1, markerstyle2, false, true,
                         pt_bins, n_bins, RL_bins, n_RLbins, include_RL0, include_RL1);
    plot_one_observable(raw_data_infile, data_correctionfactors_infile, "deltapl", markerstyle1, markerstyle2, false, true,
                         pt_bins, n_bins, RL_bins, n_RLbins, include_RL0, include_RL1);
    
    if (comp_case != 2) {
        plot_one_observable(raw_data_infile, data_correctionfactors_infile, "weights", markerstyle1, markerstyle2, false, true,
                            pt_bins, n_bins, RL_bins, n_RLbins, include_RL0, include_RL1);
    }

}






void overlay_plots_combinations() {

    // CONTROL VARIABLES HERE
    bool include_RL0 = false;
    bool include_RL1 = false;
    
    // setup variables
    bool debug = false;
    bool debug2 = false;


    // analysis variables
    const int pt_bins[] = { 20, 40, 60, 80 };
    const int n_bins = sizeof(pt_bins) / sizeof(pt_bins[0]) - 1; //3;
    
    const double RL_bins[3][8] = { { 0, 1e-2, 3e-2, 7e-2, 1.5e-1, 3e-1, 4e-1, 1 },
                            { 0, 1e-2, 2.5e-2, 4e-2, 8e-2, 2.5e-1, 4e-1, 1 },
                            { 0, 1e-2, 2.5e-2, 3e-2, 4.5e-2, 2e-1, 4e-1, 1 } };
    const int n_RLbins = sizeof(RL_bins[0]) / sizeof(RL_bins[0][0]) - 1; //gets the columns //6; //7; //5;

    if (debug2) cout << "pt_bins " << n_bins << " n_RLbins " << n_RLbins << endl;
    

    // filenames
    std::string root_indir_hic = "/software/users/blianggi/mypyjetty/storage/dEEC/rootfiles/";
    
    std::string raw_data_filename_hic = root_indir_hic + "data_firstattempt/DataHists.root"; //plots/ntuples/DataHists.root"; //FinalDataHists.root
    TFile* raw_data_infile_hic = new TFile(raw_data_filename_hic.c_str(), "READ");
    
    std::string data_correctionfactors_filename_hic = root_indir_hic + "data_correctionfactors/DataHists.root";
    TFile* data_correctionfactors_infile_hic = new TFile(data_correctionfactors_filename_hic.c_str(), "READ");

    //------------------------------------------------

    std::string root_indir_perly = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/rootfiles/";

    std::string raw_data_filename_perly = root_indir_perly + "data_firstattempt/DataHists.root"; //plots/ntuples/DataHists.root"; //FinalDataHists.root
    TFile* raw_data_infile_perly = new TFile(raw_data_filename_perly.c_str(), "READ");

    std::string binbybincorrections_filename_perly = root_indir_perly + "binbybincorrections/matched/rebinx4/DataHists_BinByBinCorr.root";
    TFile* binbybincorrections_infile_perly = new TFile(binbybincorrections_filename_perly.c_str(), "READ");


    // int num_input_files = 6;
    // std::bitset<num_input_files> pmask(0);

    // for (unsigned int i = 0; i < num_input_files; i++) {
    //     switch(i)
    //     {
    //         case kIgnore:		                    pmask[i] = true;
    //         case kAny: 			                    pmask[i] = true; 			break;
    //         case raw_data:		                    pmask[i] = true;            break;
    //         case data_correctionfactors: 			pmask[i] = true; 			break;
    //         case pythia_5TeV:                       pmask[i] = true;            break;
    //         case pythia_13TeV:                      pmask[i] = true;            break;
    //     }
    // }
    // maybe instead of pmask, make a dictionary that stores output directory, the files to run, and the marker styles?

    
    SetStyle();


    // analyze and plot!
    if (comp_case == 1) {
        plot_comparisons(raw_data_infile_hic, data_correctionfactors_infile_hic, pt_bins, n_bins, RL_bins, n_RLbins, include_RL0, include_RL1);
    } else if (comp_case == 2) {
        plot_comparisons(raw_data_infile_perly, binbybincorrections_infile_perly, pt_bins, n_bins, RL_bins, n_RLbins, include_RL0, include_RL1);
    } else if (comp_case == 3) {
        plot_comparisons(raw_data_infile_perly, , pt_bins, n_bins, RL_bins, n_RLbins, include_RL0, include_RL1);
    }
    
}