


// void cut_on_axis(THnD * hist, int axis, double minval, double maxval) {
//     hist->GetAxis(int)->SetRangeUser(minval, maxval);
// }
int colors[] = { kBlue, kRed, kGreen+2, kViolet };
int markers[] = { kFullCircle, kFullSquare, kFullDiamond, kFullStar, kFullCross };


// --- define bins ---
int pt_bins[] = { 20, 40, 60, 80, 100, 120, 150, 200 };
const int n_pt_bins = 7;
double kt_cuts[] = { 0, 0.5, 1.0, 10.0 };
int n_kt_cuts = 3;

std::string plot_filepath = "/software/users/blianggi/mypyjetty/storage/ktcuts/plots/";

struct FitRange {
    double min;
    double max;
    double sigma;
};

class EEC_Curve {
public:
    std::string name;
    std::string rl_choice; // "RL" or "ptRL"
    std::string kt_choice; // "kt" or "kappa"

    int color;
    int marker_style;
    int marker_size;

    int low_pt;
    int high_pt;
    double low_kt;
    double high_kt;
    double low_kappa;
    double high_kappa;

    // double lowx_fit;
    // double highx_fit;
    FitRange fit_range;

    double fit_param_mu;
    double fit_param_C;
    double fit_param_sg;
    
    TH1D * hist;
    TF1 * fitfunc;
    TH2D * hist2D; // this will hold the (i.e. kt vs RL) "lund plane"-like plot, but this is only stored for the EEC_Curve for the lowest pt bin !

    double peak_height_from_max;
    double peak_pos_from_max;

    EEC_Curve() {}
    EEC_Curve(std::string name_val, std::string rl_choice_val, std::string kt_choice_val, 
              int color_val, int marker_style_val) {
        name = name_val;
        rl_choice = rl_choice_val;
        kt_choice = kt_choice_val;

        color = color_val;
        marker_style = marker_style_val;
        marker_size = 1.0;
    }
};


FitRange loadFitRanges(std::string rl_add_name, std::string kt_add_name, int i, int j) {
    std::string ident_name = Form("%s_vs_%s", kt_add_name.c_str(), rl_add_name.c_str());
    
    FitRange fitrange;
    std::ifstream file("peak_ranges.csv");
    std::string line;

    // skip header - ignore the first 3 lines
    for ( int a = 0; a < 3; a++ ) std::getline(file, line);
    int a = 0;

    // format of line: curve_type, i, pt_range, min, max, min, max, min, max
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;

        std::string curve_type;
        int pt_bin; // aka "i"
        double min, max, sigma;

        std::getline(ss, token, ',');
        curve_type = token;

        std::getline(ss, token, ',');
        pt_bin = std::stoi(token);
        std::getline(ss, token, ','); // this is the pt range ... can be ignored

        if ( curve_type == ident_name && pt_bin == i ) {
            int b = 0;
            while (std::getline(ss, token, ',')) {
                min = std::stod(token);

                std::getline(ss, token, ',');
                max = std::stod(token);
                std::getline(ss, token, ',');
                sigma = std::stod(token);

                if ( b == j ) {
                    fitrange = {min, max, sigma};
                    return fitrange;
                }
                b++;
            }    
        }   
        a++;
    }

    return fitrange; // it shouldn't get to this point
}


double gaus_log(double *x, double *par) {
    double mu = par[0];
    double C  = par[1];
    double sg = par[2];

    if (x[0] <= 0) return 0.0; // avoid log(0) or negative x

    double val = C * exp( -pow(log(x[0]/mu), 2) / (2 * sg * sg) );
    return val;
}

void draw_fit_params(TF1 * f, int j, double ypos = 0.8) {
    double mu     = f->GetParameter(0);
    double mu_err = f->GetParError(0);
    double C      = f->GetParameter(1);
    double sg     = f->GetParameter(2);

    TLatex latex;
    latex.SetNDC();           // use normalized coordinates
    latex.SetTextSize(0.03);
    latex.SetTextColor(colors[j]);
    latex.DrawLatex(0.7, ypos, Form("Fit for k_{T} = %.1f - %.1f", kt_cuts[j], kt_cuts[j+1]));
    latex.SetTextColor(kBlack);
    latex.DrawLatex(0.7, ypos - 0.04, Form("#mu = %.3f #pm %.3f", mu, mu_err));
    latex.DrawLatex(0.7, ypos - 0.08, Form("C = %.3f", C));
    latex.DrawLatex(0.7, ypos - 0.12, Form("#sigma = %.3f", sg));
}


void move_stats_bar(TH2D * h2D) {
    TPaveStats *st = (TPaveStats*) h2D->GetListOfFunctions()->FindObject("stats");
    if (st) {
        st->SetX1NDC(0.2); // left x in NDC (0–1)
        st->SetX2NDC(0.4); // right x
        st->SetY1NDC(0.65); // bottom y
        st->SetY2NDC(0.85); // top y
    }
}

void set_axes_sizes(TH1D * hist, double x_title_size, double x_label_size, double y_title_size, double y_label_size) {
    hist->GetXaxis()->SetTitleSize(x_title_size); //0.10);
    hist->GetXaxis()->SetLabelSize(x_label_size); //0.08);
    hist->GetYaxis()->SetTitleSize(y_title_size); //0.09);
    hist->GetYaxis()->SetLabelSize(y_label_size); //0.08);
}


void format_1D_hist(TH1D * hist, int color, int marker, double markersize) {
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->SetMarkerStyle(marker);
    hist->SetMarkerSize(markersize);
}

void format_graph(TGraphErrors * graph, int color, int marker, double markersize, std::string title, std::string xtitle, std::string ytitle) {
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    graph->SetMarkerStyle(marker);
    graph->SetMarkerSize(markersize);

    graph->SetTitle(title.c_str());
    graph->GetXaxis()->SetTitle(xtitle.c_str());
    graph->GetYaxis()->SetTitle(ytitle.c_str());
}



// Get the histogram from THnSparse
EEC_Curve collect_hist(TFile *file, TString histname, std::string rl_add_name, std::string kt_add_name, int i, int j) {
    
    cout << "now running " << histname << " // " << rl_add_name << ", " << kt_add_name << endl;
    EEC_Curve eec_curve(Form("EEC_%s_%s", rl_add_name.c_str(), kt_add_name.c_str()), rl_add_name, kt_add_name, colors[j], markers[i]);

    THnSparse *hsparse = (THnSparse*) file->Get(histname);
    TH1D *h_jet_pt = (TH1D*) file->Get("h_jet_pt");

    // --- Get num jets ---
    h_jet_pt->GetXaxis()->SetRangeUser(pt_bins[i], pt_bins[i+1]);
    double num_jets = h_jet_pt->Integral();

    // --- Make pt cut ---
    hsparse->GetAxis(0)->SetRangeUser(pt_bins[i], pt_bins[i+1]);

    // --- Find and store 2D hist (unnormalized) ---
    if ( j == 0 ) {
        int xdim; int ydim;
        if ( rl_add_name == "RL" ) xdim = 1;
        else if ( rl_add_name == "ptRL" ) xdim = 2;
        if ( kt_add_name == "kt" ) ydim = 3;
        else if ( kt_add_name == "kappa" ) ydim = 4;

        // if (rl_add_name == "ptRL") hsparse->GetAxis(2)->SetRangeUser(0, 200);
        // if (kt_add_name == "kappa") hsparse->GetAxis(4)->SetRangeUser(0, 0.15);
        TH2D * hist2D_lund = (TH2D*) hsparse->Projection(ydim, xdim); //Project3D("zy");
        eec_curve.hist2D = hist2D_lund;
    }

    // --- Make more cuts, project, and normalize ---
    if (kt_add_name == "kt") hsparse->GetAxis(3)->SetRangeUser(kt_cuts[j], kt_cuts[j+1]);
    // else if (kt_add_name == "kappa") hsparse->GetAxis(4)->SetRangeUser(replace, replace);

    TH1D * hist_EEC_with_ktcut;
    if (rl_add_name == "RL") hist_EEC_with_ktcut = (TH1D*) hsparse->Projection(1);
    else if (rl_add_name == "ptRL") hist_EEC_with_ktcut = (TH1D*) hsparse->Projection(2);
    hist_EEC_with_ktcut->Scale(1/num_jets, "width");
    hist_EEC_with_ktcut->SetName(Form("hist_EEC_%s_%s_pt%d-%d_ktcut%.1f-%.1f", kt_add_name.c_str(), rl_add_name.c_str(), pt_bins[i], pt_bins[i+1], kt_cuts[j], kt_cuts[j+1]));
    eec_curve.hist = hist_EEC_with_ktcut;
    cout << "For pt=" << i << " number of entries = " << hist_EEC_with_ktcut->GetEntries() << endl;

    eec_curve.low_pt = pt_bins[i];
    eec_curve.high_pt = pt_bins[i+1];
    if ( kt_add_name == "kt" ) {
        eec_curve.low_kt = kt_cuts[j];
        eec_curve.high_kt = kt_cuts[j+1];
    }
    FitRange fr = loadFitRanges(rl_add_name, kt_add_name, i, j);
    eec_curve.fit_range = fr;
    // if ( rl_add_name == "RL" ) { // TODO: adjust this more????
    //     eec_curve.lowx_fit = 0.005;
    //     eec_curve.highx_fit = 1.0;
    // } else if ( rl_add_name == "ptRL" ) {
    //     eec_curve.lowx_fit = 0.01;
    //     eec_curve.highx_fit = 30.0;
    // }

    int temp_max_bin = eec_curve.hist->GetMaximumBin();
    double temp_max_pos = eec_curve.hist->GetBinCenter(temp_max_bin);
    double temp_max = eec_curve.hist->GetMaximum();
    eec_curve.fit_param_mu = temp_max_pos;
    eec_curve.fit_param_C = temp_max;
    eec_curve.fit_param_sg = fr.sigma; //0.5; // TODO: adjust this more????
    // if ( j == 0 ) eec_curve.fit_param_sg = 0.5; // TODO: adjust this more????
    // else if ( j == 1 ) eec_curve.fit_param_sg = 0.75; // TODO: adjust this more????
    // else if ( j == 2 ) eec_curve.fit_param_sg = 1.5; // TODO: adjust this more????

    
    return eec_curve;
    
}

void get_fit_for_eec_curve(EEC_Curve& eec_curve) {
    // fits - define TF1 using the external function
    // ranges: RL = 5E-3 - 1, pTRL = 2E-1 - 30
    TF1 *f = new TF1("f_gaus_log", gaus_log, eec_curve.fit_range.min, eec_curve.fit_range.max, 3); // name, function, xlow, xhigh, num_params
    f->SetParNames("mu", "C", "sigma");
    f->SetParameters(eec_curve.fit_param_mu, eec_curve.fit_param_C, eec_curve.fit_param_sg); // initial guesses

    eec_curve.hist->Fit(f, "NR");  // "R" = use TF1 range
    // Style fit
    f->SetLineColor(kOrange+7);
    f->SetLineWidth(2);

    eec_curve.fitfunc = f;
    eec_curve.fit_param_mu = f->GetParameter(0);
    eec_curve.fit_param_C = f->GetParameter(1);
    eec_curve.fit_param_sg = f->GetParameter(2);

    // f->Draw("same");
    // draw_fit_params(f, eec_curve);
}

void find_peak_pos_height_from_max(EEC_Curve& eec_curve) {
    double peak_height = eec_curve.hist->GetMaximum();
    eec_curve.peak_height_from_max = peak_height;

    int peak_pos_bin = eec_curve.hist->GetMaximumBin();
    double peak_pos = eec_curve.hist->GetBinCenter(peak_pos_bin);
    eec_curve.peak_pos_from_max = peak_pos;

    return;
}

// Find and set max for plotsmanship
void find_and_set_max(std::vector<std::vector<EEC_Curve>>& eec_curves) {
    for ( int i = 0; i < n_pt_bins; i++ ) {
        double temp_max = 0;
        for ( int j = 0; j < n_kt_cuts; j++ ) {
            temp_max = std::max(eec_curves[i][j].hist->GetMaximum(), temp_max);
        }
        // cout << "MAX! " << i << ": " << temp_max << endl;
        eec_curves[i][0].hist->SetMaximum(temp_max * 1.2);
    }
}

// Plot figures by pt bin. There should be as many of these as there are pt bins,
// for each 2D plot (lund-plane-like) and EEC with kt cuts (3 kt cut curves on each canvas).
void plot_ptbin_figs(std::vector<std::vector<EEC_Curve>> eec_curves, std::string rl_add_name, std::string kt_add_name) {
    // loop through pt bins to make plots
    for ( int i = 0; i < n_pt_bins; i++ ) {

        // Plot 2D plots (similar to Lund plane)
        TCanvas *can_2d = new TCanvas();
        gPad->SetLogx();
        // move_stats_bar(eec_curves[0].hist2D);

        eec_curves[i][0].hist2D->Draw("colz");
        can_2d->SaveAs(Form("%s/ptbins/%s_vs_%s_pt%d-%d.pdf", plot_filepath.c_str(), kt_add_name.c_str(), rl_add_name.c_str(), pt_bins[i], pt_bins[i+1]));

        // Plot EEC curves with kt cuts
        TCanvas *can_1d = new TCanvas();
        gStyle->SetOptStat(0);
        gPad->SetLogx();
        TLegend *leg_kt;
        if ( i < 2 && rl_add_name == "RL") leg_kt = new TLegend(0.15, 0.68, 0.35, 0.88);
        else leg_kt = new TLegend(0.45, 0.68, 0.65, 0.88);
        leg_kt->SetTextSize(0.035);
        leg_kt->AddEntry((TObject*)0, Form("p_{T, jet} = %d-%d", pt_bins[i], pt_bins[i+1]), "");
        for ( int j = 0; j < n_kt_cuts; j++ ) {
            // format_1D_hist(eec_curves[a-1][b].hist, eec_curves[a-1][b].color, eec_curves[a-1][b].marker_style, eec_curves[a-1][b].marker_size);
            format_1D_hist(eec_curves[i][j].hist, eec_curves[i][j].color, eec_curves[i][j].marker_style, eec_curves[i][j].marker_size);
            leg_kt->AddEntry(eec_curves[i][j].hist, Form("k_{T} = %.1f-%.1f", kt_cuts[j], kt_cuts[j+1]));
            eec_curves[i][j].hist->Draw("same");
            eec_curves[i][j].fitfunc->Draw("same");
            draw_fit_params(eec_curves[i][j].fitfunc, j, 0.85 - (0.2*j)); // using equation 0.7 - 0.2j to determine position of fit param text
        }
        leg_kt->Draw();
        can_1d->SaveAs(Form("%s/ptbins/EEC_%s_%s_pt%d-%d_ktcuts.pdf", plot_filepath.c_str(), kt_add_name.c_str(), rl_add_name.c_str(), pt_bins[i], pt_bins[i+1]));

    }
}

// Plot 1D plots of all projected EECs and ratios. These will be in one plot with 8 panels.
// Each panel will have the top be the three curves (1 curve for each kt cut), 
// and the bottom will be the ratio to the lowest kT cut.
void plot_8_panel_proj_EEC(std::vector<std::vector<EEC_Curve>> eec_curves, std::string rl_add_name, std::string kt_add_name) {

    TCanvas *can_8panel = new TCanvas("can_8panel", "Ratio plots", 800,600);
    can_8panel->Divide(4,2, 0.0, 0.0);

    // Make legend with the kt info
    TLegend * leg_kt = new TLegend(0.15, 0.3, 0.85, 0.7);
    leg_kt->SetTextSize(0.07);

    for ( int a = 1; a <= n_pt_bins; a++ ) {
        can_8panel->cd(a);

        gPad->SetLeftMargin(0.08);
        gPad->SetRightMargin(0.02);
        gPad->SetTopMargin(0.02);
        gPad->SetBottomMargin(0.08);

        // Get panel dimensions
        double left   = gPad->GetLeftMargin();
        double right  = 1 - gPad->GetRightMargin();
        double bottom = gPad->GetBottomMargin();
        double top    = 1 - gPad->GetTopMargin();

        // Create top and bottom pads inside this panel
        TPad *pad_top = new TPad(Form("pad_top_%d",a),"",0,0.3,1,1.0);
        TPad *pad_bot = new TPad(Form("pad_bot_%d",a),"",0,0.0,1,0.3);

        pad_top->SetBottomMargin(0.0); //2);  // reduce gap
        pad_bot->SetTopMargin(0.0); //5);
        pad_bot->SetBottomMargin(0.25); // leave space for x-axis labels

        // pad_top->SetTopMargin(0.02);
        // pad_top->SetBottomMargin(0.01);  // almost no gap
        pad_top->SetLeftMargin(0.15);
        pad_top->SetRightMargin(0.02);

        // pad_bot->SetTopMargin(0.01);     // shrink gap between top/bottom
        // pad_bot->SetBottomMargin(0.3);   // keep enough for x-axis labels
        pad_bot->SetLeftMargin(0.15);
        pad_bot->SetRightMargin(0.02);

        pad_top->Draw();
        pad_bot->Draw();

        std::vector<TH1D *> hratio_vec;
        double ratio_max = 0; double ratio_min = 1000;

        // plot top panel
        pad_top->cd();
        gPad->SetLogx();
        for ( int b = 0; b < n_kt_cuts; b++ ){
            set_axes_sizes(eec_curves[a-1][b].hist, 0.10, 0.08, 0.07, 0.05);
            eec_curves[a-1][b].hist->SetTitle("");
            eec_curves[a-1][b].hist->Draw("same");
            eec_curves[a-1][b].hist->SetMarkerSize(0.5);

            if ( a == 1 ) leg_kt->AddEntry(eec_curves[a-1][b].hist, Form("k_{T} = %.1f - %.1f", kt_cuts[b], kt_cuts[b+1]));

            if (b != 0) {
                TH1D * hratio = (TH1D*) eec_curves[a-1][0].hist->Clone(Form("h_ratio_%d",b));
                hratio->Divide(eec_curves[a-1][b].hist);
                format_1D_hist(hratio, eec_curves[a-1][b].color, eec_curves[a-1][b].marker_style, 0.5);
                hratio_vec.push_back(hratio);

                double temp_max = hratio->GetBinContent(hratio->GetMaximumBin());
                double temp_min = hratio->GetBinContent(hratio->GetMinimumBin());
                ratio_max = (temp_max > ratio_max) ? temp_max : ratio_max;
                ratio_min = (temp_min < ratio_min) ? temp_min : ratio_min;
            }
        }
        TLatex latex;
        latex.SetNDC();                  // use normalized coordinates in the pad
        latex.SetTextSize(0.05);
        latex.DrawLatex(0.57, 0.75, Form("%d < p_{T} < %d", pt_bins[a-1], pt_bins[a])); // 0.65, 0.7, 0.88, 0.8
        
        // plot bottom panel
        pad_bot->cd();
        gPad->SetLogx();
        gPad->SetLogy();
        cout << "DID MAX HERE " << ratio_max << " " << ratio_min << endl;
        hratio_vec[0]->SetMaximum(ratio_max);
        // hratio_vec[0]->SetMinimum(ratio_min);
        for (int b = 0; b < n_kt_cuts-1; b++){ 
            hratio_vec[b]->SetTitle("");
            hratio_vec[b]->GetYaxis()->SetTitle("k_{T}=0-0.5 / other");
            set_axes_sizes(hratio_vec[b], 0.10, 0.09, 0.09, 0.09);
            hratio_vec[b]->Draw("same");
        }
    }

    can_8panel->cd(8);
    // add legend and stuff here
    leg_kt->Draw();

    // Save!
    can_8panel->SaveAs(Form("%spanels_with_ratios_%s_%s.pdf",plot_filepath.c_str(), kt_add_name.c_str(), rl_add_name.c_str()));

}
   
// Plot the peak position and heights into one summary plot as a function of pT.
// Open circles represent just finding the maximum bin.
// Closed circles are values taken from a gaus-log fit.
void plot_peaks(std::vector<std::vector<EEC_Curve>> eec_curves, std::string rl_add_name, std::string kt_add_name) {
    TCanvas *can_peaks = new TCanvas();
    can_peaks->Divide(2);

    // TCanvas *can_peak_heights = new TCanvas();
    // TCanvas *can_peak_pos = new TCanvas();

    double pt_bins_avg[n_pt_bins] = {};
    for (int i = 0; i < n_pt_bins; i++) pt_bins_avg[i] = (pt_bins[i] + pt_bins[i+1]) / 2;

    std::vector<TGraphErrors* > graph_peak_pos;
    std::vector<TGraphErrors* > graph_peak_heights;
    std::vector<TGraphErrors* > graph_peak_pos_fit;
    std::vector<TGraphErrors* > graph_peak_heights_fit;
    double max_pos = 0;
    double max_height = 0;
    
    // --- get and plot peaks of each curve as a function of pT ---
    for (int j = 0; j < n_kt_cuts; j++) {

        std::vector<double> arr_peak_pos;
        std::vector<double> arr_peak_heights;
        std::vector<double> arr_peak_pos_fit;
        std::vector<double> arr_peak_heights_fit;

        std::vector<double> arr_pt_errs;
        std::vector<double> arr_peak_pos_errs;
        std::vector<double> arr_peak_heights_errs;
        std::vector<double> arr_peak_pos_fit_errs;
        std::vector<double> arr_peak_heights_fit_errs;
        for (int i = 0; i < n_pt_bins; i++) {
            cout << "SABING" << eec_curves[i][j].peak_pos_from_max << " AND " << eec_curves[i][j].peak_height_from_max << endl;
            arr_peak_pos.push_back(eec_curves[i][j].peak_pos_from_max);
            arr_peak_heights.push_back(eec_curves[i][j].peak_height_from_max);
            arr_peak_pos_fit.push_back(eec_curves[i][j].fit_param_mu);
            arr_peak_heights_fit.push_back(eec_curves[i][j].fit_param_C);

            arr_pt_errs.push_back(0);
            arr_peak_pos_errs.push_back(0);
            arr_peak_heights_errs.push_back(0);
            arr_peak_pos_fit_errs.push_back(eec_curves[i][j].fitfunc->GetParError(0)); // mu error
            arr_peak_heights_fit_errs.push_back(eec_curves[i][j].fitfunc->GetParError(1)); // C error
        }

        TGraphErrors *gr_peak_pos_temp = new TGraphErrors(n_pt_bins, pt_bins_avg, arr_peak_pos.data(), arr_pt_errs.data(), arr_peak_pos_errs.data());
        format_graph(gr_peak_pos_temp, colors[j], kOpenCircle, 1.0, "Peak positions", "Average p_{T} (GeV/c)", "Normalized EEC peak position");
        graph_peak_pos.push_back(gr_peak_pos_temp);

        TGraphErrors *gr_peak_pos_fit_temp = new TGraphErrors(n_pt_bins, pt_bins_avg, arr_peak_pos_fit.data(), arr_pt_errs.data(), arr_peak_pos_fit_errs.data());
        format_graph(gr_peak_pos_fit_temp, colors[j], kFullCircle, 1.0, "Peak positions", "Average p_{T} (GeV/c)", "Normalized EEC peak position");
        graph_peak_pos_fit.push_back(gr_peak_pos_fit_temp);

        double max_temp = *std::max_element(arr_peak_pos.begin(), arr_peak_pos.end());
        max_pos = ( max_temp > max_pos ) ? max_temp : max_pos;
        double max_fit_temp = *std::max_element(arr_peak_pos_fit.begin(), arr_peak_pos_fit.end());
        max_pos = ( max_fit_temp > max_pos ) ? max_fit_temp : max_pos;
        
        TGraphErrors *gr_peak_heights_temp = new TGraphErrors(n_pt_bins, pt_bins_avg, arr_peak_heights.data(), arr_pt_errs.data(), arr_peak_heights_errs.data());
        format_graph(gr_peak_heights_temp, colors[j], kOpenCircle, 1.0, "Peak heights", "Average p_{T} (GeV/c)", "Normalized EEC peak height");
        graph_peak_heights.push_back(gr_peak_heights_temp);

        TGraphErrors *gr_peak_heights_fit_temp = new TGraphErrors(n_pt_bins, pt_bins_avg, arr_peak_heights_fit.data(), arr_pt_errs.data(), arr_peak_heights_fit_errs.data());
        format_graph(gr_peak_heights_fit_temp, colors[j], kFullCircle, 1.0, "Peak heights", "Average p_{T} (GeV/c)", "Normalized EEC peak height");
        graph_peak_heights_fit.push_back(gr_peak_heights_fit_temp);

        max_temp = *std::max_element(arr_peak_heights.begin(), arr_peak_heights.end());
        max_height = ( max_temp > max_height ) ? max_temp : max_height;
        max_fit_temp = *std::max_element(arr_peak_heights_fit.begin(), arr_peak_heights_fit.end());
        max_height = ( max_fit_temp > max_height ) ? max_fit_temp : max_height;
    }

    // Set necessary min/maxes
    graph_peak_heights[0]->SetMinimum(0.0);
    graph_peak_pos[0]->SetMaximum(max_pos * 1.1);
    graph_peak_heights[0]->SetMaximum(max_height * 1.1);
    cout << "MAX PEAJK POS AND HEIGHT " << max_pos << " AND " << max_height << endl;
    // if (rl_add_name == "RL") gr_peak_pos_temp->SetMaximum(0.4);
    // if (rl_add_name == "ptRL" || kt_add_name == "kappa") {
    //     gr_peak_pos_temp->SetMaximum(100);
    //     gPad->SetLogy();
    // }
    
    // Now loop through again to plot
    for (int j = 0; j < n_kt_cuts; j++) {
        can_peaks->cd(1);
        if ( j == 0 ) graph_peak_pos[j]->Draw("AP");
        else graph_peak_pos[j]->Draw("P SAME");
        graph_peak_pos_fit[j]->Draw("P SAME");

        can_peaks->cd(2);
        if ( j == 0 ) graph_peak_heights[j]->Draw("AP");
        else graph_peak_heights[j]->Draw("P SAME");  
        graph_peak_heights_fit[j]->Draw("P SAME");

        double x, y;
        for ( int k = 0; k < n_pt_bins; k++) {
            graph_peak_pos[j]->GetPoint(k, x, y);
            std::cout << "POS Point " << k << ": x = " << x << ", y = " << y << std::endl;
            graph_peak_heights[j]->GetPoint(k, x, y);
            std::cout << "HEIGHT Point " << k << ": x = " << x << ", y = " << y << std::endl;
        }
    }

    can_peaks->cd(1);
    gPad->SetLeftMargin(0.15);
    // leg_peaks->Draw();

    can_peaks->cd(2);
    gPad->SetLeftMargin(0.15);
    // leg_peaks->Draw();
    can_peaks->SaveAs(Form("%sAll_EEC_%s_%s_peaks.pdf", plot_filepath.c_str(), kt_add_name.c_str(), rl_add_name.c_str()));


}



// Analyze one configuration at a time
// i.e. kt vs RL
void analyze_one_config(TFile * file, TString histname, std::string rl_add_name, std::string kt_add_name) {
    
    std::vector<std::vector<EEC_Curve>> all_eec_curves;

    for ( int i = 0; i < n_pt_bins; i++ ) {

        std::vector<EEC_Curve> eec_curves_per_kt;

        for ( int j = 0; j < n_kt_cuts; j++ ) {
            
            EEC_Curve eec_curve_temp = collect_hist(file, histname, rl_add_name, kt_add_name, i, j);
            eec_curves_per_kt.push_back(eec_curve_temp);
            get_fit_for_eec_curve(eec_curves_per_kt[j]);
            find_peak_pos_height_from_max(eec_curves_per_kt[j]);

        }

        all_eec_curves.push_back(eec_curves_per_kt);
    }

    size_t n_rows = all_eec_curves.size();                     // Number of rows
    size_t n_cols = all_eec_curves.empty() ? 0 : all_eec_curves[0].size();  // Number of columns (from first r
    std::cout << "Rows: " << n_rows << ", Columns: " << n_cols << std::endl;
    // Find and set max for plotsmanship
    find_and_set_max(all_eec_curves);
    
    // make 2D plots
    plot_ptbin_figs(all_eec_curves, rl_add_name, kt_add_name);

    // make 8-panel of projected EEC + ratio to blue
    plot_8_panel_proj_EEC(all_eec_curves, rl_add_name, kt_add_name);

    // put peak height and peak pos on 2-panel plot
    plot_peaks(all_eec_curves, rl_add_name, kt_add_name);

}

void analyze_dif_cuts() {
    // --- Open the ROOT file ---
    TString filename = "/software/users/blianggi/mypyjetty/storage/ktcuts/rootfiles/pythia_jse_output_ptmin20.root";
    TFile *file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return;
    }

    
    // --- Get the histograms ---
    // TString histname = "h3D_pt_vs_RL_vs_kT"; //3D - pT, RL, kT
    // TString histname_ptrl = "h3D_pt_vs_ptRL_vs_kT";
    // TString histname_kappa = "h3D_pt_vs_ptRL_vs_kappa";

    TString histname = "hND_all_pair_info"; //5D - pT, RL, <pT>RL, kT, kappa, weighted
    

    // analyze_hist(file, histname, "RL", "kt");
    // analyze_hist(file, histname_ptrl, "ptRL", "kt");
    // analyze_hist(file, histname_kappa, "ptRL", "kappa");

    // analyze_hist(file, histname_kappa, "ptRL", "kappa-ktcut0-1");
    // analyze_hist(file, histname_kappa, "ptRL", "kappa-ktcut1-10");

    // kt vs RL
    analyze_one_config(file, histname, "RL", "kt");

    // kt vs ptRL
    analyze_one_config(file, histname, "ptRL", "kt");

    
    
    
}