


// void cut_on_axis(THnD * hist, int axis, double minval, double maxval) {
//     hist->GetAxis(int)->SetRangeUser(minval, maxval);
// }
int colors[] = { kBlue, kRed, kGreen+2, kViolet };
int markers[] = { kFullCircle, kFullSquare, kFullDiamond, kFullStar, kFullCross };

std::string plot_filepath = "/software/users/blianggi/mypyjetty/storage/ktcuts/plots/";

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
    double low_rl_choice;
    double high_rl_choice;
    double low_kt_cut;
    double high_kt_cut;

    double fit_param_mu;
    double fit_param_C;
    double fit_param_sg;
    
    TH1D * hist;
    TF1 * fitfunc;

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


double gaus_log(double *x, double *par) {
    double mu = par[0];
    double C  = par[1];
    double sg = par[2];

    if (x[0] <= 0) return 0.0; // avoid log(0) or negative x

    double val = C * exp( -pow(log(x[0]/mu), 2) / (2 * sg * sg) );
    return val;
}

void draw_fit_params(TF1 * f, EEC_Curve eec_curve) {
    double mu     = f->GetParameter(0);
    double mu_err = f->GetParError(0);
    double C      = f->GetParameter(1);
    double sg     = f->GetParameter(2);

    eec_curve.fit_param_mu = mu;
    eec_curve.fit_param_C = C;
    eec_curve.fit_param_sg = sg;

    TLatex latex;
    latex.SetNDC();           // use normalized coordinates
    latex.SetTextSize(0.03);
    latex.DrawLatex(0.65, 0.6, Form("#mu = %.3f #pm %.3f", mu, mu_err));
    latex.DrawLatex(0.65, 0.55, Form("C = %.3f", C));
    latex.DrawLatex(0.65, 0.50, Form("#sigma = %.3f", sg));
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


void cut_on_th3d(TH3D * hist3D, double xmin, double xmax, double ymin, double ymax) {
    hist3D->GetXaxis()->SetRangeUser(xmin, xmax);
    hist3D->GetYaxis()->SetRangeUser(ymin, ymax);
}

void collect_all_hists(TFile *file, TString histname, std::string rl_add_name, std::string kt_add_name) {
    cout << "now running " << histname << " // " << rl_add_name << ", " << kt_add_name << endl;

    TH3D *hist3D = (TH3D*) file->Get(histname);
    TH1D *h_jet_pt = (TH1D*) file->Get("h_jet_pt");

    
}


void analyze_hist(TFile *file, TString histname, std::string rl_add_name, std::string kt_add_name) {

    // --- Make canvas ---
    TCanvas *can_all_EECs = new TCanvas();
    gPad->SetLogx();

    TLegend * leg = new TLegend();
    TLegend *leg_all_EECs = new TLegend(0.5,0.45,0.8,0.85);
    TLegend * leg_peaks = new TLegend(0.37,0.7,0.62,0.88);
    leg_peaks->SetBorderSize(1);
    leg_peaks->SetFillStyle(0); // transparent

    TCanvas *can_peak_heights = new TCanvas();
    TCanvas *can_peak_pos = new TCanvas();

    
    // --- define bins, loop through ---
    int pt_bins[] = { 20, 40, 60, 80, 100, 120, 150, 200 };
    const int n_pt_bins = 7;
    double kt_cuts[] = { 0, 0.5, 1.0, 10.0 };
    int n_kt_cuts = 3;
    double max_y_val = 0;

    // std::vector<std::vector<TH1D*>> vec_hist_EEC_with_ktcut;
    std::vector<std::vector<EEC_Curve>> all_eec_curves;

    double pt_bins_avg[n_pt_bins] = {};

    for (int i = 0; i < n_pt_bins; i++) {

        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];
        pt_bins_avg[i] = (pt_bins[i] + pt_bins[i+1]) / 2;

        std::vector<EEC_Curve> eec_curves_per_pt;

        // --- Get num jets ---
        TH1D * h_jet_pt_clone = (TH1D*) h_jet_pt->Clone("h_jet_pt_clone");
        h_jet_pt_clone->GetXaxis()->SetRangeUser(pt_min, pt_max);
        double num_jets = h_jet_pt_clone->Integral();

        // --- Plot kT vs RL ---
        TH3D * hist3D_clone = (TH3D*) hist3D->Clone("hist3D_clone");
        hist3D_clone->GetXaxis()->SetRangeUser(pt_min, pt_max);
        if (rl_add_name == "ptRL") hist3D_clone->GetYaxis()->SetRangeUser(0, 200);
        if (kt_add_name == "kappa") hist3D_clone->GetZaxis()->SetRangeUser(0, 0.15);
        TH2D * hist2D_kt_vs_RL = (TH2D*) hist3D_clone->Project3D("zy"); 
        
        TCanvas * can_2d = new TCanvas();
        can_2d->cd();
        gPad->SetLogx();
        move_stats_bar(hist2D_kt_vs_RL);
        hist2D_kt_vs_RL->Draw("colz");
        can_2d->SaveAs(Form("%s%s_vs_%s_pt%d-%d.pdf", plot_filepath.c_str(), kt_add_name.c_str(), rl_add_name.c_str(), pt_min, pt_max));

        // --- Define kT-axis cut range - now on the Y axis ---
        std::vector<TH1D*> vec_hist_EEC_with_ktcut_temp;
        TCanvas * can_EEC_with_ktcut = new TCanvas();
        gStyle->SetOptStat(0);
        gPad->SetLogx();
        for (int j = 0; j < n_kt_cuts; j++) {
            double kt_min = kt_cuts[j];
            double kt_max = kt_cuts[j+1];
            
            EEC_Curve eec_curve("", rl_add_name, kt_add_name, colors[j], markers[i]);
            eec_curve.low_pt = pt_min;
            eec_curve.high_pt = pt_max;
            eec_curve.low_kt_cut = kt_min;
            eec_curve.high_kt_cut = kt_max;

            TH2D * hist2D_kt_vs_RL_clone = (TH2D*) hist2D_kt_vs_RL->Clone("hist2D_kt_vs_RL_clone");
            hist2D_kt_vs_RL_clone->GetYaxis()->SetRangeUser(kt_min, kt_max);
            TH1D * hist_EEC_with_ktcut = (TH1D*) hist2D_kt_vs_RL_clone->ProjectionX();

            hist_EEC_with_ktcut->SetName(Form("hist_EEC_%s_%s_pt%d-%d_ktcut%.1f-%.1f", kt_add_name.c_str(), rl_add_name.c_str(), pt_min, pt_max, kt_min, kt_max));
            cout << "pushing back before " << i << " // " << j << endl;
            // vec_hist_EEC_with_ktcut_temp.push_back(hist_EEC_with_ktcut);
            // cout << "pushing back after " << i << " // " << j << endl;
            hist_EEC_with_ktcut->SetLineColor(colors[j]);
            hist_EEC_with_ktcut->SetMarkerColor(colors[j]);
            hist_EEC_with_ktcut->SetMarkerStyle(markers[i]);
            hist_EEC_with_ktcut->SetMarkerSize(1.0);
            hist_EEC_with_ktcut->Scale(1/num_jets, "width");
            if ( i == 0 ) {
                leg->AddEntry(hist_EEC_with_ktcut, Form("k_{T} = %.1f-%.1f", kt_min, kt_max));
                leg_peaks->AddEntry(hist_EEC_with_ktcut, Form("k_{T} = %.1f-%.1f", kt_min, kt_max));
                leg_all_EECs->AddEntry(hist_EEC_with_ktcut, Form("k_{T} = %.1f-%.1f", kt_min, kt_max));
            }
            max_y_val = (hist_EEC_with_ktcut->GetMaximum() > max_y_val) ? hist_EEC_with_ktcut->GetMaximum() : max_y_val;

            can_EEC_with_ktcut->cd();
            eec_curve.hist = hist_EEC_with_ktcut;
            hist_EEC_with_ktcut->Draw("SAME");   
            
            // fits - define TF1 using the external function
            // ranges: RL = 5E-3 - 1, pTRL = 2E-1 - 30
            TF1 *f;
            if (rl_add_name == "RL") f = new TF1("f_gaus_log", gaus_log, 0.005, 1, 3);
            else f = new TF1("f_gaus_log", gaus_log, 0.01, 30, 3);
            f->SetParNames("mu", "C", "sigma");
            int temp_max_bin = hist_EEC_with_ktcut->GetMaximumBin();
            double temp_max = hist_EEC_with_ktcut->GetMaximum();
            double temp_max_pos = hist_EEC_with_ktcut->GetBinCenter(temp_max_bin);
            f->SetParameters(temp_max_pos, temp_max, 0.5); // initial guesses

            hist_EEC_with_ktcut->Fit(f, "R");  // "R" = use TF1 range
            // Style fit
            f->SetLineColor(kOrange);
            f->SetLineWidth(2);
            eec_curve.fitfunc = f;
            f->Draw("same");
            draw_fit_params(f, eec_curve);

            eec_curves_per_pt.push_back(eec_curve);
            

        }
        
        all_eec_curves.push_back(eec_curves_per_pt);
        // vec_hist_EEC_with_ktcut.push_back(vec_hist_EEC_with_ktcut_temp);

        can_EEC_with_ktcut->cd();
        leg->Draw();
        can_EEC_with_ktcut->SaveAs(Form("%sEEC_%s_%s_pt%d-%d_ktcuts.pdf", plot_filepath.c_str(), kt_add_name.c_str(), rl_add_name.c_str(), pt_min, pt_max));

        // leg_all_EECs->AddEntry(NULL,"");
        leg_all_EECs->AddEntry(all_eec_curves[i][0].hist, Form("p_{T,jet} = %d-%d", pt_min, pt_max));
        
    }

    // --- plot all EECs and ratios ---
    TCanvas * can_panels = new TCanvas("can_panels", "Ratio plots", 900,600);
    can_panels->Divide(3,2);
    for (int a=1; a<=6; a++) {
        can_panels->cd(a);

        // Define legend
        TLegend * leg_pt = new TLegend(0.65, 0.7, 0.88, 0.8);
        leg_pt->SetTextSize(0.04);
        
        // get panel dimensions
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
        pad_bot->SetBottomMargin(0.35); // leave space for x-axis labels

        pad_top->Draw();
        pad_bot->Draw();

        // plot top panels
        pad_top->cd();
        gPad->SetLogx();
        for (int b = 0; b < n_kt_cuts; b++){
            set_axes_sizes(all_eec_curves[a-1][b].hist, 0.10, 0.08, 0.07, 0.05);
            all_eec_curves[a-1][b].hist->SetTitle("");
            all_eec_curves[a-1][b].hist->Draw("same");
            all_eec_curves[a-1][b].hist->SetMarkerSize(0.5);
        }
        leg_pt->AddEntry((TObject*)0, Form("%d < p_{T} < %d", pt_bins[a-1], pt_bins[a]), "");
        leg_pt->Draw();

        // plot top panels
        pad_bot->cd();
        gPad->SetLogx();
        gPad->SetLogy();
        for (int b = 1; b < n_kt_cuts; b++){
            TH1D * hratio = (TH1D*) all_eec_curves[a-1][0].hist->Clone(Form("h_ratio_%d",b));
            hratio->Divide(all_eec_curves[a-1][b].hist);

            hratio->SetTitle("");
            set_axes_sizes(hratio, 0.10, 0.09, 0.09, 0.09);
            hratio->SetMarkerColor(colors[b]);
            hratio->SetLineColor(colors[b]);
            hratio->Draw("same");
        }

    }

    can_panels->SaveAs(Form("%spanels_with_ratios_%s_%s.pdf",plot_filepath.c_str(), kt_add_name.c_str(), rl_add_name.c_str()));
    
    

    // --- plot peaks of each curve as a function of pT ---
    // --- get peaks of each curve as a function of pT ---
    // std::vector<std::vector<double>> arr_peaks;
    // std::vector<TGraph*> gr_peaks;
    for (int j = 0; j < n_kt_cuts; j++) {
        std::vector<double> arr_peak_heights_temp;
        std::vector<double> arr_peak_pos_temp;
        for (int i = 0; i < n_pt_bins; i++) {
            arr_peak_heights_temp.push_back(all_eec_curves[i][j].hist->GetMaximum());
            // cout << "max: " << i << "," << j << ": " << vec_hist_EEC_with_ktcut[i][j]->GetMaximum() << endl;

            int bin_max = all_eec_curves[i][j].hist->GetMaximumBin();
            double peak_pos = all_eec_curves[i][j].hist->GetBinCenter(bin_max);
            arr_peak_pos_temp.push_back(peak_pos);
            cout << "max: " << i << "," << j << ": " << peak_pos << endl;
        }
        // arr_peaks.push_back(arr_peak_heights_temp);
        TGraph *gr_peak_heights_temp = new TGraph(n_pt_bins, pt_bins_avg, arr_peak_heights_temp.data());
        gr_peak_heights_temp->SetMarkerColor(colors[j]);
        gr_peak_heights_temp->SetMarkerStyle(kFullCircle);
        gr_peak_heights_temp->SetTitle("Peak heights");
        gr_peak_heights_temp->GetXaxis()->SetTitle("Average p_{T}");
        gr_peak_heights_temp->GetYaxis()->SetTitle("Normalized EEC peak height");
        // gr_peaks.push_back(gr_peak_heights_temp);

        
        TGraph *gr_peak_pos_temp = new TGraph(n_pt_bins, pt_bins_avg, arr_peak_pos_temp.data());
        gr_peak_pos_temp->SetMarkerColor(colors[j]);
        gr_peak_pos_temp->SetMarkerStyle(kFullCircle);
        gr_peak_pos_temp->SetTitle("Peak positions");
        gr_peak_pos_temp->GetXaxis()->SetTitle("Average p_{T}");
        gr_peak_pos_temp->GetYaxis()->SetTitle("Normalized EEC peak position");


        can_peak_heights->cd();
        if ( j == 0 ) {
            gr_peak_heights_temp->SetMinimum(0.0);
            gr_peak_heights_temp->Draw("AP");
        }
        else gr_peak_heights_temp->Draw("P SAME");

        can_peak_pos->cd();
        if ( j == 0 ) {
            if (rl_add_name == "RL") gr_peak_pos_temp->SetMaximum(0.4);
            if (rl_add_name == "ptRL" || kt_add_name == "kappa") {
                gr_peak_pos_temp->SetMaximum(100);
                gPad->SetLogy();
            }
            gr_peak_pos_temp->Draw("AP");
        }
        else gr_peak_pos_temp->Draw("P SAME");
    }
    can_peak_heights->cd();
    leg_peaks->Draw();
    can_peak_heights->SaveAs(Form("%sAll_EEC_%s_%s_peak_heights.pdf", plot_filepath.c_str(), kt_add_name.c_str(), rl_add_name.c_str()));

    can_peak_pos->cd();
    leg_peaks->Draw();
    can_peak_pos->SaveAs(Form("%sAll_EEC_%s_%s_peak_pos.pdf", plot_filepath.c_str(), kt_add_name.c_str(), rl_add_name.c_str()));


    // --- plot all together --- and // --- plot peaks of each curve as a function of pT ---
    all_eec_curves[0][0].hist->SetMaximum(max_y_val * 1.2); // need to do this after extracting maxima
    can_all_EECs->cd();
    for (int i = 0; i < n_pt_bins; i++) {
        for (int j = 0; j < n_kt_cuts; j++) {
            all_eec_curves[i][j].hist->Draw("SAME");
        }
    }
    leg_all_EECs->Draw();
    can_all_EECs->SaveAs(Form("%sAll_EEC_%s_%s_ktcuts.pdf", plot_filepath.c_str(), kt_add_name.c_str(), rl_add_name.c_str()));
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
    TString histname = "h3D_pt_vs_RL_vs_kT"; //3D - pT, RL, kT
    TString histname_ptrl = "h3D_pt_vs_ptRL_vs_kT";
    TString histname_kappa = "h3D_pt_vs_ptRL_vs_kappa";
    

    analyze_hist(file, histname, "RL", "kt");
    analyze_hist(file, histname_ptrl, "ptRL", "kt");
    analyze_hist(file, histname_kappa, "ptRL", "kappa");

    analyze_hist(file, histname_kappa, "ptRL", "kappa-ktcut0-1");
    analyze_hist(file, histname_kappa, "ptRL", "kappa-ktcut1-10");


    collect_all_hists();
    collect_all_fits();

    
    
    
}