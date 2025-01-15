// A file that looks plots the response matrices of mc productions.



std::string outputdir = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/plots/RMs/LHC23a3";


//===========================================================================
//================================= FITTING =================================
//===========================================================================
double fit_histogram_guassianfit(TH1D * hist, double lowx, double highx) {

    /*
    // could also do:
    TF1 *linearFit = new TF1("linearFit", "pol1", 2.0, 8.0); // Fit range: [2, 8]
    hist->Fit(linearFit, "R"); // "R" = ensures the fit uses the specified range
    */
    TF1 *f1 = new TF1("f1","gausn", lowx, highx);
    // TF1 *fitfunc = new TF1("mstotal","gausn(0)", lowx, highx);

    // Perform the linear fit
    auto fitResult = hist->Fit("f1", "RS"); //, "Q"); // "pol1" = linear function, "Q" = quiet mode
    // "W": Ignore weights - set the weights of all non-zero bins to 1
    // "E": Perform better error estimation.
    // "S": The full result of the fit is returned in the TFitResultPtr (incl cov matrix)
    // 7.1.1. in https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html
    
    // Retrieve fit parameters
    TF1 *fitFunction = hist->GetFunction("f1");
    if (fitFunction) { // this is to catch the scenarios where fitfunction is empty (couldn't fit)
        double height = fitFunction->GetParameter(0);
        double mean = fitFunction->GetParameter(1);
        double sigma = fitFunction->GetParameter(2); 
        std::cout << "Fit Results: height = " << height << ", mean = " << mean  << ", sigma = " << sigma << std::endl;
    
        // auto fitResult = hist->Fit("pol1", "S");
        double chi2 = fitResult->Chi2();
        int ndf = fitResult->Ndf();
        std::cout << "Chi2/Ndf = " << chi2 / ndf << std::endl;

        return sigma;

    } else {
        return -1;
    }

    // // Draw the histogram and the fit
    // TCanvas *c1 = new TCanvas("c1", "Linear Fit", 800, 600);
    // TFile *dummy_file;
    // plot_and_save_one_histogram(c1, dummy_file, observable, hist, 3, addname);

}

//===========================================================================
//================================ FUNCTIONS ================================
//===========================================================================

void move_stat_box(TCanvas *c1, TH2D* hist2D) {
    // ge = (TGraphErrors *)c1->GetListOfPrimitives()->FindObject("Graph");
    // ge->Fit("gaus");
    // c1->Modified();
    // c1->Update();
    // TPaveText * pt = (TPaveText *)c1->GetListOfPrimitives()->FindObject("title");
    // pt->SetX1NDC(0.4);
    // pt->SetX2NDC(0.7);
    // pt->SetY1NDC(0.15);
    // pt->SetY2NDC(0.25);
    TPaveStats * ps = (TPaveStats *)hist2D->GetListOfFunctions()->FindObject("stats");
    ps->SetX1NDC(0.15);
    ps->SetX2NDC(0.55);
    c1->Modified();
    c1->Update();
}


void extract_jetpt_RM(TFile *fin) {
    std::string original_histname = Form("hResponse_JetPt_jet_pt_R0.4_1.0Scaled");
    fin->cd();
    TH2D * hist2D = (TH2D *)gDirectory->Get(original_histname.c_str());

    // save 
    TCanvas *can = new TCanvas();
    gPad->SetLogz();
    hist2D->GetZaxis()->SetRangeUser(1e-6,1e4);
    hist2D->Draw("colz");

    std::string output_name = Form("%s/RM_jetpt.pdf", outputdir.c_str());
    can->SaveAs(output_name.c_str());

}

void extract_LHC23a3_plots(TFile *fin, std::string observable) {
        
    bool logz_bool = true;
    std::string logstring = "_log"; 
    if (!logz_bool) {
        logstring = "";
    }

    fin->cd();
    std::string original_2d_histname = Form("h2D_corr_%s_JetPt_PT2050_R0.4_1.0Scaled", observable.c_str());
    TH2D * hist_2d = (TH2D *)gDirectory->Get(original_2d_histname.c_str());

    hist_2d->GetZaxis()->SetRangeUser(1e-6,1e4);

    // save 
    TCanvas *can = new TCanvas();
    if (logz_bool) gPad->SetLogz();
    hist_2d->Draw("colz same");
    // move_stat_box(can, hist_2d);
    // TPaveStats * ps = (TPaveStats *)hist_2d->GetListOfFunctions()->FindObject("stats");
    // ps->SetX1NDC(0.15);
    // TPaveStats *ps = (TPaveStats*) gPad->GetPrimitive("stats");
    // ps->SetX1NDC(0.15);
    

    std::string output_name = Form("%s/%s/RM_corr_%s_PT2050_finebins%s.pdf", outputdir.c_str(), observable.c_str(), observable.c_str(), logstring.c_str());
    can->SaveAs(output_name.c_str());

    // ==================================================================

    fin->cd();
    std::string original_4d_histname = Form("h4D_corr_%s_JetPt_R0.4_1.0Scaled", observable.c_str());
    THnSparse * thnsparse_4d = (THnSparse *)gDirectory->Get(original_4d_histname.c_str());

    // get the 2D hist for (obs det, obs truth)
    TH2D * hist2D = thnsparse_4d->Projection(2, 3);
    hist2D->GetZaxis()->SetRangeUser(1e-6,1e4);

    // save 
    TCanvas *can2 = new TCanvas();
    if (logz_bool) gPad->SetLogz();
    hist2D->Draw("colz");
    // move_stat_box(can2, hist2D);

    output_name = Form("%s/%s/RM_corr_%s_finebins%s.pdf", outputdir.c_str(), observable.c_str(), observable.c_str(), logstring.c_str());
    can2->SaveAs(output_name.c_str());


    // save - make z axis smaller
    TCanvas *can3 = new TCanvas();
    if (logz_bool) gPad->SetLogz();
    hist2D->GetZaxis()->SetRangeUser(1e-4,1e4);
    hist2D->Draw("colz");
    // move_stat_box(can3, hist2D);

    output_name = Form("%s/%s/RM_corr_%s_finebins%s_limitz.pdf", outputdir.c_str(), observable.c_str(), observable.c_str(), logstring.c_str());
    can3->SaveAs(output_name.c_str());

    // ==================================================================

    gStyle->SetOptStat(0);
    for (int i=1; i<50; i++) {

        hist2D->GetXaxis()->SetRangeUser(i - 0.2, i + 0.2);
        TH1D * hist_1 = hist2D->ProjectionY();
        hist_1->Rebin(10);
        hist_1->SetMarkerColor(kBlue);
        hist_1->SetMarkerStyle(8);
        hist_1->SetLineColor(kBlue);
        hist_1->SetMarkerSize(1);

        // fit a gaussian
        double lowx = i < 5 ? 0 : i-5;
        double highx = i > 45 ? i + 5 : 50;
        double sigma = fit_histogram_guassianfit(hist_1, lowx, highx);

        TCanvas *can3a = new TCanvas();
        // gPad->SetLogx();
        gPad->SetLogy();
        hist_1->Draw();

        if ( sigma > 0 ) {
            TLatex latex;
            latex.SetNDC();  // Use Normalized Device Coordinates (0 to 1)
            latex.SetTextSize(0.035); // Set text size
            double xtextpos = i <= 25 ? 0.7 : 0.3;
            double ytextpos = 0.8;
            latex.DrawLatex(xtextpos, ytextpos, Form("Fit #sigma = %.3f", sigma)); // (x, y, text)
        }

        can3a->SaveAs(Form("%s/%s/projections/projected_%s_jetpttruth%dGeV.pdf", outputdir.c_str(), observable.c_str(), observable.c_str(), i));

    }
    gStyle->SetOptStat(1);
    // hist2D->GetXaxis()->SetRangeUser(0.8,1.2);
    // TH1D * hist_1 = hist2D->ProjectionY();
    // hist_1->SetMarkerColor(kBlue);
    // hist_1->SetMarkerSize(2);

    // TCanvas *can3a = new TCanvas();
    // // gPad->SetLogx();
    // gPad->SetLogy();
    // hist_1->Draw();
    // can3a->SaveAs(Form("%s/%s/projections/projected_%s_jetpttruth1GeV.pdf", outputdir.c_str(), observable.c_str(), observable.c_str()));

    // hist2D->GetXaxis()->SetRangeUser(2.8,3.2);
    // TH1D * hist_2 = hist2D->ProjectionY();
    // hist_2->SetMarkerColor(kBlue);
    // hist_2->SetMarkerSize(2);

    // TCanvas *can3b = new TCanvas();
    // // gPad->SetLogx();
    // gPad->SetLogy();
    // hist_2->Draw();
    // can3b->SaveAs(Form("%s/%s/projected_%s_jetpttruth3GeV.pdf", outputdir.c_str(), observable.c_str(), observable.c_str()));

}

//===========================================================================
// 
// 
//===========================================================================
void extract_2DRM_andplot(TFile *fin, std::string observable, int ptbin, int RLbin) {

    std::string original_histname = Form("hResponse_JetPt_corr_%s_PTBIN%d_RLBIN%d_R0.4_1.0Scaled", observable.c_str(), ptbin, RLbin);
    cout << "histname" << original_histname << endl;
    // hResponse_JetPt_corr_charge_PTBIN0_RLBIN4_R0.4_1.0
    fin->cd();
    THnSparse * thnsparse = (THnSparse *)gDirectory->Get(original_histname.c_str());

    // get the 2D hist for (obs det, obs truth)
    // thnsparse->GetAxis(0)->SetRangeUser(20,40); // cut on det pt
    TH2D * hist2D = thnsparse->Projection(2, 3);
    hist2D->GetZaxis()->SetRangeUser(1e-6,1e4);

    // save 
    TCanvas *can = new TCanvas();
    gPad->SetLogz();
    hist2D->Draw("colz");

    std::string output_name = Form("%s/%s/RM_corr_%s_PTBIN%d_RLBIN%d.pdf", outputdir.c_str(), observable.c_str(), observable.c_str(), ptbin, RLbin);
    can->SaveAs(output_name.c_str());


    // //
    // hist2D->GetXaxis()->SetRangeUser(1,3);
    // TH1D * hist = hist2D->ProjectionY();

    // TCanvas *can2 = new TCanvas();
    // // gPad->SetLogx();
    // gPad->SetLogy();
    // hist->Draw();
    // can2->SaveAs(Form("%s/testing.pdf", outputdir.c_str()));

}

//===========================================================================
// This file takes in the following arguments:
// none?
// THIS IS GOOD FOR PERLMUTTER
//===========================================================================
void messing_around() { //int argc, char **argv) {

    gStyle->SetOptStat(0); // hide stats panel

    // // Reading arguments
    // std::string observable = argv[1];  // Keep argument as a string
    // int ptbin = atoi(argv[2]);     // Convert argument to an int
    // int rlbin = atoi(argv[3]);      // Convert argument to an int
   
    std::string base_slurmoutput_path = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/33905018/scaling/"; // LHC23a3
    // std::string base_slurmoutput_path = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/33818618/scaling/"; //wrong for now, LHC18b8
    std::string infile = base_slurmoutput_path + "AnalysisResultsFinal.root";

    // std::string new_filedir = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/unfolding/";
    // std::string new_filedir_addition = Form("%s/PTBIN%d/RLBIN%d/", observable, ptbin, rlbin);
    // std::string new_filename = new_filedir + new_filedir_addition + "AnalysisResults_Response.root";

    TFile* root_infile = new TFile(infile.c_str(), "READ");
    // TFile* root_outfile = new TFile(new_filename.c_str(), "RECREATE");

    extract_jetpt_RM(root_infile);
    
    std::string observable = "deltap";
    for (int i=0; i<3; i++) {
        for (int j = 0; j < 5; j++) {
            extract_2DRM_andplot(root_infile, observable, i, j);
            extract_2DRM_andplot(root_infile, "deltapt", i, j);
            extract_2DRM_andplot(root_infile, "deltapl", i, j);
            extract_2DRM_andplot(root_infile, "energyweights", i, j);
            extract_2DRM_andplot(root_infile, "charge", i, j);

        }
    }

    gStyle->SetOptStat(1);  // Enables stats panel with default settings
    extract_LHC23a3_plots(root_infile, observable);
    extract_LHC23a3_plots(root_infile, "deltapt");
    extract_LHC23a3_plots(root_infile, "deltapl");
    


}