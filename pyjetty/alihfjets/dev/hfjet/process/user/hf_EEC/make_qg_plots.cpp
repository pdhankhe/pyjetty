// ROOT macro to make quark-gluon jet plots
// Ezra Lesser (elesser@berkeley.edu)

void FormatHist(TLegend *l, TH1 *hist, TString text, int markercolor=1, int markerstyle=8) 
{
    hist->SetLineColor(markercolor);
    hist->SetMarkerColor(markercolor);
    hist->SetMarkerStyle(markerstyle);
    hist->SetMarkerSize(1.5);
    l->AddEntry(hist, text, "pl");
    return;
}

void make_qg_plots() {

//    gROOT->SetBatch(); //prevents plots from showing up
    gStyle->SetOptStat(0);

    // File containing quark vs gluon histograms
    const char infile[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/11532870/AnalysisResultsFinal.root"; //perlmutter
//    const char infile[] = "/rstorage/alice/AnalysisResults/blianggi/EEC/1296078/AnalysisResultsFinal.root"; //hiccup
    //const char infile[] = "AnalysisResults.root"; //Final.root";
    TFile* f = new TFile(infile, "READ");


    // Output directory
    //std::string outdir = "/rstorage/alice/AnalysisResults/ang/1224559/plots/";
    std::string outdir = "plots/";
    // Output file for binned results
    std::string outfile = outdir + "AnalysisResultsFinal_afteranalysis.root";
    //std::string outfile = outdir + "AnalysisResultsFinal_test.root";
    TFile* f_out = new TFile(outfile.c_str(), "RECREATE");

    // Angularity alpha value
//    std::string alpha_list[] = { "1", "1.5", "2", "3" };
//    std::string grooming_list[] = { "", "_SD_zcut02_B0" };
    std::string jetR_list[] = { "0.4" };
    for (std::string jetR : jetR_list) {
//        for (std::string alpha : alpha_list) {
//            for (std::string grooming : grooming_list) {

        // Names of histograms in the file (quark, charm, gluon)
        const std::string hq_name = "h_EEC_JetPt_quark_R" + jetR;
        const std::string hc_name = "h_EEC_JetPt_charm_R" + jetR;
        // const std::string hl_name = "h_EEC_JetPt_light_R" + jetR;
        const std::string hg_name = "h_EEC_JetPt_gluon_R" + jetR;
        const std::string hq_jet_name = "h_JetPt_quark_R" + jetR + "_jetlevel";
        const std::string hc_jet_name = "h_JetPt_charm_R" + jetR + "_jetlevel";
        // const std::string hl_jet_name = "h_JetPt_light_R" + jetR + "_jetlevel";
        const std::string hg_jet_name = "h_JetPt_gluon_R" + jetR + "_jetlevel";

        // Open histograms
        TH2* hq2D = (TH2*) f->Get(hq_name.c_str());
        TH2* hc2D = (TH2*) f->Get(hc_name.c_str());
        // TH2* hl2D = (TH2*) f->Get(hl_name.c_str());
        TH2* hg2D = (TH2*) f->Get(hg_name.c_str());
        TH1* hq1D_jet = (TH1*) f->Get(hq_jet_name.c_str());
        TH1* hc1D_jet = (TH1*) f->Get(hc_jet_name.c_str());
        // TH1* hl1D_jet = (TH1*) f->Get(hl_jet_name.c_str());
        TH1* hg1D_jet = (TH1*) f->Get(hg_jet_name.c_str());


        //const int pt_bins[] = { 10, 20, 40, 60, 80, 100, 150 };
        const int pt_bins[] = { 10, 20, 40 };
        const int n_bins = 2;
        for (int i = 0; i < n_bins; i++) {
            cout << "in pt bin" << i << endl;
            int pt_min = pt_bins[i];
            int pt_max = pt_bins[i+1];

            // Set jet pT range
            hq2D->GetXaxis()->SetRangeUser(pt_min, pt_max);
            hc2D->GetXaxis()->SetRangeUser(pt_min, pt_max);
            // hl2D->GetXaxis()->SetRangeUser(pt_min, pt_max);
            hg2D->GetXaxis()->SetRangeUser(pt_min, pt_max);
            hq1D_jet->GetXaxis()->SetRangeUser(pt_min, pt_max);
            hc1D_jet->GetXaxis()->SetRangeUser(pt_min, pt_max);
            // hl1D_jet->GetXaxis()->SetRangeUser(pt_min, pt_max);
            hg1D_jet->GetXaxis()->SetRangeUser(pt_min, pt_max);
            
//            TCanvas *c2 = new TCanvas();
//            hq2D->Draw();
//            c2->SaveAs("plots/test.pdf");

            TCanvas* c = new TCanvas();
            c->cd();
            gPad->SetLogx();

            TLegend* l = new TLegend(0.1, 0.7, 0.4, 0.9);

            // Project onto observable axis
            TH1* hq_proj = (TH1*) hq2D->ProjectionY();
            TH1* hc_proj = (TH1*) hc2D->ProjectionY();
            // TH1* hl_proj = (TH1*) hl2D->ProjectionY();
            TH1* hg_proj = (TH1*) hg2D->ProjectionY();
            
            
//            TCanvas *c2 = new TCanvas();
//            hq_proj->Draw();
//            c2->SaveAs("plots/test_proj.pdf");
                        
            // Set to appropriate name
            std::string hname = hq_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hq_proj->SetNameTitle(hname.c_str(), hname.c_str());
            hname = hc_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hc_proj->SetNameTitle(hname.c_str(), hname.c_str());
            // hname = hl_proj->GetName();
            // hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            // hl_proj->SetNameTitle(hname.c_str(), hname.c_str());
            // hname = hg_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hg_proj->SetNameTitle(hname.c_str(), hname.c_str());

            
            // Rebin
            int n_obs_bins = 50; //-1;
//            double obs_bins[51]; //[60];
            double obs_bins[] = {1.00000000e-04, 1.20226443e-04, 1.44543977e-04, 1.73780083e-04, 2.08929613e-04, 2.51188643e-04, 3.01995172e-04, 3.63078055e-04, 4.36515832e-04, 5.24807460e-04, 6.30957344e-04, 7.58577575e-04, 9.12010839e-04, 1.09647820e-03, 1.31825674e-03, 1.58489319e-03, 1.90546072e-03, 2.29086765e-03, 2.75422870e-03, 3.31131121e-03, 3.98107171e-03, 4.78630092e-03, 5.75439937e-03, 6.91830971e-03, 8.31763771e-03, 1.00000000e-02, 1.20226443e-02, 1.44543977e-02, 1.73780083e-02, 2.08929613e-02, 2.51188643e-02, 3.01995172e-02, 3.63078055e-02, 4.36515832e-02, 5.24807460e-02, 6.30957344e-02, 7.58577575e-02, 9.12010839e-02, 1.09647820e-01, 1.31825674e-01, 1.58489319e-01, 1.90546072e-01, 2.29086765e-01, 2.75422870e-01, 3.31131121e-01, 3.98107171e-01, 4.78630092e-01, 5.75439937e-01, 6.91830971e-01, 8.31763771e-01, 1.00000000e+00};
//            if (alpha == "1") {
//                n_obs_bins = 16;
//                double obs_bins_tmp[] = {0, 0.02, 0.04,  0.06,  0.08,    0.10,    0.12,     0.14,   0.16,   0.20,    0.24,       0.28,      0.32,      0.36,     0.40,   0.45,      0.5};
//                std::copy(obs_bins_tmp, obs_bins_tmp+n_obs_bins+1, obs_bins);
//            } else if (alpha == "1.5") {
//                n_obs_bins = 12;
//                double obs_bins_tmp[] = {0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.3, 0.33, 0.36};
//                std::copy(obs_bins_tmp, obs_bins_tmp+n_obs_bins+1, obs_bins);
//            } else if (alpha == "2") {
//                n_obs_bins = 15;
//                double obs_bins_tmp[] = {0, 0.02, 0.04, 0.06,  0.08,  0.10, 0.12, 0.14, 0.16, 0.18,  0.20,  0.22, 0.24, 0.26, 0.28, 0.3};
//                std::copy(obs_bins_tmp, obs_bins_tmp+n_obs_bins+1, obs_bins);
//            } else {  // alpha == "3"
//                n_obs_bins = 11;
//                double obs_bins_tmp[] = {0, 0.01, 0.02, 0.04, 0.06,  0.08,  0.10, 0.12, 0.14, 0.16, 0.18,  0.20};
//                std::copy(obs_bins_tmp, obs_bins_tmp+n_obs_bins+1, obs_bins);
//            }
            
            
            cout << "hq_proj " << hq_proj->GetNbinsX() << endl;
            cout << "hc_proj " << hc_proj->GetNbinsX() << endl;
            // cout << "hl_proj " << hl_proj->GetNbinsX() << endl;
            cout << "hg_proj " << hg_proj->GetNbinsX() << endl;
            
//            TH1* hq = (TH1*) hq_proj->Rebin(n_obs_bins, (hq_name + "rebin").c_str(), obs_bins);
//            TH1* hc = (TH1*) hc_proj->Rebin(n_obs_bins, (hc_name + "rebin").c_str(), obs_bins);
//            TH1* hg = (TH1*) hg_proj->Rebin(n_obs_bins, (hg_name + "rebin").c_str(), obs_bins);
            //unnecessary renaming for EEC...
            TH1* hq = (TH1*) hq_proj->Clone((hq_name + "rebin").c_str());
            TH1* hc = (TH1*) hc_proj->Clone((hc_name + "rebin").c_str());
            // TH1* hl = (TH1*) hl_proj->Clone((hl_name + "rebin").c_str());
            TH1* hg = (TH1*) hg_proj->Clone((hg_name + "rebin").c_str());


            // Find normalization factor
            double numjets_quark = hq1D_jet->Integral(); //GetEntries();
            double numjets_charm = hc1D_jet->Integral(); //GetEntries();
            // double numjets_light = hl1D_jet->Integral(); //GetEntries();
            double numjets_gluon = hg1D_jet->Integral(); //GetEntries();
            // double numjets_total = numjets_quark + numjets_gluon;
            
            
            // Set normalization
//            hq->Scale(1/hq->Integral(), "width");
//            hc->Scale(1/hc->Integral(), "width");
//            hg->Scale(1/hg->Integral(), "width");
            hq->Scale(1/numjets_quark, "width");
            hc->Scale(1/numjets_charm, "width");
            // hl->Scale(1/numjets_light, "width");
            hg->Scale(1/numjets_gluon, "width");
            
            cout << "num jets quark " << numjets_quark << endl;
            cout << "num jets charm " << numjets_charm << endl;
            // cout << "num jets light " << numjets_light << endl;
            cout << "num jets gluon " << numjets_gluon << endl;

            /*
            // Find plotting range -- use gluon (broadest) distribution
            int maxbin = 0;
            double somevalue = 0.000001;
            int nxbins = hg->GetXaxis()->GetNbins();
            for (int j = 0; j < nxbins; j++) {
                if (hg->GetBinContent(nxbins-j) > somevalue) {
                    maxbin = nxbins-j;
                    break;
                }
            }
            hq->GetXaxis()->SetRange(1, maxbin);
            hc->GetXaxis()->SetRange(1, maxbin);
            hg->GetXaxis()->SetRange(1, maxbin);
            */
            // Set minimum / maximum y-axis value
//            hq->SetMinimum(1e-3);
            double maxy = hg->GetMaximum();
            if (hc->GetMaximum() > maxy) maxy = hc->GetMaximum();
            // if (hl->GetMaximum() > maxy) maxy = hc->GetMaximum();
            if (hq->GetMaximum() > maxy) maxy = hq->GetMaximum();
            maxy *= 1.5;
            hq->SetMaximum(maxy);
            // cout << "maximums, hq: " << hq->GetMaximum() << " hc: " << hc->GetMaximum() << " hg: " << hg->GetMaximum() << endl;
            cout << "chosen max: " << maxy << " " << maxy/1.5 << endl;
            

            // Format histograms for plotting
            FormatHist(l, hq, "Quark-initiated jets", 1, 8);
            // hq->SetLineColor(1);
            // hq->SetMarkerColor(1);
            // hq->SetMarkerStyle(8);
            // hq->SetMarkerSize(1.5);
            // l->AddEntry(hq, "Quark-initiated jets", "pl");

            FormatHist(l, hc, "Charm-initiated jets", 2, 20);
            // hc->SetLineColor(2);
            // hc->SetMarkerColor(2);
            // hc->SetMarkerStyle(20);
            // hc->SetMarkerSize(1.5);
            // l->AddEntry(hc, "Charm-initiated jets", "pl");

            // FormatHist(l, hl, "Light-initiated jets", 3, 22);

            FormatHist(l, hg, "Gluon-initiated jets", 4, 33);
            // hg->SetLineColor(4);
            // hg->SetMarkerColor(4);
            // hg->SetMarkerStyle(33);
            // hg->SetMarkerSize(1.5);
            // l->AddEntry(hg, "Gluon-initiated jets", "pl");


            // // Fix the x-axis label
            // hq->GetXaxis()->SetTitle("R_{L}");
            // hc->GetXaxis()->SetTitle("R_{L}");
            // hg->GetXaxis()->SetTitle("R_{L}");

            // Create plot
            hq->Draw("L");
            hc->Draw("L same");
            // hl->Draw("L same");
            hg->Draw("L same");
            l->Draw("same");

            // Label plots
            TText* t = new TText(0.3, 11, "PYTHIA8");
            t->SetTextAlign(22);
            t->SetTextFont(43);
            t->SetTextSize(20);
//            t->Draw("same");
            //std::string text = "#it{#alpha} = " + alpha;
            //t->SetText(0.3, 10, text.c_str());

            std::string fname = outdir + "QG_comp_pt" + std::to_string(pt_min) + '-' + std::to_string(pt_max) + "_R" + jetR + ".pdf"; // + "_normbytype.pdf"; //"_nonorm.pdf";
            const char* fnamec = fname.c_str();
            c->SaveAs(fnamec);
            delete c;
            delete t;

            if (pt_min == 10) { //} && grooming == "") {
                // Write rebinned histograms to root file
                f_out->cd();
                hq->Write();
                hg->Write();
            }
             
        } // pT bins loop
//            } // grooming loop
//        } // alpha loop
    } // jetR loop

    f->Close();
    delete f;

    return;
}
