// ROOT macro to make quark-gluon jet plots
// Ezra Lesser (elesser@berkeley.edu)

#include<iostream>

void FormatHist(TLegend *l, TH1 *hist, TString text, int markercolor=1, int markerstyle=8) 
{
    hist->SetLineColor(markercolor);
    hist->SetMarkerColor(markercolor);
    hist->SetMarkerStyle(markerstyle);
    hist->SetMarkerSize(1.5);
    l->AddEntry(hist, text, "pl");
    return;
}

int findTopOfCurve(TH1* hist) {
    
    int val = 0;

    //for each point, find the slope
    int numbins = hist->GetNbinsX();
    int binstart = hist->FindBin(0.01);
    for (int i=binstart; i<numbins; i++) {
        double y1 = hist->GetBinContent(i);
        double y2 = hist->GetBinContent(i+1);
        // double x1 = hist->GetBinCenter(i); //not needed
        // double x2 = hist->GetBinCenter(i+1); //not needed

        double slope_num = y2-y1;
        if (slope_num < 0) {
            // double x1 = hist->GetBinCenter(i);
            val = i;
            break;
        }
    }

    return val; //return the bin number
}

TLine * drawVertLine(double x1, double y1, double y2, int color){
    auto fvertline = new TLine(x1, y1, x1, y2);
	fvertline->SetLineWidth(1);
    fvertline->SetLineColor(color);
    fvertline->SetLineStyle(2);
    return fvertline;

}


void make_qg_plots_comparison() {

//    gROOT->SetBatch(); //prevents plots from showing up
    gStyle->SetOptStat(0);
    gStyle->SetLabelOffset(0.01,"y");
    Double_t markers[10] = {kFullCircle, kFullTriangleUp, kFullDiamond, kFullSquare, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
    Double_t marker_size = 1.2;
    Double_t colors[19] = {kRed, kGreen, kBlue, kRed+1, kGreen+1, kBlue+1, kRed+2, kGreen+2, kBlue+2, kRed+3, kGreen+3, kBlue+3, kOrange, kYellow, kViolet, kOrange+1, kYellow+1, kViolet+1, kCyan+1};


    // File containing quark vs gluon histograms

    //FOR WHEN WEIGHTED/UNWEIGHTED IN SAME FILE
    const char base_path_infile[] = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/"; //perlmutter 
    // Leading pt cuts:
    // const char* run_nums[] = {"11731456", "11853285", "11853226", "11853325"};
    // const std::string labels[] = {"no leading pT cut", "leading pT cut = 4 GeV", "leading pT cut = 5.33 GeV", "leading pT cut = 6 GeV"};
    // replaceKP pairs comparison:
    const char* run_nums[] = {"11731456", "11731473", "13777236"};
    const std::string labels[] = {"charm decays off", "charm decays on", "D0->Kpi reconstruction"};
    const char infile_name[] = "/AnalysisResultsFinal.root";

    bool unweighted = false; //true = draw unweighted, false = do not draw
    bool switchON = false; //true = charm ON, false = charm OFF
    bool replaceKPON = true; //true = D0 reconstruction implemented
    // TFile* f;
    std::string add_name;
    if (replaceKPON) {
        add_name = "_replaceKP_compall.pdf"; 
    } else {
        if (switchON) {
            // f = new TFile(infile_charmON, "READ");
            add_name = "_charmdecaysON_leadingpTcutcomp.pdf"; 
        } else {
            // f = new TFile(infile_charmOFF, "READ");
            add_name = "_charmdecaysOFF_leadingpTcutcomp.pdf";
        }
    }

    // Output directory
    //std::string outdir = "/rstorage/alice/AnalysisResults/ang/1224559/plots/";
    std::string outdir = "plots/";
    // Output file for binned results
    std::string outfile = outdir + "AnalysisResultsFinal_afteranalysis.root";
    //std::string outfile = outdir + "AnalysisResultsFinal_test.root";
    TFile* f_out = new TFile(outfile.c_str(), "RECREATE");


    std::string jetR = "0.4";


    // loop over pt bins
    const int pt_bins[] = { 10, 20, 40 };
    const int n_bins = 2;
    for (int i = 0; i < n_bins; i++) {
        cout << "in pt bin" << i << endl;
        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];

        TCanvas* c = new TCanvas();
        c->cd();
        gPad->SetLogx();
        // gPad->SetLogy();

        TLegend* l = new TLegend(0.1, 0.7, 0.4, 0.9);


        // loop over the files here
        int numfiles = sizeof(run_nums)/sizeof(run_nums[0]);
        cout << "Looping over " << numfiles << " files..." << endl;
        for ( int ifile = 0; ifile < numfiles; ifile++ ) {
            
            std::string full_infile = std::string(base_path_infile) + std::string(run_nums[ifile]) + std::string(infile_name);
            const char* full_infilec = full_infile.c_str();
            cout << full_infilec << endl;
            TFile* f = new TFile(full_infilec, "READ");
            

            // Names of histograms in the file (quark, charm, gluon)
            const std::string hc_name = "h_EEC_JetPt_charm_R" + jetR;
            const std::string hl_name = "h_EEC_JetPt_light_R" + jetR;
            const std::string hg_name = "h_EEC_JetPt_gluon_R" + jetR;
            const std::string hc_unweighted_name = "h_EEC_JetPt_charm_R" + jetR + "_unweighted";
            const std::string hl_unweighted_name = "h_EEC_JetPt_light_R" + jetR + "_unweighted";
            const std::string hg_unweighted_name = "h_EEC_JetPt_gluon_R" + jetR + "_unweighted";
            const std::string hc_jet_name = "h_JetPt_charm_R" + jetR + "_jetlevel";
            const std::string hl_jet_name = "h_JetPt_light_R" + jetR + "_jetlevel";
            const std::string hg_jet_name = "h_JetPt_gluon_R" + jetR + "_jetlevel";

            // Open histograms
            TH2* hc2D = (TH2*) f->Get(hc_name.c_str());
            TH2* hl2D = (TH2*) f->Get(hl_name.c_str());
            TH2* hg2D = (TH2*) f->Get(hg_name.c_str());
            TH2* hc2D_unweighted = (TH2*) f->Get(hc_unweighted_name.c_str());//->Get(hc_unweighted_name.c_str());
            TH2* hl2D_unweighted = (TH2*) f->Get(hl_unweighted_name.c_str());//->Get(hl_unweighted_name.c_str());
            TH2* hg2D_unweighted = (TH2*) f->Get(hg_unweighted_name.c_str());//->Get(hg_unweighted_name.c_str());
            TH1* hc1D_jet = (TH1*) f->Get(hc_jet_name.c_str());
            TH1* hl1D_jet = (TH1*) f->Get(hl_jet_name.c_str());
            TH1* hg1D_jet = (TH1*) f->Get(hg_jet_name.c_str());


            // const int pt_bins[] = { 10, 20, 40 };
            // const int n_bins = 2;
            // for (int i = 0; i < n_bins; i++) {
            //     cout << "in pt bin" << i << endl;
            //     int pt_min = pt_bins[i];
            //     int pt_max = pt_bins[i+1];

            // Set jet pT range
            hc2D->GetXaxis()->SetRangeUser(pt_min, pt_max);
            hl2D->GetXaxis()->SetRangeUser(pt_min, pt_max);
            hg2D->GetXaxis()->SetRangeUser(pt_min, pt_max);
            hc2D_unweighted->GetXaxis()->SetRangeUser(pt_min, pt_max);
            hl2D_unweighted->GetXaxis()->SetRangeUser(pt_min, pt_max);
            hg2D_unweighted->GetXaxis()->SetRangeUser(pt_min, pt_max);
            hc1D_jet->GetXaxis()->SetRangeUser(pt_min, pt_max);
            hl1D_jet->GetXaxis()->SetRangeUser(pt_min, pt_max);
            hg1D_jet->GetXaxis()->SetRangeUser(pt_min, pt_max);
            
//            TCanvas *c2 = new TCanvas();
//            hq2D->Draw();
//            c2->SaveAs("plots/test.pdf");

            // Project onto observable axis
            TH1* hc = (TH1*) hc2D->ProjectionY(); //_proj = (TH1*) hc2D->ProjectionY();
            TH1* hl = (TH1*) hl2D->ProjectionY(); //_proj = (TH1*) hl2D->ProjectionY();
            TH1* hg = (TH1*) hg2D->ProjectionY(); //_proj = (TH1*) hg2D->ProjectionY();
            TH1* hc_unweighted = (TH1*) hc2D_unweighted->ProjectionY();
            TH1* hl_unweighted = (TH1*) hl2D_unweighted->ProjectionY();
            TH1* hg_unweighted = (TH1*) hg2D_unweighted->ProjectionY();

            
//            TCanvas *c2 = new TCanvas();
//            hq_proj->Draw();
//            c2->SaveAs("plots/test_proj.pdf");
                        
            // Set to appropriate name
            std::string hname = hc->GetName(); //_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hc->SetNameTitle(hname.c_str(), hname.c_str());
            hname = hl->GetName(); //_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hl->SetNameTitle(hname.c_str(), hname.c_str());
            hname = hg->GetName(); //_proj->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hg->SetNameTitle(hname.c_str(), hname.c_str());

            hname = hc_unweighted->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hc_unweighted->SetNameTitle(hname.c_str(), hname.c_str());
            hname = hl_unweighted->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hl_unweighted->SetNameTitle(hname.c_str(), hname.c_str());
            hname = hg_unweighted->GetName();
            hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
            hg_unweighted->SetNameTitle(hname.c_str(), hname.c_str());

            // Rebin
            int n_obs_bins = 50; //-1;
//            double obs_bins[51]; //[60];
            double obs_bins[] = {1.00000000e-04, 1.20226443e-04, 1.44543977e-04, 1.73780083e-04, 
                2.08929613e-04, 2.51188643e-04, 3.01995172e-04, 3.63078055e-04, 4.36515832e-04, 
                5.24807460e-04, 6.30957344e-04, 7.58577575e-04, 9.12010839e-04, 1.09647820e-03, 
                1.31825674e-03, 1.58489319e-03, 1.90546072e-03, 2.29086765e-03, 2.75422870e-03, 
                3.31131121e-03, 3.98107171e-03, 4.78630092e-03, 5.75439937e-03, 6.91830971e-03, 
                8.31763771e-03, 1.00000000e-02, 1.20226443e-02, 1.44543977e-02, 1.73780083e-02, 
                2.08929613e-02, 2.51188643e-02, 3.01995172e-02, 3.63078055e-02, 4.36515832e-02, 
                5.24807460e-02, 6.30957344e-02, 7.58577575e-02, 9.12010839e-02, 1.09647820e-01, 
                1.31825674e-01, 1.58489319e-01, 1.90546072e-01, 2.29086765e-01, 2.75422870e-01, 
                3.31131121e-01, 3.98107171e-01, 4.78630092e-01, 5.75439937e-01, 6.91830971e-01, 
                8.31763771e-01, 1.00000000e+00};
            
            // cout << "hc_proj " << hc_proj->GetNbinsX() << endl;
            // cout << "hl_proj " << hl_proj->GetNbinsX() << endl;
            // cout << "hg_proj " << hg_proj->GetNbinsX() << endl;
            

            //unnecessary renaming for EEC...
            // TH1* hc = (TH1*) hc_proj->Clone((hc_name + "rebin").c_str());
            // TH1* hl = (TH1*) hl_proj->Clone((hl_name + "rebin").c_str());
            // TH1* hg = (TH1*) hg_proj->Clone((hg_name + "rebin").c_str());


            // Find normalization factor
            double numjets_charm = hc1D_jet->Integral(); //GetEntries();
            double numjets_light = hl1D_jet->Integral(); //GetEntries();
            double numjets_gluon = hg1D_jet->Integral(); //GetEntries();
            
            
            // Set normalization
//            hq->Scale(1/hq->Integral(), "width");
//            hc->Scale(1/hc->Integral(), "width");
//            hg->Scale(1/hg->Integral(), "width");
            hc->Scale(1/numjets_charm, "width");
            hl->Scale(1/numjets_light, "width");
            hg->Scale(1/numjets_gluon, "width");
            hc_unweighted->Scale(1/numjets_charm, "width");
            hl_unweighted->Scale(1/numjets_light, "width");
            hg_unweighted->Scale(1/numjets_gluon, "width");
            // cout << "num jets quark " << numjets_quark << endl;
            cout << "num jets charm " << numjets_charm << endl;
            cout << "num jets light " << numjets_light << endl;
            cout << "num jets gluon " << numjets_gluon << endl;


            // Calculate where the slope changes sign to find the max
            cout << "top of hc" << hc->GetBinCenter(findTopOfCurve(hc)) << endl;
            cout << "top of hl" << hl->GetBinCenter(findTopOfCurve(hl)) << endl;
            cout << "top of hg" << hg->GetBinCenter(findTopOfCurve(hg)) << endl;
            cout << "top of hc_unweighted" << hc_unweighted->GetBinCenter(findTopOfCurve(hc_unweighted)) << endl;
            cout << "top of hl_unweighted" << hl_unweighted->GetBinCenter(findTopOfCurve(hl_unweighted)) << endl;
            cout << "top of hg_unweighted" << hg_unweighted->GetBinCenter(findTopOfCurve(hg_unweighted)) << endl;




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
            if (hl->GetMaximum() > maxy) maxy = hl->GetMaximum();
            if (unweighted) {
                if (hc_unweighted->GetMaximum() > maxy) maxy = hc_unweighted->GetMaximum();
                if (hl_unweighted->GetMaximum() > maxy) maxy = hl_unweighted->GetMaximum();
                if (hg_unweighted->GetMaximum() > maxy) maxy = hg_unweighted->GetMaximum();
            }
            maxy *= 1.5;
            hc->SetMaximum(maxy);
            cout << "maximums, hl: " << hl->GetMaximum() << " hc: " << hc->GetMaximum() << " hg: " << hg->GetMaximum() << endl;
            cout << "maximums, hl_unweighted: " << hl_unweighted->GetMaximum() << " hc_unweighted: " << hc_unweighted->GetMaximum() << " hg_unweighted: " << hg_unweighted->GetMaximum() << endl;
            cout << "chosen max: " << maxy << " " << maxy/1.5 << endl;
            

            // Format histograms for plotting
            cout << "labels[i] " << labels[i] << endl;
            FormatHist(l, hc, "Charm-initiated jets, " + labels[ifile], colors[3*ifile], markers[0]); //marker=20
            FormatHist(l, hl, "Light-initiated jets, " + labels[ifile], colors[3*ifile+1], markers[1]); //marker=22
            FormatHist(l, hg, "Gluon-initiated jets, " + labels[ifile], colors[3*ifile+2], markers[2] ); //marker=33

            // cout << "hists formated" << endl;
// 
            // Fix the x-axis label
            if (ifile==0){
                hc->GetXaxis()->SetTitle("R_{L}");
                hc->GetYaxis()->SetTitle("#frac{1}{N_{jet}} #times #frac{dN_{EEC}}{dR_{L}}");
            }

            // cout << "labels fixed " << endl;

            // TCanvas *c2 = new TCanvas();
            // hc->Draw("L same");
            // hl->Draw("L same");
            // hg->Draw("L same");
            // c2->SaveAs("plots/test.pdf");

            // Create plot
            hc->Draw("L same");
            hl->Draw("L same");
            hg->Draw("L same");
            // gPad->Update();

            // cout << "hists drawn" << endl;

            // draw center
            double hc_top_binpos = findTopOfCurve(hc);
            double hl_top_binpos = findTopOfCurve(hl);
            double hg_top_binpos = findTopOfCurve(hg);
            // drawVertLine(hc->GetBinCenter(hc_top_binpos), 0, hc->GetBinContent(hc_top_binpos), kRed)->Draw();
            // drawVertLine(hl->GetBinCenter(hl_top_binpos), 0, hl->GetBinContent(hl_top_binpos), kGreen)->Draw();
            // drawVertLine(hg->GetBinCenter(hg_top_binpos), 0, hg->GetBinContent(hg_top_binpos), kBlue)->Draw();
            
            // cout << "top of curve found" << endl;

            // if we want unweighted plotted
            if (unweighted){
                cout << "shouldn't be in here" << endl;
                FormatHist(l, hc_unweighted, "Charm-initiated jets, unweighted", kRed+2, 24); //fix the colors and markers on these
                FormatHist(l, hl_unweighted, "Light-initiated jets, unweighted", kGreen+3, 26);
                FormatHist(l, hg_unweighted, "Gluon-initiated jets, unweighted", kBlue+4, 27);

                hc_unweighted->Draw("L same");
                hl_unweighted->Draw("L same");
                hg_unweighted->Draw("L same");

                double hc_unweighted_top_binpos = findTopOfCurve(hc_unweighted);
                double hl_unweighted_top_binpos = findTopOfCurve(hl_unweighted);
                double hg_unweighted_top_binpos = findTopOfCurve(hg_unweighted);
                drawVertLine(hc_unweighted->GetBinCenter(hc_unweighted_top_binpos), 0, hc_unweighted->GetBinContent(hc_unweighted_top_binpos), kRed)->Draw();
                drawVertLine(hl_unweighted->GetBinCenter(hl_unweighted_top_binpos), 0, hl_unweighted->GetBinContent(hl_unweighted_top_binpos), kGreen)->Draw();
                drawVertLine(hg_unweighted->GetBinCenter(hg_unweighted_top_binpos), 0, hg_unweighted->GetBinContent(hg_unweighted_top_binpos), kBlue)->Draw();
            } 
            
            
            if (pt_min == 10) { //} && grooming == "") {
                // Write rebinned histograms to root file
                f_out->cd();
                hl->Write();
                hg->Write();
            }

            // f->Close();
            // delete f;

            // cout << "file closed" << endl;
            
             
        } // end of loop over files

        // cout << "out of files loop" << endl;

        // Draw legend
        l->Draw("same");

        // Label plots
        // TText* t = new TText(0.3, 11, "PYTHIA8");
        // t->SetTextAlign(22);
        // t->SetTextFont(43);
        // t->SetTextSize(20);
        // t->Draw("same");
        //std::string text = "#it{#alpha} = " + alpha;
        //t->SetText(0.3, 10, text.c_str());


        std::string fname = outdir + "QG_comp_pt" + std::to_string(pt_min) + '-' + std::to_string(pt_max) + "_R" + jetR + add_name; //"_charmdecaysONcomparison.pdf"; // + "_normbytype.pdf"; //"_nonorm.pdf";
        const char* fnamec = fname.c_str();
        // cout << "fnamec is " << fnamec;
        // gPad->Update();
        // cout << "gpad updated" << endl;
        c->SaveAs(fnamec);
        delete c;
        // delete t;




    } // end of pT bin loop


    return;
}
