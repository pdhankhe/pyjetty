// GOAL: Draw 3 1D histograms (from pair level) after making jet pt cuts
// run file: 'root example_for_emily.C'
// FIX!!: right now, the jet level thnsparse is 3 axes (jet pt, D0 pt, D0 y), and pair level is 4 axes (jet pt, D0 pt, D0 y, RL --> ∆jT). I will just cut on jet pt and ∆jt axes for now but this should be fixed for clarity.

double colors[16] = {kAzure+9, kTeal+9, kRed-4, kRed+1, kGreen+1, kBlue+1, kRed+2, kGreen+2, kBlue+2, kRed+3, kGreen+3, kBlue+3, kOrange+1, kViolet+1, kYellow+1, kCyan+1};

void plot_deltajt_foremily() {

    gStyle->SetOptStat(0);
    gStyle->SetTitleOffset(1.0,"y");
    
    const int num_parton_types = 3;
    std::string type_partons[num_parton_types] = { "inclusive", "light", "gluon" };
    const int num_jet_pts = 3;
    int jet_pts[num_jet_pts] = { 100, 200, 500 };

    // Initialize output file
    TFile * fout = new TFile("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/rootfiles/pythia_for_theory/PYTHIA_deltajt_hists.root", "RECREATE");

    // Loop through parton types 
    for ( int i = 0; i < num_parton_types; i++ ) {

        // Define histogram names
        std::string hname_jetlevel = Form("h_JetPt_%s_R0.4_jetlevel", type_partons[i].c_str());
        std::string hname_pairs = Form("h_corr_deltajt_JetPt_%s_R0.4", type_partons[i].c_str());

        // Make a canvas
        TCanvas * can = new TCanvas();
        can->cd();
        gPad->SetLogy();
        
        // Make a legend
        TLegend * leg = new TLegend(0.6, 0.6, 0.8, 0.8);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.04);

        // Loop through jet energies
        for ( int j = 0; j < num_jet_pts; j++ ) {

            cout << " in pt bin " << j << endl;

            // Read in your file
	        TFile * file = new TFile(Form("/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/51237766/%d/AnalysisResults.root", j+1), "READ");
            // TFile * file = new TFile("/global/cfs/cdirs/alice/blianggi/mypyjetty/analysis/testing/AnalysisResults.root", "READ");

            // Read in your histograms
            THnSparse * hsparse_jetlevel = (THnSparse*) file->Get(hname_jetlevel.c_str());
            THnSparse * hsparse_pairs = (THnSparse*) file->Get(hname_pairs.c_str());
            
            // Make cuts on pairs histogram and project
            hsparse_pairs->GetAxis(0)->SetRangeUser(jet_pts[j]-0.5,jet_pts[j]+0.5); // FIX THIS -- see above
            TH1D * hdeltajt = (TH1D * ) hsparse_pairs->Projection(3); // TH1D * hdeltajt = (TH1D * ) hsparse_jetlevel->ProjectionY();
            hdeltajt->SetDirectory(0);

            // Make cuts on jet histogram and get normalization factor & normalize
            hsparse_jetlevel->GetAxis(0)->SetRangeUser(jet_pts[j]-0.5,jet_pts[j]+0.5); // FIX THIS -- see above
            TH1D * hjetpt = (TH1D * ) hsparse_jetlevel->Projection(0);
            hjetpt->SetDirectory(0); 

            // double num_jets = hjetpt->Integral();
            // hdeltajt->Scale(1/num_jets, "width");
            double num_pairs = hdeltajt->Integral();
            hdeltajt->Scale(1/num_pairs, "width");
            
            // Format histogram
            hdeltajt->SetMarkerColorAlpha(colors[j], 1.0);
            hdeltajt->SetLineColorAlpha(colors[j], 1.0);
            hdeltajt->SetName(Form("h_corr_deltajt_%s_R0.4_pt%d", type_partons[i].c_str(), jet_pts[j]));
            hdeltajt->SetTitle(Form("#Deltaj_{T} for %s jets", type_partons[i].c_str()));
            hdeltajt->GetXaxis()->SetTitle("#Deltaj_{T}");
            hdeltajt->GetYaxis()->SetTitle("#frac{1}{N_{pairs}} #frac{dN}{d#Deltaj_{T}} #times #frac{p_{T,1}p_{T,2}}{p_{T,jet}^{2}}");
            leg->AddEntry(hdeltajt, Form("jet p_{T} = %d", jet_pts[j]), "pl");

            // Draw on canvas and save to file
            hdeltajt->Draw("SAME");
            fout->cd();
            hdeltajt->Write();

            // Close file
            file->Close();
        }

        // Add legend and save the canvas
        leg->Draw();
        can->SaveAs(Form("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC/plots/pythia_for_theory/deltajt_full_%s.pdf", type_partons[i].c_str()));    

    }

    fout->Close();
	
}