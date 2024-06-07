// ROOT macro to make quark-gluon jet plots
// Beatrice Liang-Gilman (beatrice_lg@berkeley.edu)

void SetStyle(Bool_t graypalette=true) {
  	cout << "Setting style!" << endl;
  
  	gStyle->Reset("Plain");
  	gStyle->SetOptTitle(0);
  	gStyle->SetOptStat(0);
  	if(graypalette) gStyle->SetPalette(8,0);
  	else gStyle->SetPalette(1);
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

void ProcessCanvas(TCanvas *Canvas) { 
	gStyle->SetOptStat(0);
	Canvas->SetHighLightColor(1);
	Canvas->SetFillColor(0);
	Canvas->SetBorderMode(0);
	Canvas->SetBorderSize(2);
	Canvas->SetTickx(1);
	Canvas->SetTicky(1);
	Canvas->SetFrameBorderMode(0);
	Canvas->SetFrameLineWidth(1);
 	Canvas->SetFrameBorderMode(1);
}

void FormatHist(TLegend *l, TH1 *hist, TString text, int markercolor=1, int markerstyle=8) 
{
    hist->SetLineColor(markercolor);
    hist->SetMarkerColor(markercolor);
    hist->SetMarkerStyle(markerstyle);
    hist->SetMarkerSize(1.5);
    l->AddEntry(hist, text, "pl");

	//gPad->SetTickx(); 
	//gPad->SetTicky(); 
	// h->SetLineWidth(2);
	hist->GetYaxis()->SetTitleOffset(1.05); 
	hist->GetYaxis()->SetTitleSize(0.06); //(0.042);
	hist->GetYaxis()->SetLabelSize(0.05); //(0.042);
	hist->GetYaxis()->SetLabelFont(42);
	hist->GetXaxis()->SetLabelFont(42);
	hist->GetYaxis()->SetTitleFont(42);
	hist->GetXaxis()->SetTitleFont(42);
	hist->GetXaxis()->SetTitleOffset(1.0);
	hist->GetXaxis()->SetTitleSize(0.06); //(0.042);
	hist->GetXaxis()->SetLabelSize(0.05); //(0.042);


    return;
}

//checkExtra is how many point-to-point slopes after finding a decreasing slope I want to check
int findTopOfCurve(TH1* hist, int checkExtra=1) {
    
    //for each point, find the slope from the 
    int numbins = hist->GetNbinsX();
    int binstart = hist->FindBin(0.01);
    bool falsealarm = false;

    for (int i=binstart; i<numbins; i++) {

        double y1 = hist->GetBinContent(i);
        double y2 = hist->GetBinContent(i+1);
        double slope_num = y2-y1;

        if (slope_num < 0) {
            for (int j=1; j<=checkExtra; j++) {
                y1 = hist->GetBinContent(i+j);
                y2 = hist->GetBinContent(i+1+j);
                slope_num = y2-y1;
                if (slope_num >= 0) {
                    falsealarm = true;
                    break;
                }
            }
            if (!falsealarm) return i; //return the bin number
        }
        falsealarm = false;
    }

    return 0; //return 0 if nothing found
}

// find the full width at half max
// returns a vector that gives (halfmax_height, leftPos, rightPos, full_width)
vector<double> findWidthOfCurve(TH1* hist, double topofCurve_bin) {

    vector<double> FWHM_vec;
    double topofCurve = hist->GetBinContent(topofCurve_bin);
    double halfmax_height = topofCurve/2;
    double leftPos = 0;
    double rightPos = 0;
    double full_width = 0;

    //find the left point of the half max width
    for (int i=topofCurve_bin-1; i>0.; i--) {
        double binheight = hist->GetBinContent(i);
        if (binheight <= halfmax_height) {
            leftPos = hist->GetBinCenter(i);
            break;
        }
    }
    //find the right point of the half max width
    for (int i=topofCurve_bin+1; i<hist->GetNbinsX(); i++) {
        double binheight = hist->GetBinContent(i);
        if (binheight <= halfmax_height) {
            rightPos = hist->GetBinCenter(i);
            break;
        }
    }

    full_width = rightPos - leftPos;

    FWHM_vec.push_back(halfmax_height);
    FWHM_vec.push_back(leftPos);
    FWHM_vec.push_back(rightPos);
    FWHM_vec.push_back(full_width);
    
    return FWHM_vec;
    
}

TLine * drawVertLine(double x1, double y1, double y2, int color, int linestyle=2){
    auto fvertline = new TLine(x1, y1, x1, y2);
	fvertline->SetLineWidth(1);
    fvertline->SetLineColor(color);
    fvertline->SetLineStyle(linestyle);
    return fvertline;

}

TLine * drawHoriLine(double x1, double x2, double y1, int color, int linestyle=2){
    auto fhoriline = new TLine(x1, y1, x2, y1);
	fhoriline->SetLineWidth(1);
    fhoriline->SetLineColor(color);
    fhoriline->SetLineStyle(linestyle);
    return fhoriline;

}


void make_herwig_plots_write() {

//    gROOT->SetBatch(); //prevents plots from showing up
    gStyle->SetOptStat(0);
    SetStyle();
    Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
    Double_t marker_size = 1.5;
    Double_t colors[16] = {kRed, kGreen+2, kBlue, kRed+1, kGreen+1, kBlue+1, kRed+2, kGreen+2, kBlue+2, kRed+3, kGreen+3, kBlue+3, kOrange+1, kViolet+1, kYellow+1, kCyan+1};
    Double_t marker_color[3] = {kMagenta, kMagenta+2, kMagenta+3};

    // File containing quark vs gluon histograms

     //FOR WHEN WEIGHTED/UNWEIGHTED IN SAME FILE
    // const char infile_D0_noDstar[] = "/rstorage/generators/herwig_alice/histograms/249612/216486/AnalysisResultsFinal.root"; //this is using thnsparse
    // const char infile_D0[] = "/rstorage/generators/herwig_alice/histograms/277508/271758/AnalysisResultsFinal.root"; //this is using thnsparse
    std::string infile_D0 = "/rstorage/generators/herwig_alice/scaling/343985/299990/AnalysisResultsFinal.root"; //this is using thnsparse
    std::string infile_D0_noDstar = "/rstorage/generators/herwig_alice/scaling/344485/299990/AnalysisResultsFinal.root"; //this is using thnsparse
    std::string infile_inclusive = "/rstorage/generators/herwig_alice/scaling/352629/329854/AnalysisResultsFinal.root";

    // int plot_case:
    // 0 = plot everything in a separate file, currently not set yet

    //CONTOL VARIABLES HERE
    int plot_case = 0;


    // Output directory
    std::string outdir = "plots/herwig/";
    // Output file for binned results
    std::string outfile = outdir + "AnalysisResultsFinal_herwig" + ".root"; 
    TFile* f_out = new TFile(outfile.c_str(), "RECREATE");

    int ifile = 0;
    std::string file_list[] = { infile_D0, infile_D0_noDstar, infile_inclusive };
    std::string dname_list[] = { "D0", "D0_noDstar", "inclusive" };
    for (std::string file : file_list) {
        
        TFile* f = new TFile(file.c_str(), "READ");
        std::string dname = dname_list[ifile];
        TString label1 = "";


        // Jet R value
        std::string jetR_list[] = { "0.4" };
        for (std::string jetR : jetR_list) {
            
            // different track threshold values
            std::string trkthrd_list[] = { "0.15", "0.5", "1.0" };
            int itrkthrd = 0;
            for (std::string trkthrd : trkthrd_list) {
                const std::string EEC_name = "h_jet_ENC_RL2_JetPt_Truth_R" + jetR + "_" + trkthrd + "Scaled";
                const std::string EEC_ptrl_name = "h_jet_ENC_RL2Pt_JetPt_Truth_R" + jetR + "_" + trkthrd + "Scaled";
                const std::string EEC_noweight_name = "h_jet_EEC_noweight_RL_JetPt_Truth_R" + jetR + "_" + trkthrd + "Scaled";
                const std::string jetinfo_name = "h_jet_pt_JetPt_Truth_R" + jetR + "_" + trkthrd + "Scaled";

                //-------------------------------------------------//
                // find D0 reconstruction through charm
                THnSparse* hsparse_EEC = (THnSparse*) f->Get(EEC_name.c_str());
                THnSparse* hsparse_EEC_ptrl = (THnSparse*) f->Get(EEC_ptrl_name.c_str());
                THnSparse* hsparse_EEC_noweight = (THnSparse*) f->Get(EEC_noweight_name.c_str());
                THnSparse* hsparse_jetinfo = (THnSparse*) f->Get(jetinfo_name.c_str());
                cout << "jet info name " << jetinfo_name.c_str() << endl;

                // testing - look at # jets before cuts
                cout << "numDtaggedjets from hist before cuts " << hsparse_jetinfo->Projection(0)->GetEntries() << endl;

                //==================================================//

                // save D0 pt spectrum
                TH1D *hD0_pT = hsparse_jetinfo->Projection(1); //use jet level
                hD0_pT->SetNameTitle( (dname + "_trkthrd" + trkthrd + "_D0_pt").c_str() ,  (dname + "_trkthrd" + trkthrd + "_D0_pt").c_str());
                double numjets_d0 = hD0_pT->Integral();
                hD0_pT->Scale(1/numjets_d0, "width");

                // save D0 z spectrum
                TH1D *hD0_z = hsparse_jetinfo->Projection(3); //use jet level
                hD0_z->SetNameTitle( (dname + "_trkthrd" + trkthrd + "_D0_z").c_str() ,  (dname + "_trkthrd" + trkthrd + "_D0_z").c_str());
                hD0_z->Scale(1/numjets_d0, "width");

                //==================================================//


                // pt bins
                const int pt_bins[] = { 7, 10, 15, 30 }; //{ 10, 20, 40 }; // CHANGE HERE!!
                const int d0_pt_cuts[] = { 3, 5, 5 };
                const int n_bins = 3; 
                for (int i = 0; i < n_bins; i++) {
                    cout << "in pt bin" << i << endl;
                    int pt_min = pt_bins[i];
                    int pt_max = pt_bins[i+1];

                    // define pt related variables
                    TString ptbin = TString::Format("%d #leq #it{p}_{T}^{ch. jet} < %d GeV/#it{c}, #font[122]{|}#it{#eta}_{jet}#font[122]{|} #leq 0.5", pt_min, pt_max);
                    TString ptD = TString::Format("%d #leq #it{p}_{T}^{D^{0}} < %d GeV/#it{c}, #font[122]{|}#it{y}_{D^{0}}#font[122]{|} #leq 0.8", d0_pt_cuts[i], pt_max);


                    // make a canvas for each pt range
                    TCanvas* c_EEC = new TCanvas();
                    ProcessCanvas(c_EEC);
                    // c_EEC->cd();
                    gPad->SetLogx();
                    // gPad->SetLogy();

                    TCanvas* c_EEC_ptrl = new TCanvas();
                    ProcessCanvas(c_EEC_ptrl);
                    gPad->SetLogx();

                    TCanvas* c_EEC_noweight = new TCanvas();
                    ProcessCanvas(c_EEC_noweight);
                    gPad->SetLogx();


                    TLegend* l; // = new TLegend(0.17, 0.65, 0.5, 0.85);
                    TLegend* ldummy = new TLegend(0.1797168,0.400741,0.4562155,0.8885185,"");

                    double maxy = 0;


                    // Open histograms


                    l = new TLegend(0.1797168,0.400741,0.4562155,0.8885185,""); //(0.17, 0.4, 0.5, 0.53);
                    l->SetTextSize(0.045);
                    // TLegend *leg = new TLegend(0.1797168,0.5390741,0.4562155,0.8885185,"");
                    l->AddEntry("NULL","Herwig"); //"PYTHIA 8 Monash 2013","h");
                    l->AddEntry("NULL","pp, #sqrt{#it{s}} = 13 TeV","h");
                    l->AddEntry("NULL","D^{0} #rightarrow K^{#minus} #pi^{+} and charge conj.","h");
                    l->AddEntry("NULL","in charged jets, anti-#it{k}_{T}, #it{R} = 0.4","h");
                    l->AddEntry("NULL",ptbin,"h");
                    if (ifile != 2) l->AddEntry("NULL",ptD,"h");
                    l->SetTextSize(0.037);
                    l->SetBorderSize(0);
                    // l->Draw("same");


                    // for THnSparse: make clone to work with, make cuts, get projection
                    THnSparse *hsparse_EEC_clone = (THnSparse *) hsparse_EEC->Clone("hsparse_EEC_clone");
                    THnSparse *hsparse_EEC_ptrl_clone = (THnSparse *) hsparse_EEC_ptrl->Clone("hsparse_EEC_ptrl_clone");
                    THnSparse *hsparse_EEC_noweight_clone = (THnSparse *) hsparse_EEC_noweight->Clone("hsparse_EEC_noweight_clone");
                    THnSparse *hsparse_jetinfo_clone = (THnSparse *) hsparse_jetinfo->Clone("hsparse_jetinfo_clone");

                    // get jet pT range
                    hsparse_EEC_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
                    hsparse_EEC_ptrl_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max);
                    hsparse_EEC_noweight_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max);
                    hsparse_jetinfo_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max);

                    if (ifile != 2) { // D0 cuts do not apply to inclusive
                        hsparse_EEC_clone->GetAxis(1)->SetRangeUser(d0_pt_cuts[i], pt_max); // apply cut on Dmeson pt
                        hsparse_EEC_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
                        hsparse_EEC_ptrl_clone->GetAxis(1)->SetRangeUser(d0_pt_cuts[i], pt_max);
                        hsparse_EEC_ptrl_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8);
                        hsparse_EEC_noweight_clone->GetAxis(1)->SetRangeUser(d0_pt_cuts[i], pt_max);
                        hsparse_EEC_noweight_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8);
                        hsparse_jetinfo_clone->GetAxis(1)->SetRangeUser(d0_pt_cuts[i], pt_max);
                        hsparse_jetinfo_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8);
                    }


                    // Project onto observable axis
                    TH1D *hD0_EEC_proj = hsparse_EEC_clone->Projection(4); //CALL THESE TH1*????
                    TH1D *hD0_EEC_ptrl_proj = hsparse_EEC_ptrl_clone->Projection(4);
                    TH1D *hD0_EEC_noweight_proj = hsparse_EEC_noweight_clone->Projection(4);
                    TH1D *hjetpt = hsparse_jetinfo_clone->Projection(0);
                    TH1D *hjetz_pt = hsparse_jetinfo_clone->Projection(3);
                    

                    // Set to appropriate name
                    std::string hname = hD0_EEC_proj->GetName();
                    hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_trkthrd" + trkthrd;
                    hD0_EEC_proj->SetNameTitle(hname.c_str(), hname.c_str());

                    hname = hD0_EEC_ptrl_proj->GetName();
                    hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_trkthrd" + trkthrd;
                    hD0_EEC_ptrl_proj->SetNameTitle(hname.c_str(), hname.c_str());

                    hname = hD0_EEC_noweight_proj->GetName();
                    hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_trkthrd" + trkthrd;
                    hD0_EEC_noweight_proj->SetNameTitle(hname.c_str(), hname.c_str());

                    hname = dname + "_D0_z_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_trkthrd" + trkthrd;
                    hjetz_pt->SetNameTitle(hname.c_str(), hname.c_str());


                    // Rebin or clone
                    int n_obs_bins = 13; //-1;
                    double obs_bins[14] = {0.0001, 0.00020892961308540387, 0.0004365158322401661, 0.0009120108393559096, 
                            0.0019054607179632482, 0.005754399373371567, 0.017378008287493762, 0.03630780547701014, 
                            0.07585775750291836, 0.15848931924611143, 0.3311311214825911, 0.47863009232263853, 0.6918309709189363, 1.0};

                    cout << "About to rebin" << endl;

                    std::string EEC_subname = dname + "_EEC_pt" + std::to_string(pt_min) + '-' + std::to_string(pt_max) + "_R" + jetR + "_trkthrd" + trkthrd;
                    std::string EECptrl_subname = dname + "_EECptrl_pt" + std::to_string(pt_min) + '-' + std::to_string(pt_max) + "_R" + jetR + "_trkthrd" + trkthrd;
                    std::string EECnoweight_subname = dname + "_EECnoweight_pt" + std::to_string(pt_min) + '-' + std::to_string(pt_max) + "_R" + jetR + "_trkthrd" + trkthrd;
                    

                    // TH1* hD0 = (TH1*) hD0_proj->Rebin(n_obs_bins, (hname + "rebin").c_str(), obs_bins);
                    TH1* hD0_EEC = (TH1*) hD0_EEC_proj->Clone( (EEC_subname).c_str() );
                    TH1* hD0_EEC_ptrl = (TH1*) hD0_EEC_ptrl_proj->Clone( (EECptrl_subname).c_str() );
                    TH1* hD0_EEC_noweight = (TH1*) hD0_EEC_noweight_proj->Clone( (EECnoweight_subname).c_str() );
                    cout << "Rebin done" << endl;

                    // Find normalization factor and set normalization
                    double numjets_D0 = hjetpt->Integral();
                    hD0_EEC->Scale(1/numjets_D0, "width");
                    hD0_EEC_ptrl->Scale(1/numjets_D0, "width");
                    hD0_EEC_noweight->Scale(1/numjets_D0, "width");
                    hjetz_pt->Scale(1/numjets_d0, "width");


                    //Format color and style
                    int markercolor1 = marker_color[itrkthrd]; // kMagenta; //D0, me
                    int markerstyle1 = kFullCircle;
                    label1 = "D^{0}-tagged, c-init jets";
                    

                    // Format histograms for plotting (this order needed to keep legend in order and graphs lookin good)
                    // if (plot_case != 3) {
                    //     hD0->GetXaxis()->SetTitle("#it{R}_{L}");
                    // } else {
                    //     hD0->GetXaxis()->SetTitle("p_{T}#it{R}_{L}");
                    // }
                    hD0_EEC->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
                    hD0_EEC_ptrl->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
                    hD0_EEC_noweight->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
                    
                    hD0_EEC_ptrl->GetXaxis()->SetTitle("p_{T}#it{R}_{L}");
                    
                    cout << "about to format D0" << endl;
                    FormatHist(l, hD0_EEC, label1, markercolor1, markerstyle1); //FormatHist(l, hD0, "D^{0}-tagged, c-init jets", kMagenta+3, kOpenSquare);
                    cout << "about to format 1" << endl;
                    FormatHist(ldummy, hD0_EEC_ptrl, label1, markercolor1, markerstyle1);
                    cout << "about to format 2" << endl;
                    FormatHist(ldummy, hD0_EEC_noweight, label1, markercolor1, markerstyle1);
                    // l->AddEntry("NULL","          D* decays off","h");

                    c_EEC->cd();
                    hD0_EEC->Draw("L same");
                    
                    double hD0_EEC_top_binpos = findTopOfCurve(hD0_EEC);
                    drawVertLine(hD0_EEC->GetBinCenter(hD0_EEC_top_binpos), 0, hD0_EEC->GetBinContent(hD0_EEC_top_binpos), markercolor1, 1)->Draw();
                    // vector<double> fullwidth_vec = findWidthOfCurve(hD0,  hD0_top_binpos);
                    // drawHoriLine(fullwidth_vec[1], fullwidth_vec[2], fullwidth_vec[0], kMagenta+3, 1)->Draw();

                    // draw legend
                    l->Draw("same");

                    c_EEC_ptrl->cd();
                    hD0_EEC_ptrl->Draw("L same");
                    l->Draw("same");

                    
                    c_EEC_noweight->cd();
                    hD0_EEC_noweight->Draw("L same");
                    l->Draw("same");
        



                    std::string fname_EEC = outdir + EEC_subname + ".pdf";
                    std::string fname_EECptrl = outdir + EECptrl_subname + ".pdf";
                    std::string fname_EECnoweight = outdir + EECnoweight_subname + ".pdf";
                    const char* fnamec_EEC = fname_EEC.c_str();
                    const char* fnamec_EECptrl = fname_EECptrl.c_str();
                    const char* fnamec_EECnoweight = fname_EECnoweight.c_str();

                    c_EEC->SaveAs(fnamec_EEC);
                    delete c_EEC;
                    c_EEC_ptrl->SaveAs(fnamec_EECptrl);
                    delete c_EEC_ptrl;
                    c_EEC_noweight->SaveAs(fnamec_EECnoweight);
                    delete c_EEC_noweight;

                    // if (pt_min == 10) { //} && grooming == "") {
                        // Write rebinned histograms to root file
                    f_out->cd();
                    hD0_EEC->Write();
                    hD0_EEC_ptrl->Write();
                    hD0_EEC_noweight->Write();
                    hjetz_pt->Write();
                    // }
                    
                } // pT bins loop

                itrkthrd++;
                hD0_pT->Write();
                hD0_z->Write();

            } // end track threshold loop
        
        } // jetR loop

        f->Close();
        delete f;

        ifile++;

    } // file loop

    return;
}