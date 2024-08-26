// ROOT macro to take correlation tuples and turn them into histograms
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

void ProcessCanvas(TCanvas *Canvas, bool setlogx=false, bool setlogy=false) { 
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

    if (setlogx) gPad->SetLogx();
    if (setlogy) gPad->SetLogy();
}

void FormatHist(TLegend *l, TH1 *hist, TString text, int markercolor=1, int markerstyle=8, double linealpha=1., bool drawline=false) 
{
    if (drawline) {
        // for (int k=0; k < hist->GetNbinsX();k++){
        //     hist->SetBinError(k+1, 0);
        // }
        hist->SetMarkerStyle(20);
        hist->SetMarkerColorAlpha(markercolor, 0);

        hist->SetFillStyle(0);
        hist->SetLineColorAlpha(markercolor, linealpha);
        hist->SetFillColor(markercolor);
        hist->SetLineStyle(1);
        hist->SetLineWidth(3);
    } else {
        hist->SetLineColor(markercolor);
        hist->SetMarkerColor(markercolor);
        hist->SetMarkerStyle(markerstyle);
        hist->SetMarkerSize(1.5);
    }
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

// ptrl is a boolean that says whether ptRL is being plotted (instead of RL)
// --> controls where the cutoff is to not look for the max point
//checkExtra is how many point-to-point slopes after finding a decreasing slope I want to check
int findTopOfCurve(TH1* hist, bool ptrl, int checkExtra=1) {
    
    // define the x axis start of bin search
    double startbinsearchat = ptrl ? 0.1 : 0.01;

    //for each point, find the slope from the 
    int numbins = hist->GetNbinsX();
    int binstart = hist->FindBin(startbinsearchat);
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


// ======================================================= //
//                     SOME FUNCTIONS 
// ======================================================= //

//get histogram and clone it
THnSparse * getHistAndClone(TFile *f, std::string histname) {
    THnSparse *hsparsejet = (THnSparse*) f->Get(histname.c_str());
    std::string hn = hsparsejet->GetName();
    hn += "_clone";
    THnSparse *hsparsejet_clone = (THnSparse *) hsparsejet->Clone( hn.c_str() ); //TODO: change this name!

    return hsparsejet_clone;
    
}


void applyCuts(THnSparse *hsparse, int pt_min, int pt_max, int RLaxis, double RL_min, double RL_max, 
               bool EWaxis = false) { //todo: don't need paxis anymore??

    hsparse->GetAxis(0)->SetRangeUser(pt_min, pt_max);
    // cout << "RL AXIS is " << RLaxis << " w RL min: " << RL_min << " & RL max: " << RL_max << endl;
    if (RLaxis != -1) {
        cout << " here! RL AXIS is " << RLaxis << " w RL min: " << RL_min << " & RL max: " << RL_max << endl;
        hsparse->GetAxis(RLaxis)->SetRangeUser(RL_min, RL_max);
        // hsparse->GetAxis(RLaxis)->SetRangeUser(1e-5, 1.);
        
        //hsparse->GetAxis(RLaxis)->SetRangeUser(-2, -1); //9.9e-5, 1e-4);//RL_min, RL_max);
        // cout << "UNDERFLOW BIN HAS " << hsparse->Projection(RLaxis)->GetBinContent(0) << " entries" << endl;
        // TH1D* testhist = hsparse->Projection(RLaxis);
        // for (int aa=0; aa<testhist->GetNbinsX(); aa++) {
        //     cout << testhist->GetBinLowEdge(aa) << " ";
        // }
        // cout << endl;
        
    }
    if (EWaxis) { //energy weight axis
        hsparse->GetAxis(4)->SetRangeUser(0., 0.3);
    }

    // cout << "CHECK 1A LOOKING AT NUM EMNRTIRES HERE: " << hsparse->GetEntries() << endl;

    
}

// get the observable histogram
//usually obsaxis is 3, but in new histograms it is 4. For jet level histograms, use 0.
TH1D * getObsHist(TFile *filename, std::string h_name, std::string h_jet_name, int pt_min, int pt_max, 
                  int RLaxis, double RL_min, double RL_max, std::string newhistname, int obsaxis, std::string xtitle,
                  int n_rebin_bins, double rebinbins[], int normalized=0, bool EWaxis=false) {
    
    cout << "HNAME is " << h_name << endl;
    cout << "RL AXIS is " << RLaxis << " w RL min: " << RL_min << " & RL max: " << RL_max << endl;
    
    THnSparse *hsparse = getHistAndClone(filename, h_name);
    THnSparse *hsparse_jetlevel = getHistAndClone(filename, h_jet_name);

    // cout << "CHECK 1 LOOKING AT NUM EMNRTIRES HERE: " << hsparse->GetEntries() << endl;
    // cout << "CHECK 1 LOOKING AT NUM ENTRIES: " << hsparse_jetlevel->GetEntries() << endl;

    // if (debug)
    // TH1D *testhist = hsparse->Projection(obsaxis);
    // cout << "underflow entries: " << testhist->GetBinContent(0) << " & overflow entries:" << testhist->GetBinContent(testhist->GetNbinsX()+1) << endl;
    // int sum = 0;
    // for (int aa=0; aa <= testhist->GetNbinsX()+1; aa++) {
    //     sum += testhist->GetBinContent(aa);
    // }
    // cout << " and the total sum is " << sum << endl;
    // cout << " and the integral is " << testhist->Integral() << " & " << testhist->Integral("width") << endl;

    applyCuts(hsparse, pt_min, pt_max, RLaxis, RL_min, RL_max, EWaxis);
    applyCuts(hsparse_jetlevel, pt_min, pt_max, -1, RL_min, RL_max, false); //making no cuts on RL to keep it jet level??

    
    TH1D *h_proj = hsparse->Projection(obsaxis);
    TH1D *h_proj_jetlevel = hsparse_jetlevel->Projection(0); // jet pt axis

    // cout << "CHECK 3 LOOKING AT NUM EMNRTIRES HERE: " << h_proj->GetEntries() << endl;
    // cout << "CHECK 3 LOOKING AT NUM ENTRIES: " << h_proj_jetlevel->GetEntries() << endl;
    

    std::string hname = h_proj->GetName();
    hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_RL" + std::to_string(RL_min) + "-" + std::to_string(RL_max);
    h_proj->SetNameTitle(hname.c_str(), hname.c_str());
    // cout << "HISTOGRAM NAME IS " << h_proj->GetName() << " AND TITLE " << h_proj->GetTitle() << endl;

    // allow rebin or cloning here
    TH1D* hist;
    // int n_obs_bins = sizeof(rebinbins) / sizeof(rebinbins[0]) - 1;
    // if ( n_rebin_bins != 0) {
    //     hist = (TH1D*) h_proj->Rebin(n_rebin_bins, newhistname.c_str(), rebinbins); //h_proj->GetName()
    // } else {
    //     hist = (TH1D*) h_proj->Clone(newhistname.c_str());
    // }
    hist = (TH1D*) h_proj->Clone(newhistname.c_str());
    // if (debug)
    cout << "LOOKING AT NUM BINS: " << hist->GetNbinsX() << endl;
    cout << "LOOKING AT NUM ENTRIES: " << hist->GetEntries() << endl;

    // for (int i=0; i<hist->GetNbinsX(); i++) {
    //     cout << hist->GetBinContent(i) << " ";
    // }
    // cout << endl;

    // normalize
    cout << "IS THIS NORMALIZED? " << normalized << endl;
    if (normalized == 1) {
        double numjets = hist->Integral();
        hist->Scale(1/numjets, "width");
        cout << "IN NORMALIZED1" << endl;
    } else if (normalized == 2) {
        double numjets = h_proj_jetlevel->Integral();
        hist->Scale(1/numjets, "width");
        cout << "IN NORMALIZED2" << endl;
    }
    // cout << "There are " << h_proj->GetEntries() << " pair entries in this pt bin" << endl;
    // cout << "There are " << h_proj_jetlevel->GetEntries() << " jet entries in this pt bin" << endl;

    // hist->GetXaxis()->SetTitle(xtitle.c_str()); //("#it{R}_{L}");
    // // THESE AXES LABELS ARE DEF WRONG!!
    // if (normalized != 0) hist->GetYaxis()->SetTitle( Form("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d%s}", xtitle.c_str()) );
    // else hist->GetYaxis()->SetTitle( Form("#frac{d#it{N}_{EEC}}{d%s}", xtitle.c_str()) ); 


    delete hsparse;
    delete hsparse_jetlevel;
    delete h_proj_jetlevel;
    delete h_proj;

    return hist;

}

// ======================================================= //
//                   SECONDARY FUNCTIONS 
// ======================================================= //


// open tchain - gather tchain of all TTrees of one type in a pt-hat bin
TChain * make_tchain() {

    // TChain mychain("mychain");
    TChain *mychain = new TChain("mychain");
    
    // tree1->Add(filename.c_str());
    
    return mychain;
    
}

// Add items to tchain
void add_to_tchain(TChain* mychain, std::string indir, int ibin, int jbin, std::string fname, std::string tuple_name) {
    std::string filename = indir + "/" + std::to_string(ibin) + "/" + std::to_string(jbin) + "/" + fname;
    // std::string filename = "/global/cfs/cdirs/alice/blianggi/mypyjetty/analysis/testing/AnalysisResultsTesting.root";
    mychain->Add( Form("%s/%s", filename.c_str(), tuple_name.c_str()) );
    // cout << "chain entries: " << mychain->GetEntries() << endl;
        
    // cout << "FINAL chain entries: " << mychain->GetEntries() << endl;
}

// function to make a hist of one observable in one pt and RL bin
TH1D* draw_tchain( TChain *mychain, int pt_min, int pt_max, double RL_min, double RL_max, bool weighted) {

    
    //Also make RL and otehr cuts?? maybe do this in another function??
    if (weighted) { //format is "weight *(boolean expression)"
        mychain->Draw( "obs>>htemp", Form("weights*(jet_pt > %d && jet_pt < %d && RL > %f && RL < %f)", pt_min, pt_max, RL_min, RL_max), "e" );
    } else {
        mychain->Draw( "obs>>htemp", Form("jet_pt > %d && jet_pt < %d && RL > %f && RL < %f", pt_min, pt_max, RL_min, RL_max), "e" );
    }

    // TTree *cutTree = mychain->CopyTree( Form("jet_pt > %d && jet_pt < %d && RL > %d && RL < %d", pt_min, pt_max, RL_min, RL_max) );

    TH1D * htemp = (TH1D*) gDirectory->Get("htemp");

    return htemp;
}

void save_histogram(TFile * fout, TH1D * hist, TCanvas * c_temp, std::string ptbin_name, std::string weightstr,
                    std::string ptname, std::string RLname, std::string add_name="") {
    fout->cd();
    hist->Write();

    c_temp->cd();
    hist->Draw();

    std::string outdir = "plots/thirdattempt/" + ptbin_name + "/";//"plots/test/";
    std::string fname_out = outdir + "individuals/corrhist_deltap" + weightstr + "_pt" + ptname + RLname + "_R0.4" + add_name + ".pdf";
    c_temp->SaveAs(fname_out.c_str());
    delete c_temp;


    delete hist;

}

// function to analyze ONE observable over all pt-hat, pt, and RL bins
// need to specify 2nd and onward axes sizes for a 2+D array
void analyze_observable( std::string indir, std::string fname, int n_bins, const int pt_bins[], 
                         int n_RLbins, const double RL_bins[][8], std::string observable_name,
                         TFile * f_out, std::string add_name, bool debug=false ) {

    // set up variables
    std::string pt_bin_names[n_bins]; // const std::string pt_bin_names[] = {"20-40", "40-60", "60-80"};
    for (int ipt = 0; ipt < n_bins; ipt++) {
        pt_bin_names[ipt] = std::to_string(pt_bins[ipt]) + "-" + std::to_string(pt_bins[ipt+1]);
    }

    std::string observable_tuple_name = Form("h_corr_%s_JetPt_Truth_R0.4_1.0", observable_name.c_str());
    std::string jetRname = "_R0.4"; // + jetR;
    std::string thrname = "_t1.0"; // + threshold;
    
    // start looping through all the files
    for (int ibin=1; ibin<=20; ibin++) {
        if (debug) cout << " in ibin " << ibin << endl;
        TChain *obs_chain = make_tchain();

        for (int jbin=1; jbin<=100; jbin++) {
            if (debug) cout << " in jbin " << jbin << endl;
            add_to_tchain(obs_chain, indir, ibin, jbin, fname, observable_tuple_name);

        }

        for (int i = 0; i < n_bins; i++) {
            int pt_min = pt_bins[i];
            int pt_max = pt_bins[i+1];
            std::string ptname = std::to_string(pt_min) + '-' + std::to_string(pt_max);
            if (debug) cout << " in pt bin " << i << endl;

            
            // do stuff to get one plot across all RL bins

            for (int j = 0; j < n_RLbins; j++) {
                double RL_min = RL_bins[i][j];
                double RL_max = RL_bins[i][j+1];
                std::string RLname = Form("_RL%.3f-%.3f", RL_min, RL_max); 
                std::string hist_addname = jetRname + thrname + ptname + RLname;
                if (debug) cout << " in RL bin" << j << " with " << RL_min << " - " << RL_max << endl;
                

                // // function_to_draw_tchain() - do weighted first, then unweighted
                // TH1D * htemp_weighted = draw_tchain( obs_chain, pt_min, pt_max, RL_min, RL_max, true );
                TH1D * htemp_unweighted = draw_tchain( obs_chain, pt_min, pt_max, RL_min, RL_max, false );

				// format_hist()
                TCanvas* c_temp_weighted = new TCanvas();
                ProcessCanvas(c_temp_weighted, false, true); // make these variables too?
                TCanvas* c_temp_unweighted = new TCanvas();
                ProcessCanvas(c_temp_unweighted, false, true); // make these variables too?


                // // save histograms to root file and as pdf and delete
				// save_histogram(f_out, htemp_weighted, c_temp_weighted, pt_bin_names[i], "_Weighted", ptname, RLname, add_name);
                // save_histogram(f_out, htemp_unweighted, c_temp_unweighted, pt_bin_names[i], "", ptname, RLname, add_name);

            }
        }

        //delete_tchain()
    }
}

    

// ======================================================= //
//                     MAIN FUNCTION 
// ======================================================= //

void make_histograms_from_tuples() {

    gStyle->SetOptStat(0);
    SetStyle();
    

    // Files
    // const char infile[] = "/software/users/blianggi/mypyjetty/analysis/output/100k/AnalysisResultsFinal.root"; //hiccup
    const std::string indir = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/pythiagen/histograms/29160880/ntuples/26652369"; //perlmutter, after june 2024
    const std::string infilename = "AnalysisResults.root";
    // const char infile[] = "/Volumes/NO NAME/AnalysisResultsFinal.root"; //local

    // ntuple names
    std::string jetinfo_name = "h_JETINFOjet_pt_JetPt_Truth_R0.4_1.0";
    std::string baryon_name = "h_baryon_JetPt_Truth_R0.4_1.0";
    std::string meson_name = "h_meson_JetPt_Truth_R0.4_1.0";
    std::string baryonmeson_name = "h_corr_baryonmeson_JetPt_Truth_R0.4_1.0";
    std::string deltap_name = "h_corr_deltap_JetPt_Truth_R0.4_1.0";
    std::string deltapt_name = "h_corr_deltapt_JetPt_Truth_R0.4_1.0";
    std::string energyweights_name = "h_corr_energyweights_JetPt_Truth_R0.4_1.0";
    std::string oppcharge_name = "h_corr_oppcharge_JetPt_Truth_R0.4_1.0";
    std::string samecharge_name = "h_corr_samecharge_JetPt_Truth_R0.4_1.0";

    std::string jet1D_name = "h_1Djet_pt_JetPt_Truth_R0.4_1.0;10"; // this one is a histogram

    std::string deltap_str = "deltap";


    const int pt_bins[] = { 20, 40, 60, 80 };
    const double RL_bins[3][8] = { { 0, 1e-2, 3e-2, 7e-2, 1.5e-1, 3e-1, 4e-1 }, 
                            { 0, 1e-2, 2.5e-2, 4e-2, 8e-2, 2.5e-1, 4e-1 }, 
                            { 0, 1e-2, 2.5e-2, 3e-2, 4.5e-2, 2e-1, 4e-1 } };
    const int n_bins = sizeof(pt_bins) / sizeof(pt_bins[0]) - 1; //3;
    const int n_RLbins = sizeof(RL_bins[0]) / sizeof(RL_bins[0][0]) - 1; //gets the columns //6; //7; //5;


    // Output file for binned results
    std::string root_outfile = "plots/thirdattempt/AnalysisResultsFinal.root"; // + add_name + ".root";
    TFile* f_out = new TFile(root_outfile.c_str(), "RECREATE");
    std::string add_name = "_othercorrel";



    // start analyzing different observables
    analyze_observable( indir, infilename, n_bins, pt_bins, n_RLbins, RL_bins, deltap_str, f_out, add_name );
    // analyze_observable( indir, infilename, pt_bins, RL_bins, deltapt_name, f_out);



    

    // ---------------- end ----------------



        	// tn.Draw("x>>htempx", "x>2 && x < 7");
            // TH1 * htempx = (TH1*) gDirectory->Get("htempx");
            // cout << "Num entries in htempx " << htempx->GetEntries() << endl;

            // TH2D * hist = new TH2D("hidt", "hidt", 10, 0, 10, 10, 0, 10);
            // tn.Draw("y:x>>hidt", "x>2 && x < 7", "e");





    TFile* f;

    //=============================================//
    
    // get the path for all the files
    // loop over each file
    // make histograms from each tuple, but keep it within the pt bin?
    // save into another root file?
        // only one jet R, only one track threshold
        // divide into pt or RL bins until after scaling
        // plots:
            // delta p weighted, delta p unweighted
            // delta pt weighted, delta pt unweighted
            // opp charge vs same charge quantities (weighted/unweighted??)
            // --> later, opp charge/same charge as a function of RL
            // energy weights weighted, energy weights unweighted
            // unweighted RL?? - meh
            // baryon/meson
            // later, charged EECs broken up into RL regions



    // separately, then scale all the pt hat bins
    // then merge them
    // then plot



    //=============================================//


    // CONTOL VARIABLES HERE
    int normed = 1; // 0 if unnormalized, 1 for self-normalization 
    bool weighted = false; // true if using "Weighted", false if using unweighted
    

    std::string norm_string = "";
    if (normed == 0) {
        norm_string = "_unnormalized";
    } else if (normed == 1) {
        norm_string = "_selfnorm";
    }

    // f = new TFile(infile, "READ");
    // add_name = "_othercorrel" + norm_string;
    // cout << "output name will be " << add_name << endl;

    // // plot unweighted hists
    // plot_histograms(f, add_name, normed, weighted);

    // // plot weighted hists
    // plot_histograms(f, add_name, normed, true);



    // f->Close();
    // delete f;


    return;
}
