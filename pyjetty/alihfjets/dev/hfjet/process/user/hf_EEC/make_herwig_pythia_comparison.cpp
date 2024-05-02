// ROOT macro to make quark-gluon jet plots
// this makes a comparison between herwig and pythia!!
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

void FormatHist(TLegend *l, TH1 *hist, TString text, int markercolor=1, int markerstyle=8, double markeralpha=1) 
{
    hist->SetLineColor(markercolor);
    hist->SetMarkerColorAlpha(markercolor, markeralpha);
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
int findTopOfCurve(TH1* hist, int checkExtra=1, double startsearch=0.01) {
    
    //for each point, find the slope from the 
    int numbins = hist->GetNbinsX();
    int binstart = hist->FindBin(startsearch);
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


void make_herwig_pythia_comparison() {

//    gROOT->SetBatch(); //prevents plots from showing up
    gStyle->SetOptStat(0);
    SetStyle();
    Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
    Double_t marker_size = 1.5;
    Double_t colors[16] = {kRed, kGreen+2, kBlue, kRed+1, kGreen+1, kBlue+1, kRed+2, kGreen+2, kBlue+2, kRed+3, kGreen+3, kBlue+3, kOrange+1, kViolet+1, kYellow+1, kCyan+1};


    // note: these two files include both weighted and unweighted EECs. herwig includes ptrl, pythia does not.
     //FOR WHEN WEIGHTED/UNWEIGHTED IN SAME FILE
    const char infile_herwig_D0_all[] = "/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alihfjets/dev/hfjet/process/user/hf_EEC/herwig/AnalysisResultsFinal_herwig.root";
    const char infile_pythia_D0[] = "/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alihfjets/dev/hfjet/process/user/hf_EEC/plots/final/AnalysisResultsFinalcomparison_all.root";
    const char infile_pythia_D0wDstar[] = "/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alihfjets/dev/hfjet/process/user/hf_EEC/plots/final/AnalysisResultsFinal_afteranalysis_Dstar_plotcase3.root";

    //plot cases~
    // 0 = herwig D0 vs pythia D0
    // 1 = herwig D0 vs herwig D0wDstar vs pythia D0 vs pythia D0wDstar
    
    //CONTOL VARIABLES HERE
    int plot_case = 1;
    bool logstring = false;

    TString label1 = "";

    TFile* f_herwig = new TFile(infile_herwig_D0_all, "READ"); 
    TFile* f_pythia_D0 = new TFile(infile_pythia_D0, "READ");
    TFile* f_pythia_D0wDstar = new TFile(infile_pythia_D0wDstar, "READ");


    // Output directory
    std::string outdir = "plots/final/herwig-pythia/";
    // Output file for binned results
    std::string outfile;
    if (plot_case == 0) {
        outfile = outdir + "AnalysisResultsFinal_herwig_pythia_comparison.root"; 
    } else if (plot_case == 1) {
        outfile = outdir + "AnalysisResultsFinal_herwig_pythia_D0wDstar_comparison.root"; 
    }
    TFile* f_out = new TFile(outfile.c_str(), "RECREATE");
    



    const int pt_bins[] = { 7, 10, 15, 30 }; //{ 10, 20, 40 }; // CHANGE HERE!!
    const int n_bins = 3; 
    for (int i = 0; i < n_bins; i++) {
        cout << "in pt bin" << i << endl;
        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];

        // define pt related variables
        TString ptbin = TString::Format("%d #leq #it{p}_{T}^{ch. jet} < %d GeV/#it{c}, #font[122]{|}#it{#eta}_{jet}#font[122]{|} #leq 0.5", pt_min, pt_max);
        TString ptD = TString::Format("5 #leq #it{p}_{T}^{D^{0}} < %d GeV/#it{c}, #font[122]{|}#it{y}_{D^{0}}#font[122]{|} #leq 0.8", pt_max);

        // Names of histograms in the file
        const std::string pythia_Djet_EEC_name = "hsparsejet_c_clone_proj_3_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
        const std::string pythia_DjetwithDstar_EEC_name = "hsparsejet_c_Dstar_clone_proj_3_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
        // const std::string pythia_Djet_EEC_noweight_name = "hsparsejet_c_uw_clone_proj_3_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
        
        // const std::string herwig_Djet_EEC_name = "hsparse_EEC_clone_proj_3_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_trkthrd1.0";
        // const std::string herwig_DjetwithDstar_EEC_name = "hsparse_EEC_clone_proj_3_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_trkthrd1.0";
        // // const std::string herwig_EEC_ptrl_name = "hsparse_EEC_ptrl_clone_proj_3_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_trkthrd1.0";
        // const std::string herwig_EEC_noweight_name = "hsparse_EEC_noweight_clone_proj_3_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_trkthrd1.0";

        const std::string herwig_Djet_EEC_name = "D0_noDstar_EEC_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_R0.4_trkthrd1.0";
        const std::string herwig_DjetwithDstar_EEC_name = "D0_EEC_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_R0.4_trkthrd1.0";
        
        
        // make a canvas for each pt range
        TCanvas* c = new TCanvas();
        ProcessCanvas(c);
        c->cd();
        gPad->SetLogx();
        if (logstring) {
            gPad->SetLogy();
        }
        

        TLegend* l; // = new TLegend(0.17, 0.65, 0.5, 0.85);
        TLegend* l2;
        TLegend* l3;

        double maxy = 0;


        // Open histograms


        l = new TLegend(0.1797168,0.400741,0.4562155,0.8885185,""); //(0.17, 0.4, 0.5, 0.53);
        l->SetTextSize(0.045);
        // TLegend *leg = new TLegend(0.1797168,0.5390741,0.4562155,0.8885185,"");
        l->AddEntry("NULL","PYTHIA 8 Monash 2013 and Herwig","h");
        l->AddEntry("NULL","pp, #sqrt{#it{s}} = 13 TeV","h");
        l->AddEntry("NULL","D^{0} #rightarrow K^{#minus} #pi^{+} and charge conj.","h");
        l->AddEntry("NULL","in charged jets, anti-#it{k}_{T}, #it{R} = 0.4","h");
        l->AddEntry("NULL",ptbin,"h");
        l->AddEntry("NULL",ptD,"h");
        l->SetTextSize(0.037);
        l->SetBorderSize(0);
        // l->Draw("same");

        TLegend* l_fake = new TLegend(0.,0.,0.1,0.1,"");


        TH1D* hD0_pythia = (TH1D*) f_pythia_D0->Get(pythia_Djet_EEC_name.c_str());
        TH1D* hD0wDstar_pythia = (TH1D*) f_pythia_D0wDstar->Get(pythia_DjetwithDstar_EEC_name.c_str());
        TH1D* hD0_herwig = (TH1D*) f_herwig->Get(herwig_Djet_EEC_name.c_str());
        TH1D* hD0wDstar_herwig = (TH1D*) f_herwig->Get(herwig_DjetwithDstar_EEC_name.c_str());
   


        //Format color and style
        int markercolor1 = kCyan-3; //pythia
        int markerstyle1 = kFullCircle;
        int markercolor2 = kCyan+2;//kAzure-5; //pythia D0 with Dstar
        int markerstyle2 = kOpenCircle;
        int markercolor3 = kAzure-4;//kBlue-4; //herwig
        int markerstyle3 = kFullDiamond;
        int markercolor4 = kAzure-5; //herwig D0 with Dstar
        int markerstyle4 = kOpenDiamond;
        label1 = "PYTHIA primordial D0-tagged"; //, c-init jets";
        TString label2 = "PYTHIA all D0-tagged"; //, c-init jets";
        TString label3 = "Herwig primordial D0-tagged"; //, c-init jets";
        TString label4 = "Herwig all D0-tagged"; //, c-init jets";

        

        // Format histograms for plotting (this order needed to keep legend in order and graphs lookin good)
        hD0_pythia->GetXaxis()->SetTitle("#it{R}_{L}");
        hD0_pythia->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
        hD0wDstar_pythia->GetXaxis()->SetTitle("#it{R}_{L}");
        hD0wDstar_pythia->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
        hD0_herwig->GetXaxis()->SetTitle("#it{R}_{L}");
        hD0_herwig->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
        hD0wDstar_herwig->GetXaxis()->SetTitle("#it{R}_{L}");
        hD0wDstar_herwig->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");

        cout << "about to format D0" << endl;
        
        
        FormatHist(l, hD0_pythia, label1, markercolor1, markerstyle1, 0.80); //FormatHist(l, hD0, "D^{0}-tagged, c-init jets", kMagenta+3, kOpenSquare);
        FormatHist(l, hD0_herwig, label3, markercolor3, markerstyle3, 0.80);
        if (plot_case == 1) {
            FormatHist(l, hD0wDstar_pythia, label2, markercolor2, markerstyle2, 0.80);
            FormatHist(l, hD0wDstar_herwig, label4, markercolor4, markerstyle4, 0.80);
        }

         
        // get maximum
        hD0_pythia->SetMaximum(hD0_pythia->GetMaximum()*1.5);
        
        hD0_pythia->Draw("L same");
        hD0_herwig->Draw("L same");
        if (plot_case == 1) {
            hD0wDstar_pythia->Draw("L same");
            hD0wDstar_herwig->Draw("L same");
        }
        

        l->Draw("same");
        
        
        double hD0_pythia_top_binpos = findTopOfCurve(hD0_pythia, 1, 0.08);
        drawVertLine(hD0_pythia->GetBinCenter(hD0_pythia_top_binpos), 0, hD0_pythia->GetBinContent(hD0_pythia_top_binpos), markercolor1, 1)->Draw();
        double hD0_herwig_top_binpos = findTopOfCurve(hD0_herwig, 1, 0.08);
        drawVertLine(hD0_herwig->GetBinCenter(hD0_herwig_top_binpos), 0, hD0_herwig->GetBinContent(hD0_herwig_top_binpos), markercolor3, 1)->Draw();
        if (plot_case == 1) {
            double hD0wDstar_pythia_top_binpos = findTopOfCurve(hD0wDstar_pythia, 1, 0.08);
            drawVertLine(hD0wDstar_pythia->GetBinCenter(hD0wDstar_pythia_top_binpos), 0, hD0wDstar_pythia->GetBinContent(hD0wDstar_pythia_top_binpos), markercolor2, 1)->Draw();
            double hD0wDstar_herwig_top_binpos = findTopOfCurve(hD0wDstar_herwig, 1, 0.08);
            drawVertLine(hD0wDstar_herwig->GetBinCenter(hD0wDstar_herwig_top_binpos), 0, hD0wDstar_herwig->GetBinContent(hD0wDstar_herwig_top_binpos), markercolor4, 1)->Draw();
        
            // plot line for D* too
            // if (pt_min == 15)
            hD0wDstar_pythia_top_binpos = findTopOfCurve(hD0wDstar_pythia);
            drawVertLine(hD0wDstar_pythia->GetBinCenter(hD0wDstar_pythia_top_binpos), 0, hD0wDstar_pythia->GetBinContent(hD0wDstar_pythia_top_binpos), markercolor2, 1)->Draw();
            hD0wDstar_herwig_top_binpos = findTopOfCurve(hD0wDstar_herwig);
            drawVertLine(hD0wDstar_herwig->GetBinCenter(hD0wDstar_herwig_top_binpos), 0, hD0wDstar_herwig->GetBinContent(hD0wDstar_herwig_top_binpos), markercolor4, 1)->Draw();
        
        }

        // vector<double> fullwidth_vec = findWidthOfCurve(hD0,  hD0_top_binpos);
        // drawHoriLine(fullwidth_vec[1], fullwidth_vec[2], fullwidth_vec[0], kMagenta+3, 1)->Draw();
        



        std::string fname = outdir + "herwig-pythia_comparison_pt" + std::to_string(pt_min) + '-' + std::to_string(pt_max) + "_R0.4_plotcase" + std::to_string(plot_case) + ".pdf";
        const char* fnamec = fname.c_str();
        c->SaveAs(fnamec);
        delete c;


        f_out->cd();
        hD0_pythia->Write();
        hD0wDstar_pythia->Write();
        hD0_herwig->Write();
        hD0wDstar_herwig->Write();
            
    } // pT bins loop

    f_pythia_D0->Close();
    delete f_pythia_D0;
    f_pythia_D0wDstar->Close();
    delete f_pythia_D0wDstar;
    f_herwig->Close();
    delete f_herwig;

    return;
}
