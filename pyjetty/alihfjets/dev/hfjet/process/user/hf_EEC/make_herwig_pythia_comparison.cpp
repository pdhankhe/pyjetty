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

void FormatHistwithLine(TLegend *l, TH1 *hist, TString text, int linecolor=1, int linestyle=1, double linealpha=1) 
{
    hist->SetMarkerStyle(20);
    hist->SetMarkerColorAlpha(linecolor, 0);

    hist->SetFillStyle(0);
    hist->SetLineColorAlpha(linecolor, linealpha);
    hist->SetFillColor(linecolor);
    // hist->SetLineStyle(linestyle);
    hist->SetLineWidth(3);
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
    const char infile_sherpa_D0[] = "/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alihfjets/dev/hfjet/process/user/hf_EEC/sherpa/AnalysisResultsFinal_sherpa.root";
    
    //plot cases~
    // 0 = herwig D0 vs pythia D0
    // 1 = herwig D0 vs herwig D0wDstar vs pythia D0 vs pythia D0wDstar
    // 2 = adding in sherpa D0 vs herwig D0 vs pythia D0
    
    //CONTOL VARIABLES HERE
    int plot_case = 0;
    bool logstring = false;

    TString label1 = "";

    TFile* f_herwig = new TFile(infile_herwig_D0_all, "READ"); 
    TFile* f_pythia_D0 = new TFile(infile_pythia_D0, "READ");
    TFile* f_pythia_D0wDstar = new TFile(infile_pythia_D0wDstar, "READ");
    TFile* f_sherpa = new TFile(infile_sherpa_D0, "READ");


    // Output directory
    std::string outdir = "plots/final/herwig-pythia/";
    if (plot_case == 2) {
        outdir = "plots/final/herwig-pythia-sherpa/";
    }
    // Output file for binned results
    std::string outfile;
    if (plot_case == 0) {
        outfile = outdir + "AnalysisResultsFinal_herwig_pythia_comparison.root"; 
    } else if (plot_case == 1) {
        outfile = outdir + "AnalysisResultsFinal_herwig_pythia_D0wDstar_comparison.root"; 
    } else if (plot_case == 2) {
        outfile = outdir + "AnalysisResultsFinal_herwig_pythia_sherpa_comparison.root"; 
    }
    TFile* f_out = new TFile(outfile.c_str(), "RECREATE");

    //--------------------------------------------------------//
    
    //Format color and style
    int markercolor1 = kCyan-3; //pythia
    int markerstyle1 = kFullCircle;
    int markercolor2 = kCyan+2;//kAzure-5; //pythia D0 with Dstar
    int markerstyle2 = kOpenCircle;
    int markercolor3 = kAzure-4;//kBlue-4; //herwig
    int markerstyle3 = kFullDiamond;
    int markercolor4 = kAzure-5; //herwig D0 with Dstar
    int markerstyle4 = kOpenDiamond;

    int markercolor5 = kTeal-5; //sherpa D0 c_jetpt10
    int markerstyle5 = kFullStar;
    int markercolor6 = kViolet+6; //sherpa D0 c_jetpt10_lund
    int markerstyle6 = kFullStar;
    int markercolor7 = kTeal+5; //sherpa D0 c_jetpt15
    int markerstyle7 = kFullStar;
    int markercolor8 = kViolet+8; //sherpa D0 c_jetpt15_lund
    int markerstyle8 = kFullStar;
    
    label1 = "PYTHIA primordial D0-tagged"; //, c-init jets";
    TString label2 = "PYTHIA all D0-tagged"; //, c-init jets";
    TString label3 = "Herwig primordial D0-tagged";
    TString label4 = "Herwig all D0-tagged";
    TString label5 = "Sherpa c_jetpt10"; //primordial D0-tagged";
    TString label6 = "Sherpa c_jetpt10_lund"; //primordial D0-tagged";
    TString label7 = "Sherpa c_jetpt15"; //primordial D0-tagged";
    TString label8 = "Sherpa c_jetpt15_lund"; //primordial D0-tagged";

    
    //--------------------------------------------------------//
    //look at pt bins
    const std::string herwig_Djet_pt_name = "D0_noDstar_trkthrd1.0_D0_pt";
    const std::string herwig_DjetwithDstar_pt_name = "D0_trkthrd1.0_D0_pt";
    const std::string herwig_Djet_z_name = "D0_noDstar_trkthrd1.0_D0_z";

    const std::string pythia_Djet_pt_name = "hD0_pt";
    // const std::string pythia_Djet_uw_pt_name = "hD0_uw_pt";
    const std::string pythia_Djet_z_name = "hD0z";

            
    TCanvas* c_pt = new TCanvas();
    ProcessCanvas(c_pt);
    c_pt->cd();
    gPad->SetLogy();
    gPad->SetBottomMargin(0.31);

    TPad *pad1 = new TPad("pad1","pad1",0.,0.,1.,1.);
    // pad1->SetLogy();
    // pad1->SetBottomMargin(0.5); //31);
    // pad1->Draw();

    TH1D* hD0_pt_pythia = (TH1D*) f_pythia_D0->Get(pythia_Djet_pt_name.c_str());
    // TH1D* hD0wDstar_pt_pythia = (TH1D*) f_pythia_D0wDstar->Get(pythia_Djet_pt_name.c_str());
    TH1D* hD0_pt_herwig = (TH1D*) f_herwig->Get(herwig_Djet_pt_name.c_str());
    TH1D* hD0wDstar_pt_herwig = (TH1D*) f_herwig->Get(herwig_DjetwithDstar_pt_name.c_str());

    TH1D* hD0z_pythia = (TH1D*) f_pythia_D0->Get(pythia_Djet_z_name.c_str());
    TH1D* hD0z_herwig = (TH1D*) f_herwig->Get(herwig_Djet_z_name.c_str());

    hD0_pt_pythia->GetXaxis()->SetTitle("D^{0} p_{T}");
    // hD0wDstar_pt_pythia->GetXaxis()->SetTitle("D^{0} p_{T}");
    hD0_pt_herwig->GetXaxis()->SetTitle("D^{0} p_{T}");
    hD0wDstar_pt_herwig->GetXaxis()->SetTitle("D^{0} p_{T}");

    hD0z_pythia->GetXaxis()->SetTitle("D^{0} z");
    hD0z_herwig->GetXaxis()->SetTitle("D^{0} z");

    TLegend* l_fake = new TLegend(0.5097168,0.670741,0.7362155,0.8885185,""); //not so fake anymore but whatever
    l_fake->AddEntry("NULL","in charged jets, anti-#it{k}_{T}, #it{R} = 0.4","h");
    l_fake->SetTextSize(0.037);
    l_fake->SetBorderSize(0);
    TLegend* l2 = new TLegend(0.,0.,0.1,0.1,""); //this one is fake

    FormatHist(l_fake, hD0_pt_pythia, label1, markercolor1, markerstyle1, 0.80);
    // FormatHist(l_fake, hD0wDstar_pt_pythia, label2, markercolor2, markerstyle2, 0.80);
    FormatHist(l_fake, hD0_pt_herwig, label3, markercolor3, markerstyle3, 0.80);
    // FormatHist(l_fake, hD0wDstar_pt_herwig, label4, markercolor4, markerstyle4, 0.80);

    hD0_pt_pythia->Draw("same");
    // hD0wDstar_pt_pythia->Draw("L same"); //not implemented properly yet
    hD0_pt_herwig->Draw("same");
    // hD0wDstar_pt_herwig->Draw("L same");

    //PT STUFF FIGURE OUT WHERE TO PUT THIS LATER
    // c_pt->cd();
    // pad1->cd();
    l_fake->Draw("same");
    hD0_pt_pythia->GetXaxis()->SetLabelSize(0);

    TPad *pad2 = new TPad("pad1","",0.,0.,1.,1.);
    pad2->SetTopMargin(0.71);
    pad2->SetFillColor(0);
    pad2->SetFillStyle(0);
    pad2->Draw();
    // pad2->SetLogx();
    pad2->SetLogy();
    pad2->cd();

    std::string ratio_name = "hratio_D0_pt";
    TH1D* hratio = (TH1D*) hD0_pt_herwig->Clone(ratio_name.c_str());
    hratio->Divide(hD0_pt_pythia);
    // hratio->SetMinimum(0.5);
    hratio->SetMaximum(10.);


    FormatHist(l2, hratio, "ratio", kBlack, markers[2]); //, 0.05, 0.04, 1.2, 0.035, 0.03, 1.5);
    hratio->GetYaxis()->SetTitle("#frac{Herwig}{PYTHIA}");
    hratio->GetYaxis()->SetNdivisions(5);
    hratio->GetYaxis()->SetTitleSize(0.04);
    hratio->GetYaxis()->SetTitleOffset(1.5);

    hratio->Draw();

    // // draw line at 1
    drawHoriLine(0., 200., 1., kGray+2, 1)->Draw();

    

    //PT STUFF FIGURE OUT WHERE TO PUT THIS LATER
    std::string fname_pt = outdir + "herwig-pythia_comparison_D0_pt_R0.4.pdf";
    const char* fname_ptc = fname_pt.c_str();
    c_pt->SaveAs(fname_ptc);
    delete c_pt;
    
    //--------------------------------------------------------//
    // now do z

    TCanvas* c_z = new TCanvas();
    ProcessCanvas(c_z);
    c_z->cd();
    gPad->SetLogy();
    gPad->SetBottomMargin(0.31);


    TPad *pad3 = new TPad("pad3","pad3",0.,0.,1.,1.);

    hD0z_pythia->GetXaxis()->SetTitle("D^{0} z");
    hD0z_herwig->GetXaxis()->SetTitle("D^{0} z");

    FormatHist(l2, hD0z_pythia, label1, markercolor1, markerstyle1, 0.80);
    FormatHist(l2, hD0z_herwig, label3, markercolor3, markerstyle3, 0.80);

    hD0z_pythia->Draw("same");
    hD0z_herwig->Draw("same");

    l_fake->Draw("same");
    hD0z_pythia->GetXaxis()->SetLabelSize(0);

    TPad *pad4 = new TPad("pad3","",0.,0.,1.,1.);
    pad4->SetTopMargin(0.71);
    pad4->SetFillColor(0);
    pad4->SetFillStyle(0);
    pad4->Draw();
    // pad4->SetLogx();
    pad4->SetLogy();
    pad4->cd();

    ratio_name = "hratio_D0z";
    TH1D* hratioz = (TH1D*) hD0z_herwig->Clone(ratio_name.c_str());
    hratioz->Divide(hD0z_pythia);
    // hratioz->SetMinimum(0.5);
    hratioz->SetMaximum(10.);


    FormatHist(l2, hratioz, "ratio", kBlack, markers[2]); //, 0.05, 0.04, 1.2, 0.035, 0.03, 1.5);
    hratioz->GetYaxis()->SetTitle("#frac{Herwig}{PYTHIA}");
    hratioz->GetYaxis()->SetNdivisions(5);
    hratioz->GetYaxis()->SetTitleSize(0.04);
    hratioz->GetYaxis()->SetTitleOffset(1.5);

    hratioz->Draw();

    // // draw line at 1
    drawHoriLine(0., 200., 1., kGray+2, 1)->Draw();

    

    //PT STUFF FIGURE OUT WHERE TO PUT THIS LATER
    std::string fname_z = outdir + "herwig-pythia_comparison_D0_z_R0.4.pdf";
    const char* fname_zc = fname_z.c_str();
    c_z->SaveAs(fname_zc);
    delete c_z;



    //--------------------------------------------------------//



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

        // Names of histograms in the file
        const std::string pythia_Djet_EEC_name = "hsparsejet_c_clone_proj_4_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max); //NOTE: after thnsparse axis change -> this projection # turned to a 4
        const std::string pythia_DjetwithDstar_EEC_name = "hsparsejet_c_Dstar_clone_proj_3_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
        // const std::string pythia_Djet_EEC_noweight_name = "hsparsejet_c_uw_clone_proj_3_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
        
        // const std::string herwig_Djet_EEC_name = "hsparse_EEC_clone_proj_3_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_trkthrd1.0";
        // const std::string herwig_DjetwithDstar_EEC_name = "hsparse_EEC_clone_proj_3_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_trkthrd1.0";
        // // const std::string herwig_EEC_ptrl_name = "hsparse_EEC_ptrl_clone_proj_3_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_trkthrd1.0";
        // const std::string herwig_EEC_noweight_name = "hsparse_EEC_noweight_clone_proj_3_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_trkthrd1.0";

        const std::string herwig_Djet_EEC_name = "D0_noDstar_EEC_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_R0.4_trkthrd1.0";
        const std::string herwig_DjetwithDstar_EEC_name = "D0_EEC_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_R0.4_trkthrd1.0";
        
        const std::string sherpa_c_jetpt10_EEC_name = "c_jetpt10_EEC_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_R0.4_trkthrd1.0";
        const std::string sherpa_c_jetpt10_lund_EEC_name = "c_jetpt10_lund_EEC_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_R0.4_trkthrd1.0";
        const std::string sherpa_c_jetpt15_EEC_name = "c_jetpt15_EEC_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_R0.4_trkthrd1.0";
        const std::string sherpa_c_jetpt15_lund_EEC_name = "c_jetpt15_lund_EEC_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_R0.4_trkthrd1.0";
        
        const std::string pythia_Djet_z_name = "hD0z_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
        const std::string herwig_Djet_z_perpt_name = "D0_noDstar_D0_z_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_trkthrd1.0";
        const std::string herwig_DjetwithDstar_z_perpt_name = "D0_D0_z_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max) + "_trkthrd1.0";
        
        
        // make a canvas for each pt range
        TCanvas* c = new TCanvas();
        ProcessCanvas(c);
        c->cd();
        gPad->SetLogx();
        if (logstring) {
            gPad->SetLogy();
        }

        TCanvas* c_z2 = new TCanvas();
        ProcessCanvas(c_z2);
        c_z2->cd();
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetBottomMargin(0.31);
        TPad *pad1a = new TPad("pad1a","pad1a",0.,0.,1.,1.);
        

        TLegend* l; // = new TLegend(0.17, 0.65, 0.5, 0.85);

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

        TH1D* hc_jetpt10_sherpa = (TH1D*) f_sherpa->Get(sherpa_c_jetpt10_EEC_name.c_str());
        TH1D* hc_jetpt10_lund_sherpa = (TH1D*) f_sherpa->Get(sherpa_c_jetpt10_lund_EEC_name.c_str());
        TH1D* hc_jetpt15_sherpa = (TH1D*) f_sherpa->Get(sherpa_c_jetpt15_EEC_name.c_str());
        TH1D* hc_jetpt15_lund_sherpa = (TH1D*) f_sherpa->Get(sherpa_c_jetpt15_lund_EEC_name.c_str());

        TH1D* hD0_z_pythia = (TH1D*) f_pythia_D0->Get(pythia_Djet_z_name.c_str());
        TH1D* hD0_z_herwig = (TH1D*) f_herwig->Get(herwig_Djet_z_perpt_name.c_str()); 
        TH1D* hD0wDstar_z_herwig = (TH1D*) f_herwig->Get(herwig_DjetwithDstar_z_perpt_name.c_str()); //not used

        // Format histograms for plotting (this order needed to keep legend in order and graphs lookin good)
        hD0_pythia->GetXaxis()->SetTitle("#it{R}_{L}");
        hD0_pythia->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
        hD0wDstar_pythia->GetXaxis()->SetTitle("#it{R}_{L}");
        hD0wDstar_pythia->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
        hD0_herwig->GetXaxis()->SetTitle("#it{R}_{L}");
        hD0_herwig->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
        hD0wDstar_herwig->GetXaxis()->SetTitle("#it{R}_{L}");
        hD0wDstar_herwig->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
        // not bothering formatting sherpa

        hD0_z_pythia->GetXaxis()->SetTitle("D^{0} z");
        
        cout << "about to format D0" << endl;
        
        
        FormatHist(l, hD0_pythia, label1, markercolor1, markerstyle1, 0.80); //FormatHist(l, hD0, "D^{0}-tagged, c-init jets", kMagenta+3, kOpenSquare);
        FormatHist(l, hD0_herwig, label3, markercolor3, markerstyle3, 0.80);
        if (plot_case == 1) {
            FormatHist(l, hD0wDstar_pythia, label2, markercolor2, markerstyle2, 0.80);
            FormatHist(l, hD0wDstar_herwig, label4, markercolor4, markerstyle4, 0.80);
        }
        if (plot_case == 2) {
            for (int j=0; j < hc_jetpt10_sherpa->GetNbinsX();j++){
                hc_jetpt10_sherpa->SetBinError(j+1, 0);
                hc_jetpt10_lund_sherpa->SetBinError(j+1, 0);
                hc_jetpt15_sherpa->SetBinError(j+1, 0);
                hc_jetpt15_lund_sherpa->SetBinError(j+1, 0);
            }

            FormatHistwithLine(l, hc_jetpt10_sherpa, label5, markercolor5, markerstyle5, 0.80);
            FormatHistwithLine(l, hc_jetpt10_lund_sherpa, label6, markercolor6, markerstyle6, 0.80);
            FormatHistwithLine(l, hc_jetpt15_sherpa, label7, markercolor7, markerstyle7, 0.80);
            FormatHistwithLine(l, hc_jetpt15_lund_sherpa, label8, markercolor8, markerstyle8, 0.80);
            
        }

         
        // get maximum
        hD0_pythia->SetMaximum(hD0_pythia->GetMaximum()*1.5);
        
        c->cd();
        hD0_pythia->Draw("L same");
        hD0_herwig->Draw("L same");
        if (plot_case == 1) {
            hD0wDstar_pythia->Draw("L same");
            hD0wDstar_herwig->Draw("L same");
        } else if (plot_case == 2) {
            hc_jetpt10_sherpa->Draw("L same");
            hc_jetpt10_lund_sherpa->Draw("L same");
            hc_jetpt15_sherpa->Draw("L same");
            hc_jetpt15_lund_sherpa->Draw("L same");
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
        
        } else if (plot_case == 2) {
            double hc_jetpt10_sherpa_top_binpos = findTopOfCurve(hc_jetpt10_sherpa, 1, 0.08);
            drawVertLine(hc_jetpt10_sherpa->GetBinCenter(hc_jetpt10_sherpa_top_binpos), 0, hc_jetpt10_sherpa->GetBinContent(hc_jetpt10_sherpa_top_binpos), markercolor5, 1)->Draw();
            double hc_jetpt10_lund_sherpa_top_binpos = findTopOfCurve(hc_jetpt10_lund_sherpa, 1, 0.08);
            drawVertLine(hc_jetpt10_lund_sherpa->GetBinCenter(hc_jetpt10_lund_sherpa_top_binpos), 0, hc_jetpt10_lund_sherpa->GetBinContent(hc_jetpt10_lund_sherpa_top_binpos), markercolor6, 1)->Draw();
            double hc_jetpt15_sherpa_top_binpos = findTopOfCurve(hc_jetpt15_sherpa, 1, 0.08);
            drawVertLine(hc_jetpt15_sherpa->GetBinCenter(hc_jetpt15_sherpa_top_binpos), 0, hc_jetpt15_sherpa->GetBinContent(hc_jetpt15_sherpa_top_binpos), markercolor7, 1)->Draw();
            double hc_jetpt15_lund_sherpa_top_binpos = findTopOfCurve(hc_jetpt15_lund_sherpa, 1, 0.08);
            drawVertLine(hc_jetpt15_lund_sherpa->GetBinCenter(hc_jetpt15_lund_sherpa_top_binpos), 0, hc_jetpt15_lund_sherpa->GetBinContent(hc_jetpt15_lund_sherpa_top_binpos), markercolor8, 1)->Draw();
        
        }

        // vector<double> fullwidth_vec = findWidthOfCurve(hD0,  hD0_top_binpos);
        // drawHoriLine(fullwidth_vec[1], fullwidth_vec[2], fullwidth_vec[0], kMagenta+3, 1)->Draw();
        



        std::string fname = outdir + "herwig-pythia_comparison_pt" + std::to_string(pt_min) + '-' + std::to_string(pt_max) + "_R0.4_plotcase" + std::to_string(plot_case) + ".pdf";
        const char* fnamec = fname.c_str();
        c->SaveAs(fnamec);
        delete c;

        // save D0 z by pt range
        c_z2->cd();
        hD0_z_pythia->Draw("same");
        hD0_z_herwig->Draw("same");

        FormatHist(l2, hD0_z_pythia, label1, markercolor1, markerstyle1, 0.80);
        FormatHist(l2, hD0_z_herwig, label3, markercolor3, markerstyle3, 0.80);

        l_fake->Draw("same");
        hD0_z_pythia->GetXaxis()->SetLabelSize(0);

        TPad *pad2a = new TPad("pad1a","",0.,0.,1.,1.);
        pad2a->SetTopMargin(0.71);
        pad2a->SetFillColor(0);
        pad2a->SetFillStyle(0);
        pad2a->Draw();
        // pad2a->SetLogx();
        pad2a->SetLogy();
        pad2a->cd();

        ratio_name = "hratio_D0_z";
        TH1D* hratio_z = (TH1D*) hD0_z_herwig->Clone(ratio_name.c_str());
        hratio_z->Divide(hD0_z_pythia);
        // hratio_z->SetMinimum(0.5);
        hratio_z->SetMaximum(10.);


        FormatHist(l2, hratio_z, "ratio", kBlack, markers[2]); //, 0.05, 0.04, 1.2, 0.035, 0.03, 1.5);
        hratio_z->GetYaxis()->SetTitle("#frac{Herwig}{PYTHIA}");
        hratio_z->GetYaxis()->SetNdivisions(5);
        hratio_z->GetYaxis()->SetTitleSize(0.04);
        hratio_z->GetYaxis()->SetTitleOffset(1.5);

        hratio_z->Draw();

        // // draw line at 1
        drawHoriLine(0., 1., 1., kGray+2, 1)->Draw();

    

        //PT STUFF FIGURE OUT WHERE TO PUT THIS LATER
        std::string fname_z2 = outdir + "herwig-pythia_comparison_z_pt" + std::to_string(pt_min) + '-' + std::to_string(pt_max) + "_R0.4_plotcase" + std::to_string(plot_case) + ".pdf";
        const char* fname_z2c = fname_z2.c_str();
        c_z2->SaveAs(fname_z2c);
        delete c_z2;






        f_out->cd();
        hD0_pythia->Write();
        hD0wDstar_pythia->Write();
        hD0_herwig->Write();
        hD0wDstar_herwig->Write();

        hc_jetpt10_sherpa->Write();
        hc_jetpt10_lund_sherpa->Write();
        hc_jetpt15_sherpa->Write();
        hc_jetpt15_lund_sherpa->Write();

        hD0_z_pythia->Write();
        hD0_z_herwig->Write();

            
    } // pT bins loop

    f_out->cd();
    hD0_pt_pythia->Write();
    // hD0wDstar_pt_pythia->Write();
    hD0_pt_herwig->Write();
    hD0wDstar_pt_herwig->Write();

    hD0z_pythia->Write();
    hD0z_herwig->Write();





    f_pythia_D0->Close();
    delete f_pythia_D0;
    f_pythia_D0wDstar->Close();
    delete f_pythia_D0wDstar;
    f_herwig->Close();
    delete f_herwig;

    return;
}
