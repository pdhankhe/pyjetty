
#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <unordered_map>
#include <vector>

#include <TAttFill.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TLegend.h>
#include <TTree.h>
#endif

// #include "/software/users/blianggi/template.C"
#include "TMath.h"
#include "iostream"




// =============================================================================================================
// --------------------------------------------- template functions --------------------------------------------
// =============================================================================================================


Double_t Errorsum(Double_t dx, Double_t dy) {
    return sqrt((dx)*(dx)+(dy)*(dy));
}


Double_t ErrorCalUnCorr(Double_t x, Double_t y, Double_t dx, Double_t dy) {
    return (x/y)*sqrt((dx/x)*(dx/x)+(dy/y)*(dy/y));
}

void LoadLibs() {
  	gSystem->Load("libCore.so");
  	gSystem->Load("libGeom.so");
  	gSystem->Load("libPhysics.so");
  	gSystem->Load("libVMC");
  	gSystem->Load("libTree");
  	gSystem->Load("libMinuit");
  	/*
  	gSystem->Load("libSTEERBase");
  	gSystem->Load("libESD");
  	gSystem->Load("libAOD");
  	gSystem->Load("libANALYSIS");
  	gSystem->Load("libANALYSISalice");
  	gSystem->Load("libCORRFW");
  	gSystem->Load("libPWGTools");
	*/
}

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
  	gStyle->SetLabelOffset(0.01,"y");
  	gStyle->SetLabelOffset(0.01,"x");
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

void ProcessHisto(TH1 *h, Double_t size=1.2, Int_t col=1, Int_t style=20) { 
	//gPad->SetTickx(); 
	//gPad->SetTicky(); 
	h->SetMarkerSize(size); 
	h->SetMarkerColor(col);
	h->SetLineWidth(2);
	h->SetLineColor(col);
	h->SetMarkerStyle(style); 
	h->GetYaxis()->SetTitleOffset(1.); 
	h->GetYaxis()->SetTitleSize(0.042);
	h->GetYaxis()->SetLabelSize(0.042);
	h->GetYaxis()->SetLabelFont(42);
	h->GetXaxis()->SetLabelFont(42);
	h->GetYaxis()->SetTitleFont(42);
	h->GetXaxis()->SetTitleFont(42);
	h->GetXaxis()->SetTitleOffset(1.0);
	h->GetXaxis()->SetTitleSize(0.042);
	h->GetXaxis()->SetLabelSize(0.042);
}

void ProcessLegend(TLegend *leg) {
	leg->SetFillColor(kWhite);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetTextSize(0.043);
	leg->Draw();
}

void ProcessLtxText(TLatex *ltx, Double_t size=0.04, Int_t col=1) {
	ltx->SetTextColor(col);
	ltx->SetTextSize(size);
	ltx->SetTextFont(42);
	ltx->SetNDC();
	ltx->Draw();
}

//---------------------------------------------------------------------------//

void makeCutsOnThNSparse(THnSparse* hsparse, int pt_min, int pt_max) {
    hsparse->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
    hsparse->GetAxis(1)->SetRangeUser(5., pt_max); // apply cut on Dmeson pt
    hsparse->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
        
}

TLine * drawHoriLine(double x1, double x2, double y1, int color, int linestyle=2){
    auto fhoriline = new TLine(x1, y1, x2, y1);
	fhoriline->SetLineWidth(1);
    fhoriline->SetLineColor(color);
    fhoriline->SetLineStyle(linestyle);
    return fhoriline;

}



// =============================================================================================================
// ------------------------------------------------- MAIN METHOD -----------------------------------------------
// =============================================================================================================

void PairEfficiency_PYTHIA() {
	
	// some styling setup
	LoadLibs();
	SetStyle();
	Double_t markers[7] = {kFullCircle, kFullSquare, kFullStar, kFullDiamond, kOpenCircle, kOpenSquare, kOpenDiamond};
    Double_t marker_size = 1.2;
    Double_t colors[10] = {kGreen+1, kRed+2, kGreen+2, kBlue+2, kOrange+2, kViolet+1, kYellow+1,kRed+1, kOrange+1, kCyan+1};	

    bool using_hiccup = true;
    // pt?

	// input and output files 
    TString path_to_inputfiles;
    TString infile_D0_truth; //const char infile_D0_truth[];
    TString infile_D0_det; //const char infile_D0_det[];
    TString infile_D0wDstar_truth; //const char infile_D0wDstar_truth[];
    TString infile_D0wDstar_det; //const char infile_D0wDstar_det[];

    if (using_hiccup) { //hiccup
        path_to_inputfiles = "/rstorage/alice/AnalysisResults/blianggi/EEC/";
        infile_D0_truth = "from_perlmutter/18063568/AnalysisResultsFinal.root";
        infile_D0_det = "60901/AnalysisResultsFinal.root";
        infile_D0wDstar_truth = "from_perlmutter/18044499/AnalysisResultsFinal.root";
        infile_D0wDstar_det = "60902/AnalysisResultsFinal.root";
    } else { //perlmutter
        path_to_inputfiles = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/";
        infile_D0_truth = "17651853/AnalysisResultsFinal.root";
        infile_D0_det = "17651853/AnalysisResultsFinal.root";
        infile_D0wDstar_truth = "17651853/AnalysisResultsFinal.root";
        infile_D0wDstar_det = "17651853/AnalysisResultsFinal.root";
    }
    

	
	// TString inputfile = "1252842/AnalysisResultsFinal.root";
    // const char* outputfile = fname.c_str();
	// TString outputfile = "plots/pair_efficiency.pdf";
	
	// TString test_inputfile = "/software/users/blianggi/pyjetty/pyjetty/alihfjets/dev/hfjet/process/user/hf_EEC/TestOutput/AnalysisResults.root";
	// TString test_outputfile = "./test_pair_efficiency.pdf";

    
	
	// Open files
	TFile *file_D0_truth = TFile::Open(path_to_inputfiles + infile_D0_truth, "read");
    TFile *file_D0_det = TFile::Open(path_to_inputfiles + infile_D0_det, "read");
    TFile *file_D0wDstar_truth = TFile::Open(path_to_inputfiles + infile_D0wDstar_truth, "read");
    TFile *file_D0wDstar_det = TFile::Open(path_to_inputfiles + infile_D0wDstar_det, "read");
    
    // Get histograms
    std::string jetR = "0.4";
    const std::string hc_name = "h_EEC_JetPt_charm_R" + jetR;
    const std::string hc_jet_name = "h_JetPt_charm_R" + jetR + "_jetlevel";
    const int pt_bins[] = { 10, 15, 30 }; //{ 10, 20, 40 };
    const int n_bins = 2;
    

    // Get histograms
    THnSparse* hsparsejet_c_D0_truth = (THnSparse*) file_D0_truth->Get(hc_name.c_str());
    THnSparse* hsparsejet_c_jetlevel_D0_truth = (THnSparse*) file_D0_truth->Get(hc_jet_name.c_str());  
    THnSparse* hsparsejet_c_D0_det = (THnSparse*) file_D0_det->Get(hc_name.c_str());
    THnSparse* hsparsejet_c_jetlevel_D0_det = (THnSparse*) file_D0_det->Get(hc_jet_name.c_str());
    
    THnSparse* hsparsejet_c_D0wDstar_truth = (THnSparse*) file_D0wDstar_truth->Get(hc_name.c_str());
    THnSparse* hsparsejet_c_jetlevel_D0wDstar_truth = (THnSparse*) file_D0wDstar_truth->Get(hc_jet_name.c_str());
    THnSparse* hsparsejet_c_D0wDstar_det = (THnSparse*) file_D0wDstar_det->Get(hc_name.c_str());
    THnSparse* hsparsejet_c_jetlevel_D0wDstar_det = (THnSparse*) file_D0wDstar_det->Get(hc_jet_name.c_str());
    
    // for THnSparse: make clone to work with, make cuts, get projection
    THnSparse *hsparsejet_c_D0_truth_clone = (THnSparse *) hsparsejet_c_D0_truth->Clone("hsparsejet_c_D0_truth_clone");
    THnSparse *hsparsejet_c_jetlevel_D0_truth_clone = (THnSparse *) hsparsejet_c_jetlevel_D0_truth->Clone("hsparsejet_c_jetlevel_D0_truth_clone");
    THnSparse *hsparsejet_c_D0_det_clone = (THnSparse *) hsparsejet_c_D0_det->Clone("hsparsejet_c_D0_det_clone");
    THnSparse *hsparsejet_c_jetlevel_D0_det_clone = (THnSparse *) hsparsejet_c_jetlevel_D0_det->Clone("hsparsejet_c_jetlevel_D0_det_clone");

    THnSparse *hsparsejet_c_D0wDstar_truth_clone = (THnSparse *) hsparsejet_c_D0wDstar_truth->Clone("hsparsejet_c_D0wDstar_truth_clone");
    THnSparse *hsparsejet_c_jetlevel_D0wDstar_truth_clone = (THnSparse *) hsparsejet_c_jetlevel_D0wDstar_truth->Clone("hsparsejet_c_jetlevel_D0wDstar_truth_clone");
    THnSparse *hsparsejet_c_D0wDstar_det_clone = (THnSparse *) hsparsejet_c_D0wDstar_det->Clone("hsparsejet_c_D0wDstar_det_clone");
    THnSparse *hsparsejet_c_jetlevel_D0wDstar_det_clone = (THnSparse *) hsparsejet_c_jetlevel_D0wDstar_det->Clone("hsparsejet_c_jetlevel_D0wDstar_det_clone");


    //ieworijwoij
    for (int i = 0; i < n_bins; i++) {
        cout << "in pt bin" << i << endl;
        int pt_min = pt_bins[i];
        int pt_max = pt_bins[i+1];
        TString ptbin = TString::Format("%d #leq #it{p}_{T}^{ch. jet} < %d GeV/#it{c}, #font[122]{|}#it{#eta}_{jet}#font[122]{|} #leq 0.5", pt_min, pt_max);
        TString ptD = TString::Format("5 #leq #it{p}_{T}^{D^{0}} < %d GeV/#it{c}, #font[122]{|}#it{y}_{D^{0}}#font[122]{|} #leq 0.8", pt_max);
        // TString ptbin = pt_min.c_str() + " #leq #it{p}_{T}^{ch. jet} < " + pt_max.c_str() + " GeV/#it{c}, #font[122]{|}#it{#eta}_{jet}#font[122]{|} #leq 0.5";
        // TString ptD = "5 #leq #it{p}_{T}^{D^{0}} < " + pt_max.c_str() + " GeV/#it{c}, #font[122]{|}#it{y}_{D^{0}}#font[122]{|} #leq 0.8";
        
        TString outdir = "plots/";
        TString outputfile = outdir + TString::Format("pair_efficiency_pt%d-%d_R%s.pdf", pt_min, pt_max, jetR.c_str()); 
    

        // make cuts on jet pT, D^0 pT, D^0 rapidity
        makeCutsOnThNSparse(hsparsejet_c_D0_truth_clone, pt_min, pt_max);
        makeCutsOnThNSparse(hsparsejet_c_jetlevel_D0_truth_clone, pt_min, pt_max);
        makeCutsOnThNSparse(hsparsejet_c_D0_det_clone, pt_min, pt_max);
        makeCutsOnThNSparse(hsparsejet_c_jetlevel_D0_det_clone, pt_min, pt_max);

        makeCutsOnThNSparse(hsparsejet_c_D0wDstar_truth_clone, pt_min, pt_max);
        makeCutsOnThNSparse(hsparsejet_c_jetlevel_D0wDstar_truth_clone, pt_min, pt_max);
        makeCutsOnThNSparse(hsparsejet_c_D0wDstar_det_clone, pt_min, pt_max);
        makeCutsOnThNSparse(hsparsejet_c_jetlevel_D0wDstar_det_clone, pt_min, pt_max);

        
        
        // hsparsejet_c_D0_truth_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
        // hsparsejet_c_D0_truth_clone->GetAxis(1)->SetRangeUser(5., pt_max); // apply cut on Dmeson pt
        // hsparsejet_c_D0_truth_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
        // hsparsejet_c_jetlevel_D0_truth_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
        // hsparsejet_c_jetlevel_D0_truth_clone->GetAxis(1)->SetRangeUser(5., pt_max); // apply cut on Dmeson pt
        // hsparsejet_c_jetlevel_D0_truth_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
        
        // hsparsejet_c_D0_truth_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
        // hsparsejet_c_D0_truth_clone->GetAxis(1)->SetRangeUser(5., pt_max); // apply cut on Dmeson pt
        // hsparsejet_c_D0_truth_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
        // hsparsejet_c_jetlevel_D0_truth_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
        // hsparsejet_c_jetlevel_D0_truth_clone->GetAxis(1)->SetRangeUser(5., pt_max); // apply cut on Dmeson pt
        // hsparsejet_c_jetlevel_D0_truth_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
        
        // hsparsejet_c_D0_truth_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
        // hsparsejet_c_D0_truth_clone->GetAxis(1)->SetRangeUser(5., pt_max); // apply cut on Dmeson pt
        // hsparsejet_c_D0_truth_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
        // hsparsejet_c_jetlevel_D0_truth_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
        // hsparsejet_c_jetlevel_D0_truth_clone->GetAxis(1)->SetRangeUser(5., pt_max); // apply cut on Dmeson pt
        // hsparsejet_c_jetlevel_D0_truth_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
        
        // hsparsejet_c_D0_truth_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
        // hsparsejet_c_D0_truth_clone->GetAxis(1)->SetRangeUser(5., pt_max); // apply cut on Dmeson pt
        // hsparsejet_c_D0_truth_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
        // hsparsejet_c_jetlevel_D0_truth_clone->GetAxis(0)->SetRangeUser(pt_min, pt_max); // apply cut on jet pt
        // hsparsejet_c_jetlevel_D0_truth_clone->GetAxis(1)->SetRangeUser(5., pt_max); // apply cut on Dmeson pt
        // hsparsejet_c_jetlevel_D0_truth_clone->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
        
        // Project onto observable axis
        TH1D *hD0_truth = hsparsejet_c_D0_truth_clone->Projection(3); 
        TH1D *hD0_truth_jet_pt = hsparsejet_c_jetlevel_D0_truth_clone->Projection(0);
        TH1D *hD0_det = hsparsejet_c_D0_det_clone->Projection(3); 
        TH1D *hD0_det_jet_pt = hsparsejet_c_jetlevel_D0_det_clone->Projection(0);

        TH1D *hD0wDstar_truth = hsparsejet_c_D0wDstar_truth_clone->Projection(3); 
        TH1D *hD0wDstar_truth_jet_pt = hsparsejet_c_jetlevel_D0wDstar_truth_clone->Projection(0);
        TH1D *hD0wDstar_det = hsparsejet_c_D0wDstar_det_clone->Projection(3); 
        TH1D *hD0wDstar_det_jet_pt = hsparsejet_c_jetlevel_D0wDstar_det_clone->Projection(0);

        // Set to appropriate name
        // std::string hname = hD0->GetName();
        // hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
        // hD0->SetNameTitle(hname.c_str(), hname.c_str());

        // Find normalization factor
        double numjets_D0_truth = hD0_truth_jet_pt->Integral();
        double numjets_D0_det = hD0_det_jet_pt->Integral();
        double numjets_D0wDstar_truth = hD0wDstar_truth_jet_pt->Integral();
        double numjets_D0wDstar_det = hD0wDstar_det_jet_pt->Integral();

        // Normalize
        hD0_truth->Scale(1/numjets_D0_truth, "width");
        hD0_det->Scale(1/numjets_D0_det, "width");
        hD0wDstar_truth->Scale(1/numjets_D0wDstar_truth, "width");
        hD0wDstar_det->Scale(1/numjets_D0wDstar_det, "width");

        // Rebin?
        // int n_obs_bins = 5; //50; //-1;
        // double obs_bins[] = {1,2,3,4,5};


        // Take the ratio
        TH1D *hratio_pairs_truth = (TH1D *) hD0_truth->Clone("hratio_pairs_truth");
        TH1D *hratio_pairs_det = (TH1D *) hD0_det->Clone("hratio_pairs_det");
        hratio_pairs_truth->Divide(hD0wDstar_truth);
        hratio_pairs_det->Divide(hD0wDstar_det);


        // Format plot
        // int markercolor1 = kGreen-5; //D0, me
        // int markerstyle1 = kFullCircle;
        // int markercolor2 = 0;
        // int markerstyle2 = 0;

        hratio_pairs_truth->SetTitle("Truth Level Pair Efficiency");
        hratio_pairs_det->SetTitle("Detector Level Pair Efficiency");
        hratio_pairs_truth->GetXaxis()->SetTitle("#it{R}_{L}"); //"rec R_{L} / gen R_{L}");
        hratio_pairs_det->GetXaxis()->SetTitle("#it{R}_{L}"); //"rec R_{L} / gen R_{L}");
        hratio_pairs_truth->GetYaxis()->SetTitle(" c->D^{0} / (c->D^{0} + c->D^{*}) "); //"rec R_{L} / gen R_{L}");
        hratio_pairs_det->GetYaxis()->SetTitle(" c->D^{0} / (c->D^{0} + c->D^{*}) "); //"rec R_{L} / gen R_{L}");
        // EEC y axis label: "#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}"
        
        ProcessHisto(hratio_pairs_truth, marker_size, kRed+2, markers[0]);
        ProcessHisto(hratio_pairs_det, marker_size, kBlue+2, markers[4]);

        // set y axis limits
        hratio_pairs_truth->SetMaximum(2.);
        hratio_pairs_det->SetMaximum(2.);
        

        // make and format the canvas
        TCanvas *c1 = new TCanvas("c1", "c1");
        ProcessCanvas(c1);
        gPad->SetLogx();

        // plot the ratio
        hratio_pairs_truth->Draw();
        hratio_pairs_det->Draw("same");

        drawHoriLine(1e-4, 1, 1, kBlack, 3)->Draw();


        // Draw legends

        TLegend* l = new TLegend(0.1797168,0.550741,0.5062155,0.8885185,""); //(0.17, 0.4, 0.5, 0.53);
        l->SetTextSize(0.045);
        l->SetFillColor(kWhite);
        // TLegend *leg = new TLegend(0.1797168,0.5390741,0.4562155,0.8885185,"");
        l->AddEntry("NULL","PYTHIA 8 Monash 2013","h");
        l->AddEntry("NULL","pp, #sqrt{#it{s}} = 13 TeV","h");
        l->AddEntry("NULL","D^{0} #rightarrow K^{#minus} #pi^{+} and charge conj.","h");
        l->AddEntry("NULL","in charged jets, anti-#it{k}_{T}, #it{R} = 0.4","h");
        l->AddEntry("NULL",ptbin,"h");
        l->AddEntry("NULL",ptD,"h");
        l->SetTextSize(0.037);
        l->SetBorderSize(0);
        l->Draw("same");

        auto legend_labels = new TLegend(0.7,0.65,0.85,0.85);
        legend_labels->AddEntry(hratio_pairs_truth,"truth level","lep");
        legend_labels->AddEntry(hratio_pairs_det,"detector level","lep"); //"f" = boxes
        legend_labels->Draw("same");

        c1->SaveAs(outputfile);
        delete c1;


    }

  
    
    
	// // make a legend
	// TLegend *leg= new TLegend( 0.49, 0.7, 0.73, 0.95, NULL,"brNDC");
    // 	leg->AddEntry("NULL","pp, #sqrt{#it{s}} = 13.02 TeV","h");
    // 	leg->AddEntry("NULL","D^{0}-tagged ch. jets, anti-k_{T}, #it{R} = 0.4","h");
    // 	leg->AddEntry("NULL","10 < #it{p}_{T}^{ch jet} < 20 GeV/#it{c}, #font[122]{|}#it{#eta}_{jet}#font[122]{|} #leq 0.5","h");
    // 	leg->AddEntry("NULL","5 < #it{p}_{T}^{D^{0}} < 20 GeV/#it{c}, #font[122]{|}#it{y}_{D^{0}}#font[122]{|} #leq 0.8","h");
    // 	ProcessLegend(leg);
	// leg->Draw("same");
    
    
    
    
    
    
    // // Plot other stuff - be more specific later
    // auto legend_genvsreco = new TLegend(0.2,0.5,0.35,0.7);
    // legend_genvsreco->AddEntry(hist_jet_gen_prompt_norm_RL,"generated","lep");
    // legend_genvsreco->AddEntry(hist_jet_reco_prompt_norm_RL,"reconstructed","lep"); //"f" = boxes
    
    
    // TCanvas *c2 = new TCanvas("c2","c2");
    // ProcessCanvas(c2);
    // gPad->SetLogx();
    
    // ProcessHisto(hist_jet_gen_prompt_norm_RL, marker_size, kRed+2, markers[0]);
    // ProcessHisto(hist_jet_reco_prompt_norm_RL, marker_size, kBlue+2, markers[0]);
    // hist_jet_reco_prompt_norm_RL->SetTitle("Prompt Pairs");
    // hist_jet_gen_prompt_norm_RL->Draw();
    // hist_jet_reco_prompt_norm_RL->Draw("same");
    
    // TLegend *leg_prompt = new TLegend( 0.2,0.7,0.35,0.8, NULL,"brNDC");
    // leg_prompt->AddEntry("NULL","Prompt jets","h");
    // ProcessLegend(leg_prompt);
    // leg_prompt->Draw("same");
    // legend_genvsreco->Draw("same");
    
    
    // TCanvas *c3 = new TCanvas("c3","c3");
    // ProcessCanvas(c3);
    // gPad->SetLogx();
    
    // ProcessHisto(hist_jet_gen_nonprompt_norm_RL, marker_size, kRed+2, markers[0]);
    // ProcessHisto(hist_jet_reco_nonprompt_norm_RL, marker_size, kBlue+2, markers[0]);
    // hist_jet_reco_nonprompt_norm_RL->SetTitle("Non-Prompt Pairs");
    // hist_jet_gen_nonprompt_norm_RL->Draw();
    // hist_jet_reco_nonprompt_norm_RL->Draw("same");
    
    // TLegend *leg_nonprompt = new TLegend( 0.2,0.7,0.35,0.8, NULL,"brNDC");
    // leg_nonprompt->AddEntry("NULL","Non-prompt jets","h");
    // ProcessLegend(leg_nonprompt);
    // leg_nonprompt->Draw("same");
    // legend_genvsreco->Draw("same");
    
    
    
    
    
    
    file_D0_truth->Close();
    delete file_D0_truth;
    file_D0_det->Close();
    delete file_D0_det;

    file_D0wDstar_truth->Close();
    delete file_D0wDstar_truth;
    file_D0wDstar_det->Close();
    delete file_D0wDstar_det;


}

