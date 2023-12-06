
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

	// input and output files 
    TString path_to_inputfiles = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/";
    const char infile_D0_truth[] = "17651853/AnalysisResultsFinal.root";
    const char infile_D0_det[] = "17651853/AnalysisResultsFinal.root";
    const char infile_D0wDstar_truth[] = "17651853/AnalysisResultsFinal.root";
    const char infile_D0wDstar_det[] = "17651853/AnalysisResultsFinal.root";

	
	// TString inputfile = "1252842/AnalysisResultsFinal.root";
	TString outputfile = "plots/pair_efficiency.pdf";
	
	// TString test_inputfile = "/software/users/blianggi/pyjetty/pyjetty/alihfjets/dev/hfjet/process/user/hf_EEC/TestOutput/AnalysisResults.root";
	// TString test_outputfile = "./test_pair_efficiency.pdf";

    
	
	// Open file
	TFile *file0 = TFile::Open(path_to_inputfiles + inputfile, "read");
//    TFile *file0 = TFile::Open(test_inputfile, "read");
    THnSparse *hsparsejet_gen = (THnSparse *) file0->Get("hsparsejet_gen");
	THnSparse *hsparsejet_reco = (THnSparse *) file0->Get("hsparsejet_reco");
    THnSparse *hsparsejetcandinfo_gen = (THnSparse *) file0->Get("hsparsejetcand_gen");
    THnSparse *hsparsejetcandinfo_reco = (THnSparse *) file0->Get("hsparsejetcand_reco");


	// make clones to work with
    THnSparse *hsparsejet_gen_prompt_clone = (THnSparse *) hsparsejet_gen->Clone("hsparsejet_gen_prompt_clone");
    THnSparse *hsparsejet_gen_nonprompt_clone = (THnSparse *) hsparsejet_gen->Clone("hsparsejet_gen_nonprompt_clone");
    THnSparse *hsparsejet_reco_prompt_clone = (THnSparse *) hsparsejet_reco->Clone("hsparsejet_reco_prompt_clone");
    THnSparse *hsparsejet_reco_nonprompt_clone = (THnSparse *) hsparsejet_reco->Clone("hsparsejet_reco_nonprompt_clone");
    
    THnSparse *hsparsejetcandinfo_gen_prompt_clone = (THnSparse *) hsparsejetcandinfo_gen->Clone("hsparsejetcandinfo_gen_prompt_clone");
    THnSparse *hsparsejetcandinfo_gen_nonprompt_clone = (THnSparse *) hsparsejetcandinfo_gen->Clone("hsparsejetcandinfo_gen_nonprompt_clone");
    THnSparse *hsparsejetcandinfo_reco_prompt_clone = (THnSparse *) hsparsejetcandinfo_reco->Clone("hsparsejetcandinfo_reco_prompt_clone");
    THnSparse *hsparsejetcandinfo_reco_nonprompt_clone = (THnSparse *) hsparsejetcandinfo_reco->Clone("hsparsejetcandinfo_reco_nonprompt_clone");
    
    
    //make cuts on jet pT, D^0 pT, cand
    hsparsejet_gen_prompt_clone->GetAxis(0)->SetRangeUser(10.,20); // apply cut on jet pt
    hsparsejet_gen_prompt_clone->GetAxis(1)->SetRangeUser(5.,20.); // apply cut on Dmeson pt
    hsparsejet_gen_prompt_clone->GetAxis(3)->SetRangeUser(1.,2.); // apply cut on cand - prompt is 1-2
    hsparsejetcandinfo_gen_prompt_clone->GetAxis(0)->SetRangeUser(10.,20);
    hsparsejetcandinfo_gen_prompt_clone->GetAxis(1)->SetRangeUser(5.,20.);
    hsparsejetcandinfo_gen_prompt_clone->GetAxis(3)->SetRangeUser(1.,2.);
    
    hsparsejet_gen_nonprompt_clone->GetAxis(0)->SetRangeUser(10.,20);
    hsparsejet_gen_nonprompt_clone->GetAxis(1)->SetRangeUser(5.,20.);
    hsparsejet_gen_nonprompt_clone->GetAxis(3)->SetRangeUser(2.,3.); //nonprompt is 2-3
    hsparsejetcandinfo_gen_nonprompt_clone->GetAxis(0)->SetRangeUser(10.,20);
    hsparsejetcandinfo_gen_nonprompt_clone->GetAxis(1)->SetRangeUser(5.,20.);
    hsparsejetcandinfo_gen_nonprompt_clone->GetAxis(3)->SetRangeUser(2.,3.);
    
    hsparsejet_reco_prompt_clone->GetAxis(0)->SetRangeUser(10.,20);
    hsparsejet_reco_prompt_clone->GetAxis(1)->SetRangeUser(5.,20.);
    hsparsejet_reco_prompt_clone->GetAxis(3)->SetRangeUser(1.,2.);
    hsparsejetcandinfo_reco_prompt_clone->GetAxis(0)->SetRangeUser(10.,20);
    hsparsejetcandinfo_reco_prompt_clone->GetAxis(1)->SetRangeUser(5.,20.);
    hsparsejetcandinfo_reco_prompt_clone->GetAxis(3)->SetRangeUser(1.,2.);
    
    hsparsejet_reco_nonprompt_clone->GetAxis(0)->SetRangeUser(10.,20);
    hsparsejet_reco_nonprompt_clone->GetAxis(1)->SetRangeUser(5.,20.);
    hsparsejet_reco_nonprompt_clone->GetAxis(3)->SetRangeUser(2.,3.);
    hsparsejetcandinfo_reco_nonprompt_clone->GetAxis(0)->SetRangeUser(10.,20);
    hsparsejetcandinfo_reco_nonprompt_clone->GetAxis(1)->SetRangeUser(5.,20.);
    hsparsejetcandinfo_reco_nonprompt_clone->GetAxis(3)->SetRangeUser(2.,3.);


	// access the histograms for gen & reco, prompt & nonprompt
	TH1D *hist_jet_gen_prompt_RL = hsparsejet_gen_prompt_clone->Projection(4); //technically don't need to do this separately
    TH1D *hist_jet_gen_nonprompt_RL = hsparsejet_gen_nonprompt_clone->Projection(4);
	TH1D *hist_jet_reco_prompt_RL = hsparsejet_reco_prompt_clone->Projection(4);
    TH1D *hist_jet_reco_nonprompt_RL = hsparsejet_reco_nonprompt_clone->Projection(4);
    
    TH1D *hist_jet_gen_prompt_cand = hsparsejetcandinfo_gen_prompt_clone->Projection(3); // or take from the jet cand one thats per jet??
    TH1D *hist_jet_gen_nonprompt_cand = hsparsejetcandinfo_gen_nonprompt_clone->Projection(3);
    TH1D *hist_jet_reco_prompt_cand = hsparsejetcandinfo_reco_prompt_clone->Projection(3);
    TH1D *hist_jet_reco_nonprompt_cand = hsparsejetcandinfo_reco_nonprompt_clone->Projection(3);
    
    
    //testing
    cout << "num bins " << hist_jet_gen_prompt_cand->GetNbinsX() << endl;
    cout << "num bins " << hist_jet_gen_nonprompt_cand->GetNbinsX() << endl;
    cout << "num bins " << hist_jet_reco_prompt_cand->GetNbinsX() << endl;
    cout << "num bins " << hist_jet_reco_nonprompt_cand->GetNbinsX() << endl;
    
    //get the number of prompt and nonprompt jets - need to increase the number +1 in GetBinContent bc of the way it is binned - should fix later.
    double numPromptJets_gen = hist_jet_gen_prompt_cand->GetBinContent(1);
    double numNonPromptJets_gen = hist_jet_gen_nonprompt_cand->GetBinContent(1);
    double numPromptJets_reco = hist_jet_reco_prompt_cand->GetBinContent(1);
    double numNonPromptJets_reco = hist_jet_reco_nonprompt_cand->GetBinContent(1);
    /*cout << "num reflection gen jets: " << hsparsejetcandinfo_gen->Projection(3)->GetBinContent(1) << endl;
    cout << "num reflection reco jets: " << hsparsejetcandinfo_reco->Projection(3)->GetBinContent(1) << endl;
    cout << "num prompt gen jets: " << hsparsejetcandinfo_gen->Projection(3)->GetBinContent(2) << endl;
    cout << "num nonprompt reco jets: " << hsparsejetcandinfo_reco->Projection(3)->GetBinContent(3) << endl;
    cout << "num prompt gen jets: " << hsparsejetcandinfo_gen->Projection(3)->GetBinContent(2) << endl;
    cout << "num nonprompt reco jets: " << hsparsejetcandinfo_reco->Projection(3)->GetBinContent(3) << endl;
    */
    cout << "num prompt gen jets: " << numPromptJets_gen << endl;
    cout << "num NON-prompt gen jets: " << numNonPromptJets_gen << endl;
    cout << "num prompt reco jets: " << numPromptJets_reco << endl;
    cout << "num NON-prompt reco jets: " << numNonPromptJets_reco << endl;
    
    // normalize
    TH1D *hist_jet_gen_prompt_norm_RL = (TH1D *) hist_jet_gen_prompt_RL->Clone("hist_jet_gen_prompt_norm_RL");
    TH1D *hist_jet_gen_nonprompt_norm_RL = (TH1D *) hist_jet_gen_nonprompt_RL->Clone("hist_jet_gen_nonprompt_norm_RL");
    TH1D *hist_jet_reco_prompt_norm_RL = (TH1D *) hist_jet_reco_prompt_RL->Clone("hist_jet_reco_prompt_norm_RL");
    TH1D *hist_jet_reco_nonprompt_norm_RL = (TH1D *) hist_jet_reco_nonprompt_RL->Clone("hist_jet_reco_nonprompt_norm_RL");
 
    
    hist_jet_gen_prompt_norm_RL->Scale(1 / numPromptJets_gen); // IS THIS RIGHT IDK
    hist_jet_gen_nonprompt_norm_RL->Scale(1 / numNonPromptJets_gen);
    hist_jet_reco_prompt_norm_RL->Scale(1 / numPromptJets_reco);
    hist_jet_reco_nonprompt_norm_RL->Scale(1 / numNonPromptJets_reco);
    cout << "CHECKPOINT 1" << endl;
    /*
    TCanvas *c1a = new TCanvas("c1a");
    hist_jet_gen_prompt_norm_RL->Draw();
    TCanvas *c2a = new TCanvas("c2a");
    hist_jet_gen_nonprompt_norm_RL->Draw();
    TCanvas *c3a = new TCanvas("c3a");
    hist_jet_reco_prompt_norm_RL->Draw();
    TCanvas *c4a = new TCanvas("c4a");
    hist_jet_reco_nonprompt_norm_RL->Draw();
    */
    

	// divide the rec/gen
//	TH1D *hratio_pairs = (TH1D *) hsparsejet_reco_RL->Clone("hratio_pairs"); //as a record of what works, delete later
//	hratio_pairs->Divide(hsparsejet_gen_RL);
    TH1D *hratio_pairs_prompt = (TH1D *) hist_jet_reco_prompt_RL->Clone("hratio_pairs_prompt");
    TH1D *hratio_pairs_nonprompt = (TH1D *) hist_jet_reco_nonprompt_RL->Clone("hratio_pairs_nonprompt");
    hratio_pairs_prompt->Divide(hist_jet_gen_prompt_RL);
    hratio_pairs_nonprompt->Divide(hist_jet_gen_nonprompt_RL);
    cout << "CHECKPOINT 2" << endl;

	// format the plot
//	hratio_pairs->SetTitle("Pair Efficiency"); //hratio_pairs->CenterTitle(); //again, what worked
//	hratio_pairs->GetYaxis()->SetTitle("rec R_{L} / gen R_{L}");
//	ProcessHisto(hratio_pairs, marker_size, kRed+2, markers[0]);
    hratio_pairs_prompt->SetTitle("Prompt Pair Efficiency");
    hratio_pairs_nonprompt->SetTitle("Non-prompt Pair Efficiency");
    hratio_pairs_prompt->GetYaxis()->SetTitle("rec R_{L} / gen R_{L}");
    hratio_pairs_nonprompt->GetYaxis()->SetTitle("rec R_{L} / gen R_{L}");
    ProcessHisto(hratio_pairs_prompt, marker_size, kRed+2, markers[0]);
    ProcessHisto(hratio_pairs_nonprompt, marker_size, kBlue+2, markers[4]);

    
    // make and format the canvas
    TCanvas *c1 = new TCanvas("c1", "c1");
    ProcessCanvas(c1);
    gPad->SetLogx();
    
    // plot the ratio
    hratio_pairs_prompt->Draw();
    hratio_pairs_nonprompt->Draw("same");
    
	// make a legend
	TLegend *leg= new TLegend( 0.49, 0.7, 0.73, 0.95, NULL,"brNDC");
    	leg->AddEntry("NULL","pp, #sqrt{#it{s}} = 13.02 TeV","h");
    	leg->AddEntry("NULL","D^{0}-tagged ch. jets, anti-k_{T}, #it{R} = 0.4","h");
    	leg->AddEntry("NULL","10 < #it{p}_{T}^{ch jet} < 20 GeV/#it{c}, #font[122]{|}#it{#eta}_{jet}#font[122]{|} #leq 0.5","h");
    	leg->AddEntry("NULL","5 < #it{p}_{T}^{D^{0}} < 20 GeV/#it{c}, #font[122]{|}#it{y}_{D^{0}}#font[122]{|} #leq 0.8","h");
    	ProcessLegend(leg);
	leg->Draw("same");
    
    auto legend_pnp = new TLegend(0.7,0.15,0.85,0.35);
    legend_pnp->AddEntry(hratio_pairs_prompt,"prompt","lep");
    legend_pnp->AddEntry(hratio_pairs_nonprompt,"nonprompt","lep"); //"f" = boxes
    legend_pnp->Draw("same");
    
    
    
    
    // Plot other stuff - be more specific later
    auto legend_genvsreco = new TLegend(0.2,0.5,0.35,0.7);
    legend_genvsreco->AddEntry(hist_jet_gen_prompt_norm_RL,"generated","lep");
    legend_genvsreco->AddEntry(hist_jet_reco_prompt_norm_RL,"reconstructed","lep"); //"f" = boxes
    
    
    TCanvas *c2 = new TCanvas("c2","c2");
    ProcessCanvas(c2);
    gPad->SetLogx();
    
    ProcessHisto(hist_jet_gen_prompt_norm_RL, marker_size, kRed+2, markers[0]);
    ProcessHisto(hist_jet_reco_prompt_norm_RL, marker_size, kBlue+2, markers[0]);
    hist_jet_reco_prompt_norm_RL->SetTitle("Prompt Pairs");
    hist_jet_gen_prompt_norm_RL->Draw();
    hist_jet_reco_prompt_norm_RL->Draw("same");
    
    TLegend *leg_prompt = new TLegend( 0.2,0.7,0.35,0.8, NULL,"brNDC");
    leg_prompt->AddEntry("NULL","Prompt jets","h");
    ProcessLegend(leg_prompt);
    leg_prompt->Draw("same");
    legend_genvsreco->Draw("same");
    
    
    TCanvas *c3 = new TCanvas("c3","c3");
    ProcessCanvas(c3);
    gPad->SetLogx();
    
    ProcessHisto(hist_jet_gen_nonprompt_norm_RL, marker_size, kRed+2, markers[0]);
    ProcessHisto(hist_jet_reco_nonprompt_norm_RL, marker_size, kBlue+2, markers[0]);
    hist_jet_reco_nonprompt_norm_RL->SetTitle("Non-Prompt Pairs");
    hist_jet_gen_nonprompt_norm_RL->Draw();
    hist_jet_reco_nonprompt_norm_RL->Draw("same");
    
    TLegend *leg_nonprompt = new TLegend( 0.2,0.7,0.35,0.8, NULL,"brNDC");
    leg_nonprompt->AddEntry("NULL","Non-prompt jets","h");
    ProcessLegend(leg_nonprompt);
    leg_nonprompt->Draw("same");
    legend_genvsreco->Draw("same");
    
    
    
    
    
    // save!
    c1->Print(outputfile + "(");
    c2->Print(outputfile);
    c3->Print(outputfile + ")");


}

