// ROOT macro to make D0 when using softQCD
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
   gStyle->SetPadLeftMargin(0.175);
  // gStyle->SetPadRightMargin(0.1);
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


void FormatHist(TLegend *l, TH1 *hist, TString text, int markercolor=1, int markerstyle=8, double markeralpha=1.,
              double xtitlesize=0.06, double xlabelsize=0.05, double xoffset=1.0,
              double ytitlesize=0.06, double ylabelsize=0.05, double yoffset=1.05, double markersize=1.5) 
              // double xtitlesize=0.04, double xlabelsize=0.04, double xoffset=1.2,
              // double ytitlesize=0.04, double ylabelsize=0.04, double yoffset=1.0) 
{
  hist->SetLineColor(markercolor);
  hist->SetMarkerColorAlpha(markercolor, markeralpha);
  hist->SetMarkerStyle(markerstyle);
  hist->SetMarkerSize(markersize);
  l->AddEntry(hist, text, "pl");

  //gPad->SetTickx(); 
  //gPad->SetTicky(); 
  // h->SetLineWidth(2);
  hist->GetYaxis()->SetTitleOffset(yoffset); //(1.05); 
  hist->GetYaxis()->SetTitleSize(ytitlesize); //(0.042); //the axis number labels
  hist->GetYaxis()->SetLabelSize(ylabelsize); //(0.042);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleOffset(xoffset);
  hist->GetXaxis()->SetTitleSize(xtitlesize); //(0.042);
  hist->GetXaxis()->SetLabelSize(xlabelsize); //(0.042);


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



//get histogram and clone it
THnSparse * getHistAndClone(TFile *f, std::string histname) {
  THnSparse *hsparsejet = (THnSparse*) f->Get(histname.c_str());
  THnSparse *hsparsejet_clone = (THnSparse *) hsparsejet->Clone("hsparsejet_c_clone"); //TODO: change this name!

  return hsparsejet_clone;
  
}


void applyCuts(THnSparse *hsparse, int pt_min, int pt_max, int d0_pt_cut, bool d0cuts=false) {
  hsparse->GetAxis(0)->SetRangeUser(pt_min, pt_max);
  if (d0cuts) {
      hsparse->GetAxis(1)->SetRangeUser(d0_pt_cut, pt_max); //d0_pt_cuts[i], pt_max); // apply cut on Dmeson pt
      hsparse->GetAxis(2)->SetRangeUser(-0.8, 0.8); // apply cut on Dmeson rapidity
  }
}

// get the observable histogram
//usually obsaxis is 3, but in new histograms it is 4. For jet level histograms, use 0.
TH1D * getObsHist(TFile *filename, std::string h_name, std::string h_jet_name, 
                int pt_min, int pt_max, int d0_pt_cut, std::string newhistname, 
                bool d0cuts=false, int obsaxis=3, bool ptrl=false, bool norm=true) {
  
  THnSparse *hsparse = getHistAndClone(filename, h_name);
  THnSparse *hsparse_jetlevel = getHistAndClone(filename, h_jet_name);

  applyCuts(hsparse, pt_min, pt_max, d0_pt_cut, d0cuts);
  applyCuts(hsparse_jetlevel, pt_min, pt_max, d0_pt_cut, d0cuts);
  
  TH1D *h_proj = hsparse->Projection(obsaxis);
  TH1D *h_proj_jetlevel = hsparse_jetlevel->Projection(0); // jet pt axis

  cout << "HNAME " << h_name << " num jets: " << h_proj_jetlevel->GetEntries() << endl;

  std::string hname = h_proj->GetName();
  hname += "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);
  h_proj->SetNameTitle(hname.c_str(), hname.c_str());

  // allow rebin or cloning here
  TH1D* hist = (TH1D*) h_proj->Clone(newhistname.c_str());

  //normalize
  double numjets = h_proj_jetlevel->Integral();
  if (norm) {
      cout << "getting normalized by # of jets here!!" << endl;
      hist->Scale(1/numjets, "width");
  }
  // cout << "There are " << h_proj->GetEntries() << " pair entries in this pt bin" << endl;
  // cout << "There are " << h_proj_jetlevel->GetEntries() << " jet entries in this pt bin" << endl;

  if (ptrl) {
      hist->GetXaxis()->SetTitle("#it{p}_{T}#it{R}_{L}");
  } else {
      hist->GetXaxis()->SetTitle("#it{R}_{L}");
  }
  // hist->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #times #frac{d#it{N}_{EEC}}{d#it{R}_{L}}");
  hist->GetYaxis()->SetTitle("#Sigma_{EEC}(#it{R}_{L})");

  return hist;
}

double drawWithLineAtMax(TH1D *hist, TLegend *l, std::string labeltext, int markercolor, int markerstyle, double lowx,
                  bool ptrl=false, bool removexaxis=false, double markeralpha=1, double xoffset=1.0, bool verbosity=false) {
  // FormatHist(l, hist, "gluon-init jets", markercolor_g, markerstyle_g);
  FormatHist(l, hist, labeltext, markercolor, markerstyle, markeralpha,
             0.06, 0.05, xoffset, 0.06, 0.05, 1.05);
  hist->GetXaxis()->SetRangeUser(lowx, 1.);
  hist->Draw("L same");
  if (removexaxis) hist->GetXaxis()->SetLabelSize(0);
  
  double hist_top_binpos = findTopOfCurve(hist, ptrl);
  drawVertLine(hist->GetBinCenter(hist_top_binpos), 0, hist->GetBinContent(hist_top_binpos), markercolor, 1)->Draw();

  if (verbosity) {
      cout << hist->GetName() << " top bincenter " << hist->GetBinCenter(hist_top_binpos) << ", " << hist->GetName() << " top binpos" << hist->GetBinContent(hist_top_binpos) << endl;
  }   

  return hist->GetBinCenter(hist_top_binpos);
}

void drawNoLine(TH1D *hist, TLegend *l, std::string labeltext, int markercolor, int markerstyle, double lowx,
                float markersize, bool removexaxis=false, double markeralpha=1, double xoffset=1.0, bool verbosity=false) {
  // FormatHist(l, hist, "gluon-init jets", markercolor_g, markerstyle_g);
  FormatHist(l, hist, labeltext, markercolor, markerstyle, markeralpha,
             0.06, 0.05, xoffset, 0.06, 0.05, 1.05, markersize);
  hist->GetXaxis()->SetRangeUser(lowx, 1.);
  hist->Draw("L same");
  if (removexaxis) hist->GetXaxis()->SetLabelSize(0);
  
  if (verbosity) {
      cout << "no prints here!" << endl;
  }   

  return;
}

double getPeak(TH1D *hist, bool ptrl=false, bool verbosity=false) {
  
  double hist_top_binpos = findTopOfCurve(hist, ptrl);
  if (verbosity) {
      cout << hist->GetName() << " top bincenter " << hist->GetBinCenter(hist_top_binpos) << ", " << hist->GetName() << " top binpos" << hist->GetBinContent(hist_top_binpos) << endl;
  }   

  return hist->GetBinCenter(hist_top_binpos);
}

double getPeakErr(TH1D *hist, bool ptrl=false, bool verbosity=false) {
  
  double hist_top_binpos = findTopOfCurve(hist, ptrl);
  if (verbosity) {
      cout << hist->GetName() << " top bincenter " << hist->GetBinCenter(hist_top_binpos) << ", " << hist->GetName() << " top binpos" << hist->GetBinContent(hist_top_binpos) << endl;
  }   

  return hist->GetBinWidth(hist_top_binpos)/2;
}

void plotGraph(TMultiGraph *mg, TLegend *l, int n_bins, double *pt_x, double *peaks_y, double *err_y, TString text, int markercolor, int markerstyle, double markeralpha=1.0) {
  // TGraph *g = new TGraph(n_bins,pt_x,peaks_y);
  double err_x[n_bins];
  for (int i=0; i<n_bins; i++) {
      err_x[i] = 0;
  }
  TGraphErrors *g = new TGraphErrors(n_bins,pt_x,peaks_y, err_x, err_y);
  g->SetMarkerColorAlpha(markercolor, markeralpha);
  g->SetMarkerStyle(markerstyle);
  g->SetMarkerSize(1.5);
  g->SetLineColorAlpha(markercolor, markeralpha);
  // if (markerstyle == kFullDiamond) g->SetMarkerSize(1.75);
  if (markercolor == kRed && markerstyle == kOpenCircle) g->SetMarkerSize(1);

  l->AddEntry(g, text, "pl");

  g->GetXaxis()->SetTitle("#it{p}_{T}");
  g->GetYaxis()->SetTitle("#it{R}_{L} peak position");
  g->SetMaximum(5.);

  if (text == "#frac{b-init full jets}{l-init full jets}") {
      cout << "here!!!!!" << endl;
      g->SetMaximum(5.);
  }
  
 
  // g->Draw("ap SAME");
  // return g;
  mg->Add(g);

}

// TPad * makeTopPad(TH1D *hdummy, double ymax, double botmarg=0.31) {
TPad * makeTopPad(double botmarg=0.31) {
  TPad *pad1 = new TPad("pad1","pad1",0.,0.,1.,1.);
  pad1->SetLogx();
  pad1->SetFillColor(0);
  pad1->SetFillStyle(0);
  pad1->SetTopMargin(0.025);
  pad1->SetBottomMargin(botmarg);
  pad1->Draw();
  pad1->cd();

  // hdummy->SetMaximum(ymax);
  // hdummy->Draw();

  return pad1;
}

TPad * makeBottomPad(double topmarg=0.71, double botmarg=0.975) {
  TPad *pad2 = new TPad("pad1","",0.,0.,1.,1.);
  pad2->SetTopMargin(topmarg); //0.71);
  pad2->SetBottomMargin(botmarg);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);
  pad2->Draw();
  pad2->SetLogx();
  // pad2->SetGridy();
  pad2->cd();

  return pad2;
}

TPad * plotRatio(TH1D *h1, TH1D *h2, std::string ratio_name, TLegend *l, int markercolor, int markerstyle, std::string yaxislabel, 
             double ymin=0.5, double ymax=1.5, double topmarg=0.71, double botmarg=0.975, double lowx=1e-3, std::string legendlabel="ratio",
             bool removexaxis=false, bool drawlineatone=false, bool justratio=false, int ndiv=5) {   

  TPad *pad2;
  if (!justratio) pad2 = makeBottomPad(topmarg, botmarg);

  TH1D* hratio = (TH1D*) h1->Clone(ratio_name.c_str());
  hratio->Divide(h2);
  hratio->SetMinimum(ymin);
  hratio->SetMaximum(ymax);
  
  //xtitlesize, xlabelsize, xoffset, ytitlesize, ylabelsize, yoffset
  FormatHist(l, hratio, legendlabel, markercolor, markerstyle, 1.0, 0.05, 0.04, 1.2, 0.035, 0.03, 1.5);
  hratio->GetYaxis()->SetNdivisions(5);
  if (removexaxis) hratio->GetXaxis()->SetLabelSize(0);
  else hratio->GetXaxis()->SetTitleOffset(0.95);
  
  hratio->GetXaxis()->SetRangeUser(lowx,1);

  hratio->Draw();

  hratio->GetYaxis()->SetTitleSize(0);

  if (!justratio) {
      double ypos = ((1-topmarg) + botmarg) / 2;
      TLatex *t = new TLatex(0.075,ypos,yaxislabel.c_str());
      t->SetTextAlign(22);
      t->SetTextColor(kBlack);
      t->SetTextFont(43);
      t->SetTextSize(14);
      // t->SetTextAngle(45);
      t->SetNDC(kTRUE);
      t->Draw();

      // draw line at 1
      
      if (drawlineatone) drawHoriLine(lowx, 1., 1., kGray+2, 1)->Draw();
      // drawHoriLine(1e-3, 1., 0.9, kGray+2, 6)->Draw();
      // drawHoriLine(1e-3, 1., 1.1, kGray+2, 6)->Draw();
  }

  return pad2;

}

// this function specifically plotting 2 ratios - one that is h1/h3, and one that is h2/h3
TPad * plotRatio2(TH1D *h1, TH1D *h2, TH1D *h3, std::string ratio_name1, std::string ratio_name2, 
             int markercolor1, int markerstyle1, double markersize1, double markeralpha1,
             int markercolor2, int markerstyle2, double markersize2, double markeralpha2, 
             std::string yaxislabel, double ymin=0.5, double ymax=1.5, 
             double topmarg=0.71, double botmarg=0.975, double lowx=1e-3, 
             std::string legendlabel1="ratio", std::string legendlabel2="ratio", 
             bool removexaxis=false, bool justratio=false, int ndiv=5) {   

  TPad *pad2;
  if (!justratio) pad2 = makeBottomPad(topmarg, botmarg);

  TLegend* l = new TLegend(0.50,0.300741,0.80,0.3585185,"");

  TH1D* hratio = (TH1D*) h1->Clone(ratio_name1.c_str());
  hratio->Divide(h3);
  hratio->SetMinimum(ymin);
  hratio->SetMaximum(ymax);

  TH1D* hratio2 = (TH1D*) h2->Clone(ratio_name2.c_str());
  hratio2->Divide(h3);
  
  //xtitlesize, xlabelsize, xoffset, ytitlesize, ylabelsize, yoffset                
  FormatHist(l, hratio, legendlabel1, markercolor1, markerstyle1, markeralpha1, 0.05, 0.04, 1.2, 0.035, 0.03, 1.05, markersize1);
  FormatHist(l, hratio2, legendlabel2, markercolor2, markerstyle2, markeralpha2, 0.05, 0.04, 1.2, 0.035, 0.03, 1.05, markersize2);
  hratio->GetYaxis()->SetNdivisions(5);
  hratio2->GetYaxis()->SetNdivisions(5);
  if (removexaxis) hratio->GetXaxis()->SetLabelSize(0);
  else hratio->GetXaxis()->SetTitleOffset(0.95);
  
  hratio->GetXaxis()->SetRangeUser(lowx,1);
  hratio2->GetXaxis()->SetRangeUser(lowx,1);

  hratio->Draw("same");
  hratio2->Draw("same");

  // hratio->GetYaxis()->SetTitleSize(0);
  hratio2->GetYaxis()->SetTitleSize(0);
  hratio->GetYaxis()->SetTitle(yaxislabel.c_str());
  hratio->GetYaxis()->SetTitleSize(0.06);
  hratio->GetYaxis()->SetTitleOffset(0.9);

  /*if (!justratio) {
      double ypos = ((1-topmarg) + botmarg) / 2;
      TLatex *t = new TLatex(0.075,ypos,yaxislabel.c_str());
      t->SetTextAlign(22);
      t->SetTextColor(kBlack);
      t->SetTextFont(43);
      t->SetTextSize(14);
      // t->SetTextAngle(45);
      t->SetNDC(kTRUE);
      t->Draw();

      // draw line at 1
      
      // drawHoriLine(1e-3, 1., 0.9, kGray+2, 6)->Draw();
      // drawHoriLine(1e-3, 1., 1.1, kGray+2, 6)->Draw();
  }
  */

  drawHoriLine(lowx, 1., 1., kGray+2, 1)->Draw();
  drawVertLine(0.4, ymin, ymax, kBlack, 2)->Draw();
  l->Draw("same");

  return pad2;

}


//===========================================================================//

void make_D0_softqcd_plots() {

//    gROOT->SetBatch(); //prevents plots from showing up
  gStyle->SetOptStat(0);
  SetStyle();
  Double_t markers[10] = {kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullStar, kOpenCircle, kOpenTriangleUp, kOpenDiamond, kOpenSquare, kOpenStar};
  Double_t marker_size = 1.5;
  Double_t colors[16] = {kRed, kGreen+2, kBlue, kRed+1, kGreen+1, kBlue+1, kRed+2, kGreen+2, kBlue+2, kRed+3, kGreen+3, kBlue+3, kOrange+1, kViolet+1, kYellow+1, kCyan+1};

  //=====================================================================//
  // File containing quark vs gluon histograms
  // the c/b enhanced means those files will only have c/b jets
  // most are charged jets unless otherwise specified

  const TString basedir = "/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/"; //"/Volumes/WORK USB/EEC/slurmfiles/"; //"/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/EEC/";
  const TString root_filename = "/AnalysisResultsFinal.root";

  const TString infile_D0_softqcd = basedir + "37158355" + root_filename;
  // const TString infile_c_enhanced_D0 = basedir + "25536401" + root_filename; //--replaceKP 1 --chinitscat 3

  
  // open the files
  TFile *f_D0_softqcd = new TFile(infile_D0_softqcd, "READ");
  // TFile *f_c_enhanced_D0 = new TFile(infile_c_enhanced_D0, "READ");

  // std::vector<TFile*> files;
  std::string add_name;

  add_name = ".pdf";
  cout << "output name will be " << add_name << endl;

  // Output directory
  std::string outdir= "plots/final/softqcd/"; //testing/";//"plots/test/";

  // Output file for binned results
  std::string outfile = outdir + "AnalysisResultsFinal_softqcd.root";
  TFile* f_out = new TFile(outfile.c_str(), "RECREATE");


  // Jet r value
  std::string jetR_list[] = { "0.4" };
  for (std::string jetR : jetR_list) {

      // Names of histograms in the file (quark, charm, gluon)
      const std::string hc_name = "h_EEC_JetPt_charm_R" + jetR;
      const std::string hl_name = "h_EEC_JetPt_light_R" + jetR;
      const std::string hg_name = "h_EEC_JetPt_gluon_R" + jetR;
      const std::string hi_name = "h_EEC_JetPt_inclusive_R" + jetR;
      const std::string hb_name = "h_EEC_JetPt_beauty_R" + jetR;
      const std::string hc_jet_name = "h_JetPt_charm_R" + jetR + "_jetlevel";
      const std::string hl_jet_name = "h_JetPt_light_R" + jetR + "_jetlevel";
      const std::string hg_jet_name = "h_JetPt_gluon_R" + jetR + "_jetlevel";
      const std::string hi_jet_name = "h_JetPt_inclusive_R" + jetR + "_jetlevel";
      const std::string hb_jet_name = "h_JetPt_beauty_R" + jetR + "_jetlevel";

      //-------------------------------------------------//
      //TODO: IF THNSPARSE CLONES ARE BEFORE THE PT BIN LOOP, MOVE THEM INTO IT!!
      //-------------------------------------------------//
      // Make canvases that save across pt bins
      TCanvas* c_D0 = new TCanvas();
      ProcessCanvas(c_D0);
      c_D0->cd();
      gPad->SetLogx();
      // gPad->SetLogy();


      //-------------------------------------------------//

      //const int pt_bins[] = { 10, 20, 40, 60, 80, 100, 150 };
      // const int pt_bins[] = { 20, 40 }; 
      const int pt_bins[] = {7, 10, 15, 30, 50, 70, 100, 150, 200 }; //, 100, 150 }; //{ 10, 20, 40 };
      const int d0_pt_cuts[] = { 3, 5, 5, 5, 5, 5, 5, 5 }; //, 5, 5 };
      const int n_bins = 4; //7;
      // double pt_x[n_bins] = { 20. }; //
      double pt_x[n_bins] = { 10., 15., 30., 50. }; //, 70., 100., 150., 200.0 };
      

      //Format color and style
      int markercolor_c = kRed; //charm
      int markerstyle_c = kFullCircle;
      int markercolor_g = kViolet+2; //gluon
      int markerstyle_g = 33; //diamond
      int markercolor_l = kGreen+2; //light
      int markerstyle_l = 21; //square
      int markercolor_i = kMagenta-6; //inclusive
      int markerstyle_i = 29; //star
      int markercolor_D0 = kOrange+7; //kGreen-5; //D0
      int markerstyle_D0 = kFullCircle;
      int markercolor_b = kAzure; //beauty
      int markerstyle_b = 34;

      int markerstyle_c_ch = kOpenCircle; //charged jets when being compared to full
      int markerfill_c_chn = 3944; //charged jets with neutral hadrons

      int markerstyle_l_ch = kOpenSquare; //charged jets when being compared to full
      int markerfill_b_chn = 3944; //charged jets with neutral hadrons

      
      for (int i = 0; i < n_bins; i++) {
          cout << "in pt bin" << i << endl;
          int pt_min = pt_bins[i];
          int pt_max = pt_bins[i+1];


          // define pt related variables
          TString ptbin = TString::Format("%d #leq #it{p}_{T}^{ch. jet} < %d GeV/#it{c}, #font[122]{|}#it{#eta}_{jet}#font[122]{|} #leq 0.5", pt_min, pt_max);
          TString ptD = TString::Format("%d #leq #it{p}_{T}^{D^{0}} < %d GeV/#it{c}, #font[122]{|}#it{y}_{D^{0}}#font[122]{|} #leq 0.8", d0_pt_cuts[i], pt_max);
          // if (plot_case == 3 or plot_case == 11 or plot_case == 12 or plot_case == 5) {
          //     ptbin = TString::Format("%d #leq #it{p}_{T}^{full jet} < %d GeV/#it{c}, #font[122]{|}#it{#eta}_{jet}#font[122]{|} #leq 0.5", pt_min, pt_max);
          // } else if (plot_case == 13 or plot_case == 21) {
          //     ptD = TString::Format("%d #leq #it{p}_{T}^{D^{0}} < %d GeV/#it{c}, #font[122]{|}#it{y}_{D^{0}}#font[122]{|} #leq 0.8", 5, pt_max);
          // }

          // make a canvas for each pt range
          TCanvas* c = new TCanvas();
          ProcessCanvas(c);
          c->cd();
          gPad->SetLogx();


          TLegend* l; // = new TLegend(0.17, 0.65, 0.5, 0.85);
          TLegend* l2 = new TLegend(0.1797168,0.400741,0.4562155,0.8885185,""); //dummy legend

          double maxy = 0;


          // Open histograms


          // TODO: figure out plot_case 13 and 21!!
          l = new TLegend(0.1957168,0.760741,0.462155,0.9505185,"");
          l->SetTextSize(0.028); //l->SetTextSize(0.045);
          l->SetBorderSize(0);
          l->SetFillStyle(0);
          l->AddEntry("NULL","PYTHIA 8 Monash 2013","h");
          l->AddEntry("NULL","pp, #sqrt{#it{s}} = 13 TeV","h");
          // l->AddEntry("NULL","anti-#it{k}_{T}, #it{R} = 0.4","h");
          // else {
          l->AddEntry("NULL","D^{0} #rightarrow K^{#minus} #pi^{+} and charge conj.","h");
          l->AddEntry("NULL","anti-#it{k}_{T} ch. jets, #it{R} = 0.4","h");
          // }
          l->AddEntry("NULL",ptbin,"h");
          
          // l->Draw("same");

          TLegend *l_main2 = new TLegend(0.1957168,0.550741,0.462155,0.755,""); //(0.17, 0.4, 0.5, 0.53);
          // l->SetTextSize(0.045);
          l_main2->SetTextSize(0.028);
          l_main2->SetBorderSize(0);
          l_main2->SetFillStyle(0);
          l_main2->AddEntry("NULL",ptD,"h");
          


          std::string pt_name = "_pt" + std::to_string(pt_min) + "-" + std::to_string(pt_max);


          //format: getObsHist(TFile *filename, std::string h_name, std::string h_jet_name, int pt_min, int pt_max, int d0_pt_cut, std::string newhistname, bool d0cuts=false, int obsaxis=3)
          TH1D *h_D0_softqcd = getObsHist(f_D0_softqcd, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_D0_softqcd" + pt_name, false, 4);
          
          // TH1D *h_c_enhanced_D0 = getObsHist(f_c_enhanced_D0, hc_name, hc_jet_name, pt_min, pt_max, d0_pt_cuts[i], "h_c_enhanced_D0" + pt_name, true, 4);


          

          // TCanvas *c1 = new TCanvas();
          // gPad->SetLogx();
          double lowx = 0.0;
          FormatHist(l, h_D0_softqcd, "softQCD", kRed, kFullDiamond, 1.0,
             0.06, 0.05, 1.0, 0.06, 0.05, 1.05, 1.5);
          h_D0_softqcd->Draw();
          // drawNoLine(h_D0_softqcd, l, "softQCD", kRed, kFullDiamond, lowx, 1.5, false, 1.0);

          // draw legend
          l->Draw("same");
          l_main2->Draw("same");


          std::string fname = outdir + "QG_comp" + pt_name + "_R" + jetR + add_name; //"_charmdecaysONcomparison.pdf"; // + "_normbytype.pdf"; //"_nonorm.pdf";
          const char* fnamec = fname.c_str();
          c->SaveAs(fnamec);
          delete c;

          f_out->cd();
          
          h_D0_softqcd->Write();
          // h_c_enhanced_D0->Write();
          

          delete h_D0_softqcd;
          // delete h_c_enhanced_D0;
          
           
      } // pT bins loop



  } // jetR loop

  

  // f_c_enhanced_D0->Close();

  // delete f_c_enhanced_D0;
  


  return;
}
