#include <stdio.h>
#include "fitfunctions.h"


TF1 * fit_histogram_linearfit(TH1D * hist, std::string observable, std::string addname, double fitmax) {

    /*
    // could also do:
    TF1 *linearFit = new TF1("linearFit", "pol1", 2.0, 8.0); // Fit range: [2, 8]
    hist->Fit(linearFit, "R"); // "R" = ensures the fit uses the specified range
    */

    // Perform the linear fit
    // hist->Fit("pol1", "R"); //, "Q"); // "pol1" = linear function, "Q" = quiet mode
    // TF1 *fitFunction = hist->GetFunction("pol1");
    // auto fitResult = hist->Fit("pol1", "S");

    TF1 *fitFunction = new TF1("fitFunction", "pol1", 0.0, fitmax); // Fit range: [2, 8]
    fitFunction->SetParameter(1,0.005); // set initial slope parameter to 0.005
    auto fitResult = hist->Fit(fitFunction, "SR");
    // "W": Ignore weights - set the weights of all non-zero bins to 1
    // "E": Perform better error estimation.
    // "S": The full result of the fit is returned in the TFitResultPtr (incl cov matrix)
    // 7.1.1. in https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html
    
    // Retrieve fit parameters
    double slope = fitFunction->GetParameter(1);  // Slope of the line
    double intercept = fitFunction->GetParameter(0); // Intercept of the line
    double chi2 = fitResult->Chi2(); //or fitFunction->GetChisquare();
    int ndf = fitResult->Ndf(); //or fitFunction->GetNDF();
    double pvalue = fitFunction->GetProb();

    // std::cout << "Fit Results: slope = " << slope << ", intercept = " << intercept << std::endl;
    // std::cout << "Chi2/Ndf = " << chi2 / ndf << std::endl;

    

    // Create a TPaveText (x1, y1, x2, y2 in NDC coordinates)
    TPaveText* pavetext = new TPaveText(0.2, 0.6, 0.5, 0.8, "NDC"); // Coordinates in normalized device space
    // pavetext->SetFillColor(0);         // Set background color (0 for transparent)
    // pavetext->SetTextColor(1);         // Set text color
    // pavetext->SetTextSize(0.03);       // Set text size
    // pavetext->SetBorderSize(1);        // Set border size

    // Add multiple lines
    pavetext->AddText(Form("Fit results: y = %.3fx + %.3f", slope, intercept));
    pavetext->AddText(Form("#Chi^{2}/Ndf = %.4f", chi2 / ndf));
    pavetext->AddText(Form("p-value = %.4f", pvalue));

    // Draw the histogram and the fit
    TCanvas *c1 = new TCanvas("c1", "Linear Fit", 800, 600);
    TFile *dummy_file;
    plot_and_save_one_histogram(c1, dummy_file, observable, hist, 4, addname, pavetext);

    return fitFunction;
}

TF1 * fit_histogram_quadfit(TH1D * hist, std::string observable, std::string addname, double fitmax) {

    TF1 *fitFunction = new TF1("fitFunction", "pol2", 0.0, fitmax); // Fit range: [2, 8]
    auto fitResult = hist->Fit(fitFunction, "SR");

    // Retrieve fit parameters
    double p0 = fitFunction->GetParameter(0); // a in ax^2 + bx + c
    double p1 = fitFunction->GetParameter(1); // b in ax^2 + bx + c
    double p2 = fitFunction->GetParameter(2); // b in ax^2 + bx + c
    double chi2 = fitResult->Chi2(); //or fitFunction->GetChisquare();
    int ndf = fitResult->Ndf(); //or fitFunction->GetNDF();
    double pvalue = fitFunction->GetProb();

    // Create a TPaveText (x1, y1, x2, y2 in NDC coordinates)
    TPaveText* pavetext = new TPaveText(0.2, 0.6, 0.5, 0.8, "NDC"); // Coordinates in normalized device space

    // Add multiple lines
    pavetext->AddText(Form("Fit results: y = %.5fx^{2} + %.3fx + %.3f", p2, p1, p0));
    pavetext->AddText(Form("#Chi^{2}/Ndf = %.4f", chi2 / ndf));
    pavetext->AddText(Form("p-value = %.4f", pvalue));

    // Draw the histogram and the fit
    TCanvas *c1 = new TCanvas("c1", "Quadratic Fit", 800, 600);
    TFile *dummy_file;
    plot_and_save_one_histogram(c1, dummy_file, observable, hist, 5, addname, pavetext);

    return fitFunction;

}

TF1 * fit_histogram_expofit(TH1D * hist, std::string observable, std::string addname, double fitmax) {

    TF1 *fitFunction = new TF1("fitFunction", "expo", 0.0, fitmax); // expo: f(x) = exp(p0+p1*x)
    auto fitResult = hist->Fit(fitFunction, "SR");

    // Retrieve fit parameters
    double p0 = fitFunction->GetParameter(0);  // a in ax^2 + bx + c
    double p1 = fitFunction->GetParameter(1); // b in ax^2 + bx + c
    double chi2 = fitResult->Chi2(); //or fitFunction->GetChisquare();
    int ndf = fitResult->Ndf(); //or fitFunction->GetNDF();
    double pvalue = fitFunction->GetProb();

    // Create a TPaveText (x1, y1, x2, y2 in NDC coordinates)
    TPaveText* pavetext = new TPaveText(0.2, 0.6, 0.5, 0.8, "NDC"); // Coordinates in normalized device space

    // Add multiple lines
    pavetext->AddText(Form("Fit results: y = exp(%.3fx + %.3f)", p1, p0));
    pavetext->AddText(Form("#Chi^{2}/Ndf = %.4f", chi2 / ndf));
    pavetext->AddText(Form("p-value = %.4f", pvalue));

    // Draw the histogram and the fit
    TCanvas *c1 = new TCanvas("c1", "Exponential Fit", 800, 600);
    TFile *dummy_file;
    plot_and_save_one_histogram(c1, dummy_file, observable, hist, 6, addname, pavetext);

    return fitFunction;
}

void fitfunctions()
{
   printf("Hello World\n");
}