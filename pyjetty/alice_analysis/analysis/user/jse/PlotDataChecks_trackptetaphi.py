#!/usr/bin/env python3
"""Data-check plots (track pt eta phi) for the JSE analysis."""

import os
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)  # Turned off stats box for clean plots

# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------
SITE = "perlmutter"          # "perlmutter" or "hiccup"

period = "LHC26_thin_medium"
file_num = "260-264" #"263" #260-264
# period = "LHC24_ppref"
# file_num = "69"
input_filename = f"/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/rootfiles/data_checks/{period}/{file_num}/HistsDataCheck.root"
output_dir = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots/data_checks/track_ptetaphi"
output_ext = f"{period}_{file_num}"
os.makedirs(output_dir, exist_ok=True)

f = ROOT.TFile.Open(input_filename, "READ")
h3d = f.Get("h_track_pt_eta_phi")
h_pt = f.Get("h_pt")

# Canvas setup helper
def create_canvas():
    c = ROOT.TCanvas("c", "c", 800, 600)
    c.SetGrid()
    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.12)
    c.SetBottomMargin(0.12)
    return c

def normalize_hist(h, y_title):
    """Normalizes a 1D histogram by total counts N_trk and bin width dx."""
    n_trk = h.Integral()
    if n_trk > 0:
        h.Scale(1.0, "width")  # Divide by bin width dx
        h.Scale(1.0 / n_trk)   # Divide by total track yield N_trk
    h.SetYTitle(y_title)
    h.SetMarkerStyle(20)
    h.SetMarkerSize(0.9)
    h.SetMarkerColor(ROOT.kBlack)
    h.SetLineColor(ROOT.kBlack)

# Define ROOT LaTeX strings as raw plain strings to avoid Python f-string escaping bugs
y_title_pt = r"#frac{1}{N_{trk}} #frac{dN_{trk}}{dp_{T}} [(GeV/c)^{-1}]"
y_title_eta = r"#frac{1}{N_{trk}} #frac{dN_{trk}}{d#eta}"
y_title_phi = r"#frac{1}{N_{trk}} #frac{dN_{trk}}{d#phi} [rad^{-1}]"
y_title_2d = r"#frac{1}{N_{trk}} #frac{d^{2}N_{trk}}{d#eta d#phi}"

# ---------------------------------------------------------------------------
# 1. Standalone h_pt
# ---------------------------------------------------------------------------
c = create_canvas()
c.SetLogy()
h_pt_norm = h_pt.Clone("h_pt_standalone_norm")
normalize_hist(h_pt_norm, y_title_pt)
h_pt_norm.SetTitle("Standalone p_{T} Spectrum;p_{T} [GeV/c];" + y_title_pt)
h_pt_norm.Draw("E1 P")
c.SaveAs(os.path.join(output_dir, f"{output_ext}_h_pt_standalone.pdf"))

# ---------------------------------------------------------------------------
# 2. Full Spectrums
# ---------------------------------------------------------------------------
# Full pT
c = create_canvas()
c.SetLogy()
h_pt_full = h3d.ProjectionX("h_pt_full")
normalize_hist(h_pt_full, y_title_pt)
h_pt_full.SetTitle("Full p_{T} Spectrum;p_{T} [GeV/c];" + y_title_pt)
h_pt_full.Draw("E1 P")
c.SaveAs(os.path.join(output_dir, f"{output_ext}_full_pt.pdf"))

# Full Eta
c = create_canvas()
c.SetLogy(False)
h_eta_full = h3d.ProjectionY("h_eta_full")
normalize_hist(h_eta_full, y_title_eta)
h_eta_full.SetTitle("Full #eta Spectrum;#eta;" + y_title_eta)
h_eta_full.Draw("E1 P")
c.SaveAs(os.path.join(output_dir, f"{output_ext}_full_eta.pdf"))

# Full Phi
c = create_canvas()
h_phi_full = h3d.ProjectionZ("h_phi_full")
normalize_hist(h_phi_full, y_title_phi)
h_phi_full.SetTitle("Full #phi Spectrum;#phi [rad];" + y_title_phi)
h_phi_full.Draw("E1 P")
c.SaveAs(os.path.join(output_dir, f"{output_ext}_full_phi.pdf"))

# ---------------------------------------------------------------------------
# 3. Binned Projections in pT
# ---------------------------------------------------------------------------
pt_bins = [0, 2, 5, 10, 20, 40, 60, 80, 100, 150, 200]
pt_axis = h3d.GetXaxis()

# Collections to store histograms for uniform Y-axis scaling
eta_hists = []
phi_hists = []
phi_rebin_hists = []
phi_eta_hists = []
bin_tags = []
title_suffixes = []

# Pass 1: Extract, normalize, and measure max Y values
max_y_eta = 0.0
max_y_phi = 0.0
max_y_phi_rebin = 0.0

for i in range(len(pt_bins) - 1):
    bin_min = pt_axis.FindBin(pt_bins[i])
    bin_max = pt_axis.FindBin(pt_bins[i+1]) - 1
    if bin_max < bin_min:
        bin_max = bin_min
        
    tag = f"pt_{pt_bins[i]}_{pt_bins[i+1]}"
    title_suffix = f" ({pt_bins[i]} #leq p_{{T}} < {pt_bins[i+1]} GeV/c)"
    
    bin_tags.append(tag)
    title_suffixes.append(title_suffix)
    
    # Eta projection
    h_eta = h3d.ProjectionY(f"h_eta_{tag}", bin_min, bin_max, 0, -1)
    normalize_hist(h_eta, y_title_eta)
    max_y_eta = max(max_y_eta, h_eta.GetMaximum())
    eta_hists.append(h_eta)
    
    # Phi projection
    h_phi = h3d.ProjectionZ(f"h_phi_{tag}", bin_min, bin_max, 0, -1)
    normalize_hist(h_phi, y_title_phi)
    max_y_phi = max(max_y_phi, h_phi.GetMaximum())
    phi_hists.append(h_phi)
    
    # Phi projection (rebinned x8)
    h_phi_rebin = h3d.ProjectionZ(f"h_phi_rebinned_raw_{tag}", bin_min, bin_max, 0, -1)
    h_phi_rebin.Rebin(8)
    normalize_hist(h_phi_rebin, y_title_phi)
    max_y_phi_rebin = max(max_y_phi_rebin, h_phi_rebin.GetMaximum())
    phi_rebin_hists.append(h_phi_rebin)
    
    # 2D projection
    pt_axis.SetRange(bin_min, bin_max)
    h_phi_eta = h3d.Project3D(f"zy_h_phi_eta_{tag}")
    n_2d = h_phi_eta.Integral()
    if n_2d > 0:
        h_phi_eta.Scale(1.0, "width")
        h_phi_eta.Scale(1.0 / n_2d)
    phi_eta_hists.append(h_phi_eta)
    pt_axis.SetRange(0, 0)

# Pass 2: Set global Y limits and draw/save
for i in range(len(eta_hists)):
    tag = bin_tags[i]
    title_suffix = title_suffixes[i]
    
    # 3a. Eta Plot
    h_eta = eta_hists[i]
    h_eta.SetMaximum(1.2 * max_y_eta)
    h_eta.SetMinimum(0.0)
    h_eta.SetTitle("#eta Spectrum" + title_suffix + ";#eta;" + y_title_eta)
    c = create_canvas()
    h_eta.Draw("E1 P")
    c.SaveAs(os.path.join(output_dir, f"{output_ext}_eta_{tag}.pdf"))
    
    # 3b. Phi Plot
    h_phi = phi_hists[i]
    h_phi.SetMaximum(1.2 * max_y_phi)
    h_phi.SetMinimum(0.0)
    h_phi.SetTitle("#phi Spectrum" + title_suffix + ";#phi [rad];" + y_title_phi)
    c = create_canvas()
    h_phi.Draw("E1 P")
    c.SaveAs(os.path.join(output_dir, f"{output_ext}_phi_{tag}.pdf"))
    
    # 3c. Rebinned Phi Plot
    h_phi_rebin = phi_rebin_hists[i]
    h_phi_rebin.SetMaximum(1.2 * max_y_phi_rebin)
    h_phi_rebin.SetMinimum(0.0)
    h_phi_rebin.SetTitle("#phi Spectrum (Rebinned x4)" + title_suffix + ";#phi [rad];" + y_title_phi)
    c_rebin = create_canvas()
    h_phi_rebin.Draw("E1 P")
    c_rebin.SaveAs(os.path.join(output_dir, f"{output_ext}_phi_rebinned_{tag}.pdf"))
    
    # 3d. 2D Phi vs Eta Plot
    h_phi_eta = phi_eta_hists[i]
    h_phi_eta.SetTitle("#phi vs #eta" + title_suffix + ";#eta;#phi [rad];" + y_title_2d)
    c = create_canvas()
    h_phi_eta.Draw("COLZ")
    c.SaveAs(os.path.join(output_dir, f"{output_ext}_phi_vs_eta_{tag}.pdf"))

f.Close()
print("All plots successfully generated with uniform Y-axis limits.")