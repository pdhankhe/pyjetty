import ROOT

jobid="53567700"
add_ext = False

# f = ROOT.TFile(f"/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data_checks/{jobid}/MergedHistsDataCheck.root")
f = ROOT.TFile(f"/global/cfs/projectdirs/alice/blianggi/mypyjetty/analysis/testing/HistsDataCheck.root")

if add_ext:
    f = ROOT.TFile(f"/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data_checks/{jobid}/1/HistsDataCheckExt.root")
    f_2 = ROOT.TFile(f"/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data_checks/{jobid}/2/HistsDataCheckExt.root")
    f_luisa = ROOT.TFile(f"/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data_checks/{jobid}/1/trackPt_QAhist.root")
    f_luisa_2 = ROOT.TFile(f"/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data_checks/{jobid}/2/trackPt_QAhist_file2.root")
f.h_cuts.Print("all")     # see total/sel8/global-track counts

c = ROOT.TCanvas("c", "", 1500, 450)
c.Divide(3, 1)
c.cd(1); ROOT.gPad.SetLogy(); f.h_pt.Draw()
c.cd(2); f.h_eta.Draw()
c.cd(3); f.h_phi.Draw()
if add_ext:
    c.SaveAs("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots/data_checks/merged_1_ext.pdf")
else:
    c.SaveAs("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots/data_checks/merged_1.pdf")


# ============================================================
# Two-panel canvas: all-event track pT (left) vs
# jet-constituent / triggered-event track pT (right)
# ============================================================
ROOT.gStyle.SetOptStat(0)

# y-projection of the 2D over jet pT >= 8 GeV/c
xlo = f.h_trackpt_jetpt.GetXaxis().FindBin(8.0 + 1e-6)
xhi = f.h_trackpt_jetpt.GetNbinsX()
h_injet = f.h_trackpt_jetpt.ProjectionY("h_injet", xlo, xhi)
h_injet.SetDirectory(0)   # detach so it survives the file going out of scope

c_pt = ROOT.TCanvas("c_pt", "", 1200, 600)
c_pt.Divide(2, 1)

# ----- Left panel: all tracks in all events -----
c_pt.cd(1)
ROOT.gPad.SetLogy()
ROOT.gPad.SetLeftMargin(0.13)

f.h_pt.SetLineColor(ROOT.kBlack)
f.h_pt.SetLineWidth(2)
f.h_pt.SetTitle(";p_{T} (GeV/c);counts")
f.h_pt.Draw("HIST")

legL = ROOT.TLegend(0.40, 0.78, 0.88, 0.88)
legL.SetBorderSize(0)
legL.SetFillStyle(0)
legL.AddEntry(f.h_pt, "All tracks in all events", "l")
legL.Draw()

# ----- Right panel: jet constituents vs triggered-event tracks (+ ratio) -----
c_pt.cd(2)
pad_R = ROOT.gPad

# main spectra pad (top)
p_top = ROOT.TPad("p_top_R", "", 0.0, 0.32, 1.0, 1.0)
p_top.SetLogy()
p_top.SetLeftMargin(0.13)
p_top.SetBottomMargin(0.02)
p_top.Draw()
p_top.cd()

h_injet.SetLineColor(ROOT.kRed + 1)
h_injet.SetLineWidth(2)
f.h_pt_jettrig.SetLineColor(ROOT.kBlue + 1)
f.h_pt_jettrig.SetLineWidth(2)

# draw the taller one first so nothing gets clipped
if f.h_pt_jettrig.GetMaximum() >= h_injet.GetMaximum():
    lead, sub = f.h_pt_jettrig, h_injet
else:
    lead, sub = h_injet, f.h_pt_jettrig

lead.SetStats(0)
lead.SetTitle(";;counts")
lead.GetYaxis().SetTitleSize(0.05)
lead.GetYaxis().SetTitleOffset(1.1)
lead.GetXaxis().SetLabelSize(0)          # hide x labels on top pad
lead.Draw("HIST")
sub.Draw("HIST SAME")

legR = ROOT.TLegend(0.30, 0.75, 0.88, 0.88)
legR.SetBorderSize(0)
legR.SetFillStyle(0)
legR.AddEntry(h_injet,         "Tracks in all jets > 8 GeV/c", "l")
legR.AddEntry(f.h_pt_jettrig,  "All tracks in events with a jet > 8 GeV/c", "l")
legR.Draw()

# ratio pad (bottom)
pad_R.cd()
p_bot = ROOT.TPad("p_bot_R", "", 0.0, 0.0, 1.0, 0.32)
p_bot.SetLeftMargin(0.13)
p_bot.SetTopMargin(0.02)
p_bot.SetBottomMargin(0.3)
p_bot.Draw()
p_bot.cd()

h_ratio = h_injet.Clone("h_ratio_R")
h_ratio.SetStats(0)
h_ratio.SetLineColor(ROOT.kBlack)
h_ratio.SetLineWidth(2)
h_ratio.Divide(f.h_pt_jettrig)           # in-jet / all-event-tracks

h_ratio.SetTitle(";p_{T} (GeV/c);in-jet / all")
h_ratio.GetYaxis().SetNdivisions(505)
h_ratio.GetYaxis().SetTitleSize(0.11)
h_ratio.GetYaxis().SetTitleOffset(0.45)
h_ratio.GetYaxis().SetLabelSize(0.09)
h_ratio.GetXaxis().SetTitleSize(0.11)
h_ratio.GetXaxis().SetTitleOffset(1.0)
h_ratio.GetXaxis().SetLabelSize(0.09)
# h_ratio.GetYaxis().SetRangeUser(0, 1.1)   # uncomment to fix ratio range
h_ratio.Draw("HIST")

# reference line at 1
line = ROOT.TLine(h_ratio.GetXaxis().GetXmin(), 1.0, h_ratio.GetXaxis().GetXmax(), 1.0)
line.SetLineStyle(2)
line.SetLineColor(ROOT.kGray + 2)
line.Draw()

c_pt.Update()
c_pt.SaveAs("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots/data_checks/pt_panels_1.pdf")


# === ETA ===
# ============================================================
# Two-panel canvas for eta (mirrors the pT version)
# ============================================================
h_injet_eta = f.h_tracketa_jetpt.ProjectionY("h_injet_eta", xlo, xhi)
h_injet_eta.SetDirectory(0)

c_eta = ROOT.TCanvas("c_eta", "", 1200, 600)
c_eta.Divide(2, 1)

# ----- Left: all tracks in all events -----
c_eta.cd(1)
ROOT.gPad.SetLeftMargin(0.13)
f.h_eta.SetLineColor(ROOT.kBlack)
f.h_eta.SetLineWidth(2)
f.h_eta.SetTitle(";#eta;counts")
f.h_eta.Draw("HIST")

legL_eta = ROOT.TLegend(0.30, 0.18, 0.78, 0.28)
legL_eta.SetBorderSize(0); legL_eta.SetFillStyle(0)
legL_eta.AddEntry(f.h_eta, "All tracks in all events", "l")
legL_eta.Draw()

# ----- Right: jet constituents vs triggered-event tracks -----
c_eta.cd(2)
ROOT.gPad.SetLeftMargin(0.13)
h_injet_eta.SetLineColor(ROOT.kRed + 1)
h_injet_eta.SetLineWidth(2)
h_injet_eta.SetTitle(";#eta;counts")

f.h_eta_jettrig.SetLineColor(ROOT.kBlue + 1)
f.h_eta_jettrig.SetLineWidth(2)

if f.h_eta_jettrig.GetMaximum() >= h_injet_eta.GetMaximum():
    f.h_eta_jettrig.SetTitle(";#eta;counts")
    f.h_eta_jettrig.Draw("HIST")
    h_injet_eta.Draw("HIST SAME")
else:
    h_injet_eta.Draw("HIST")
    f.h_eta_jettrig.Draw("HIST SAME")

legR_eta = ROOT.TLegend(0.25, 0.18, 0.80, 0.32)
legR_eta.SetBorderSize(0); legR_eta.SetFillStyle(0)
legR_eta.AddEntry(h_injet_eta,      "Tracks in all jets > 8 GeV/c", "l")
legR_eta.AddEntry(f.h_eta_jettrig,  "All tracks in events with a jet > 8 GeV/c", "l")
legR_eta.Draw()

c_eta.Update()
c_eta.SaveAs("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots/data_checks/eta_panels_1.pdf")


# ============================================================
# track eta vs track pt: fiducial constituents vs edge (>8 GeV) constituents
# ============================================================
c_2d = ROOT.TCanvas("c_2d", "", 1300, 600)
c_2d.Divide(2, 1)

c_2d.cd(1)
ROOT.gPad.SetRightMargin(0.15); ROOT.gPad.SetLogz()
f.h_tracketa_trackpt_fid.SetTitle("Constituents of fiducial jets;p_{T}^{track} [GeV/c];#eta^{track}")
f.h_tracketa_trackpt_fid.Draw("COLZ")

c_2d.cd(2)
ROOT.gPad.SetRightMargin(0.15); ROOT.gPad.SetLogz()
f.h_tracketa_trackpt_edge.SetTitle("Constituents of jets>8 GeV outside fiducial;p_{T}^{track} [GeV/c];#eta^{track}")
f.h_tracketa_trackpt_edge.Draw("COLZ")

c_2d.SaveAs("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots/data_checks/tracketa_trackpt_panels_1.pdf")


if add_ext:
    # Dedicated canvas that plots ONLY the variable-binned h_pt_ext
    c_ext = ROOT.TCanvas("c_ext", "", 600, 500)
    c_ext.SetLogy()
    # Divide by bin width so the variable bins are shown as a proper density
    h_pt_ext_norm = f.h_pt_ext.Clone("h_pt_ext_norm")
    h_pt_ext_norm.SetDirectory(0)          # detach from file so it survives
    h_pt_ext_norm.Scale(1.0, "width")
    h_pt_ext_norm.GetYaxis().SetTitle("counts / (GeV/c)")
    h_pt_ext_norm.Draw("hist e0")
    c_ext.SaveAs("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots/data_checks/h_pt_ext_1.pdf")

    # Dedicated canvas to compare me+Luisa
    c_comp = ROOT.TCanvas("c_comp", "", 600, 500)
    c_comp.SetLogy()

    f.h_pt.SetMarkerColorAlpha(ROOT.kBlue, 0.75)
    f.h_pt.SetLineColorAlpha(ROOT.kBlue, 0.75)
    f.h_pt.Draw("hist e0")
    f.h_pt_float.SetMarkerColorAlpha(ROOT.kGreen+2, 0.75)
    f.h_pt_float.SetLineColorAlpha(ROOT.kGreen+2, 0.75)
    f.h_pt_float.SetLineStyle(ROOT.kDashed)
    f.h_pt_float.Draw("hist e0 same")
    f_luisa.hQA_trackPt.SetMarkerColorAlpha(ROOT.kRed, 0.75)
    f_luisa.hQA_trackPt.SetLineColorAlpha(ROOT.kRed, 0.75)
    f_luisa.hQA_trackPt.SetLineStyle(3)
    f_luisa.hQA_trackPt.Draw("hist e0 same")

    # Legend: (x1, y1, x2, y2) in normalized pad coordinates (0-1)
    leg = ROOT.TLegend(0.6, 0.75, 0.85, 0.85)
    leg.SetBorderSize(0)          # no box border
    leg.SetFillStyle(0)           # transparent background
    leg.AddEntry(f.h_pt,              "Beatrice double",    "l")
    leg.AddEntry(f.h_pt_float,        "Beatrice float",    "l")
    leg.AddEntry(f_luisa.hQA_trackPt, "Luisa", "l")
    leg.Draw()

    c_comp.SaveAs("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots/data_checks/h_pt_comp_1.pdf")
    
    # Dedicated canvas to compare me+Luisa
    c_comp_2 = ROOT.TCanvas("c_comp_2", "", 600, 500)
    c_comp_2.SetLogy()

    f_2.h_pt.SetMarkerColorAlpha(ROOT.kBlue, 0.75)
    f_2.h_pt.SetLineColorAlpha(ROOT.kBlue, 0.75)
    f_2.h_pt.Draw("hist e0")
    f_luisa_2.hQA_trackPt.SetMarkerColorAlpha(ROOT.kRed, 0.75)
    f_luisa_2.hQA_trackPt.SetLineColorAlpha(ROOT.kRed, 0.75)
    f_luisa_2.hQA_trackPt.SetLineStyle(ROOT.kDashed)
    f_luisa_2.hQA_trackPt.Draw("hist e0 same")

    # Legend: (x1, y1, x2, y2) in normalized pad coordinates (0-1)
    leg_2 = ROOT.TLegend(0.6, 0.75, 0.85, 0.85)
    leg_2.SetBorderSize(0)          # no box border
    leg_2.SetFillStyle(0)           # transparent background
    leg_2.AddEntry(f_2.h_pt,              "Beatrice file 2",    "l")
    leg_2.AddEntry(f_luisa_2.hQA_trackPt, "Luisa file 2", "l")
    leg_2.Draw()

    c_comp_2.SaveAs("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots/data_checks/h_pt_comp_2.pdf")