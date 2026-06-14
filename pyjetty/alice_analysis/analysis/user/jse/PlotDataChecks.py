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
ylo = f.h_trackpt_jetpt.GetYaxis().FindBin(8.0 + 1e-6)
yhi = f.h_trackpt_jetpt.GetNbinsY()
h_injet = f.h_trackpt_jetpt.ProjectionY("h_injet", ylo, yhi)
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

# ----- Right panel: jet constituents vs triggered-event tracks -----
c_pt.cd(2)
ROOT.gPad.SetLogy()
ROOT.gPad.SetLeftMargin(0.13)

h_injet.SetLineColor(ROOT.kRed + 1)
h_injet.SetLineWidth(2)
h_injet.SetTitle(";p_{T} (GeV/c);counts")

f.h_pt_jettrig.SetLineColor(ROOT.kBlue + 1)
f.h_pt_jettrig.SetLineWidth(2)

# draw the taller one first so nothing gets clipped
if f.h_pt_jettrig.GetMaximum() >= h_injet.GetMaximum():
    f.h_pt_jettrig.SetTitle(";p_{T} (GeV/c);counts")
    f.h_pt_jettrig.Draw("HIST")
    h_injet.Draw("HIST SAME")
else:
    h_injet.Draw("HIST")
    f.h_pt_jettrig.Draw("HIST SAME")

legR = ROOT.TLegend(0.30, 0.75, 0.88, 0.88)
legR.SetBorderSize(0)
legR.SetFillStyle(0)
legR.AddEntry(h_injet,         "Tracks in all jets > 8 GeV/c", "l")
legR.AddEntry(f.h_pt_jettrig,  "All tracks in events with a jet > 8 GeV/c", "l")
legR.Draw()

c_pt.Update()
c_pt.SaveAs("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots/data_checks/pt_panels_1.pdf")


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