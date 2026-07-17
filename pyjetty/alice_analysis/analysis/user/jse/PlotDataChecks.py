import ROOT

jobid="53567700"
add_ext = False

# # PERLMUTTER FILEPATHS
input_base = f"/global/cfs/projectdirs/alice/blianggi/mypyjetty/analysis/testing" # perlmutter
extra_input_base = f"/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data_checks/{jobid}" # perlmutter
output_base = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots/data_checks" # perlmutter

# # HICCUP FILEPATHS
# input_base = "/software/users/blianggi/mypyjetty/analysis/testing" # hiccup
# extra_input_base = "" # no extra_input_base on hiccup
# output_base = "/software/users/blianggi/mypyjetty/storage/jse/plots/data_checks" # hiccup

# f = ROOT.TFile(f"/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data_checks/{jobid}/MergedHistsDataCheck.root")
f = ROOT.TFile(f"{input_base}/HistsDataCheck.root")
print("using file:", f.GetName())

if add_ext:
    f = ROOT.TFile(f"{extra_input_base}/1/HistsDataCheckExt.root")
    f_2 = ROOT.TFile(f"{extra_input_base}/2/HistsDataCheckExt.root")
    f_luisa = ROOT.TFile(f"{extra_input_base}/1/trackPt_QAhist.root")
    f_luisa_2 = ROOT.TFile(f"{extra_input_base}/2/trackPt_QAhist_file2.root")
f.h_cuts.Print("all")     # see total/sel8/global-track counts

c = ROOT.TCanvas("c", "", 1500, 450)
c.Divide(3, 1)
c.cd(1); ROOT.gPad.SetLogy(); f.h_pt.Draw()
c.cd(2); f.h_eta.Draw()
c.cd(3); f.h_phi.Draw()
if add_ext:
    c.SaveAs(f"{output_base}/merged_1_ext.pdf")
else:
    c.SaveAs(f"{output_base}/merged_1.pdf")


# ============================================================
# Two-panel canvas: all-event track pT (left) vs
# jet-constituent / triggered-event track pT (right)
# ============================================================
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadGridX(1)
ROOT.gStyle.SetPadGridY(1)

# y-projection of the 2D over jet pT >= 8 GeV/c
xlo = f.h_trackpt_jetpt.GetXaxis().FindBin(8.0 + 1e-6)
xhi = f.h_trackpt_jetpt.GetXaxis().FindBin(20.0+1e-6) #GetNbinsX()
print("Getnbinsx():", f.h_trackpt_jetpt.GetXaxis().GetNbins())
h_injet = f.h_trackpt_jetpt.ProjectionY("h_injet", xlo, xhi)
h_injet.SetDirectory(0)   # detach so it survives the file going out of scope

# ----- build the edge-jet (>8 GeV/c) track-pt spectrum -----  # h_trackpt_jetpt_edge axes: x = jet pT, y = track pT
xax = f.h_trackpt_jetpt_edge.GetXaxis()
bin_lo = xax.FindBin(8.0 + 1e-6)     # first bin with jet pT > 8
bin_hi = xax.FindBin(20.0 + 1e-6) #xax.GetNbins() + 1                  # include overflow
h_injet_edge = f.h_trackpt_jetpt_edge.ProjectionY("h_injet_edge", bin_lo, bin_hi)

# h_pt_R06_trig: all tracks in events with an R=0.6 jet > 8 GeV/c
f.h_pt_R06_trig.SetLineColor(ROOT.kOrange + 7)
f.h_pt_R06_trig.SetLineWidth(2)

c_pt = ROOT.TCanvas("c_pt", "", 1200, 600)
c_pt.Divide(2, 1)

# ----- Left panel: all tracks in all events -----
c_pt.cd(1)
ROOT.gPad.SetLogy()
ROOT.gPad.SetLeftMargin(0.13)

f.h_pt.SetLineColor(ROOT.kBlack)
f.h_pt.SetLineWidth(1)
f.h_pt.SetTitle(";p_{T} (GeV/c);counts")
f.h_pt.Draw("HIST")
f.h_pt_R06_trig_fid.SetLineColor(ROOT.kOrange + 2)
f.h_pt_R06_trig_fid.SetLineWidth(1)
f.h_pt_R06_trig_fid.Draw("HIST SAME")

legL = ROOT.TLegend(0.40, 0.82, 0.88, 0.88)
legL.SetBorderSize(0)
legL.SetFillStyle(0)
legL.AddEntry(f.h_pt, "All tracks in all events", "l")
legL.AddEntry(f.h_pt_R06_trig_fid, "Tracks in events with R=0.6 jet > 8 GeV/c, jet |eta| < 0.3", "l")
legL.Draw()

# ----- Right side: jet constituents vs triggered-event tracks (+ ratio) -----
# NOTE: draw directly on c_pt (the full canvas), spanning the right half.
# Do NOT use c_pt.cd(2) here.
c_pt.cd()

# main spectra pad (top), right half of canvas
p_top = ROOT.TPad("p_top_R", "", 0.5, 0.32, 1.0, 1.0)
p_top.SetLogy()
p_top.SetLeftMargin(0.13)
p_top.SetBottomMargin(0.02)
p_top.Draw()
p_top.cd()

h_injet.SetLineColor(ROOT.kRed + 1)
h_injet.SetLineWidth(2)
f.h_pt_jettrig.SetLineColor(ROOT.kBlue + 1)
f.h_pt_jettrig.SetLineWidth(2)
h_injet_edge.SetLineColor(ROOT.kGreen + 2)
h_injet_edge.SetLineWidth(2)

# optional histogram: all tracks in events with a jet > 8 GeV/c with no #eta restriction
if hasattr(f, 'h_pt_jettrig_noetarestr'):
    f.h_pt_jettrig_noetarestr.SetLineColor(ROOT.kMagenta + 2)
    f.h_pt_jettrig_noetarestr.SetLineWidth(2)

h_injet_sum = h_injet.Clone("h_injet_sum")
h_injet_sum.SetDirectory(0)
h_injet_sum.SetLineColor(ROOT.kGray + 1)
h_injet_sum.SetLineWidth(2)
h_injet_sum.SetLineStyle(7)
h_injet_sum.Add(h_injet_edge)

# # ADDING CHECKS HISTOGRAMS HERE
# f.hCHECK_trackpt_jetsabove8gev_fid.SetLineColorAlpha(ROOT.kRed - 7, 0.8)
# f.hCHECK_trackpt_jetsabove8gev_fid.SetLineWidth(1)
# f.hCHECK_trackpt_jetsabove8gev_fid.SetLineStyle(8)
# f.hCHECK_trackpt_jetsabove8gev_edge.SetLineColorAlpha(ROOT.kGreen - 3, 0.8)
# f.hCHECK_trackpt_jetsabove8gev_edge.SetLineWidth(1)
# f.hCHECK_trackpt_jetsabove8gev_edge.SetLineStyle(8)
# f.hCHECK_trackpt_jetsabove8gev_all.SetLineColorAlpha(ROOT.kYellow, 0.8)
# f.hCHECK_trackpt_jetsabove8gev_all.SetLineWidth(1)
# f.hCHECK_trackpt_jetsabove8gev_all.SetLineStyle(8)
# f.hCHECK_trackpt_evt_withjetsabove8gev_all.SetLineColorAlpha(ROOT.kMagenta - 4, 0.8)
# f.hCHECK_trackpt_evt_withjetsabove8gev_all.SetLineWidth(1)
# f.hCHECK_trackpt_evt_withjetsabove8gev_all.SetLineStyle(8)

# pick the tallest for drawing first
max_list = [h_injet.GetMaximum(), f.h_pt_jettrig.GetMaximum(), h_injet_edge.GetMaximum(),
            h_injet_sum.GetMaximum(), f.h_pt_R06_trig.GetMaximum()]
if hasattr(f, 'h_pt_jettrig_noetarestr'):
    max_list.append(f.h_pt_jettrig_noetarestr.GetMaximum())
hmax = max(max_list)
if hmax == f.h_pt_jettrig.GetMaximum():
    lead = f.h_pt_jettrig
elif hmax == h_injet.GetMaximum():
    lead = h_injet
elif hmax == h_injet_edge.GetMaximum():
    lead = h_injet_edge
elif hmax == f.h_pt_R06_trig.GetMaximum():
    lead = f.h_pt_R06_trig
elif hasattr(f, 'h_pt_jettrig_noetarestr') and hmax == f.h_pt_jettrig_noetarestr.GetMaximum():
    lead = f.h_pt_jettrig_noetarestr
else:
    lead = h_injet_sum

lead.SetStats(0)
lead.SetTitle(";;counts")
lead.GetYaxis().SetTitleSize(0.05)
lead.GetYaxis().SetTitleOffset(1.1)
lead.GetXaxis().SetLabelSize(0)
lead.Draw("HIST")
draw_list = [f.h_pt_jettrig, h_injet, h_injet_edge, h_injet_sum, f.h_pt_R06_trig]
# draw_list.append(f.hCHECK_trackpt_jetsabove8gev_fid)
# draw_list.append(f.hCHECK_trackpt_jetsabove8gev_edge)
# draw_list.append(f.hCHECK_trackpt_jetsabove8gev_all)
# draw_list.append(f.hCHECK_trackpt_evt_withjetsabove8gev_all)
if hasattr(f, 'h_pt_jettrig_noetarestr'):
    draw_list.insert(1, f.h_pt_jettrig_noetarestr)
for h in draw_list:
    if h is not lead:
        h.Draw("HIST SAME")

legR = ROOT.TLegend(0.30, 0.65, 0.88, 0.88)
legR.SetBorderSize(0)
legR.SetFillStyle(0)
legR.AddEntry(f.h_pt_jettrig,  "All tracks in events with a jet > 8 GeV/c", "l")
if hasattr(f, 'h_pt_jettrig_noetarestr'):
    legR.AddEntry(f.h_pt_jettrig_noetarestr,
                 "All tracks in events with a jet > 8 GeV/c with no #eta restriction",
                 "l")
legR.AddEntry(f.h_pt_R06_trig, "All tracks in events with an R=0.6 jet > 8 GeV/c", "l")
legR.AddEntry(h_injet,         "Tracks in all jets > 8 GeV/c, |#eta| #leq 0.5", "l")
legR.AddEntry(h_injet_edge,    "Tracks in all jets > 8 GeV/c, |#eta| > 0.5", "l")
legR.AddEntry(h_injet_sum,     "Sum of in-jet and edge tracks", "l")
# if hasattr(f, 'hCHECK_trackpt_jetsabove8gev_fid'):
#     legR.AddEntry(f.hCHECK_trackpt_jetsabove8gev_fid, "CHECK Tracks in all jets > 8 GeV/c, |#eta| #leq 0.5", "l")
# if hasattr(f, 'hCHECK_trackpt_jetsabove8gev_edge'):
#     legR.AddEntry(f.hCHECK_trackpt_jetsabove8gev_edge, "CHECK Tracks in all jets > 8 GeV/c, |#eta| > 0.5", "l")
# if hasattr(f, 'hCHECK_trackpt_jetsabove8gev_all'):
#     legR.AddEntry(f.hCHECK_trackpt_jetsabove8gev_all, "CHECK Tracks in all jets > 8 GeV/c, all eta", "l")
# if hasattr(f, 'hCHECK_trackpt_evt_withjetsabove8gev_all'):
#     legR.AddEntry(f.hCHECK_trackpt_evt_withjetsabove8gev_all, "CHECK Tracks in events with a jet > 8 GeV/c, all eta", "l")
legR.Draw()

# ratio pad (bottom), right half of canvas
c_pt.cd()
p_bot = ROOT.TPad("p_bot_R", "", 0.5, 0.0, 1.0, 0.32)
p_bot.SetLeftMargin(0.13)
p_bot.SetTopMargin(0.02)
p_bot.SetBottomMargin(0.3)
p_bot.Draw()
p_bot.cd()

h_ratio = h_injet.Clone("h_ratio_R")
h_ratio.SetStats(0)
h_ratio.SetLineColor(ROOT.kRed + 1)
h_ratio.SetLineWidth(2)
h_ratio.Divide(f.h_pt_jettrig_noetarestr)            # in-jet (fiducial) / all

h_ratio_edge = h_injet_edge.Clone("h_ratio_edge_R")
h_ratio_edge.SetStats(0)
h_ratio_edge.SetLineColor(ROOT.kGreen + 2)
h_ratio_edge.SetLineWidth(2)
h_ratio_edge.Divide(f.h_pt_jettrig_noetarestr)       # in-jet (edge) / all

h_ratio.SetTitle(";p_{T} (GeV/c);in-jet / all (purple)")
h_ratio.GetYaxis().SetNdivisions(505)
h_ratio.GetYaxis().SetTitleSize(0.11)
h_ratio.GetYaxis().SetTitleOffset(0.45)
h_ratio.GetYaxis().SetLabelSize(0.09)
h_ratio.GetXaxis().SetTitleSize(0.11)
h_ratio.GetXaxis().SetTitleOffset(1.0)
h_ratio.GetXaxis().SetLabelSize(0.09)
# h_ratio.GetYaxis().SetRangeUser(0, 1.1)
h_ratio.Draw("HIST")
h_ratio_edge.Draw("HIST SAME")

line = ROOT.TLine(h_ratio.GetXaxis().GetXmin(), 1.0, h_ratio.GetXaxis().GetXmax(), 1.0)
line.SetLineStyle(2)
line.SetLineColor(ROOT.kGray + 2)
line.Draw()

c_pt.Update()
c_pt.SaveAs(f"{output_base}/pt_panels_1.pdf")


# ============================================================
# Jet pT comparison: R=0.4 vs R=0.6 vs R=0.4-in-R=0.6-events
# with ratio panel (R=0.6 / R=0.4)
# ============================================================
def make_jetpt_panels(h_R04_fine, h_R06, h_R04_in_R06, outname,
                      r06_label="R = 0.6 jets",
                      r04_in_label="R = 0.4 jets in R = 0.6 events",
                      h_extra=None, extra_label=None, tag=""):
    """Two-panel jet pT figure: spectra (top) + R=0.6/R=0.4 ratio (bottom).
    `h_extra` is an optional additional spectrum drawn on the top pad, and
    (if present) a second R=0.6/extra ratio is drawn on the bottom pad,
    colored to match the extra curve.
    `tag` keeps ROOT object names unique between calls."""
    c = ROOT.TCanvas(f"c_jetpt{tag}", "", 700, 700)

    # rebin R=0.4 (1000 bins) to match R=0.6 hists (100 bins); both span 0-100
    h_R04 = h_R04_fine.Clone(f"h_jet_pt_rb{tag}")
    h_R04.SetDirectory(0)
    h_R04.Rebin(10)

    # top pad: spectra
    p_top = ROOT.TPad(f"p_jet_top{tag}", "", 0.0, 0.32, 1.0, 1.0)
    p_top.SetLogy()
    p_top.SetLeftMargin(0.13)
    p_top.SetBottomMargin(0.02)
    p_top.Draw()
    p_top.cd()

    EXTRA_COLOR = ROOT.kMagenta + 1
    h_R04.SetLineColor(ROOT.kBlue + 1);  h_R04.SetLineWidth(2)
    h_R06.SetLineColor(ROOT.kRed + 1);   h_R06.SetLineWidth(2)
    h_R04_in_R06.SetLineColor(ROOT.kGreen + 2); h_R04_in_R06.SetLineWidth(2)

    # build the list of spectra to draw, appending the optional extra
    hists = [h_R04, h_R06, h_R04_in_R06]
    if h_extra is not None:
        h_extra.SetLineColor(EXTRA_COLOR)
        h_extra.SetLineWidth(2)
        hists.append(h_extra)

    # draw tallest first
    lead = max(hists, key=lambda h: h.GetMaximum())
    lead.SetStats(0)
    lead.SetTitle(";;counts")
    lead.GetYaxis().SetTitleSize(0.05)
    lead.GetYaxis().SetTitleOffset(1.1)
    lead.GetXaxis().SetLabelSize(0)
    lead.Draw("HIST")
    for h in hists:
        if h is not lead:
            h.Draw("HIST SAME")

    leg = ROOT.TLegend(0.40, 0.68, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(h_R06,        r06_label, "l")
    leg.AddEntry(h_R04,        "R = 0.4 jets, |#eta| < 0.5", "l")
    leg.AddEntry(h_R04_in_R06, r04_in_label, "l")
    if h_extra is not None:
        leg.AddEntry(h_extra, extra_label, "l")
    leg.Draw()

    # bottom pad: ratio R=0.6 / R=0.4 (+ R=0.6 / extra if present)
    c.cd()
    p_bot = ROOT.TPad(f"p_jet_bot{tag}", "", 0.0, 0.0, 1.0, 0.32)
    p_bot.SetLeftMargin(0.13)
    p_bot.SetTopMargin(0.02)
    p_bot.SetBottomMargin(0.3)
    p_bot.SetGridy()
    p_bot.Draw()
    p_bot.cd()

    h_ratio = h_R04.Clone(f"h_jet_ratio{tag}")
    h_ratio.SetDirectory(0)
    h_ratio.SetStats(0)
    h_ratio.SetLineColor(ROOT.kBlue)
    h_ratio.SetLineWidth(2)
    h_ratio.Divide(h_R06)

    # second ratio vs. the extra curve, colored to match it
    h_ratio_extra = None
    if h_extra is not None:
        h_ratio_extra = h_extra.Clone(f"h_jet_ratio_extra{tag}")
        h_ratio_extra.SetDirectory(0)
        h_ratio_extra.SetStats(0)
        h_ratio_extra.SetLineColor(EXTRA_COLOR)
        h_ratio_extra.SetLineWidth(2)
        h_ratio_extra.Divide(h_R06)

    ylabel = "R=0.4 / R=0.6" if h_extra is None else "R=0.4 / R=0.6 ratios"
    h_ratio.SetTitle(f";p_{{T,jet}} (GeV/c);{ylabel}")
    h_ratio.GetYaxis().SetNdivisions(505)
    h_ratio.GetYaxis().SetTitleSize(0.11)
    h_ratio.GetYaxis().SetTitleOffset(0.45)
    h_ratio.GetYaxis().SetLabelSize(0.09)
    h_ratio.GetXaxis().SetTitleSize(0.11)
    h_ratio.GetXaxis().SetTitleOffset(1.0)
    h_ratio.GetXaxis().SetLabelSize(0.09)

    # set a common y-range so both ratios are visible
    if h_ratio_extra is not None:
        ymin = min(h_ratio.GetMinimum(0.0), h_ratio_extra.GetMinimum(0.0))
        ymax = max(h_ratio.GetMaximum(),    h_ratio_extra.GetMaximum())
        pad = 0.05 * (ymax - ymin) if ymax > ymin else 0.1
        h_ratio.GetYaxis().SetRangeUser(ymin - pad, ymax + pad)

    h_ratio.Draw("HIST")
    if h_ratio_extra is not None:
        h_ratio_extra.Draw("HIST SAME")

    line = ROOT.TLine(h_ratio.GetXaxis().GetXmin(), 1.0,
                      h_ratio.GetXaxis().GetXmax(), 1.0)
    line.SetLineStyle(2)
    line.SetLineColor(ROOT.kGray + 2)
    line.Draw()

    c.Update()
    c.SaveAs(outname)

    # keep references alive so ROOT doesn't garbage-collect the pads/objects
    return c, p_top, p_bot, leg, h_ratio, h_ratio_extra, line, h_R04, h_extra


# non-fiducial
_keep1 = make_jetpt_panels(
    f.h_jet_pt, f.h_jet_pt_R06, f.h_jet_pt_R04_with_R06,
    f"{output_base}/jet_pt_panels_1.pdf",
)

# fiducial (|eta| < 0.3)
_keep2 = make_jetpt_panels(
    f.h_jet_pt, f.h_jet_pt_R06_fid, f.h_jet_pt_R04_with_R04_fid,
    f"{output_base}/jet_pt_panels_1_fid.pdf",
    r06_label="R = 0.6 jets, |#eta| < 0.3",
    r04_in_label="R = 0.4 jets, |#eta| < 0.5 in R = 0.6 events (|#eta| < 0.3)",
    h_extra=f.h_jet_pt_R04_with_R06_fid,
    extra_label="R = 0.4, |#eta| < 0.3 jets in R = 0.6 events (|#eta| < 0.3)",
    tag="_fid",
)


# === ETA ===
# ============================================================
# Two-panel canvas for eta (mirrors the pT version)
# ============================================================
h_injet_eta = f.h_tracketa_jetpt.ProjectionY("h_injet_eta", xlo, xhi)
h_injet_eta.SetDirectory(0)

h_injet_eta_edge = f.h_tracketa_jetpt_edge.ProjectionY("h_injet_eta_edge", bin_lo, bin_hi)

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
h_injet_eta_edge.SetLineColor(ROOT.kGreen + 2)
h_injet_eta_edge.SetLineWidth(2)

hmax = max(h_injet_eta.GetMaximum(), f.h_eta_jettrig.GetMaximum(), h_injet_eta_edge.GetMaximum())
if hmax == f.h_eta_jettrig.GetMaximum():
    lead = f.h_eta_jettrig
elif hmax == h_injet_eta.GetMaximum():
    lead = h_injet_eta
else:
    lead = h_injet_eta_edge

lead.SetStats(0)
lead.SetTitle(";;counts")
lead.GetYaxis().SetTitleSize(0.05)
lead.GetYaxis().SetTitleOffset(1.1)
lead.GetXaxis().SetLabelSize(0)
lead.Draw("HIST")
for h in (f.h_eta_jettrig, h_injet_eta, h_injet_eta_edge):
    if h is not lead:
        h.Draw("HIST SAME")

legR_eta = ROOT.TLegend(0.25, 0.18, 0.80, 0.32)
legR_eta.SetBorderSize(0); legR_eta.SetFillStyle(0)
legR_eta.AddEntry(f.h_eta_jettrig,  "All tracks in events with a jet > 8 GeV/c", "l")
legR_eta.AddEntry(h_injet_eta,      "Tracks in all jets > 8 GeV/c, |#eta| #leq 0.5", "l")
legR_eta.AddEntry(h_injet_eta_edge,    "Tracks in all jets > 8 GeV/c, |#eta| > 0.5", "l")
legR_eta.Draw()

c_eta.Update()
c_eta.SaveAs(f"{output_base}/eta_panels_1.pdf")


# ============================================================
# track eta vs track pt: fiducial / triggered-fiducial / edge (>8 GeV) constituents
# ============================================================

# ---- the three 2D histograms ----
h_fid     = f.h_tracketa_trackpt_fid
h_trigfid = f.h_tracketa_trackpt_jettrig_fid
h_edge    = f.h_tracketa_trackpt_edge

# for log z, the minimum must be > 0; use the smallest nonzero bin content
def min_nonzero(h):
    m = None
    for bx in range(1, h.GetNbinsX() + 1):
        for by in range(1, h.GetNbinsY() + 1):
            c = h.GetBinContent(bx, by)
            if c > 0 and (m is None or c < m):
                m = c
    return m if m is not None else 1.0

# ---- common z-range across all three histograms ----
zmin = min(min_nonzero(h_fid), min_nonzero(h_trigfid), min_nonzero(h_edge))
zmax = max(h_fid.GetMaximum(), h_trigfid.GetMaximum(), h_edge.GetMaximum())

for h in (h_fid, h_trigfid, h_edge):
    h.SetMinimum(zmin)
    h.SetMaximum(zmax)
    h.SetStats(0)

c_2d = ROOT.TCanvas("c_2d", "", 1800, 600)

# ---- left: 2D fiducial (no color bar) ----
p_left = ROOT.TPad("p_left", "", 0.0, 0.0, 0.34, 1.0)
p_left.SetRightMargin(0.02); p_left.SetLogz()
p_left.Draw()
p_left.cd()
h_fid.SetTitle("Constituents of fiducial jets, p_{T,jet} > 5 GeV/c;p_{T}^{track} [GeV/c];#eta^{track}")
h_fid.Draw("COL")
c_2d.cd()

# ---- middle: 2D triggered fiducial (no color bar) ----
p_mid = ROOT.TPad("p_mid", "", 0.34, 0.0, 0.66, 1.0)
p_mid.SetRightMargin(0.02); p_mid.SetLogz()
p_mid.Draw()
p_mid.cd()
h_trigfid.SetTitle("Constituents of fiducial jets, p_{T,jet} > 8 GeV/c;p_{T}^{track} [GeV/c];#eta^{track}")
h_trigfid.Draw("COL")
c_2d.cd()

# ---- right: 2D edge (with shared color bar) ----
p_right = ROOT.TPad("p_right", "", 0.66, 0.0, 1.0, 1.0)
p_right.SetRightMargin(0.15); p_right.SetLogz()
p_right.Draw()
p_right.cd()
h_edge.SetTitle("Constituents of edge jets, p_{T,jet} > 8 GeV/c;p_{T}^{track} [GeV/c];#eta^{track}")
h_edge.Draw("COLZ")
c_2d.cd()

c_2d.SaveAs(f"{output_base}/tracketa_trackpt_panels_1.pdf")

# c_2d.cd(1)
# ROOT.gPad.SetRightMargin(0.15); ROOT.gPad.SetLogz()
# f.h_tracketa_trackpt_fid.SetTitle("Constituents of fiducial jets;p_{T}^{track} [GeV/c];#eta^{track}")
# f.h_tracketa_trackpt_fid.Draw("COLZ")

# c_2d.cd(2)
# ROOT.gPad.SetRightMargin(0.15); ROOT.gPad.SetLogz()
# f.h_tracketa_trackpt_edge.SetTitle("Constituents of jets>8 GeV outside fiducial;p_{T}^{track} [GeV/c];#eta^{track}")
# f.h_tracketa_trackpt_edge.Draw("COLZ")



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
    c_ext.SaveAs(f"{output_base}/h_pt_ext_1.pdf")

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

    c_comp.SaveAs(f"{output_base}/h_pt_comp_1.pdf")
    
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

    c_comp_2.SaveAs(f"{output_base}/h_pt_comp_2.pdf")