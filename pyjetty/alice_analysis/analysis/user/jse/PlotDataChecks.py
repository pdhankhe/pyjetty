#!/usr/bin/env python3
"""Data-check plots (track/jet kinematics + inclusive EEC) for the JSE analysis."""

import math
import ROOT

# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------
SITE = "perlmutter"          # "perlmutter" or "hiccup"
JOBID = "53567700"

PATHS = {
    "perlmutter": dict(
        input_base="/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data_checks/56458992",
        extra_input_base=f"/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data_checks/{JOBID}",
        output_base="/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots/data_checks",
        full_eec_file="/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data/57696990/AnalysisResultsMerged_ungroomedbins.root" #"/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data/55778272/AnalysisResultsMerged.root",
    ),
    "hiccup": dict(
        input_base="/rstorage/alice/AnalysisResults/blianggi/jse/data_checks/1839116",
        extra_input_base="",
        output_base="/software/users/blianggi/mypyjetty/storage/jse/plots/data_checks",
        full_eec_file="/rstorage/alice/AnalysisResults/blianggi/jse/data/55778272/AnalysisResultsMerged.root",
    ),
}
CFG = PATHS[SITE]
OUT = CFG["output_base"]

JETPT_LO, JETPT_HI = 8.0, 20.0        # jet-pT window for the y-projections
ALICE_LABEL = ("ALICE pp #sqrt{s} = 5.36 TeV",
               "LHC24 ppref pass 1, JE derived",
               "anti-k_{T} R = 0.4, |#eta_{jet}| < 0.5")

# EEC section
GROOM_TAG = "sd0.1"                   # "sd0.1" or "maxkt" (hist_full is ungroomed anyway)
USE_PTRL = False                      # False -> hist_full_...  True -> hist_full_ptRL_...
PT_BINS = [(10, 20), (20, 40), (40, 60), (60, 80),
           (80, 100), (100, 120), (120, 150), (150, 200)]
COUNTER_NAME = "counters_jetpt%d_%d_%s"
BIN_NJETS, BIN_SUMPT = 2, 3
COLORS = [ROOT.kBlack, ROOT.kBlue, ROOT.kOrange + 7, ROOT.kGreen + 2,
          ROOT.kRed + 1, ROOT.kMagenta + 1, ROOT.kCyan + 2, ROOT.kAzure + 1,
          ROOT.kYellow + 1, ROOT.kGray + 1]
MARKERS = [ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kFullDiamond,
           ROOT.kFullStar, ROOT.kFullTriangleUp, ROOT.kFullTriangleDown,
           ROOT.kFullCross, ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kFullDiamond]

_keep = []                            # references so ROOT doesn't collect pads/legends


# ---------------------------------------------------------------------------
# generic ROOT helpers
# ---------------------------------------------------------------------------
def style(h, color=None, width=2, linestyle=None, marker=None, msize=None,
          title=None, stats=False):
    """Apply the usual line/marker cosmetics and return the histogram."""
    if color is not None:
        h.SetLineColor(color)
        if marker is not None:
            h.SetMarkerColor(color)
    if width is not None:
        h.SetLineWidth(width)
    if linestyle is not None:
        h.SetLineStyle(linestyle)
    if marker is not None:
        h.SetMarkerStyle(marker)
    if msize is not None:
        h.SetMarkerSize(msize)
    if title is not None:
        h.SetTitle(title)
    h.SetStats(stats)
    return h


def make_legend(box, entries, header=None, textsize=None):
    """entries: iterable of (hist, label) or (hist, label, draw_option)."""
    leg = ROOT.TLegend(*box)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    if textsize is not None:
        leg.SetTextSize(textsize)
    if header is not None:
        leg.SetHeader(header)
    for entry in entries:
        leg.AddEntry(entry[0], entry[1], entry[2] if len(entry) > 2 else "l")
    leg.Draw()
    _keep.append(leg)
    return leg


def draw_stack(hists, ytitle="counts", hide_xlabels=True):
    """Draw the tallest histogram first (it owns the axes), the rest with SAME."""
    lead = max(hists, key=lambda h: h.GetMaximum())
    lead.SetStats(0)
    lead.SetTitle(";;%s" % ytitle)
    lead.GetYaxis().SetTitleSize(0.05)
    lead.GetYaxis().SetTitleOffset(1.1)
    if hide_xlabels:
        lead.GetXaxis().SetLabelSize(0)
    lead.Draw("HIST")
    for h in hists:
        if h is not lead:
            h.Draw("HIST SAME")
    return lead


def split_pads(canvas, tag, x1=0.0, x2=1.0, split=0.32, logy=False, gridy=False):
    """Top (spectra) / bottom (ratio) pad pair spanning [x1, x2] of `canvas`."""
    canvas.cd()
    top = ROOT.TPad("p_top%s" % tag, "", x1, split, x2, 1.0)
    top.SetLeftMargin(0.13)
    top.SetBottomMargin(0.02)
    top.SetLogy(logy)
    top.Draw()

    canvas.cd()
    bot = ROOT.TPad("p_bot%s" % tag, "", x1, 0.0, x2, split)
    bot.SetLeftMargin(0.13)
    bot.SetTopMargin(0.02)
    bot.SetBottomMargin(0.3)
    if gridy:
        bot.SetGridy()
    bot.Draw()

    _keep.extend([top, bot])
    return top, bot


def make_ratio(num, den, name, color):
    r = num.Clone(name)
    r.SetDirectory(0)
    style(r, color)
    r.Divide(den)
    _keep.append(r)
    return r


def style_ratio_axes(h, xtitle, ytitle):
    h.SetTitle(";%s;%s" % (xtitle, ytitle))
    ya, xa = h.GetYaxis(), h.GetXaxis()
    ya.SetNdivisions(505)
    ya.SetTitleSize(0.11); ya.SetTitleOffset(0.45); ya.SetLabelSize(0.09)
    xa.SetTitleSize(0.11); xa.SetTitleOffset(1.0); xa.SetLabelSize(0.09)


def draw_unity_line(h):
    line = ROOT.TLine(h.GetXaxis().GetXmin(), 1.0, h.GetXaxis().GetXmax(), 1.0)
    line.SetLineStyle(2)
    line.SetLineColor(ROOT.kGray + 2)
    line.Draw()
    _keep.append(line)
    return line


def draw_alice_label(x=0.17, ytop=0.84, lines=ALICE_LABEL, dy=0.04):
    """Draw the standard label block on the *current* pad."""
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(0.035)
    lat.SetTextAlign(11)
    for i, text in enumerate(lines):
        lat.DrawLatex(x, ytop - i * dy, text)
    _keep.append(lat)
    return lat


def framed_canvas(name, xrange_, yrange, xtitle, ytitle,
                  size=(900, 700), logx=False, logy=False, grid=True):
    c = ROOT.TCanvas(name, "", *size)
    c.SetLeftMargin(0.13)
    c.SetBottomMargin(0.13)
    c.SetRightMargin(0.05)
    c.SetTicks(1, 1)
    c.SetGrid(grid, grid)
    c.SetLogx(logx)
    c.SetLogy(logy)
    frame = c.DrawFrame(xrange_[0], yrange[0], xrange_[1], yrange[1])
    frame.GetXaxis().SetTitle(xtitle)
    frame.GetYaxis().SetTitle(ytitle)
    frame.GetXaxis().SetTitleSize(0.045)
    frame.GetYaxis().SetTitleSize(0.045)
    frame.GetYaxis().SetTitleOffset(1.35)
    _keep.extend([c, frame])
    return c, frame


def project_jetpt_window(h2d, name, lo=JETPT_LO, hi=JETPT_HI):
    """Y-projection of a (jet pT, track X) 2D over the [lo, hi] jet-pT window."""
    ax, eps = h2d.GetXaxis(), 1e-6
    h = h2d.ProjectionY(name, ax.FindBin(lo + eps), ax.FindBin(hi + eps))
    h.SetDirectory(0)                 # survive the file going out of scope
    return h


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def plot_track_qa(f, num_acc_events):
    """merged_1.pdf (raw) and merged_1_divbyevents.pdf (per-event normalized)."""
    f.h_phi.SetTitle("global track #varphi")
    f.h_phi.GetXaxis().SetTitle("#varphi")

    c = ROOT.TCanvas("c_qa", "", 1500, 450)
    c.Divide(3, 1)
    for i, (h, logy) in enumerate([(f.h_pt, True), (f.h_eta, False), (f.h_phi, False)], 1):
        c.cd(i)
        ROOT.gPad.SetLogy(logy)
        h.Draw()
    c.SaveAs("%s/merged_1.pdf" % OUT)

    panels = [
        (f.h_pt,  "#frac{1}{N_{events}}#frac{dN_{tracks}}{dp_{T}} (GeV/c)^{-1}", True,  None),
        (f.h_eta, "#frac{1}{N_{events}}#frac{dN_{tracks}}{d#eta}",               False, None),
        (f.h_phi, "#frac{1}{N_{events}}#frac{dN_{tracks}}{d#varphi}",            False, 0.0),
    ]
    c_alt = ROOT.TCanvas("c_alt", "", 1500, 450)
    c_alt.Divide(3, 1)
    for i, (h, ytitle, logy, ymin) in enumerate(panels, 1):
        c_alt.cd(i)
        ROOT.gPad.SetLeftMargin(0.16)
        ROOT.gPad.SetLogy(logy)
        h.SetStats(0)
        h.GetYaxis().SetTitleOffset(1.8)
        h.Scale(1.0 / num_acc_events)
        h.GetYaxis().SetTitle(ytitle)
        if ymin is not None:
            h.SetMinimum(ymin)
        h.Draw()
    c_alt.SaveAs("%s/merged_1_divbyevents.pdf" % OUT)
    _keep.extend([c, c_alt])


def plot_pt_panels(f):
    """pt_panels_1.pdf: all-event track pT (left) vs in-jet / triggered pT (right)."""
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadGridX(1)
    ROOT.gStyle.SetPadGridY(1)

    print("Getnbinsx():", f.h_trackpt_jetpt.GetXaxis().GetNbins())
    h_injet = project_jetpt_window(f.h_trackpt_jetpt, "h_injet")
    h_injet_edge = project_jetpt_window(f.h_trackpt_jetpt_edge, "h_injet_edge")

    h_sum = h_injet.Clone("h_injet_sum")
    h_sum.SetDirectory(0)
    h_sum.Add(h_injet_edge)

    style(h_injet, ROOT.kRed + 1)
    style(h_injet_edge, ROOT.kGreen + 2)
    style(h_sum, ROOT.kGray + 1, linestyle=7)
    style(f.h_pt_jettrig, ROOT.kBlue + 1)
    style(f.h_pt_R06_trig, ROOT.kOrange + 7)
    h_noeta = getattr(f, "h_pt_jettrig_noetarestr", None)
    if h_noeta:
        style(h_noeta, ROOT.kMagenta + 2)

    c_pt = ROOT.TCanvas("c_pt", "", 1200, 600)
    c_pt.Divide(2, 1)

    # ----- left: all tracks in all events -----
    c_pt.cd(1)
    ROOT.gPad.SetLogy()
    ROOT.gPad.SetLeftMargin(0.13)
    style(f.h_pt, ROOT.kBlack, width=1, title=";p_{T} (GeV/c);counts").Draw("HIST")
    style(f.h_pt_R06_trig_fid, ROOT.kOrange + 2, width=1).Draw("HIST SAME")
    make_legend((0.40, 0.82, 0.88, 0.88), [
        (f.h_pt, "All tracks in all events"),
        (f.h_pt_R06_trig_fid, "Tracks in events with R=0.6 jet > 8 GeV/c, jet |eta| < 0.3"),
    ])

    # ----- right half of the canvas: spectra + ratio (drawn on c_pt itself) -----
    p_top, p_bot = split_pads(c_pt, "_R", x1=0.5, logy=True)

    p_top.cd()
    hists = [f.h_pt_jettrig] + ([h_noeta] if h_noeta else []) + \
            [h_injet, h_injet_edge, h_sum, f.h_pt_R06_trig]
    draw_stack(hists)
    entries = [(f.h_pt_jettrig, "All tracks in events with a jet > 8 GeV/c")]
    if h_noeta:
        entries.append((h_noeta, "All tracks in events with a jet > 8 GeV/c with no #eta restriction"))
    entries += [
        (f.h_pt_R06_trig, "All tracks in events with an R=0.6 jet > 8 GeV/c"),
        (h_injet,         "Tracks in all jets > 8 GeV/c, |#eta| #leq 0.5"),
        (h_injet_edge,    "Tracks in all jets > 8 GeV/c, |#eta| > 0.5"),
        (h_sum,           "Sum of in-jet and edge tracks"),
    ]
    make_legend((0.30, 0.65, 0.88, 0.88), entries)

    p_bot.cd()
    denom = f.h_pt_jettrig_noetarestr
    h_ratio = make_ratio(h_injet, denom, "h_ratio_R", ROOT.kRed + 1)
    h_ratio_edge = make_ratio(h_injet_edge, denom, "h_ratio_edge_R", ROOT.kGreen + 2)
    style_ratio_axes(h_ratio, "p_{T} (GeV/c)", "in-jet / all (purple)")
    h_ratio.Draw("HIST")
    h_ratio_edge.Draw("HIST SAME")
    draw_unity_line(h_ratio)

    c_pt.Update()
    c_pt.SaveAs("%s/pt_panels_1.pdf" % OUT)
    _keep.extend([c_pt, h_injet, h_injet_edge, h_sum])
    return h_ratio


def make_jetpt_panels(h_R04_fine, h_R06, h_R04_in_R06, outname,
                      r06_label="R = 0.6 jets",
                      r04_in_label="R = 0.4 jets in R = 0.6 events",
                      h_extra=None, extra_label=None, tag=""):
    """Jet pT spectra (top) + R=0.4/R=0.6 ratio(s) (bottom).

    `h_extra` is an optional extra spectrum; if given, a second ratio vs. R=0.6
    is drawn in the same colour. `tag` keeps ROOT object names unique.
    """
    EXTRA_COLOR = ROOT.kMagenta + 1
    c = ROOT.TCanvas("c_jetpt%s" % tag, "", 700, 700)

    # rebin R=0.4 (1000 bins) to match the R=0.6 hists (100 bins); both span 0-100
    h_R04 = h_R04_fine.Clone("h_jet_pt_rb%s" % tag)
    h_R04.SetDirectory(0)
    h_R04.Rebin(10)

    p_top, p_bot = split_pads(c, "_jet%s" % tag, logy=True, gridy=True)

    p_top.cd()
    hists = [style(h_R04, ROOT.kBlue + 1),
             style(h_R06, ROOT.kRed + 1),
             style(h_R04_in_R06, ROOT.kGreen + 2)]
    if h_extra is not None:
        hists.append(style(h_extra, EXTRA_COLOR))
    draw_stack(hists)

    entries = [(h_R06, r06_label),
               (h_R04, "R = 0.4 jets, |#eta| < 0.5"),
               (h_R04_in_R06, r04_in_label)]
    if h_extra is not None:
        entries.append((h_extra, extra_label))
    make_legend((0.40, 0.68, 0.88, 0.88), entries)
    draw_alice_label()

    p_bot.cd()
    h_ratio = make_ratio(h_R04, h_R06, "h_jet_ratio%s" % tag, ROOT.kBlue)
    h_ratio_extra = (make_ratio(h_extra, h_R06, "h_jet_ratio_extra%s" % tag, EXTRA_COLOR)
                     if h_extra is not None else None)

    ylabel = "R=0.4 / R=0.6" + ("" if h_extra is None else " ratios")
    style_ratio_axes(h_ratio, "p_{T,jet} (GeV/c)", ylabel)
    if h_ratio_extra is not None:                       # common y-range for both
        ymin = min(h_ratio.GetMinimum(0.0), h_ratio_extra.GetMinimum(0.0))
        ymax = max(h_ratio.GetMaximum(), h_ratio_extra.GetMaximum())
        pad = 0.05 * (ymax - ymin) if ymax > ymin else 0.1
        h_ratio.GetYaxis().SetRangeUser(ymin - pad, ymax + pad)
    h_ratio.Draw("HIST")
    if h_ratio_extra is not None:
        h_ratio_extra.Draw("HIST SAME")
    draw_unity_line(h_ratio)

    c.Update()
    c.SaveAs(outname)
    _keep.extend([c, h_R04])
    return c


def plot_jetpt_figures(f):
    make_jetpt_panels(f.h_jet_pt, f.h_jet_pt_R06, f.h_jet_pt_R04_with_R06,
                      "%s/jet_pt_panels_1.pdf" % OUT)
    make_jetpt_panels(
        f.h_jet_pt, f.h_jet_pt_R06_fid, f.h_jet_pt_R04_with_R04_fid,
        "%s/jet_pt_panels_1_fid.pdf" % OUT,
        r06_label="R = 0.6 jets, |#eta| < 0.3",
        r04_in_label="R = 0.4 jets, |#eta| < 0.5 in R = 0.6 events (|#eta| < 0.3)",
        h_extra=f.h_jet_pt_R04_with_R06_fid,
        extra_label="R = 0.4, |#eta| < 0.3 jets in R = 0.6 events (|#eta| < 0.3)",
        tag="_fid",
    )


def plot_tracks_in_R_jets(f, R=0.6):
    Rval = "06" if R == 0.6 else "04"
    if R == 0.6:
        h = getattr(f, "h_trackpt_R06_trig", None)
        h_fid = getattr(f, "h_trackpt_R06_trig_fid", None)
    elif R == 0.4:
        # h = getattr(f, "h_trackpt_jetpt", None)
        # h_fid = getattr(f, "h_trackpt_jetpt_edge", None)
        h = project_jetpt_window(f.h_trackpt_jetpt, "h_injet")
        h_fid = project_jetpt_window(f.h_trackpt_jetpt_edge, "h_injet_edge")

    if not (h and h_fid):
        print("warning: h_trackpt_R%s_trig or h_trackpt_R%s_trig_fid not found in input file" % (Rval, Rval))
        return

    h = style(h.Clone("h_trackpt_R%s_trig_plot" % Rval), ROOT.kBlue + 1)
    h.SetDirectory(0)
    style(h_fid, ROOT.kGreen + 2)

    ymax = h.GetMaximum()
    c, _ = framed_canvas("c_tracks_in_R%s_jets" % Rval,
                         (h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax()),
                         (1e4, ymax * 2 if ymax > 0 else 1.0),
                         "p_{T, track}", "counts", logy=True)
    print("BOUNDS HERE!!! ", h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax(), "y:", ymax * 2 if ymax > 0 else 1.0)
    h.Draw("HIST SAME")
    h_fid.Draw("HIST SAME")
    make_legend((0.45, 0.62, 0.92, 0.72), [
        (h, "tracks in R=%.1f jets > 8 GeV" % R),
        (h_fid, "tracks in R=%.1f jets > 8 GeV (fiducial)" % R),
    ], textsize=0.032)
    draw_alice_label(x=0.465)
    c.RedrawAxis()
    c.SaveAs("%s/tracks_in_R%s_jets.pdf" % (OUT, Rval))
    _keep.append(h)


def plot_jet_eta(f):
    h = getattr(f, "h_jet_eta", None)
    if not h:
        print("warning: h_jet_eta not found in input file")
        return

    h = style(h.Clone("h_jet_eta_plot"), ROOT.kBlue + 1, title=";#eta_{jet};counts")
    h.SetDirectory(0)
    ymax = h.GetMaximum()
    c, _ = framed_canvas("c_jet_eta",
                         (h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax()),
                         (0.0, ymax * 1.20 if ymax > 0 else 1.0),
                         "#eta_{jet}", "counts")
    h.Draw("HIST SAME")
    make_legend((0.62, 0.72, 0.92, 0.88), [(h, "jet #eta spectrum")], textsize=0.032)
    draw_alice_label()
    c.RedrawAxis()
    c.SaveAs("%s/jet_eta_spectrum.pdf" % OUT)
    _keep.append(h)


def plot_eta_panels(f):
    """eta_panels_1.pdf -- the eta analogue of the pT two-panel figure."""
    h_injet = project_jetpt_window(f.h_tracketa_jetpt, "h_injet_eta")
    h_injet_edge = project_jetpt_window(f.h_tracketa_jetpt_edge, "h_injet_eta_edge")

    c = ROOT.TCanvas("c_eta", "", 1200, 600)
    c.Divide(2, 1)

    c.cd(1)
    ROOT.gPad.SetLeftMargin(0.13)
    style(f.h_eta, ROOT.kBlack, title=";#eta;counts").Draw("HIST")
    make_legend((0.30, 0.18, 0.78, 0.28), [(f.h_eta, "All tracks in all events")])

    c.cd(2)
    ROOT.gPad.SetLeftMargin(0.13)
    style(h_injet, ROOT.kRed + 1, title=";#eta;counts")
    style(h_injet_edge, ROOT.kGreen + 2)
    style(f.h_eta_jettrig, ROOT.kBlue + 1)
    draw_stack([f.h_eta_jettrig, h_injet, h_injet_edge])
    make_legend((0.25, 0.18, 0.80, 0.32), [
        (f.h_eta_jettrig, "All tracks in events with a jet > 8 GeV/c"),
        (h_injet,         "Tracks in all jets > 8 GeV/c, |#eta| #leq 0.5"),
        (h_injet_edge,    "Tracks in all jets > 8 GeV/c, |#eta| > 0.5"),
    ])

    c.Update()
    c.SaveAs("%s/eta_panels_1.pdf" % OUT)
    _keep.extend([c, h_injet, h_injet_edge])


def min_nonzero(h):
    """Smallest strictly-positive bin content of a 2D hist (for log z)."""
    contents = (h.GetBinContent(bx, by)
                for bx in range(1, h.GetNbinsX() + 1)
                for by in range(1, h.GetNbinsY() + 1))
    positive = [c for c in contents if c > 0]
    return min(positive) if positive else 1.0


def plot_tracketa_trackpt(f):
    """tracketa_trackpt_panels_1.pdf -- three 2D maps sharing one z-range."""
    panels = [
        (f.h_tracketa_trackpt_fid,
         "Constituents of fiducial jets, p_{T,jet} > 5 GeV/c", (0.00, 0.34), 0.02, "COL"),
        (f.h_tracketa_trackpt_jettrig_fid,
         "Constituents of fiducial jets, p_{T,jet} > 8 GeV/c", (0.34, 0.66), 0.02, "COL"),
        (f.h_tracketa_trackpt_edge,
         "Constituents of edge jets, p_{T,jet} > 8 GeV/c",     (0.66, 1.00), 0.15, "COLZ"),
    ]
    hists = [p[0] for p in panels]
    zmin = min(min_nonzero(h) for h in hists)
    zmax = max(h.GetMaximum() for h in hists)

    c = ROOT.TCanvas("c_2d", "", 1800, 600)
    for i, (h, title, (x1, x2), rmargin, opt) in enumerate(panels):
        c.cd()
        pad = ROOT.TPad("p_2d_%d" % i, "", x1, 0.0, x2, 1.0)
        pad.SetRightMargin(rmargin)
        pad.SetLogz()
        pad.Draw()
        pad.cd()
        h.SetMinimum(zmin)
        h.SetMaximum(zmax)
        h.SetStats(0)
        h.SetTitle("%s;p_{T}^{track} [GeV/c];#eta^{track}" % title)
        h.Draw(opt)
        _keep.append(pad)
    c.SaveAs("%s/tracketa_trackpt_panels_1.pdf" % OUT)
    _keep.append(c)


# ---------------------------------------------------------------------------
# inclusive EEC across jet pT bins
# ---------------------------------------------------------------------------
_counter_cache = {}


def get_norm_factors(f, lo, hi, cut=GROOM_TAG, numjets_bin=BIN_NJETS):
    """(num_jets, <jet pT>) from the counters hist, or (None, None) if unusable."""
    key = (lo, hi, cut)
    if key in _counter_cache:
        return _counter_cache[key]

    cname = COUNTER_NAME % (lo, hi, cut)
    hc = f.Get(cname)
    result = (None, None)
    if not hc:
        print("WARNING: %s not found" % cname)
    elif hc.GetNbinsX() < max(numjets_bin, BIN_SUMPT):
        print("WARNING: %s has only %d bins" % (cname, hc.GetNbinsX()))
    elif hc.GetBinContent(numjets_bin) <= 0:
        print("WARNING: %s has zero jets in bin %d" % (cname, numjets_bin))
    else:
        njets = hc.GetBinContent(numjets_bin)
        result = (njets, hc.GetBinContent(BIN_SUMPT) / njets)

    _counter_cache[key] = result
    return result


def normalize_eec(hist, f, lo, hi, cut=GROOM_TAG, ptRL=False, numjets_bin=BIN_NJETS):
    """Scale an EEC hist in memory: 1/N_jets with "width" (plus log<pt>/<pt> for ptRL)."""
    njets, avg_pt = get_norm_factors(f, lo, hi, cut, numjets_bin=numjets_bin)
    if hist is None or njets is None:
        return None
    if ptRL:
        hist.Scale(math.log(avg_pt) / avg_pt)
    hist.Scale(1.0 / njets, "width")
    return hist


def load_eec_hists(f):
    """Fetch, normalize and style the full-EEC hist for each jet pT bin."""
    # stem = "hist_full_ptRL" if USE_PTRL else "hist_full"
    stem = "QA_hist_full"
    hists = []
    for i, (lo, hi) in enumerate(PT_BINS):
        # hname = "%s_jetpt%d_%d_%s" % (stem, lo, hi, GROOM_TAG)
        hname = "%s_pt%d_%d_wwjetpt" % (stem, lo, hi)
        h = f.Get(hname)
        if not h:
            print("WARNING: missing %s -- skipping" % hname)
            continue

        h = h.Clone("clone_%s" % hname)
        h.SetDirectory(0)
        INCL_BIN_NJETS=1
        if normalize_eec(h, f, lo, hi, GROOM_TAG, ptRL=USE_PTRL, numjets_bin=INCL_BIN_NJETS) is None:
            print("WARNING: no norm factors for %d-%d -- skipping" % (lo, hi))
            continue

        njets, avgpt = get_norm_factors(f, lo, hi, GROOM_TAG)
        print("%s: N_jets = %.0f, <p_T> = %.2f GeV/c" % (hname, njets, avgpt))
        style(h, COLORS[i % len(COLORS)], marker=MARKERS[i % len(MARKERS)], msize=1.1)
        hists.append((h, lo, hi))
    return hists


def plot_full_eec(filename):
    f = ROOT.TFile.Open(filename, "READ")
    if not f or f.IsZombie():
        raise IOError("could not open %s" % filename)

    hists = load_eec_hists(f)
    if not hists:
        raise RuntimeError("no histograms found")

    values = [h.GetBinContent(b) for h, _, _ in hists
              for b in range(1, h.GetNbinsX() + 1)]
    ymin, ymax = min(values), max(values)
    span = ymax - ymin
    ylo = 0.0 if ymin >= 0 else ymin - 0.10 * span   # anchor at 0 for a positive dist.
    yhi = ymax + 0.35 * span                         # headroom for legend/latex

    ax = hists[0][0].GetXaxis()
    xlow = ax.GetXmin() if ax.GetXmin() > 0 else ax.GetBinLowEdge(2)

    xtitle = "p_{T}R_{L}" if USE_PTRL else "R_{L}"
    c, frame = framed_canvas("c_eec", (xlow, ax.GetXmax()), (ylo, yhi), xtitle,
                             "#frac{1}{N_{jet}} #frac{dN_{EEC}}{d%s}" % xtitle,
                             logx=True, grid=False)
    frame.GetXaxis().SetMoreLogLabels()

    for h, _, _ in hists:
        h.Draw("PE same")
    make_legend((0.62, 0.57, 0.92, 0.87),
                [(h, "%d < p_{T,jet} < %d GeV/c" % (lo, hi), "lp") for h, lo, hi in hists],
                header="full EEC", textsize=0.032)
    draw_alice_label()

    c.RedrawAxis()
    outname = "%s/full_eec_all_ptbins.pdf" % OUT
    c.SaveAs(outname)
    f.Close()
    print("wrote %s" % outname)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    f = ROOT.TFile("%s/HistsDataCheckMerged.root" % CFG["input_base"])
    print("using file:", f.GetName())

    f.h_cuts.Print("all")                     # total / sel8 / global-track counts
    num_acc_events = f.h_cuts.GetBinContent(
        f.h_cuts.GetXaxis().FindBin("events_sel8_rct"))
    print("num_acc_events:", num_acc_events)

    # NOTE: plot_track_qa scales h_pt/h_eta/h_phi by 1/N_events in place, so the
    # figures below inherit that normalization (as in the original script).
    plot_track_qa(f, num_acc_events)
    plot_pt_panels(f)
    plot_jetpt_figures(f)
    plot_tracks_in_R_jets(f, 0.6)
    plot_tracks_in_R_jets(f, 0.4)
    plot_jet_eta(f)
    plot_eta_panels(f)
    plot_tracketa_trackpt(f)

    plot_full_eec(CFG["full_eec_file"])


if __name__ == "__main__":
    main()