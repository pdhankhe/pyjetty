#!/usr/bin/env python3
"""
Plot response matrices produced by make_rms.py.

Reads the ROOT file and produces:
  - 2D jet pt response (det vs part)
  - 2D groomed jet pt response (det vs part)
  - For each 6D THnSparse (AA, AB, BB, rad):
        * jet pt   det vs part   (axes 0,1)  -> 2D
        * RL       det vs part   (axes 2,3)  -> 2D
        * weight   det vs part   (axes 4,5)  -> 2D
        * 1D projections (det and part overlaid) for jet pt, RL, weight

    Usage:
        python PlotRMCrosschecks.py # -i crosschecks.root
"""

import os
import argparse
import ROOT
import math

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)          # annotation replaces the pad title
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(255)
ROOT.gStyle.SetPadGridX(1)
ROOT.gStyle.SetPadGridY(1)
ROOT.TH1.AddDirectory(False)

EEC_LABELS = ["AA", "AB", "BB", "rad"]

# text block drawn on every plot
INFO_LINES = [
    "MC pp #sqrt{s} = 5.36 TeV",
    "LHC26c5, JJ",
    "anti-#it{k}_{T} R = 0.4, |#eta_{jet}| < 0.5",
]

# 6D axis layout (must match make_rms.py):
#   0: pt_det   1: pt_part   2: RL_det   3: RL_part   4: w_det   5: w_part
AX_PT_DET,  AX_PT_PART  = 0, 1
AX_RL_DET,  AX_RL_PART  = 2, 3
AX_W_DET,   AX_W_PART   = 4, 5

# (name, det axis, part axis, axis title, logx on 1D, log x&y on 2D, logy on 1D)
PROJECTIONS = [
    ("jetpt",  AX_PT_DET, AX_PT_PART, "#it{p}_{T,jet} (GeV/#it{c})", True,  False, True),
    ("RL",     AX_RL_DET, AX_RL_PART, "#it{R}_{L}",                  True,  True,  False),
    ("weight", AX_W_DET,  AX_W_PART,  "weight",                      True,  True,  True),
]

_CANVAS_COUNTER = [0]


def _new_canvas(w=800, h=700):
    """Unique canvas name so ROOT never warns about replacing an old one."""
    _CANVAS_COUNTER[0] += 1
    return ROOT.TCanvas(f"c{_CANVAS_COUNTER[0]}", "c", w, h)


def _log_ok(axis):
    """A log axis needs a strictly positive lower edge."""
    return axis.GetXmin() > 0


def draw_text(lines, x=0.16, y=0.86, size=0.033, dy=0.045, font=42,
              color=ROOT.kBlack):
    """Draw lines of text in NDC coords. Returns the TLatex; keep the returned
    reference alive until the canvas is saved."""
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(font)
    latex.SetTextSize(size)
    latex.SetTextColor(color)
    latex.SetTextAlign(13)          # left / top
    for i, line in enumerate(lines):
        latex.DrawLatex(x, y - i * dy, line)
    return latex


def draw_hori_line(x1, x2, y1, color, linestyle, linewidth=1):
    line = ROOT.TLine(x1, y1, x2, y1)
    line.SetLineWidth(linewidth)
    line.SetLineColor(color)
    line.SetLineStyle(linestyle)
    line.Draw("SAME")
    return line


def draw_1d_overlay(h_det, h_part, xtitle, outpath,
                    logx=False, logy=False, lines=None, extra_lines=None, matchvsall=False, showmorexlab=False):
    """Overlay det and part 1D projections and save."""
    c = _new_canvas(800, 650)
    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.06)

    if logx and _log_ok(h_det.GetXaxis()):
        c.SetLogx()
        if showmorexlab:
            h_det.GetXaxis().SetMoreLogLabels()  # <--- Forces 10, 20, 50, 100, 200... labels
            h_det.GetXaxis().SetNoExponent()
    if logy:
        c.SetLogy()

    h_det.SetLineColor(ROOT.kAzure + 1)
    h_det.SetLineWidth(2)
    h_det.SetMarkerColor(ROOT.kAzure + 1)
    h_det.SetMarkerStyle(20)
    h_det.SetMarkerSize(0.8)

    h_part.SetLineColor(ROOT.kRed + 1)
    h_part.SetLineWidth(2)
    h_part.SetMarkerColor(ROOT.kRed + 1)
    h_part.SetMarkerStyle(24)
    h_part.SetMarkerSize(0.8)

    ymax = max(h_det.GetMaximum(), h_part.GetMaximum())

    if logy and ymax > 0:
        # log scale needs a strictly positive range; base the minimum on the
        # smallest nonzero bin content across both histograms
        mins = [h.GetMinimum(0.0) for h in (h_det, h_part)]
        mins = [m for m in mins if m > 0]
        ymin = min(mins) * 0.5 if mins else 1e-3
        ymax_draw = 1e3 * ymax          # headroom for the text block
        if ymax_draw <= ymin:           # guard degenerate case
            ymax_draw = ymin * 100
        h_det.SetMinimum(ymin)
        h_det.SetMaximum(ymax_draw)
    else:
        if logy:                        # requested log but nothing positive
            c.SetLogy(0)                # turn it back off, avoids painter errors
        h_det.SetMinimum(0.0)
        h_det.SetMaximum(1.6 * ymax if ymax > 0 else 1.0)

    h_det.GetXaxis().SetTitle(xtitle)
    h_det.GetYaxis().SetTitle("counts")
    h_det.GetYaxis().SetTitleOffset(1.4)

    h_det.Draw("E1")
    h_part.Draw("E1 SAME")

    leg = ROOT.TLegend(0.68, 0.74, 0.92, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.033)
    if matchvsall:
        leg.AddEntry(h_det, "matched", "lep")
        leg.AddEntry(h_part, "all", "lep")
    else:
        leg.AddEntry(h_det, "detector", "lep")
        leg.AddEntry(h_part, "particle", "lep")
    leg.Draw()

    block = list(lines) if lines else []
    if extra_lines:
        block += list(extra_lines)
    txt = draw_text(block) if block else None

    c.RedrawAxis()
    c.SaveAs(outpath)
    c.Close()


def process_sparse(hs, lab, grooming_str, outdir, lines=None, consolidated_dict=None):
    """Make 2D det-vs-part and 1D det/part overlays for each observable."""
    for name, ax_det, ax_part, xtitle, logx1d, log2d, logy1d in PROJECTIONS:
        tag = [f"{lab}, {name} response ({grooming_str} jetpt bins)"]

        # ---- 2D det vs part ----
        # THnSparse::Projection(y, x) -> TH2 with x = second arg, y = first arg
        h2 = hs.Projection(ax_part, ax_det)          # x = det, y = part
        h2.SetName(f"{lab}_{name}_2D")
        h2.SetTitle("")
        h2.GetXaxis().SetTitle(f"{xtitle} (det)")
        h2.GetYaxis().SetTitle(f"{xtitle} (part)")

        # Store for consolidated plotting across labels
        if consolidated_dict is not None:
            if name not in consolidated_dict:
                consolidated_dict[name] = {}
            consolidated_dict[name][lab] = h2

        draw_th2(
            h2,
            os.path.join(outdir, f"resp2D_{grooming_str}_{lab}_{name}.pdf"),
            logx=log2d, logy=log2d, logz=True,
            lines=lines, extra_lines=tag,
        )

        # ---- 1D projections (det and part) ----
        h_det = hs.Projection(ax_det)
        h_det.SetName(f"{lab}_{name}_det")
        h_det.SetTitle("")
        h_part = hs.Projection(ax_part)
        h_part.SetName(f"{lab}_{name}_part")
        h_part.SetTitle("")

        draw_1d_overlay(
            h_det, h_part, xtitle,
            os.path.join(outdir, f"proj1D_{lab}_{name}.pdf"),
            logx=logx1d, logy=logy1d,
            lines=lines, extra_lines=[f"{lab}, {name}"],
        )



def plot_consolidated_2d(consolidated_h2s, eec_labels, projections, grooming_str, outdir, lines=None):
    """Generate 4-panel consolidated 2D response plots with aligned Z scales for each observable."""
    for name, _, _, xtitle, _, log2d, _ in projections:
        if name not in consolidated_h2s or not consolidated_h2s[name]:
            continue

        c_cons = ROOT.TCanvas(f"c_resp2D_{grooming_str}_{name}_consolidated", f"Consolidated 2D Responses - {name}", 1600, 1400)
        c_cons.Divide(2, 2)

        # Determine global z-range across the labels for scale alignment
        max_z = max((h2.GetMaximum() for h2 in consolidated_h2s[name].values()), default=-1)
        min_z = min((h2.GetMinimum() for h2 in consolidated_h2s[name].values()), default=-1)

        txt_refs = []  # <--- Keeps TLatex instances alive for ALL pads until canvas saves

        for i, lab in enumerate(eec_labels, 1):
            if lab not in consolidated_h2s[name]:
                continue
                
            pad = c_cons.cd(i)
            pad.SetLogx(log2d)
            pad.SetLogy(log2d)
            pad.SetLogz(True)
            pad.SetRightMargin(0.15)

            h2 = consolidated_h2s[name][lab]
            if max_z > 0:
                h2.SetMaximum(max_z * 1.05)
            if min_z > 0:
                h2.SetMinimum(min_z / 1.05)
                
            h2.SetTitle(f"{lab} - {name}")
            h2.Draw("COLZ")

            if i == 1:
                block = list(lines) if lines else []
                block += [f"{lab} {name} pairs"]
            else:
                block = [f"{lab} {name} pairs"]
                
            txt = draw_text(block) if block else None
            if txt:
                txt_refs.append(txt)  # Store reference in list

        c_cons.SaveAs(os.path.join(outdir, f"consolidated_resp2D_{grooming_str}_{name}.pdf"))
        c_cons.Close()
        
        
def PlotGiven1D(hist, outpath, xtitle="", ytitle="", lines=None, extra_lines=None, ytext=0.86, logx=False, logy=False, horiline=False, yscale01=False, showmorexlab=False):
    # hist = f.Get(hist_name)
    c = _new_canvas(800, 650)

    if logx and _log_ok(hist.GetXaxis()):
        c.SetLogx()
        if showmorexlab:
            hist.GetXaxis().SetMoreLogLabels()  # <--- Forces 10, 20, 50, 100, 200... labels
            hist.GetXaxis().SetNoExponent()
    if logy:
        c.SetLogy()

    if yscale01: #set the y axis from 0 to 1
        hist.SetMinimum(0.0)
        hist.SetMaximum(1.0)
    
    if xtitle != "":
        hist.GetXaxis().SetTitle(xtitle)
    if ytitle != "":
        hist.GetYaxis().SetTitle(ytitle)

    hist.Draw("E1")
    
    block = list(lines) if lines else []
    if extra_lines:
        block += list(extra_lines)
    txt = draw_text(block, y=ytext) if block else None #x=0.16, y=0.86, size=0.033, dy=0.045

    
    line = None                      # keep alive until SaveAs
    if horiline:
        line = draw_hori_line(hist.GetXaxis().GetXmin(),
                              hist.GetXaxis().GetXmax(),
                              0, ROOT.kGray + 2, 2, linewidth=2)

    c.RedrawAxis()
    c.SaveAs(outpath)
    c.Close()

    # leg = ROOT.TLegend(0.68, 0.74, 0.92, 0.88)
    # leg.SetBorderSize(0)
    # leg.SetFillStyle(0)
    # leg.SetTextFont(42)
    # leg.SetTextSize(0.033)
    # leg.AddEntry(h_det, "detector", "lep")
    # leg.AddEntry(h_part, "particle", "lep")
    # leg.Draw()


def PlotGiven2D(hist2D, outpath, xtitle="", ytitle="", lines=None, extra_lines=None, logx=False, logy=False, logz=False, zscale01=False):
    #, logx=False, logy=False, logz=True, zmin=None, zmax=None,

    c = _new_canvas(800, 650) #(800, 700)
    c.SetRightMargin(0.15)
    c.SetLeftMargin(0.13)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.06)

    hist2D.GetXaxis().SetTitleOffset(1.1)
    hist2D.GetYaxis().SetTitleOffset(1.3)
    
    if logx and _log_ok(hist2D.GetXaxis()):
        c.SetLogx()
    if logy and _log_ok(hist2D.GetYaxis()):
        c.SetLogy()
    if logz and hist2D.GetEntries() > 0:
        c.SetLogz()
        # if zmin is None:
        #     zmin = hist2D.GetMinimum(0.0)      # smallest bin content strictly > 0
        #     if zmin <= 0:                 # everything zero/negative
        #         zmin = 1e-3
    # if zmin is not None:
    #     hist2D.SetMinimum(zmin)
    # if zmax is not None:
    #     hist2D.SetMaximum(zmax)
        
    if zscale01: #set the z axis from 0 to 1
        hist2D.SetMinimum(0.0)
        hist2D.SetMaximum(1.0)
    
    if xtitle != "":
        hist.GetXaxis().SetTitle(xtitle)
    if ytitle != "":
        hist.GetYaxis().SetTitle(ytitle)

    hist2D.Draw("COLZ")

    block = list(lines) if lines else []
    if extra_lines:
        block += list(extra_lines)
    txt = draw_text(block) if block else None

    c.RedrawAxis()
    c.SaveAs(outpath)
    c.Close()



def plot_QA(f, grooming_str, outdir):
    '''
    Residual plots are 2D: xaxis: part, yaxis: (part-det)/part
    '''

    residual_dir = os.path.join(outdir, "residuals")
    os.makedirs(residual_dir, exist_ok=True)

    # First, jet pt
    residual_jetpt = f.Get(f"res_jetpt_{grooming_str}")
    CalculateAndPlotResAndScale(residual_jetpt, "jetpt", "", f"jetpt", residual_dir, logx=True)

    # Then, RL & weights
    residuals_rl = {}
    residuals_w = {}
    for lab in EEC_LABELS:
        res_rl = f.Get(f"res_rl_{lab}")
        res_w = f.Get(f"res_w_{lab}")
        CalculateAndPlotResAndScale(res_rl, "rl", lab, f"rl_{lab}", residual_dir, logx=True)
        CalculateAndPlotResAndScale(res_w, "weight", lab, f"weight_{lab}", residual_dir, logx=True)
        # residuals_rl[lab] = res_rl
        # residuals_w[lab] = res_w

    # Finally, track efficiencies/purities/resolution
    trk_eff_num = f.Get("trk_eff_num")
    trk_eff_den = f.Get("trk_eff_den")
    trk_pur_num = f.Get("trk_pur_num")
    trk_pur_den = f.Get("trk_pur_den")



    trk_res_pt = f.Get("trk_res_pt")
    CalculateAndPlotResAndScale(trk_res_pt, "trackpt", "", f"trackpt", residual_dir, logx=True)


def CalculateAndPlotResAndScale(h2D, name, lab, hname, outdir,
                                logx=False, min_entries=100,
                                fit_range=(-1.0, 1.0),
                                jes_range=(-0.5, 0.5),
                                jer_range=(0.0, 0.5),
                                draw_slices=True, lines=None):
    """
    Input: 2D residual plot, xaxis: part, yaxis: (part-det)/part
    Returns (hJER, hJES):
        hJER = sigma of the Gaussian fit in each particle-pT slice
        hJES = mean  of the Gaussian fit in each particle-pT slice
    """

    if not h2D or h2D.GetEntries() == 0:
        print(f"WARNING: skipping {hname}, histogram missing or empty")
        return None, None

    base = h2D.GetName()
    os.makedirs(outdir, exist_ok=True)

    slice_pdf = None
    if draw_slices:
        slice_dir = os.path.join(outdir, "slices")
        os.makedirs(slice_dir, exist_ok=True)
        slice_pdf = os.path.join(slice_dir, f"slices_{hname}.pdf")

    hMean, hSigma = _fit_slices_with_plots(
        h2D, min_entries=min_entries, fit_range=fit_range,
        slice_pdf=slice_pdf, lines=lines,
    )

    hJES = hMean.Clone(f"hJES_{hname}")
    hJER = hSigma.Clone(f"hJER_{hname}")

    os.makedirs(outdir, exist_ok=True)
    PlotGiven1D(hJER, os.path.join(outdir, f"JER_{hname}.pdf"),
                logx=logx, logy=False, horiline=False)
    PlotGiven1D(hJES, os.path.join(outdir, f"JES_{hname}.pdf"),
                logx=logx, logy=False, horiline=False)



def _fit_slices_with_plots(h2D, min_entries=100, fit_range=(-1.0, 1.0),
                           niter=4, nsig=2.0, slice_pdf=None,
                           ncols=4, nrows=3, logy=True, lines=None):
    """
    Fit a Gaussian to every y-projection of h2D, iterating the fit window to
    +-nsig*sigma around the current mean.

    Returns (hMean, hSigma) with fit parameters and their errors, and
    optionally writes a multi-page PDF showing each slice with its fit.
    """
    base = h2D.GetName()

    hMean = h2D.ProjectionX(f"{base}_fitmean")
    hMean.Reset()
    hMean.SetDirectory(0)
    hSigma = hMean.Clone(f"{base}_fitsigma")
    hSigma.SetDirectory(0)

    keep = []            # ROOT objects must stay alive until Print()
    c = None
    ipad = 0
    npads = ncols * nrows

    if slice_pdf:
        os.makedirs(os.path.dirname(slice_pdf) or ".", exist_ok=True)
        c = _new_canvas(400 * ncols, 340 * nrows)
        c.Print(slice_pdf + "[")     # open multi-page document

    n_ok = n_skip = n_fail = 0

    for i in range(1, h2D.GetNbinsX() + 1):
        p = h2D.ProjectionY(f"_py_{base}_{i}", i, i)
        p.SetDirectory(0)
        p.SetTitle("")

        lo_edge = h2D.GetXaxis().GetBinLowEdge(i)
        hi_edge = h2D.GetXaxis().GetBinUpEdge(i)

        if p.GetEntries() < min_entries:
            n_skip += 1
            continue

        m, s = p.GetMean(), p.GetRMS()
        if s <= 0:
            n_fail += 1
            continue

        g = ROOT.TF1(f"_g_{base}_{i}", "gaus", *fit_range)
        g.SetParameters(p.GetMaximum(), m, s)
        ok = False
        for _ in range(niter):
            lo = max(m - nsig * s, fit_range[0])
            hi = min(m + nsig * s, fit_range[1])
            if hi <= lo:
                break
            if int(p.Fit(g, "QNRS", "", lo, hi)) != 0:
                break
            m, s = g.GetParameter(1), abs(g.GetParameter(2))
            ok = True

        if not ok:
            n_fail += 1
        else:
            n_ok += 1
            hMean.SetBinContent(i, g.GetParameter(1))
            hMean.SetBinError(i, g.GetParError(1))
            hSigma.SetBinContent(i, abs(g.GetParameter(2)))
            hSigma.SetBinError(i, g.GetParError(2))

        # ---------------- draw this slice ----------------
        if not slice_pdf:
            continue

        ipad += 1
        pad = c.cd(ipad if ipad <= npads else 1)
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.04)
        pad.SetTopMargin(0.06)
        pad.SetBottomMargin(0.13)

        p.SetMarkerStyle(20)
        p.SetMarkerSize(0.6)
        p.SetLineColor(ROOT.kBlack)
        p.GetXaxis().SetTitle(h2D.GetYaxis().GetTitle() or "(part-det)/part")
        p.GetXaxis().SetTitleSize(0.05)
        p.GetXaxis().SetLabelSize(0.045)
        p.GetYaxis().SetTitle("counts")
        p.GetYaxis().SetTitleSize(0.05)
        p.GetYaxis().SetTitleOffset(1.3)
        p.GetYaxis().SetLabelSize(0.045)

        if logy and p.GetMaximum() > 0:
            pad.SetLogy()
            p.SetMinimum(0.5)
            p.SetMaximum(20 * p.GetMaximum())
        else:
            p.SetMinimum(0.0)
            p.SetMaximum(1.7 * p.GetMaximum())

        p.Draw("E1")

        if ok:
            gd = g.Clone(f"_gd_{base}_{i}")
            gd.SetRange(max(m - nsig * s, fit_range[0]),
                        min(m + nsig * s, fit_range[1]))
            gd.SetLineColor(ROOT.kRed + 1)
            gd.SetLineWidth(2)
            gd.Draw("SAME")
            keep.append(gd)

            # dashed extrapolation over the full range, to expose the tails
            ge = g.Clone(f"_ge_{base}_{i}")
            ge.SetRange(*fit_range)
            ge.SetLineColor(ROOT.kRed + 1)
            ge.SetLineStyle(2)
            ge.SetLineWidth(1)
            ge.Draw("SAME")
            keep.append(ge)

        # zero reference
        zl = ROOT.TLine(0.0, p.GetMinimum(), 0.0, p.GetMaximum())
        zl.SetLineColor(ROOT.kGray + 2)
        zl.SetLineStyle(3)
        zl.Draw("SAME")
        keep.append(zl)

        ndf = g.GetNDF()
        chi2ndf = g.GetChisquare() / ndf if ndf > 0 else float("nan")
        info = [
            f"bin {i}: [{lo_edge:.3g}, {hi_edge:.3g}]",
            f"N = {p.GetEntries():.0f}",
            f"#mu = {g.GetParameter(1):+.4f} #pm {g.GetParError(1):.4f}",
            f"#sigma = {abs(g.GetParameter(2)):.4f} #pm {g.GetParError(2):.4f}",
            f"#chi^{{2}}/ndf = {chi2ndf:.2f}",
            f"RMS = {p.GetRMS():.4f}",
        ]
        if not ok:
            info.append("#color[2]{FIT FAILED}")
        keep.append(draw_text(info, x=0.60, y=0.90, size=0.045, dy=0.058))

        if ipad == npads:
            c.Print(slice_pdf)
            c.Clear()
            c.Divide(ncols, nrows)
            ipad = 0
            keep = []

    if slice_pdf:
        if ipad > 0:
            c.Print(slice_pdf)
        c.Print(slice_pdf + "]")     # close document
        c.Close()

    print(f"    [{base}] slices: {n_ok} fitted, {n_skip} below "
          f"{min_entries} entries, {n_fail} failed")

    return hMean, hSigma


def GetHist(f, histname):
    print("histname", histname)
    hist = f.Get(histname) #f.histname
    return hist

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input",
                    default="/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rm_crosscheck/58272568/crosschecks_merged.root",
                    help="input ROOT file (from make_rms.py)")
    ap.add_argument("-o", "--outdir",
                    default="/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots/response_matrices/crosschecks",
                    help="output directory for plots")
    args = ap.parse_args()

    
    os.makedirs(args.outdir, exist_ok=True)
    # os.makedirs(os.path.join(updated_outdir, "effpur"), exist_ok=True)

    f = ROOT.TFile.Open(args.input)
    if not f or f.IsZombie():
        raise RuntimeError(f"could not open {args.input}")

    
    hnevents = f.Get("hNevents")
    print("There are a total of ", hnevents.GetBinContent(hnevents.FindBin(0)), "events in anch mc.")

    suffix = ["det", "part"]

    # 2D histograms
    # the brackets get replaced with det/part
    # histogram name, x title, y title, logx, logy, logz
    hist2D_list = [
        ("h_pair_pt_{}",   "#it{p}_{T,i} [GeV/#it{c}]", "#it{p}_{T,j} [GeV/#it{c}]", False, False, True),
        ("h_pair_pid_{}",  "particle from subjet A",    "particle from subjet B",    False, False, True),
    ]
    for h2D_str, xtitle, ytitle, logx, logy, logz in hist2D_list:
        for suf in suffix:
            h2D_name = h2D_str.format(suf)
            print(h2D_name)
            h2D = f.Get(h2D_name)
            outpath = os.path.join(args.outdir, f"{h2D_name}.pdf")
            PlotGiven2D(h2D, outpath, logx=logx, logy=logy, logz=logz) #, xtitle="", ytitle="", lines=None, extra_lines=None)


    # 1D histograms
    hist1D_list = [ "h_n_const_orig_{}", "h_n_const_groomed_{}", "h_n_const_removed_{}", "h_jet_pt_orig_{}", "h_jet_pt_groomed_{}", "h_pair_dr_{}", "h_pair_mass_{}" ]
    # histogram name, x title, y title logx, logy
    hist1D_list = [
        ("h_n_const_orig_{}",     "N_{constituents}",                 "counts",  False, True),
        ("h_n_const_groomed_{}",  "N_{groomed constituents}",         "counts",  False, True),
        ("h_n_const_removed_{}",  "N_{removed}",                      "counts",  False, True),
        ("h_jet_pt_orig_{}",      "#it{p}_{T,jet} [GeV/#it{c}]",      "counts",  False, True),
        ("h_jet_pt_groomed_{}",   "#it{p}_{T,gr. jet} [GeV/#it{c}]",  "counts",  False, True),
        ("h_pair_dr_{}",          "#Delta R",                         "counts",  False, True),
        ("h_pair_mass_{}",        "Pair Mass [GeV]",                  "counts",  False, True),
    ]
    for h_str, xtitle, ytitle, logx, logy in hist1D_list:
        hists = []
        for suf in suffix:
            h_name = h_str.format(suf)
            h = f.Get(h_name)
            hists.append(h)
            
            # PlotGiven1D(h, outpath) #, xtitle="", ytitle="")
        output_name = h_str.format("part_vs_det")
        outpath = os.path.join(args.outdir, f"{output_name}.pdf")
        draw_1d_overlay(hists[0], hists[1], xtitle, outpath, logx=logx, logy=logy) #, lines=None, extra_lines)


if __name__ == "__main__":
    main()