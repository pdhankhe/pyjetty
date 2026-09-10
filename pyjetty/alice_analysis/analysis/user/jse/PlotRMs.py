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
        python plot_rms.py -i response.root -o plots
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

def draw_th2(h, outpath, logx=False, logy=False, logz=True,
             zmin=None, zmax=None, lines=None, extra_lines=None):
    """Draw a TH2 as a colz plot and save."""
    c = _new_canvas(800, 700)
    c.SetRightMargin(0.15)
    c.SetLeftMargin(0.13)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.06)

    if logx and _log_ok(h.GetXaxis()):
        c.SetLogx()
    if logy and _log_ok(h.GetYaxis()):
        c.SetLogy()
    if logz and h.GetEntries() > 0:
        c.SetLogz()
        if zmin is None:
            zmin = h.GetMinimum(0.0)      # smallest bin content strictly > 0
            if zmin <= 0:                 # everything zero/negative
                zmin = 1e-3
    if zmin is not None:
        h.SetMinimum(zmin)
    if zmax is not None:
        h.SetMaximum(zmax)

    h.GetXaxis().SetTitleOffset(1.1)
    h.GetYaxis().SetTitleOffset(1.3)
    h.Draw("COLZ")

    block = list(lines) if lines else []
    if extra_lines:
        block += list(extra_lines)
    txt = draw_text(block) if block else None   # local ref keeps it alive

    c.RedrawAxis()
    c.SaveAs(outpath)
    c.Close()


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

        if name in ("RL", "weight"):
            CalculateResidual(h_part, h_det, name, lab, f"hresidual_{name}_{lab}", outdir, logx=logx1d)
        # CalculateJER(h2, name, lab, f"hjer_{name}_{lab}", outdir, logx=logx1d)
        # CalculateJES(h2, name, lab, f"hjes_{name}_{lab}", outdir, logx=logx1d)


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


def PlotGiven2D(hist2D, outpath, xtitle="", ytitle="", lines=None, extra_lines=None, zscale01=False):
    c = _new_canvas(800, 650)
    
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


def CalculateResidual(hpart, hdet, name, lab, hname, outdir, logx=False, correlated=False, plot=True):
    """
    Bin-by-bin residual (hpart - hdet)/hpart.

    correlated=False assumes hpart and hdet are statistically independent.
    correlated=True is appropriate when both are filled from the same
    matched jets, in which case the naive propagation overestimates the
    uncertainty; there we only propagate the det-level error.
    RETURNS: the TH2D (i.e. truth jetpt in x-axis, (jetpt_truth - jetpt_det) / jetpt_truth in y-axis) for JER/JES calculation
    """
    if hpart.GetNbinsX() != hdet.GetNbinsX():
        raise ValueError("CalculateResidual: incompatible binning "
                         f"({hpart.GetNbinsX()} vs {hdet.GetNbinsX()})")

    # Make 1D histogram version
    hresidual = hpart.Clone(hname)
    hresidual.SetDirectory(0)
    hresidual.Reset()
    hresidual.SetTitle(hname)
    hresidual.GetYaxis().SetTitle(f"{lab} (part #minus det)/part")


    for ibin in range(hpart.GetNcells()):   # includes under/overflow
        p  = hpart.GetBinContent(ibin)
        d  = hdet.GetBinContent(ibin)
        ep = hpart.GetBinError(ibin)
        ed = hdet.GetBinError(ibin)

        if p == 0.:
            hresidual.SetBinContent(ibin, 0.)
            hresidual.SetBinError(ibin, 0.)
            continue

        r = (p - d) / p

        # r = 1 - d/p  ->  dr/dd = -1/p,  dr/dp = d/p^2
        if correlated:
            err = abs(ed / p)
        else:
            err = math.sqrt((ed / p)**2 + (d * ep / (p * p))**2)

        # Fill 1D
        hresidual.SetBinContent(ibin, r)
        hresidual.SetBinError(ibin, err)

    
    if plot:
        PlotGiven1D(hresidual, os.path.join(outdir, f"residual_{lab}_{name}.pdf"), logx=logx, logy=False, horiline=True)

    # hresidual.ResetStats()
    return hresidual


# def CalculateJER(h2D, name, lab, hname, outdir, logx=False):
#     """
#     Calculate Jet Energy Resolution: sigma(pT_det) / mean(pT_det) per pT_truth bin using TProfile.
#     h2D axes: X = det-level (reco), Y = truth-level (gen)
#     """

#     # 1. Profile X (pT_det) along Y (pT_truth)
#     # Option "s" calculates spread/std-dev per bin instead of error-on-mean
#     prof = h2D.ProfileY(f"prof_{hname}", 1, -1, "s")
#     prof.SetDirectory(0)

#     # 2. Extract binning array to match truth pT (Y-axis of h2D)
#     y_axis = h2D.GetYaxis()
#     n_bins = y_axis.GetNbins()
#     y_bins = y_axis.GetXbins().GetArray()

#     # Create local 1D histograms detached from ROOT gDirectory
#     h_jer = ROOT.TH1D(f"h_jer_{name}", f"JER - {lab};p_{{T, truth}} [GeV];#sigma(p_{{T, det}}) / #langle p_{{T, det}} #rangle", n_bins, y_bins)
#     h_jer.SetDirectory(0)

#     # 3. Calculate JES and JER per truth bin
#     for i in range(1, n_bins + 1):
#         pt_truth = y_axis.GetBinCenter(i)
#         entries = prof.GetBinEntries(i)

#         if entries < 10:  # Skip bins with low statistics
#             continue

#         mean = prof.GetBinContent(i)      # <pT_det>
#         sigma = prof.GetBinError(i)       # std-dev(pT_det) from "s" option
#         mean_err = sigma / math.sqrt(entries) if entries > 0 else 0.0

#         if mean <= 0 or pt_truth <= 0:
#             continue

#         # JER = sigma(pT_det) / <pT_det>
#         jer = sigma / mean
#         # Error propagation for sigma / mean assuming normal population
#         sigma_err = sigma / math.sqrt(2.0 * entries) if entries > 1 else 0.0
#         jer_err = jer * math.sqrt((sigma_err / sigma)**2 + (mean_err / mean)**2)
        
#         h_jer.SetBinContent(i, jer)
#         h_jer.SetBinError(i, jer_err)

#     # 4. Parametrize resolution with standard N/S/C curve
#     x_min = h_jer.GetXaxis().GetXmin()
#     x_max = h_jer.GetXaxis().GetXmax()

#     jer_fit = ROOT.TF1(f"jer_fit_{name}", "sqrt(([0]/x)^2 + ([1]/sqrt(x))^2 + [2]^2)", x_min, x_max)
#     jer_fit.SetParNames("N (Noise)", "S (Stochastic)", "C (Constant)")
#     jer_fit.SetParameters(1.0, 0.8, 0.05)

#     h_jer.Fit(jer_fit, "RQ0")

#     # 5. Plot outputs
#     PlotGiven1D(h_jer, os.path.join(outdir, f"JER_{lab}_{name}.pdf"), logx=logx, logy=False, horiline=False)

#     return h_jer, jer_fit


# def CalculateJES(h2D, name, lab, hname, outdir, logx=False):
#     """
#     Calculate Jet Energy Scale: mean(pT_det) / pT_truth_center per pT_truth bin using TProfile.
#     h2 axes: X = det-level (reco), Y = truth-level (gen)
#     """

#     # 1. Profile X (pT_det) along Y (pT_truth)
#     # Option "s" calculates spread/std-dev per bin instead of error-on-mean
#     prof = h2D.ProfileY(f"prof_{hname}", 1, -1, "s")
#     prof.SetDirectory(0)

#     # 2. Extract binning array to match truth pT (Y-axis of h2D)
#     y_axis = h2D.GetYaxis()
#     n_bins = y_axis.GetNbins()
#     y_bins = y_axis.GetXbins().GetArray()

#     h_jes = ROOT.TH1D(f"h_jes_{name}", f"JES - {lab};p_{{T, truth}} [GeV];#langle p_{{T, det}} #rangle / p_{{T, truth}}", n_bins, y_bins)
#     h_jes.SetDirectory(0)

#     # 3. Calculate JES and JER per truth bin
#     for i in range(1, n_bins + 1):
#         pt_truth = y_axis.GetBinCenter(i)
#         entries = prof.GetBinEntries(i)

#         if entries < 10:  # Skip bins with low statistics
#             continue

#         mean = prof.GetBinContent(i)      # <pT_det>
#         sigma = prof.GetBinError(i)       # std-dev(pT_det) from "s" option
#         mean_err = sigma / math.sqrt(entries) if entries > 0 else 0.0

#         if mean <= 0 or pt_truth <= 0:
#             continue

#         # JES = <pT_det> / pT_truth
#         jes = mean / pt_truth
#         jes_err = mean_err / pt_truth
#         h_jes.SetBinContent(i, jes)
#         h_jes.SetBinError(i, jes_err)

    
#     # 5. Plot outputs
#     PlotGiven1D(h_jes, os.path.join(outdir, f"JES_{lab}_{name}.pdf"), logx=logx, logy=False, horiline=True)

#     return h_jes


def GetHist(f, histname):
    print("histname", histname)
    hist = f.Get(histname) #f.histname
    return hist

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input",
                    # default="/rstorage/alice/AnalysisResults/blianggi/jse/rms/1836481/reponse_merged.root",
                    default="/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/57911278/response_{}_merged_partial.root",
                    # default="/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/1836481/reponse_merged.root",
                    # default="/global/cfs/cdirs/alice/blianggi/mypyjetty/analysis/testing/response.root",
                    help="input ROOT file (from make_rms.py)")
    ap.add_argument("-o", "--outdir",
                    # default="/software/users/blianggi/mypyjetty/storage/jse/plots/response_matrices",
                    default="/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots/response_matrices",
                    help="output directory for plots")
    ap.add_argument("--groomed-bins", default=True, help="True if looking at groomed bins, False if looking at ungroomed bins")
    args = ap.parse_args()

    # Groomed vs ungroomed bins
    gr_bins_bool = bool(args.groomed_bins)
    grooming_str = "groomed" if gr_bins_bool else "ungroomed"
    
    updated_outdir = os.path.join(args.outdir, f"{grooming_str}_bins")
    os.makedirs(updated_outdir, exist_ok=True)
    os.makedirs(os.path.join(updated_outdir, "effpur"), exist_ok=True)

    # fill the {} placeholder only if one is present, so an explicit -i still works
    infile = args.input.format(grooming_str) if "{}" in args.input else args.input
    f = ROOT.TFile.Open(infile)
    if not f or f.IsZombie():
        raise RuntimeError(f"could not open {infile}")
    

    # -----------------------------------------------------------------------
    # 2D jet pt and groomed jet pt responses
    # -----------------------------------------------------------------------
    h_jetpt = f.Get(f"resp_jetpt_{grooming_str}")
    if h_jetpt:

        jetptstr = "#it{p}_{T,g}" if gr_bins_bool else "it{p}_{T}"
        ext_line = ["groomed jet #it{p}_{T} response"] if gr_bins_bool else ["jet #it{p}_{T} response"]
        
        h_jetpt.SetTitle("")
        h_jetpt.GetXaxis().SetTitle(f"{jetptstr}^{{det}} (GeV/#it{{c}})")
        h_jetpt.GetYaxis().SetTitle(f"{jetptstr}^{{part}} (GeV/#it{{c}})")
        draw_th2(h_jetpt, os.path.join(updated_outdir, "resp2D_jetpt.pdf"),
                 logx=False, logy=False, logz=True,
                 lines=INFO_LINES, extra_lines=ext_line)
        
        # CalculateResidual(h_jetpt.ProjectionY(), h_jetpt.ProjectionX(), "jetpt", "", f"hresidual_{grooming_str}_jetpt", updated_outdir, logx=logx1d) #part, det
        # CalculateJER(h_jetpt, "jetpt", "", f"hjer_{grooming_str}_jetpt", updated_outdir, logx=True)
        # CalculateJES(h_jetpt, "jetpt", "", f"hjes_{grooming_str}_jetpt", updated_outdir, logx=True)
    else:
        print(f"WARNING: resp_jetpt_{grooming_str} not found")

    
    # Plot efficiency/purity plots
    draw_1d_overlay(GetHist(f, f"jet_match_gen_eff_num_{grooming_str}"), GetHist(f, f"jet_all_gen_eff_den_{grooming_str}"), "#it{p}_{T}^{tr} (GeV/#it{c})", os.path.join(updated_outdir, f"effpur/jet_gen_{grooming_str}_matchedvsall.pdf"),
                    lines=INFO_LINES, extra_lines=["truth jet p_{T}"], logx=True, logy=True, matchvsall=True, showmorexlab=True)
    PlotGiven1D(GetHist(f, f"jet_efficiency_{grooming_str}_new"), os.path.join(updated_outdir, f"effpur/jet_efficiency_{grooming_str}_new.pdf"), 
                xtitle="#it{p}_{T}^{tr} (GeV/#it{c})", ytitle="Efficiency = matched gen/all gen", lines=INFO_LINES, extra_lines=["truth jet p_{T}"], ytext=0.3, logx=True, yscale01=True, showmorexlab=True)
    draw_1d_overlay(GetHist(f, f"jet_match_rec_pur_num_{grooming_str}"), GetHist(f, f"jet_all_rec_pur_den_{grooming_str}"), "#it{p}_{T}^{det} (GeV/#it{c})", os.path.join(updated_outdir, f"effpur/jet_rec_{grooming_str}_matchedvsall.pdf"),
                    lines=INFO_LINES, extra_lines=["detector jet p_{T}"], logx=True, logy=True, matchvsall=True, showmorexlab=True)
    PlotGiven1D(GetHist(f, f"jet_purity_{grooming_str}_new"), os.path.join(updated_outdir, f"effpur/jet_purity_{grooming_str}_new.pdf"), 
                xtitle="#it{p}_{T}^{det} (GeV/#it{c})", ytitle="Purity = matched rec/all rec", lines=INFO_LINES, extra_lines=["detector jet p_{T}"], ytext=0.3, logx=True, yscale01=True, showmorexlab=True)

    if gr_bins_bool:

        PlotGiven2D(f.lund_matched_gen, os.path.join(updated_outdir, "effpur/lund_matched_gen.pdf"), lines=INFO_LINES, extra_lines=["matched truth splittings"])
        PlotGiven2D(f.lund_all_gen, os.path.join(updated_outdir, "effpur/lund_all_gen.pdf"), lines=INFO_LINES, extra_lines=["all truth splittings"])
        PlotGiven2D(f.lund_matched_rec, os.path.join(updated_outdir, "effpur/lund_matched_rec.pdf"), lines=INFO_LINES, extra_lines=["all detector splittings"])
        PlotGiven2D(f.lund_all_rec, os.path.join(updated_outdir, "effpur/lund_all_rec.pdf"), lines=INFO_LINES, extra_lines=["matched detector splittings"])

        PlotGiven2D(f.lund_split_efficiency_new, os.path.join(updated_outdir, "effpur/lund_split_efficiency_new.pdf"), lines=INFO_LINES, extra_lines=["truth splittings"], zscale01=True)
        PlotGiven2D(f.lund_split_purity_new, os.path.join(updated_outdir, "effpur/lund_split_purity_new.pdf"), lines=INFO_LINES, extra_lines=["detector splittings"], zscale01=True)

        for lab in EEC_LABELS:

            draw_1d_overlay(GetHist(f, f"pair_match_gen_eff_num_{lab}"), GetHist(f, f"pair_all_gen_eff_den_{lab}"), "#it{R}_{L}^{tr}", os.path.join(updated_outdir, f"effpur/pair_gen_{lab}_matchedvsall.pdf"),
                            lines=INFO_LINES, extra_lines=[f"truth {lab} pairs"], logx=True, logy=True, matchvsall=True)
            PlotGiven1D(GetHist(f, f"pair_efficiency_{lab}_new"), os.path.join(updated_outdir, f"effpur/pair_efficiency_{lab}_new.pdf"), 
                        xtitle="#it{R}_{L}^{tr}", ytitle="Efficiency = matched gen/all gen", 
                        lines=INFO_LINES, extra_lines=[f"truth {lab} pairs"], ytext=0.3, 
                        logx=True, yscale01=True)
            draw_1d_overlay(GetHist(f, f"pair_match_rec_pur_num_{lab}"), GetHist(f, f"pair_all_rec_pur_den_{lab}"), "#it{R}_{L}^{det}", os.path.join(updated_outdir, f"effpur/pair_rec_{lab}_matchedvsall.pdf"),
                            lines=INFO_LINES, extra_lines=[f"detector {lab} pairs"], logx=True, logy=True, matchvsall=True)
            PlotGiven1D(GetHist(f, f"pair_purity_{lab}_new"), os.path.join(updated_outdir, f"effpur/pair_purity_{lab}_new.pdf"), 
                        xtitle="#it{R}_{L}^{det}", ytitle="Purity = matched rec/all rec", 
                        lines=INFO_LINES, extra_lines=[f"detector {lab} pairs"], ytext=0.3, 
                        logx=True, yscale01=True)

    


    # -----------------------------------------------------------------------
    # 6D THnSparse projections
    # -----------------------------------------------------------------------
    # Dictionary to collect histograms: {name: {lab: h2}}
    consolidated_h2s = {}
    
    for lab in EEC_LABELS:
        hs = f.Get(f"resp6_{grooming_str}_{lab}")
        if not hs:
            print(f"WARNING: resp6_{grooming_str}_{lab} not found")
            continue
        process_sparse(hs, lab, grooming_str, updated_outdir, lines=INFO_LINES, consolidated_dict=consolidated_h2s)
    
    # Call the consolidated plot function after looping through all labels
    plot_consolidated_2d(consolidated_h2s, EEC_LABELS, PROJECTIONS, grooming_str, updated_outdir, lines=INFO_LINES)

    f.Close()
    print(f"Wrote plots to {updated_outdir}/")


if __name__ == "__main__":
    main()