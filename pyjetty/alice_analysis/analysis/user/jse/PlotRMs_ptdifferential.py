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

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)          # annotation replaces the pad title
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(255)
ROOT.TH1.AddDirectory(False)

EEC_LABELS = ["AA", "AB", "BB", "rad"]
PT_BINS = [10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 150.0, 200.0, 500.0]

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



def process_sparse(hs, lab, grooming_str, ptmin, ptmax, outdir, lines=None, weight_linbins=False):
    """Make 2D det-vs-part and 1D det/part overlays for each observable."""
    for name, ax_det, ax_part, xtitle, logx1d, log2d, logy1d in PROJECTIONS:
        if weight_linbins and name != "weight":
            continue

        grstr = "gr." if grooming_str == "groomed" else ""
        tag = [f"{lab}, {name} response ({grooming_str} jetpt bins)", f"{grstr} pt={ptmin}-{ptmax}"]

        # ---- 2D det vs part ----
        hs.GetAxis(AX_PT_DET).SetRangeUser(ptmin, ptmax)
        hs.GetAxis(AX_PT_PART).SetRangeUser(ptmin, ptmax)

        # THnSparse::Projection(y, x) -> TH2 with x = second arg, y = first arg
        h2 = hs.Projection(ax_part, ax_det)          # x = det, y = part
        h2.SetName(f"{lab}_{name}_2D_pt{ptmin}-{ptmax}")
        h2.SetTitle("")
        h2.GetXaxis().SetTitle(f"{xtitle} (det)")
        h2.GetYaxis().SetTitle(f"{xtitle} (part)")
        draw_th2(
            h2,
            os.path.join(outdir, f"resp2D_{grooming_str}_{name}_{lab}_pt{ptmin}-{ptmax}.pdf"),
            logx=log2d, logy=log2d, logz=True,
            lines=lines, extra_lines=tag,
        )

        # ---- 1D projections (det and part) ----
        h_det = hs.Projection(ax_det)
        h_det.SetName(f"{lab}_{name}_pt{ptmin}-{ptmax}_det")
        h_det.SetTitle("")
        h_part = hs.Projection(ax_part)
        h_part.SetName(f"{lab}_{name}_pt{ptmin}-{ptmax}_part")
        h_part.SetTitle("")

        draw_1d_overlay(
            h_det, h_part, xtitle,
            os.path.join(outdir, f"proj1D_{lab}_{name}_pt{ptmin}-{ptmax}.pdf"),
            logx=logx1d, logy=logy1d,
            lines=lines, extra_lines=[f"{lab}, {name}", f"{grstr} pt={ptmin}-{ptmax}"],
        )

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


def GetHist(f, histname):
    hist = f.Get(histname)
    return hist
    
# Only really need to run this script with groomed bins??
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input-resp",
                    default="/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/57911278/response_{}_merged_partial.root",
                    # default="/global/cfs/cdirs/alice/blianggi/mypyjetty/analysis/testing/response_ptdifferential.root",
                    help="input ROOT file (from make_rms.py)")
    ap.add_argument("-ipt", "--input-resp-ptdif",
                    default="/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/57911278/response_ptdifferential_merged.root",
                    # default="/global/cfs/cdirs/alice/blianggi/mypyjetty/analysis/testing/response_ptdifferential.root",
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
    grstr = "gr." if grooming_str == "groomed" else ""
    
    pt_outdir = os.path.join(args.outdir, f"{grooming_str}_bins/ptdifferential")
    os.makedirs(pt_outdir, exist_ok=True)
    os.makedirs(os.path.join(pt_outdir, "rms"), exist_ok=True)
    os.makedirs(os.path.join(pt_outdir, "splittings"), exist_ok=True)
    os.makedirs(os.path.join(pt_outdir, "pairs"), exist_ok=True)
    os.makedirs(os.path.join(pt_outdir, "rms/weight_linbins"), exist_ok=True)

    # fill the {} placeholder only if one is present, so an explicit -i still works
    infile = args.input_resp.format(grooming_str) if "{}" in args.input_resp else args.input_resp
    f = ROOT.TFile.Open(infile)
    if not f or f.IsZombie():
        raise RuntimeError(f"could not open {infile}")
    f_pt = ROOT.TFile.Open(args.input_resp_ptdif)
    if not f_pt or f_pt.IsZombie():
        raise RuntimeError(f"could not open {args.input_resp_ptdif}")

    # Alternative file: linear bins for energy weights, just for pt differential in weight distribution (not purity/eff)
    lin_weights_filepath = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/57259276/response_merged_partial.root" # not correct for high pt jets
    f_lin_weights = ROOT.TFile.Open(lin_weights_filepath)

    # FIRST make cuts on 3D matrices
    # -----------------------------------------------------------------------
    # 6D THnSparse projections
    # -----------------------------------------------------------------------
    for lab in EEC_LABELS:
        for ptbin in range(len(PT_BINS)-1):
            ptmin = PT_BINS[ptbin]
            ptmax = PT_BINS[ptbin+1]

            hs = f.Get(f"resp6_{grooming_str}_{lab}")
            if not hs:
                print(f"WARNING: resp6_{grooming_str}_{lab} not found")
                continue
            process_sparse(hs, lab, grooming_str, ptmin, ptmax, os.path.join(pt_outdir,"rms"), lines=INFO_LINES)

    f.Close()
    print(f"Wrote plots to {args.outdir}/")

    #1B: look at linear energy weights
    for lab in EEC_LABELS:
        for ptbin in range(len(PT_BINS)-1):
            ptmin = PT_BINS[ptbin]
            ptmax = PT_BINS[ptbin+1]

            hs = f_lin_weights.Get(f"resp6_{lab}")
            if not hs:
                print(f"WARNING: resp6_{lab} not found")
                continue
            process_sparse(hs, lab, grooming_str, ptmin, ptmax, os.path.join(pt_outdir,"rms/weight_linbins"), lines=INFO_LINES, weight_linbins=True)

    f_lin_weights.Close()


    
    # SECOND look at pt differential splitting/purity plots

    # Plot efficiency/purity plots
    for ptbin in range(len(PT_BINS)-1):
        ptmin = PT_BINS[ptbin]
        ptmax = PT_BINS[ptbin+1]
        
        PlotGiven2D(GetHist(f_pt, f"h_lund_matched_gen_pt{ptmin}-{ptmax}"), os.path.join(pt_outdir, f"splittings/lund_matched_gen_pt{ptmin}-{ptmax}.pdf"),
                    lines=INFO_LINES, extra_lines=[f"truth splittings", f"{grstr} pt={ptmin}-{ptmax}"])
        PlotGiven2D(GetHist(f_pt, f"h_lund_all_gen_pt{ptmin}-{ptmax}"), os.path.join(pt_outdir, f"splittings/lund_all_gen_pt{ptmin}-{ptmax}.pdf"),
                    lines=INFO_LINES, extra_lines=[f"truth splittings", f"{grstr} pt={ptmin}-{ptmax}"])
        PlotGiven2D(GetHist(f_pt, f"h_lund_matched_rec_pt{ptmin}-{ptmax}"), os.path.join(pt_outdir, f"splittings/lund_matched_rec_pt{ptmin}-{ptmax}.pdf"),
                    lines=INFO_LINES, extra_lines=[f"detector splittings", f"{grstr} pt={ptmin}-{ptmax}"])
        PlotGiven2D(GetHist(f_pt, f"h_lund_all_rec_pt{ptmin}-{ptmax}"), os.path.join(pt_outdir, f"splittings/lund_all_rec_pt{ptmin}-{ptmax}.pdf"),
                    lines=INFO_LINES, extra_lines=[f"detector splittings", f"{grstr} pt={ptmin}-{ptmax}"])

        PlotGiven2D(GetHist(f_pt, f"h_lund_split_efficiency_pt{ptmin}-{ptmax}_new"), os.path.join(pt_outdir, f"splittings/lund_split_efficiency_pt{ptmin}-{ptmax}_new.pdf"),
                    lines=INFO_LINES, extra_lines=[f"truth splittings", f"{grstr} pt={ptmin}-{ptmax}"], zscale01=True)
        PlotGiven2D(GetHist(f_pt, f"h_lund_split_purity_pt{ptmin}-{ptmax}_new"), os.path.join(pt_outdir, f"splittings/lund_split_purity_pt{ptmin}-{ptmax}_new.pdf"),
                    lines=INFO_LINES, extra_lines=[f"detector splittings", f"{grstr} pt={ptmin}-{ptmax}"], zscale01=True)
 
        
        for lab in EEC_LABELS:
            draw_1d_overlay(GetHist(f_pt, f"pair_match_gen_eff_num_pt{ptmin}-{ptmax}_{lab}"), GetHist(f_pt, f"pair_all_gen_eff_den_pt{ptmin}-{ptmax}_{lab}"), "#it{R}_{L}^{tr}", os.path.join(pt_outdir, f"pairs/pair_gen_{lab}_pt{ptmin}-{ptmax}_matchedvsall.pdf"),
                        lines=INFO_LINES, extra_lines=[f"truth {lab} pairs", f"{grstr} pt={ptmin}-{ptmax}"], logx=True, logy=True, matchvsall=True)
            PlotGiven1D(GetHist(f_pt, f"pair_efficiency_{lab}_pt{ptmin}-{ptmax}_new"), os.path.join(pt_outdir, f"pairs/pair_efficiency_{lab}_pt{ptmin}-{ptmax}_new.pdf"), 
                        xtitle="#it{R}_{L}^{tr}", ytitle="Efficiency = matched gen/all gen", 
                        lines=INFO_LINES, extra_lines=[f"truth {lab} pairs", f"{grstr} pt={ptmin}-{ptmax}"], ytext=0.4,
                        logx=True, yscale01=True)
            draw_1d_overlay(GetHist(f_pt, f"pair_match_rec_pur_num_pt{ptmin}-{ptmax}_{lab}"), GetHist(f_pt, f"pair_all_rec_pur_den_pt{ptmin}-{ptmax}_{lab}"), "#it{R}_{L}^{det}", os.path.join(pt_outdir, f"pairs/pair_rec_{lab}_pt{ptmin}-{ptmax}_matchedvsall.pdf"),
                        lines=INFO_LINES, extra_lines=[f"detector {lab} pairs", f"{grstr} pt={ptmin}-{ptmax}"], logx=True, logy=True, matchvsall=True)
            PlotGiven1D(GetHist(f_pt, f"pair_purity_{lab}_pt{ptmin}-{ptmax}_new"), os.path.join(pt_outdir, f"pairs/pair_purity_{lab}_pt{ptmin}-{ptmax}_new.pdf"), 
                        xtitle="#it{R}_{L}^{det}", ytitle="Purity = matched rec/all rec", 
                        lines=INFO_LINES, extra_lines=[f"detector {lab} pairs", f"{grstr} pt={ptmin}-{ptmax}"], ytext=0.4,
                        logx=True, yscale01=True)

    f_pt.Close()
    


if __name__ == "__main__":
    main()