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
    ("weight", AX_W_DET,  AX_W_PART,  "weight",                      False, False, True),
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
                    logx=False, logy=False, lines=None, extra_lines=None, matchvsall=False):
    """Overlay det and part 1D projections and save."""
    c = _new_canvas(800, 650)
    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.06)

    if logx and _log_ok(h_det.GetXaxis()):
        c.SetLogx()
    if logy:
        c.SetLogy()

    # # normalize to unit area for shape comparison (guard against empty)
    # for h in (h_det, h_part):
    #     integral = h.Integral()
    #     if integral > 0:
    #         h.Scale(1.0 / integral)

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


def process_sparse(hs, lab, outdir, lines=None):
    """Make 2D det-vs-part and 1D det/part overlays for each observable."""
    for name, ax_det, ax_part, xtitle, logx1d, log2d, logy1d in PROJECTIONS:
        tag = [f"{lab}, {name} response"]

        # ---- 2D det vs part ----
        # THnSparse::Projection(y, x) -> TH2 with x = second arg, y = first arg
        h2 = hs.Projection(ax_part, ax_det)          # x = det, y = part
        h2.SetName(f"{lab}_{name}_2D")
        h2.SetTitle("")
        h2.GetXaxis().SetTitle(f"{xtitle} (det)")
        h2.GetYaxis().SetTitle(f"{xtitle} (part)")
        draw_th2(
            h2,
            os.path.join(outdir, f"resp2D_{lab}_{name}.pdf"),
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

def PlotGiven1D(hist, outpath, logx=False, logy=False):
    # hist = f.Get(hist_name)
    c = _new_canvas(800, 650)
    # c.SetLeftMargin(0.13)
    # c.SetRightMargin(0.05)
    # c.SetBottomMargin(0.12)
    # c.SetTopMargin(0.06)
    if logx and _log_ok(hist.GetXaxis()):
        c.SetLogx()
    if logy:
        c.SetLogy()


    hist.Draw("E1")
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


def PlotGiven2D(hist2D, outpath):
    c = _new_canvas(800, 650)
    # c.SetLeftMargin(0.13)
    # c.SetRightMargin(0.05)
    # c.SetBottomMargin(0.12)
    # c.SetTopMargin(0.06)

    hist2D.Draw("COLZ")
    c.SaveAs(outpath)
    c.Close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input",
                    # default="/rstorage/alice/AnalysisResults/blianggi/jse/rms/1836481/reponse_merged.root",
                    default="/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/57259276/response_merged_partial.root",
                    # default="/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/1836481/reponse_merged.root",
                    # default="/global/cfs/cdirs/alice/blianggi/mypyjetty/analysis/testing/response.root",
                    help="input ROOT file (from make_rms.py)")
    ap.add_argument("-o", "--outdir",
                    # default="/software/users/blianggi/mypyjetty/storage/jse/plots/response_matrices",
                    default="/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots/response_matrices",
                    help="output directory for plots")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    f = ROOT.TFile.Open(args.input)
    if not f or f.IsZombie():
        raise RuntimeError(f"could not open {args.input}")

    # -----------------------------------------------------------------------
    # 2D jet pt and groomed jet pt responses
    # -----------------------------------------------------------------------
    h_jetpt = f.Get("resp_jetpt")
    if h_jetpt:
        h_jetpt.SetTitle("")
        h_jetpt.GetXaxis().SetTitle("#it{p}_{T}^{det} (GeV/#it{c})")
        h_jetpt.GetYaxis().SetTitle("#it{p}_{T}^{part} (GeV/#it{c})")
        draw_th2(h_jetpt, os.path.join(args.outdir, "resp2D_jetpt.pdf"),
                 logx=False, logy=False, logz=True,
                 lines=INFO_LINES, extra_lines=["jet #it{p}_{T} response"])
    else:
        print("WARNING: resp_jetpt not found")

    h_groomed = f.Get("resp_groomed_jetpt")
    if h_groomed:
        h_groomed.SetTitle("")
        h_groomed.GetXaxis().SetTitle("#it{p}_{T,g}^{det} (GeV/#it{c})")
        h_groomed.GetYaxis().SetTitle("#it{p}_{T,g}^{part} (GeV/#it{c})")
        draw_th2(h_groomed, os.path.join(args.outdir, "resp2D_groomed_jetpt.pdf"),
                 logx=False, logy=False, logz=True,
                 lines=INFO_LINES, extra_lines=["groomed jet #it{p}_{T} response"])
    else:
        print("WARNING: resp_groomed_jetpt not found")

    # Plot efficiency/purity plots
    draw_1d_overlay(f.jet_match_gen_eff_num_ungroomed, f.jet_all_gen_eff_den_ungroomed, "#it{p}_{T} (GeV/#it{c})", os.path.join(args.outdir, "effpur/jet_gen_ungroomed_matchedvsall.pdf"),logx=True, logy=True)
    PlotGiven1D(f.jet_efficiency_ungroomed_new, os.path.join(args.outdir, "effpur/jet_efficiency_ungroomed_new.pdf"), logx=True)
    draw_1d_overlay(f.jet_match_rec_pur_num_ungroomed, f.jet_all_rec_pur_den_ungroomed, "#it{p}_{T} (GeV/#it{c})", os.path.join(args.outdir, "effpur/jet_rec_ungroomed_matchedvsall.pdf"),logx=True, logy=True)
    PlotGiven1D(f.jet_purity_ungroomed_new, os.path.join(args.outdir, "effpur/jet_purity_ungroomed_new.pdf"), logx=True)

    draw_1d_overlay(f.jet_match_gen_eff_num_groomed, f.jet_all_gen_eff_den_groomed, "#it{p}_{T,g} (GeV/#it{c})", os.path.join(args.outdir, "effpur/jet_gen_groomed_matchedvsall.pdf"),logx=True, logy=True)
    PlotGiven1D(f.jet_efficiency_groomed_new, os.path.join(args.outdir, "effpur/jet_efficiency_groomed_new.pdf"), logx=True)
    draw_1d_overlay(f.jet_match_rec_pur_num_groomed, f.jet_all_rec_pur_den_groomed, "#it{p}_{T} (GeV/#it{c})", os.path.join(args.outdir, "effpur/jet_rec_groomed_matchedvsall.pdf"),logx=True, logy=True)
    PlotGiven1D(f.jet_purity_groomed_new, os.path.join(args.outdir, "effpur/jet_purity_groomed_new.pdf"), logx=True)

    draw_1d_overlay(f.pair_match_gen_eff_num_AA, f.pair_all_gen_eff_den_AA, "#it{R}_{L}", os.path.join(args.outdir, "effpur/pair_gen_AA_matchedvsall.pdf"),logx=True, logy=True)
    PlotGiven1D(f.pair_efficiency_AA_new, os.path.join(args.outdir, "effpur/pair_efficiency_AA_new.pdf"))
    draw_1d_overlay(f.pair_match_rec_pur_num_AA, f.pair_all_rec_pur_den_AA, "#it{p}_{T}", os.path.join(args.outdir, "effpur/pair_rec_AA_matchedvsall.pdf"),logx=True, logy=True)
    PlotGiven1D(f.pair_purity_AA_new, os.path.join(args.outdir, "effpur/pair_purity_AA_new.pdf"))

    draw_1d_overlay(f.pair_match_gen_eff_num_AB, f.pair_all_gen_eff_den_AB, "#it{R}_{L}", os.path.join(args.outdir, "effpur/pair_gen_AB_matchedvsall.pdf"),logx=True, logy=True)
    PlotGiven1D(f.pair_efficiency_AB_new, os.path.join(args.outdir, "effpur/pair_efficiency_AB_new.pdf"))
    draw_1d_overlay(f.pair_match_rec_pur_num_AB, f.pair_all_rec_pur_den_AB, "#it{p}_{T}", os.path.join(args.outdir, "effpur/pair_rec_AB_matchedvsall.pdf"),logx=True, logy=True)
    PlotGiven1D(f.pair_purity_AB_new, os.path.join(args.outdir, "effpur/pair_purity_AB_new.pdf"))

    draw_1d_overlay(f.pair_match_gen_eff_num_BB, f.pair_all_gen_eff_den_BB, "#it{R}_{L}", os.path.join(args.outdir, "effpur/pair_gen_BB_matchedvsall.pdf"),logx=True, logy=True)
    PlotGiven1D(f.pair_efficiency_BB_new, os.path.join(args.outdir, "effpur/pair_efficiency_BB_new.pdf"))
    draw_1d_overlay(f.pair_match_rec_pur_num_BB, f.pair_all_rec_pur_den_BB, "#it{p}_{T}", os.path.join(args.outdir, "effpur/pair_rec_BB_matchedvsall.pdf"),logx=True, logy=True)
    PlotGiven1D(f.pair_purity_BB_new, os.path.join(args.outdir, "effpur/pair_purity_BB_new.pdf"))

    draw_1d_overlay(f.pair_match_gen_eff_num_rad, f.pair_all_gen_eff_den_rad, "#it{R}_{L}", os.path.join(args.outdir, "effpur/pair_gen_rad_matchedvsall.pdf"),logx=True, logy=True)
    PlotGiven1D(f.pair_efficiency_rad_new, os.path.join(args.outdir, "effpur/pair_efficiency_rad_new.pdf"))
    draw_1d_overlay(f.pair_match_rec_pur_num_rad, f.pair_all_rec_pur_den_rad, "#it{p}_{T}", os.path.join(args.outdir, "effpur/pair_rec_rad_matchedvsall.pdf"),logx=True, logy=True)
    PlotGiven1D(f.pair_purity_rad_new, os.path.join(args.outdir, "effpur/pair_purity_rad_new.pdf"))

    PlotGiven2D(f.lund_split_efficiency_new, os.path.join(args.outdir, "effpur/lund_split_efficiency_new.pdf"))
    PlotGiven2D(f.lund_split_purity_new, os.path.join(args.outdir, "effpur/lund_split_purity_new.pdf"))



    # -----------------------------------------------------------------------
    # 6D THnSparse projections
    # -----------------------------------------------------------------------
    for lab in EEC_LABELS:
        hs = f.Get(f"resp6_{lab}")
        if not hs:
            print(f"WARNING: resp6_{lab} not found")
            continue
        process_sparse(hs, lab, args.outdir, lines=INFO_LINES)

    f.Close()
    print(f"Wrote plots to {args.outdir}/")


if __name__ == "__main__":
    main()