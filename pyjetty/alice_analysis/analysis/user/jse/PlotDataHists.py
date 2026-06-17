# Plot DATA curves for JSE analysis (stage 3)
#
# Single data source, single AnalysisResults.root, jet-pt RANGE bins.

import os
import ROOT
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from enum import IntEnum


class Color(IntEnum):
    BLUE   = ROOT.TColor.GetColor("#1f77b4")
    GREEN  = ROOT.TColor.GetColor("#2ca02c")
    RED    = ROOT.TColor.GetColor("#d62728")
    ORANGE = ROOT.TColor.GetColor("#ff7f0e")


class PlotDataCurves:
    def __init__(self):
        ROOT.gROOT.SetBatch(True)
        ROOT.gStyle.SetLegendBorderSize(0)
        ROOT.gStyle.SetLegendFillColor(0)
        ROOT.gStyle.SetPadGridX(1)
        ROOT.gStyle.SetPadGridY(1)

        self.crosscheck = True

        # jet pT RANGE bins: 10-80 GeV in 7 bins of 10 GeV -> (lo, hi) tuples
        self.target_jet_pts = [(lo, lo + 10) for lo in range(10, 80, 10)]
        # -> [(10,20),(20,30),(30,40),(40,50),(50,60),(60,70),(70,80)]

        self.cut_modes = [("sd", 0.1), ("maxkt", None)]

        # *** single input file containing ALL slices ***
        # self.rootfile_path = ("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/rootfiles/AnalysisResultsCombined.root")
        self.rootfile_path = ("/global/cfs/cdirs/alice/blianggi/mypyjetty/analysis/testing/AnalysisResults.root")
        self.data_rootfile = None

        self.cut_mode = ""
        self.z_cut = None
        self.base_plot_dir = ("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots") #data

        self._persistent_canvases = []
        self._canvas_counter = 0

        # den_weight key -> token used in histogram names
        self.WEIGHT_TOKEN = {"jet": "wwjetpt", "rad": "wwradpt"}

    # -------------------------------------------------------------------------
    # pt-range helpers
    # -------------------------------------------------------------------------

    @staticmethod
    def _ptrange_token(jetpt):
        """jetpt is a (lo, hi) tuple -> 'jetpt70_80'."""
        lo, hi = jetpt
        return f"jetpt{lo}_{hi}"

    @staticmethod
    def _ptrange_label(jetpt):
        lo, hi = jetpt
        return f"{lo}-{hi}"

    @staticmethod
    def _ptrange_tag(jetpt):
        """For filenames: 'jetpt70_80'."""
        lo, hi = jetpt
        return f"jetpt{lo}_{hi}"

    # -------------------------------------------------------------------------
    # Name building
    # -------------------------------------------------------------------------

    def get_cut_suffix(self, z_cut):
        return f"sd{z_cut}" if self.cut_mode == "sd" else "maxkt"

    def get_cut_label(self, z_cut):
        return f"SD z_{{cut}} = {z_cut}" if self.cut_mode == "sd" else "maxkt selection"

    def get_passed_label(self):
        return "SD" if self.cut_mode == "sd" else "max k_{T}"

    def _ptRL(self, prefix, ptRL):
        """Return 'prefix_ptRL_' or 'prefix_' for name construction."""
        return f"{prefix}_ptRL_" if ptRL else f"{prefix}_"

    # -------------------------------------------------------------------------
    # File handling (single file, opened once)
    # -------------------------------------------------------------------------

    def open_rootfile(self):
        if self.data_rootfile and not self.data_rootfile.IsZombie():
            return
        self.data_rootfile = ROOT.TFile.Open(self.rootfile_path)
        if not self.data_rootfile or self.data_rootfile.IsZombie():
            raise RuntimeError(f"Could not open {self.rootfile_path}")

    # -------------------------------------------------------------------------
    # Output paths and canvases
    # -------------------------------------------------------------------------

    def get_output_dir(self, subdir, z_cut, plot_type=None):
        parts = [self.base_plot_dir, subdir, self.get_cut_suffix(z_cut)]
        if plot_type:
            parts.extend(plot_type.split('/'))
        path = os.path.join(*parts)
        os.makedirs(path, exist_ok=True)
        return path

    def make_canvas(self, base_name, tag="", z_cut="", den_weight="",
                    title=None, w=800, h=600):
        self._canvas_counter += 1
        unique_name = (
            f"{base_name}_{tag}_{self.get_cut_suffix(z_cut)}"
            f"_{den_weight}_{self._canvas_counter}"
        )
        canvas = ROOT.TCanvas(unique_name, title or base_name, w, h)
        self._persistent_canvases.append(canvas)
        return canvas

    def _save_canvas(self, canvas, subdir, z_cut, filename, plot_type=None):
        canvas.SaveAs(os.path.join(self.get_output_dir(subdir, z_cut, plot_type), filename))

    # -------------------------------------------------------------------------
    # Formatting helpers
    # -------------------------------------------------------------------------

    def FormatHist(self, hist, color, linestyle, markerstyle=0, coloralpha=1):
        if hist is None:
            return
        hist.SetLineColor(color)
        if coloralpha != 1:
            hist.SetLineColorAlpha(color, coloralpha)
        hist.SetLineStyle(linestyle)
        hist.SetLineWidth(2)
        if markerstyle > 0:
            hist.SetLineStyle(1)
            hist.SetMarkerStyle(markerstyle)
            hist.SetMarkerColor(color)

    def draw_hori_line(self, x1, x2, y1, color, linestyle, linewidth=1):
        line = ROOT.TLine(x1, y1, x2, y1)
        line.SetLineWidth(linewidth)
        line.SetLineColor(color)
        line.SetLineStyle(linestyle)
        line.Draw("SAME")
        return line

    def MakeEventLeg(self, jetpt_label, z_cut, den_weight="",
                     x1=0.15, y1=0.7, x2=0.40, y2=0.88):
        leg = ROOT.TLegend(x1, y1, x2, y2)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetMargin(0)
        leg.AddEntry(ROOT.nullptr, "pp data, R = 0.4 jets", "")
        leg.AddEntry(ROOT.nullptr, f"jet p_{{T}} = {jetpt_label} GeV/c", "")
        leg.AddEntry(ROOT.nullptr, self.get_cut_label(z_cut), "")
        if den_weight:
            leg.AddEntry(ROOT.nullptr,
                         f"using weight p_{{T,1}}p_{{T,2}} / p_{{T,{den_weight}}}^{{2}}", "")
        return leg

    # -------------------------------------------------------------------------
    # Histogram retrieval
    # -------------------------------------------------------------------------

    def _get(self, name):
        h = self.data_rootfile.Get(name)
        if h:
            h.SetDirectory(0)
            return h
        return None

    def GetSubjetEECHists(self, jetpt, z_cut, den_weight, ptRL=False):
        cut = self.get_cut_suffix(z_cut)
        pr = self._ptrange_token(jetpt)
        w = self.WEIGHT_TOKEN[den_weight]

        # hist_full has NO weight token; has a ptRL variant
        full_name = f"{self._ptRL('hist_full', ptRL)}{pr}_{cut}"
        rad_name  = f"{self._ptRL('hist_rad', ptRL)}{pr}_{cut}_{w}"
        aa_name   = f"{self._ptRL('hist_AA', ptRL)}{pr}_{cut}_{w}"
        bb_name   = f"{self._ptRL('hist_BB', ptRL)}{pr}_{cut}_{w}"
        ab_name   = f"{self._ptRL('hist_AB', ptRL)}{pr}_{cut}_{w}"

        return (self._get(full_name), self._get(rad_name), self._get(aa_name),
                self._get(bb_name), self._get(ab_name))

    def GetCABHist(self, jetpt, z_cut, den_weight, ptRL=False):
        cut = self.get_cut_suffix(z_cut)
        pr = self._ptrange_token(jetpt)
        w = self.WEIGHT_TOKEN[den_weight]
        name = f"{self._ptRL('CAB', ptRL)}{pr}_{cut}_{w}"
        return self._get(name)

    def GetRGHist(self, jetpt, z_cut):
        cut = self.get_cut_suffix(z_cut)
        pr = self._ptrange_token(jetpt)
        return self._get(f"hist_rg_{pr}_{cut}")

    def GetRadPropHist(self, radprop, jetpt, z_cut):
        cut = self.get_cut_suffix(z_cut)
        pr = self._ptrange_token(jetpt)
        return self._get(f"radiator_{radprop}_{pr}_{cut}")

    def GetPassedJetPtHist(self, jetpt, z_cut):
        """hist_jetpt_jetpt..._cut  (post-cut). 'all' variant is pre-cut."""
        cut = self.get_cut_suffix(z_cut)
        pr = self._ptrange_token(jetpt)
        return self._get(f"hist_jetpt_{pr}_{cut}")

    # -------------------------------------------------------------------------
    # Drawing utilities
    # -------------------------------------------------------------------------

    def _draw_eec_set(self, hists, first=True):
        draw_opt = "HIST" if first else "HIST SAME"
        for h in hists:
            h.Draw(draw_opt)
            draw_opt = "HIST SAME"

    def _draw_eec_canvas(self, canvas, hists_primary, ev_leg=None, legend=None,
                         crosscheck_hist=None):
        canvas.cd()
        all_hists = list(hists_primary)
        if crosscheck_hist:
            all_hists.append(crosscheck_hist)
        global_max = max(h.GetMaximum() for h in all_hists if h)
        hists_primary[0].GetYaxis().SetRangeUser(0, global_max * 1.4)

        self._draw_eec_set(hists_primary, first=True)
        if crosscheck_hist:
            crosscheck_hist.Draw("HIST SAME")
        for leg in (ev_leg, legend):
            if leg:
                leg.Draw()

    # -------------------------------------------------------------------------
    # Plot: basic subjet EEC
    # -------------------------------------------------------------------------

    def plot_basic(self, jetpt, z_cut, den_weight):
        cut = self.get_cut_suffix(z_cut)
        tag = self._ptrange_tag(jetpt)
        label = self._ptrange_label(jetpt)

        for ptRL, subdir_suffix, canvas_name in [
            (False, "subjet_eec",      "canvas_basic"),
            (True,  "subjet_eec/ptRL", "canvas_ptRL"),
        ]:
            hists = self.GetSubjetEECHists(jetpt, z_cut, den_weight, ptRL=ptRL)
            hist_full, hist_rad, hist_AA, hist_BB, hist_AB = hists
            if not all([hist_rad, hist_AA, hist_BB, hist_AB]):
                continue

            canvas = self.make_canvas(canvas_name, tag=tag,
                                      z_cut=z_cut, den_weight=den_weight)
            canvas.SetLogx()

            self.FormatHist(hist_full, ROOT.kGray + 1, ROOT.kSolid)
            self.FormatHist(hist_rad, ROOT.kBlack,     ROOT.kSolid)
            self.FormatHist(hist_AA,  Color.BLUE,      ROOT.kSolid)
            self.FormatHist(hist_BB,  Color.ORANGE,    ROOT.kSolid)
            self.FormatHist(hist_AB,  Color.GREEN,     ROOT.kSolid)

            crosscheck = None
            if not ptRL and self.crosscheck:
                crosscheck = hist_AA.Clone(
                    f"hist_crosscheck_data_{tag}_{cut}_{den_weight}")
                crosscheck.SetDirectory(0)
                crosscheck.Add(hist_BB)
                crosscheck.Add(hist_AB)
                self.FormatHist(crosscheck, ROOT.kGray + 2, ROOT.kDashed)

            ev_leg = self.MakeEventLeg(label, z_cut, den_weight)
            legend = ROOT.TLegend(0.68, 0.5, 0.88, 0.68)
            if hist_full is not None:
                legend.AddEntry(hist_full,
                                f"all jets that passed {self.get_passed_label()}", "l")
            legend.AddEntry(hist_rad, "radiator", "l")
            if crosscheck:
                legend.AddEntry(crosscheck, "cross-check #Sigma_{subjets} = rad.", "l")
            legend.AddEntry(hist_AA, "AxA", "l")
            legend.AddEntry(hist_BB, "BxB", "l")
            legend.AddEntry(hist_AB, "AxB", "l")

            hists_to_draw = [h for h in hists if h is not None]
            self._draw_eec_canvas(canvas, hists_to_draw, ev_leg=ev_leg,
                                  legend=legend, crosscheck_hist=crosscheck)

            prefix = "subjet_eec_ptRL" if ptRL else "subjet_eec"
            output_name = (f"{prefix}_comparison_data_{tag}"
                           f"_R0.4_{cut}_ww{den_weight}.pdf")
            self._save_canvas(canvas, "data", z_cut, output_name,
                              plot_type=subdir_suffix)

    # -------------------------------------------------------------------------
    # Plot: R_g overlaid with AA, BB, AB (SD only)
    # -------------------------------------------------------------------------

    def plot_rg(self, jetpt, z_cut, den_weight):
        if self.cut_mode != "sd":
            return
        cut = self.get_cut_suffix(z_cut)
        tag = self._ptrange_tag(jetpt)
        label = self._ptrange_label(jetpt)

        hist_rg = self.GetRGHist(jetpt, z_cut)
        if not hist_rg:
            print(f"WARNING: R_g histogram not found for data {tag} {cut}")
            return

        _, _, hist_AA, hist_BB, hist_AB = self.GetSubjetEECHists(
            jetpt, z_cut, den_weight)
        if not all([hist_AA, hist_BB, hist_AB]):
            print(f"WARNING: AA/BB/AB not found for data {tag} {cut}")
            return

        # ----- Canvas 1: R_g vs all EEC components -----
        canvas = self.make_canvas("canvas_rg", tag=tag, z_cut=z_cut,
                                  den_weight=den_weight,
                                  title="R_{g} vs EEC components")
        canvas.SetLogx()
        canvas.cd()

        hist_AA.SetLineColor(Color.BLUE)
        hist_BB.SetLineColor(Color.ORANGE)
        hist_AB.SetLineColor(Color.GREEN)
        hist_rg.SetLineColor(ROOT.kRed + 1)
        hist_rg.SetLineStyle(ROOT.kDashed)

        y_max = max(hist_AA.GetMaximum(), hist_BB.GetMaximum(),
                    hist_AB.GetMaximum(), hist_rg.GetMaximum()) * 1.3
        hist_AA.SetMaximum(y_max)
        hist_AA.SetMinimum(0)
        hist_AA.GetXaxis().SetTitle("R_{L} or R_{g}")
        hist_AA.GetYaxis().SetTitle("(1/N_{jets}) dN/d(R_{L} or R_{g})")

        hist_AA.Draw("HIST")
        hist_BB.Draw("HIST SAME")
        hist_AB.Draw("HIST SAME")
        hist_rg.Draw("HIST SAME")

        ev_leg = self.MakeEventLeg(label, z_cut, den_weight)
        legend = ROOT.TLegend(0.55, 0.65, 0.88, 0.88)
        legend.AddEntry(hist_AA, "AxA (EEC)", "l")
        legend.AddEntry(hist_BB, "BxB (EEC)", "l")
        legend.AddEntry(hist_AB, "AxB (EEC)", "l")
        legend.AddEntry(hist_rg, "R_{g} = #DeltaR_{AB}/R", "l")
        ev_leg.Draw()
        legend.Draw()

        output_name = (f"rg_vs_eec_data_{tag}_R0.4_{cut}_ww{den_weight}.pdf")
        self._save_canvas(canvas, "data", z_cut, output_name, plot_type="rg")

        # ----- Canvas 2: R_g vs AxB self-normalized with ratio panel -----
        hist_rg_norm = hist_rg.Clone(f"hist_rg_norm_data_{tag}")
        hist_AB_norm = hist_AB.Clone(f"hist_AB_norm_data_{tag}")
        hist_rg_norm.SetDirectory(0)
        hist_AB_norm.SetDirectory(0)
        if hist_rg_norm.Integral() > 0:
            hist_rg_norm.Scale(1.0 / hist_rg_norm.Integral())
        if hist_AB_norm.Integral() > 0:
            hist_AB_norm.Scale(1.0 / hist_AB_norm.Integral())

        canvas2 = self.make_canvas("canvas_rg_vs_AxB", tag=tag, z_cut=z_cut,
                                   den_weight=den_weight,
                                   title="R_{g} vs AxB (self-normalized)")

        pad1 = ROOT.TPad("pad1_rgAxB", "pad1_rgAxB", 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0.02)
        pad1.SetLogx()
        pad1.Draw()
        canvas2.cd()
        pad2 = ROOT.TPad("pad2_rgAxB", "pad2_rgAxB", 0, 0.0, 1, 0.3)
        pad2.SetTopMargin(0.02)
        pad2.SetBottomMargin(0.35)
        pad2.SetLogx()
        pad2.SetGridy()
        pad2.Draw()

        pad1.cd()
        hist_rg_norm.SetLineColor(ROOT.kRed + 1)
        hist_rg_norm.SetLineStyle(ROOT.kDashed)
        hist_AB_norm.SetLineColor(Color.GREEN)
        hist_AB_norm.SetLineStyle(ROOT.kSolid)
        y_max2 = max(hist_rg_norm.GetMaximum(), hist_AB_norm.GetMaximum()) * 1.3
        hist_AB_norm.SetMaximum(y_max2)
        hist_AB_norm.SetMinimum(0)
        hist_AB_norm.GetYaxis().SetTitle("Self-normalized")
        hist_AB_norm.GetYaxis().SetTitleSize(0.05)
        hist_AB_norm.GetYaxis().SetTitleOffset(0.9)
        hist_AB_norm.GetXaxis().SetLabelSize(0)
        hist_AB_norm.Draw("HIST")
        hist_rg_norm.Draw("HIST SAME")

        ev_leg2 = self.MakeEventLeg(label, z_cut, den_weight)
        legend2 = ROOT.TLegend(0.65, 0.74, 0.88, 0.88)
        legend2.AddEntry(hist_AB_norm, "AxB (EEC)", "l")
        legend2.AddEntry(hist_rg_norm, "R_{g} = #DeltaR_{AB}/R", "l")
        ev_leg2.Draw()
        legend2.Draw()

        pad2.cd()
        hist_ratio = hist_AB_norm.Clone(f"hist_ratio_AxB_rg_data_{tag}")
        hist_ratio.SetDirectory(0)
        hist_ratio.Divide(hist_rg_norm)
        hist_ratio.SetLineColor(ROOT.kBlack)
        hist_ratio.SetMarkerStyle(20)
        hist_ratio.SetMarkerSize(0.7)
        hist_ratio.GetYaxis().SetTitle("AxB / R_{g}")
        hist_ratio.GetYaxis().SetNdivisions(505)
        hist_ratio.GetYaxis().SetTitleSize(0.11)
        hist_ratio.GetYaxis().SetTitleOffset(0.4)
        hist_ratio.GetYaxis().SetLabelSize(0.09)
        hist_ratio.GetXaxis().SetTitle("R_{L} or R_{g}")
        hist_ratio.GetXaxis().SetTitleSize(0.12)
        hist_ratio.GetXaxis().SetTitleOffset(1.0)
        hist_ratio.GetXaxis().SetLabelSize(0.09)
        hist_ratio.SetMinimum(0.0)
        hist_ratio.SetMaximum(2.0)
        hist_ratio.Draw("EP")
        # corrected argument order: (x1, x2, y, color, linestyle)
        self.draw_hori_line(hist_ratio.GetXaxis().GetXmin(),
                            hist_ratio.GetXaxis().GetXmax(),
                            1.0, ROOT.kGray + 2, ROOT.kDashed)

        canvas2.cd()
        output_name2 = (f"rg_vs_AxB_data_{tag}_R0.4_{cut}_ww{den_weight}.pdf")
        self._save_canvas(canvas2, "data", z_cut, output_name2, plot_type="rg")

    # -------------------------------------------------------------------------
    # Plot: C_AB (regular + ptRL)
    # -------------------------------------------------------------------------

    def _plot_CAB_impl(self, jetpt, z_cut, ptRL=False):
        cut = self.get_cut_suffix(z_cut)
        tag = self._ptrange_tag(jetpt)
        label = self._ptrange_label(jetpt)

        canvas_name = "can_CAB_ptRL" if ptRL else "can_CAB"
        title       = "C_{AB}_ptRL"  if ptRL else "C_{AB}"
        subdir      = "CAB/ptRL"     if ptRL else "CAB"
        file_prefix = "CAB_ptRL"     if ptRL else "CAB"

        can_CAB = self.make_canvas(canvas_name, tag=tag, z_cut=z_cut, title=title)
        can_CAB.SetLogx()
        can_CAB.SetLogy()

        hist_wwjetpt = self.GetCABHist(jetpt, z_cut, "jet", ptRL=ptRL)
        hist_wwradpt = self.GetCABHist(jetpt, z_cut, "rad", ptRL=ptRL)
        if not (hist_wwjetpt and hist_wwradpt):
            print(f"WARNING: CAB hist(s) not found for data {tag} {cut} (ptRL={ptRL})")
            return

        self.FormatHist(hist_wwjetpt, Color.BLUE,   ROOT.kSolid)
        self.FormatHist(hist_wwradpt, Color.ORANGE, ROOT.kSolid)

        ev_leg = self.MakeEventLeg(label, z_cut)
        legend = ROOT.TLegend(0.2, 0.5, 0.4, 0.65)
        legend.AddEntry(hist_wwjetpt, "weight uses jet pt", "l")
        legend.AddEntry(hist_wwradpt, "weight uses rad pt", "l")

        can_CAB.cd()
        hist_wwjetpt.Draw("HIST")
        hist_wwradpt.Draw("HIST SAME")
        ev_leg.Draw()
        legend.Draw()
        self.draw_hori_line(1e-3, 1, 1, ROOT.kGray + 3, 9)

        output_name = (f"{file_prefix}_data_{tag}_R0.4_{cut}.pdf")
        self._save_canvas(can_CAB, "data", z_cut, output_name, plot_type=subdir)

    def plot_CAB(self, jetpt, z_cut):
        self._plot_CAB_impl(jetpt, z_cut, ptRL=False)

    def plot_CAB_ptRL(self, jetpt, z_cut):
        self._plot_CAB_impl(jetpt, z_cut, ptRL=True)

    # -------------------------------------------------------------------------
    # Plot: radiator property
    # -------------------------------------------------------------------------

    def plot_rad_prop(self, radprop, jetpt, z_cut):
        cut = self.get_cut_suffix(z_cut)
        tag = self._ptrange_tag(jetpt)
        label = self._ptrange_label(jetpt)

        hist = self.GetRadPropHist(radprop, jetpt, z_cut)
        if not hist:
            print(f"WARNING: radiator_{radprop} hist not found for data {tag} {cut}")
            return

        can = self.make_canvas("can_radprop", tag=tag, z_cut=z_cut, den_weight=radprop)
        self.FormatHist(hist, ROOT.kBlack, ROOT.kSolid)
        can.cd()
        hist.Draw("HIST")
        ev_leg = self.MakeEventLeg(label, z_cut)
        ev_leg.Draw()
        output_name = (f"radiator_{radprop}_data_{tag}_R0.4_{cut}.pdf")
        self._save_canvas(can, "data", z_cut, output_name, "radiator_prop")

    # -------------------------------------------------------------------------
    # Plot: EEC components overlaid across all jet pT
    # -------------------------------------------------------------------------

    def plot_basic_acrosspt(self, z_cut, den_weight):
        cut = self.get_cut_suffix(z_cut)
        n_pt = len(self.target_jet_pts)
        # build colour/marker cycles long enough for 7 bins
        base_colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kOrange + 7, ROOT.kGreen + 2,
                       ROOT.kRed + 1, ROOT.kMagenta + 1, ROOT.kCyan + 2]
        base_markers = [ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kFullDiamond,
                        ROOT.kFullStar, ROOT.kFullTriangleUp, ROOT.kFullTriangleDown,
                        ROOT.kFullCross]
        colors  = [base_colors[i % len(base_colors)] for i in range(n_pt)]
        markers = [base_markers[i % len(base_markers)] for i in range(n_pt)]

        panel_labels = ["radiator", "AxA", "BxB", "AxB"]

        canvas = self.make_canvas("canvas_acrosspt", z_cut=z_cut,
                                  den_weight=den_weight, w=1200, h=1000)
        canvas.Divide(2, 2, 0.005, 0.005)

        ev_leg = self.MakeEventLeg("10-80", z_cut, den_weight)
        legend_pt = ROOT.TLegend(0.4, 0.6, 0.77, 0.88)
        legend_pt.SetBorderSize(0)

        dummy_graphs = []
        for jetpt, marker, color in zip(self.target_jet_pts, markers, colors):
            d = ROOT.TGraph(1)
            d.SetMarkerStyle(marker)
            d.SetMarkerColor(color)
            d.SetLineColor(color)
            legend_pt.AddEntry(d, f"jet p_{{T}} = {self._ptrange_label(jetpt)} GeV/c", "p")
            dummy_graphs.append(d)

        # persistent[ijetpt] = [rad, AA, BB, AB] or None
        persistent = []
        for ijetpt, jetpt in enumerate(self.target_jet_pts):
            hists = self.GetSubjetEECHists(jetpt, z_cut, den_weight)
            _, hr, hAA, hBB, hAB = hists
            if not all([hr, hAA, hBB, hAB]):
                persistent.append(None)
                continue
            tag = self._ptrange_tag(jetpt)
            comps = []
            for h, name in zip([hr, hAA, hBB, hAB], ["rad", "AA", "BB", "AB"]):
                c = h.Clone(f"acrosspt_{name}_data_{tag}_{cut}_ww{den_weight}")
                c.SetDirectory(0)
                if c.Integral() > 0:
                    c.Scale(1.0 / c.Integral())
                self.FormatHist(c, colors[ijetpt], ijetpt + 1, markers[ijetpt])
                comps.append(c)
            persistent.append(comps)

        latex_labels = []
        for ipanel, label in enumerate(panel_labels):
            pad = canvas.cd(ipanel + 1)
            pad.SetLogx()
            pad.SetLeftMargin(0.14)
            pad.SetBottomMargin(0.14)

            ymax = 0.0
            for ijetpt in range(n_pt):
                if persistent[ijetpt] is None:
                    continue
                ymax = max(ymax, persistent[ijetpt][ipanel].GetMaximum())

            first = True
            for ijetpt in range(n_pt):
                if persistent[ijetpt] is None:
                    continue
                h = persistent[ijetpt][ipanel]
                if first:
                    h.SetMaximum(ymax * 1.3)
                    h.SetMinimum(0)
                    h.GetYaxis().SetTitle("self-normalized")
                    h.GetYaxis().SetTitleSize(0.06)
                    h.GetYaxis().SetTitleOffset(0.9)
                    h.GetXaxis().SetTitleSize(0.06)
                    h.GetXaxis().SetTitle("R_{L}")
                h.Draw("PE" if first else "PE SAME")
                first = False

            latex = ROOT.TLatex()
            latex.SetNDC()
            latex.SetTextSize(0.07)
            latex.DrawLatex(0.18, 0.84, label)
            latex_labels.append(latex)

            if ipanel == 0:
                ev_leg.Draw()
                legend_pt.Draw()

        output_name = (f"subjet_eec_acrosspt_data_alljetpt_R0.4_{cut}_ww{den_weight}.pdf")
        self._save_canvas(canvas, "data", z_cut, output_name, "subjet_eec")

        # ---- C_AB overlaid across pT ----
        can_CAB = self.make_canvas("can_CAB_acrosspt", z_cut=z_cut,
                                   den_weight=den_weight, title="C_{AB}")
        can_CAB.SetLogx()
        can_CAB.SetLogy()
        can_CAB.cd()
        legend_CAB = ROOT.TLegend(0.38, 0.12, 0.58, 0.32)
        legend_CAB.SetBorderSize(0)
        persistent_CAB = []
        first_CAB = True
        for ijetpt, jetpt in enumerate(self.target_jet_pts):
            h = self.GetCABHist(jetpt, z_cut, den_weight)
            if not h:
                continue
            tag = self._ptrange_tag(jetpt)
            c = h.Clone(f"acrosspt_CAB_data_{tag}_{cut}_ww{den_weight}")
            c.SetDirectory(0)
            self.FormatHist(c, colors[ijetpt], ijetpt + 1, markers[ijetpt])
            legend_CAB.AddEntry(c, f"jet p_{{T}} = {self._ptrange_label(jetpt)} GeV/c", "pe")
            persistent_CAB.append(c)
            c.Draw("PE" if first_CAB else "PE SAME")
            first_CAB = False

        if persistent_CAB:
            ev_leg.Draw()
            legend_CAB.Draw()
            self.draw_hori_line(1e-3, 1, 1, ROOT.kGray + 3, 9)

        output_name = (f"CAB_acrosspt_data_alljetpt_R0.4_{cut}_ww{den_weight}.pdf")
        self._save_canvas(can_CAB, "data", z_cut, output_name, plot_type="CAB")

    # -------------------------------------------------------------------------
    # MPV (Most Probable Value) machinery
    # -------------------------------------------------------------------------

    @staticmethod
    def gaus_log(x, mu, C, sg):
        return C * np.exp(-(np.log(x / mu)) ** 2 / (2 * sg * sg))

    @staticmethod
    def _select_fit_window(xs, ys, yerrs, n_side=3, n_min=2):
        n = len(ys)
        if n == 0:
            return xs, ys, yerrs, None
        peak = int(np.argmax(ys))
        left_avail = peak
        right_avail = n - 1 - peak
        left = min(n_side, left_avail)
        right = min(n_side, right_avail)
        target_total = 2 * n_side + 1
        deficit = target_total - (left + right + 1)
        if deficit > 0:
            extra_left = min(deficit, left_avail - left)
            left += extra_left
            deficit -= extra_left
            extra_right = min(deficit, right_avail - right)
            right += extra_right
        lo = peak - left
        hi = peak + right + 1
        return xs[lo:hi], ys[lo:hi], yerrs[lo:hi], peak

    def _fit_mpv(self, hist, fit_range=None, save_path=None, plot_title=None):
        if hist is None:
            return None, None
        nb = hist.GetNbinsX()
        xs, ys, yerrs = [], [], []
        for ib in range(1, nb + 1):
            x = hist.GetBinCenter(ib)
            y = hist.GetBinContent(ib)
            e = hist.GetBinError(ib)
            if x <= 0 or y <= 0:
                continue
            if fit_range is not None and (x < fit_range[0] or x > fit_range[1]):
                continue
            xs.append(x); ys.append(y); yerrs.append(e if e > 0 else 1.0)

        if len(xs) < 5:
            return None, None
        xs = np.array(xs); ys = np.array(ys); yerrs = np.array(yerrs)
        xs_full, ys_full, yerrs_full = xs, ys, yerrs
        xs, ys, yerrs, _ = self._select_fit_window(xs, ys, yerrs, n_side=3, n_min=2)

        if len(xs) < 3:
            print(f"  Not enough points around peak to fit ({len(xs)} found)")
            if save_path is not None:
                try:
                    self._save_fit_diagnostic(xs_full, ys_full, yerrs_full, None,
                                              save_path, plot_title, None, None, p0=None)
                except Exception:
                    pass
            return None, None

        mu0 = xs[np.argmax(ys)]
        C0  = ys.max()
        half_max = ys.max() / 2.0
        above = xs[ys >= half_max]
        if len(above) >= 2:
            sg0 = (np.log(above.max()) - np.log(above.min())) / 2.355
            sg0 = max(sg0, 0.05)
        else:
            sg0 = 0.5
        p0 = [mu0, C0, sg0]

        mu, mu_err, popt = None, None, None
        try:
            popt, pcov = curve_fit(self.gaus_log, xs, ys, p0=p0,
                                   sigma=yerrs, absolute_sigma=False, maxfev=5000)
            mu_fit, _, _ = popt
            if not np.isfinite(mu_fit) or mu_fit <= 0:
                popt = None
            else:
                mu = float(mu_fit)
                mu_err = float(np.sqrt(pcov[0, 0])) if pcov is not None else 0.0
        except Exception as ex:
            print(f"  MPV fit failed: {ex}")
            popt = None

        if save_path is not None:
            try:
                self._save_fit_diagnostic(xs_full, ys_full, yerrs_full, popt,
                                          save_path, plot_title, mu, mu_err,
                                          p0=p0, fit_xs=xs)
            except Exception as ex:
                print(f"  Failed to save fit diagnostic: {ex}")
        return mu, mu_err

    def _save_fit_diagnostic(self, xs, ys, yerrs, popt, save_path, title,
                             mu, mu_err, p0=None, fit_xs=None):
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.errorbar(xs, ys, yerr=yerrs, fmt="o", markersize=4, color="black",
                    capsize=2, label="data", zorder=2)
        if fit_xs is not None and len(fit_xs):
            mask = np.isin(xs, fit_xs)
            ax.errorbar(xs[mask], ys[mask], yerr=yerrs[mask], fmt="o", markersize=6,
                        mfc="none", mec="C2", mew=1.5, capsize=2,
                        label="fit window", zorder=2.2)
            grid_lo, grid_hi = fit_xs.min(), fit_xs.max()
        else:
            grid_lo, grid_hi = xs.min(), xs.max()
        x_grid = np.logspace(np.log10(grid_lo), np.log10(grid_hi), 400)
        if p0 is not None:
            ax.plot(x_grid, self.gaus_log(x_grid, *p0), "--", color="C0",
                    linewidth=1.5, alpha=0.8,
                    label=(f"initial guess\n$\\mu_0$={p0[0]:.4g}, "
                           f"$C_0$={p0[1]:.4g}, $\\sigma_0$={p0[2]:.4g}"), zorder=2.5)
        if popt is not None:
            ax.plot(x_grid, self.gaus_log(x_grid, *popt), "-", color="C3", linewidth=2,
                    label=(f"gaus_log fit\n$\\mu$={popt[0]:.4g}, "
                           f"$C$={popt[1]:.4g}, $\\sigma$={popt[2]:.4g}"), zorder=3)
            ax.axvline(popt[0], color="C3", linestyle="--", alpha=0.6,
                       label=f"peak $\\mu$={popt[0]:.4g} $\\pm$ {mu_err:.2g}")
        else:
            ax.text(0.5, 0.5, "FIT FAILED", transform=ax.transAxes, ha="center",
                    va="center", fontsize=20, color="red", alpha=0.5)
        if p0 is not None:
            ax.axvline(p0[0], color="C0", linestyle=":", alpha=0.5,
                       label=f"seed $\\mu_0$={p0[0]:.4g}")
        ax.set_xscale("log")
        ax.set_xlabel("x")
        ax.set_ylabel("entries")
        if title:
            ax.set_title(title, fontsize=10)
        ax.legend(fontsize=8, loc="best")
        ax.grid(True, which="both", alpha=0.3)
        fig.tight_layout()
        fig.savefig(save_path)
        plt.close(fig)

    def _init_mpv_storage(self):
        # self.mpv_data[case_key][component][jetpt_mid] = (mu, mu_err)
        self.mpv_data = {}

    def _fit_diag_path(self, case_key, component, tag):
        return os.path.join(self.base_plot_dir, "mpv_summary", "fits", case_key,
                            f"fit_{case_key}_{component}_{tag}.pdf")

    def _store_mpv(self, case_key, component, jetpt_mid, mu, mu_err):
        if mu is None:
            return
        self.mpv_data.setdefault(case_key, {}).setdefault(component, {})[jetpt_mid] = (mu, mu_err)

    def _collect_mpv_basic(self, jetpt, z_cut, den_weight):
        cut = self.get_cut_suffix(z_cut)
        tag = self._ptrange_tag(jetpt)
        lo, hi = jetpt
        jetpt_mid = 0.5 * (lo + hi)   # x-position for the summary plot
        for ptRL, t in [(False, "basic"), (True, "basic_ptRL")]:
            hists = self.GetSubjetEECHists(jetpt, z_cut, den_weight, ptRL=ptRL)
            case_key = f"{t}_data_{cut}_ww{den_weight}"
            for comp_name, h in zip(["full", "rad", "AA", "BB", "AB"], hists):
                if h is None:
                    continue
                save_path = self._fit_diag_path(
                    f"phys_data_{cut}_ww{den_weight}" + ("_ptRL" if ptRL else ""),
                    comp_name, tag)
                title = (f"{t} | data | jet p_T={self._ptrange_label(jetpt)} | "
                         f"{cut} | ww{den_weight} | {comp_name}")
                mu, mu_err = self._fit_mpv(h, save_path=save_path, plot_title=title)
                self._store_mpv(case_key, comp_name, jetpt_mid, mu, mu_err)

    def _plot_mpv_summary(self, case_key, title, xlabel, ylabel, outdir, filename,
                          component_styles=None):
        if case_key not in self.mpv_data:
            return
        data = self.mpv_data[case_key]
        if not data:
            return
        fig, ax = plt.subplots(figsize=(7, 5))
        default_colors = ["k", "C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"]
        default_markers = ["o", "s", "D", "^", "v", "P", "X", "*", "h", "<"]
        for i, (comp, jetpt_dict) in enumerate(sorted(data.items())):
            if not jetpt_dict:
                continue
            pts  = sorted(jetpt_dict.keys())
            mus  = [jetpt_dict[p][0] for p in pts]
            errs = [jetpt_dict[p][1] for p in pts]
            if component_styles and comp in component_styles:
                style = component_styles[comp]
            else:
                style = dict(color=default_colors[i % len(default_colors)],
                             marker=default_markers[i % len(default_markers)],
                             linestyle="-", label=comp)
            ax.errorbar(pts, mus, yerr=errs, capsize=3, markersize=7, **style)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=11)
        ax.set_xscale("log")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8, loc="best", framealpha=0.9)
        fig.tight_layout()
        os.makedirs(outdir, exist_ok=True)
        fig.savefig(os.path.join(outdir, filename))
        plt.close(fig)

    def plot_all_mpv_summaries(self):
        base = os.path.join(self.base_plot_dir, "mpv_summary")
        basic_styles = {
            "full": dict(color="gray",  marker="o", linestyle="-", label="full (passed cut)"),
            "rad":  dict(color="black", marker="s", linestyle="-", label="radiator"),
            "AA":   dict(color="C0",    marker="D", linestyle="-", label="AxA"),
            "BB":   dict(color="C1",    marker="^", linestyle="-", label="BxB"),
            "AB":   dict(color="C2",    marker="v", linestyle="-", label="AxB"),
        }
        for case_key in sorted(self.mpv_data.keys()):
            parts = case_key.split("_")
            # basic_data_{cut}_ww{den}  or  basic_ptRL_data_{cut}_ww{den}
            ptRL = (parts[1] == "ptRL")
            offset = 2 if ptRL else 1
            cut_suffix = parts[offset + 1]
            den_token  = parts[offset + 2]
            obs = "p_TR_L" if ptRL else "R_L"
            ylabel = f"Peak position in {obs}"
            title = (f"MPV summary — data, {cut_suffix}, {den_token}"
                     + (" (ptRL)" if ptRL else ""))
            outdir = os.path.join(base, "data", cut_suffix, "ptRL" if ptRL else "RL")
            fname = (f"mpv_basic{'_ptRL' if ptRL else ''}_data"
                     f"_{cut_suffix}_{den_token}.pdf")
            self._plot_mpv_summary(case_key, title, "Jet p_T [GeV/c]", ylabel,
                                   outdir, fname, component_styles=basic_styles)
        print(f"MPV summary plots saved under: {base}")

    # -------------------------------------------------------------------------
    # Main loop
    # -------------------------------------------------------------------------

    def plot(self):
        self.open_rootfile()
        self._init_mpv_storage()

        for i, jetpt in enumerate(self.target_jet_pts):
            for cut_mode, z_cut in self.cut_modes:
                self.cut_mode = cut_mode
                self.z_cut = z_cut
                print(f"Processing {cut_mode} mode, z_cut={z_cut}, "
                      f"jetpt={self._ptrange_label(jetpt)}...")

                for den_weight in ["jet", "rad"]:
                    self.plot_basic(jetpt, z_cut, den_weight)
                    self._collect_mpv_basic(jetpt, z_cut, den_weight)
                    if self.cut_mode == "sd":
                        self.plot_rg(jetpt, z_cut, den_weight)

                    if i == 0:
                        self.plot_basic_acrosspt(z_cut, den_weight)

                self.plot_CAB(jetpt, z_cut)
                self.plot_CAB_ptRL(jetpt, z_cut)
                self.plot_rad_prop("pt",   jetpt, z_cut)
                self.plot_rad_prop("lnkt", jetpt, z_cut)

        print("\nGenerating MPV summary plots...")
        self.plot_all_mpv_summaries()


if __name__ == "__main__":
    plotter = PlotDataCurves()
    plotter.plot()