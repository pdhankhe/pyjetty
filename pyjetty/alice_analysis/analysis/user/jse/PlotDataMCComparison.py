# Plot overlapping DATA / MC (PYTHIA + Herwig) curves for JSE analysis,
# with an MC/data ratio panel underneath each plot.
#
# Data: ranges (single AnalysisResults.root), MC: single-value jetpt files per generator.
# Data drawn as points; inclusive PYTHIA and inclusive Herwig drawn as two line styles.

import os
import ROOT
import math
import numpy as np
from enum import IntEnum

import mpv_plotting as mpv

class Color(IntEnum):
    BLUE   = ROOT.TColor.GetColor("#1f77b4")
    GREEN  = ROOT.TColor.GetColor("#2ca02c")
    RED    = ROOT.TColor.GetColor("#d62728")
    ORANGE = ROOT.TColor.GetColor("#ff7f0e")
    MAGENTA = ROOT.TColor.GetColor("#e377c2")
    CYAN    = ROOT.TColor.GetColor("#17becf")


class PlotDataMCComparison:
    def __init__(self):
        ROOT.gROOT.SetBatch(True)
        ROOT.gStyle.SetLegendBorderSize(0)
        ROOT.gStyle.SetLegendFillColor(0)
        ROOT.gStyle.SetPadGridX(1)
        ROOT.gStyle.SetPadGridY(1)

        # ------------------------------------------------------------------
        # Flexible data-range <-> MC-value pairing.
        # Each entry: (data_jetpt_range_tuple, mc_jetpt_value)
        # Add more rows here later for additional bins.
        # ------------------------------------------------------------------
        self.pt_pairs = [
            ((10, 20), None),
            ((20, 40), None),
            ((40, 60), None),
            ((50, 60), ["50", "50_55"]),
            ((60, 80), None),
            ((80, 100), None),
            ((100, 120), ["100", "100_110"]),
            ((120, 150), None),
            ((150, 200), None)
        ]

        self.cut_modes   = [("sd", 0.1)]  # , ("maxkt", None)]
        self.den_weights = ["jet", "rad"]

        # inclusive-only MC for the overlay
        self.mc_parton = "inclusive"
        # (gen_string, display_label, ROOT_linestyle)
        self.mc_generators = [
            ("pythia", "PYTHIA8",  ROOT.kSolid),
            ("herwig", "Herwig7",  ROOT.kDashed),
        ]

        # ------------------------------------------------------------------
        # Files
        # ------------------------------------------------------------------
        self.data_rootfile_path = (
            # "/global/cfs/cdirs/alice/blianggi/mypyjetty/analysis/testing/"
            # "AnalysisResults.root"
            "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data/55778272/AnalysisResultsMerged.root"

        )
        self.mc_rootfile_template = (
            "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/rootfiles/"
            "jse_preliminary_curves_{gen}_jetpt{jetpt}.root"
        )

        self.mc_rootfile_template_narrow = (
            "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/rootfiles/narrowerbins1.1/"
            "jse_preliminary_curves_{gen}_jetpt{jetpt}.root"
        )

        self.base_plot_dir = (
            "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots/"
            "data/mc_comparison"
        )

        # ------------------------------------------------------------------
        # MC histogram name templates (confirmed by user).
        #   full: (parton, jetpt, cut)
        #   rad/AA/BB/AB/CAB: (parton, jetpt, cut, den_weight)
        # ------------------------------------------------------------------
        self.MC_FULL = "hist_full_{}_jetpt{}_{}"
        self.MC_RAD  = "hist_rad_{}_jetpt{}_{}_ww{}pt"
        self.MC_AA   = "hist_AA_{}_jetpt{}_{}_ww{}pt"
        self.MC_BB   = "hist_BB_{}_jetpt{}_{}_ww{}pt"
        self.MC_AB   = "hist_AB_{}_jetpt{}_{}_ww{}pt"
        self.MC_CAB  = "CAB_{}_jetpt{}_{}_ww{}pt"

        # den_weight key -> token used in DATA histogram names
        self.WEIGHT_TOKEN = {"jet": "wwjetpt", "rad": "wwradpt"}

        # ratio-panel y-axis range
        self.ratio_ymin = 0.0
        self.ratio_ymax = 2.0

        self.cut_mode = ""
        self.z_cut = None

        self.data_rootfile = None
        self.mc_rootfiles = {}   # (gen, mc_jetpt) -> TFile

        self._persistent = []    # keep canvases, pads, lines, hists alive
        self._canvas_counter = 0
        self._counter_cache = {}   # (pr, cut) -> (num_jets, avg_jet_pt)

    def _get_norm_factors(self, jetpt, z_cut):
        """Return (num_jets_passed_cut, avg_jet_pt) for this slice, cached.

        counters bin 1 = num_jets
        counters bin 2 = num_jets_passed_cut
        counters bin 3 = sum_jetpt_passed_cut
        avg_jet_pt = sum_jetpt_passed_cut / num_jets_passed_cut  (correct global average)
        Returns (None, None) if counters missing or num_jets <= 0.
        """
        cut = self.get_cut_suffix(z_cut)
        pr = self._data_ptrange_token(jetpt)
        key = (pr, cut)
        if key in self._counter_cache:
            return self._counter_cache[key]

        counters = self._get_from(self.data_rootfile, f"counters_{pr}_{cut}")
        if not counters:
            print(f"WARNING: counters_{pr}_{cut} not found")
            self._counter_cache[key] = (None, None)
            return (None, None)

        num_jets = counters.GetBinContent(2)
        sum_jetpt = counters.GetBinContent(3)
        if num_jets <= 0:
            self._counter_cache[key] = (None, None)
            return (None, None)

        avg_jet_pt = sum_jetpt / num_jets
        self._counter_cache[key] = (num_jets, avg_jet_pt)
        return (num_jets, avg_jet_pt)

    def _normalize_eec(self, hist, jetpt, z_cut, ptRL=False):
        """Scale a fetched EEC/ptRL hist in memory (does not touch the file).

        Regular:  1/num_jets with "width"
        ptRL:     log(<pt>)/<pt> content factor, then 1/num_jets with "width"
        """
        if hist is None:
            return None
        num_jets, avg_jet_pt = self._get_norm_factors(jetpt, z_cut)
        if num_jets is None:
            return hist
        if ptRL:
            hist.Scale(math.log(avg_jet_pt) / avg_jet_pt)
        hist.Scale(1.0 / num_jets, "width")
        return hist

    def get_cab_data(self, jetpt, z_cut, den_weight, ptRL=False):
        """Recalculate C_{AB} for data using the method from PlotDataHists.py.
        C_{AB} = AB / sqrt(AA * BB)
        """
        cut = self.get_cut_suffix(z_cut)
        pr = self._data_ptrange_token(jetpt)
        w = self.WEIGHT_TOKEN[den_weight]

        # Get components
        h_AA = self.get_data_component("AA", jetpt, z_cut, den_weight, ptRL=ptRL)
        h_BB = self.get_data_component("BB", jetpt, z_cut, den_weight, ptRL=ptRL)
        h_AB = self.get_data_component("AB", jetpt, z_cut, den_weight, ptRL=ptRL)

        if not all([h_AA, h_BB, h_AB]):
            return None

        name = f"{self._data_ptRL('CAB_recalc', ptRL)}{pr}_{cut}_{w}"
        hist_CAB = h_AB.Clone(name)
        hist_CAB.SetTitle("C_{AB};R_{L};AxB / #sqrt{AxA #times BxB}")
        hist_CAB.SetDirectory(0)
        hist_CAB.Reset()

        for i in range(1, hist_CAB.GetNbinsX() + 1):
            aa, bb, ab = h_AA.GetBinContent(i), h_BB.GetBinContent(i), h_AB.GetBinContent(i)
            sig_aa, sig_bb, sig_ab = h_AA.GetBinError(i), h_BB.GetBinError(i), h_AB.GetBinError(i)
            if aa > 0 and bb > 0:
                denom = math.sqrt(aa * bb)
                result = ab / denom
                term_ab = sig_ab / denom
                term_aa = (result * sig_aa) / (2 * aa)
                term_bb = (result * sig_bb) / (2 * bb)
                hist_CAB.SetBinContent(i, result)
                hist_CAB.SetBinError(i, np.sqrt(term_ab**2 + term_aa**2 + term_bb**2))
            else:
                hist_CAB.SetBinContent(i, 0)
                hist_CAB.SetBinError(i, 0)
        return hist_CAB

    # =====================================================================
    # Naming helpers
    # =====================================================================

    def get_cut_suffix(self, z_cut):
        return f"sd{z_cut}" if self.cut_mode == "sd" else "maxkt"

    def get_cut_label(self, z_cut):
        return f"SD z_{{cut}} = {z_cut}" if self.cut_mode == "sd" else "maxkt selection"

    def _mc_path(self, gen, mc_jetpt):
        """Full MC rootfile path. Range-style tokens (with '_', e.g. '50_55')
        live in the narrower-bins folder; single-value tokens use the normal one."""
        tmpl = (self.mc_rootfile_template_narrow if "_" in str(mc_jetpt)
                else self.mc_rootfile_template)
        return tmpl.format(gen=gen, jetpt=mc_jetpt)
        
    @staticmethod
    def _data_ptrange_token(jetpt):
        lo, hi = jetpt
        return f"jetpt{lo}_{hi}"

    @staticmethod
    def _data_ptrange_label(jetpt):
        lo, hi = jetpt
        return f"{lo}-{hi}"

    @staticmethod
    def _data_ptrange_tag(jetpt):
        lo, hi = jetpt
        return f"jetpt{lo}_{hi}"

    @staticmethod
    def _mc_pt_tag(mc_jetpt):
        return f"jetpt{mc_jetpt}"

    @staticmethod
    def _data_ptRL(prefix, ptRL):
        return f"{prefix}_ptRL_" if ptRL else f"{prefix}_"

    # =====================================================================
    # File handling
    # =====================================================================

    def open_files(self):
        self.data_rootfile = ROOT.TFile.Open(self.data_rootfile_path)
        if not self.data_rootfile or self.data_rootfile.IsZombie():
            raise RuntimeError(f"Could not open DATA file {self.data_rootfile_path}")

        # flatten: each pt_pair's MC part may be a single token, a list, or None
        mc_jetpts = set()
        for _, mc in self.pt_pairs:
            if mc is None:
                continue
            if isinstance(mc, (list, tuple)):
                mc_jetpts.update(str(m) for m in mc)
            else:
                mc_jetpts.add(str(mc))
        mc_jetpts = sorted(mc_jetpts)
        print("analyzing mc jet pts:", mc_jetpts)

        for gen, _, _ in self.mc_generators:
            for mc_jetpt in mc_jetpts:
                path = self._mc_path(gen, mc_jetpt)
                print("Opening MC file:", path)
                f = ROOT.TFile.Open(path)
                if not f or f.IsZombie():
                    print(f"WARNING: could not open MC file {path}")
                    self.mc_rootfiles[(gen, mc_jetpt)] = None
                else:
                    self.mc_rootfiles[(gen, mc_jetpt)] = f

    @staticmethod
    def _get_from(tfile, name):
        if tfile is None:
            return None
        h = tfile.Get(name)
        if h:
            h.SetDirectory(0)
            return h
        return None

    # =====================================================================
    # DATA histogram retrieval  (mirrors PlotDataCurves naming)
    # =====================================================================

    def get_data_component(self, comp, jetpt, z_cut, den_weight, ptRL=False):
        cut = self.get_cut_suffix(z_cut)
        pr = self._data_ptrange_token(jetpt)
        w = self.WEIGHT_TOKEN[den_weight]
        if comp == "rad":
            name = f"{self._data_ptRL('hist_rad', ptRL)}{pr}_{cut}_{w}"
        elif comp == "AA":
            name = f"{self._data_ptRL('hist_AA', ptRL)}{pr}_{cut}_{w}"
        elif comp == "BB":
            name = f"{self._data_ptRL('hist_BB', ptRL)}{pr}_{cut}_{w}"
        elif comp == "AB":
            name = f"{self._data_ptRL('hist_AB', ptRL)}{pr}_{cut}_{w}"
        elif comp == "CAB":
            name = f"{self._data_ptRL('CAB', ptRL)}{pr}_{cut}_{w}"
        else:
            raise ValueError(f"Unknown component {comp}")
        return self._get_from(self.data_rootfile, name)

    # =====================================================================
    # MC histogram retrieval
    # =====================================================================

    def get_mc_component(self, comp, gen, mc_jetpt, z_cut, den_weight, ptRL=False):
        cut = self.get_cut_suffix(z_cut)
        tfile = self.mc_rootfiles.get((gen, mc_jetpt))
        print("Getting mc rootfile:", "gen:", gen, "mc_jetpt:", mc_jetpt, "-->", tfile)
        pt = mc_jetpt

        base_map = {
            "rad": self.MC_RAD, "AA": self.MC_AA,
            "BB": self.MC_BB, "AB": self.MC_AB, "CAB": self.MC_CAB,
        }
        tmpl = base_map[comp]  # all requested comps are weighted
        # Insert ptRL modifier into the leading base part of the template.
        # e.g. "hist_AA_{}_jetpt{}_{}_ww{}pt" -> "hist_AA_ptRL_{}_jetpt{}_{}_ww{}pt"
        if ptRL:
            head, rest = tmpl.split("_{}", 1)
            tmpl = f"{head}_ptRL_{{}}{rest}"
        name = tmpl.format(self.mc_parton, pt, cut, den_weight)
        return self._get_from(tfile, name)

    # =====================================================================
    # Canvas / output helpers
    # =====================================================================

    def make_canvas(self, base_name, tag="", z_cut="", den_weight="",
                    title=None, w=800, h=700):
        self._canvas_counter += 1
        unique = (f"{base_name}_{tag}_{self.get_cut_suffix(z_cut)}"
                  f"_{den_weight}_{self._canvas_counter}")
        c = ROOT.TCanvas(unique, title or base_name, w, h)
        self._persistent.append(c)
        return c

    def get_output_dir(self, plot_type=None):
        parts = [self.base_plot_dir, self.get_cut_suffix(self.z_cut)]
        if plot_type:
            parts.extend(plot_type.split('/'))
        path = os.path.join(*parts)
        os.makedirs(path, exist_ok=True)
        return path

    def _save_canvas(self, canvas, filename, plot_type=None):
        canvas.SaveAs(os.path.join(self.get_output_dir(plot_type), filename))

    def FormatHistPoints(self, hist, color, markerstyle=20):
        if hist is None:
            return
        hist.SetLineColor(color)
        hist.SetMarkerColor(color)
        hist.SetMarkerStyle(markerstyle)
        hist.SetMarkerSize(0.9)
        hist.SetLineWidth(1)

    def FormatHistLine(self, hist, color, linestyle):
        if hist is None:
            return
        hist.SetLineColor(color)
        hist.SetLineStyle(linestyle)
        hist.SetLineWidth(2)
        hist.SetMarkerStyle(0)

    def MakeEventLeg(self, jetpt_label, z_cut, den_weight="",
                     x1=0.15, y1=0.70, x2=0.45, y2=0.88):
        leg = ROOT.TLegend(x1, y1, x2, y2)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetMargin(0)
        leg.AddEntry(ROOT.nullptr, "pp data vs MC, R = 0.4 jets", "")
        leg.AddEntry(ROOT.nullptr, f"jet p_{{T}} = {jetpt_label} GeV/c", "")
        leg.AddEntry(ROOT.nullptr, self.get_cut_label(z_cut), "")
        if den_weight:
            leg.AddEntry(ROOT.nullptr,
                         f"using weight p_{{T,1}}p_{{T,2}} / p_{{T,{den_weight}}}^{{2}}", "")
        return leg

    # =====================================================================
    # Color families for the 2x2 quadrant canvas
    #   data = darkest shade, MC generators = lighter variants
    # =====================================================================

    def _color_family(self, comp):
        """Return (data_color, [mc_color_0, mc_color_1]) for a component."""
        families = {
            # rad -> black / grays
            "rad": (ROOT.kBlack,
                    [ROOT.TColor.GetColor("#555555"),
                     ROOT.TColor.GetColor("#999999")]),
            # AA -> blues
            "AA":  (ROOT.TColor.GetColor("#08306b"),
                    [ROOT.TColor.GetColor("#2171b5"),
                     ROOT.TColor.GetColor("#6baed6")]),
            # BB -> oranges
            "BB":  (ROOT.TColor.GetColor("#8c2d04"),
                    [ROOT.TColor.GetColor("#ec7014"),
                     ROOT.TColor.GetColor("#fec44f")]),
            # AB -> greens
            "AB":  (ROOT.TColor.GetColor("#00441b"),
                    [ROOT.TColor.GetColor("#238b45"),
                     ROOT.TColor.GetColor("#74c476")]),
        }
        return families[comp]

    def _draw_component_in_pad(self, canvas, cell, comp, data_jetpt, mc_jetpt,
                               z_cut, den_weight, ptRL=False):
        """Draw one component's data+MC overlay + ratio into a cell region.

        `cell` = (x_lo, y_lo, x_hi, y_hi) in canvas NDC.
        Returns True if anything was drawn.
        """
        cut = self.get_cut_suffix(z_cut)
        data_tag = self._data_ptrange_tag(data_jetpt)
        suffix = "_ptRL" if ptRL else ""
        obs = "p_{T}R_{L}" if ptRL else "R_{L}"

        data_color, mc_colors = self._color_family(comp)

        # normalize MC bins (single token, list, or None)
        if mc_jetpt is None:
            mc_bins = []
        elif isinstance(mc_jetpt, (list, tuple)):
            mc_bins = list(mc_jetpt)
        else:
            mc_bins = [mc_jetpt]

        # ---- retrieve ----
        if comp == "CAB":
            h_data = self.get_cab_data(data_jetpt, z_cut, den_weight, ptRL=ptRL)
        else:
            h_data = self.get_data_component(comp, data_jetpt, z_cut, den_weight, ptRL=ptRL)
        mc_raw = []  # (hist, gen_label, mcb, linestyle)
        for gen, gen_label, ls in self.mc_generators:
            for mcb in mc_bins:
                h = self.get_mc_component(comp, gen, mcb, z_cut, den_weight, ptRL=ptRL)
                if h is not None:
                    mc_raw.append((h, gen_label, mcb, ls))

        if h_data is None and not mc_raw:
            return False

        # ---- style clones ----
        if h_data is not None:
            h_data = h_data.Clone(f"q_data_{comp}_{data_tag}_{cut}_ww{den_weight}{suffix}")
            h_data.SetDirectory(0)
            if comp != "CAB":
                h_data = self._normalize_eec(h_data, data_jetpt, z_cut, ptRL=ptRL)
            self.FormatHistPoints(h_data, data_color, markerstyle=20)
            self._persistent.append(h_data)

        styled_mc = []  # (hist, legend_label, color)
        gen_order = [gl for _, gl, _ in self.mc_generators]
        line_styles = [1, 2, 3, 7, 9]
        for (h, gen_label, mcb, ls) in mc_raw:
            gi = gen_order.index(gen_label)
            bi = mc_bins.index(mcb)
            color = mc_colors[gi] if gi < len(mc_colors) else mc_colors[-1]
            style = line_styles[bi % len(line_styles)]
            hc = h.Clone(f"q_mc_{gen_label}_{mcb}_{comp}_{data_tag}_{cut}"
                         f"_ww{den_weight}{suffix}")
            hc.SetDirectory(0)
            self.FormatHistLine(hc, color, style)
            leg_label = (f"{gen_label} (incl., jetpt{mcb})"
                         if len(mc_bins) > 1 else f"{gen_label} (incl.)")
            styled_mc.append((hc, leg_label, color))
            self._persistent.append(hc)

        drawn = ([h_data] if h_data is not None else []) + [h for h, _, _ in styled_mc]

        # ---- build main + ratio sub-pads inside the cell ----
        x_lo, y_lo, x_hi, y_hi = cell
        y_split = y_lo + 0.30 * (y_hi - y_lo)   # ratio occupies bottom 30% of the cell

        self._canvas_counter += 1
        tag = f"{comp}_{data_tag}_{cut}_ww{den_weight}{suffix}_{self._canvas_counter}"

        canvas.cd()
        pad_main = ROOT.TPad(f"qmain_{tag}", "main", x_lo, y_split, x_hi, y_hi)
        pad_main.SetLeftMargin(0.16)
        pad_main.SetRightMargin(0.04)
        pad_main.SetTopMargin(0.06)
        pad_main.SetBottomMargin(0.02)
        pad_main.SetLogx()
        if comp == "CAB":
            pad_main.SetLogy()
        pad_main.Draw()
        self._persistent.append(pad_main)

        canvas.cd()
        pad_ratio = ROOT.TPad(f"qratio_{tag}", "ratio", x_lo, y_lo, x_hi, y_split)
        pad_ratio.SetLeftMargin(0.16)
        pad_ratio.SetRightMargin(0.04)
        pad_ratio.SetTopMargin(0.02)
        pad_ratio.SetBottomMargin(0.35)
        pad_ratio.SetLogx()
        pad_ratio.SetGridy()
        pad_ratio.Draw()
        self._persistent.append(pad_ratio)

        # ------------------- main pad -------------------
        pad_main.cd()
        y_max = max((h.GetMaximum() for h in drawn if h), default=1.0)
        first_hist = drawn[0]
        if comp == "CAB":
            positive_mins = [h.GetMinimum(0) for h in drawn if h]
            first_hist.SetMaximum(y_max * 3.0)
            first_hist.SetMinimum(max(1e-6, min(positive_mins, default=1e-6)))
        else:
            first_hist.SetMaximum(y_max * 1.5)
            first_hist.SetMinimum(0)
        ytitle = "C_{AB}" if comp == "CAB" else f"(1/N_{{jets}}) dN/d{obs}"
        first_hist.GetYaxis().SetTitle(ytitle)
        first_hist.GetYaxis().SetTitleSize(0.06)
        first_hist.GetYaxis().SetTitleOffset(1.1)
        first_hist.GetYaxis().SetLabelSize(0.05)
        first_hist.GetXaxis().SetLabelSize(0)   # hide x labels on main pad

        first = True
        if h_data is not None:
            h_data.Draw("PE")
            first = False
        for h, _, _ in styled_mc:
            h.Draw("HIST SAME" if not first else "HIST")
            first = False
        if h_data is not None:
            h_data.Draw("PE SAME")

        if comp == "CAB":
            line = ROOT.TLine(first_hist.GetXaxis().GetXmin(), 1.0,
                              first_hist.GetXaxis().GetXmax(), 1.0)
            line.SetLineColor(ROOT.kGray + 2)
            line.SetLineStyle(ROOT.kDashed)
            line.Draw("SAME")
            self._persistent.append(line)

        # per-quadrant legend
        comp_title = {"rad": "radiator", "AA": "AA", "BB": "BB",
                      "AB": "AB", "CAB": "C_{AB}"}.get(comp, comp)
        leg = ROOT.TLegend(0.55, 0.66, 0.94, 0.92)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetTextSize(0.055)
        leg.SetHeader(comp_title)
        if h_data is not None:
            leg.AddEntry(h_data, "data", "pe")
        for h, leg_label, _ in styled_mc:
            leg.AddEntry(h, leg_label, "l")
        leg.Draw()
        self._persistent.append(leg)

        # ------------------- ratio pad (MC / data) -------------------
        pad_ratio.cd()
        ratio_hists = []
        if h_data is not None:
            for h, leg_label, color in styled_mc:
                r = h.Clone(f"qratio_{comp}_{data_tag}_{cut}"
                            f"_ww{den_weight}{suffix}_{len(ratio_hists)}")
                r.SetDirectory(0)
                r.Divide(h_data)   # MC / data
                r.SetLineColor(color)
                r.SetLineStyle(h.GetLineStyle())
                r.SetLineWidth(2)
                r.SetMarkerStyle(0)
                ratio_hists.append(r)
                self._persistent.append(r)

        if ratio_hists:
            r0 = ratio_hists[0]
            r0.SetMinimum(self.ratio_ymin)
            r0.SetMaximum(self.ratio_ymax)
            r0.GetYaxis().SetTitle("MC / data")
            r0.GetYaxis().SetNdivisions(505)
            r0.GetYaxis().SetTitleSize(0.13)
            r0.GetYaxis().SetTitleOffset(0.5)
            r0.GetYaxis().SetLabelSize(0.11)
            r0.GetXaxis().SetTitle(obs)
            r0.GetXaxis().SetTitleSize(0.14)
            r0.GetXaxis().SetTitleOffset(1.0)
            r0.GetXaxis().SetLabelSize(0.11)
            r0.Draw("HIST")
            for r in ratio_hists[1:]:
                r.Draw("HIST SAME")

            unity = ROOT.TLine(r0.GetXaxis().GetXmin(), 1.0,
                               r0.GetXaxis().GetXmax(), 1.0)
            unity.SetLineColor(ROOT.kGray + 2)
            unity.SetLineStyle(ROOT.kDashed)
            unity.Draw("SAME")
            self._persistent.append(unity)
        else:
            txt = ROOT.TLatex()
            txt.SetNDC()
            txt.SetTextSize(0.15)
            txt.DrawLatex(0.25, 0.5, "no data for ratio")
            self._persistent.append(txt)

        return True

    # =====================================================================
    # 2x2 quadrant canvas: rad / AA / BB / AB
    # =====================================================================

    def plot_quadrants(self, data_jetpt, mc_jetpt, z_cut, den_weight, ptRL=False):
        cut = self.get_cut_suffix(z_cut)
        data_tag = self._data_ptrange_tag(data_jetpt)
        data_label = self._data_ptrange_label(data_jetpt)
        suffix = "_ptRL" if ptRL else ""

        canvas = self.make_canvas("can_quadrants", tag=data_tag, z_cut=z_cut,
                                  den_weight=den_weight,
                                  title="rad/AA/BB/AB data vs MC",
                                  w=1000, h=1000)

        # Cell regions in canvas NDC: (x_lo, y_lo, x_hi, y_hi).
        # Leave a small strip at top (y up to 0.97) for the overall header.
        #   top-left  = rad      top-right  = AA
        #   bot-left  = BB       bot-right  = AB
        top_y = 0.97
        cells = {
            "rad": (0.00, 0.50, 0.50, top_y),
            "AA":  (0.50, 0.50, 1.00, top_y),
            "BB":  (0.00, 0.00, 0.50, 0.50),
            "AB":  (0.50, 0.00, 1.00, 0.50),
        }

        any_drawn = False
        for comp, cell in cells.items():
            drawn = self._draw_component_in_pad(canvas, cell, comp,
                                                data_jetpt, mc_jetpt,
                                                z_cut, den_weight, ptRL=ptRL)
            any_drawn = any_drawn or drawn

        if not any_drawn:
            print(f"WARNING: no data or MC for quadrant plot "
                  f"(ptRL={ptRL}, ww{den_weight}, cut={cut}) - skipping")
            return

        # Overall header across the top
        canvas.cd()
        header = ROOT.TLatex()
        header.SetNDC()
        header.SetTextSize(0.020)
        header.SetTextAlign(22)
        header.DrawLatex(0.5, 0.985,
                         f"pp data vs MC, R = 0.4, jet p_{{T}} = {data_label} GeV/c, "
                         f"{self.get_cut_label(z_cut)}, ww{den_weight}")
        self._persistent.append(header)

        fname = (f"datamc_quadrants{suffix}_{data_tag}"
                 f"_vs_mc{self._mc_pt_tag(mc_jetpt)}_R0.4_{cut}_ww{den_weight}.pdf")
        self._save_canvas(canvas, fname, plot_type="quadrants")
        
    # =====================================================================
    # Core overlay for one component (with MC/data ratio panel)
    # =====================================================================

    def plot_component(self, comp, data_jetpt, mc_jetpt, z_cut, den_weight,
                       ptRL=False):
        cut = self.get_cut_suffix(z_cut)
        data_tag = self._data_ptrange_tag(data_jetpt)
        data_label = self._data_ptrange_label(data_jetpt)
        suffix = "_ptRL" if ptRL else ""

        # normalize MC bins to a list (single token, list, or None all OK)
        if mc_jetpt is None:
            mc_bins = []
        elif isinstance(mc_jetpt, (list, tuple)):
            mc_bins = list(mc_jetpt)
        else:
            mc_bins = [mc_jetpt]
        print("MC BINS", mc_bins)

        # ---- retrieve ----
        if comp == "CAB":
            h_data = self.get_cab_data(data_jetpt, z_cut, den_weight, ptRL=ptRL)
        else:
            h_data = self.get_data_component(comp, data_jetpt, z_cut, den_weight, ptRL=ptRL)
        mc_raw = []  # (hist, gen_label, mcb, linestyle)
        for gen, gen_label, ls in self.mc_generators:
            for mcb in mc_bins:
                h = self.get_mc_component(comp, gen, mcb, z_cut, den_weight, ptRL=ptRL)
                # TO DO!!! for some reason the narrower bins are registering here as none...
                if h is not None:
                    print("adding to mc_raw:", "gen_label:", gen_label, "mcb:", mcb, "ls:", ls)
                    mc_raw.append((h, gen_label, mcb, ls))

        if h_data is None and not mc_raw:
            print(f"WARNING: no data or MC for {comp} "
                  f"(ptRL={ptRL}, ww{den_weight}, cut={cut}) - skipping")
            return

        # ---- style (clone to avoid mutating cached hists) ----
        if h_data is not None:
            h_data = h_data.Clone(f"data_{comp}_{data_tag}_{cut}_ww{den_weight}{suffix}")
            h_data.SetDirectory(0)
            if comp != "CAB":
                h_data = self._normalize_eec(h_data, data_jetpt, z_cut, ptRL=ptRL)
            self.FormatHistPoints(h_data, ROOT.kBlack, markerstyle=20)
            self._persistent.append(h_data)

        # a color per generator, a line style per mc-bin within that generator
        gen_colors = [Color.RED, Color.BLUE, Color.GREEN, Color.MAGENTA,
                      Color.CYAN, Color.ORANGE]
        line_styles = [1, 2, 3, 7, 9]  # solid, dashed, dotted, ...

        # map each generator label / bin to a stable index for color/style
        gen_order = [gl for _, gl, _ in self.mc_generators]
        bin_order = mc_bins

        styled_mc = []  # (hist, legend_label, color)
        for (h, gen_label, mcb, ls) in mc_raw:
            gi = gen_order.index(gen_label)
            bi = bin_order.index(mcb)
            color = gen_colors[gi % len(gen_colors)]
            style = line_styles[bi % len(line_styles)]
            hc = h.Clone(f"mc_{gen_label}_{mcb}_{comp}_{data_tag}_{cut}"
                         f"_ww{den_weight}{suffix}")
            hc.SetDirectory(0)
            self.FormatHistLine(hc, color, style)
            # legend shows bin only when there's more than one, to avoid clutter
            leg_label = (f"{gen_label} (incl., jetpt{mcb})"
                         if len(mc_bins) > 1 else f"{gen_label} (incl.)")
            styled_mc.append((hc, leg_label, color))
            self._persistent.append(hc)
            print("gen_label:", gen_label, "leg_label:", leg_label, "color:", color)

        drawn_main = ([h_data] if h_data is not None else []) + [h for h, _, _ in styled_mc]

        # =================================================================
        # Canvas with two pads: main (top) + ratio (bottom)
        # =================================================================
        cname = f"can_{comp}{suffix}"
        canvas = self.make_canvas(cname, tag=data_tag, z_cut=z_cut,
                                  den_weight=den_weight, title=f"{comp} data/MC")

        pad1 = ROOT.TPad(f"pad1_{cname}_{self._canvas_counter}", "top", 0, 0.30, 1, 1.0)
        pad1.SetBottomMargin(0.02)
        pad1.SetLogx()
        if comp == "CAB":
            pad1.SetLogy()
        pad1.Draw()
        self._persistent.append(pad1)

        canvas.cd()
        pad2 = ROOT.TPad(f"pad2_{cname}_{self._canvas_counter}", "bot", 0, 0.0, 1, 0.30)
        pad2.SetTopMargin(0.02)
        pad2.SetBottomMargin(0.35)
        pad2.SetLogx()
        pad2.SetGridy()
        pad2.Draw()
        self._persistent.append(pad2)

        # ------------------- top pad -------------------
        pad1.cd()
        obs = "p_{T}R_{L}" if ptRL else "R_{L}"
        y_max = max((h.GetMaximum() for h in drawn_main if h), default=1.0)
        first_hist = drawn_main[0]
        if comp == "CAB":
            positive_mins = [h.GetMinimum(0) for h in drawn_main if h]
            first_hist.SetMaximum(y_max * 3.0)
            first_hist.SetMinimum(max(1e-6, min(positive_mins, default=1e-6)))
        else:
            first_hist.SetMaximum(y_max * 1.4)
            first_hist.SetMinimum(0)
        ytitle = "C_{AB}" if comp == "CAB" else f"(1/N_{{jets}}) dN/d{obs}"
        first_hist.GetYaxis().SetTitle(ytitle)
        first_hist.GetYaxis().SetTitleSize(0.05)
        first_hist.GetYaxis().SetTitleOffset(0.9)
        first_hist.GetXaxis().SetLabelSize(0)  # hide x labels on top pad

        first = True
        if h_data is not None:
            h_data.Draw("PE")
            first = False
        for h, _, _ in styled_mc:
            h.Draw("HIST SAME" if not first else "HIST")
            first = False
        if h_data is not None:
            h_data.Draw("PE SAME")  # data on top

        if comp == "CAB":
            line = ROOT.TLine(first_hist.GetXaxis().GetXmin(), 1.0,
                              first_hist.GetXaxis().GetXmax(), 1.0)
            line.SetLineColor(ROOT.kGray + 2)
            line.SetLineStyle(ROOT.kDashed)
            line.Draw("SAME")
            self._persistent.append(line)

        ev_leg = self.MakeEventLeg(data_label, z_cut, den_weight)
        if comp == "CAB":
            legend = ROOT.TLegend(0.47, 0.62, 0.73, 0.88)
        else:
            legend = ROOT.TLegend(0.62, 0.62, 0.88, 0.88)
        if h_data is not None:
            legend.AddEntry(h_data, "data", "pe")
        for h, leg_label, _ in styled_mc:
            legend.AddEntry(h, leg_label, "l") #f"{leg_label} (incl.)", "l")
        ev_leg.Draw()
        legend.Draw()
        self._persistent.extend([ev_leg, legend])

        # ------------------- bottom pad (MC / data) -------------------
        pad2.cd()
        ratio_hists = []
        if h_data is not None:
            for h, leg_label, color in styled_mc:
                r = h.Clone(f"ratio_{comp}_{data_tag}_{cut}"
                            f"_ww{den_weight}{suffix}_{len(ratio_hists)}")
                r.SetDirectory(0)
                r.Divide(h_data)   # MC / data
                r.SetLineColor(color)
                r.SetLineStyle(h.GetLineStyle())
                r.SetLineWidth(2)
                r.SetMarkerStyle(0)
                ratio_hists.append(r)
                self._persistent.append(r)

        if ratio_hists:
            r0 = ratio_hists[0]
            r0.SetMinimum(self.ratio_ymin)
            r0.SetMaximum(self.ratio_ymax)
            r0.GetYaxis().SetTitle("MC / data")
            r0.GetYaxis().SetNdivisions(505)
            r0.GetYaxis().SetTitleSize(0.11)
            r0.GetYaxis().SetTitleOffset(0.4)
            r0.GetYaxis().SetLabelSize(0.09)
            r0.GetXaxis().SetTitle(obs)
            r0.GetXaxis().SetTitleSize(0.12)
            r0.GetXaxis().SetTitleOffset(1.0)
            r0.GetXaxis().SetLabelSize(0.09)
            r0.Draw("HIST")
            for r in ratio_hists[1:]:
                r.Draw("HIST SAME")

            unity = ROOT.TLine(r0.GetXaxis().GetXmin(), 1.0,
                               r0.GetXaxis().GetXmax(), 1.0)
            unity.SetLineColor(ROOT.kGray + 2)
            unity.SetLineStyle(ROOT.kDashed)
            unity.Draw("SAME")
            self._persistent.append(unity)
        else:
            # no data -> nothing to divide by; leave a note
            txt = ROOT.TLatex()
            txt.SetNDC()
            txt.SetTextSize(0.15)
            txt.DrawLatex(0.2, 0.5, "no data for ratio")
            self._persistent.append(txt)

        # ---- save ----
        canvas.cd()
        subdir = f"{comp}/ptRL" if ptRL else comp
        mc_tag = "_".join(self._mc_pt_tag(m) for m in mc_bins) if mc_bins else "none"
        fname = (f"datamc_{comp}{suffix}_{data_tag}"
                 f"_vs_mc{mc_tag}_R0.4_{cut}_ww{den_weight}.pdf")
        self._save_canvas(canvas, fname, plot_type=subdir)

    # =====================================================================
    # MPV (Most Probable Value) summaries: data vs PYTHIA vs Herwig
    # =====================================================================

    def _init_mpv_storage(self):
        """Initialize MPV storage.

        Structure: self.mpv_data[case_key][source][jetpt_x] = (mu, mu_err)
          case_key = f"{comp}_{cut_suffix}_ww{den_weight}{'_ptRL' if ptRL else ''}"
          source   in {"data", "pythia", "herwig"}
          jetpt_x  = numeric x used on the summary plot (see _mpv_x below)
        """
        self.mpv_data = {}
        self._mpv_fit_cache = {}   # avoid re-fitting the same histogram

    @staticmethod
    def _mpv_x_data(data_jetpt):
        """Numeric x-position for a DATA range: use the bin center."""
        lo, hi = data_jetpt
        return 0.5 * (lo + hi)

    @staticmethod
    def _mpv_x_mc(mc_jetpt):
        """Numeric x-position for an MC single-value jetpt token.

        mc_jetpt may be "50", "100", or a range token like "50_55".
        """
        s = str(mc_jetpt)
        if "_" in s:
            lo, hi = s.split("_", 1)
            return 0.5 * (float(lo) + float(hi))
        return float(s)

    def _mpv_diag_path(self, source, comp, cut_suffix, den_weight, ptRL, xlabel):
        """Path for a per-histogram fit diagnostic PDF."""
        suffix = "_ptRL" if ptRL else ""
        return os.path.join(
            self.base_plot_dir, cut_suffix, "mpv_summary", "fits",
            f"{comp}{suffix}_ww{den_weight}",
            f"fit_{source}_{comp}{suffix}_ww{den_weight}_x{xlabel}.pdf"
        )

    def _fit_hist_cached(self, hist, cache_key, save_path=None, title=None):
        """Fit a histogram for its MPV, caching by cache_key so we never
        fit the same underlying histogram twice."""
        if cache_key in self._mpv_fit_cache:
            return self._mpv_fit_cache[cache_key]
        mu, mu_err = mpv.fit_mpv(hist, save_path=save_path, plot_title=title)
        self._mpv_fit_cache[cache_key] = (mu, mu_err)
        return mu, mu_err

    def collect_mpv(self, comp, data_jetpt, mc_jetpt, z_cut, den_weight,
                    ptRL=False):
        """Fit data + MC for one (comp, pt-pair, cut, weight, ptRL) and store
        the peak positions.

        mc_jetpt may be a single token or a list of tokens; each MC bin is
        plotted as its own point with a small horizontal offset for visibility.
        """
        cut = self.get_cut_suffix(z_cut)
        suffix = "_ptRL" if ptRL else ""
        case_key = f"{comp}_{cut}_ww{den_weight}{suffix}"
        obs = "p_TR_L" if ptRL else "R_L"

        # normalize mc_jetpt to a list (accept single value or None)
        if mc_jetpt is None:
            mc_bins = []
        elif isinstance(mc_jetpt, (list, tuple)):
            mc_bins = list(mc_jetpt)
        else:
            mc_bins = [mc_jetpt]

        # ---- DATA ----
        x_data = self._mpv_x_data(data_jetpt)
        h_data = self.get_data_component(comp, data_jetpt, z_cut, den_weight,
                                         ptRL=ptRL)
        if h_data is not None:
            ck = ("data", comp, self._data_ptrange_token(data_jetpt), cut,
                  den_weight, ptRL)
            path = self._mpv_diag_path("data", comp, cut, den_weight, ptRL,
                                       self._data_ptrange_token(data_jetpt))
            title = (f"data | {comp} | {self._data_ptrange_label(data_jetpt)} | "
                     f"{cut} | ww{den_weight} | {obs}")
            mu, mu_err = self._fit_hist_cached(h_data, ck, save_path=path,
                                               title=title)
            if mu is not None:
                (self.mpv_data.setdefault(case_key, {})
                     .setdefault("data", {})[x_data]) = (mu, mu_err)

        # ---- MC generators (one point per generator x mc-bin) ----
        # small horizontal offset per MC point so overlapping bins are visible
        MPV_X_OFFSET = 0.6
        n_gen = len(self.mc_generators)
        for gi, (gen, gen_label, _) in enumerate(self.mc_generators):
            gen_offset = MPV_X_OFFSET * (gi - (n_gen - 1) / 2.0)
            for mcb in mc_bins:
                h_mc = self.get_mc_component(comp, gen, mcb, z_cut, den_weight,
                                             ptRL=ptRL)
                if h_mc is None:
                    continue
                x_mc = self._mpv_x_mc(mcb) + gen_offset
                ck = (gen, comp, str(mcb), cut, den_weight, ptRL)
                path = self._mpv_diag_path(gen, comp, cut, den_weight, ptRL,
                                           str(mcb))
                title = (f"{gen_label} | {comp} | jetpt{mcb} | "
                         f"{cut} | ww{den_weight} | {obs}")
                mu, mu_err = self._fit_hist_cached(h_mc, ck, save_path=path,
                                                   title=title)
                if mu is not None:
                    (self.mpv_data.setdefault(case_key, {})
                         .setdefault(gen, {})[x_mc]) = (mu, mu_err)

    def _mpv_styles(self):
        """matplotlib styles for the three sources on the summary plot."""
        return {
            "data":   dict(color="black", marker="o", linestyle="none",
                           label="data"),
            "pythia": dict(color="C3", marker="s", linestyle="-",
                           label="PYTHIA8"),
            "herwig": dict(color="C0", marker="^", linestyle="--",
                           label="Herwig7"),
        }

    def plot_all_mpv_summaries(self):
        """Produce peak-position vs jet-pT summaries (data vs PYTHIA vs Herwig)
        for every collected (component, cut, weight, ptRL) case."""
        styles = self._mpv_styles()

        for case_key in sorted(self.mpv_data.keys()):
            data = self.mpv_data[case_key]
            if not data:
                continue

            # parse: "{comp}_{cut}_ww{den}[_ptRL]"
            ptRL = case_key.endswith("_ptRL")
            core = case_key[:-5] if ptRL else case_key
            # core = "{comp}_{cut}_ww{den}"
            comp, cut_suffix, ww_token = core.split("_", 2)
            den_weight = ww_token[2:]          # strip leading "ww"
            obs = "p_{T}R_{L}" if ptRL else "R_{L}"

            title = (f"MPV summary — {comp}, data vs PYTHIA vs Herwig, "
                     f"{cut_suffix}, ww{den_weight}"
                     + (" (ptRL)" if ptRL else ""))
            outdir = os.path.join(self.base_plot_dir, cut_suffix,
                                  "mpv_summary", comp,
                                  "ptRL" if ptRL else "RL")
            fname = (f"mpv_summary_{comp}{'_ptRL' if ptRL else ''}"
                     f"_{cut_suffix}_ww{den_weight}.pdf")

            mpv.plot_mpv_summary(
                data, title, "Jet p_T [GeV/c]",
                f"Peak position in {obs}",
                outdir, fname,
                component_styles=styles,
            )

        print(f"MPV summary plots saved under: "
              f"{os.path.join(self.base_plot_dir, self.get_cut_suffix(self.z_cut), 'mpv_summary')}")

    def plot_mpv(self):
        """Standalone entry point: open files, collect MPVs, make summaries.

        Only inclusive MC is used (self.mc_parton == 'inclusive'), and data is
        inherently inclusive here, so all summaries are the inclusive case.
        """
        self.open_files()
        self._init_mpv_storage()

        components = ["rad", "AA", "BB", "AB"] #, "CAB"]

        for data_jetpt, mc_jetpt in self.pt_pairs:
            for cut_mode, z_cut in self.cut_modes:
                self.cut_mode = cut_mode
                self.z_cut = z_cut
                for den_weight in self.den_weights:
                    for comp in components:
                        for ptRL in (False, True):
                            self.collect_mpv(comp, data_jetpt, mc_jetpt,
                                             z_cut, den_weight, ptRL=ptRL)

        self.plot_all_mpv_summaries()
        print("\nDone with MPV summaries.")
        
    # =====================================================================
    # Main loop
    # =====================================================================

    def plot(self):
        self.open_files()

        components = ["rad", "AA", "BB", "AB", "CAB"]

        for data_jetpt, mc_jetpt in self.pt_pairs:
            for cut_mode, z_cut in self.cut_modes:
                self.cut_mode = cut_mode
                self.z_cut = z_cut
                print(f"Processing {cut_mode} z_cut={z_cut} | "
                      f"data {self._data_ptrange_label(data_jetpt)} "
                      f"<-> MC {mc_jetpt} ...")
                for den_weight in self.den_weights:
                    # individual components (with ratio panels)
                    for comp in components:
                        for ptRL in (False, True):
                            self.plot_component(comp, data_jetpt, mc_jetpt,
                                                z_cut, den_weight, ptRL=ptRL)
                    # combined 2x2 quadrant canvas (rad/AA/BB/AB)
                    for ptRL in (False, True):
                        self.plot_quadrants(data_jetpt, mc_jetpt,
                                            z_cut, den_weight, ptRL=ptRL)

        print(f"\nDone. Plots saved under: {self.base_plot_dir}")

if __name__ == "__main__":
    plotter = PlotDataMCComparison()
    plotter.plot()
    plotter.plot_mpv()    # new MPV summaries