# Plot DATA curves for JSE analysis (stage 3)
#
# Single data source, single AnalysisResults.root, jet-pt RANGE bins.

import os
import ROOT
import math
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
    def __init__(self, groomed_binning=False):
        ROOT.gROOT.SetBatch(True)
        ROOT.gStyle.SetLegendBorderSize(0)
        ROOT.gStyle.SetLegendFillColor(0)
        ROOT.gStyle.SetPadGridX(1)
        ROOT.gStyle.SetPadGridY(1)

        ROOT.gStyle.SetOptStat(0)

        self.crosscheck = False #True

        # -------------------------------------------------------------
        # Binning mode.
        #   False -> slices in ungroomed jet pT, names use  'jetpt{lo}_{hi}'
        #   True  -> slices in groomed/radiator pT, names use 'gjetpt{lo}_{hi}'
        # -------------------------------------------------------------
        self.groomed_binning = groomed_binning

        # jet pT RANGE bins (available in file:
        #   [(10,20),(20,40),(40,60),(60,80),(80,100),(100,120),(120,150),(150,200),(50,60)])
        self.target_jet_pts_ungroomed = [
            (10, 20), (20, 40), (40, 60), (60, 80), (80, 100), (100, 120), (120, 150), (150, 200)
            # (60, 80), (80, 100), (100, 120), (120, 150), (150, 200)
        ]
        # groomed pT is strictly below the ungroomed pT of the same jet, so the
        # useful slices sit lower; adjust to whatever the writer actually filled.
        self.target_jet_pts_groomed = [
            (10, 20), (20, 40), (40, 60), (60, 80), (80, 100), (100, 120), (120, 150), (150, 200)
        ]
        self.target_jet_pts = (self.target_jet_pts_groomed if self.groomed_binning
                               else self.target_jet_pts_ungroomed)

        # self.cut_modes = [("sd", 0.1), ("maxkt", None)]
        self.cut_modes = [("sd", 0.1)] #, ("sd", 0.2)]

        print("self.groomed_binning", self.groomed_binning)

        # perlmutter
        if self.groomed_binning:
            print("Using groomed binning file")
            self.rootfile_path = ("/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data/57540479/AnalysisResultsMerged_groomedbins.root") # Groomed, hiccup
            self.zcut2_rootfile_path = ("/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data/57174964/AnalysisResultsMerged_groomedbins.root")
        else:
            print("Using ungroomed binning file")
            self.rootfile_path = ("/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data/57174964/AnalysisResultsMerged_ungroomedbins.root") # Ungroomed, perlmutter
            self.zcut2_rootfile_path = ("/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data/57174964/AnalysisResultsMerged_ungroomedbins.root") # small bins
        

        # self.rootfile_path = ("/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data/55778272/AnalysisResultsMerged.root")
        # # self.rootfile_path = ("/global/cfs/cdirs/alice/blianggi/mypyjetty/analysis/testing/AnalysisResults.root")
        # self.zcut2_rootfile_path = ("/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data/56140697/AnalysisResultsMerged.root")
        
        # perlmutter smaller bins
        # self.rootfile_path = ("/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data/56300667/AnalysisResultsMerged.root") # perlmutter, smaller bins
        # self.zcut2_rootfile_path = ("/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data/56300667/AnalysisResultsMerged.root") # perlmutter, smaller bins
        
        # hiccup
        # if self.groomed_binning:
        #     print("Using groomed binning file")
        #     self.rootfile_path = ("/rstorage/alice/AnalysisResults/blianggi/jse/data/1839690/AnalysisResultsMerged.root") # Groomed, hiccup
        #     self.zcut2_rootfile_path = ("/rstorage/alice/AnalysisResults/blianggi/jse/data/1839690/AnalysisResultsMerged.root")
        # else:
        #     print("Using ungroomed binning file")
        #     self.rootfile_path = ("/rstorage/alice/AnalysisResults/blianggi/jse/data/56300667/AnalysisResultsMerged.root") # Ungroomed, hiccup
        #     self.zcut2_rootfile_path = ("/rstorage/alice/AnalysisResults/blianggi/jse/data/56300667/AnalysisResultsMerged.root")
        self.data_rootfile = None
        self.current_rootfile_path = None

        self.cut_mode = ""
        self.z_cut = None
        self.base_plot_dir = ("/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots")
        # self.base_plot_dir = ("/software/users/blianggi/mypyjetty/storage/jse/plots")

        self.outputforvassu = "/global/cfs/cdirs/alice/blianggi/mypyjetty/jse/InclusivePassedSDEECs.root"
        self.outputforvassufile = ROOT.TFile.Open(self.outputforvassu, "RECREATE")

        self._persistent_canvases = []
        self._canvas_counter = 0

        self._counter_cache = {}   # (pr, cut) -> (num_jets, avg_jet_pt)

        self.WEIGHT_TOKEN = {"jet": "wwjetpt", "rad": "wwradpt", "null": "wwnullpt"}

        self._keep = []          # anything that must outlive the drawing call

    # -------------------------------------------------------------------------
    # pt-range helpers
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # pt-range helpers
    # -------------------------------------------------------------------------

    @property
    def binning_tag(self):
        """Slice-token prefix: 'gjetpt' for groomed binning, 'jetpt' otherwise."""
        return "gjetpt" if self.groomed_binning else "jetpt"

    @property
    def binning_dir(self):
        """Subdirectory so the two binnings never overwrite each other."""
        return "groomed_bins" if self.groomed_binning else "ungroomed_bins"

    @property
    def binning_desc(self):
        """How to describe the binning variable in legends."""
        return "gr. jet p_{T}" if self.groomed_binning else "jet p_{T}"

    def _ptrange_token(self, jetpt):
        """jetpt is a (lo, hi) tuple -> 'jetpt70_80' or 'gjetpt70_80'."""
        lo, hi = jetpt
        return f"{self.binning_tag}{lo}_{hi}"

    @staticmethod
    def _ptrange_label(jetpt):
        lo, hi = jetpt
        return f"{lo}-{hi}"

    def _ptrange_tag(self, jetpt):
        """For filenames; same token as the histogram names."""
        return self._ptrange_token(jetpt)

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
    # Normalization factors (read-only)
    # -------------------------------------------------------------------------

    def _get_norm_factors(self, jetpt, z_cut):
        """Return (num_jets_passed_cut, avg_jet_pt) for this slice, cached.

        #counters bin 1 = num_jets
        counters bin 2 = num_jets_passed_cut
        counters bin 3 = sum_jetpt_passed_cut
        avg_jet_pt = sum_jetpt_passed_cut / num_jets_passed_cut  (correct global average)
        Returns (None, None) if counters missing or num_jets <= 0.
        """
        cut = self.get_cut_suffix(z_cut)
        pr = self._ptrange_token(jetpt)
        key = (pr, cut)
        if key in self._counter_cache:
            return self._counter_cache[key]

        counters = self._get(f"counters_{pr}_{cut}")
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

    def _normalize_eec(self, hist, jetpt, z_cut, ptRL=False, save=False):
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

        # hist.Scale(1.0 / num_jets, "width")

        oldname = hist.GetName()
        self.outputforvassufile.cd()
        hist.Scale(1.0, "width")
        hist.SetName(f"{oldname}_justdividedbybinwidth")
        if "full" in oldname and not ptRL and save:
            hist.Write()
            print("wrote!!")
        if ptRL:
            print("histogram maximum", hist.GetName(), hist.GetMaximum())
        hist.Scale(1.0 / num_jets)
        hist.SetName(f"{oldname}")
        if "full" in oldname and not ptRL and save:
            hist.Write()
        return hist


    # -------------------------------------------------------------------------
    # File handling (single file, opened once)
    # -------------------------------------------------------------------------

    def open_rootfile(self, z_cut=None):
        # Determine which file we should be using
        target_path = self.zcut2_rootfile_path if z_cut == 0.2 else self.rootfile_path

        # If already open and it's the correct file, just return
        if self.data_rootfile and not self.data_rootfile.IsZombie():
            if self.current_rootfile_path == target_path:
                return
            # Otherwise, we need to switch files
            self.data_rootfile.Close()

        # Open the target file
        self.data_rootfile = ROOT.TFile.Open(target_path)
        if not self.data_rootfile or self.data_rootfile.IsZombie():
            raise RuntimeError(f"Could not open {target_path}")

        self.current_rootfile_path = target_path

    # -------------------------------------------------------------------------
    # Output paths and canvases
    # -------------------------------------------------------------------------

    def get_output_dir(self, subdir, z_cut, plot_type=None):
        parts = [self.base_plot_dir, subdir, self.binning_dir,
                 self.get_cut_suffix(z_cut)]
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
        self._keep.append(line)
        return line

    def draw_rl_low_region(self, hist=None, x_max=0.01, color=ROOT.kGray + 1,
                        alpha=0.35, redraw=None):
        """Shade the region R_L < x_max on the current pad."""
        pad = ROOT.gPad
        if not pad:
            return None

        # frame coordinates are only valid once the pad has been laid out
        pad.Modified()
        pad.Update()

        x_lo, x_hi = pad.GetUxmin(), pad.GetUxmax()
        y_lo, y_hi = pad.GetUymin(), pad.GetUymax()
        if pad.GetLogx():
            x_lo, x_hi = 10.0 ** x_lo, 10.0 ** x_hi
        if pad.GetLogy():
            y_lo, y_hi = 10.0 ** y_lo, 10.0 ** y_hi

        x_right = min(x_max, x_hi)
        if not np.isfinite(x_lo) or not np.isfinite(x_right) or x_right <= x_lo:
            return None
        if not (np.isfinite(y_lo) and np.isfinite(y_hi)) or y_hi <= y_lo:
            return None

        box = ROOT.TBox(x_lo, y_lo, x_right, y_hi)
        box.SetFillColorAlpha(color, alpha)
        box.SetLineColor(color)
        box.SetLineWidth(0)
        box.Draw()

        # optionally put the curves back on top of the band
        if redraw:
            for h, opt in redraw:
                h.Draw(opt + " SAME")

        pad.RedrawAxis()
        pad.Modified()
        pad.Update()
        self._keep.append(box)      # <-- the actual fix
        return box

    def MakeEventLeg(self, jetpt_label, z_cut, den_weight="",
                     x1=0.15, y1=0.7, x2=0.40, y2=0.88):
        leg = ROOT.TLegend(x1, y1, x2, y2)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetMargin(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.04)
        leg.AddEntry(ROOT.nullptr, "ALICE WIP", "")
        leg.AddEntry(ROOT.nullptr, "pp data, R = 0.4 jets", "")
        leg.AddEntry(ROOT.nullptr, f"{self.binning_desc} = {jetpt_label} GeV/c", "")
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

    def GetSubjetEECHists(self, jetpt, z_cut, den_weight, ptRL=False, norm_by_jets=True, save=False):
        cut = self.get_cut_suffix(z_cut)
        pr = self._ptrange_token(jetpt)
        w = self.WEIGHT_TOKEN[den_weight]

        # hist_full has NO weight token; has a ptRL variant
        full_name = f"{self._ptRL('hist_full', ptRL)}{pr}_{cut}"
        rad_name  = f"{self._ptRL('hist_rad', ptRL)}{pr}_{cut}_{w}"
        aa_name   = f"{self._ptRL('hist_AA', ptRL)}{pr}_{cut}_{w}"
        bb_name   = f"{self._ptRL('hist_BB', ptRL)}{pr}_{cut}_{w}"
        ab_name   = f"{self._ptRL('hist_AB', ptRL)}{pr}_{cut}_{w}"

        hists = (self._get(full_name), self._get(rad_name), self._get(aa_name),
                 self._get(bb_name), self._get(ab_name))
        if norm_by_jets:
            print("about to normalize function", jetpt, z_cut, den_weight, ptRL, norm_by_jets)
            if den_weight == "rad":
                save = False #only save when looking at full eec
            final_hists = tuple(self._normalize_eec(h, jetpt, z_cut, ptRL=ptRL, save=save) for h in hists)
        return final_hists if norm_by_jets else hists

    # Can't use this anymore because CAB was calculated with unnormalized histograms... though maybe this shouldn't make a difference?
    # def GetCABHist(self, jetpt, z_cut, den_weight, ptRL=False):
        # cut = self.get_cut_suffix(z_cut)
        # pr = self._ptrange_token(jetpt)
        # w = self.WEIGHT_TOKEN[den_weight]
        # name = f"{self._ptRL('CAB', ptRL)}{pr}_{cut}_{w}"
    #     return self._get(name)

    def GetCABHist(self, jetpt, z_cut, den_weight, hist_AA, hist_BB, hist_AB, ptRL=False):
        cut = self.get_cut_suffix(z_cut)
        pr = self._ptrange_token(jetpt)
        w = self.WEIGHT_TOKEN[den_weight]
        name = f"{self._ptRL('CAB', ptRL)}{pr}_{cut}_{w}"
        print("name", name)

        hist_CAB = hist_AB.Clone(name)
        hist_CAB.SetTitle("C_{AB};R_{L};AxB / #sqrt{AxA #times BxB}")
        hist_CAB.Reset()
        for i in range(1, hist_CAB.GetNbinsX() + 1):
            aa, bb, ab = hist_AA.GetBinContent(i), hist_BB.GetBinContent(i), hist_AB.GetBinContent(i)
            sig_aa, sig_bb, sig_ab = hist_AA.GetBinError(i), hist_BB.GetBinError(i), hist_AB.GetBinError(i)
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
        draw_opt = "EP" if first else "EP SAME" #"HIST" if first else "HIST SAME"
        for h in hists:
            h.Draw(draw_opt)
            draw_opt = "EP SAME" #"HIST SAME"

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
        self.draw_rl_low_region(hists_primary[0]) 
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
            print("subjeteec1")
            hists = self.GetSubjetEECHists(jetpt, z_cut, den_weight, ptRL=ptRL, save=True)
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

            ev_leg = self.MakeEventLeg(label, z_cut, den_weight, x1=0.15, y1=0.65, x2=0.40, y2=0.88)
            legend = ROOT.TLegend(0.68, 0.65, 0.88, 0.88)
            if hist_full is not None and den_weight == "jet":
                legend.AddEntry(hist_full,
                                f"all jets that passed {self.get_passed_label()}", "l")
            legend.AddEntry(hist_rad, "radiator", "l")
            if crosscheck:
                legend.AddEntry(crosscheck, "cross-check #Sigma_{subjets} = rad.", "l")
            legend.AddEntry(hist_AA, "AxA", "l")
            legend.AddEntry(hist_BB, "BxB", "l")
            legend.AddEntry(hist_AB, "AxB", "l")

            excluded = {None}
            if den_weight == "rad":
                excluded.add(hist_full)
            hists_to_draw = [h for h in hists if h not in excluded]
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

        print("subjeteec2")
        _, _, hist_AA, hist_BB, hist_AB = self.GetSubjetEECHists(
            jetpt, z_cut, den_weight, False, False)
        if not all([hist_AA, hist_BB, hist_AB]):
            print(f"WARNING: AA/BB/AB not found for data {tag} {cut}")
            return

        # # ----- Canvas 1: R_g vs all EEC components -----
        # canvas = self.make_canvas("canvas_rg", tag=tag, z_cut=z_cut,
        #                           den_weight=den_weight,
        #                           title="R_{g} vs EEC components")
        # canvas.SetLogx()
        # canvas.cd()

        # hist_AA.SetLineColor(Color.BLUE)
        # hist_BB.SetLineColor(Color.ORANGE)
        # hist_AB.SetLineColor(Color.GREEN)
        # hist_rg.SetLineColor(ROOT.kRed + 1)
        # hist_rg.SetLineStyle(ROOT.kDashed)

        # y_max = max(hist_AA.GetMaximum(), hist_BB.GetMaximum(),
        #             hist_AB.GetMaximum(), hist_rg.GetMaximum()) * 1.3
        # hist_AA.SetMaximum(y_max)
        # hist_AA.SetMinimum(0)
        # hist_AA.GetXaxis().SetTitle("R_{L} or R_{g}")
        # hist_AA.GetYaxis().SetTitle("(1/N_{jets}) dN/d(R_{L} or R_{g})")

        # hist_AA.Draw("HIST")
        # hist_BB.Draw("HIST SAME")
        # hist_AB.Draw("HIST SAME")
        # hist_rg.Draw("HIST SAME")

        # ev_leg = self.MakeEventLeg(label, z_cut, den_weight)
        # legend = ROOT.TLegend(0.55, 0.65, 0.88, 0.88)
        # legend.AddEntry(hist_AA, "AxA (EEC)", "l")
        # legend.AddEntry(hist_BB, "BxB (EEC)", "l")
        # legend.AddEntry(hist_AB, "AxB (EEC)", "l")
        # legend.AddEntry(hist_rg, "R_{g} = #DeltaR_{AB}/R", "l")
        # ev_leg.Draw()
        # legend.Draw()

        # output_name = (f"rg_vs_eec_data_{tag}_R0.4_{cut}_ww{den_weight}.pdf")
        # self._save_canvas(canvas, "data", z_cut, output_name, plot_type="rg")

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
        self.draw_rl_low_region(hist_AB_norm)

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
        self.draw_rl_low_region(hist_ratio)
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

        print("subjeteec34")
        _, _, hist_AA_jet, hist_BB_jet, hist_AB_jet = self.GetSubjetEECHists(jetpt, z_cut, "jet", ptRL=ptRL)
        _, _, hist_AA_rad, hist_BB_rad, hist_AB_rad = self.GetSubjetEECHists(jetpt, z_cut, "rad", ptRL=ptRL)
        hist_wwjetpt = self.GetCABHist(jetpt, z_cut, "jet", hist_AA_jet, hist_BB_jet, hist_AB_jet, ptRL=ptRL)
        hist_wwradpt = self.GetCABHist(jetpt, z_cut, "rad", hist_AA_rad, hist_BB_rad, hist_AB_rad, ptRL=ptRL)
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
        self.draw_rl_low_region(hist_wwjetpt)
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
        # Expanded colour/marker cycles
        base_colors = [
            ROOT.kBlack, ROOT.kBlue, ROOT.kOrange + 7, ROOT.kGreen + 2,
            ROOT.kRed + 1, ROOT.kMagenta + 1, ROOT.kCyan + 2, ROOT.kAzure + 1,
            ROOT.kYellow + 1, ROOT.kGray + 1
        ]
        base_markers = [
            ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kFullDiamond,
            ROOT.kFullStar, ROOT.kFullTriangleUp, ROOT.kFullTriangleDown,
            ROOT.kFullCross, ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kFullDiamond
        ]
        colors  = [base_colors[i % len(base_colors)] for i in range(n_pt)]
        markers = [base_markers[i % len(base_markers)] for i in range(n_pt)]

        panel_labels = ["radiator", "AxA", "BxB", "AxB"]

        canvas = self.make_canvas("canvas_acrosspt", z_cut=z_cut,
                                  den_weight=den_weight, w=1200, h=1000)
        canvas.Divide(2, 2, 0.005, 0.005)

        pt_lo = min(lo for lo, _ in self.target_jet_pts)
        pt_hi = max(hi for _, hi in self.target_jet_pts)
        ev_leg = self.MakeEventLeg(f"{pt_lo}-{pt_hi}", z_cut, den_weight)
        legend_pt = ROOT.TLegend(0.6, 0.6, 0.88, 0.88) # Restored vertical dimensions
        legend_pt.SetBorderSize(0)

        dummy_graphs = []
        for jetpt, marker, color in zip(self.target_jet_pts, markers, colors):
            d = ROOT.TGraph(1)
            d.SetMarkerStyle(marker)
            d.SetMarkerColor(color)
            d.SetLineColor(color)
            # legend_pt.AddEntry(d, f"jet p_{{T}} = {self._ptrange_label(jetpt)} GeV/c", "p")
            legend_pt.AddEntry(d, f"{self.binning_desc} = {self._ptrange_label(jetpt)} GeV/c", "p")
            dummy_graphs.append(d)

        # persistent[ijetpt] = [rad, AA, BB, AB] or None
        persistent = []
        for ijetpt, jetpt in enumerate(self.target_jet_pts):
            print("subjeteec5")
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
                    h.GetYaxis().SetTitle("(1/N_{jets}) dN/dR_{L}")
                    h.GetYaxis().SetTitleSize(0.06)
                    h.GetYaxis().SetTitleOffset(0.9)
                    h.GetXaxis().SetTitleSize(0.06)
                    h.GetXaxis().SetTitle("R_{L}")
                h.Draw("PE" if first else "PE SAME")
                self.draw_rl_low_region(h)
                first = False

            latex = ROOT.TLatex()
            latex.SetNDC()
            latex.SetTextSize(0.07)
            latex.DrawLatex(0.18, 0.84, label)
            latex_labels.append(latex)

            if ipanel == 0:
                ev_leg.Draw()
                legend_pt.Draw()

        output_name = (f"subjet_eec_acrosspt_data_all{self.binning_tag}"f"_R0.4_{cut}_ww{den_weight}.pdf")
        self._save_canvas(canvas, "data", z_cut, output_name, "subjet_eec")

        # ---- C_AB overlaid across pT ----
        can_CAB = self.make_canvas("can_CAB_acrosspt", z_cut=z_cut,
                                   den_weight=den_weight, title="C_{AB}")
        can_CAB.SetLogx()
        can_CAB.SetLogy()
        can_CAB.cd()
        legend_CAB = ROOT.TLegend(0.65, 0.12, 0.88, 0.32) # Moved to bottom right
        legend_CAB.SetBorderSize(0)
        persistent_CAB = []
        first_CAB = True
        for ijetpt, jetpt in enumerate(self.target_jet_pts):
            print("subjeteec6")
            hists = self.GetSubjetEECHists(jetpt, z_cut, den_weight)
            _, hr, hAA, hBB, hAB = hists
            h = self.GetCABHist(jetpt, z_cut, den_weight, hAA, hBB, hAB)
            if not h:
                continue
            tag = self._ptrange_tag(jetpt)
            c = h.Clone(f"acrosspt_CAB_data_{tag}_{cut}_ww{den_weight}")
            c.SetDirectory(0)
            self.FormatHist(c, colors[ijetpt], ijetpt + 1, markers[ijetpt])

            # Extract slope for RL = 0.2 to 1
            slope_str = ""
            try:
                # Extract points
                bins = h.GetNbinsX()
                xs, ys = [], []
                for b in range(1, bins + 1):
                    xb = h.GetBinCenter(b)
                    if 0.2 <= xb <= 1.0:
                        xs.append(xb)
                        ys.append(h.GetBinContent(b))

                if len(xs) >= 2:
                    # Fit straight line y = mx + c
                    coeffs = np.polyfit(xs, ys, 1)
                    slope = coeffs[0]
                    slope_str = f" (slope={slope:.2e})"
            except Exception as e:
                print(f"Fit failed for {tag}: {e}")

            # legend_CAB.AddEntry(c, f"jet p_{{T}} = {self._ptrange_label(jetpt)} GeV/c{slope_str}", "pe")
            legend_CAB.AddEntry(d, f"{self.binning_desc} = {self._ptrange_label(jetpt)} GeV/c", "p") #TODO: fix?
            persistent_CAB.append(c)
            c.Draw("PE" if first_CAB else "PE SAME")
            self.draw_rl_low_region(c)
            first_CAB = False

        if persistent_CAB:
            ev_leg.Draw()
            legend_CAB.Draw()
            self.draw_hori_line(1e-3, 1, 1, ROOT.kGray + 3, 9)

        output_name = (f"CAB_acrosspt_data_all{self.binning_tag}"f"_R0.4_{cut}_ww{den_weight}.pdf")
        self._save_canvas(can_CAB, "data", z_cut, output_name, plot_type="CAB")

        # ---- C_AB overlaid across pT (with Fit Lines) ----
        can_CAB_fits = self.make_canvas("can_CAB_acrosspt_fits", z_cut=z_cut,
                                        den_weight=den_weight, title="C_{AB} fits")
        can_CAB_fits.SetLogx()
        can_CAB_fits.SetLogy()
        can_CAB_fits.cd()
        legend_CAB_fits = ROOT.TLegend(0.65, 0.12, 0.88, 0.32)
        legend_CAB_fits.SetBorderSize(0)
        persistent_CAB_fits = []
        first_CAB_fits = True
        for ijetpt, jetpt in enumerate(self.target_jet_pts):
            print("subjeteec7")
            hists = self.GetSubjetEECHists(jetpt, z_cut, den_weight)
            _, hr, hAA, hBB, hAB = hists
            h = self.GetCABHist(jetpt, z_cut, den_weight, hAA, hBB, hAB)
            if not h:
                continue
            tag = self._ptrange_tag(jetpt)
            c_fits = h.Clone(f"acrosspt_CAB_fits_data_{tag}_{cut}_ww{den_weight}")
            c_fits.SetDirectory(0)
            self.FormatHist(c_fits, colors[ijetpt], ijetpt + 1, markers[ijetpt])

            slope_str_fits = ""
            try:
                bins = h.GetNbinsX()
                xs, ys = [], []
                for b in range(1, bins + 1):
                    xb = h.GetBinCenter(b)
                    yb = h.GetBinContent(b)
                    if 0.2 <= xb <= 1.0 and yb > 0:
                        xs.append(xb)
                        ys.append(yb)

                if len(xs) >= 2:
                    # Fit log(y) = m*log(x) + b  =>  y = exp(b) * x^m
                    # This looks linear in log(x) and log(y)
                    coeffs = np.polyfit(np.log(xs), np.log(ys), 1)
                    m, b_int = coeffs[0], coeffs[1]
                    A = np.exp(b_int)
                    slope_str_fits = f" (slope={m:.3f})"

                    # Draw power-law fit line: y = [0] * x^[1]
                    fit_func = ROOT.TF1(f"fit_{tag}_{ijetpt}", "[0]*TMath.Power(x, [1])", 0.2, 1.0)
                    fit_func.SetParameters(A, m)
                    fit_func.SetLineColor(colors[ijetpt])
                    fit_func.SetLineWidth(2)
                    fit_func.Draw("SAME")
            except Exception as e:
                print(f"Fit failed for {tag} (fits version): {e}")

            # legend_CAB_fits.AddEntry(c_fits, f"jet p_{{T}} = {self._ptrange_label(jetpt)} GeV/c{slope_str_fits}", "pe")
            legend_CAB_fits.AddEntry(d, f"{self.binning_desc} = {self._ptrange_label(jetpt)} GeV/c", "p") #TODO: fix?persistent_CAB_fits.append(c_fits)
            c_fits.Draw("PE" if first_CAB_fits else "PE SAME")
            self.draw_rl_low_region(c_fits)
            first_CAB_fits = False

        if persistent_CAB_fits:
            ev_leg.Draw()
            legend_CAB_fits.Draw()
            self.draw_hori_line(1e-3, 1, 1, ROOT.kGray + 3, 9)

        output_name_fits = (f"CAB_acrosspt_fits_data_all{self.binning_tag}"f"_R0.4_{cut}_ww{den_weight}.pdf")
        self._save_canvas(can_CAB_fits, "data", z_cut, output_name_fits, plot_type="CAB")

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
        self.mpv_meta = {}   # case_key -> dict(cut=, den=, ptRL=, binning=)

    def _mpv_case_key(self, cut, den_weight, ptRL):
        return (f"basic_{self.binning_tag}_{cut}_ww{den_weight}"
                + ("_ptRL" if ptRL else ""))

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
        for ptRL in (False, True):
            print("subjeteec8")
            hists = self.GetSubjetEECHists(jetpt, z_cut, den_weight, ptRL=ptRL)
            case_key = self._mpv_case_key(cut, den_weight, ptRL)
            self.mpv_meta[case_key] = dict(cut=cut, den=den_weight, ptRL=ptRL,
                                           binning=self.binning_tag,
                                           groomed=self.groomed_binning)
            diag_key = (f"phys_data_{self.binning_tag}_{cut}_ww{den_weight}"
                        + ("_ptRL" if ptRL else ""))
            for comp_name, h in zip(["full", "rad", "AA", "BB", "AB"], hists):
                if h is None:
                    continue
                save_path = self._fit_diag_path(diag_key, comp_name, tag)
                title = (f"{'basic_ptRL' if ptRL else 'basic'} | data | "
                         f"{self.binning_desc}={self._ptrange_label(jetpt)} | "
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
            meta = self.mpv_meta.get(case_key)
            if meta is None:
                print(f"WARNING: no metadata for MPV case {case_key}, skipping")
                continue
            cut_suffix = meta["cut"]
            den_token  = f"ww{meta['den']}"
            ptRL       = meta["ptRL"]
            binning    = meta["binning"]
            obs = "p_TR_L" if ptRL else "R_L"
            ylabel = f"Peak position in {obs}"
            xlabel = ("Groomed jet p_T [GeV/c]" if meta["groomed"]
                      else "Jet p_T [GeV/c]")
            title = (f"MPV summary — data, {binning}, {cut_suffix}, {den_token}"
                     + (" (ptRL)" if ptRL else ""))
            outdir = os.path.join(base, "data", self.binning_dir, cut_suffix,
                                  "ptRL" if ptRL else "RL")
            fname = (f"mpv_basic{'_ptRL' if ptRL else ''}_data_{binning}"
                     f"_{cut_suffix}_{den_token}.pdf")
            self._plot_mpv_summary(case_key, title, xlabel, ylabel,
                                   outdir, fname, component_styles=basic_styles)
        print(f"MPV summary plots saved under: {base}")

    # -------------------------------------------------------------------------
    # Main loop
    # -------------------------------------------------------------------------

    def plot(self):
        # self.normalize_all()
        self._init_mpv_storage()

        for i, jetpt in enumerate(self.target_jet_pts):
            for cut_mode, z_cut in self.cut_modes:
                self.cut_mode = cut_mode
                self.z_cut = z_cut
                self.open_rootfile(z_cut)
                print(f"Processing {cut_mode} mode, z_cut={z_cut}, "
                      f"jetpt={self._ptrange_label(jetpt)}...")

                for den_weight in ["jet", "rad"]: #, "null"]:
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
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--binning", choices=["ungroomed", "groomed", "both"],
                    default="groomed",
                    help="which pT slicing to read: 'jetpt...' or 'gjetpt...' hists")
    args = ap.parse_args()

    modes = ([False, True] if args.binning == "both"
             else [args.binning == "groomed"])
    for groomed in modes:
        print(f"\n=== binning: {'groomed' if groomed else 'ungroomed'} ===")
        PlotDataCurves(groomed_binning=groomed).plot()