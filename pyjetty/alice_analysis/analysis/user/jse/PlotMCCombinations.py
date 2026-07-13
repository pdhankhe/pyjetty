# Plot multiplicity / counting curves for JSE analysis
# - N particles in subjet A vs B  (2D; ProjectionX/Y overlaid, plus 2D comparisons)
# - N ungroomed / N groomed
# - combAA, combBB, combAB, combTotal
# Comparisons: PYTHIA vs HERWIG, across pt bins, parton types, and cut modes.

import os
import ROOT
from enum import IntEnum


class Color(IntEnum):
    BLUE   = ROOT.TColor.GetColor("#1f77b4")  # matplotlib C0
    GREEN  = ROOT.TColor.GetColor("#2ca02c")  # matplotlib C2
    RED    = ROOT.TColor.GetColor("#d62728")  # matplotlib C3
    ORANGE = ROOT.TColor.GetColor("#ff7f0e")  # matplotlib C1


class PlotCounting:
    def __init__(self):
        ROOT.gROOT.SetBatch(True)          # must be first — headless nodes
        ROOT.gStyle.SetLegendBorderSize(0)
        ROOT.gStyle.SetLegendFillColor(0)
        ROOT.gStyle.SetPadGridX(1)
        ROOT.gStyle.SetPadGridY(1)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPalette(ROOT.kBird)

        # ---- generator toggle ----
        self.include_pythia = True   # flip to True once pythia is up to date
        self.generators = (["pythia", "herwig"] if self.include_pythia
                           else ["herwig"])

        self.target_jet_pts = [50, 100, 200, 500]
        self.partontypes    = ["inclusive", "quark", "gluon"]
        self.cut_modes      = [("sd", 0.1)] #, ("maxkt", None)]

        self.input_base = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/rootfiles" # perlmutter
        # self.input_base = "/software/users/blianggi/mypyjetty/storage/jse/rootfiles" # hiccup
        self.rootfile_template = f"{self.input_base}/jse_preliminary_curves_{{gen}}_jetpt{{jetpt}}.root"
        self.current_jetpt   = None
        self.pythia_rootfile = None
        self.herwig_rootfile = None

        self.cut_mode = ""
        # self.base_plot_dir = "/software/users/blianggi/mypyjetty/storage/jse/plots" # hiccup
        self.base_plot_dir = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots" # perlmutter



        # Histogram name templates:  hist_<NAME>_<partontype>_jetpt<jetpt>_<cut_suffix>
        self.NA_NB_TEMPLATE     = "hist_nA_nB_{}_jetpt{}_{}"          # 2D
        self.NUNGROOMED_TEMPLATE = "hist_nTotalUnGroomed_{}_jetpt{}_{}"
        self.NGROOMED_TEMPLATE   = "hist_nTotalGroomed_{}_jetpt{}_{}"
        self.COMBAA_TEMPLATE    = "hist_combAA_{}_jetpt{}_{}"
        self.COMBBB_TEMPLATE    = "hist_combBB_{}_jetpt{}_{}"
        self.COMBAB_TEMPLATE    = "hist_combAB_{}_jetpt{}_{}"
        self.COMBTOTAL_TEMPLATE = "hist_combTotal_{}_jetpt{}_{}"

        self._persistent_canvases = []
        self._canvas_counter      = 0

    # -------------------------------------------------------------------------
    # File / canvas / formatting helpers  (carried over from old code)
    # -------------------------------------------------------------------------

    def get_file(self, gen):
        return self.pythia_rootfile if gen == "pythia" else self.herwig_rootfile

    def set_current_rootfiles(self, jetpt):
        if self.current_jetpt == jetpt:
            return
        for attr in ("pythia_rootfile", "herwig_rootfile"):
            f = getattr(self, attr, None)
            if f:
                try:
                    f.Close()
                except Exception:
                    pass
        # only open pythia if it's enabled
        if self.include_pythia:
            self.pythia_rootfile = ROOT.TFile.Open(
                self.rootfile_template.format(gen="pythia", jetpt=jetpt))
        else:
            self.pythia_rootfile = None
        self.herwig_rootfile = ROOT.TFile.Open(
            self.rootfile_template.format(gen="herwig", jetpt=jetpt))
        self.current_jetpt = jetpt

    def get_cut_suffix(self, z_cut):
        return f"sd{z_cut}" if self.cut_mode == "sd" else "maxkt"

    def get_cut_label(self, z_cut):
        return f"SD z_{{cut}} = {z_cut}" if self.cut_mode == "sd" else "maxkt selection"

    def get_output_dir(self, subdir, z_cut, plot_type=None):
        # Top-level "counting" folder, then subdir / cut_suffix / plot_type
        parts = [self.base_plot_dir, "counting", subdir, self.get_cut_suffix(z_cut)]
        if plot_type:
            parts.extend(plot_type.split('/'))
        path = os.path.join(*parts)
        os.makedirs(path, exist_ok=True)
        return path

    def make_canvas(self, base_name, gen="", partontype="", target_jetpt="",
                    z_cut="", title=None, w=800, h=600):
        self._canvas_counter += 1
        unique_name = (
            f"{base_name}_{gen}_{partontype}_{target_jetpt}"
            f"_{self.get_cut_suffix(z_cut)}_{self._canvas_counter}"
        )
        canvas = ROOT.TCanvas(unique_name, title or base_name, w, h)
        self._persistent_canvases.append(canvas)
        return canvas

    def FormatHist(self, hist, color, linestyle, markerstyle=0, coloralpha=1):
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

    def MakeEventLeg(self, gen, parton_type, jet_pt, z_cut,
                     x1=0.15, y1=0.7, x2=0.40, y2=0.88):
        event_leg = ROOT.TLegend(x1, y1, x2, y2)
        event_leg.SetBorderSize(0)
        event_leg.SetFillColor(0)
        event_leg.SetMargin(0)
        event_leg.AddEntry(ROOT.nullptr, f"{gen} {parton_type} jets", "")
        event_leg.AddEntry(ROOT.nullptr, f"jet p_{{T}} = {jet_pt} GeV/c", "")
        event_leg.AddEntry(ROOT.nullptr, self.get_cut_label(z_cut), "")
        return event_leg

    def _save_canvas(self, canvas, subdir, z_cut, filename, plot_type=None):
        canvas.SaveAs(os.path.join(
            self.get_output_dir(subdir, z_cut, plot_type), filename))

    # -------------------------------------------------------------------------
    # Histogram retrieval
    # -------------------------------------------------------------------------

    def _get(self, file, template, partontype, target_jetpt, z_cut):
        if file is None:
            return None
        cut_suffix = self.get_cut_suffix(z_cut)
        h = file.Get(template.format(partontype, target_jetpt, cut_suffix))
        if h:
            h.SetDirectory(0)   # detach so it survives file.Close()
        return h

    def GetNA_NB(self, file, partontype, target_jetpt, z_cut):
        """Return the 2D N_A vs N_B histogram."""
        return self._get(file, self.NA_NB_TEMPLATE, partontype, target_jetpt, z_cut)

    def GetCountingHists(self, file, partontype, target_jetpt, z_cut):
        """Return the 1D counting histograms as a dict."""
        return {
            "ungroomed":  self._get(file, self.NUNGROOMED_TEMPLATE, partontype, target_jetpt, z_cut),
            "groomed":    self._get(file, self.NGROOMED_TEMPLATE,   partontype, target_jetpt, z_cut),
            "combAA":     self._get(file, self.COMBAA_TEMPLATE,     partontype, target_jetpt, z_cut),
            "combBB":     self._get(file, self.COMBBB_TEMPLATE,     partontype, target_jetpt, z_cut),
            "combAB":     self._get(file, self.COMBAB_TEMPLATE,     partontype, target_jetpt, z_cut),
            "combTotal":  self._get(file, self.COMBTOTAL_TEMPLATE,  partontype, target_jetpt, z_cut),
        }
    
    def _rebin_to_match(self, hist, ref):
        """Rebin `hist` to match `ref`'s bin WIDTH.

        Assumes both start at the same low edge with uniform bins, and that
        ref's bin width is an integer multiple of hist's bin width.
        (Here: combAA/AB/BB have width 18 over [-0.5, comb_max-0.5];
         combTotal has width 72 over [-0.5, 4*comb_max-0.5], so factor = 4.)

        Returns a new (rebinned) histogram; the original is untouched.
        """
        if hist is None or ref is None:
            return hist

        w_hist = hist.GetXaxis().GetBinWidth(1)
        w_ref  = ref.GetXaxis().GetBinWidth(1)

        if w_hist <= 0:
            return hist

        factor = round(w_ref / w_hist)
        if factor < 1:
            factor = 1

        if factor == 1:
            return hist  # already matching width, nothing to do

        # Rebin(group) requires the factor to divide the number of bins evenly.
        if hist.GetNbinsX() % factor != 0:
            print(f"[WARN] _rebin_to_match: {hist.GetName()} has "
                  f"{hist.GetNbinsX()} bins, not divisible by factor {factor}; "
                  f"leaving histogram unchanged.")
            return hist

        # Clone first so the caller's original is preserved.
        new_name = f"{hist.GetName()}_rebin{factor}"
        h_clone = hist.Clone(new_name)
        h_clone.Rebin(factor)   # integer Rebin: no edge array, no edge-mismatch warnings
        return h_clone

    # -------------------------------------------------------------------------
    # Generic 1D PYTHIA vs HERWIG overlay (with ratio panel)
    # -------------------------------------------------------------------------

    def _plot_1d_hvp(self, hist_p, hist_h, partontype, target_jetpt, z_cut,
                     xtitle, ytitle, legend_title, file_tag, plot_type,
                     normalize=True, logy=False):
        if hist_p is None or hist_h is None:
            print(f"WARNING: missing {file_tag} for {partontype} "
                  f"jetpt{target_jetpt} {self.get_cut_suffix(z_cut)}")
            return

        cut_suffix = self.get_cut_suffix(z_cut)

        # clone so normalization doesn't mutate the originals
        hp = hist_p.Clone(f"{hist_p.GetName()}_p_clone")
        hh = hist_h.Clone(f"{hist_h.GetName()}_h_clone")
        hp.SetDirectory(0)
        hh.SetDirectory(0)

        if normalize:
            if hp.Integral() > 0:
                hp.Scale(1.0 / hp.Integral())
            if hh.Integral() > 0:
                hh.Scale(1.0 / hh.Integral())
            ytitle = "self-normalized " + ytitle

        self.FormatHist(hp, Color.BLUE, ROOT.kSolid)
        self.FormatHist(hh, Color.RED,  ROOT.kSolid)

        canvas = self.make_canvas("can_count_hvp", partontype=partontype,
                                  target_jetpt=target_jetpt, z_cut=z_cut,
                                  title=legend_title)
        if logy:
            # logy lives on the top pad below
            pass

        # ----- two-panel: top distribution, bottom ratio -----
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0.02)
        if logy:
            pad1.SetLogy()
        pad1.Draw()

        canvas.cd()
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, 0.3)
        pad2.SetTopMargin(0.02)
        pad2.SetBottomMargin(0.35)
        pad2.SetGridy()
        pad2.Draw()

        # ----- top -----
        pad1.cd()
        y_max = max(hp.GetMaximum(), hh.GetMaximum())
        if logy:
            hp.SetMaximum(y_max * 5)
            hp.SetMinimum(max(1e-6, min(
                v for v in [hp.GetMinimum(0), hh.GetMinimum(0)] if v > 0) * 0.5)
                if y_max > 0 else 1e-6)
        else:
            hp.SetMaximum(y_max * 1.4)
            hp.SetMinimum(0)
        hp.GetYaxis().SetTitle(ytitle)
        hp.GetYaxis().SetTitleSize(0.05)
        hp.GetYaxis().SetTitleOffset(0.9)
        hp.GetXaxis().SetLabelSize(0)

        hp.Draw("HIST")
        hh.Draw("HIST SAME")

        ev_leg = self.MakeEventLeg("PYTHIA and HERWIG", partontype,
                                   target_jetpt, z_cut)
        legend = ROOT.TLegend(0.6, 0.7, 0.88, 0.88)
        legend.AddEntry(hp, f"PYTHIA (mean: {hist_p.GetMean():.2f})", "l")
        legend.AddEntry(hh, f"HERWIG (mean: {hist_h.GetMean():.2f})", "l")
        ev_leg.Draw()
        legend.Draw()

        # ----- bottom (ratio pythia/herwig) -----
        pad2.cd()
        ratio = hp.Clone(f"ratio_{file_tag}_{partontype}_{target_jetpt}_{cut_suffix}")
        ratio.SetDirectory(0)
        ratio.Divide(hh)
        ratio.SetLineColor(ROOT.kBlack)
        ratio.SetMarkerStyle(20)
        ratio.SetMarkerSize(0.7)
        ratio.GetYaxis().SetTitle("PYTHIA / HERWIG")
        ratio.GetYaxis().SetNdivisions(505)
        ratio.GetYaxis().SetTitleSize(0.11)
        ratio.GetYaxis().SetTitleOffset(0.4)
        ratio.GetYaxis().SetLabelSize(0.09)
        ratio.GetXaxis().SetTitle(xtitle)
        ratio.GetXaxis().SetTitleSize(0.12)
        ratio.GetXaxis().SetTitleOffset(1.0)
        ratio.GetXaxis().SetLabelSize(0.09)
        ratio.SetMinimum(0.0)
        ratio.SetMaximum(2.0)
        ratio.Draw("EP")

        line = ROOT.TLine(ratio.GetXaxis().GetXmin(), 1.0,
                          ratio.GetXaxis().GetXmax(), 1.0)
        line.SetLineColor(ROOT.kGray + 2)
        line.SetLineStyle(ROOT.kDashed)
        line.Draw("SAME")

        canvas.cd()
        output_name = (f"{file_tag}_PYTHIA_VS_HERWIG_{partontype}"
                       f"_jetpt{target_jetpt}_R0.4_{cut_suffix}.pdf")
        self._save_canvas(canvas, "pythia_vs_herwig", z_cut, output_name,
                          plot_type=plot_type)

    # -------------------------------------------------------------------------
    # Generic 1D single-generator distribution (no ratio panel)
    # -------------------------------------------------------------------------

    def _plot_1d_single(self, hist, gen, partontype, target_jetpt, z_cut,
                        xtitle, ytitle, file_tag, plot_type,
                        normalize=True, logy=False):
        if hist is None:
            print(f"WARNING: missing {file_tag} for {gen} {partontype} "
                  f"jetpt{target_jetpt} {self.get_cut_suffix(z_cut)}")
            return

        cut_suffix = self.get_cut_suffix(z_cut)

        h = hist.Clone(f"{hist.GetName()}_single_clone")
        h.SetDirectory(0)

        if normalize:
            if h.Integral() > 0:
                h.Scale(1.0 / h.Integral())
            ytitle = "self-normalized " + ytitle

        color = Color.BLUE if gen == "pythia" else Color.RED
        self.FormatHist(h, color, ROOT.kSolid)

        canvas = self.make_canvas("can_count_single", gen=gen,
                                  partontype=partontype,
                                  target_jetpt=target_jetpt, z_cut=z_cut,
                                  title=file_tag)
        canvas.cd()
        if logy:
            canvas.SetLogy()

        y_max = h.GetMaximum()
        if logy:
            h.SetMaximum(y_max * 5)
        else:
            h.SetMaximum(y_max * 1.4)
            h.SetMinimum(0)
        h.GetXaxis().SetTitle(xtitle)
        h.GetYaxis().SetTitle(ytitle)
        h.Draw("HIST")

        ev_leg = self.MakeEventLeg(gen.upper(), partontype, target_jetpt, z_cut)
        legend = ROOT.TLegend(0.6, 0.78, 0.88, 0.88)
        legend.AddEntry(h, f"{gen.upper()} (mean: {hist.GetMean():.2f})", "l")
        ev_leg.Draw()
        legend.Draw()

        output_name = (f"{file_tag}_{gen.upper()}_{partontype}"
                       f"_jetpt{target_jetpt}_R0.4_{cut_suffix}.pdf")
        self._save_canvas(canvas, gen, z_cut, output_name, plot_type=plot_type)

    # -------------------------------------------------------------------------
    # Plot: N_A vs N_B projections (overlaid X and Y), PYTHIA vs HERWIG
    # -------------------------------------------------------------------------

    def plot_nA_nB_projections(self, partontype, target_jetpt, z_cut):
        self.set_current_rootfiles(target_jetpt)
        cut_suffix = self.get_cut_suffix(z_cut)

        h2_p = (self.GetNA_NB(self.pythia_rootfile, partontype, target_jetpt, z_cut)
                if self.include_pythia else None)
        h2_h = self.GetNA_NB(self.herwig_rootfile, partontype, target_jetpt, z_cut)

        if h2_h is None:
            print(f"WARNING: missing herwig nA_nB for {partontype} "
                  f"jetpt{target_jetpt} {cut_suffix}")
            return
        if self.include_pythia and h2_p is None:
            print(f"WARNING: missing pythia nA_nB for {partontype} "
                  f"jetpt{target_jetpt} {cut_suffix}")
            return

        # ProjectionX = N in A, ProjectionY = N in B
        projX_h = h2_h.ProjectionX(f"projX_h_{partontype}_{target_jetpt}_{cut_suffix}")
        projY_h = h2_h.ProjectionY(f"projY_h_{partontype}_{target_jetpt}_{cut_suffix}")
        projX_h.SetDirectory(0)
        projY_h.SetDirectory(0)

        projX_p = projY_p = None
        if self.include_pythia:
            projX_p = h2_p.ProjectionX(f"projX_p_{partontype}_{target_jetpt}_{cut_suffix}")
            projY_p = h2_p.ProjectionY(f"projY_p_{partontype}_{target_jetpt}_{cut_suffix}")
            projX_p.SetDirectory(0)
            projY_p.SetDirectory(0)

        # --- self-normalize each projection ---
        projections = [p for p in (projX_p, projY_p, projX_h, projY_h) if p is not None]
        # Include this block if self-normalization is desired
        for h in projections:
            if h.Integral() > 0:
                h.Scale(1.0 / h.Integral())

        if self.include_pythia:
            self.FormatHist(projX_p, Color.BLUE, ROOT.kSolid)
            self.FormatHist(projY_p, Color.BLUE, ROOT.kDashed)
        self.FormatHist(projX_h, Color.RED, ROOT.kSolid)
        self.FormatHist(projY_h, Color.RED, ROOT.kDashed)

        canvas = self.make_canvas("can_nAnB_proj", partontype=partontype,
                                  target_jetpt=target_jetpt, z_cut=z_cut,
                                  title="N_A and N_B projections")
        canvas.cd()

        y_max = max(h.GetMaximum() for h in projections)
        first_hist = projX_h if not self.include_pythia else projX_p
        first_hist.SetMaximum(y_max * 1.4)
        first_hist.SetMinimum(0)
        first_hist.GetXaxis().SetTitle("N particles in subjet")
        first_hist.GetYaxis().SetTitle("self-normalized counts") # "counts") #

        first = True
        for h in projections:
            h.Draw("HIST" if first else "HIST SAME")
            first = False

        gen_label = "PYTHIA and HERWIG" if self.include_pythia else "HERWIG"
        ev_leg = self.MakeEventLeg(gen_label, partontype, target_jetpt, z_cut)
        legend = ROOT.TLegend(0.55, 0.62, 0.88, 0.88)
        if self.include_pythia:
            legend.AddEntry(projX_p, "PYTHIA N in A", "l")
            legend.AddEntry(projY_p, "PYTHIA N in B", "l")
        legend.AddEntry(projX_h, "HERWIG N in A", "l")
        legend.AddEntry(projY_h, "HERWIG N in B", "l")
        ev_leg.Draw()
        legend.Draw()

        gen_tag = "PYTHIA_VS_HERWIG" if self.include_pythia else "HERWIG"
        output_name = (f"nA_nB_projections_{gen_tag}_{partontype}"
                       f"_jetpt{target_jetpt}_R0.4_{cut_suffix}.pdf")
        subdir = "pythia_vs_herwig" if self.include_pythia else "herwig"
        self._save_canvas(canvas, subdir, z_cut, output_name,
                          plot_type="nA_nB/projections")

    # -------------------------------------------------------------------------
    # Plot: N_A vs N_B 2D correlation, PYTHIA vs HERWIG comparison
    # (only meaningful when both generators are present)
    # -------------------------------------------------------------------------

    def plot_nA_nB_2d(self, partontype, target_jetpt, z_cut):
        self.set_current_rootfiles(target_jetpt)
        cut_suffix = self.get_cut_suffix(z_cut)

        if not self.include_pythia:
            # no pythia -> just draw the herwig 2D on its own (raw counts, logz)
            h2_h = self.GetNA_NB(self.herwig_rootfile, partontype, target_jetpt, z_cut)
            if h2_h is None:
                return
            hh = h2_h.Clone(f"h2h_only_{partontype}_{target_jetpt}_{cut_suffix}")
            hh.SetDirectory(0)

            canvas = self.make_canvas("can_nAnB_2d_herwig", partontype=partontype,
                                      target_jetpt=target_jetpt, z_cut=z_cut,
                                      title="N_A vs N_B 2D (HERWIG)")
            canvas.cd()
            ROOT.gPad.SetRightMargin(0.15)
            ROOT.gPad.SetLogz()
            hh.SetTitle(f"HERWIG {partontype} (p_{{T}}={target_jetpt});N in A;N in B")
            hh.Draw("COLZ")

            output_name = (f"nA_nB_2d_HERWIG_{partontype}"
                           f"_jetpt{target_jetpt}_R0.4_{cut_suffix}.pdf")
            self._save_canvas(canvas, "herwig", z_cut, output_name,
                              plot_type="nA_nB/2d")
            return

        h2_p = self.GetNA_NB(self.pythia_rootfile, partontype, target_jetpt, z_cut)
        h2_h = self.GetNA_NB(self.herwig_rootfile, partontype, target_jetpt, z_cut)
        if h2_p is None or h2_h is None:
            return

        # raw counts (no normalization)
        hp = h2_p.Clone(f"h2p_{partontype}_{target_jetpt}_{cut_suffix}")
        hh = h2_h.Clone(f"h2h_{partontype}_{target_jetpt}_{cut_suffix}")
        hp.SetDirectory(0)
        hh.SetDirectory(0)

        # ----- (a) side-by-side COLZ -----
        canvas = self.make_canvas("can_nAnB_2d", partontype=partontype,
                                  target_jetpt=target_jetpt, z_cut=z_cut,
                                  title="N_A vs N_B 2D", w=1200, h=550)
        canvas.Divide(2, 1, 0.005, 0.005)

        canvas.cd(1)
        ROOT.gPad.SetRightMargin(0.15)
        ROOT.gPad.SetLogz()
        hp.SetTitle(f"PYTHIA {partontype} (p_{{T}}={target_jetpt});N in A;N in B")
        hp.Draw("COLZ")

        canvas.cd(2)
        ROOT.gPad.SetRightMargin(0.15)
        ROOT.gPad.SetLogz()
        hh.SetTitle(f"HERWIG {partontype} (p_{{T}}={target_jetpt});N in A;N in B")
        hh.Draw("COLZ")

        output_name = (f"nA_nB_2d_sidebyside_PYTHIA_VS_HERWIG_{partontype}"
                       f"_jetpt{target_jetpt}_R0.4_{cut_suffix}.pdf")
        self._save_canvas(canvas, "pythia_vs_herwig", z_cut, output_name,
                          plot_type="nA_nB/2d")

        # ----- (b) overlay: PYTHIA COLZ + HERWIG contour lines -----
        canvas2 = self.make_canvas("can_nAnB_2d_overlay", partontype=partontype,
                                   target_jetpt=target_jetpt, z_cut=z_cut,
                                   title="N_A vs N_B overlay")
        canvas2.cd()
        ROOT.gPad.SetRightMargin(0.15)
        ROOT.gPad.SetLogz()
        hp.SetTitle(f"{partontype} (p_{{T}}={target_jetpt}): "
                    f"PYTHIA (color) vs HERWIG (contours);N in A;N in B")
        hp.Draw("COLZ")
        hh.SetContour(6)
        hh.SetLineColor(ROOT.kBlack)
        hh.Draw("CONT3 SAME")

        leg = ROOT.TLegend(0.15, 0.78, 0.45, 0.88)
        leg.AddEntry(hh, "HERWIG (contours)", "l")
        leg.Draw()

        output_name2 = (f"nA_nB_2d_overlay_PYTHIA_VS_HERWIG_{partontype}"
                        f"_jetpt{target_jetpt}_R0.4_{cut_suffix}.pdf")
        self._save_canvas(canvas2, "pythia_vs_herwig", z_cut, output_name2,
                          plot_type="nA_nB/2d")

        # ----- (c) ratio: PYTHIA / HERWIG (2D, raw counts) -----
        canvas3 = self.make_canvas("can_nAnB_2d_ratio", partontype=partontype,
                                   target_jetpt=target_jetpt, z_cut=z_cut,
                                   title="N_A vs N_B ratio")
        canvas3.cd()
        ROOT.gPad.SetRightMargin(0.15)
        ratio = hp.Clone(f"ratio2d_{partontype}_{target_jetpt}_{cut_suffix}")
        ratio.SetDirectory(0)
        ratio.Divide(hh)
        ratio.SetTitle(f"{partontype} (p_{{T}}={target_jetpt}): "
                       f"PYTHIA / HERWIG;N in A;N in B")
        ratio.SetMinimum(0.0)
        ratio.SetMaximum(2.0)
        ratio.Draw("COLZ")

        output_name3 = (f"nA_nB_2d_ratio_PYTHIA_VS_HERWIG_{partontype}"
                        f"_jetpt{target_jetpt}_R0.4_{cut_suffix}.pdf")
        self._save_canvas(canvas3, "pythia_vs_herwig", z_cut, output_name3,
                          plot_type="nA_nB/2d")

    # -------------------------------------------------------------------------
    # Plot: 1D counting histograms, PYTHIA vs HERWIG  (or single-gen fallback)
    # -------------------------------------------------------------------------

    def plot_counting_hvp(self, partontype, target_jetpt, z_cut):
        self.set_current_rootfiles(target_jetpt)

        # if pythia is off, plot herwig distributions on their own
        if not self.include_pythia:
            hh = self.GetCountingHists(self.herwig_rootfile, partontype,
                                       target_jetpt, z_cut)
            self._plot_1d_single(
                hh["ungroomed"], "herwig", partontype, target_jetpt, z_cut,
                xtitle="N particles (ungroomed jet)", ytitle="counts",
                file_tag="nUngroomed", plot_type="multiplicity")
            self._plot_1d_single(
                hh["groomed"], "herwig", partontype, target_jetpt, z_cut,
                xtitle="N particles (groomed jet)", ytitle="counts",
                file_tag="nGroomed", plot_type="multiplicity")
            # self._plot_1d_single(
            #     hh["combTotal"], "herwig", partontype, target_jetpt, z_cut,
            #     xtitle="N_{total} combinations", ytitle="counts",
            #     file_tag="combTotal", plot_type="combinations", logy=True)
            return

        hp = self.GetCountingHists(self.pythia_rootfile, partontype, target_jetpt, z_cut)
        hh = self.GetCountingHists(self.herwig_rootfile, partontype, target_jetpt, z_cut)

        # multiplicity distributions
        self._plot_1d_hvp(
            hp["ungroomed"], hh["ungroomed"], partontype, target_jetpt, z_cut,
            xtitle="N particles (ungroomed jet)",
            ytitle="counts", legend_title="N ungroomed",
            file_tag="nUngroomed", plot_type="multiplicity")

        self._plot_1d_hvp(
            hp["groomed"], hh["groomed"], partontype, target_jetpt, z_cut,
            xtitle="N particles (groomed jet)",
            ytitle="counts", legend_title="N groomed",
            file_tag="nGroomed", plot_type="multiplicity")

        # # only the total combination distribution (individual AA/AB/BB removed)
        # self._plot_1d_hvp(
        #     hp["combTotal"], hh["combTotal"], partontype, target_jetpt, z_cut,
        #     xtitle="N_{total} combinations",
        #     ytitle="counts", legend_title="Total combinations",
        #     file_tag="combTotal", plot_type="combinations", logy=True)

    # -------------------------------------------------------------------------
    # Plot: all comb components overlaid (single generator, one canvas)
    # -------------------------------------------------------------------------

    def plot_comb_components(self, gen, partontype, target_jetpt, z_cut):
        self.set_current_rootfiles(target_jetpt)
        file = self.get_file(gen)
        cut_suffix = self.get_cut_suffix(z_cut)

        hd = self.GetCountingHists(file, partontype, target_jetpt, z_cut)

        # use combTotal as the reference binning and rebin AA/AB/BB to match
        ref = hd["combTotal"]
        combAA = self._rebin_to_match(hd["combAA"], ref)
        combBB = self._rebin_to_match(hd["combBB"], ref)
        combAB = self._rebin_to_match(hd["combAB"], ref)

        comps = [("combAA", combAA, Color.BLUE),
                 ("combBB", combBB, Color.ORANGE),
                 ("combAB", combAB, Color.GREEN),
                 ("combTotal", ref, ROOT.kBlack)]
        comps = [(name, h, c) for name, h, c in comps if h is not None]
        if not comps:
            return

        canvas = self.make_canvas("can_comb_comp", gen=gen, partontype=partontype,
                                  target_jetpt=target_jetpt, z_cut=z_cut,
                                  title="combination components")
        canvas.SetLogy()
        canvas.cd()

        y_max = max(h.GetMaximum() for _, h, _ in comps)
        first = True
        legend = ROOT.TLegend(0.6, 0.65, 0.88, 0.88)
        for name, h, color in comps:
            self.FormatHist(h, color, ROOT.kSolid)
            if first:
                h.SetMaximum(y_max * 5)
                h.GetXaxis().SetTitle("N combinations")
                h.GetYaxis().SetTitle("counts")
                h.Draw("HIST")
                first = False
            else:
                h.Draw("HIST SAME")
            legend.AddEntry(h, name, "l")

        ev_leg = self.MakeEventLeg(gen.upper(), partontype, target_jetpt, z_cut)
        ev_leg.Draw()
        legend.Draw()

        output_name = (f"comb_components_{gen}_{partontype}"
                       f"_jetpt{target_jetpt}_R0.4_{cut_suffix}.pdf")
        self._save_canvas(canvas, gen, z_cut, output_name,
                          plot_type="combinations")

    # -------------------------------------------------------------------------
    # Plot: all comb components overlaid, PYTHIA vs HERWIG on one canvas
    # (color = component, line style = generator; raw counts, logy)
    # -------------------------------------------------------------------------

    def plot_comb_components_hvp(self, partontype, target_jetpt, z_cut):
        self.set_current_rootfiles(target_jetpt)
        cut_suffix = self.get_cut_suffix(z_cut)

        if not self.include_pythia:
            # nothing to compare against; fall back to single-gen per generator
            for gen in self.generators:
                self.plot_comb_components(gen, partontype, target_jetpt, z_cut)
            return

        hp = self.GetCountingHists(self.pythia_rootfile, partontype, target_jetpt, z_cut)
        hh = self.GetCountingHists(self.herwig_rootfile, partontype, target_jetpt, z_cut)

        # rebin AA/AB/BB to combTotal's binning, per generator
        ref_p = hp["combTotal"]
        ref_h = hh["combTotal"]

        comb_p = {
            "combAA":    self._rebin_to_match(hp["combAA"], ref_p),
            "combBB":    self._rebin_to_match(hp["combBB"], ref_p),
            "combAB":    self._rebin_to_match(hp["combAB"], ref_p),
            "combTotal": ref_p,
        }
        comb_h = {
            "combAA":    self._rebin_to_match(hh["combAA"], ref_h),
            "combBB":    self._rebin_to_match(hh["combBB"], ref_h),
            "combAB":    self._rebin_to_match(hh["combAB"], ref_h),
            "combTotal": ref_h,
        }

        # color per component, line style per generator
        comp_colors = {
            "combAA":    Color.BLUE,
            "combBB":    Color.ORANGE,
            "combAB":    Color.GREEN,
            "combTotal": ROOT.kBlack,
        }
        comp_order = ["combAA", "combBB", "combAB", "combTotal"]

        # collect drawable (name, hist, color, style, gen_label) tuples
        curves = []
        for name in comp_order:
            cp = comb_p.get(name)
            ch = comb_h.get(name)
            if cp is not None:
                curves.append((name, cp, comp_colors[name], ROOT.kSolid,  "PYTHIA"))
            if ch is not None:
                curves.append((name, ch, comp_colors[name], ROOT.kDashed, "HERWIG"))

        curves = [c for c in curves if c[1] is not None]
        if not curves:
            return

        canvas = self.make_canvas("can_comb_comp_hvp", partontype=partontype,
                                  target_jetpt=target_jetpt, z_cut=z_cut,
                                  title="combination components PYTHIA vs HERWIG")
        canvas.SetLogy()
        canvas.cd()

        y_max = max(h.GetMaximum() for _, h, _, _, _ in curves)

        legend = ROOT.TLegend(0.58, 0.55, 0.88, 0.88)
        first = True
        for name, h, color, style, gen_label in curves:
            self.FormatHist(h, color, style)
            if first:
                h.SetMaximum(y_max * 5)
                h.GetXaxis().SetTitle("N combinations")
                h.GetYaxis().SetTitle("counts")
                h.Draw("HIST")
                first = False
            else:
                h.Draw("HIST SAME")
            legend.AddEntry(h, f"{gen_label} {name}", "l")

        ev_leg = self.MakeEventLeg("PYTHIA and HERWIG", partontype,
                                   target_jetpt, z_cut)
        ev_leg.Draw()
        legend.Draw()

        output_name = (f"comb_components_PYTHIA_VS_HERWIG_{partontype}"
                       f"_jetpt{target_jetpt}_R0.4_{cut_suffix}.pdf")
        self._save_canvas(canvas, "pythia_vs_herwig", z_cut, output_name,
                          plot_type="combinations")
    
    
    # -------------------------------------------------------------------------
    # Main plot loop
    # -------------------------------------------------------------------------

    def plot(self):
        for target_jetpt in self.target_jet_pts:
            for cut_mode, z_cut in self.cut_modes:
                self.cut_mode = cut_mode
                self.z_cut    = z_cut
                print(f"Processing {cut_mode} mode, z_cut={z_cut}, "
                      f"jetpt={target_jetpt}...")

                for partontype in self.partontypes:
                    # --- N_A vs N_B ---
                    self.plot_nA_nB_projections(partontype, target_jetpt, z_cut)
                    self.plot_nA_nB_2d(partontype, target_jetpt, z_cut)

                    # --- 1D counting: PYTHIA vs HERWIG (or single-gen) ---
                    self.plot_counting_hvp(partontype, target_jetpt, z_cut)

                    # --- comb components overlay, per enabled generator ---
                    for gen in self.generators:
                        self.plot_comb_components(gen, partontype, target_jetpt, z_cut)
                    
                    # --- comb components overlay, PYTHIA vs HERWIG ---
                    self.plot_comb_components_hvp(partontype, target_jetpt, z_cut)


# =============================================================================
# Entry point
# =============================================================================

if __name__ == "__main__":
    plotter = PlotCounting()
    plotter.plot()