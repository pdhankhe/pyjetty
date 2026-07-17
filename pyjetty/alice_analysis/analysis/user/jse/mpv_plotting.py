# mpv_plotting.py
"""Reusable MPV (Most Probable Value) fitting and plotting utilities."""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def gaus_log(x, mu, C, sg):
    """Log-normal distribution function."""
    return C * np.exp(-(np.log(x / mu)) ** 2 / (2 * sg * sg))


def select_fit_window(xs, ys, yerrs, n_side=3, n_min=2):
    """Restrict (xs, ys, yerrs) to a window around the maximum of ys."""
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


def fit_mpv(hist, fit_range=None, save_path=None, plot_title=None):
    """Extract the peak position (mu) from a ROOT histogram by fitting gaus_log.

    Returns (mu, mu_err) or (None, None) if the fit fails.
    """
    if hist is None:
        return None, None

    xaxis = hist.GetXaxis()
    xmin, xmax = xaxis.GetXmin(), xaxis.GetXmax()
    if xmin <= 0:
        lo = None
        for ib in range(1, hist.GetNbinsX() + 1):
            edge = hist.GetBinLowEdge(ib)
            if edge > 0:
                lo = edge
                break
        xmin = lo if lo is not None else hist.GetBinCenter(1)
    xrange = (xmin, xmax)

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
    xs, ys, yerrs, peak_idx = select_fit_window(xs, ys, yerrs, n_side=3, n_min=2)

    if len(xs) < 3:
        print(f"  Not enough points around peak to fit ({len(xs)} found)")
        if save_path is not None:
            try:
                save_fit_diagnostic(xs_full, ys_full, yerrs_full,
                                    None, save_path, plot_title,
                                    None, None, p0=None, xrange=xrange)
            except Exception:
                pass
        return None, None

    mu0 = xs[np.argmax(ys)]
    C0 = ys.max()
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
        popt, pcov = curve_fit(
            gaus_log, xs, ys,
            p0=p0, sigma=yerrs, absolute_sigma=False, maxfev=5000,
        )
        mu_fit, C_fit, sg_fit = popt
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
            save_fit_diagnostic(
                xs_full, ys_full, yerrs_full,
                popt, save_path, plot_title, mu, mu_err,
                p0=p0, fit_xs=xs, xrange=xrange,
            )
        except Exception as ex:
            print(f"  Failed to save fit diagnostic: {ex}")

    return mu, mu_err


def save_fit_diagnostic(xs, ys, yerrs, popt, save_path,
                        title, mu, mu_err, p0=None, fit_xs=None, xrange=None):
    """Save a matplotlib plot showing data points, initial-guess curve, and fitted curve."""
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 5))

    ax.errorbar(xs, ys, yerr=yerrs, fmt="o", markersize=4,
                color="black", capsize=2, label="data", zorder=2)

    if fit_xs is not None and len(fit_xs):
        mask = np.isin(xs, fit_xs)
        ax.errorbar(xs[mask], ys[mask], yerr=yerrs[mask], fmt="o",
                    markersize=6, mfc="none", mec="C2", mew=1.5,
                    capsize=2, label="fit window", zorder=2.2)
        grid_lo, grid_hi = fit_xs.min(), fit_xs.max()
    else:
        grid_lo, grid_hi = xs.min(), xs.max()

    x_grid = np.logspace(np.log10(grid_lo), np.log10(grid_hi), 400)

    if p0 is not None:
        y_init = gaus_log(x_grid, *p0)
        ax.plot(x_grid, y_init, "--", color="C0", linewidth=1.5, alpha=0.8,
                label=(f"initial guess\n"
                       f"$\\mu_0$ = {p0[0]:.4g}, $C_0$ = {p0[1]:.4g}, "
                       f"$\\sigma_0$ = {p0[2]:.4g}"), zorder=2.5)

    if popt is not None:
        y_fit = gaus_log(x_grid, *popt)
        ax.plot(x_grid, y_fit, "-", color="C3", linewidth=2,
                label=(f"gaus_log fit\n"
                       f"$\\mu$ = {popt[0]:.4g}, $C$ = {popt[1]:.4g}, "
                       f"$\\sigma$ = {popt[2]:.4g}"), zorder=3)
        ax.axvline(popt[0], color="C3", linestyle="--", alpha=0.6,
                   label=f"peak $\\mu$ = {popt[0]:.4g} $\\pm$ {mu_err:.2g}")
    else:
        ax.text(0.5, 0.5, "FIT FAILED", transform=ax.transAxes,
                ha="center", va="center", fontsize=20, color="red", alpha=0.5)

    if p0 is not None:
        ax.axvline(p0[0], color="C0", linestyle=":", alpha=0.5,
                   label=f"seed $\\mu_0$ = {p0[0]:.4g}")

    ax.set_xscale("log")
    ax.set_xlabel("x")
    ax.set_ylabel("entries")
    if title:
        ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, which="both", alpha=0.3)
    if xrange is not None:
        ax.set_xlim(*xrange)
    fig.tight_layout()
    fig.savefig(save_path)
    plt.close(fig)


def plot_mpv_summary(data, title, xlabel, ylabel, outdir, filename,
                     component_styles=None, ylim=None):
    """Generic MPV summary plot: y = peak position vs x = jet pT.

    `data` is a dict {component: {jetpt: (mu, mu_err)}}.
    """
    if not data:
        return

    fig, ax = plt.subplots(figsize=(7, 5))

    default_colors = ["k", "C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"]
    default_markers = ["o", "s", "D", "^", "v", "P", "X", "*", "h", "<"]

    for i, (comp, jetpt_dict) in enumerate(sorted(data.items())):
        if not jetpt_dict:
            continue
        pts = sorted(jetpt_dict.keys())
        mus = [jetpt_dict[p][0] for p in pts]
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
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="best", framealpha=0.9)
    if ylim is not None:
        ax.set_ylim(*ylim)
    fig.tight_layout()

    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, filename))
    plt.close(fig)


def plot_mpv_ratio(num_data, den_data, comp_labels, styles,
                   title, xlabel, ylabel, outdir, filename, ylim=None):
    """Plot ratio of MPV (num/den) vs jet pT for each component."""
    fig, ax = plt.subplots(figsize=(7, 5))
    plotted = False

    for comp in comp_labels:
        nd = num_data.get(comp, {})
        dd = den_data.get(comp, {})
        if not nd or not dd:
            continue
        pts = sorted(set(nd.keys()) & set(dd.keys()))
        if not pts:
            continue

        ratios, rerrs = [], []
        for p in pts:
            n, ne = nd[p]
            d, de = dd[p]
            if n is None or d is None or d == 0:
                continue
            r = n / d
            rel = 0.0
            if ne and n:
                rel += (ne / n) ** 2
            if de and d:
                rel += (de / d) ** 2
            ratios.append(r)
            rerrs.append(abs(r) * np.sqrt(rel))
        if not ratios:
            continue

        style = styles.get(comp, dict(marker="o", linestyle="-", label=comp))
        ax.errorbar(pts[:len(ratios)], ratios, yerr=rerrs, capsize=3,
                    markersize=7, **style)
        plotted = True

    if not plotted:
        plt.close(fig)
        return

    ax.axhline(1.0, color="gray", linestyle=":", alpha=0.7)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11)
    ax.set_xscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="best", framealpha=0.9)
    if ylim is not None:
        ax.set_ylim(*ylim)
    fig.tight_layout()
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, filename))
    plt.close(fig)