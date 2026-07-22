import os
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
GENERATORS = ["pythia", "herwig"]
target_jet_pts = [50, 100, 200, 500]
PTHAT_VALUES = [40, 80, 160, 400]

PYTHIA_BASE = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/pythia_otf/55555648"
HERWIG_BASE = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/herwiggen/tree_gen/55293842"

OUTPUT_DIR  = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots/"
OUTPUT_NAME = "qg_fraction_vs_pt.pdf"

LABEL_FRAC = 0.8       # pthat label = LABEL_FRAC * ptbin
N_BINS_PER_SEGMENT = 5 # Approximate total bins = len(PTHAT_VALUES) * N_BINS_PER_SEGMENT

def make_path(gen, ptbin):
    base = PYTHIA_BASE if gen == "pythia" else HERWIG_BASE
    return f"{base}/{ptbin}gev/JetsForAnalysisCombined.parquet"

def get_pt_edges(lo, hi, pthat_values):
    """Create log-spaced edges that explicitly include pthat boundaries."""
    anchors = [lo]
    for p in pthat_values:
        if lo < p < hi:
            anchors.append(p)
    anchors.append(hi)

    edges = []
    for i in range(len(anchors) - 1):
        # Create bins between each anchor point
        seg = np.logspace(np.log10(anchors[i]), np.log10(anchors[i+1]), N_BINS_PER_SEGMENT + 1)
        edges.extend(seg[:-1])
    edges.append(hi)
    return np.array(edges)

# ---------------------------------------------------------------------------
# Stream the file in batches; dedup within each batch only.
# ---------------------------------------------------------------------------
def accumulate_counts(gen, ptbin, pt_edges):
    path = make_path(gen, ptbin)
    if not os.path.exists(path):
        print(f"  [MISSING] {path}")
        return None

    cols = pq.ParquetFile(path).schema.names
    if "parton_pid" not in cols:
        print(f"  [SKIP] {gen} {ptbin}gev has no parton_pid column")
        return None

    n_pt = len(pt_edges) - 1
    nq = np.zeros(n_pt, dtype=np.int64)
    ng = np.zeros(n_pt, dtype=np.int64)

    pf = pq.ParquetFile(path)
    for batch in pf.iter_batches(
        batch_size=2_000_000,
        columns=["event_id", "jet_id", "jet_pt", "parton_pid"],
    ):
        df = batch.to_pandas().drop_duplicates(subset=["event_id", "jet_id"])

        apid = df["parton_pid"].abs().values
        is_q = (apid >= 1) & (apid <= 6)
        is_g = (apid == 21)

        bin_idx = np.digitize(df["jet_pt"].values, pt_edges) - 1
        valid = (bin_idx >= 0) & (bin_idx < n_pt)

        nq += np.bincount(bin_idx[valid & is_q], minlength=n_pt)[:n_pt]
        ng += np.bincount(bin_idx[valid & is_g], minlength=n_pt)[:n_pt]

        del df, batch

    return nq, ng

# ---------------------------------------------------------------------------
# Cheap jet_pt range from parquet column statistics
# ---------------------------------------------------------------------------
def jet_pt_range(gen, ptbin):
    path = make_path(gen, ptbin)
    if not os.path.exists(path):
        return None
    pf = pq.ParquetFile(path)
    lo, hi = np.inf, -np.inf
    col_idx = pf.schema.names.index("jet_pt")
    for rg in range(pf.num_row_groups):
        stats = pf.metadata.row_group(rg).column(col_idx).statistics
        if stats is not None and stats.has_min_max:
            lo = min(lo, stats.min)
            hi = max(hi, stats.max)
    if not np.isfinite(lo) or not np.isfinite(hi):
        return None
    return lo, hi

def find_crossing(x, y):
    """Find x where y crosses 0.5 using linear interpolation."""
    for i in range(len(y) - 1):
        if (y[i] < 0.5 <= y[i+1]) or (y[i] > 0.5 >= y[i+1]):
            return x[i] + (0.5 - y[i]) * (x[i+1] - x[i]) / (y[i+1] - y[i])
    return None

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
results = {}

for ptbin in target_jet_pts:
    ranges = [r for r in (jet_pt_range(g, ptbin) for g in GENERATORS)
              if r is not None]
    if not ranges:
        continue
    lo = min(r[0] for r in ranges)
    hi = max(r[1] for r in ranges)
    if lo <= 0:
        lo = min(r[1] for r in ranges) * 1e-3

    pt_edges = get_pt_edges(lo, hi, PTHAT_VALUES)
    centers = np.sqrt(pt_edges[:-1] * pt_edges[1:])

    for gen in GENERATORS:
        counts = accumulate_counts(gen, ptbin, pt_edges)
        if counts is None:
            continue
        nq, ng = counts
        ntot = nq + ng
        with np.errstate(invalid="ignore", divide="ignore"):
            qfrac = np.where(ntot > 0, nq / ntot, np.nan)
            gfrac = np.where(ntot > 0, ng / ntot, np.nan)
            qerr = np.where(ntot > 0,
                            np.sqrt(qfrac * (1 - qfrac) / ntot), np.nan)
            gerr = np.where(ntot > 0,
                            np.sqrt(gfrac * (1 - gfrac) / ntot), np.nan)
        results[(gen, ptbin)] = (centers, qfrac, gfrac, qerr, gerr)
        print(f"{gen:7s} {ptbin}gev: {int(ntot.sum())} q/g jets")

# Find crossing points
print("\n--- 50% Crossing Points (jet pT) ---")
for ptbin in target_jet_pts:
    pthat = LABEL_FRAC * ptbin
    for gen in GENERATORS:
        if (gen, ptbin) in results:
            centers, _, gfrac, _, _ = results[(gen, ptbin)]
            cross = find_crossing(centers, gfrac)
            print(f"{gen:7s} pthat={pthat:g} GeV: {cross:.2f} GeV" if cross else f"{gen:7s} pthat={pthat:g} GeV: Not found")

# ---------------------------------------------------------------------------
# Plot 1: qg_fraction_vs_pt.pdf
# ---------------------------------------------------------------------------
n_bins = len(target_jet_pts)
ncols = 2
nrows = int(np.ceil(n_bins / ncols))

fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4.5 * nrows),
                         squeeze=False)
axes = axes.flatten()

style = {
    ("pythia", "quark"): dict(color="tab:blue", marker="o", ls="-"),
    ("pythia", "gluon"): dict(color="tab:blue", marker="s", ls="--", markerfacecolor="none"),
    ("herwig", "quark"): dict(color="tab:red",  marker="o", ls="-"),
    ("herwig", "gluon"): dict(color="tab:red",  marker="s", ls="--", markerfacecolor="none"),
}

for ax, ptbin in zip(axes, target_jet_pts):
    pthat_label = LABEL_FRAC * ptbin
    present = [g for g in GENERATORS if (g, ptbin) in results]
    if not present:
        ax.set_title(rf"$\hat{{p}}_T$ = {pthat_label:g} GeV (no data)")
        continue

    ax.axvspan(0, pthat_label, color='grey', alpha=0.3, zorder=0)

    for gen in present:
        centers, qfrac, gfrac, qerr, gerr = results[(gen, ptbin)]
        ax.errorbar(centers, qfrac, yerr=qerr,
                    label=f"{gen} quark", **style[(gen, "quark")],
                    markersize=4, capsize=2)
        ax.errorbar(centers, gfrac, yerr=gerr,
                    label=f"{gen} gluon", **style[(gen, "gluon")],
                    markersize=4, capsize=2)

    ax.set_title(rf"$\hat{{p}}_T$ = {pthat_label:g} GeV")
    ax.set_xlabel(r"jet $p_T$ [GeV]")
    ax.set_xscale("log")
    ax.set_ylabel("flavor fraction")
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(ncol=2, fontsize=8)

for ax in axes[n_bins:]:
    ax.set_visible(False)

fig.suptitle("Quark / Gluon Jet Fraction vs. jet " r"$p_T$"
             "  (anti-kT R = 0.4) (q-frac = $N_q/(N_q+N_g)$)", y=1.0)
fig.tight_layout()

os.makedirs(OUTPUT_DIR, exist_ok=True)
out_path = os.path.join(OUTPUT_DIR, OUTPUT_NAME)
fig.savefig(out_path, format="pdf", bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Plot 2: Gluon curves combined
# ---------------------------------------------------------------------------
fig_gluon, ax_gluon = plt.subplots(figsize=(8, 6))
colors = {"pythia": "tab:blue", "herwig": "tab:red"}
markers = {50: "o", 100: "s", 200: "^", 500: "d"}
linestyles = {50: "-", 100: "--", 200: "-.", 500: ":"}

for gen in GENERATORS:
    # Exclude pthat=400 (ptbin 500)
    for ptbin in [50, 100, 200]:
        if (gen, ptbin) in results:
            centers, _, gfrac, _, gerr = results[(gen, ptbin)]
            pthat = LABEL_FRAC * ptbin

            # Remove points where jet pt <= pthat
            mask = centers > pthat
            c_filt = centers[mask]
            g_filt = gfrac[mask]
            e_filt = gerr[mask]

            ax_gluon.errorbar(c_filt, g_filt, yerr=e_filt,
                              color=colors[gen], marker=markers[ptbin],
                              ls=linestyles[ptbin],
                              label=f"{gen} $\hat{{p}}_T$={pthat:g}",
                              markersize=4, capsize=2, alpha=0.7)

ax_gluon.axhline(0.5, color='black', linestyle=':', alpha=0.5, label='50% fraction')
ax_gluon.set_title("Combined Gluon Fractions")
ax_gluon.set_xlabel(r"jet $p_T$ [GeV]")
ax_gluon.set_ylabel("Gluon fraction")
ax_gluon.set_xscale("log")
ax_gluon.set_ylim(0, 1)
ax_gluon.grid(True, alpha=0.3, which="both")
ax_gluon.legend(ncol=2, fontsize=8)
fig_gluon.savefig(os.path.join(OUTPUT_DIR, "gluon_fractions_combined.pdf"), bbox_inches="tight")
plt.close(fig_gluon)

# ---------------------------------------------------------------------------
# Plot 3: Ratio of pthat=80 to pthat=40
# ---------------------------------------------------------------------------
fig_ratio, ax_ratio = plt.subplots(figsize=(8, 6))
for gen in GENERATORS:
    if (gen, 50) in results and (gen, 100) in results:
        # pthat=40 (ptbin 50), pthat=80 (ptbin 100)
        x40, _, g40, _, _ = results[(gen, 50)]
        x80, _, g80, _, _ = results[(gen, 100)]

        g40_interp = np.interp(x80, x40, g40)
        ratio = g80 / g40_interp

        ax_ratio.plot(x80, ratio, label=gen, color=colors[gen], marker='o', markersize=4)

ax_ratio.set_title(r"Ratio of Gluon Fraction ($\hat{p}_T=80$ / $\hat{p}_T=40$)")
ax_ratio.set_xlabel(r"jet $p_T$ [GeV]")
ax_ratio.set_ylabel("Ratio")
ax_ratio.set_xscale("log")
ax_ratio.grid(True, alpha=0.3, which="both")
ax_ratio.legend()
fig_ratio.savefig(os.path.join(OUTPUT_DIR, "gluon_ratio_80_40.pdf"), bbox_inches="tight")
plt.close(fig_ratio)

# ---------------------------------------------------------------------------
# Plot 4: Combined data plot (Stitched)
# ---------------------------------------------------------------------------
stitch_config = [
    (50, 40, 80),
    (100, 80, 160),
    (200, 160, 400),
    (500, 400, 500),
]

fig_stitch, ax_stitch = plt.subplots(figsize=(8, 6))
for gen in GENERATORS:
    all_x, all_q, all_g = [], [], []
    for ptbin, pt_min, pt_max in stitch_config:
        if (gen, ptbin) in results:
            centers, qfrac, gfrac, _, _ = results[(gen, ptbin)]
            mask = (centers >= pt_min) & (centers <= pt_max)
            all_x.extend(centers[mask])
            all_q.extend(qfrac[mask])
            all_g.extend(gfrac[mask])

    all_x = np.array(all_x)
    all_q = np.array(all_q)
    all_g = np.array(all_g)

    idx = np.argsort(all_x)
    ax_stitch.plot(all_x[idx], all_q[idx], label=f"{gen} quark", color=style[(gen, "quark")]["color"], marker='o', markersize=3)
    ax_stitch.plot(all_x[idx], all_g[idx], label=f"{gen} gluon", color=style[(gen, "gluon")]["color"], marker='s', markersize=3)

ax_stitch.set_title("Stitched Quark/Gluon Fractions")
ax_stitch.set_xlabel(r"jet $p_T$ [GeV]")
ax_stitch.set_ylabel("fraction")
ax_stitch.set_xscale("log")
ax_stitch.set_ylim(0, 1)
ax_stitch.grid(True, alpha=0.3, which="both")
ax_stitch.legend()
fig_stitch.savefig(os.path.join(OUTPUT_DIR, "composite_qg_fraction.pdf"), bbox_inches="tight")
plt.close(fig_stitch)

print(f"Saved all plots to {OUTPUT_DIR}")
