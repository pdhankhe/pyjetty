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

PYTHIA_BASE = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/pythia_otf/55555648"
HERWIG_BASE = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/herwiggen/tree_gen/55293842"

OUTPUT_DIR  = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots/"
OUTPUT_NAME = "qg_fraction_vs_pt.pdf"

LABEL_FRAC = 0.8       # pthat label = LABEL_FRAC * ptbin
N_PT_BINS = 20

def make_path(gen, ptbin):
    base = PYTHIA_BASE if gen == "pythia" else HERWIG_BASE
    return f"{base}/{ptbin}gev/JetsForAnalysisCombined.parquet"

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
    # guard against non-positive lower edge (log requires > 0)
    if lo <= 0:
        lo = min(r[1] for r in ranges) * 1e-3  # small positive fallback
    pt_edges = np.logspace(np.log10(lo), np.log10(hi), N_PT_BINS + 1)
    centers = np.sqrt(pt_edges[:-1] * pt_edges[1:])  # geometric centers

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

# ---------------------------------------------------------------------------
# Plot
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
    ax.set_xscale("log")          # <-- add this line
    ax.set_ylabel("flavor fraction")
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3, which="both")   # show minor gridlines too
    ax.legend(ncol=2, fontsize=8)

for ax in axes[n_bins:]:
    ax.set_visible(False)

fig.suptitle("Quark / Gluon Jet Fraction vs. jet " r"$p_T$"
             "  (anti-kT R = 0.4) (q-frac = $N_q/(N_q+N_g)$)", y=1.0)
fig.tight_layout()

os.makedirs(OUTPUT_DIR, exist_ok=True)
out_path = os.path.join(OUTPUT_DIR, OUTPUT_NAME)
fig.savefig(out_path, format="pdf", bbox_inches="tight")
plt.show()
print(f"Saved plot to {out_path}")