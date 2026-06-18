import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# Which generators to include. Add "pythia" back when you have the file.
GENERATORS = ["herwig"]

# Fill these in with your actual pt-hat bin values (the {..}gev directory names)
target_jet_pts = [50, 100, 200, 500]   # <-- EDIT to match your bins

PYTHIA_BASE = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/pythia_otf/53423546"
HERWIG_BASE = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/herwiggen/tree_gen/54380351"
# HERWIG_BASE = "/rstorage/generators/herwig_alice/tree_gen/1006458"

# Output: full path is OUTPUT_DIR / OUTPUT_NAME
OUTPUT_DIR  = "/software/users/blianggi/mypyjetty/storage/jse/plots/"
OUTPUT_DIR  = "/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots/"
OUTPUT_NAME = "qg_fraction_vs_pt.pdf"

def make_path(gen, ptbin):
    base = PYTHIA_BASE if gen == "pythia" else HERWIG_BASE
    return f"{base}/{ptbin}gev/JetsForAnalysisCombined.parquet"

# jet_pt histogram binning used to draw the fraction curves
N_PT_BINS = 20      # number of bins across the jet_pt range in each panel

# ---------------------------------------------------------------------------
# Flavor classification
# ---------------------------------------------------------------------------
def classify_flavor(pid):
    apid = abs(int(pid))
    if 1 <= apid <= 6:
        return "quark"
    elif apid == 21:
        return "gluon"
    else:
        return "unidentified"

# ---------------------------------------------------------------------------
# Load file -> per-jet records (quark/gluon only)
# ---------------------------------------------------------------------------
def load_jet_level(gen, ptbin):
    path = make_path(gen, ptbin)
    if not os.path.exists(path):
        print(f"  [warning] missing file: {path}")
        return None

    df = pd.read_parquet(
        path, columns=["event_id", "jet_id", "jet_pt", "parton_pid"]
    )
    # one row per jet (jet_pt / parton_pid constant within a jet)
    jets = df.drop_duplicates(subset=["event_id", "jet_id"]).copy()
    jets["flavor"] = jets["parton_pid"].apply(classify_flavor)

    # drop unidentified jets entirely, per your definition
    jets = jets[jets["flavor"].isin(["quark", "gluon"])].copy()
    return jets

# ---------------------------------------------------------------------------
# Compute quark/gluon fraction vs jet_pt for one jet sample
# ---------------------------------------------------------------------------
def fraction_vs_pt(jets, pt_edges):
    centers = 0.5 * (pt_edges[:-1] + pt_edges[1:])
    bin_idx = np.digitize(jets["jet_pt"].values, pt_edges) - 1

    is_quark = (jets["flavor"] == "quark").values
    is_gluon = (jets["flavor"] == "gluon").values

    qfrac = np.full(len(centers), np.nan)
    gfrac = np.full(len(centers), np.nan)
    qerr  = np.full(len(centers), np.nan)
    gerr  = np.full(len(centers), np.nan)

    for b in range(len(centers)):
        mask = bin_idx == b
        n = mask.sum()
        if n == 0:
            continue
        nq = is_quark[mask].sum()
        ng = is_gluon[mask].sum()       # nq + ng == n (only q/g remain)
        qfrac[b] = nq / n
        gfrac[b] = ng / n
        # binomial error on the fraction
        qerr[b] = np.sqrt(qfrac[b] * (1 - qfrac[b]) / n)
        gerr[b] = np.sqrt(gfrac[b] * (1 - gfrac[b]) / n)

    return centers, qfrac, gfrac, qerr, gerr

# ---------------------------------------------------------------------------
# Load everything
# ---------------------------------------------------------------------------
data = {}   # data[(gen, ptbin)] = jets dataframe
for ptbin in target_jet_pts:
    for gen in GENERATORS:
        jets = load_jet_level(gen, ptbin)
        if jets is not None and len(jets):
            data[(gen, ptbin)] = jets
            print(f"{gen:7s} {ptbin}gev: {len(jets)} q/g jets")

# ---------------------------------------------------------------------------
# Plot: one panel per pt-hat bin
# ---------------------------------------------------------------------------
n_bins = len(target_jet_pts)
ncols = 2
nrows = int(np.ceil(n_bins / ncols))

fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4.5 * nrows),
                         squeeze=False)
axes = axes.flatten()

style = {
    ("pythia", "quark"): dict(color="tab:blue", marker="o", ls="-"),
    ("pythia", "gluon"): dict(color="tab:blue", marker="s", ls="--"),
    ("herwig", "quark"): dict(color="tab:red",  marker="o", ls="-"),
    ("herwig", "gluon"): dict(color="tab:red",  marker="s", ls="--"),
}

for ax, ptbin in zip(axes, target_jet_pts):
    present = [data[(g, ptbin)] for g in GENERATORS if (g, ptbin) in data]
    if not present:
        ax.set_title(f"{ptbin} GeV (no data)")
        continue

    all_pt = np.concatenate([j["jet_pt"].values for j in present])
    lo, hi = np.percentile(all_pt, [1, 99])   # trim extreme tails
    pt_edges = np.linspace(lo, hi, N_PT_BINS + 1)

    for gen in GENERATORS:
        if (gen, ptbin) not in data:
            continue
        centers, qfrac, gfrac, qerr, gerr = fraction_vs_pt(
            data[(gen, ptbin)], pt_edges
        )
        ax.errorbar(centers, qfrac, yerr=qerr,
                    label=f"{gen} quark", **style[(gen, "quark")],
                    markersize=4, capsize=2)
        ax.errorbar(centers, gfrac, yerr=gerr,
                    label=f"{gen} gluon", **style[(gen, "gluon")],
                    markersize=4, capsize=2)

    ax.set_title(rf"$\hat{{p}}_T$ bin: {ptbin} GeV")
    ax.set_xlabel(r"jet $p_T$ [GeV]")
    ax.set_ylabel("flavor fraction")
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)
    ax.legend(ncol=2, fontsize=8)

# hide any unused panels
for ax in axes[n_bins:]:
    ax.set_visible(False)

fig.suptitle("Quark / Gluon Jet Fraction vs. jet " r"$p_T$"
             "  (q-frac = $N_q/(N_q+N_g)$)", y=1.0)
fig.tight_layout()

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
os.makedirs(OUTPUT_DIR, exist_ok=True)
out_path = os.path.join(OUTPUT_DIR, OUTPUT_NAME)
fig.savefig(out_path, format="pdf", bbox_inches="tight")
plt.show()
print(f"Saved plot to {out_path}")