import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import json
import uproot

# ─────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────

target_jet_pts = [50, 100, 200, 500]

pythia_jobid = "small_test" #"53423546"
herwig_jobid = "small_test" #"54380351"

# Paths
pythia_base   = f"/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/pythia_otf/{pythia_jobid}"
slurm_base    = f"{pythia_base}/output"
herwig_base   = f"/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/herwiggen/tree_gen/{herwig_jobid}"
herwig_hepmc_base = f"/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/herwiggen/hepmc/{herwig_jobid}"
output_dir    = f"/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/plots/pythia_vs_herwig/small_test"

# File index ranges for each pt bin in the Pythia slurm outputs
slurm_ranges = {
    50:  (1,   1),
    100: (2, 2),
    200: (3, 3),
    500: (4, 4),
}

BIN_WIDTH = 1.0  # GeV

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# ─────────────────────────────────────────────
# DATA LOADING UTILS
# ─────────────────────────────────────────────

def get_root_hist_data(file_path, hist_name):
    """Extracts values and bin edges from a ROOT histogram."""
    try:
        with uproot.open(file_path) as file:
            hist = file[hist_name]
            return hist.values(), hist.axis().edges()
    except Exception as e:
        print(f"Error reading {hist_name} from {file_path}: {e}")
        return None, None

def get_jet_pt_distribution(parquet_path):
    """
    Loads jet data from parquet.
    Assuming each row is a constituent, and 'jet_pt' is the pT of the jet it belongs to.
    """
    try:
        df = pd.read_parquet(parquet_path, columns=['jet_id', 'jet_pt'])
        jet_pts = df.groupby('jet_id')['jet_pt'].first().to_numpy()
        return jet_pts
    except Exception as e:
        print(f"Error reading parquet {parquet_path}: {e}")
        return None

def parse_pythia_slurm_file(filepath):
    """Extract sigmaGen (mb) and N_accepted from a single Pythia slurm output file."""
    sigma = None
    n_accepted = None
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        sum_match = re.search(
            r'\|\s+sum\s+\|\s+(\d+)\s+(\d+)\s+(\d+)\s+\|\s+([\d.e+\-]+)\s+([\d.e+\-]+)',
            content
        )
        # print("pythia slurm:", sum_match)
        if sum_match:
            n_accepted = int(sum_match.group(3))
            sigma = float(sum_match.group(4))  # mb
            print("n_accepted", n_accepted, "sigma", sigma)
    except Exception as e:
        print(f"Error parsing Pythia slurm file {filepath}: {e}")
    return sigma, n_accepted

def get_pythia_scale_factor(target_jet_pt):
    """Parse slurm files for a given target_jet_pt bin. Returns scale_f in mb/event."""
    idx_start, idx_end = slurm_ranges[target_jet_pt]
    all_files = sorted(glob.glob(os.path.join(slurm_base, f"slurm-{pythia_jobid}_*.out")))

    selected_files = [
        f for f in all_files
        if idx_start <= int(re.search(rf'slurm-{pythia_jobid}_(\d+)\.out', f).group(1)) <= idx_end
    ]

    if not selected_files:
        print(f"Warning: No Pythia slurm files found for target_jet_pt={target_jet_pt}")
        return None

    sigma_values, n_accepted_values = [], []
    for filepath in selected_files:
        sigma, n_accepted = parse_pythia_slurm_file(filepath)
        if sigma is not None:
            sigma_values.append(sigma)
            n_accepted_values.append(n_accepted)

    if not sigma_values:
        return None

    sigma_mean = np.mean(sigma_values)
    n_total = int(np.sum(n_accepted_values))
    return sigma_mean / n_total  # mb/event

def get_herwig_scale_factor(target_jet_pt):
    """
    Parse Herwig .out file for a given target_jet_pt bin.
    Expects a file like: .../50gev/LHC_5360_MPI_jse-S1.out
    Returns scale factor in nb/event.
    """
    file_path = os.path.join(herwig_hepmc_base, f"{target_jet_pt}gev", "LHC_5360_MPI_jse-S1.out")
    try:
        with open(file_path, 'r') as f:
            content = f.read()
        # Look for the "Total (from generated events)" line
        # Format: Total (from generated events):            10000        10000      9.44(5)e+03
        match = re.search(r'Total \(from generated events\):\s+(\d+)\s+\d+\s+([\d.e+\- \(\)]+)', content)
        if match:
            n_events = int(match.group(1))
            cs_str = match.group(2).strip()
            # Remove uncertainty (5) -> 9.44e+03
            cs_clean = re.sub(r'\(\d+\)', '', cs_str)
            return float(cs_clean) / n_events # nb/event
    except Exception as e:
        print(f"Error parsing Herwig scale factor from {file_path}: {e}")
    return None

# ─────────────────────────────────────────────
# MAIN ANALYSIS & PLOTTING
# ─────────────────────────────────────────────

def run_comparison():
    # 1. Collect Data
    # ---------------------------------------------------------
    all_data = {}
    for target_pt in target_jet_pts:
        print(f"\n>>> Collecting data for target pT = {target_pt} GeV")
        py_dir = os.path.join(pythia_base, f"{target_pt}gev")
        hw_dir = os.path.join(herwig_base,  f"{target_pt}gev")

        # Data holders
        data = {
            'py_particle_pt': None, 'py_particle_pt_edges': None,
            'hw_particle_pt': None, 'hw_particle_pt_edges': None,
            'py_particle_eta': None, 'py_particle_eta_edges': None,
            'hw_particle_eta': None, 'hw_particle_eta_edges': None,
            'py_jets': None, 'hw_jets': None,
            'py_scale': None, 'hw_scale': None,
        }

        # Root data
        py_vals, py_edges = get_root_hist_data(os.path.join(py_dir, "AnalysisResults.root"), "hparticlepT")
        hw_vals, hw_edges = get_root_hist_data(os.path.join(hw_dir, "AnalysisResults.root"), "hparticlepT")
        data['py_particle_pt'], data['py_particle_pt_edges'] = py_vals, py_edges
        data['hw_particle_pt'], data['hw_particle_pt_edges'] = hw_vals, hw_edges

        py_vals, py_edges = get_root_hist_data(os.path.join(py_dir, "AnalysisResults.root"), "hparticleEta")
        hw_vals, hw_edges = get_root_hist_data(os.path.join(hw_dir, "AnalysisResults.root"), "hparticleEta")
        data['py_particle_eta'], data['py_particle_eta_edges'] = py_vals, py_edges
        data['hw_particle_eta'], data['hw_particle_eta_edges'] = hw_vals, hw_edges

        # Parquet data
        data['py_jets'] = get_jet_pt_distribution(os.path.join(py_dir, "JetsForAnalysis.parquet"))
        data['hw_jets'] = get_jet_pt_distribution(os.path.join(hw_dir, "JetsForAnalysis.parquet"))

        # Scale factors
        data['py_scale'] = get_pythia_scale_factor(target_pt)
        data['hw_scale'] = get_herwig_scale_factor(target_pt)

        all_data[target_pt] = data

    # 2. Plotting Functions
    # ---------------------------------------------------------
    def plot_combined(metric_key, title, xlabel, ylabel, is_xsec=False, is_ratio=False):
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()

        for i, target_pt in enumerate(target_jet_pts):
            ax = axes[i]
            d = all_data[target_pt]

            if is_xsec:
                # Calculate dsigma/dpT
                bin_edges = np.linspace(0, target_pt*2, 51)
                bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
                current_bin_width = (target_pt*2) / 50
                py_counts, _ = np.histogram(d['py_jets'], bins=bin_edges)
                hw_counts, _ = np.histogram(d['hw_jets'], bins=bin_edges)
                py_val = d['py_scale'] * 1e9 * py_counts / current_bin_width
                hw_val = d['hw_scale'] * 1e3 * hw_counts / current_bin_width
                edges = bin_edges[:-1]


                ax.step(edges, py_val, color='blue', label='Pythia', where='post')
                ax.step(edges, hw_val, color='red', label='Herwig', where='post')
                ax.set_yscale('log')
            elif is_ratio:
                # Calculate ratio using 50 bins
                bin_edges = np.linspace(0, target_pt*2, 51)
                bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
                current_bin_width = (target_pt*2) / 50
                py_counts, _ = np.histogram(d['py_jets'], bins=bin_edges)
                hw_counts, _ = np.histogram(d['hw_jets'], bins=bin_edges)
                py_xsec = d['py_scale'] * 1e9 * py_counts / current_bin_width
                hw_xsec = d['hw_scale'] * 1e3 * hw_counts / current_bin_width
                ratio = np.divide(hw_xsec, py_xsec, out=np.zeros_like(hw_xsec), where=py_xsec!=0)


                ax.step(bin_edges[:-1], ratio, color='purple', where='post', linewidth=2)
                ax.axhline(1.0, color='black', linestyle='--')
                if np.nanmax(ratio) > 0:
                    ax.set_ylim(0, 2.0) #max(2.0, np.nanmax(ratio) * 1.1))
            else:
                # Handle Particle distributions (ROOT)
                py_vals, py_edges = d[f'py_{metric_key}'], d[f'py_{metric_key}_edges']
                hw_vals, hw_edges = d[f'hw_{metric_key}'], d[f'hw_{metric_key}_edges']
                if py_vals is not None and hw_vals is not None:
                    ax.step(py_edges[:-1], py_vals, color='blue', label='Pythia', where='mid')
                    ax.step(hw_edges[:-1], hw_vals, color='red', label='Herwig', where='mid')
                    ax.set_yscale('log')

            ax.set_title(f"Target pT = {target_pt} GeV")
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.legend()
            ax.grid(True, alpha=0.3)

        plt.suptitle(title, fontsize=16)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        filename = f"{metric_key}_{'xsec_ratio' if is_ratio else 'xsec' if is_xsec else 'dist'}_all.pdf"
        plt.savefig(os.path.join(output_dir, filename))
        plt.close()
        print(f"Saved combined plot: {filename}")

    # Execute plots
    plot_combined("particle_pt", "Generated Particle pT Comparison", "pT [GeV]", "Counts")
    plot_combined("particle_eta", "Generated Particle Eta Comparison", "Eta", "Counts")

    # Custom handler for jet pT distributions
    def plot_jet_dist_all():
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()
        for i, target_pt in enumerate(target_jet_pts):
            ax = axes[i]
            d = all_data[target_pt]
            if d['py_jets'] is not None and d['hw_jets'] is not None:
                ax.hist(d['py_jets'], bins=50, range=(0, target_pt*2), color='blue',
                         label='Pythia', histtype='step') #, density=True)
                ax.hist(d['hw_jets'], bins=50, range=(0, target_pt*2), color='red',
                         label='Herwig', histtype='step') #, density=True)
                ax.set_title(f"Target pT = {target_pt} GeV")
                ax.set_xlabel("Jet pT [GeV]")
                ax.set_ylabel("Counts") #"Normalized Counts")
                ax.legend()
                ax.grid(True, alpha=0.3)
        plt.suptitle("Jet pT Distribution Comparison", fontsize=16)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(os.path.join(output_dir, "jet_pt_dist_all.pdf"))
        plt.close()

    plot_jet_dist_all()
    plot_combined("jet_pt", "Jet pT Cross Section Comparison", "Jet pT [GeV]", "d$\sigma$/dpT [pb/GeV]", is_xsec=True)
    plot_combined("jet_pt", "Jet pT XSec Ratio (Herwig/Pythia)", "Jet pT [GeV]", "Ratio", is_ratio=True)

if __name__ == "__main__":
    run_comparison()
