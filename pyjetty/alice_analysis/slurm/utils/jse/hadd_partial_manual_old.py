import ROOT
import argparse
import os
import gc

# python3 -u hadd_partial_manual_old.py --dir /global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/57259276 --grooming groomed --nbatches 4 --budget-mb 500 --keep-tmp
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)   # objects are detached from the input file

GLOBAL_HIST_LIST = [
    'resp_jetpt',
    'resp6_AA', 'resp6_AB', 'resp6_BB', 'resp6_rad',
    'jet_match_gen_eff_num_groomed', 'jet_all_gen_eff_den_groomed', 'jet_match_rec_pur_num_groomed', 'jet_all_rec_pur_den_groomed',
    'lund_matched_gen', 'lund_matched_rec', 'lund_all_gen', 'lund_all_rec',
    'pair_match_gen_eff_num_full_ungroomed', 'pair_all_gen_eff_den_full_ungroomed', 'pair_match_rec_pur_num_full_ungroomed', 'pair_all_rec_pur_den_full_ungroomed',
    'pair_match_gen_eff_num_AA', 'pair_all_gen_eff_den_AA', 'pair_match_rec_pur_num_AA', 'pair_all_rec_pur_den_AA',
    'pair_match_gen_eff_num_AB', 'pair_all_gen_eff_den_AB', 'pair_match_rec_pur_num_AB', 'pair_all_rec_pur_den_AB',
    'pair_match_gen_eff_num_BB', 'pair_all_gen_eff_den_BB', 'pair_match_rec_pur_num_BB', 'pair_all_rec_pur_den_BB',
    'pair_match_gen_eff_num_rad', 'pair_all_gen_eff_den_rad', 'pair_match_rec_pur_num_rad', 'pair_all_rec_pur_den_rad',
    'summary_efficiencies'
]


default_subpath_file = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/anchmc_subpath_filelist.txt"


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------
def rss_mb():
    cur = hwm = 0.0
    try:
        with open("/proc/self/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    cur = float(line.split()[1]) / 1024.0
                elif line.startswith("VmHWM:"):
                    hwm = float(line.split()[1]) / 1024.0
    except IOError:
        pass
    return cur, hwm


def print_mem(tag):
    cur, hwm = rss_mb()
    print("[mem] {:<50s} rss={:8.1f} MB  peak={:8.1f} MB".format(tag, cur, hwm), flush=True)


def free(obj):
    """Actually release the C++ memory. This is what was missing."""
    if obj is None:
        return
    try:
        ROOT.SetOwnership(obj, False)
        obj.Delete()
    except Exception:
        pass


def open_read(path):
    if not os.path.exists(path):
        return None
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        if f:
            f.Close()
        return None
    return f


def close_read(f):
    try:
        f.Close()
        ROOT.gROOT.GetListOfFiles().Remove(f)
    except Exception:
        pass


def read_subpaths(filelist):
    out = []
    with open(filelist) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                out.append(line.strip('/'))
    return out


def full_name(h_name, grooming_str):
    return h_name.format(grooming_str) if "{}" in h_name else h_name


def obj_size_mb(obj):
    try:
        if isinstance(obj, ROOT.TH1):
            n = obj.GetNcells()
            per = 8 if 'D' in obj.ClassName() else 4
            s = n * per
            if obj.GetSumw2N() > 0:
                s += n * 8
            return s / 1024.0 ** 2
        if hasattr(obj, "Hresponse"):
            return 1.4 * obj_size_mb(obj.Hresponse())
    except Exception:
        pass
    return 0.0


def scan_and_group(sample_path, names, budget_mb):
    """Measure each object once, then bin the names into groups under budget_mb."""
    f = open_read(sample_path)
    if f is None:
        raise RuntimeError("cannot open sample file " + sample_path)
    sizes = {}
    print("\n--- object sizes in {} ---".format(sample_path))
    for n in names:
        o = f.Get(n)
        if not o:
            sizes[n] = 0.0
            print("  {:<45s} MISSING".format(n))
            continue
        sizes[n] = obj_size_mb(o)
        print("  {:<45s} {:>18s} {:8.1f} MB".format(n, o.ClassName(), sizes[n]))
        free(o)
    close_read(f)
    gc.collect()
    print("  {:<45s} {:>18s} {:8.1f} MB  <-- one full set".format("TOTAL", "", sum(sizes.values())))

    groups, cur, cur_size = [], [], 0.0
    for n in names:
        s = sizes.get(n, 0.0)
        if s >= budget_mb:
            groups.append([n])
            continue
        if cur and cur_size + s > budget_mb:
            groups.append(cur)
            cur, cur_size = [], 0.0
        cur.append(n)
        cur_size += s
    if cur:
        groups.append(cur)

    print("\nmerge groups (budget {:.0f} MB):".format(budget_mb))
    for g in groups:
        print("  {:8.1f} MB : {}".format(sum(sizes.get(n, 0.0) for n in g), ", ".join(g)))
    print("estimated peak above baseline ~{:.0f} MB\n".format(
        2 * max(sum(sizes.get(n, 0.0) for n in g) for g in groups)))
    return groups


# ----------------------------------------------------------------------
# core merge: one group of objects at a time
# ----------------------------------------------------------------------
def merge_group(file_paths, names, fout, missing, label=""):
    acc = {}
    for i, fpath in enumerate(file_paths):
        fin = open_read(fpath)
        if fin is None:
            missing.add(fpath)
            continue
        for n in names:
            obj = fin.Get(n)
            if not obj:
                continue
            if n in acc:
                try:
                    acc[n].Add(obj)
                except Exception as e:
                    print("  cannot Add {}: {}".format(n, e))
                free(obj)                      # <-- release every object read
            else:
                if isinstance(obj, ROOT.TH1):
                    obj.SetDirectory(0)
                ROOT.SetOwnership(obj, False)
                acc[n] = obj                   # adopt as accumulator, do not Clone
        close_read(fin)
        del fin
        if (i + 1) % 500 == 0:
            gc.collect()
            print_mem("  {}{} files".format(label, i + 1))

    fout.cd()
    for n in list(acc.keys()):
        o = acc.pop(n)
        o.Write(n)
        free(o)
    acc.clear()
    gc.collect()


def merge_stage(file_paths, groups, out_path, missing, label):
    fout = ROOT.TFile(out_path, "RECREATE")
    for ig, names in enumerate(groups):
        merge_group(file_paths, names, fout, missing, label="{}g{} ".format(label, ig + 1))
        print_mem("{}group {}/{} written".format(label, ig + 1, len(groups)))
    fout.Close()
    del fout
    gc.collect()
    return out_path


def split_into_batches(items, nbatches):
    nbatches = max(1, min(nbatches, len(items)))
    n, start, chunks = len(items), 0, []
    for i in range(nbatches):
        size = n // nbatches + (1 if i < n % nbatches else 0)
        chunks.append(items[start:start + size])
        start += size
    return chunks


# ----------------------------------------------------------------------
# ratios, computed at the very end from the finished file
# ----------------------------------------------------------------------
def write_ratios(path, grooming_str):
    g = grooming_str
    recipes = [(f'jet_match_gen_eff_num_{g}', f'jet_all_gen_eff_den_{g}', f'jet_efficiency_{g}_new'),
               (f'jet_match_rec_pur_num_{g}', f'jet_all_rec_pur_den_{g}', f'jet_purity_{g}_new')]
    if g == "groomed":
        recipes += [('lund_matched_gen', 'lund_all_gen', 'lund_split_efficiency_new'),
                    ('lund_matched_rec', 'lund_all_rec', 'lund_split_purity_new')]
        for tag in ['full_ungroomed', 'AA', 'AB', 'BB', 'rad']:
            recipes += [(f'pair_match_gen_eff_num_{tag}', f'pair_all_gen_eff_den_{tag}', f'pair_efficiency_{tag}_new'),
                        (f'pair_match_rec_pur_num_{tag}', f'pair_all_rec_pur_den_{tag}', f'pair_purity_{tag}_new')]

    f = ROOT.TFile.Open(path, "UPDATE")
    for num, den, newname in recipes:
        hn, hd = f.Get(num), f.Get(den)
        if not hn or not hd:
            print("  ratio {}: missing {}".format(newname, num if not hn else den))
            continue
        r = hn.Clone(newname)
        r.SetDirectory(0)
        r.Divide(hd)
        f.cd()
        r.Write(newname, ROOT.TObject.kOverwrite)
        free(r)
        free(hn)
        free(hd)
    f.Close()
    gc.collect()


# ----------------------------------------------------------------------
def finishing_merge(dir_path, grooming_str, nbatches=2, budget_mb=500.0, keep_tmp=False):
    print("LOADING ROOUNFOLD")
    roo = False
    while not roo:
        roo = ROOT.gSystem.Load("libRooUnfold")
    roo2 = False
    while not roo2:
        roo2 = ROOT.gSystem.Load("libRooUnfold.so")
    print_mem("after libRooUnfold load")

    subpaths = read_subpaths(default_subpath_file)
    file_paths = ["{}/{}/response.root".format(dir_path, sp) for sp in subpaths]

    raw_list = GLOBAL_HIST_LIST
    names = [full_name(h, grooming_str) for h in raw_list]

    groups = scan_and_group(file_paths[0], names, budget_mb)
    missing = set()

    batches = split_into_batches(file_paths, nbatches)
    print("merging {} files for '{}' in {} batches of sizes {}".format(
        len(file_paths), grooming_str, len(batches), [len(b) for b in batches]), flush=True)

    partials = []
    tmpdir = "{}/tmp".format(dir_path)
    os.makedirs(tmpdir, exist_ok=True)
    for ib, batch in enumerate(batches):
        tmp = "{}/tmp/tmp_{}_batch{}.root".format(dir_path, grooming_str, ib)
        if keep_tmp and os.path.exists(tmp):
            print("=== batch {}/{} already exists, skipping ===".format(ib + 1, len(batches)))
            partials.append(tmp)
            continue
        print("=== batch {}/{} ({} files) ===".format(ib + 1, len(batches), len(batch)), flush=True)
        partials.append(merge_stage(batch, groups, tmp, missing, "b{} ".format(ib + 1)))

    out_path = "{}/response_merged_partial.root".format(dir_path, grooming_str)
    print("=== combining {} partial files ===".format(len(partials)), flush=True)
    merge_stage(partials, groups, out_path, missing, "final ")

    print("=== computing efficiencies and purities ===", flush=True)
    write_ratios(out_path, grooming_str)

    if not keep_tmp:
        for p in partials:
            try:
                os.remove(p)
            except OSError:
                pass
    if not keep_tmp and os.path.isdir(tmpdir) and not os.listdir(tmpdir):
        os.rmdir(tmpdir)

    nonexistent.extend(sorted(missing))
    print("written to {}".format(out_path))
    print("{} files were missing or unreadable".format(len(missing)))
    print_mem("done with " + grooming_str)
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', default=None)
    parser.add_argument('--nbatches', type=int, default=2)
    parser.add_argument('--budget-mb', type=float, default=500.0,
                        help="max MB of accumulators held at once; peak is ~2x this")
    parser.add_argument('--grooming', default="both", choices=["groomed", "ungroomed", "both"])
    parser.add_argument('--keep-tmp', action='store_true',
                        help="keep per-batch files and reuse them on a re-run")
    flags = parser.parse_args()

    nonexistent = []
    if flags.grooming in ("groomed", "both"):
        finishing_merge(flags.dir, "groomed", flags.nbatches, flags.budget_mb, flags.keep_tmp)
    if flags.grooming in ("ungroomed", "both"):
        finishing_merge(flags.dir, "ungroomed", flags.nbatches, flags.budget_mb, flags.keep_tmp)