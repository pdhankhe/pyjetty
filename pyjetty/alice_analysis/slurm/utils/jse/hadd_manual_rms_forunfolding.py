import ROOT
import argparse
import os
import gc
import shutil

# python3 -u hadd_manual_rms_forunfolding.py \
#     --dir /global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/57485566 \
#     --batch-size 200 --group-budget-mb 500 --resume --mem-limit-mb 6000

#     --filelist /global/.../rms/anchmc_subpath_filelist.txt \

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

hist_list = ['roounfold_response_1D_{}', 'jet_match_rec_pur_num_{}', 'jet_all_rec_pur_den_{}',
             'jet_match_gen_eff_num_{}', 'jet_all_gen_eff_den_{}', 
             '{}_full_ungroomed_roounfold_response', '{}_AA_roounfold_response',
             '{}_AB_roounfold_response', '{}_BB_roounfold_response', '{}_rad_roounfold_response', 
             '{}_full_ungroomed_reco', '{}_full_ungroomed_reco_unmatched', '{}_AA_reco', '{}_AA_reco_unmatched', 
             '{}_AB_reco', '{}_AB_reco_unmatched', '{}_BB_reco', '{}_BB_reco_unmatched',
             '{}_rad_reco', '{}_rad_reco_unmatched']

default_subpath_file = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/anchmc_subpath_filelist.txt"


# ----------------------------------------------------------------------
# memory instrumentation
# ----------------------------------------------------------------------
def rss_mb():
    try:
        cur = hwm = 0.0
        with open("/proc/self/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    cur = float(line.split()[1]) / 1024.0
                elif line.startswith("VmHWM:"):
                    hwm = float(line.split()[1]) / 1024.0
        return cur, hwm
    except IOError:
        import resource
        m = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
        return m, m


def print_mem(tag):
    cur, hwm = rss_mb()
    print("[mem] {:<45s} rss={:8.1f} MB  peak={:8.1f} MB".format(tag, cur, hwm), flush=True)


def obj_size_mb(obj):
    try:
        if isinstance(obj, ROOT.TH1):
            n = obj.GetNcells()
            cname = obj.ClassName()
            per = 8 if ('D' in cname) else (4 if ('F' in cname) else 8)
            size = n * per
            if obj.GetSumw2N() > 0:
                size += n * 8
            return size / 1024.0 ** 2
        if hasattr(obj, "Hresponse"):
            return 1.4 * obj_size_mb(obj.Hresponse())
    except Exception:
        pass
    return 0.0


# ----------------------------------------------------------------------
# cleanup helpers
# ----------------------------------------------------------------------
def safe_delete(obj):
    if obj is None:
        return
    try:
        ROOT.SetOwnership(obj, False)
        obj.Delete()
    except Exception:
        pass
    del obj


def close_file(fin):
    if fin is None:
        return
    try:
        fin.Close()
        ROOT.gROOT.GetListOfFiles().Remove(fin)
    except Exception:
        pass
    del fin


def open_file(fpath):
    if not os.path.exists(fpath):
        return None
    f = ROOT.TFile.Open(fpath, "READ")
    if not f or f.IsZombie():
        if f:
            close_file(f)
        return None
    return f


def sweep():
    gc.collect()


# ----------------------------------------------------------------------
# size scan and grouping
# ----------------------------------------------------------------------
def scan_sizes(sample_file, grooming_str, quiet=False):
    sizes = {}
    f = open_file(sample_file)
    if f is None:
        raise RuntimeError("cannot open sample file {}".format(sample_file))
    if not quiet:
        print("\n--- object sizes in {} ---".format(sample_file))
    total = 0.0
    for h in hist_list:
        hname = h.format(grooming_str)
        o = f.Get(hname)
        if not o:
            sizes[hname] = 0.0
            if not quiet:
                print("  {:<40s} MISSING".format(hname))
            continue
        s = obj_size_mb(o)
        sizes[hname] = s
        total += s
        if not quiet:
            print("  {:<40s} {:>18s}  {:8.1f} MB".format(hname, o.ClassName(), s))
        safe_delete(o)
    close_file(f)
    sweep()
    if not quiet:
        print("  {:<40s} {:>18s}  {:8.1f} MB  <-- one full set".format("TOTAL", "", total))
        print("---------------------------------------\n")
    return sizes


def build_groups(grooming_str, sizes, budget_mb):
    """Objects are merged group by group. Peak RAM ~ 2x the largest group."""
    groups, cur, cur_size = [], [], 0.0
    for h in hist_list:
        hname = h.format(grooming_str)
        s = sizes.get(hname, 0.0)
        if s >= budget_mb:
            groups.append([hname])
            continue
        if cur and cur_size + s > budget_mb:
            groups.append(cur)
            cur, cur_size = [], 0.0
        cur.append(hname)
        cur_size += s
    if cur:
        groups.append(cur)
    print("merge groups (budget {:.0f} MB):".format(budget_mb))
    for g in groups:
        print("   {:8.1f} MB : {}".format(sum(sizes.get(hname, 0.0) for hname in g), ", ".join(g)))
    worst = max(sum(sizes.get(hname, 0.0) for hname in g) for g in groups)
    print("   estimated peak above baseline: ~{:.0f} MB\n".format(2 * worst))
    return groups


# ----------------------------------------------------------------------
# merging
# ----------------------------------------------------------------------
def output_is_complete(path, grooming_str):
    if not os.path.exists(path):
        return False
    f = open_file(path)
    if f is None:
        return False
    keys = set(k.GetName() for k in f.GetListOfKeys())
    close_file(f)
    return all(h.format(grooming_str) in keys for h in hist_list)


def merge_batch(file_paths, out_path, groups, bad_files, verbose_mem=False):
    """One group of objects at a time; memory is independent of len(file_paths)."""
    good = []
    for fpath in file_paths:
        f = open_file(fpath)
        if f is None:
            bad_files.append(fpath)
            continue
        good.append(fpath)
        close_file(f)

    if not good:
        return None, 0

    fout = ROOT.TFile(out_path, "RECREATE")
    n_written = 0

    for ig, group in enumerate(groups):
        acc = {}
        for fpath in good:
            fin = open_file(fpath)
            if fin is None:
                continue
            for h_name in group:
                obj = fin.Get(h_name)
                if not obj:
                    continue
                if h_name in acc:
                    acc[h_name].Add(obj)
                else:
                    c = obj.Clone(h_name)
                    ROOT.SetOwnership(c, False)
                    if isinstance(c, ROOT.TH1):
                        c.SetDirectory(0)
                    acc[h_name] = c
                safe_delete(obj)
            close_file(fin)

        fout.cd()
        for h_name in list(acc.keys()):
            acc[h_name].Write(h_name)
            safe_delete(acc.pop(h_name))
            n_written += 1
        acc.clear()
        sweep()
        if verbose_mem:
            print_mem("    group {}/{} done".format(ig + 1, len(groups)))

    fout.Close()
    del fout
    sweep()
    return (out_path, len(good)) if n_written else (None, 0)


def hierarchical_merge(indir, subpaths, grooming_str, batch_size, tmpdir, keep_tmp,
                       groups, verbose_mem, mem_limit, resume):
    os.makedirs(tmpdir, exist_ok=True)

    current = [os.path.join(indir, sp, "response.root") for sp in subpaths]
    bad_files = []
    level = 0
    previous_level_files = []

    print("starting from {} input files, batch size {}".format(len(current), batch_size))
    print_mem("start")

    while len(current) > 1:
        next_level = []
        n_batches = (len(current) + batch_size - 1) // batch_size
        print("\n=== level {}: merging {} files in {} batches ===".format(level, len(current), n_batches))

        for ib in range(n_batches):
            batch = current[ib * batch_size:(ib + 1) * batch_size]
            out_path = os.path.join(tmpdir, "level{}_{:05d}.root".format(level, ib))

            if resume and output_is_complete(out_path, grooming_str):
                print("  batch {:>4d}/{:<4d}  already done, skipping".format(ib + 1, n_batches), flush=True)
                next_level.append(out_path)
                continue

            res, n_used = merge_batch(batch, out_path, groups, bad_files, verbose_mem)
            print("  batch {:>4d}/{:<4d}  {:>4d}/{:<4d} files -> {}".format(
                ib + 1, n_batches, n_used, len(batch),
                os.path.basename(out_path) if res else "EMPTY"), flush=True)
            if res:
                next_level.append(res)

            sweep()
            print_mem("after level {} batch {}".format(level, ib + 1))
            cur, _ = rss_mb()
            if mem_limit and cur > mem_limit:
                raise MemoryError("RSS {:.0f} MB exceeds --mem-limit-mb {:.0f}".format(cur, mem_limit))

        if not keep_tmp:
            for f in previous_level_files:
                if f not in next_level:
                    try:
                        os.remove(f)
                    except OSError:
                        pass

        if not next_level:
            raise RuntimeError("level {} produced no output".format(level))

        previous_level_files = next_level
        current = next_level
        level += 1

    final_path = os.path.join(indir, f"response_{grooming_str}_forunfolding_merged.root")
    if os.path.dirname(os.path.abspath(current[0])) == os.path.abspath(tmpdir):
        shutil.move(current[0], final_path)
    else:
        shutil.copy(current[0], final_path)

    if not keep_tmp and os.path.isdir(tmpdir) and not os.listdir(tmpdir):
        os.rmdir(tmpdir)

    print("\nwritten to {}".format(final_path))
    print_mem("final")
    return final_path, bad_files


def load_roounfold():
    print("LOADING ROOUNFOLD")
    roo = False
    while roo is False:
        roo = ROOT.gSystem.Load("libRooUnfold")
    print(roo)
    roo2 = False
    while roo2 is False:
        roo2 = ROOT.gSystem.Load("libRooUnfold.so")
    print(roo2)


def read_subpaths(filelist):
    out = []
    with open(filelist) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                out.append(line.strip('/'))
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', default=None)
    parser.add_argument('--filelist', default=default_subpath_file)
    parser.add_argument('--groomingstr', default="groomed", help="choose either 'groomed' or 'ungroomed'")
    parser.add_argument('--batch-size', type=int, default=200,
                        help="files per mini-merge; does NOT affect memory, only I/O")
    parser.add_argument('--group-budget-mb', type=float, default=500.0,
                        help="max MB of accumulators held at once; peak RAM is ~2x this")
    parser.add_argument('--tmpdir', default=None)
    parser.add_argument('--keep-tmp', action='store_true')
    parser.add_argument('--resume', action='store_true',
                        help="skip intermediate files that already exist and look complete")
    parser.add_argument('--verbose-mem', action='store_true')
    parser.add_argument('--mem-limit-mb', type=float, default=0.0)
    parser.add_argument('--scan-only', action='store_true')
    flags = parser.parse_args()

    load_roounfold()
    subpaths = read_subpaths(flags.filelist)
    print_mem("after libRooUnfold load")

    sizes = scan_sizes(os.path.join(flags.dir, subpaths[0], "response.root"), flags.groomingstr)
    if flags.scan_only:
        raise SystemExit(0)
    groups = build_groups(flags.groomingstr, sizes, flags.group_budget_mb)

    tmpdir = flags.tmpdir if flags.tmpdir else os.path.join(flags.dir, "tmp_merge")
    final_path, bad_files = hierarchical_merge(flags.dir, subpaths, flags.groomingstr, flags.batch_size, tmpdir,
                                               flags.keep_tmp, groups, flags.verbose_mem,
                                               flags.mem_limit_mb, flags.resume)

    print("\n{} files missing / unreadable / empty".format(len(bad_files)))
    if bad_files:
        log = os.path.join(flags.dir, "missing_files.txt")
        with open(log, "w") as f:
            for b in bad_files:
                f.write(b + "\n")
        print("list written to {}".format(log))

        
        
# import ROOT
# import argparse
# import os
# import shutil

# # python3 -u hadd_manual_rms_forunfolding.py \
# #     --dir /global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/57485566 \
# #     --batch-size 15

# #     --filelist /global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/anchmc_subpath_filelist.txt \

# ROOT.gROOT.SetBatch(True)
# # Detach histograms from the TFile they are read from, so they survive fin.Close()
# ROOT.TH1.AddDirectory(False)

# hist_list = ['roounfold_response_1D', 'jet_match_rec_pur_num_groomed', 'jet_all_rec_pur_den_groomed',
#              'jet_match_gen_eff_num_groomed', 'jet_all_gen_eff_den_groomed', 'AA_roounfold_response',
#              'AB_roounfold_response', 'BB_roounfold_response', 'rad_roounfold_response', 'AA_reco',
#              'AA_reco_unmatched', 'AB_reco', 'AB_reco_unmatched', 'BB_reco', 'BB_reco_unmatched',
#              'rad_reco', 'rad_reco_unmatched']

# default_subpath_file = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/anchmc_subpath_filelist.txt"


# def load_roounfold():
#     print("LOADING ROOUNFOLD")
#     roo = False
#     while roo is False:
#         roo = ROOT.gSystem.Load("libRooUnfold")
#     print(roo)
#     roo2 = False
#     while roo2 is False:
#         roo2 = ROOT.gSystem.Load("libRooUnfold.so")
#     print(roo2)


# def read_subpaths(filelist):
#     """Read lines like '2/559385/338' and return them as a clean list."""
#     subpaths = []
#     with open(filelist) as f:
#         for line in f:
#             line = line.strip()
#             if not line or line.startswith('#'):
#                 continue
#             subpaths.append(line.strip('/'))
#     return subpaths


# def detach(obj):
#     """Make sure the object is owned by Python/us and not by the input file."""
#     if isinstance(obj, ROOT.TH1):
#         obj.SetDirectory(0)
#     ROOT.SetOwnership(obj, True)
#     return obj


# def merge_batch(file_paths, out_path, bad_files):
#     """Merge one batch of files into out_path. Returns (out_path or None, n_files_used)."""
#     merged = {}
#     n_used = 0

#     for fpath in file_paths:
#         fin = None
#         try:
#             if not os.path.exists(fpath):
#                 bad_files.append(fpath)
#                 continue

#             fin = ROOT.TFile.Open(fpath)
#             if not fin or fin.IsZombie():
#                 bad_files.append(fpath)
#                 continue

#             found_something = False
#             for h_name in hist_list:
#                 obj = fin.Get(h_name)
#                 if not obj:
#                     continue
#                 if h_name in merged:
#                     merged[h_name].Add(obj)
#                 else:
#                     merged[h_name] = detach(obj.Clone(h_name))
#                 found_something = True

#             if found_something:
#                 n_used += 1
#             else:
#                 bad_files.append(fpath)

#         except Exception as e:
#             print("  skipping {} ({})".format(fpath, e))
#             bad_files.append(fpath)
#         finally:
#             if fin:
#                 fin.Close()

#     if not merged:
#         return None, 0

#     fout = ROOT.TFile(out_path, "RECREATE")
#     fout.cd()
#     for h_name, obj in merged.items():
#         obj.Write(h_name)
#     fout.Close()

#     merged.clear()
#     return out_path, n_used


# def hierarchical_merge(indir, subpaths, batch_size, tmpdir, keep_tmp):
#     """Merge in successive levels: files -> mini-merges -> merges of mini-merges -> ... -> 1 file."""
#     os.makedirs(tmpdir, exist_ok=True)

#     current = [os.path.join(indir, sp, "response.root") for sp in subpaths]
#     bad_files = []
#     level = 0
#     previous_level_files = []   # intermediates that can be deleted once consumed

#     print("starting from {} input files, batch size {}".format(len(current), batch_size))

#     while len(current) > 1:
#         next_level = []
#         n_batches = (len(current) + batch_size - 1) // batch_size
#         print("\n=== level {}: merging {} files in {} batches ===".format(level, len(current), n_batches))

#         for ib in range(n_batches):
#             batch = current[ib * batch_size:(ib + 1) * batch_size]
#             out_path = os.path.join(tmpdir, "level{}_{:05d}.root".format(level, ib))
#             res, n_used = merge_batch(batch, out_path, bad_files)
#             print("  batch {}/{}: {} / {} files merged -> {}".format(
#                 ib + 1, n_batches, n_used, len(batch), os.path.basename(out_path) if res else "EMPTY"))
#             if res:
#                 next_level.append(res)

#         # clean up the intermediates we just consumed
#         if not keep_tmp:
#             for f in previous_level_files:
#                 try:
#                     os.remove(f)
#                 except OSError:
#                     pass

#         if not next_level:
#             raise RuntimeError("level {} produced no output - nothing could be merged".format(level))

#         previous_level_files = next_level
#         current = next_level
#         level += 1

#     if not current:
#         raise RuntimeError("no files to merge")

#     final_path = os.path.join(indir, "response_forunfolding_merged.root")
#     if os.path.dirname(os.path.abspath(current[0])) == os.path.abspath(tmpdir):
#         shutil.move(current[0], final_path)
#     else:
#         # only one input file to begin with
#         shutil.copy(current[0], final_path)

#     if not keep_tmp and os.path.isdir(tmpdir) and not os.listdir(tmpdir):
#         os.rmdir(tmpdir)

#     print("\nwritten to {}".format(final_path))
#     print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
#     return final_path, bad_files


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--dir', default=None, help="base directory containing the subpaths")
#     parser.add_argument('--batch-size', type=int, default=20, help="how many files per mini-merge")
#     parser.add_argument('--tmpdir', default=None, help="where to put intermediate merges (default <dir>/tmp_merge)")
#     parser.add_argument('--keep-tmp', action='store_true', help="do not delete intermediate merge files")
#     parser.add_argument('--filelist', default='/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/anchmc_subpath_filelist.txt', help="text file with one subpath per line")
#     flags = parser.parse_args()

#     load_roounfold()

#     subpaths = read_subpaths(flags.filelist)
#     tmpdir = flags.tmpdir if flags.tmpdir else os.path.join(flags.dir, "tmp_merge")

#     final_path, bad_files = hierarchical_merge(flags.dir, subpaths, flags.batch_size, tmpdir, flags.keep_tmp)

#     print("\n{} files were missing / unreadable / empty".format(len(bad_files)))
#     if bad_files:
#         log = os.path.join(flags.dir, "missing_files.txt")
#         with open(log, "w") as f:
#             for b in bad_files:
#                 f.write(b + "\n")
#         print("list written to {}".format(log))