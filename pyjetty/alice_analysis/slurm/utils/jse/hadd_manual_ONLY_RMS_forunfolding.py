import ROOT
import argparse
import os
import gc
import shutil

# python3 -u hadd_manual_rms_forunfolding.py \
#     --dir /global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/58919321 \
#     --batch-size 20 --group-budget-mb 500 --resume --mem-limit-mb 6000
#
#     optional:
#     --filelist /global/.../rms/anchmc_subpath_filelist.txt
#     --roounfold-include /dir/containing/RooUnfoldResponse.h

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

hist_list = ['roounfold_response_1D_{}', 'roounfold_response_1D_{}_oldversion', 
             '{}_full_ungroomed_roounfold_response', '{}_AA_roounfold_response',
             '{}_AB_roounfold_response', '{}_BB_roounfold_response', '{}_rad_roounfold_response']

default_subpath_file = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/anchmc_subpath_filelist.txt"

COMPLETE_MARKER = "__merge_complete__"
SKIP_FILE_NAME = "skip_pairs.txt"


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


# ----------------------------------------------------------------------
# RooUnfold loading (library + headers so PyROOT knows the methods)
# ----------------------------------------------------------------------
_RU_CLASS = None


def _loaded_lib_path(name="libRooUnfold"):
    """Path of the libRooUnfold actually mapped into this process."""
    try:
        with open("/proc/self/maps") as f:
            for line in f:
                parts = line.split()
                if not parts:
                    continue
                p = parts[-1]
                if os.path.basename(p).startswith(name) and ".so" in p:
                    return os.path.realpath(p)
    except IOError:
        pass
    try:
        p = str(ROOT.gSystem.DynamicPathName(name, True))
        if p:
            return os.path.realpath(p)
    except Exception:
        pass
    return None


def _find_header_dir(libpath, extra_dirs):
    header = "RooUnfoldResponse.h"
    cands = [d for d in extra_dirs if d]
    cands += [d for d in os.environ.get("ROOT_INCLUDE_PATH", "").split(":") if d]
    if libpath:
        libdir = os.path.dirname(libpath)
        prefix = os.path.dirname(libdir)
        cands += [libdir,
                  os.path.join(prefix, "include"),
                  os.path.join(prefix, "include", "RooUnfold"),
                  os.path.join(prefix, "src"),       # RooUnfold 1.x keeps headers in src/
                  prefix]
    for d in cands:
        if os.path.isfile(os.path.join(d, header)):
            return d, cands
    # last resort: bounded search under the install prefix
    if libpath:
        prefix = os.path.dirname(os.path.dirname(libpath))
        base = prefix.rstrip("/").count("/")
        for root, dirs, files in os.walk(prefix):
            if header in files:
                return root, cands
            if root.count("/") - base >= 4:
                dirs[:] = []
    return None, cands


def load_roounfold(include_dir=None):
    global _RU_CLASS
    print("LOADING ROOUNFOLD")
    rc = ROOT.gSystem.Load("libRooUnfold")
    print("  gSystem.Load('libRooUnfold') -> {}".format(rc))
    if rc < 0:
        raise RuntimeError("could not load libRooUnfold (is its lib dir in LD_LIBRARY_PATH?)")

    lib = _loaded_lib_path()
    print("  library : {}".format(lib))

    hdr_dir, tried = _find_header_dir(lib, [include_dir])
    if hdr_dir is None:
        raise RuntimeError("could not find RooUnfoldResponse.h; searched:\n    " + "\n    ".join(tried) +
                           "\n  pass --roounfold-include /dir/containing/RooUnfoldResponse.h")
    print("  headers : {}".format(hdr_dir))

    ROOT.gInterpreter.AddIncludePath(hdr_dir)
    if not ROOT.gInterpreter.Declare('#include "RooUnfoldResponse.h"'):
        raise RuntimeError("interpreter failed to parse RooUnfoldResponse.h from {}".format(hdr_dir))

    cls = ROOT.RooUnfoldResponse
    missing = [m for m in ("Hresponse", "Hmeasured", "Htruth", "Add") if not hasattr(cls, m)]
    if missing:
        raise RuntimeError("RooUnfoldResponse is still missing methods {} after loading headers".format(missing))
    _RU_CLASS = cls
    print("  RooUnfoldResponse bound OK")


def get_ru_class():
    if _RU_CLASS is None:
        raise RuntimeError("load_roounfold() must be called first")
    return _RU_CLASS


# ----------------------------------------------------------------------
# object helpers
# ----------------------------------------------------------------------
def is_response(obj):
    """Detect by C++ class name, not by Python attributes (the proxy may be a bare TObject)."""
    try:
        return (not isinstance(obj, ROOT.TH1)) and obj.ClassName().startswith("RooUnfoldResponse")
    except Exception:
        return False


def take_ownership(obj):
    """Give Python ownership of an object read from a file and return a correctly typed proxy.
    ALWAYS use the return value: for responses it may be a different proxy object."""
    if isinstance(obj, ROOT.TH1):
        obj.SetDirectory(0)
        ROOT.SetOwnership(obj, True)
        return obj
    if is_response(obj):
        if hasattr(obj, "Hresponse"):
            ROOT.SetOwnership(obj, True)
            return obj
        # re-bind the same C++ pointer with the proper RooUnfoldResponse type;
        # transfer ownership to the new proxy so the object is deleted exactly once
        cast = ROOT.BindObject(ROOT.addressof(obj), get_ru_class())
        ROOT.SetOwnership(obj, False)
        ROOT.SetOwnership(cast, True)
        return cast
    ROOT.SetOwnership(obj, True)
    return obj


def hist_mb(h):
    if not h:
        return 0.0
    n = h.GetNcells()
    cname = h.ClassName()
    per = 4 if cname.endswith('F') else (2 if cname.endswith('S') else (1 if cname.endswith('C') else 8))
    size = n * per
    if h.GetSumw2N() > 0:
        size += n * 8
    return size / 1024.0 ** 2


def obj_size_mb(obj):
    if isinstance(obj, ROOT.TH1):
        return hist_mb(obj)
    if is_response(obj):
        total = 0.0
        for m in ("Hresponse", "Hmeasured", "Htruth", "Hfakes"):
            f = getattr(obj, m, None)
            if f is not None:
                total += hist_mb(f())
        return total
    return 0.0


def _axis_sig(ax):
    return (ax.GetNbins(), round(ax.GetXmin(), 9), round(ax.GetXmax(), 9))


def _hist_sig(h):
    if not h:
        return None
    return (h.GetDimension(),) + _axis_sig(h.GetXaxis()) + _axis_sig(h.GetYaxis()) + _axis_sig(h.GetZaxis())


def obj_signature(obj):
    """Binning fingerprint; objects with different signatures are not added together."""
    if isinstance(obj, ROOT.TH1):
        return ("TH1", _hist_sig(obj))
    if is_response(obj):
        hr, hm, ht = obj.Hresponse(), obj.Hmeasured(), obj.Htruth()
        if not hr or not hm or not ht:
            raise RuntimeError("RooUnfoldResponse has null internal histogram(s)")
        return ("RU", _hist_sig(hm), _hist_sig(ht), _hist_sig(hr))
    raise TypeError("unsupported class {}".format(obj.ClassName()))


def add_into(acc_obj, obj):
    if isinstance(acc_obj, ROOT.TH1):
        if not acc_obj.Add(obj):
            raise RuntimeError("TH1::Add returned false")
    elif is_response(acc_obj):
        acc_obj.Add(obj)
    else:
        raise TypeError("don't know how to add {}".format(acc_obj.ClassName()))


def short(e):
    s = str(e).strip().splitlines()
    return s[0][:200] if s else type(e).__name__


# ----------------------------------------------------------------------
# file helpers
# ----------------------------------------------------------------------
def close_file(fin):
    if fin is None:
        return
    try:
        fin.Close()
        ROOT.gROOT.GetListOfFiles().Remove(fin)
    except Exception:
        pass


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
        o = take_ownership(o)
        s = obj_size_mb(o)
        sizes[hname] = s
        total += s
        if not quiet:
            print("  {:<40s} {:>28s}  {:8.1f} MB".format(hname, o.ClassName(), s))
        del o
    close_file(f)
    sweep()
    if not quiet:
        print("  {:<40s} {:>28s}  {:8.1f} MB  <-- one full set".format("TOTAL", "", total))
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
        print("   {:8.1f} MB : {}".format(sum(sizes.get(n, 0.0) for n in g), ", ".join(g)))
    worst = max(sum(sizes.get(n, 0.0) for n in g) for g in groups)
    print("   estimated peak above baseline: ~{:.0f} MB\n".format(2 * worst))
    return groups


# ----------------------------------------------------------------------
# merging
# ----------------------------------------------------------------------
def output_is_complete(path):
    f = open_file(path)
    if f is None:
        return False
    ok = bool(f.GetListOfKeys().FindObject(COMPLETE_MARKER))
    close_file(f)
    return ok


def merge_batch(file_paths, out_path, groups, bad_files, problems, skip_pairs, level, verbose_mem=False):
    """One group of objects at a time; each (file, histogram) is merged independently."""
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
    if not fout or fout.IsZombie():
        raise RuntimeError("cannot create {}".format(out_path))
    n_written = 0

    # breadcrumb: records the (file, histogram) currently being touched.
    # flush() hands it to the kernel, so it survives a segfault of this process.
    breadcrumb_path = out_path + ".breadcrumb"
    bf = open(breadcrumb_path, "w")

    def crumb(text):
        bf.seek(0)
        bf.truncate()
        bf.write(text)
        bf.flush()

    try:
        for ig, group in enumerate(groups):
            acc, sigs = {}, {}

            for fpath in good:
                if verbose_mem:
                    print_mem("    pre  {}".format(fpath[-45:]))

                fin = open_file(fpath)
                if fin is None:
                    problems.append((level, fpath, "*", "could not reopen file"))
                    continue

                try:
                    for h_name in group:
                        if (fpath, h_name) in skip_pairs:
                            continue
                        crumb(fpath + "\t" + h_name)
                        obj = None
                        try:
                            obj = fin.Get(h_name)
                            if not obj:
                                obj = None
                                if level == 0:
                                    problems.append((level, fpath, h_name, "missing"))
                                continue
                            obj = take_ownership(obj)
                            sig = obj_signature(obj)

                            if h_name not in acc:
                                # first good copy becomes the accumulator (no clone needed)
                                acc[h_name] = obj
                                sigs[h_name] = sig
                            else:
                                if sig != sigs[h_name]:
                                    raise RuntimeError("binning mismatch: {} vs accumulator {}".format(
                                        sig, sigs[h_name]))
                                add_into(acc[h_name], obj)
                        except Exception as e:
                            msg = short(e)
                            problems.append((level, fpath, h_name, msg))
                            print("  !!! skipping {} from {}: {}".format(h_name, fpath, msg), flush=True)
                        finally:
                            obj = None  # drops our reference; non-accumulator copies are freed here
                finally:
                    close_file(fin)

                if verbose_mem:
                    print_mem("    post")

            fout.cd()
            for h_name in group:
                if h_name not in acc:
                    continue
                o = acc.pop(h_name)
                try:
                    o.Write(h_name)
                    n_written += 1
                except Exception as e:
                    problems.append((level, out_path, h_name, "write failed: " + short(e)))
                    print("  !!! could not write {}: {}".format(h_name, short(e)), flush=True)
                del o
            acc.clear()
            sweep()
            if verbose_mem:
                print_mem("    group {}/{} done".format(ig + 1, len(groups)))

        if n_written:
            fout.cd()
            ROOT.TNamed(COMPLETE_MARKER, str(len(good))).Write()
    finally:
        bf.close()
        if os.path.exists(breadcrumb_path):
            os.remove(breadcrumb_path)
        fout.Close()
        sweep()

    if not n_written:
        try:
            os.remove(out_path)
        except OSError:
            pass
        return None, 0
    return out_path, len(good)


def load_skip_pairs(tmpdir):
    pairs = set()
    path = os.path.join(tmpdir, SKIP_FILE_NAME)
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                line = line.rstrip("\n")
                if "\t" in line:
                    fp, hn = line.split("\t", 1)
                    pairs.add((fp, hn))
    if pairs:
        print("loaded {} (file, histogram) pairs to skip from previous crashes".format(len(pairs)))
    return pairs


def check_crash_breadcrumb(out_path, skip_pairs, tmpdir):
    """A leftover breadcrumb means the previous run segfaulted on that (file, histogram)."""
    bc = out_path + ".breadcrumb"
    if not os.path.exists(bc):
        return
    with open(bc) as f:
        txt = f.read().strip()
    os.remove(bc)
    if "\t" not in txt:
        return
    fp, hn = txt.split("\t", 1)
    print("  !!! previous run crashed on {} in {} -> will skip this pair".format(hn, fp), flush=True)
    skip_pairs.add((fp, hn))
    with open(os.path.join(tmpdir, SKIP_FILE_NAME), "a") as f:
        f.write(fp + "\t" + hn + "\n")


def hierarchical_merge(indir, subpaths, grooming_str, batch_size, tmpdir, keep_tmp,
                       groups, verbose_mem, mem_limit, resume):
    os.makedirs(tmpdir, exist_ok=True)

    current = [os.path.join(indir, sp, "response.root") for sp in subpaths]
    bad_files, problems = [], []
    skip_pairs = load_skip_pairs(tmpdir)
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

            check_crash_breadcrumb(out_path, skip_pairs, tmpdir)

            if resume and output_is_complete(out_path):
                print("  batch {:>4d}/{:<4d}  already done, skipping".format(ib + 1, n_batches), flush=True)
                next_level.append(out_path)
                continue

            res, n_used = merge_batch(batch, out_path, groups, bad_files, problems,
                                      skip_pairs, level, verbose_mem)
            print("  batch {:>4d}/{:<4d}  {:>4d}/{:<4d} files -> {}".format(
                ib + 1, n_batches, n_used, len(batch),
                os.path.basename(out_path) if res else "EMPTY"), flush=True)
            if res:
                next_level.append(res)

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

    final_path = os.path.join(indir, "response_{}_forunfolding_merged.root".format(grooming_str))
    if os.path.dirname(os.path.abspath(current[0])) == os.path.abspath(tmpdir):
        shutil.move(current[0], final_path)
    else:
        shutil.copy(current[0], final_path)

    if not keep_tmp and os.path.isdir(tmpdir):
        leftovers = [x for x in os.listdir(tmpdir) if x != SKIP_FILE_NAME]
        if not leftovers:
            shutil.rmtree(tmpdir, ignore_errors=True)

    print("\nwritten to {}".format(final_path))
    print_mem("final")
    return final_path, bad_files, problems


# ----------------------------------------------------------------------
# misc
# ----------------------------------------------------------------------
def read_subpaths(filelist):
    out = []
    with open(filelist) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                out.append(line.strip('/'))
    return out


def check_final(final_path, grooming_str):
    f = open_file(final_path)
    if f is None:
        print("could not open final file for checking")
        return
    keys = set(k.GetName() for k in f.GetListOfKeys())
    close_file(f)
    missing = [h.format(grooming_str) for h in hist_list if h.format(grooming_str) not in keys]
    if missing:
        print("\nWARNING: {} histograms absent from final output:".format(len(missing)))
        for m in missing:
            print("   " + m)
    else:
        print("\nall {} histograms present in final output".format(len(hist_list)))


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', required=True)
    parser.add_argument('--filelist', default=default_subpath_file)
    parser.add_argument('--groomingstr', default="groomed", help="choose either 'groomed' or 'ungroomed'")
    parser.add_argument('--batch-size', type=int, default=200,
                        help="files per mini-merge; does NOT affect memory, only I/O")
    parser.add_argument('--group-budget-mb', type=float, default=500.0,
                        help="max MB of accumulators held at once; peak RAM is ~2x this")
    parser.add_argument('--tmpdir', default=None)
    parser.add_argument('--keep-tmp', action='store_true')
    parser.add_argument('--resume', action='store_true',
                        help="skip intermediate files that were completed by a previous run")
    parser.add_argument('--verbose-mem', action='store_true')
    parser.add_argument('--mem-limit-mb', type=float, default=0.0)
    parser.add_argument('--scan-only', action='store_true')
    parser.add_argument('--roounfold-include', default=None,
                        help="directory containing RooUnfoldResponse.h (auto-detected if omitted)")
    flags = parser.parse_args()

    load_roounfold(flags.roounfold_include)
    subpaths = read_subpaths(flags.filelist)
    print_mem("after libRooUnfold load")

    sizes = scan_sizes(os.path.join(flags.dir, subpaths[0], "response.root"), flags.groomingstr)

    ru_names = [h.format(flags.groomingstr) for h in hist_list if "roounfold_response" in h]
    zero = [n for n in ru_names if sizes.get(n, 0.0) == 0.0]
    if zero:
        print("WARNING: these responses were missing or unreadable in the sample file:")
        for n in zero:
            print("   " + n)

    if flags.scan_only:
        raise SystemExit(0)
    groups = build_groups(flags.groomingstr, sizes, flags.group_budget_mb)

    tmpdir = flags.tmpdir if flags.tmpdir else os.path.join(flags.dir, "tmp_merge")
    final_path, bad_files, problems = hierarchical_merge(
        flags.dir, subpaths, flags.groomingstr, flags.batch_size, tmpdir,
        flags.keep_tmp, groups, flags.verbose_mem, flags.mem_limit_mb, flags.resume)

    print("\n{} files missing / unreadable".format(len(bad_files)))
    if bad_files:
        log = os.path.join(flags.dir, "missing_files.txt")
        with open(log, "w") as f:
            for b in bad_files:
                f.write(b + "\n")
        print("list written to {}".format(log))

    if problems:
        counts = {}
        for _, _, hn, _ in problems:
            counts[hn] = counts.get(hn, 0) + 1
        print("\n{} skipped (file, histogram) pairs:".format(len(problems)))
        for hn, n in sorted(counts.items(), key=lambda x: -x[1]):
            print("   {:<45s} {:6d}".format(hn, n))
        log = os.path.join(flags.dir, "merge_problems_{}.txt".format(flags.groomingstr))
        with open(log, "w") as f:
            for lvl, fp, hn, msg in problems:
                f.write("{}\t{}\t{}\t{}\n".format(lvl, fp, hn, msg))
        print("details written to {}".format(log))

    check_final(final_path, flags.groomingstr)