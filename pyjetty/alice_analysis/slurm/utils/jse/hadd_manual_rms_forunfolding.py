import ROOT
import argparse
import os
import gc
import json
import math
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

NEVENTS_NAME = "hNevents"

hist_list = [NEVENTS_NAME,
             'roounfold_response_1D_{}', 'roounfold_response_1D_{}_oldversion', 'jet_match_rec_pur_num_{}', 'jet_all_rec_pur_den_{}',
             'jet_match_gen_eff_num_{}', 'jet_all_gen_eff_den_{}',
             '{}_full_ungroomed_roounfold_response', '{}_AA_roounfold_response',
             '{}_AB_roounfold_response', '{}_BB_roounfold_response', '{}_rad_roounfold_response',
             '{}_full_ungroomed_reco_matched', '{}_full_ungroomed_reco_all', '{}_AA_reco_matched', '{}_AA_reco_all',
             '{}_AB_reco_matched', '{}_AB_reco_all', '{}_BB_reco_matched', '{}_BB_reco_all',
             '{}_rad_reco_matched', '{}_rad_reco_all',
             '{}_full_ungroomed_gen_matched', '{}_full_ungroomed_gen_all', '{}_AA_gen_matched', '{}_AA_gen_all',
             '{}_AB_gen_matched', '{}_AB_gen_all', '{}_BB_gen_matched', '{}_BB_gen_all',
             '{}_rad_gen_matched', '{}_rad_gen_all']

default_subpath_file = "/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/anchmc_subpath_filelist.txt"

COMPLETE_MARKER = "__merge_complete__"
LEDGER_NAME = "__event_ledger__"
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


def check_mem(mem_limit, where):
    if not mem_limit:
        return
    cur, _ = rss_mb()
    if cur > mem_limit:
        raise MemoryError("RSS {:.0f} MB exceeds --mem-limit-mb {:.0f} ({})".format(cur, mem_limit, where))


# ----------------------------------------------------------------------
# RooUnfold loading (library + headers so PyROOT knows the methods)
# ----------------------------------------------------------------------
_RU_CLASS = None

_CPP_HELPERS = r'''
#include "RooUnfoldResponse.h"
#include "TH1.h"
namespace mergehelp {
  // Delete the internal histograms of a RooUnfoldResponse.
  // The response object itself must NOT be deleted afterwards (its pointers now dangle);
  // it is a few hundred bytes and is intentionally leaked. This is safe whether or not
  // RooUnfold's destructor owns the histograms.
  inline void free_response_hists(RooUnfoldResponse* r) {
    if (!r) return;
    const TH1* hs[4] = { r->Hmeasured(), r->Hfakes(), r->Htruth(), r->Hresponse() };
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < i; ++j)
        if (hs[i] && hs[j] == hs[i]) hs[i] = nullptr;
    for (int i = 0; i < 4; ++i)
      if (hs[i]) delete const_cast<TH1*>(hs[i]);
  }
}
'''


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
                  os.path.join(prefix, "src"),
                  prefix]
    for d in cands:
        if os.path.isfile(os.path.join(d, header)):
            return d, cands
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

    if not ROOT.gInterpreter.Declare(_CPP_HELPERS):
        raise RuntimeError("failed to compile C++ cleanup helper")

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
    """Return a correctly typed proxy for an object read from a file.
    ALWAYS use the return value.
    - TH1: Python owns it and deletes it when the last reference goes away.
    - RooUnfoldResponse: Python never deletes it; call free_obj() to release its memory."""
    if isinstance(obj, ROOT.TH1):
        obj.SetDirectory(0)
        ROOT.SetOwnership(obj, True)
        return obj
    if is_response(obj):
        ROOT.SetOwnership(obj, False)
        if hasattr(obj, "Hresponse"):
            return obj
        cast = ROOT.BindObject(ROOT.addressof(obj), get_ru_class())
        ROOT.SetOwnership(cast, False)
        return cast
    ROOT.SetOwnership(obj, True)
    return obj


def free_obj(obj):
    """Release the memory of an object obtained through take_ownership()."""
    if obj is None:
        return
    if is_response(obj):
        try:
            ROOT.mergehelp.free_response_hists(obj)
        except Exception as e:
            print("  !!! could not free response internals: {}".format(short(e)), flush=True)
    # TH1s are Python-owned and are deleted when the caller drops its reference


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
# event accounting
#
# ledger = {
#   "t0":      total level-0 events (sum of hNevents) of all readable input files in this subtree
#   "events":  {hist_name: events whose contents are included in that histogram}
#   "losses":  [[level, file, hist_name, n_events, reason], ...]   (cumulative over the subtree)
#   "unknown": [[level, file, reason], ...]   files whose event count could not be determined
# }
# invariant for every histogram h:  events[h] + sum(losses for h) == t0
# ----------------------------------------------------------------------
def near(a, b):
    return abs(a - b) <= 1e-6 * max(1.0, abs(b))


def read_nevents(fin):
    """Number of events stored in hNevents at FindBin(0). Returns (value, None) or (None, reason)."""
    o = fin.Get(NEVENTS_NAME)
    if not o:
        return None, "{} missing".format(NEVENTS_NAME)
    o = take_ownership(o)
    if not isinstance(o, ROOT.TH1):
        cname = o.ClassName()
        free_obj(o)
        return None, "{} is a {}, not a histogram".format(NEVENTS_NAME, cname)
    v = float(o.GetBinContent(o.GetXaxis().FindFixBin(0.0)))
    del o
    if math.isnan(v) or math.isinf(v) or v < 0:
        return None, "{} has invalid content {}".format(NEVENTS_NAME, v)
    return v, None


def read_ledger(fin):
    o = fin.Get(LEDGER_NAME)
    if not o:
        return None
    ROOT.SetOwnership(o, True)
    try:
        return json.loads(str(o.GetTitle()))
    except ValueError:
        return None


def write_ledger(ledger):
    n = ROOT.TNamed(LEDGER_NAME, json.dumps(ledger, separators=(",", ":")))
    n.Write()


def empty_ledger(names):
    return {"t0": 0.0, "events": {h: 0.0 for h in names}, "losses": [], "unknown": []}


def losses_by_hist(ledger):
    out = {}
    for _, _, h, n, _ in ledger["losses"]:
        out[h] = out.get(h, 0.0) + n
    return out


def closure_errors(ledger, names):
    lost = losses_by_hist(ledger)
    errs = []
    for h in names:
        inc = ledger["events"].get(h, 0.0)
        if not near(inc + lost.get(h, 0.0), ledger["t0"]):
            errs.append("closure failed for {}: included {:.0f} + lost {:.0f} != total {:.0f}".format(
                h, inc, lost.get(h, 0.0), ledger["t0"]))
    return errs


def verify_output(out_path, ledger, names, inputs=None):
    """Read hNevents back from a written file and compare with the expected sum of its inputs."""
    errs = []
    f = open_file(out_path)
    if f is None:
        return ["cannot reopen {} for verification".format(out_path)]
    try:
        r, why = read_nevents(f)
        stored = read_ledger(f)
    finally:
        close_file(f)

    exp = ledger["events"].get(NEVENTS_NAME, 0.0)
    if r is None:
        if exp:
            errs.append("{} not readable in output ({}) but {:.0f} events expected".format(NEVENTS_NAME, why, exp))
    elif not near(r, exp):
        errs.append("{} in output = {:.0f}, expected sum of merged inputs = {:.0f} (diff {:+.0f})".format(
            NEVENTS_NAME, r, exp, r - exp))
        if inputs:
            for fpath, led in inputs:
                errs.append("    input {:>14.0f} events  {}".format(led["events"].get(NEVENTS_NAME, 0.0), fpath))
    if stored is None:
        errs.append("event ledger could not be read back from output")
    errs += closure_errors(ledger, names)
    return errs


def print_batch_losses(batch_losses):
    if not batch_losses:
        return
    per_h = {}
    for _, fp, h, n, _ in batch_losses:
        tot, files = per_h.get(h, (0.0, set()))
        files.add(fp)
        per_h[h] = (tot + n, files)
    print("    events NOT included in this batch:", flush=True)
    for h, (tot, files) in sorted(per_h.items()):
        print("      {:<45s} {:>14.0f} events from {} file(s)".format(h, tot, len(files)), flush=True)


def level_check(level, ledgers, orphan):
    t0 = orphan["t0"] + sum(l["t0"] for l in ledgers)
    nev = sum(l["events"].get(NEVENTS_NAME, 0.0) for l in ledgers)
    lost = sum(n for l in ledgers + [orphan] for (_, _, h, n, _) in l["losses"] if h == NEVENTS_NAME)
    ok = near(nev + lost, t0)
    print("\n  level {} event check: {:.0f} events in readable inputs, {:.0f} in {} of {} outputs, "
          "{:.0f} lost  ->  {}".format(level, t0, nev, NEVENTS_NAME, len(ledgers), lost,
                                       "OK" if ok else "MISMATCH ({:+.0f})".format(nev + lost - t0)), flush=True)
    return ok


# ----------------------------------------------------------------------
# size scan and grouping
# ----------------------------------------------------------------------
def scan_sizes(sample_file, grooming_str, quiet=False):
    sizes = {}
    print_mem("before size scan")
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
            print("  {:<40s} {:>28s}  {:8.2f} MB".format(hname, o.ClassName(), s))
        free_obj(o)
        del o
    n, why = read_nevents(f)
    close_file(f)
    sweep()
    if not quiet:
        print("  {:<40s} {:>28s}  {:8.1f} MB  <-- one full set".format("TOTAL", "", total))
        print("  events in this file ({} at FindBin(0)): {}".format(
            NEVENTS_NAME, "{:.0f}".format(n) if n is not None else "UNREADABLE: " + why))
        print("---------------------------------------\n")
    print_mem("after size scan (should be close to 'before')")
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
def load_complete_ledger(path):
    """Ledger of a finished intermediate file, or None if the file is absent/incomplete/old-format."""
    f = open_file(path)
    if f is None:
        return None
    try:
        if not f.GetListOfKeys().FindObject(COMPLETE_MARKER):
            return None
        return read_ledger(f)
    finally:
        close_file(f)


def merge_batch(file_paths, out_path, groups, names, bad_files, problems, skip_pairs, level,
                known_ledgers, mem_limit=0.0, verbose_mem=False):
    """One group of objects at a time; each (file, histogram) is merged independently.
    Returns (out_path or None, n_inputs_used, ledger, losses_in_this_batch, verification_errors)."""
    inputs, inherited_losses, unknown, batch_losses = [], [], [], []
    t0 = 0.0

    def lose(fpath, h, n, reason):
        if n:
            batch_losses.append([level, fpath, h, n, reason])

    # ---- pre-pass: which inputs are usable, and how many events does each carry ----
    for fpath in file_paths:
        f = open_file(fpath)
        if f is None:
            bad_files.append(fpath)
            if level == 0:
                unknown.append([level, fpath, "cannot open file"])
            else:
                kl = known_ledgers.get(fpath)
                if kl is None:
                    unknown.append([level, fpath, "cannot open intermediate file; its ledger is unknown"])
                else:
                    t0 += kl["t0"]
                    inherited_losses.extend(kl["losses"])
                    unknown.extend(kl["unknown"])
                    for h in names:
                        lose(fpath, h, kl["events"].get(h, 0.0), "intermediate file unreadable")
            continue
        try:
            if level == 0:
                n, why = read_nevents(f)
                if n is None:
                    msg = why + " -> whole file skipped (its histograms could not be normalised)"
                    unknown.append([level, fpath, msg])
                    problems.append((level, fpath, "*", msg))
                    print("  !!! {}: {}".format(fpath, msg), flush=True)
                    continue
                led = {"t0": n, "events": {h: n for h in names}, "losses": [], "unknown": []}
            else:
                led = read_ledger(f) or known_ledgers.get(fpath)
                if led is None:
                    msg = "intermediate file has no event ledger -> skipped"
                    unknown.append([level, fpath, msg])
                    problems.append((level, fpath, "*", msg))
                    print("  !!! {}: {}".format(fpath, msg), flush=True)
                    continue
        finally:
            close_file(f)
        inputs.append((fpath, led))
        t0 += led["t0"]
        inherited_losses.extend(led["losses"])
        unknown.extend(led["unknown"])

    ledger = {"t0": t0, "events": {h: 0.0 for h in names}, "losses": [], "unknown": unknown}

    def finish_ledger():
        ledger["losses"] = inherited_losses + batch_losses

    if not inputs:
        finish_ledger()
        return None, 0, ledger, batch_losses, []

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

    acc, sigs = {}, {}
    try:
        for ig, group in enumerate(groups):
            acc, sigs = {}, {}

            for fpath, led in inputs:
                if verbose_mem:
                    print_mem("    pre  {}".format(fpath[-45:]))

                fin = open_file(fpath)
                if fin is None:
                    problems.append((level, fpath, "*", "could not reopen file"))
                    for h_name in group:
                        lose(fpath, h_name, led["events"].get(h_name, 0.0), "could not reopen file")
                    continue

                try:
                    for h_name in group:
                        n_ev = led["events"].get(h_name, 0.0)
                        if (fpath, h_name) in skip_pairs:
                            lose(fpath, h_name, n_ev, "skipped: crashed here in a previous run")
                            continue
                        crumb(fpath + "\t" + h_name)
                        obj = None
                        kept = False
                        try:
                            obj = fin.Get(h_name)
                            if not obj:
                                obj = None
                                if level == 0 or n_ev:
                                    problems.append((level, fpath, h_name, "missing"))
                                lose(fpath, h_name, n_ev, "missing")
                                continue
                            obj = take_ownership(obj)
                            sig = obj_signature(obj)

                            if h_name not in acc:
                                # first good copy becomes the accumulator (no clone needed)
                                acc[h_name] = obj
                                sigs[h_name] = sig
                                kept = True
                            else:
                                if sig != sigs[h_name]:
                                    raise RuntimeError("binning mismatch: {} vs accumulator {}".format(
                                        sig, sigs[h_name]))
                                add_into(acc[h_name], obj)
                            ledger["events"][h_name] += n_ev
                        except Exception as e:
                            msg = short(e)
                            problems.append((level, fpath, h_name, msg))
                            lose(fpath, h_name, n_ev, msg)
                            print("  !!! skipping {} from {}: {}".format(h_name, fpath, msg), flush=True)
                        finally:
                            if obj is not None and not kept:
                                free_obj(obj)
                            obj = None
                finally:
                    close_file(fin)

                if verbose_mem:
                    print_mem("    post")
                check_mem(mem_limit, "level {} after {}".format(level, fpath))

            fout.cd()
            for h_name in group:
                if h_name not in acc:
                    continue
                o = acc.pop(h_name)
                try:
                    if o.Write(h_name) <= 0:
                        raise RuntimeError("Write returned 0 bytes")
                    n_written += 1
                except Exception as e:
                    msg = "write failed: " + short(e)
                    problems.append((level, out_path, h_name, msg))
                    lose(out_path, h_name, ledger["events"][h_name], msg)
                    ledger["events"][h_name] = 0.0
                    print("  !!! could not write {}: {}".format(h_name, short(e)), flush=True)
                finally:
                    free_obj(o)
                del o
            acc.clear()
            sweep()
            if verbose_mem:
                print_mem("    group {}/{} done".format(ig + 1, len(groups)))

        finish_ledger()
        if n_written:
            fout.cd()
            write_ledger(ledger)
            ROOT.TNamed(COMPLETE_MARKER, str(len(inputs))).Write()
    finally:
        # free any accumulators left over if we bailed out mid-group
        for o in acc.values():
            free_obj(o)
        acc.clear()
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
        return None, 0, ledger, batch_losses, []

    errs = verify_output(out_path, ledger, names, inputs)
    return out_path, len(inputs), ledger, batch_losses, errs


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


def report_check_errors(out_path, errs, check_failures):
    if not errs:
        return
    print("  !!! EVENT CHECK FAILED for {}:".format(os.path.basename(out_path)), flush=True)
    for e in errs:
        print("      " + e, flush=True)
    check_failures.append((out_path, errs))


def hierarchical_merge(indir, subpaths, zcut_str, grooming_str, names, batch_size, tmpdir, keep_tmp,
                       groups, verbose_mem, mem_limit, resume):
    os.makedirs(tmpdir, exist_ok=True)

    current = [os.path.join(indir, sp, "response{}.root".format(zcut_str)) for sp in subpaths]
    bad_files, problems, check_failures = [], [], []
    skip_pairs = load_skip_pairs(tmpdir)
    known_ledgers = {}
    orphan = empty_ledger(names)   # batches that produced no output file (events all lost)
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

            if resume:
                led = load_complete_ledger(out_path)
                if led is not None:
                    known_ledgers[out_path] = led
                    print("  batch {:>4d}/{:<4d}  already done, skipping".format(ib + 1, n_batches), flush=True)
                    report_check_errors(out_path, verify_output(out_path, led, names), check_failures)
                    next_level.append(out_path)
                    continue

            res, n_used, led, batch_losses, errs = merge_batch(
                batch, out_path, groups, names, bad_files, problems, skip_pairs, level,
                known_ledgers, mem_limit, verbose_mem)
            status = "OK" if (res and not errs) else ("CHECK FAILED" if res else "")
            print("  batch {:>4d}/{:<4d}  {:>4d}/{:<4d} files -> {}   {:.0f} events  {}".format(
                ib + 1, n_batches, n_used, len(batch),
                os.path.basename(out_path) if res else "EMPTY",
                led["events"].get(NEVENTS_NAME, 0.0), status), flush=True)
            print_batch_losses(batch_losses)
            report_check_errors(out_path, errs, check_failures)

            if res:
                known_ledgers[res] = led
                next_level.append(res)
            else:
                orphan["t0"] += led["t0"]
                orphan["losses"] += led["losses"]
                orphan["unknown"] += led["unknown"]

            print_mem("after level {} batch {}".format(level, ib + 1))
            check_mem(mem_limit, "after level {} batch {}".format(level, ib + 1))

        if not next_level:
            raise RuntimeError("level {} produced no output".format(level))

        if not level_check(level, [known_ledgers[p] for p in next_level], orphan):
            check_failures.append(("level {}".format(level), ["level total mismatch"]))

        if not keep_tmp:
            for f in previous_level_files:
                if f not in next_level:
                    try:
                        os.remove(f)
                    except OSError:
                        pass
                known_ledgers.pop(f, None)

        previous_level_files = next_level
        current = next_level
        level += 1

    final_path = os.path.join(indir, "response_{}_forunfolding_merged.root".format(grooming_str))
    if os.path.dirname(os.path.abspath(current[0])) == os.path.abspath(tmpdir):
        shutil.move(current[0], final_path)
    else:
        shutil.copy(current[0], final_path)

    final_ledger = None
    f = open_file(final_path)
    if f is not None:
        final_ledger = read_ledger(f)
        close_file(f)
    if final_ledger is not None:
        final_ledger["t0"] += orphan["t0"]
        final_ledger["losses"] += orphan["losses"]
        final_ledger["unknown"] += orphan["unknown"]

    if not keep_tmp and os.path.isdir(tmpdir):
        leftovers = [x for x in os.listdir(tmpdir) if x != SKIP_FILE_NAME]
        if not leftovers:
            shutil.rmtree(tmpdir, ignore_errors=True)

    print("\nwritten to {}".format(final_path))
    print_mem("final")
    return final_path, bad_files, problems, final_ledger, check_failures


# ----------------------------------------------------------------------
# reporting
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


def event_report(final_path, ledger, names, check_failures, outdir, grooming_str):
    print("\n================ EVENT ACCOUNTING ================")
    if ledger is None:
        print("no event ledger in the final file (was there only one input file?)")
        return

    f = open_file(final_path)
    r, why = read_nevents(f) if f is not None else (None, "cannot open final file")
    close_file(f)

    t0 = ledger["t0"]
    exp = ledger["events"].get(NEVENTS_NAME, 0.0)
    print("events in readable input files (sum of {} over level-0 files): {:.0f}".format(NEVENTS_NAME, t0))
    print("input files skipped because their event count is unknown:      {}".format(len(ledger["unknown"])))
    if r is None:
        print("final {}: UNREADABLE ({})".format(NEVENTS_NAME, why))
    else:
        print("final {}: {:.0f}   (expected {:.0f})  {}".format(
            NEVENTS_NAME, r, exp, "OK" if near(r, exp) else "MISMATCH ({:+.0f})".format(r - exp)))

    lost = losses_by_hist(ledger)
    print("\n  {:<45s} {:>14s} {:>14s} {:>9s}".format("histogram", "events incl.", "missing", "missing%"))
    for h in names:
        inc = ledger["events"].get(h, 0.0)
        miss = t0 - inc
        flag = "" if near(miss, lost.get(h, 0.0)) else "   <-- inconsistent with loss log"
        print("  {:<45s} {:>14.0f} {:>14.0f} {:>8.4f}%{}".format(
            h, inc, miss, 100.0 * miss / t0 if t0 else 0.0, flag))

    errs = closure_errors(ledger, names)
    if check_failures or errs:
        print("\n!!! {} event check failure(s) during merging, {} closure error(s) in final ledger".format(
            len(check_failures), len(errs)))
        for e in errs:
            print("    " + e)
    else:
        print("\nall per-batch and per-level event checks passed")

    log = os.path.join(outdir, "event_losses_{}.txt".format(grooming_str))
    with open(log, "w") as fl:
        fl.write("# events NOT included, per (file, histogram)\n")
        fl.write("# level\tfile\thistogram\tn_events\treason\n")
        for lvl, fp, h, n, reason in sorted(ledger["losses"], key=lambda x: (x[1], x[2])):
            fl.write("{}\t{}\t{}\t{:.0f}\t{}\n".format(lvl, fp, h, n, reason))
        fl.write("\n# files whose event count could not be determined (not included anywhere)\n")
        fl.write("# level\tfile\treason\n")
        for lvl, fp, reason in ledger["unknown"]:
            fl.write("{}\t{}\t{}\n".format(lvl, fp, reason))
        if check_failures:
            fl.write("\n# event check failures\n")
            for where, es in check_failures:
                for e in es:
                    fl.write("{}\t{}\n".format(where, e))
    n_files = len(set(x[1] for x in ledger["losses"]))
    print("{} loss entries from {} file(s) written to {}".format(len(ledger["losses"]), n_files, log))


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
    parser.add_argument('--zcut', type=float, default=0.1, help="which zcut do you want to look at (for naming)")
    parser.add_argument('--keep-tmp', action='store_true')
    parser.add_argument('--resume', action='store_true',
                        help="skip intermediate files that were completed by a previous run")
    parser.add_argument('--verbose-mem', action='store_true')
    parser.add_argument('--mem-limit-mb', type=float, default=0.0)
    parser.add_argument('--scan-only', action='store_true')
    parser.add_argument('--roounfold-include', default=None,
                        help="directory containing RooUnfoldResponse.h (auto-detected if omitted)")
    flags = parser.parse_args()

    names = [h.format(flags.groomingstr) for h in hist_list]

    load_roounfold(flags.roounfold_include)
    subpaths = read_subpaths(flags.filelist)
    print_mem("after libRooUnfold load")

    zcut_str = ""
    if flags.zcut == 0.1:
        zcut_str = f"_zcut0p1"
    elif flags.zcut == 0.2:
        zcut_str = f"_zcut0p2"

    sizes = scan_sizes(os.path.join(flags.dir, subpaths[0], "response{}.root".format(zcut_str)), flags.groomingstr)

    ru_names = [n for n in names if "roounfold_response" in n]
    zero = [n for n in ru_names if sizes.get(n, 0.0) == 0.0]
    if zero:
        print("WARNING: these responses were missing or unreadable in the sample file:")
        for n in zero:
            print("   " + n)

    if flags.scan_only:
        raise SystemExit(0)
    groups = build_groups(flags.groomingstr, sizes, flags.group_budget_mb)

    tmpdir = flags.tmpdir if flags.tmpdir else os.path.join(flags.dir, "tmp_merge")
    final_path, bad_files, problems, final_ledger, check_failures = hierarchical_merge(
        flags.dir, subpaths, zcut_str, flags.groomingstr, names, flags.batch_size, tmpdir,
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
    event_report(final_path, final_ledger, names, check_failures, flags.dir, flags.groomingstr)