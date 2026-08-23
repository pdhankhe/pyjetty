import ROOT
import argparse
import os
import glob
import tempfile
import shutil

'''
HOW TO RUN
python3 hadd_manual.py --dir /global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/57259276  --batch 50
'''

# Load RooUnfold dynamic libraries safely before reading files
def load_roounfold():
    res1 = ROOT.gSystem.Load("libRooUnfold")
    res2 = ROOT.gSystem.Load("libRooUnfold.so")
    # gSystem.Load returns 0 on success, 1 if already loaded, -1 on error
    if res1 < 0 and res2 < 0:
        print("Warning: Could not load libRooUnfold via gSystem.Load.")

def detach_object(obj):
    """Detaches ROOT objects from file directory memory to avoid garbage collection crashes."""
    if hasattr(obj, "SetDirectory"):
        obj.SetDirectory(0)
    ROOT.SetOwnership(obj, True) # Hand ownership to Python GC

def merge_batch(files, batch_id, temp_dir):
    """Merges a small batch of files into a temporary root file to save memory."""
    if not files:
        return None

    first_file = ROOT.TFile.Open(files[0], "READ")
    if not first_file or first_file.IsZombie():
        print(f"Error opening file: {files[0]}")
        return None

    merged_objects = {}
    keys = first_file.GetListOfKeys()
    
    for key in keys:
        obj = key.ReadObj()
        if obj and hasattr(obj, "Add"):
            name = obj.GetName()
            cloned = obj.Clone()
            detach_object(cloned)
            merged_objects[name] = cloned

    first_file.Close()

    # Merge remaining files in this batch
    for f_path in files[1:]:
        f = ROOT.TFile.Open(f_path, "READ")
        if not f or f.IsZombie():
            print(f"Warning: Could not open {f_path}, skipping.")
            continue

        for name, merged_obj in list(merged_objects.items()):
            obj = f.Get(name)
            
            # Check for CPyCppyy_NoneType or missing objects
            if obj is None:
                print(f"Warning: Object '{name}' not found in {f_path}")
                continue
                
            if hasattr(obj, "Add"):
                merged_obj.Add(obj)
            else:
                print(f"Warning: Object '{name}' in {f_path} has no 'Add' method.")
                
        f.Close()

    # Write batch results to a temporary file on disk
    batch_file_path = os.path.join(temp_dir, f"batch_{batch_id}.root")
    out_file = ROOT.TFile(batch_file_path, "RECREATE")
    out_file.cd()
    
    for obj in merged_objects.values():
        obj.Write()
        
    out_file.Close()

    # Clear objects from memory explicitly
    merged_objects.clear()
    return batch_file_path

def merge_all_hists(root_dir, output_filename="response_merged.root", batch_size=50):
    load_roounfold()

    search_pattern = os.path.join(root_dir, "**", "response.root")
    files = glob.glob(search_pattern, recursive=True)

    if not files:
        print(f"No files found matching pattern: {search_pattern}")
        return

    print(f"Found {len(files)} files to merge. Processing in batches of {batch_size}...")

    temp_dir = tempfile.mkdtemp(prefix="hadd_batches_")
    batch_files = []

    try:
        # Stage 1: Merge into batches
        for i in range(0, len(files), batch_size):
            batch = files[i : i + batch_size]
            batch_id = i // batch_size
            print(f"Processing batch {batch_id} ({i // batch_size + 1}/{(len(files)-1)//batch_size + 1})...")
            batch_path = merge_batch(batch, batch_id, temp_dir)
            if batch_path:
                batch_files.append(batch_path)

        if not batch_files:
            print("No data merged into batches.")
            return

        # Stage 2: Merge intermediate batch files into final output
        print(f"Merging {len(batch_files)} batch files into final output...")

        first_batch = ROOT.TFile.Open(batch_files[0], "READ")
        if not first_batch or first_batch.IsZombie():
            print("Could not open first batch file.")
            return

        final_objects = {}
        for key in first_batch.GetListOfKeys():
            obj = key.ReadObj()
            if obj and hasattr(obj, "Add"):
                cloned = obj.Clone()
                detach_object(cloned)
                final_objects[obj.GetName()] = cloned

        first_batch.Close()

        for b_path in batch_files[1:]:
            f = ROOT.TFile.Open(b_path, "READ")
            if not f or f.IsZombie():
                continue
            for name, final_obj in final_objects.items():
                obj = f.Get(name)
                if obj is not None and hasattr(obj, "Add"):
                    final_obj.Add(obj)
            f.Close()

        # Final Write
        full_output_path = os.path.join(root_dir, output_filename)
        out_file = ROOT.TFile(full_output_path, "RECREATE")
        out_file.cd()
        for obj in final_objects.values():
            obj.Write()
        out_file.Close()

        print(f"Successfully merged all histograms into {full_output_path}")

    finally:
        shutil.rmtree(temp_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', required=True, help="Root directory to search for response.root files")
    parser.add_argument('--out', default="response_merged.root", help="Output filename")
    parser.add_argument('--batch', type=int, default=50, help="Number of files per batch (default: 50)")
    args = parser.parse_args()

    merge_all_hists(args.dir, args.out, args.batch)