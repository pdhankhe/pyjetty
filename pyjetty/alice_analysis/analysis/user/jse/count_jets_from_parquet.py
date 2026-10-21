from pathlib import Path
import pandas as pd

look_at_data=False
look_at_anchmc=True

# IF DATA:
if look_at_data:
    # Target directory search path
    base_dir = Path("/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data/58416144/")
    target_filename = "DataJetsForAnalysis.parquet"

    # Find all matching parquet files recursively
    file_paths = list(base_dir.rglob(target_filename))

    print(f"Found {len(file_paths)} parquet files. Counting unique jets...")

    total_jets = 0
    total_constituents = 0

    for file_path in file_paths:
        # Read only the 'event_id' and 'jet_id' columns to save memory & time
        df = pd.read_parquet(file_path, columns=["event_id", "jet_id"])
        
        # Count unique jets in this file using (event_id, jet_id)
        n_jets = df[["event_id", "jet_id"]].drop_duplicates().shape[0]
        
        total_jets += n_jets
        total_constituents += len(df)
        
        print(f"File: {file_path.relative_to(base_dir)} | Jets: {n_jets:,} | Constituents: {len(df):,}")

    print("-" * 60)
    print(f"TOTAL UNIQUE JETS:         {total_jets:,}")
    print(f"TOTAL CONSTITUENTS (ROWS): {total_constituents:,}")


# IF ANCH MC:
if look_at_anchmc:
    # Target directory search path
    base_dir = Path("/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/58919321")
    file_paths = list(base_dir.rglob("jets_out.parquet"))
    print("len filepaths", len(file_paths))

    total_stats = {
        "det": {"all": 0, "matched": 0},
        "part": {"all": 0, "matched": 0}
    }

    for file_path in file_paths:
        # Read only the level and matching columns (1 row = 1 jet)
        df = pd.read_parquet(file_path, columns=["event", "level", "jet_index", "is_matched"])
        
        for lvl in ["det", "part"]:
            subset = df[df["level"] == lvl]
        
            # Deduplicate using (event, jet_index) for this specific level
            unique_all = subset[["event", "jet_index"]].drop_duplicates().shape[0]
            
            # Unique matched jets for this level
            matched_subset = subset[subset["is_matched"] == True]
            unique_matched = matched_subset[["event", "jet_index"]].drop_duplicates().shape[0]
            
            total_stats[lvl]["all"] += unique_all
            total_stats[lvl]["matched"] += unique_matched

        print(f"Processed: {file_path.relative_to(base_dir)}")

    print("=" * 60)
    print("DETECTOR LEVEL ('det'):")
    print(f"  Total Unique Jets:   {total_stats['det']['all']:,}")
    print(f"  Matched Unique Jets: {total_stats['det']['matched']:,}")
    if total_stats['det']['all'] > 0:
        print(f"  Match Efficiency:    {total_stats['det']['matched'] / total_stats['det']['all']:.2%}")

    print("-" * 60)
    print("PARTICLE LEVEL ('part'):")
    print(f"  Total Unique Jets:   {total_stats['part']['all']:,}")
    print(f"  Matched Unique Jets: {total_stats['part']['matched']:,}")
    if total_stats['part']['all'] > 0:
        print(f"  Match Efficiency:    {total_stats['part']['matched'] / total_stats['part']['all']:.2%}")
    print("=" * 60)
