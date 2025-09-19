import pandas as pd
import glob

pkl_files = sorted(glob.glob("/fast2/def-jacquesp/Projet_Pangenome/3_partial_cooccurrence_matrices_PKL/result_*.pkl"))
merged_df = pd.DataFrame()
for i, file in enumerate(pkl_files):
    df = pd.read_pickle(file)
    merged_df = pd.concat([merged_df, df])
    
    if (i + 1) % 50 == 0:  # Optional: save intermediate results every 50 files
        print(f"{i+1} files merged...")

print(f"Finished, writing pickle")
merged_df.to_pickle("/fast2/def-jacquesp/Projet_Pangenome/3_partial_cooccurrence_matrices_PKL/cooccurrence_matrix.pkl")
#merged_df.to_csv("merged_results.csv")  # For inspection or cross-tool use TOO LONG
