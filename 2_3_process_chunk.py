import pickle
import ijson
from pathlib import Path
import sys
import pandas as pd
from multiprocessing import Pool

def split_list(lst, n):
    k, m = divmod(len(lst), n)
    return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n)]

def process_files_batch(file_batch):
    results = []
    for jf in file_batch:
        shortname = jf.name.rsplit('.', 2)[0]
        seen = set()
        with open(jf, 'rb') as f:
            for gene in ijson.items(f, 'features.item.gene'):
                if gene and gene not in seen:
                    results.append((shortname, gene))
                    seen.add(gene)
    return results

# Get chunk ID from argument
chunk_id = int(sys.argv[1])
outpath = '/nfs3_ib/nfs-ip34/fast2/def-jacquesp/Projet_Pangenome/'
chunk_path = f'{outpath}/2_chunked_annotations_json_PKL/chunk_{chunk_id:04d}.pkl'
with open(chunk_path, 'rb') as f:
    file_batch = pickle.load(f)

# Run processing
pools = 24
with Pool(pools) as pool:
    results = list(pool.imap_unordered(process_files_batch, split_list(file_batch,pools)))
flat = [item for sublist in results for item in sublist]
df = pd.DataFrame(flat, columns=['genome', 'gene'])
matrix = pd.crosstab(df['genome'], df['gene']).astype('uint8')

# Save output
output_matrix_path = f'{outpath}/3_partial_cooccurrence_matrices_PKL/result_{chunk_id:04d}.pkl'
Path(f'{outpath}/3_partial_cooccurrence_matrices_PKL').mkdir(exist_ok=True)
with open(output_matrix_path, 'wb') as f:
    pickle.dump(matrix, f)
print(f'Done chunk {chunk_id}, {len(matrix)} entries.')
