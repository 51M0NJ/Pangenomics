from pathlib import Path
import pickle

def chunk(lst, size):
    for i in range(0, len(lst), size):
        yield lst[i:i + size]

files = list(Path('/nfs3_ib/nfs-ip34/fast2/def-jacquesp/Projet_Pangenome/1_sample_annotations_JSON/').glob('*.bakta.json'))
chunks_path = '/nfs3_ib/nfs-ip34/fast2/def-jacquesp/Projet_Pangenome/2_chunked_annotations_json_PKL'
chunk_size = 1000
file_chunks = list(chunk(files, chunk_size))
# Save chunk lists for external sbatch jobs
Path(chunks_path).mkdir(exist_ok=True)
for i, c in enumerate(file_chunks):
    with open(f'{chunks_path}/chunk_{i:04d}.pkl', 'wb+') as f:
        pickle.dump(c, f)
print(f'Saved {len(file_chunks)} chunk files.')
