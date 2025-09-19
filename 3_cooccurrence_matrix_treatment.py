# Master code to run in background to filter the matrix and produce some figures - takes ~3 days. 
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import warnings, random ,hdbscan
import sys, pickle
from datetime import datetime
from collections import Counter as cnt
from sklearn.metrics import pairwise_distances
import umap
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
warnings.filterwarnings("ignore", category=FutureWarning, module="seaborn")
warnings.filterwarnings("ignore", category=FutureWarning, module="sklearn")
figpath='/home/def-jacquesp/jeanneau/9-Pangenomics/'
with open('treatment.log','w') as logfile:logfile.write(f'[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}] loading matrix from pickle...\n')
df = pd.read_pickle('/fast2/def-jacquesp/Projet_Pangenome/cooccurrence_matrix.pkl')
with open('treatment.log','a') as logfile: logfile.write(f'[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}] matrix loaded {df.shape} \n')
nb_of_sample_expected=int(sys.argv[1])
variability_threshold=int(sys.argv[2])
keio_genes=list(pd.read_csv('/home/def-jacquesp/jeanneau/ReferenceData/Annotations/BW25113_features_supplemented.bed',sep='\t', names=(list('abcdef')))['d'])
df_keiogenes = df[[x for x in df.columns if x in keio_genes]]
with open('treatment.log','a') as logfile: logfile.write(f'[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}]  - kept only keio genes {df_keiogenes.shape}\n')
df_keiogenes_subsampled = df_keiogenes.loc[
    random.sample(list(df_keiogenes.index), min(int(nb_of_sample_expected*1.5),df.shape[0]))
    ].fillna(0).astype('uint8')
with open('treatment.log','a') as logfile: logfile.write(f'[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}]  - drawn first subset of samples {df_keiogenes_subsampled.shape}\n')
X = df_keiogenes_subsampled.values
row_names = df_keiogenes_subsampled.index
row_sums = X.sum(axis=1)
mask = row_sums > variability_threshold # control variability?
X_filtered = X[mask]
names_filtered = row_names[mask]
df_keiogenes_subsampled = pd.DataFrame(X_filtered, index=names_filtered, columns=df_keiogenes_subsampled.columns)
with open('treatment.log','a') as logfile: logfile.write(f'[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}]  - kept sample of >={variability_threshold} genes {df_keiogenes_subsampled.shape}\n')
df_keiogenes_subsampled = df_keiogenes_subsampled[[c for c in df_keiogenes_subsampled.columns if df_keiogenes_subsampled[c].sum()>0]]
included_genes = list(df_keiogenes_subsampled.columns)
with open('treatment.log','a') as logfile: logfile.write(f'[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}]  - kept non-empty genes {df_keiogenes_subsampled.shape}\n')
included_samples = random.sample(list(df_keiogenes_subsampled.index), min(nb_of_sample_expected, len(df_keiogenes_subsampled)))
df_keiogenes_subsampled = df_keiogenes_subsampled.loc[included_samples].fillna(0).astype(int)
included_samples = list(df_keiogenes_subsampled.index)
with open('treatment.log','a') as logfile: logfile.write(f'[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}]  - drawn final subset of samples {df_keiogenes_subsampled.shape}\n')
strain_labels = df_keiogenes_subsampled.index.to_numpy()
gene_labels = df_keiogenes_subsampled.columns.to_numpy()
arr_keiogenes = np.array(df_keiogenes_subsampled, dtype=bool)
arr_keiogenes, unique_idx = np.unique(arr_keiogenes, axis=0, return_index=True)
strain_labels = strain_labels[unique_idx]
with open("strain_labels.pkl", "wb") as f: pickle.dump(strain_labels, f)
with open("gene_labels.pkl", "wb") as f: pickle.dump(gene_labels, f)
with open("arr_keiogenes.pkl", "wb") as f: pickle.dump(arr_keiogenes, f) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<CHECKPOINT
with open('treatment.log','a') as logfile: logfile.write(f'[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}] Matricification, redundancy removal and Pickling done {arr_keiogenes.shape}\n')
strains, genes = arr_keiogenes.shape
frequencies = sorted(list(arr_keiogenes.sum(axis=0)))
plt.figure(figsize=(10,4))
sns.histplot(frequencies,bins=25)
plt.xlabel(f'Number of strains in which genes occur ({strains} strains)')
plt.ylabel(f'Gene occurrence ({genes} genes)')
plt.savefig(figpath+'Fig1_GeneOccurrenceHist.png')
plt.close()
with open('treatment.log','a') as logfile: logfile.write(f'[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}] Fig1_GeneOccurrenceHist written\n')
np.infty = np.inf
# Initialize UMAP with appropriate distance for binary vectors
reducer = umap.UMAP(
    metric='hamming',         # Better suited for binary data
    n_neighbors=4,           # Controls local structure preservation
    min_dist=0.5,             # Controls spread of clusters
    n_components=2,           # 2D projection
    random_state=42
)
# Fit and transform the data
embedding = reducer.fit_transform(arr_keiogenes)
with open("embedding.pkl", "wb") as f: pickle.dump(embedding, f) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<CHECKPOINT
with open('treatment.log','a') as logfile: logfile.write(f'[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}] UMAP embedding completed\n')
clusterer = hdbscan.HDBSCAN(
    metric='hamming',
    min_cluster_size=2,              # ← allow very small clusters
#    min_samples=2,                   # ← minimum density requirement
    cluster_selection_epsilon=0.005,   # ← allow larger radius for cluster separation
#    cluster_selection_method='leaf', # ← finer granularity (splits parent clusters)
    gen_min_span_tree=True,
    core_dist_n_jobs=72
)
clusterer.fit(arr_keiogenes)
with open("clusterer.pkl", "wb") as f: pickle.dump(clusterer, f) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<CHECKPOINT
with open('treatment.log','a') as logfile: logfile.write(f'[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}] HDBSCAN clustering completed\n')
labels = clusterer.labels_
n_clusters = labels.max() + 1
# Identify noise points (label == -1)
noise_mask = labels == -1
n_noise = noise_mask.sum()
# Assign new unique cluster IDs to each noise point
labels_fixed = labels.copy()
labels_fixed[noise_mask] = np.arange(n_clusters, n_clusters + n_noise)
sns.scatterplot(
    sorted(np.log10(list(cnt(labels_fixed).values())),reverse=True),
    linewidth=0.1)
plt.xlabel('Cluster label')
plt.ylabel('log10(Cluster size)')
_=plt.text(20000,3,f"{n_clusters} clusters",fontsize=28)
_=plt.text(20000,2.2,f"{n_noise} noises",fontsize=28)
plt.tight_layout()
plt.savefig(figpath+'Fig2_ClusterSizes.png')
plt.close()
with open('treatment.log','a') as logfile: logfile.write(f'[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}] Fig2_ClusterSizes written\n')
embedding_wlabel = np.column_stack((embedding, labels_fixed))
plt.figure(figsize=(6, 6))
sns.scatterplot(x=embedding_wlabel[:, 0], y=embedding_wlabel[:, 1], s=5, alpha=0.5, hue=embedding_wlabel[:, -1], palette='tab20', legend=False)
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
plt.tight_layout()
plt.savefig(figpath+'Fig3_UMAPembedding_HDBSCANcolored.png',dpi=300)
plt.close()
with open('treatment.log','a') as logfile: logfile.write(f'[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}] Fig3_UMAPembedding_HDBSCANcolored written\n')
representatives = []
cluster_representant_dict = {}
for cluster_id in np.unique(labels_fixed):
    #if cluster_id == -1:
    #    continue  # Skip noise
    cluster_mask = labels_fixed == cluster_id
    cluster_points = arr_keiogenes[cluster_mask]
    # Compute pairwise Hamming distances within the cluster
    dists = pairwise_distances(cluster_points, metric='hamming')
    # Find the medoid: index with smallest total distance
    medoid_local_idx = np.argmin(dists.sum(axis=1))
    # Map back to global index
    global_idx = np.where(cluster_mask)[0][medoid_local_idx]
    representatives.append(global_idx)
    cluster_representant_dict[cluster_id]=strain_labels[global_idx]
clustered_arr_keiogenes = arr_keiogenes[representatives]
with open("cluster_representant_dict.pkl", "wb") as f: pickle.dump(cluster_representant_dict, f) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<CHECKPOINT
with open("clustered_arr_keiogenes.pkl", "wb") as f: pickle.dump(clustered_arr_keiogenes, f) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<CHECKPOINT
with open('treatment.log','a') as logfile: logfile.write(f'[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}] Representant isolated {clustered_arr_keiogenes.shape}\n')
# this step is reputed UNFEASIBLE on large object ~70k x 4k
#sys.setrecursionlimit(100000)
#row_linkage = linkage(pdist(clustered_arr_keiogenes, metric="hamming"), method="average")
#col_linkage = linkage(pdist(clustered_arr_keiogenes.T, metric="hamming"), method="average")
#g = sns.clustermap(clustered_arr_keiogenes,
#                   row_linkage=row_linkage,
#                   col_linkage=col_linkage,
#                   figsize=(12, 6),
#                   cbar_pos=None)
#g.ax_heatmap.set_xlabel("Genes")
#g.ax_heatmap.set_ylabel("Strains")
#plt.savefig(figpath + 'Fig4_RepresentantClustermap.png', dpi=300, bbox_inches='tight')
#plt.close()
#with open('treatment.log','a') as logfile: logfile.write(f'[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}] Fig4_RepresentantClustermap written\n JOB DONE')
