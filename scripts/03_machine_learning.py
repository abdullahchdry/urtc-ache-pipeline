#!/usr/bin/env python3
"""
03_umap_clustering.py

Description:
    Run UMAP embedding and K-means clustering on ADMET-filtered compounds,
    select cluster medoids, compute clustering metrics, and save outputs.

Usage:
    python scripts/03_umap_clustering.py

Requirements:
    - pandas
    - numpy
    - scikit-learn
    - umap-learn
    - matplotlib

Data structure:
    data/processed/       # contains admet_filtered.csv and stores clustering outputs (CSVs)
    figures/              # contains plots (UMAP + clusters)

Before running:
    Ensure previous filtering step has produced data/processed/admet_filtered.csv.
"""

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from umap import UMAP
import matplotlib.pyplot as plt

# Define project directories
BASE_DIR = Path(__file__).resolve().parent.parent
PROCESSED_DIR = BASE_DIR / "data" / "processed"
FIGURES_DIR = BASE_DIR / "figures"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# File paths
INPUT_FILE      = PROCESSED_DIR / "admet_filtered.csv"
OUTPUT_ALL      = PROCESSED_DIR / "admet_umap_clusters.csv"
OUTPUT_MEDOIDS  = PROCESSED_DIR / "cluster_medoids.csv"
OUTPUT_METRICS  = PROCESSED_DIR / "clustering_metrics.csv"
OUTPUT_SIZES    = PROCESSED_DIR / "cluster_sizes.csv"
OUTPUT_PLOT     = FIGURES_DIR / "umap_clusters.png"

# Parameters
FEATURES = [
    'logS','logP','TPSA','nRot','BBB','MW','nHD','nHA',
    'logVDss','cl-plasma','t0.5','PPB','Fsp3'
]
N_CLUSTERS = 12

# Logging setup
logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger(__name__)


def main():
    # Verify input exists
    if not INPUT_FILE.exists():
        log.error(f"Error: filtered data not found at {INPUT_FILE}")
        sys.exit(1)

    # Load data and standardize features
    df = pd.read_csv(INPUT_FILE)
    X = StandardScaler().fit_transform(df[FEATURES])
    log.info(f"Loaded {len(df)} compounds for clustering")

    # UMAP embedding
    umap_model = UMAP(n_components=2, n_neighbors=15, min_dist=0.1, random_state=42)
    coords = umap_model.fit_transform(X)
    df[['UMAP1','UMAP2']] = coords
    df.to_csv(OUTPUT_ALL, index=False)
    log.info(f"UMAP coordinates saved to {OUTPUT_ALL}")

    # K-means clustering
    km = KMeans(n_clusters=N_CLUSTERS, random_state=42)
    df['cluster'] = km.fit_predict(coords)

    # Compute clustering metrics
    sil_score = silhouette_score(coords, df['cluster'])
    sizes = df['cluster'].value_counts().sort_index()
    metrics = {
        'silhouette_score': sil_score,
        'mean_cluster_size': sizes.mean(),
        'min_cluster_size': sizes.min(),
        'max_cluster_size': sizes.max()
    }
    pd.DataFrame([metrics]).to_csv(OUTPUT_METRICS, index=False)
    log.info(f"Clustering metrics saved to {OUTPUT_METRICS}")

    # Save cluster sizes
    cluster_sizes_df = pd.DataFrame({
        'cluster': sizes.index + 1,
        'size': sizes.values
    })
    cluster_sizes_df.to_csv(OUTPUT_SIZES, index=False)
    log.info(f"Cluster sizes saved to {OUTPUT_SIZES}")

    # Select medoids for each cluster
    medoids = []
    centers = km.cluster_centers_
    for cid in range(N_CLUSTERS):
        subset = df[df['cluster'] == cid]
        dists = np.linalg.norm(subset[['UMAP1','UMAP2']].values - centers[cid], axis=1)
        medoid = subset.iloc[np.argmin(dists)].copy()
        medoid['MedoidID'] = cid + 1
        medoids.append(medoid)
    medoids_df = pd.DataFrame(medoids)
    medoids_df.to_csv(OUTPUT_MEDOIDS, index=False)
    log.info(f"Cluster medoids saved to {OUTPUT_MEDOIDS}")

    # Plot UMAP clusters with medoids
    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(
        df['UMAP1'], df['UMAP2'],
        c=df['cluster'], cmap='tab20',
        s=50, alpha=0.8
    )
    plt.scatter(
        medoids_df['UMAP1'], medoids_df['UMAP2'],
        facecolors='none', edgecolors='k',
        s=150, marker='X', linewidth=1.5
    )

    # Build legend
    handles, _ = scatter.legend_elements(prop="colors", num=N_CLUSTERS)
    labels = [f"Cluster {i}" for i in range(1, N_CLUSTERS+1)]
    medoid_handle = plt.Line2D(
        [], [], marker='X', linestyle='None',
        markersize=10, markeredgewidth=1.5,
        markeredgecolor='k', markerfacecolor='none'
    )
    handles.append(medoid_handle)
    labels.append('Medoids')
    plt.legend(
        handles, labels,
        title="Clusters",
        bbox_to_anchor=(1.02, 1), loc='upper left',
        fontsize=8, title_fontsize=10
    )

    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.title('UMAP Projection with K-Means Clusters')
    plt.tight_layout(rect=(0, 0, 0.85, 1))
    plt.savefig(OUTPUT_PLOT, dpi=300)
    plt.close()
    log.info(f"UMAP cluster plot saved to {OUTPUT_PLOT}")

if __name__ == "__main__":
    main()
