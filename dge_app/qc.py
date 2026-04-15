from __future__ import annotations

import numpy as np
import pandas as pd


def build_sample_qc_table(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    gene_id_column: str,
    sample_id_column: str,
    sample_columns: list[str],
) -> pd.DataFrame:
    count_matrix = counts_df.set_index(gene_id_column)[sample_columns]

    qc = pd.DataFrame(
        {
            sample_id_column: sample_columns,
            "library_size": count_matrix.sum(axis=0).values,
            "detected_genes": (count_matrix > 0).sum(axis=0).values,
            "median_count": count_matrix.median(axis=0).values,
            "zero_fraction": ((count_matrix == 0).sum(axis=0) / count_matrix.shape[0]).values,
        }
    )
    return metadata_df.merge(qc, on=sample_id_column, how="left")


def compute_log_cpm(counts_df: pd.DataFrame, gene_id_column: str, sample_columns: list[str]) -> pd.DataFrame:
    count_matrix = counts_df.set_index(gene_id_column)[sample_columns].astype(float)
    lib_sizes = count_matrix.sum(axis=0)
    cpm = (count_matrix / lib_sizes) * 1_000_000
    log_cpm = np.log2(cpm + 1.0)
    return log_cpm


def compute_pca(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    gene_id_column: str,
    sample_id_column: str,
    sample_columns: list[str],
) -> tuple[pd.DataFrame, list[float]]:
    log_cpm = compute_log_cpm(counts_df, gene_id_column, sample_columns)
    matrix = log_cpm.transpose()
    centered = matrix - matrix.mean(axis=0)

    u, s, _ = np.linalg.svd(centered.to_numpy(), full_matrices=False)
    total_variance = np.square(s).sum()
    explained = (
        [float((value**2 / total_variance) * 100.0) for value in s[:2]]
        if total_variance > 0
        else [0.0, 0.0]
    )
    coords = u[:, :2] * s[:2]
    pca_df = pd.DataFrame(
        {
            sample_id_column: matrix.index,
            "PC1": coords[:, 0] if coords.shape[1] > 0 else 0.0,
            "PC2": coords[:, 1] if coords.shape[1] > 1 else 0.0,
        }
    )
    return metadata_df.merge(pca_df, on=sample_id_column, how="left"), explained


def compute_sample_distance_matrix(
    counts_df: pd.DataFrame,
    gene_id_column: str,
    sample_columns: list[str],
) -> pd.DataFrame:
    log_cpm = compute_log_cpm(counts_df, gene_id_column, sample_columns).transpose()
    samples = list(log_cpm.index)
    matrix = log_cpm.to_numpy()
    diffs = matrix[:, None, :] - matrix[None, :, :]
    distances = np.sqrt(np.square(diffs).sum(axis=2))
    return pd.DataFrame(distances, index=samples, columns=samples)
