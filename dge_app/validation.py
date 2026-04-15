from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import pandas as pd


@dataclass
class PreparedInputs:
    counts_df: pd.DataFrame
    metadata_df: pd.DataFrame
    sample_columns: list[str]
    design_factor: str
    reference_level: str
    contrast_level: str
    covariate_columns: list[str] = field(default_factory=list)


def coerce_numeric_counts(
    counts_df: pd.DataFrame,
    gene_id_column: str,
    sample_columns: list[str],
) -> tuple[pd.DataFrame, list[str]]:
    """
    Validate and coerce the counts matrix.

    Returns
    -------
    tuple[pd.DataFrame, list[str]]
        Cleaned counts DataFrame and a list of duplicate gene IDs that were removed.
    """
    if gene_id_column not in counts_df.columns:
        raise ValueError(
            f"Counts matrix is missing the gene identifier column '{gene_id_column}'."
        )

    if not sample_columns:
        raise ValueError("Counts matrix must include at least one sample count column.")

    missing_columns = [col for col in sample_columns if col not in counts_df.columns]
    if missing_columns:
        example = ", ".join(f"'{c}'" for c in missing_columns[:8])
        extra = f" … ({len(missing_columns) - 8} more)" if len(missing_columns) > 8 else ""
        raise ValueError(
            f"The following sample columns were not found in the counts matrix: {example}{extra}"
        )

    coerced = counts_df.copy()
    coerced[sample_columns] = coerced[sample_columns].apply(pd.to_numeric, errors="coerce")
    if coerced[sample_columns].isna().any().any():
        raise ValueError(
            "Counts matrix contains non-numeric values in the selected sample columns. "
            "Ensure you have selected only raw count columns and not accidental text columns."
        )
    if (coerced[sample_columns] < 0).any().any():
        raise ValueError("Counts matrix cannot contain negative values.")

    # Track duplicates before dropping
    dup_mask = coerced.duplicated(subset=[gene_id_column], keep="first")
    duplicate_genes: list[str] = coerced.loc[dup_mask, gene_id_column].astype(str).tolist()
    coerced = coerced.drop_duplicates(subset=[gene_id_column]).reset_index(drop=True)
    return coerced, duplicate_genes


def parse_annotation_text(annotation_df: pd.DataFrame, gene_id_column: str) -> pd.DataFrame:
    if gene_id_column not in annotation_df.columns:
        raise ValueError(
            f"Annotation table is missing the gene identifier column '{gene_id_column}'."
        )
    return annotation_df.drop_duplicates(subset=[gene_id_column]).copy()


def prepare_inputs(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    gene_id_column: str,
    sample_id_column: str,
    sample_columns: list[str],
    design_factor: str,
    reference_level: str,
    contrast_level: str,
    covariate_columns: list[str] | None = None,
) -> PreparedInputs:
    if covariate_columns is None:
        covariate_columns = []

    if sample_id_column not in metadata_df.columns:
        raise ValueError(
            f"Metadata table is missing the sample identifier column '{sample_id_column}'."
        )
    if design_factor not in metadata_df.columns:
        raise ValueError(
            f"Metadata table is missing the design column '{design_factor}'."
        )
    for cov in covariate_columns:
        if cov not in metadata_df.columns:
            raise ValueError(
                f"Metadata table is missing the covariate column '{cov}'."
            )

    metadata_sample_ids = set(metadata_df[sample_id_column].astype(str))
    sample_columns = [str(col) for col in sample_columns]
    counts_df = counts_df.copy()
    counts_df.columns = [str(col) for col in counts_df.columns]
    metadata_df = metadata_df.copy()
    metadata_df[sample_id_column] = metadata_df[sample_id_column].astype(str)

    missing_samples = sorted(set(sample_columns) - metadata_sample_ids)
    if missing_samples:
        raise ValueError(
            "Metadata table is missing sample annotations for: "
            + ", ".join(missing_samples[:12])
            + ("…" if len(missing_samples) > 12 else "")
        )

    metadata_prepared = metadata_df[metadata_df[sample_id_column].isin(sample_columns)].copy()
    metadata_prepared = metadata_prepared.drop_duplicates(subset=[sample_id_column])
    metadata_prepared = (
        metadata_prepared.set_index(sample_id_column).loc[sample_columns].reset_index()
    )

    levels = set(metadata_prepared[design_factor].astype(str))
    if reference_level not in levels:
        raise ValueError(
            f"Reference group '{reference_level}' was not found in '{design_factor}'. "
            f"Available levels: {sorted(levels)}"
        )
    if contrast_level not in levels:
        raise ValueError(
            f"Comparison group '{contrast_level}' was not found in '{design_factor}'. "
            f"Available levels: {sorted(levels)}"
        )
    if reference_level == contrast_level:
        raise ValueError("Reference and comparison groups must be different.")

    if metadata_prepared[design_factor].nunique() < 2:
        raise ValueError("The design column must contain at least two groups.")

    return PreparedInputs(
        counts_df=counts_df,
        metadata_df=metadata_prepared,
        sample_columns=sample_columns,
        design_factor=design_factor,
        reference_level=reference_level,
        contrast_level=contrast_level,
        covariate_columns=list(covariate_columns),
    )
