from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import pandas as pd

R_IMPORT_ERROR_MESSAGE = (
    "The analysis backend requires `rpy2` plus an R installation with "
    "`DESeq2`, `edgeR`, and `limma` available. See the README for setup steps."
)


class DGEAnalysisError(RuntimeError):
    """Raised when the DGE backend cannot run or returns invalid output."""


@dataclass
class DGEResult:
    method: str
    results: pd.DataFrame
    normalized_counts: pd.DataFrame
    summary: dict
    session_info: str


@dataclass
class BackendStatus:
    available: bool
    summary: str
    details: list[str]


METHOD_PACKAGE_MAP = {
    "DESeq2": ["DESeq2", "tibble"],
    "edgeR": ["edgeR", "tibble"],
    "limma-voom": ["edgeR", "limma", "tibble"],
}

# Module-level cache — prevents re-initialising R on every Streamlit rerun,
# which was causing the "promise already under evaluation" crash and SIGINT warnings.
_BACKEND_STATUS_CACHE: dict[str, BackendStatus] | None = None


def list_available_methods() -> list[str]:
    return list(METHOD_PACKAGE_MAP.keys())


def get_backend_status(force_refresh: bool = False) -> dict[str, BackendStatus]:
    global _BACKEND_STATUS_CACHE
    if _BACKEND_STATUS_CACHE is not None and not force_refresh:
        return _BACKEND_STATUS_CACHE

    methods = list_available_methods()
    try:
        robjects, pandas2ri, importr = _import_r_dependencies()
    except DGEAnalysisError as exc:
        message = str(exc)
        _BACKEND_STATUS_CACHE = {
            m: BackendStatus(
                available=False,
                summary="Python-to-R bridge (rpy2) is unavailable.",
                details=[message],
            )
            for m in methods
        }
        return _BACKEND_STATUS_CACHE

    try:
        with _r_conversion_context(robjects, pandas2ri):
            importr("base")
            importr("utils")
            importr("stats")
    except Exception as exc:
        message = "R is installed but could not be loaded through rpy2."
        _BACKEND_STATUS_CACHE = {
            m: BackendStatus(available=False, summary=message, details=[message, str(exc)])
            for m in methods
        }
        return _BACKEND_STATUS_CACHE

    statuses: dict[str, BackendStatus] = {}
    for method, package_names in METHOD_PACKAGE_MAP.items():
        missing: list[str] = []
        for pkg in package_names:
            try:
                with _r_conversion_context(robjects, pandas2ri):
                    importr(pkg)
            except Exception:
                missing.append(pkg)

        if missing:
            statuses[method] = BackendStatus(
                available=False,
                summary=f"Missing R packages: {', '.join(missing)}",
                details=[
                    f"Install required packages for {method}: {', '.join(package_names)}.",
                    "Run `Rscript .\\install_r_packages.R` in the environment used by Streamlit.",
                ],
            )
        else:
            statuses[method] = BackendStatus(
                available=True,
                summary="Ready",
                details=[f"{method} backend is available."],
            )

    _BACKEND_STATUS_CACHE = statuses
    return statuses


def _build_r_script(method: str) -> str:
    """
    Build a self-contained R function string.

    Arguments passed from Python:
        counts_df, metadata_df, gene_id_col, sample_id_col,
        design_factor, reference_level, contrast_level,
        covariate_columns_str (comma-separated), alpha_value
    """
    method_key = method.lower()

    if method_key == "deseq2":
        body = """
            suppressPackageStartupMessages(library(DESeq2))
            dds <- DESeqDataSetFromMatrix(
              countData = counts_mat,
              colData   = metadata_df,
              design    = design_formula
            )
            dds[[design_factor]] <- relevel(as.factor(dds[[design_factor]]), ref = reference_level)
            deseq_error <- NULL
            dds <- tryCatch(
              DESeq(dds),
              error = function(e) {{ deseq_error <<- conditionMessage(e); dds }}
            )
            if (!is.null(deseq_error)) {{
              if (grepl("all gene-wise dispersion estimates are within 2 orders", deseq_error, fixed = TRUE)) {{
                dds <- estimateSizeFactors(dds)
                dds <- estimateDispersionsGeneEst(dds)
                dispersions(dds) <- mcols(dds)$dispGeneEst
                dds <- nbinomWaldTest(dds)
              }} else {{
                stop(deseq_error)
              }}
            }}
            result_tbl        <- as.data.frame(results(dds,
                                   contrast  = c(design_factor, contrast_level, reference_level),
                                   alpha     = alpha_value))
            result_tbl$base_mean        <- result_tbl$baseMean
            result_tbl$log2_fold_change <- result_tbl$log2FoldChange
            result_tbl$statistic        <- result_tbl$stat
            result_tbl$p_value          <- result_tbl$pvalue
            result_tbl$adjusted_p_value <- result_tbl$padj
            normalized_counts <- as.data.frame(counts(dds, normalized = TRUE))
        """
    elif method_key == "edger":
        body = """
            suppressPackageStartupMessages(library(edgeR))
            metadata_df[[design_factor]] <- relevel(as.factor(metadata_df[[design_factor]]), ref = reference_level)
            design_matrix <- model.matrix(design_formula, data = metadata_df)
            coef_name     <- tail(colnames(design_matrix), 1)
            y  <- DGEList(counts = counts_mat)
            y  <- calcNormFactors(y)
            y  <- estimateDisp(y, design_matrix)
            fit <- glmQLFit(y, design_matrix)
            qlf <- glmQLFTest(fit, coef = coef_name)
            result_tbl <- topTags(qlf, n = nrow(counts_mat), sort.by = "none")$table
            result_tbl$base_mean        <- rowMeans(cpm(y, log = FALSE))
            result_tbl$log2_fold_change <- result_tbl$logFC
            result_tbl$statistic        <- result_tbl$F
            result_tbl$p_value          <- result_tbl$PValue
            result_tbl$adjusted_p_value <- result_tbl$FDR
            normalized_counts <- as.data.frame(cpm(y, log = FALSE, normalized.lib.sizes = TRUE))
        """
    elif method_key == "limma-voom":
        body = """
            suppressPackageStartupMessages(library(edgeR))
            suppressPackageStartupMessages(library(limma))
            metadata_df[[design_factor]] <- relevel(as.factor(metadata_df[[design_factor]]), ref = reference_level)
            design_matrix <- model.matrix(design_formula, data = metadata_df)
            coef_name     <- tail(colnames(design_matrix), 1)
            y   <- DGEList(counts = counts_mat)
            y   <- calcNormFactors(y)
            v   <- voom(y, design = design_matrix, plot = FALSE)
            fit <- lmFit(v, design_matrix)
            fit <- eBayes(fit)
            result_tbl <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
            result_tbl$base_mean        <- rowMeans(v$E)
            result_tbl$log2_fold_change <- result_tbl$logFC
            result_tbl$statistic        <- result_tbl$t
            result_tbl$p_value          <- result_tbl$P.Value
            result_tbl$adjusted_p_value <- result_tbl$adj.P.Val
            normalized_counts <- as.data.frame(v$E)
        """
    else:
        raise DGEAnalysisError(f"Unsupported method '{method}'.")

    return f"""
        function(counts_df, metadata_df,
                 gene_id_col, sample_id_col,
                 design_factor, reference_level, contrast_level,
                 covariate_columns_str, alpha_value) {{
          suppressWarnings({{

            # ── Counts matrix ──────────────────────────────────────────────
            gene_ids  <- as.character(counts_df[[gene_id_col]])
            count_cols <- colnames(counts_df)[colnames(counts_df) != gene_id_col]
            counts_mat <- as.matrix(counts_df[, count_cols, drop = FALSE])
            rownames(counts_mat) <- gene_ids
            mode(counts_mat) <- "numeric"

            # ── Metadata ───────────────────────────────────────────────────
            metadata_df[[sample_id_col]] <- as.character(metadata_df[[sample_id_col]])
            rownames(metadata_df) <- metadata_df[[sample_id_col]]
            metadata_df <- metadata_df[, colnames(metadata_df) != sample_id_col, drop = FALSE]

            # ── Design formula (supports multiple covariates) ──────────────
            covariate_vec <- strsplit(covariate_columns_str, ",")[[1]]
            covariate_vec <- covariate_vec[nchar(trimws(covariate_vec)) > 0]
            formula_terms  <- c(covariate_vec, design_factor)
            design_formula <- as.formula(paste("~", paste(formula_terms, collapse = " + ")))

            {body}

            result_tbl        <- tibble::rownames_to_column(result_tbl,        var = "gene_id")
            normalized_counts <- tibble::rownames_to_column(normalized_counts, var = "gene_id")
            list(
              results           = result_tbl,
              normalized_counts = normalized_counts,
              session_info      = paste(capture.output(sessionInfo()), collapse = "\\n")
            )
          }})
        }}
    """


def _import_r_dependencies():
    try:
        from rpy2 import robjects
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.packages import importr
    except Exception as exc:
        raise DGEAnalysisError(R_IMPORT_ERROR_MESSAGE) from exc
    return robjects, pandas2ri, importr


def _r_conversion_context(robjects, pandas2ri):
    return (robjects.default_converter + pandas2ri.converter).context()


def run_dge_analysis(
    method: str,
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    gene_id_column: str,
    sample_id_column: str,
    sample_columns: list[str],
    design_factor: str,
    reference_level: str,
    contrast_level: str,
    covariate_columns: list[str],
    alpha: float,
    lfc_cutoff: float,
    annotation_df: Optional[pd.DataFrame] = None,
) -> DGEResult:
    # Subset and reorder to ensure column positions are deterministic
    counts_df = counts_df[[gene_id_column] + sample_columns]
    if sample_id_column != metadata_df.columns[0]:
        reordered = [sample_id_column] + [c for c in metadata_df.columns if c != sample_id_column]
        metadata_df = metadata_df[reordered]

    robjects, pandas2ri, importr = _import_r_dependencies()

    try:
        with _r_conversion_context(robjects, pandas2ri):
            importr("base")
            importr("utils")
            importr("stats")
            importr("tibble")
            if method == "DESeq2":
                importr("DESeq2")
            elif method == "edgeR":
                importr("edgeR")
            elif method == "limma-voom":
                importr("edgeR")
                importr("limma")
    except Exception as exc:
        raise DGEAnalysisError(
            f"Required R package import failed while preparing {method}. {R_IMPORT_ERROR_MESSAGE}"
        ) from exc

    with _r_conversion_context(robjects, pandas2ri):
        try:
            r_function = robjects.r(_build_r_script(method))
            result = r_function(
                counts_df,
                metadata_df,
                gene_id_column,
                sample_id_column,
                design_factor,
                reference_level,
                contrast_level,
                ",".join(covariate_columns) if covariate_columns else "",
                alpha,
            )
            raw_results           = _extract_named_result(result, "results")
            raw_normalized_counts = _extract_named_result(result, "normalized_counts")
            raw_session_info      = _extract_named_result(result, "session_info")

            results_df        = _to_pandas_dataframe(raw_results, pandas2ri)
            normalized_counts = _to_pandas_dataframe(raw_normalized_counts, pandas2ri)
            if isinstance(raw_session_info, str):
                session_info = raw_session_info
            elif hasattr(raw_session_info, "__getitem__"):
                session_info = str(raw_session_info[0])
            else:
                session_info = str(raw_session_info)
        except Exception as exc:
            raise DGEAnalysisError(
                f"{method} failed to run. Check that your count matrix uses raw integer-like counts "
                f"and that your metadata columns and contrast are valid. Backend error: {exc}"
            ) from exc

    results_df = _tidy_results(results_df, gene_id_column, annotation_df, alpha, lfc_cutoff)
    normalized_counts = normalized_counts.rename(
        columns={"gene_id": gene_id_column}
    ).set_index(gene_id_column)

    return DGEResult(
        method=method,
        results=results_df,
        normalized_counts=normalized_counts,
        summary=_summarize_results(results_df),
        session_info=session_info,
    )


def _tidy_results(
    results_df: pd.DataFrame,
    gene_id_column: str,
    annotation_df: Optional[pd.DataFrame],
    alpha: float,
    lfc_cutoff: float,
) -> pd.DataFrame:
    results_df = results_df.rename(columns={"gene_id": gene_id_column})

    keep_cols = [gene_id_column, "base_mean", "log2_fold_change", "statistic",
                 "p_value", "adjusted_p_value"]
    existing = [c for c in keep_cols if c in results_df.columns]
    tidy = results_df[existing].copy()

    tidy["p_value"]          = pd.to_numeric(tidy["p_value"],          errors="coerce")
    tidy["adjusted_p_value"] = pd.to_numeric(tidy["adjusted_p_value"], errors="coerce")
    tidy["log2_fold_change"] = pd.to_numeric(tidy["log2_fold_change"], errors="coerce")
    tidy["base_mean"]        = pd.to_numeric(tidy["base_mean"],        errors="coerce")
    if "statistic" in tidy.columns:
        tidy["statistic"] = pd.to_numeric(tidy["statistic"], errors="coerce")

    # Track which genes have NaN adjusted p-value BEFORE filling (DESeq2 independent filtering)
    originally_nan_padj = tidy["adjusted_p_value"].isna()

    tidy["is_significant"] = (
        tidy["adjusted_p_value"].fillna(1.0).le(alpha)
        & tidy["log2_fold_change"].abs().ge(lfc_cutoff)
    )

    # Assign direction — "Filtered" reserved for DESeq2-filtered genes (NaN padj)
    tidy["direction"] = "Not significant"
    tidy.loc[originally_nan_padj, "direction"] = "Filtered"
    tidy.loc[tidy["is_significant"] & tidy["log2_fold_change"].gt(0), "direction"] = "Upregulated"
    tidy.loc[tidy["is_significant"] & tidy["log2_fold_change"].lt(0), "direction"] = "Downregulated"

    tidy = tidy.sort_values(["adjusted_p_value", "p_value", "log2_fold_change"], na_position="last")

    if annotation_df is not None:
        tidy = tidy.merge(annotation_df, on=gene_id_column, how="left")

    return tidy.reset_index(drop=True)


def apply_thresholds(res_df: pd.DataFrame, alpha: float, lfc_cutoff: float) -> pd.DataFrame:
    """Re-evaluate is_significant and direction live from current thresholds.

    Use this in the UI to update display results without re-running the analysis.
    """
    df = res_df.copy()
    originally_nan = res_df["adjusted_p_value"].isna()
    df["is_significant"] = (
        df["adjusted_p_value"].fillna(1.0).le(alpha)
        & df["log2_fold_change"].abs().ge(lfc_cutoff)
    )
    df["direction"] = "Not significant"
    df.loc[originally_nan, "direction"] = "Filtered"
    df.loc[df["is_significant"] & df["log2_fold_change"].gt(0), "direction"] = "Upregulated"
    df.loc[df["is_significant"] & df["log2_fold_change"].lt(0), "direction"] = "Downregulated"
    return df


def _extract_named_result(result, key: str):
    if hasattr(result, "rx2"):
        return result.rx2(key)
    try:
        return result[key]
    except Exception:
        pass
    names_attr  = getattr(result, "names", None)
    names_value = names_attr() if callable(names_attr) else names_attr
    names       = list(names_value or [])
    if key in names:
        return result[names.index(key)]
    raise DGEAnalysisError(f"Expected analysis output key '{key}' was not found.")


def _to_pandas_dataframe(value, pandas2ri) -> pd.DataFrame:
    if isinstance(value, pd.DataFrame):
        return value.copy()
    converted = pandas2ri.rpy2py(value)
    if isinstance(converted, pd.DataFrame):
        return converted
    return pd.DataFrame(converted)


def _summarize_results(results_df: pd.DataFrame) -> dict:
    sig = results_df[results_df["is_significant"]]
    return {
        "genes_tested":     int(results_df.shape[0]),
        "significant_genes": int(sig.shape[0]),
        "upregulated":      int((sig["direction"] == "Upregulated").sum()),
        "downregulated":    int((sig["direction"] == "Downregulated").sum()),
        "filtered":         int((results_df["direction"] == "Filtered").sum()),
    }
