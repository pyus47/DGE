from __future__ import annotations

import json
import os
import subprocess
import sys
from datetime import datetime
from itertools import combinations
from typing import Any

import pandas as pd
import streamlit as st
from scipy.stats import hypergeom

from dge_app.analysis import (
    DGEAnalysisError,
    apply_thresholds,
    get_backend_status,
    list_available_methods,
    run_dge_analysis,
)
from dge_app.plots import (
    make_detected_genes_plot,
    make_gene_expression_plot,
    make_library_size_plot,
    make_lfc_scatter,
    make_log2fc_correlation_plot,
    make_ma_plot,
    make_method_comparison_table,
    make_normalization_boxplot,
    make_overlap_bar,
    make_pca_plot,
    make_pvalue_histogram,
    make_sample_distance_heatmap,
    make_top_gene_heatmap,
    make_volcano_plot,
)
from dge_app.qc import (
    build_sample_qc_table,
    compute_log_cpm,
    compute_pca,
    compute_sample_distance_matrix,
)
from dge_app.report import build_pdf_report
from dge_app.validation import coerce_numeric_counts, parse_annotation_text, prepare_inputs

st.set_page_config(page_title="DGE", page_icon="🧬", layout="wide")

local_r_path = os.path.join(os.getcwd(), "R_libs")
os.environ['R_LIBS_USER'] = local_r_path

# Function to run the R installation script
def verify_r_environment():
    # We check for a 'flag' file to avoid running the 10-minute install on every refresh
    if not os.path.exists("r_packages_installed.flag"):
        with st.spinner("Installing Bioconductor packages (DESeq2, edgeR)... This may take 10-15 minutes."):
            try:
                # This runs your .R file exactly like you would in a terminal
                subprocess.run(["Rscript", "install_r_packages.R"], check=True)
                # Create the flag file
                with open("r_packages_installed.flag", "w") as f:
                    f.write("Success")
                st.success("R environment ready!")
            except subprocess.CalledProcessError as e:
                st.error(f"R installation failed: {e}")
                st.stop()

verify_r_environment()


# ─────────────────────────────────────────────────────────────────────────────
# Session-state initialisation
# ─────────────────────────────────────────────────────────────────────────────

def _init_state() -> None:
    defaults: dict[str, Any] = {
        "counts_df":      None,
        "metadata_df":    None,
        "annotation_df":  None,
        "results":        {},
        "run_log":        [],
        "analysis_params": {"alpha": 0.05, "lfc_cutoff": 1.0, "top_n": 50},
    }
    for key, val in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = val


# ─────────────────────────────────────────────────────────────────────────────
# Cached helpers
# ─────────────────────────────────────────────────────────────────────────────

@st.cache_resource(show_spinner="Checking R backend…")
def _get_backend_status_cached() -> dict:
    """Cached once per server session — prevents repeated R initialisation crashes."""
    return get_backend_status()


@st.cache_data(show_spinner=False)
def _load_example_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    genes  = [f"Gene_{i:04d}" for i in range(1, 301)]
    counts = pd.DataFrame({
        "gene_id": genes,
        "ctrl_1":  [120 + i * 2 for i in range(300)],
        "ctrl_2":  [115 + i * 2 for i in range(300)],
        "ctrl_3":  [118 + i * 2 for i in range(300)],
        "treat_1": [122 + i * 2 if i < 220 else 220 + i * 3 for i in range(300)],
        "treat_2": [119 + i * 2 if i < 220 else 212 + i * 3 for i in range(300)],
        "treat_3": [121 + i * 2 if i < 220 else 216 + i * 3 for i in range(300)],
    })
    metadata = pd.DataFrame({
        "sample":    ["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3"],
        "condition": ["control"] * 3 + ["treated"] * 3,
        "batch":     ["A", "A", "B", "A", "B", "B"],
    })
    return counts, metadata


@st.cache_data(show_spinner=False)
def _compute_qc_cached(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    gene_id_column: str,
    sample_id_column: str,
    sample_columns: tuple[str, ...],
):
    sc = list(sample_columns)
    qc_df       = build_sample_qc_table(counts_df, metadata_df, gene_id_column, sample_id_column, sc)
    pca_df, exp = compute_pca(counts_df, metadata_df, gene_id_column, sample_id_column, sc)
    dist_df     = compute_sample_distance_matrix(counts_df, gene_id_column, sc)
    log_cpm     = compute_log_cpm(counts_df, gene_id_column, sc)
    return qc_df, pca_df, exp, dist_df, log_cpm


# ─────────────────────────────────────────────────────────────────────────────
# General helpers
# ─────────────────────────────────────────────────────────────────────────────

def _read_table(f) -> pd.DataFrame:
    name = f.name.lower()
    if name.endswith(".xlsx") or name.endswith(".xls"):
        return pd.read_excel(f)
    if name.endswith(".csv"):
        return pd.read_csv(f)
    return pd.read_csv(f, sep=None, engine="python")


def _download_csv(df: pd.DataFrame, label: str, filename: str) -> None:
    st.download_button(
        label=label,
        data=df.to_csv(index=False).encode("utf-8"),
        file_name=filename,
        mime="text/csv",
    )


def _live_results(method: str) -> pd.DataFrame:
    """Return the stored results for *method* with thresholds re-applied from current UI state."""
    result = st.session_state.results[method]
    alpha  = st.session_state.get("res_alpha",  st.session_state.analysis_params["alpha"])
    lfc    = st.session_state.get("res_lfc",    st.session_state.analysis_params["lfc_cutoff"])
    return apply_thresholds(result.results, alpha, lfc)


def _live_summary(res_df: pd.DataFrame) -> dict:
    sig = res_df[res_df["is_significant"]]
    return {
        "genes_tested":      int(res_df.shape[0]),
        "significant_genes": int(sig.shape[0]),
        "upregulated":       int((sig["direction"] == "Upregulated").sum()),
        "downregulated":     int((sig["direction"] == "Downregulated").sum()),
        "filtered":          int((res_df["direction"] == "Filtered").sum()),
    }


def _formula_str(design_factor: str, covariates: list[str]) -> str:
    if covariates:
        return f"~ {' + '.join(covariates)} + {design_factor}"
    return f"~ {design_factor}"


def _parse_gmt_text(text: str) -> dict[str, set[str]]:
    gene_sets: dict[str, set[str]] = {}
    for line in text.splitlines():
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        name  = parts[0]
        genes = {p for p in parts[2:] if p}
        if genes:
            gene_sets[name] = genes
    return gene_sets


def _run_enrichment(sig_genes: set[str], bg: set[str], gene_sets: dict[str, set[str]]) -> pd.DataFrame:
    rows = []
    M, N = len(bg), len(sig_genes)
    if M == 0 or N == 0:
        return pd.DataFrame()
    for name, members in gene_sets.items():
        members_bg = members & bg
        if not members_bg:
            continue
        n = len(members_bg)
        k = len(sig_genes & members_bg)
        if k == 0:
            continue
        rows.append({"term": name, "overlap": k, "set_size": n,
                     "p_value": float(hypergeom.sf(k - 1, M, n, N))})
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows).sort_values("p_value").reset_index(drop=True)
    out["rank"]           = range(1, len(out) + 1)
    out["adjusted_p_value"] = (out["p_value"] * len(out)).clip(upper=1.0)
    return out


def _save_run(params: dict, methods: list[str]) -> None:
    st.session_state.run_log.append({
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "methods":   methods,
        "params":    params,
        "python":    sys.version,
    })


def _summary_lines() -> list[str]:
    lines  = ["DGE Studio — Analysis Summary"]
    params = st.session_state.analysis_params
    alpha  = params.get("alpha", 0.05)
    lfc    = params.get("lfc_cutoff", 1.0)
    if st.session_state.results:
        for method, result in st.session_state.results.items():
            rdf = apply_thresholds(result.results, alpha, lfc)
            s   = _live_summary(rdf)
            lines.append(
                f"{method}: tested={s['genes_tested']}, sig={s['significant_genes']}, "
                f"up={s['upregulated']}, down={s['downregulated']}, filtered={s['filtered']}"
            )
    lines.append(f"Run log entries: {len(st.session_state.run_log)}")
    return lines


# ─────────────────────────────────────────────────────────────────────────────
# Tab 1 — Setup & Run
# ─────────────────────────────────────────────────────────────────────────────

def _render_setup(backend_status: dict) -> None:
    st.markdown(
        "Upload your count matrix and sample metadata, configure the experimental design, "
        "select analysis methods, then click **▶ Run Analysis**. "
        "You can run all three methods (DESeq2, edgeR, limma-voom) and compare their results "
        "in the **🔁 Method Comparison** tab."
    )

    data_loaded = (
        st.session_state.counts_df is not None
        and st.session_state.metadata_df is not None
    )

    # ── Section 1: Data ───────────────────────────────────────────────────────
    with st.expander("📁 Data Input", expanded=not data_loaded):
        if st.button("Load bundled example data (300 genes × 6 samples, 2 conditions)"):
            counts_df, metadata_df = _load_example_inputs()
            st.session_state.counts_df     = counts_df
            st.session_state.metadata_df   = metadata_df
            st.session_state.annotation_df = None
            st.success("Example dataset loaded.")
            st.rerun()

        st.divider()
        c1, c2, c3 = st.columns(3)
        with c1:
            count_file = st.file_uploader(
                "Count matrix", type=["csv", "tsv", "txt", "xlsx"],
                help="Rows = genes, columns = samples. First column should be gene IDs.",
                key="fu_counts",
            )
        with c2:
            meta_file = st.file_uploader(
                "Sample metadata", type=["csv", "tsv", "txt", "xlsx"],
                help="One row per sample. Must have a column matching count matrix column names.",
                key="fu_metadata",
            )
        with c3:
            annot_file = st.file_uploader(
                "Gene annotation (optional)", type=["csv", "tsv", "txt", "xlsx"],
                help="Optional table with gene symbols or descriptions.",
                key="fu_annotation",
            )

        if st.button("Apply uploaded files"):
            if count_file is None or meta_file is None:
                st.error("Both count matrix and sample metadata are required.")
            else:
                st.session_state.counts_df     = _read_table(count_file)
                st.session_state.metadata_df   = _read_table(meta_file)
                st.session_state.annotation_df = _read_table(annot_file) if annot_file else None
                st.success("Files loaded.")
                st.rerun()

        if data_loaded:
            p1, p2 = st.columns(2)
            with p1:
                counts = st.session_state.counts_df
                st.caption(f"**Count matrix** — {counts.shape[0]:,} genes × {counts.shape[1]} columns")
                st.dataframe(counts.head(5), use_container_width=True, hide_index=True)
            with p2:
                meta = st.session_state.metadata_df
                st.caption(f"**Sample metadata** — {meta.shape[0]} samples × {meta.shape[1]} columns")
                st.dataframe(meta.head(5), use_container_width=True, hide_index=True)

    # ── Section 2: Column Mapping & Design ────────────────────────────────────
    with st.expander("⚙️ Column Mapping & Experimental Design", expanded=data_loaded):
        if not data_loaded:
            st.info("Load data first.")
        else:
            counts_df   = st.session_state.counts_df.copy()
            metadata_df = st.session_state.metadata_df.copy()
            counts_df.columns   = counts_df.columns.map(str)
            metadata_df.columns = metadata_df.columns.map(str)

            counts_cols = list(counts_df.columns)
            meta_cols   = list(metadata_df.columns)

            st.markdown("**Column Mapping**")
            cm1, cm2, cm3 = st.columns(3)

            with cm1:
                default_gene    = "gene_id" if "gene_id" in counts_cols else counts_cols[0]
                prev_gene       = st.session_state.get("w_gene_id_column", default_gene)
                gene_idx        = counts_cols.index(prev_gene) if prev_gene in counts_cols else 0
                gene_id_col     = st.selectbox(
                    "Gene ID column", counts_cols, index=gene_idx,
                    key="w_gene_id_column",
                    help="Column in the count matrix that holds gene identifiers.",
                )

            with cm2:
                default_sample  = "sample" if "sample" in meta_cols else meta_cols[0]
                prev_sample     = st.session_state.get("w_sample_id_column", default_sample)
                sample_idx      = meta_cols.index(prev_sample) if prev_sample in meta_cols else 0
                sample_id_col   = st.selectbox(
                    "Sample ID column (metadata)", meta_cols, index=sample_idx,
                    key="w_sample_id_column",
                    help="Column in metadata whose values match count matrix column names.",
                )

            with cm3:
                metadata_ids    = set(metadata_df[sample_id_col].astype(str))
                candidates      = [c for c in counts_cols if c != gene_id_col]
                auto_selected   = [c for c in candidates if c in metadata_ids] or candidates
                prev_sc         = st.session_state.get("w_sample_columns", auto_selected)
                valid_sc        = [c for c in prev_sc if c in candidates]
                sample_cols     = st.multiselect(
                    "Sample count columns", options=candidates,
                    default=valid_sc if valid_sc else auto_selected,
                    key="w_sample_columns",
                    help="Which count matrix columns represent samples to include.",
                )

            if sample_cols:
                in_counts_not_meta = [s for s in sample_cols if s not in metadata_ids]
                in_meta_not_counts = [s for s in metadata_ids if s not in sample_cols]
                if in_counts_not_meta:
                    names = "`, `".join(in_counts_not_meta[:5])
                    st.warning(
                        f"⚠️ {len(in_counts_not_meta)} count column(s) have no metadata: `{names}`"
                        + ("…" if len(in_counts_not_meta) > 5 else "")
                    )
                if in_meta_not_counts:
                    names = "`, `".join(list(in_meta_not_counts)[:5])
                    st.info(
                        f"ℹ️ {len(in_meta_not_counts)} metadata sample(s) excluded: `{names}`"
                        + ("…" if len(in_meta_not_counts) > 5 else "")
                    )
                if not in_counts_not_meta:
                    st.success(f"✅ All {len(sample_cols)} selected samples match metadata.")

            st.divider()
            st.markdown("**Experimental Design**")

            design_candidates = [c for c in meta_cols if c != sample_id_col]
            if not design_candidates:
                st.error("Metadata must have at least one column beyond the sample ID column.")
            else:
                dd1, dd2 = st.columns(2)
                with dd1:
                    default_factor = "condition" if "condition" in design_candidates else design_candidates[0]
                    prev_factor    = st.session_state.get("w_design_factor", default_factor)
                    fac_idx        = (design_candidates.index(prev_factor)
                                      if prev_factor in design_candidates else 0)
                    design_factor  = st.selectbox(
                        "Primary design factor", design_candidates, index=fac_idx,
                        key="w_design_factor",
                        help="The main biological variable under test (e.g. treatment condition).",
                    )
                with dd2:
                    cov_options    = [c for c in design_candidates if c != design_factor]
                    prev_covs      = st.session_state.get("w_covariate_columns", [])
                    valid_covs     = [c for c in prev_covs if c in cov_options]
                    covariate_cols = st.multiselect(
                        "Covariates / Batch factors (optional)", options=cov_options,
                        default=valid_covs,
                        key="w_covariate_columns",
                        help="Additional variables to control for (e.g. batch, sex). Multi-factor designs supported.",
                    )

                if design_factor in metadata_df.columns:
                    levels = sorted(metadata_df[design_factor].dropna().astype(str).unique().tolist())
                    if len(levels) >= 2:
                        lc1, lc2 = st.columns(2)
                        with lc1:
                            default_ref = "control" if "control" in levels else levels[0]
                            prev_ref    = st.session_state.get("w_reference_level", default_ref)
                            ref_idx     = levels.index(prev_ref) if prev_ref in levels else 0
                            ref_level   = st.selectbox(
                                "Reference group", levels, index=ref_idx,
                                key="w_reference_level",
                                help="Baseline / control group (denominator of fold change).",
                            )
                        with lc2:
                            ctr_opts    = [l for l in levels if l != ref_level]
                            default_ctr = [ctr_opts[0]] if ctr_opts else []
                            prev_ctrs   = st.session_state.get("w_contrast_levels", default_ctr)
                            valid_ctrs  = [l for l in prev_ctrs if l in ctr_opts]
                            
                            st.checkbox("Select all comparison groups", key="w_all_contrasts",
                                        on_change=lambda: st.session_state.update({
                                            "w_contrast_levels": ctr_opts if st.session_state.w_all_contrasts else []
                                        }))
                            
                            ctr_levels = st.multiselect(
                                "Comparison group(s)", ctr_opts, 
                                default=valid_ctrs if valid_ctrs else default_ctr,
                                key="w_contrast_levels",
                                help="Numerator group(s) for fold change.",
                            )

                        formula = _formula_str(design_factor, covariate_cols)
                        msg = f"Design formula : {formula}"
                        if ctr_levels:
                            msg += f"\nContrasts      : {design_factor}  ({', '.join(ctr_levels)})  vs  {ref_level}"
                        st.code(msg, language=None)

                        dm_cols    = covariate_cols + [design_factor]
                        dm_preview = metadata_df[
                            [sample_id_col] + [c for c in dm_cols if c in metadata_df.columns]
                        ]
                        st.caption("Design preview (metadata subset)")
                        st.dataframe(dm_preview.head(8), use_container_width=True, hide_index=True)
                    else:
                        st.error(f"'{design_factor}' has fewer than 2 levels: {levels}")

    # ── Section 3: Analysis Settings ──────────────────────────────────────────
    with st.expander("🔬 Analysis Settings", expanded=False):
        runnable = [m for m in list_available_methods() if backend_status[m].available]

        bcols = st.columns(len(list_available_methods()))
        for i, m in enumerate(list_available_methods()):
            with bcols[i]:
                status = backend_status[m]
                icon   = "🟢" if status.available else "🔴"
                st.metric(f"{icon} {m}", "Ready" if status.available else "Unavailable")
                if not status.available:
                    with st.expander("Details"):
                        for d in status.details:
                            st.caption(d)

        st.divider()
        if not runnable:
            st.error("No analysis methods available. Install R and the required packages.")
        else:
            methods = st.multiselect(
                "Methods to run", options=runnable, default=runnable,
                key="w_methods",
                help="Run all three for the Method Comparison tab.",
            )
            sc1, sc2, sc3 = st.columns(3)
            with sc1:
                alpha = st.slider(
                    "FDR threshold (α)", 0.001, 0.2,
                    float(st.session_state.analysis_params["alpha"]), 0.001,
                    key="w_alpha",
                    help="Genes with adjusted p-value ≤ α are called significant.",
                )
            with sc2:
                lfc_cutoff = st.slider(
                    "|log₂FC| threshold", 0.0, 5.0,
                    float(st.session_state.analysis_params["lfc_cutoff"]), 0.1,
                    key="w_lfc_cutoff",
                    help="Minimum absolute fold change. Set to 0 to use p-value only.",
                )
            with sc3:
                top_n = st.slider(
                    "Top N genes (heatmap/table)", 10, 200,
                    int(st.session_state.analysis_params["top_n"]), 10,
                    key="w_top_n",
                    help="Maximum significant genes shown in heatmaps and tables.",
                )

    # ── Run button ─────────────────────────────────────────────────────────────
    st.divider()

    can_run = (
        data_loaded
        and bool(st.session_state.get("w_sample_columns"))
        and bool(st.session_state.get("w_design_factor"))
        and bool(st.session_state.get("w_methods"))
    )

    if not can_run:
        missing = []
        if not data_loaded:
            missing.append("data files")
        if not st.session_state.get("w_sample_columns"):
            missing.append("sample columns")
        if not st.session_state.get("w_design_factor"):
            missing.append("design factor")
        if not st.session_state.get("w_methods"):
            missing.append("analysis methods")
        st.info(f"Complete setup before running: {', '.join(missing)}.")

    btn_col, _ = st.columns([1, 4])
    with btn_col:
        run_clicked = st.button(
            "▶ Run Analysis", type="primary",
            disabled=not can_run,
            use_container_width=True,
        )

    if run_clicked:
        gene_id_col   = st.session_state.w_gene_id_column
        sample_id_col = st.session_state.w_sample_id_column
        sample_cols   = list(st.session_state.w_sample_columns)
        design_factor = st.session_state.w_design_factor
        ref_level     = st.session_state.w_reference_level
        ctr_level     = st.session_state.w_contrast_level
        cov_cols      = list(st.session_state.get("w_covariate_columns", []))
        methods       = list(st.session_state.get("w_methods", []))
        alpha         = float(st.session_state.get("w_alpha", 0.05))
        lfc_cutoff    = float(st.session_state.get("w_lfc_cutoff", 1.0))
        top_n         = int(st.session_state.get("w_top_n", 50))
        params        = {"alpha": alpha, "lfc_cutoff": lfc_cutoff, "top_n": top_n}

        # Validate inputs
        try:
            counts_c, dup_genes = coerce_numeric_counts(
                st.session_state.counts_df.copy(), gene_id_col, sample_cols
            )
            if dup_genes:
                st.warning(
                    f"⚠️ {len(dup_genes)} duplicate gene ID(s) removed (first occurrence kept). "
                    f"Examples: {', '.join(dup_genes[:5])}"
                )
            prepared = prepare_inputs(
                counts_df=counts_c,
                metadata_df=st.session_state.metadata_df.copy(),
                gene_id_column=gene_id_col,
                sample_id_column=sample_id_col,
                sample_columns=sample_cols,
                design_factor=design_factor,
                reference_level=ref_level,
                contrast_level=ctr_level,
                covariate_columns=cov_cols,
            )
        except ValueError as exc:
            st.error(str(exc))
            st.stop()

        annotations = None
        if st.session_state.annotation_df is not None:
            try:
                annotations = parse_annotation_text(st.session_state.annotation_df, gene_id_col)
            except ValueError as exc:
                st.warning(f"Annotation file skipped: {exc}")

        run_progress = st.progress(0.0, text="Starting…")
        new_results: dict = {}

        for i, method in enumerate(methods):
            run_progress.progress(i / len(methods), text=f"Running {method}…")
            try:
                result = run_dge_analysis(
                    method=method,
                    counts_df=prepared.counts_df,
                    metadata_df=prepared.metadata_df,
                    gene_id_column=gene_id_col,
                    sample_id_column=sample_id_col,
                    sample_columns=prepared.sample_columns,
                    design_factor=prepared.design_factor,
                    reference_level=prepared.reference_level,
                    contrast_level=prepared.contrast_level,
                    covariate_columns=prepared.covariate_columns,
                    alpha=alpha,
                    lfc_cutoff=lfc_cutoff,
                    annotation_df=annotations,
                )
                new_results[method] = result
                run_progress.progress((i + 1) / len(methods), text=f"{method} ✓")
            except DGEAnalysisError as exc:
                st.error(f"**{method}** failed: {exc}")

        run_progress.empty()

        if new_results:
            st.session_state.results       = new_results
            st.session_state.analysis_params = params
            _save_run(params, list(new_results.keys()))
            ok = len(new_results)
            st.success(
                f"✅ Analysis complete — {ok}/{len(methods)} method(s) succeeded. "
                "View results in the **📊 Results** or **🔁 Method Comparison** tab."
            )


# ─────────────────────────────────────────────────────────────────────────────
# Tab 2 — Quality Control
# ─────────────────────────────────────────────────────────────────────────────

def _render_qc() -> None:
    if st.session_state.counts_df is None or st.session_state.metadata_df is None:
        st.info("Load data in the **📂 Setup & Run** tab first.")
        return

    gene_id_col   = st.session_state.get("w_gene_id_column")
    sample_id_col = st.session_state.get("w_sample_id_column")
    sample_cols   = st.session_state.get("w_sample_columns", [])
    design_factor = st.session_state.get("w_design_factor", "")
    cov_cols      = st.session_state.get("w_covariate_columns", [])

    if not gene_id_col or not sample_id_col or not sample_cols:
        st.info("Complete **Column Mapping** in the Setup tab first.")
        return

    try:
        qc_df, pca_df, explained, dist_df, log_cpm = _compute_qc_cached(
            st.session_state.counts_df,
            st.session_state.metadata_df,
            gene_id_col,
            sample_id_col,
            tuple(sample_cols),
        )
    except Exception as exc:
        st.error(f"QC computation failed: {exc}")
        return

    st.write(
        "Inspect sample quality **before** running differential expression analysis. "
        "Look for outlier samples, unexpected batch effects, or library size imbalances."
    )

    # ── QC Plot Controls ─────────────────────────────────────────────────────
    st.markdown("---")
    qc_c1, qc_c2, qc_c3 = st.columns([1, 1, 2])
    all_meta_cols = list(st.session_state.metadata_df.columns)
    with qc_c1:
        color_idx = all_meta_cols.index(design_factor) if design_factor in all_meta_cols else 0
        qc_color_by = st.selectbox("Color by", all_meta_cols, index=color_idx, key="qc_color_by")
    with qc_c2:
        shape_options = [None] + [c for c in all_meta_cols if c != sample_id_col]
        qc_shape_by = st.selectbox("Shape by", shape_options, index=0, key="qc_shape_by")
    st.markdown("---")

    c1, c2 = st.columns(2)
    with c1:
        st.plotly_chart(make_library_size_plot(qc_df, sample_id_col, qc_color_by),
                        use_container_width=True)
        with st.expander("💡 Library Size"):
            st.write(
                "Total reads per sample. Large differences (>5×) between samples may "
                "indicate library preparation problems. Each method normalises for this automatically."
            )
    with c2:
        st.plotly_chart(make_detected_genes_plot(qc_df, sample_id_col, qc_color_by),
                        use_container_width=True)
        with st.expander("💡 Detected Genes"):
            st.write(
                "Number of genes with at least 1 read. Samples from the same condition "
                "should have similar detection rates. A marked outlier may be a failed library."
            )

    c3, c4 = st.columns(2)
    with c3:
        st.plotly_chart(make_pca_plot(pca_df, sample_id_col, qc_color_by, explained, qc_shape_by),
                        use_container_width=True)
        with st.expander("💡 PCA"):
            st.write(
                "Principal Component Analysis on log-CPM values. "
                "**Samples from the same condition should cluster together.** "
                "If groups overlap, the differential signal may be weak. "
                "Scattered within-group samples suggest batch effects or outliers."
            )
    with c4:
        st.plotly_chart(make_sample_distance_heatmap(dist_df), use_container_width=True)
        with st.expander("💡 Sample Distance Heatmap"):
            st.write(
                "Euclidean distances between samples in log-CPM space. "
                "Darker = more similar. Samples in the same condition should have "
                "small mutual distances. A clear block structure validates condition-driven expression."
            )

    c5, c6 = st.columns(2)
    with c5:
        st.plotly_chart(make_normalization_boxplot(log_cpm, "Expression Distribution (log₂ CPM)"),
                        use_container_width=True)
        with st.expander("💡 Expression Distribution"):
            st.write(
                "log₂ CPM distributions per sample. Medians and spreads should be "
                "roughly similar across samples. A drastically different distribution "
                "may indicate a problematic sample."
            )
    with c6:
        st.caption("QC metrics table")
        st.dataframe(qc_df, use_container_width=True, hide_index=True)


# ─────────────────────────────────────────────────────────────────────────────
# Tab 3 — Results (one sub-tab per method)
# ─────────────────────────────────────────────────────────────────────────────

def _render_results() -> None:
    if not st.session_state.results:
        st.info(
            "No results yet. Configure your analysis in **📂 Setup & Run** and click **▶ Run Analysis**. "
            "You can run all three methods to compare them in the **🔁 Method Comparison** tab."
        )
        return

    # Live threshold sliders — results re-computed immediately, no re-run needed
    st.caption("Adjust thresholds to update results live (no re-running required).")
    params  = st.session_state.analysis_params
    tc1, tc2, tc3 = st.columns(3)
    with tc1:
        alpha = st.slider("FDR threshold (α)", 0.001, 0.2,
                          float(st.session_state.get("res_alpha", params["alpha"])),
                          0.001, key="res_alpha")
    with tc2:
        lfc = st.slider("|log₂FC| threshold", 0.0, 5.0,
                        float(st.session_state.get("res_lfc", params["lfc_cutoff"])),
                        0.1, key="res_lfc")
    with tc3:
        top_n = st.slider("Top N genes", 10, 200,
                          int(st.session_state.get("res_top_n", params["top_n"])),
                          10, key="res_top_n")

    gene_id_col   = st.session_state.get("w_gene_id_column", "gene_id")
    sample_id_col = st.session_state.get("w_sample_id_column", "sample")
    design_factor = st.session_state.get("w_design_factor", "condition")

    method_names = list(st.session_state.results.keys())
    method_tabs  = st.tabs(method_names)

    for tab, method in zip(method_tabs, method_names):
        with tab:
            result = st.session_state.results[method]
            res_df = apply_thresholds(result.results, alpha, lfc)
            s      = _live_summary(res_df)

            # Metric row
            m1, m2, m3, m4, m5 = st.columns(5)
            m1.metric("Genes tested",    f"{s['genes_tested']:,}")
            m2.metric("✅ Significant",  f"{s['significant_genes']:,}")
            m3.metric("⬆ Upregulated",  f"{s['upregulated']:,}")
            m4.metric("⬇ Downregulated", f"{s['downregulated']:,}")
            m5.metric("🚫 Filtered",     f"{s['filtered']:,}",
                      help="Genes with NA adjusted p-value (DESeq2 independent filtering).")

            # Interpretation card
            n_sig, n_up, n_down = s["significant_genes"], s["upregulated"], s["downregulated"]
            n_test = s["genes_tested"]
            pct    = f"{100 * n_sig / n_test:.1f}%" if n_test > 0 else "0%"

            if n_sig == 0:
                interp = (
                    f"No genes pass the current criteria (FDR ≤ {alpha}, |log₂FC| ≥ {lfc}). "
                    "Try relaxing the thresholds above, or check the QC tab for data quality issues."
                )
            else:
                up_row   = res_df[res_df["direction"] == "Upregulated"]
                down_row = res_df[res_df["direction"] == "Downregulated"]
                top_up   = (
                    f"Top upregulated: **{up_row.iloc[0][gene_id_col]}** "
                    f"(log₂FC = {up_row.iloc[0]['log2_fold_change']:.2f}, "
                    f"FDR = {up_row.iloc[0]['adjusted_p_value']:.2e}). "
                    if not up_row.empty else ""
                )
                top_down = (
                    f"Top downregulated: **{down_row.iloc[0][gene_id_col]}** "
                    f"(log₂FC = {down_row.iloc[0]['log2_fold_change']:.2f}, "
                    f"FDR = {down_row.iloc[0]['adjusted_p_value']:.2e})."
                    if not down_row.empty else ""
                )
                interp = (
                    f"**{method}** identified **{n_sig} significant genes** ({pct} of {n_test:,} tested): "
                    f"**{n_up} upregulated**, **{n_down} downregulated** "
                    f"(FDR ≤ {alpha}, |log₂FC| ≥ {lfc}). {top_up}{top_down}"
                )

            with st.expander("💡 Result Interpretation", expanded=True):
                st.write(interp)

            st.divider()

            # Plots
            p1, p2 = st.columns(2)
            with p1:
                st.plotly_chart(
                    make_volcano_plot(res_df, alpha=alpha, lfc_cutoff=lfc, gene_column=gene_id_col),
                    use_container_width=True,
                )
                with st.expander("💡 Volcano plot"):
                    st.write(
                        "X-axis: log₂ fold change (positive = higher in comparison group). "
                        "Y-axis: −log₁₀(adjusted p-value). Top-right/top-left = most significant with "
                        "largest effect. Dashed lines are your current thresholds. "
                        "The 10 most significant genes are labelled."
                    )
            with p2:
                st.plotly_chart(
                    make_ma_plot(res_df, gene_column=gene_id_col),
                    use_container_width=True,
                )
                with st.expander("💡 MA plot"):
                    st.write(
                        "X-axis: mean expression (log scale). Y-axis: log₂ fold change. "
                        "Fold changes of genes at the far left (low count) are less reliable. "
                        "A good pattern shows most genes near log₂FC = 0 (the horizontal dashed line)."
                    )

            p3, p4 = st.columns(2)
            with p3:
                st.plotly_chart(make_pvalue_histogram(res_df), use_container_width=True)
                with st.expander("💡 P-value distribution"):
                    st.write(
                        "A peak near 0 with a roughly flat right tail is ideal — it indicates real "
                        "differential signal. A peak near 1, or a U-shape, suggests the method "
                        "or the design might be misspecified."
                    )
            with p4:
                st.plotly_chart(
                    make_top_gene_heatmap(result.normalized_counts, res_df,
                                          gene_column=gene_id_col, top_n=top_n),
                    use_container_width=True,
                )
                with st.expander("💡 Top DE Genes Heatmap"):
                    st.write(
                        f"Only statistically significant genes are included (up to {top_n}). "
                        "Expression values are row z-scored for comparability across genes. "
                        "Rows and columns are hierarchically clustered. "
                        "A clean block pattern confirms consistent differential expression."
                    )

            # Gene-level explorer
            st.markdown("---")
            st.markdown("**Gene-level exploration**")
            ge1, ge2 = st.columns([1, 2])

            with ge1:
                # Searches FULL result set, not top_n slice — fixes Flaw 8
                search = st.text_input(
                    "Search all genes", placeholder="Type gene ID…",
                    key=f"search_{method}",
                    help="Searches all tested genes, not just the top-N table.",
                )
                disp_cols = [c for c in [gene_id_col, "base_mean", "log2_fold_change",
                                          "adjusted_p_value", "direction", "is_significant"]
                             if c in res_df.columns]
                if search:
                    mask = res_df[gene_id_col].astype(str).str.contains(search, case=False, na=False)
                    show = res_df[mask][disp_cols]
                    st.caption(f"{len(show)} matching genes")
                    st.dataframe(show.head(top_n), use_container_width=True, hide_index=True)
                else:
                    sig_show = res_df[res_df["is_significant"]][disp_cols]
                    st.caption(f"Top {min(top_n, len(sig_show))} significant genes")
                    st.dataframe(sig_show.head(top_n), use_container_width=True, hide_index=True)

            with ge2:
                all_genes = res_df[gene_id_col].astype(str).tolist()
                if all_genes and st.session_state.metadata_df is not None:
                    sel_gene = st.selectbox(
                        "Expression boxplot for gene",
                        options=all_genes[:5000],
                        key=f"genebox_{method}",
                    )
                    st.plotly_chart(
                        make_gene_expression_plot(
                            result.normalized_counts,
                            st.session_state.metadata_df,
                            sample_id_column=sample_id_col,
                            design_factor=design_factor,
                            gene_id=sel_gene,
                        ),
                        use_container_width=True,
                    )
                    with st.expander("💡 Expression boxplot"):
                        st.write(
                            "Normalized expression per sample, split by the design factor. "
                            "Look for a clear group separation. Large within-group variance "
                            "may reduce power to detect differential expression."
                        )

            # Enrichment
            st.markdown("---")
            st.markdown("**Functional Enrichment (ORA)**")
            st.caption("Upload a GMT gene set file to run over-representation analysis on significant genes.")
            gmt_file = st.file_uploader("Upload GMT file", type=["gmt"], key=f"gmt_{method}")
            if gmt_file is not None:
                gene_sets = _parse_gmt_text(gmt_file.getvalue().decode("utf-8", errors="ignore"))
                sig_genes = set(res_df[res_df["is_significant"]][gene_id_col].astype(str))
                background = set(res_df[gene_id_col].astype(str))
                enrich_df  = _run_enrichment(sig_genes, background, gene_sets)
                if enrich_df.empty:
                    st.info("No enriched terms found at the current threshold.")
                else:
                    st.dataframe(enrich_df.head(25), use_container_width=True, hide_index=True)
                    _download_csv(enrich_df, "Download enrichment table",
                                  f"{method.lower()}_enrichment.csv")

            # Downloads
            st.markdown("---")
            dc1, dc2, dc3 = st.columns(3)
            with dc1:
                _download_csv(res_df, "All results (CSV)", f"{method.lower()}_full.csv")
            with dc2:
                sig_only = res_df[res_df["is_significant"]]
                if not sig_only.empty:
                    _download_csv(sig_only, "Significant genes (CSV)",
                                  f"{method.lower()}_significant.csv")
            with dc3:
                _download_csv(result.normalized_counts.reset_index(),
                              "Normalized counts (CSV)",
                              f"{method.lower()}_normalized.csv")


# ─────────────────────────────────────────────────────────────────────────────
# Tab 4 — Method Comparison
# ─────────────────────────────────────────────────────────────────────────────

def _render_comparison() -> None:
    if not st.session_state.results:
        st.info("Run at least one analysis method in **📂 Setup & Run** first.")
        return

    results      = st.session_state.results
    method_names = list(results.keys())
    gene_id_col  = st.session_state.get("w_gene_id_column", "gene_id")

    # Use live thresholds from Results tab (or fall back to stored params)
    alpha = float(st.session_state.get("res_alpha",
                  st.session_state.analysis_params.get("alpha", 0.05)))
    lfc   = float(st.session_state.get("res_lfc",
                  st.session_state.analysis_params.get("lfc_cutoff", 1.0)))

    thresholded = {m: apply_thresholds(r.results, alpha, lfc) for m, r in results.items()}
    sig_sets    = {
        m: set(df[df["is_significant"]][gene_id_col].astype(str))
        for m, df in thresholded.items()
    }

    if len(method_names) == 1:
        st.warning(
            f"Only **{method_names[0]}** results are available. "
            "Go to **📂 Setup & Run**, add more methods, and re-run to compare results."
        )
        return

    # ── Overview cards ─────────────────────────────────────────────────────────
    st.markdown("### Method Overview")
    metric_cols = st.columns(len(method_names))
    for col, method in zip(metric_cols, method_names):
        df    = thresholded[method]
        n_sig = int(df["is_significant"].sum())
        n_up  = int((df["direction"] == "Upregulated").sum())
        n_dn  = int((df["direction"] == "Downregulated").sum())
        with col:
            st.metric(method, f"{n_sig:,} sig. genes")
            st.caption(f"⬆ {n_up:,} up  |  ⬇ {n_dn:,} down")

    # Consensus summary
    all_agree  = set.intersection(*sig_sets.values())
    all_union  = set.union(*sig_sets.values())
    any_two    = set()
    for m1, m2 in combinations(method_names, 2):
        any_two |= sig_sets[m1] & sig_sets[m2]
    unique_to_one = len(all_union) - len(any_two)

    cc1, cc2, cc3 = st.columns(3)
    cc1.metric(f"All {len(method_names)} methods agree", f"{len(all_agree):,} genes")
    cc2.metric("At least 2 methods agree",                f"{len(any_two):,} genes")
    cc3.metric("Unique to 1 method",                      f"{unique_to_one:,} genes")

    if all_union:
        pct = 100 * len(all_agree) / len(all_union)
        with st.expander("💡 What does this mean?", expanded=True):
            st.write(
                f"**{len(all_agree):,} genes ({pct:.0f}% of all significant genes) are identified by "
                f"all {len(method_names)} methods.** High consensus gives greater confidence — these "
                "genes are robust to the choice of statistical method. "
                f"**{unique_to_one:,} genes are unique to a single method** — inspect these to understand "
                "where methods diverge (often due to differences in normalisation, dispersion estimation, "
                "or how lowly expressed genes are handled)."
            )

    st.divider()

    # ── Gene × Method table ────────────────────────────────────────────────────
    st.markdown("### Gene × Method Agreement")
    st.write(
        "Every gene significant in at least one method is listed. "
        "✅ = method called it significant; ❌ = not significant."
    )

    comp_df = make_method_comparison_table(thresholded, gene_id_col)

    if comp_df.empty:
        st.info("No significant genes found by any method. Try relaxing the thresholds in the Results tab.")
    else:
        gene_filter = st.text_input(
            "Filter by gene ID", placeholder="Type to search…", key="comp_search"
        )
        show_df = comp_df
        if gene_filter:
            mask    = comp_df[gene_id_col].astype(str).str.contains(gene_filter, case=False, na=False)
            show_df = comp_df[mask]

        pretty_df = show_df.copy()
        for m in method_names:
            if m in pretty_df.columns:
                pretty_df[m] = pretty_df[m].map({True: "✅", False: "❌"})

        st.dataframe(pretty_df, use_container_width=True, hide_index=True)
        st.caption(f"{len(show_df):,} of {len(comp_df):,} genes shown")
        _download_csv(comp_df, "Download comparison table (CSV)", "method_comparison.csv")

    st.divider()

    # ── Overlap bar ────────────────────────────────────────────────────────────
    st.markdown("### Gene Set Overlap (Exclusive Counts)")
    overlap_counts: dict[str, int] = {}
    for r in range(1, len(method_names) + 1):
        for group in combinations(method_names, r):
            group_set = sig_sets[group[0]].copy()
            for m in group[1:]:
                group_set &= sig_sets[m]
            others = set(method_names) - set(group)
            exclusive = group_set.copy()
            for om in others:
                exclusive -= sig_sets[om]
            overlap_counts[" & ".join(group)] = len(exclusive)

    st.plotly_chart(make_overlap_bar(overlap_counts), use_container_width=True)
    with st.expander("💡 Reading the overlap chart"):
        st.write(
            "Each bar shows genes **exclusive** to that set or combination. "
            "'DESeq2 & edgeR' = genes found by both DESeq2 and edgeR but **not** limma-voom. "
            "A tall 'all methods' bar indicates strong consensus."
        )

    st.divider()

    # ── Pairwise LFC scatter ───────────────────────────────────────────────────
    st.markdown("### Pairwise log₂FC Comparison")
    st.write(
        "Compare fold-change estimates between methods. Points on the diagonal line indicate "
        "agreement on effect size. Colour shows which methods called the gene significant."
    )

    sc1, sc2 = st.columns(2)
    with sc1:
        method_a = st.selectbox("Method A (x-axis)", method_names, index=0, key="lfc_a")
    with sc2:
        mb_opts  = [m for m in method_names if m != method_a]
        method_b = st.selectbox("Method B (y-axis)", mb_opts, index=0, key="lfc_b")

    st.plotly_chart(
        make_lfc_scatter(thresholded, gene_id_col, method_a, method_b),
        use_container_width=True,
    )
    with st.expander("💡 Interpreting the scatter"):
        st.write(
            "**Spearman ρ** measures rank agreement of fold-change estimates. "
            "ρ ≥ 0.95 = very high concordance. "
            "**Green** (Both significant) = highest confidence. "
            "**Red/Blue** (Only one method) = borderline genes or dispersion estimation differences — "
            "worth a closer look."
        )

    if len(method_names) >= 3:
        st.divider()
        st.markdown("### log₂FC Correlation Across All Methods")
        st.plotly_chart(
            make_log2fc_correlation_plot(thresholded, gene_id_col),
            use_container_width=True,
        )
        with st.expander("💡 Correlation heatmap"):
            st.write(
                "Spearman correlations of log₂FC across all pairs of methods. "
                "Values close to 1.0 indicate high agreement. "
                "All three methods typically exceed ρ = 0.95 on clean data. "
                "Lower values may indicate differences in handling of lowly expressed genes."
            )

    st.divider()

    # ── Download & Report ──────────────────────────────────────────────────────
    st.markdown("### Download & Report")
    dl1, dl2, dl3 = st.columns(3)
    with dl1:
        if not comp_df.empty:
            _download_csv(comp_df, "Comparison table (CSV)", "method_comparison.csv")
    with dl2:
        payload = {
            "gene_id_column":    st.session_state.get("w_gene_id_column"),
            "design_factor":     st.session_state.get("w_design_factor"),
            "reference_level":   st.session_state.get("w_reference_level"),
            "contrast_level":    st.session_state.get("w_contrast_level"),
            "covariate_columns": st.session_state.get("w_covariate_columns", []),
            "analysis_params":   st.session_state.analysis_params,
            "run_log":           st.session_state.run_log,
        }
        st.download_button(
            "Session config (JSON)",
            data=json.dumps(payload, indent=2).encode("utf-8"),
            file_name="dge_session.json",
            mime="application/json",
        )
    with dl3:
        st.download_button(
            "PDF report",
            data=build_pdf_report(_summary_lines()),
            file_name="dge_report.pdf",
            mime="application/pdf",
        )


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

_init_state()

backend_status = _get_backend_status_cached()

st.title("🧬 DGE Studio")
st.caption(
    "Rapid differential gene expression analysis with DESeq2, edgeR, and limma-voom. "
    "Run all three methods and compare their results gene by gene."
)

tab_setup, tab_qc, tab_results, tab_compare = st.tabs([
    "📂 Setup & Run",
    "🔬 Quality Control",
    "📊 Results",
    "🔁 Method Comparison",
])

with tab_setup:
    _render_setup(backend_status)

with tab_qc:
    _render_qc()

with tab_results:
    _render_results()

with tab_compare:
    _render_comparison()