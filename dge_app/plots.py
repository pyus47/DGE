import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.stats import spearmanr

def make_library_size_plot(qc_df: pd.DataFrame, sample_id_col: str, qc_color_by: str) -> go.Figure:
    fig = px.bar(
        qc_df, x=sample_id_col, y="library_size", color=qc_color_by,
        title="Library Size per Sample"
    )
    fig.update_layout(xaxis_title="Sample", yaxis_title="Total Counts (Library Size)", margin=dict(t=50, l=0, r=0, b=0))
    return fig

def make_detected_genes_plot(qc_df: pd.DataFrame, sample_id_col: str, qc_color_by: str) -> go.Figure:
    fig = px.bar(
        qc_df, x=sample_id_col, y="detected_genes", color=qc_color_by,
        title="Detected Genes (>0 reads)"
    )
    fig.update_layout(xaxis_title="Sample", yaxis_title="Number of Genes", margin=dict(t=50, l=0, r=0, b=0))
    return fig

def make_pca_plot(pca_df: pd.DataFrame, sample_id_col: str, qc_color_by: str, explained: list, qc_shape_by: str | None) -> go.Figure:
    fig = px.scatter(
        pca_df, x="PC1", y="PC2", color=qc_color_by, symbol=qc_shape_by, hover_data=[sample_id_col],
        title="PCA (log₂ CPM)", size_max=12
    )
    fig.update_traces(marker=dict(size=12, line=dict(width=1, color="DarkSlateGrey")))
    fig.update_layout(
        xaxis_title=f"PC1 ({explained[0]:.1f}%)",
        yaxis_title=f"PC2 ({explained[1]:.1f}%)",
        margin=dict(t=50, l=0, r=0, b=0)
    )
    return fig

def make_sample_distance_heatmap(dist_df: pd.DataFrame) -> go.Figure:
    fig = px.imshow(
        dist_df, text_auto=".1f", color_continuous_scale="Plasma",
        title="Sample Distances"
    )
    fig.update_layout(xaxis_title="Sample", yaxis_title="Sample", margin=dict(t=50, l=0, r=0, b=0))
    return fig

def make_normalization_boxplot(log_cpm: pd.DataFrame, title: str = "Expression Distribution (log₂ CPM)") -> go.Figure:
    fig = px.box(
        log_cpm, title=title, points="outliers"
    )
    fig.update_layout(xaxis_title="Sample", yaxis_title="log₂ CPM", margin=dict(t=50, l=0, r=0, b=0), showlegend=False)
    return fig

def make_volcano_plot(res_df: pd.DataFrame, alpha: float = 0.05, lfc_cutoff: float = 1.0, gene_column: str = "gene_id") -> go.Figure:
    df = res_df.copy()
    df["-log10_padj"] = -np.log10(df["adjusted_p_value"].fillna(1.0).replace(0, 1e-300))
    
    color_map = {
        "Upregulated": "#EF553B",
        "Downregulated": "#00CC96",
        "Not significant": "#555555",
        "Filtered": "#333333"
    }
    
    fig = px.scatter(
        df, x="log2_fold_change", y="-log10_padj", color="direction",
        hover_data=[gene_column, "p_value", "adjusted_p_value", "base_mean"],
        color_discrete_map=color_map,
        title="Volcano Plot"
    )
    
    fig.add_hline(y=-np.log10(alpha), line_dash="dash", line_color="gray", annotation_text=f"FDR={alpha}")
    fig.add_vline(x=lfc_cutoff, line_dash="dash", line_color="gray")
    fig.add_vline(x=-lfc_cutoff, line_dash="dash", line_color="gray")
    
    top_10 = df[df["is_significant"]].sort_values("adjusted_p_value").head(10)
    for i, row in top_10.iterrows():
        fig.add_annotation(
            x=row["log2_fold_change"], y=row["-log10_padj"],
            text=row[gene_column], showarrow=True, arrowhead=1, ax=20, ay=-20,
            font=dict(size=10)
        )
        
    fig.update_layout(xaxis_title="log₂ Fold Change", yaxis_title="−log₁₀(adjusted p-value)", margin=dict(t=50, l=0, r=0, b=0))
    return fig

def make_ma_plot(res_df: pd.DataFrame, gene_column: str = "gene_id") -> go.Figure:
    df = res_df.copy()
    color_map = {
        "Upregulated": "#EF553B",
        "Downregulated": "#00CC96",
        "Not significant": "#555555",
        "Filtered": "#333333"
    }
    fig = px.scatter(
        df, x="base_mean", y="log2_fold_change", color="direction",
        hover_data=[gene_column, "adjusted_p_value"],
        log_x=True, color_discrete_map=color_map,
        title="MA Plot"
    )
    fig.add_hline(y=0, line_dash="dash", line_color="white")
    fig.update_layout(xaxis_title="Base Mean Expression", yaxis_title="log₂ Fold Change", margin=dict(t=50, l=0, r=0, b=0))
    return fig

def make_pvalue_histogram(res_df: pd.DataFrame) -> go.Figure:
    df = res_df.dropna(subset=["p_value"])
    fig = px.histogram(
        df, x="p_value", nbins=50, title="P-value Distribution"
    )
    fig.update_layout(xaxis_title="p-value", yaxis_title="Frequency", margin=dict(t=50, l=0, r=0, b=0))
    return fig

def make_top_gene_heatmap(normalized_counts: pd.DataFrame, res_df: pd.DataFrame, gene_column: str = "gene_id", top_n: int = 50) -> go.Figure:
    sig_genes = res_df[res_df["is_significant"]].sort_values("adjusted_p_value").head(top_n)[gene_column]
    if sig_genes.empty:
        fig = go.Figure()
        fig.update_layout(title="No significant genes to display")
        return fig
        
    plot_data = normalized_counts.loc[sig_genes].astype(float)
    z_scored = plot_data.apply(lambda row: (row - row.mean()) / (row.std() or 1), axis=1)
    
    fig = px.imshow(
        z_scored, color_continuous_scale="RdBu_r", zmin=-3, zmax=3,
        title=f"Top {len(sig_genes)} Significant Genes (Row Z-score)", aspect="auto"
    )
    fig.update_layout(yaxis_title="Gene", xaxis_title="Sample", margin=dict(t=50, l=0, r=0, b=0))
    return fig

def make_gene_expression_plot(normalized_counts: pd.DataFrame, metadata_df: pd.DataFrame, sample_id_column: str, design_factor: str, gene_id: str) -> go.Figure:
    if gene_id not in normalized_counts.index:
        fig = go.Figure()
        fig.update_layout(title=f"Gene {gene_id} not found")
        return fig
        
    expr = normalized_counts.loc[gene_id].to_frame(name="normalized_expression")
    expr.index.name = sample_id_column
    df = expr.reset_index().merge(metadata_df[[sample_id_column, design_factor]], on=sample_id_column, how="left")
    
    fig = px.box(
        df, x=design_factor, y="normalized_expression", points="all",
        hover_data=[sample_id_column], color=design_factor,
        title=f"Expression: {gene_id}"
    )
    fig.update_layout(xaxis_title=design_factor.capitalize(), yaxis_title="Normalized Expression", margin=dict(t=50, l=0, r=0, b=0))
    return fig

def make_method_comparison_table(thresholded: dict, gene_id_col: str) -> pd.DataFrame:
    dfs = []
    for method, df in thresholded.items():
        sig_df = df[df["is_significant"]][[gene_id_col]].copy()
        sig_df[method] = True
        dfs.append(sig_df)
        
    if not dfs:
        return pd.DataFrame()
        
    all_sig_genes = pd.concat([d[[gene_id_col]] for d in dfs]).drop_duplicates()
    merged = all_sig_genes
    for df in dfs:
        merged = merged.merge(df, on=gene_id_col, how="left")
        
    for method in thresholded.keys():
        if method in merged.columns:
            merged[method] = merged[method].fillna(False)
            
    return merged

def make_overlap_bar(overlap_counts: dict) -> go.Figure:
    df = pd.DataFrame(list(overlap_counts.items()), columns=["Combination", "Exclusive Count"])
    fig = px.bar(
        df, x="Combination", y="Exclusive Count", text="Exclusive Count",
        title="Gene Set Overlap", color="Combination"
    )
    fig.update_traces(textposition='outside')
    fig.update_layout(xaxis_title="Method Combination", yaxis_title="Number of Exclusive Genes", showlegend=False, margin=dict(t=50, l=0, r=0, b=0))
    return fig

def make_lfc_scatter(thresholded: dict, gene_id_col: str, method_a: str, method_b: str) -> go.Figure:
    df_a = thresholded[method_a][[gene_id_col, "log2_fold_change", "is_significant"]].copy()
    df_b = thresholded[method_b][[gene_id_col, "log2_fold_change", "is_significant"]].copy()
    
    merged = df_a.merge(df_b, on=gene_id_col, suffixes=(f"_{method_a}", f"_{method_b}"))
    if merged.empty:
        fig = go.Figure()
        fig.update_layout(title="No overlapping genes found")
        return fig
        
    def get_category(row):
        a_sig = row[f"is_significant_{method_a}"]
        b_sig = row[f"is_significant_{method_b}"]
        if a_sig and b_sig:
            return "Both Significant"
        elif a_sig:
            return f"Only {method_a}"
        elif b_sig:
            return f"Only {method_b}"
        else:
            return "Neither"
            
    merged["Significance"] = merged.apply(get_category, axis=1)
    
    mask = merged[f"log2_fold_change_{method_a}"].notna() & merged[f"log2_fold_change_{method_b}"].notna()
    if mask.sum() > 1:
        corr, _ = spearmanr(merged.loc[mask, f"log2_fold_change_{method_a}"], merged.loc[mask, f"log2_fold_change_{method_b}"])
    else:
        corr = 0.0
    
    color_map = {
        "Both Significant": "#00CC96",
        f"Only {method_a}": "#EF553B",
        f"Only {method_b}": "#AB63FA",
        "Neither": "#444444"
    }
    
    fig = px.scatter(
        merged, x=f"log2_fold_change_{method_a}", y=f"log2_fold_change_{method_b}",
        color="Significance", hover_data=[gene_id_col], color_discrete_map=color_map,
        title=f"log₂FC Comparison (Spearman ρ = {corr:.3f})"
    )
    
    min_val = merged[[f"log2_fold_change_{method_a}", f"log2_fold_change_{method_b}"]].min().min()
    max_val = merged[[f"log2_fold_change_{method_a}", f"log2_fold_change_{method_b}"]].max().max()
    if pd.notna(min_val) and pd.notna(max_val):
        fig.add_shape(
            type="line", x0=min_val, y0=min_val, x1=max_val, y1=max_val,
            line=dict(color="white", dash="dash")
        )
    
    fig.update_layout(xaxis_title=f"{method_a} log₂FC", yaxis_title=f"{method_b} log₂FC", margin=dict(t=50, l=0, r=0, b=0))
    return fig

def make_log2fc_correlation_plot(thresholded: dict, gene_id_col: str) -> go.Figure:
    methods = list(thresholded.keys())
    merged = None
    for m in methods:
        lfc_df = thresholded[m][[gene_id_col, "log2_fold_change"]].rename(columns={"log2_fold_change": m})
        if merged is None:
            merged = lfc_df
        else:
            merged = merged.merge(lfc_df, on=gene_id_col, how="inner")
            
    if merged is None or merged.empty:
        fig = go.Figure()
        fig.update_layout(title="No overlapping genes")
        return fig
        
    corr_matrix = merged[methods].corr(method="spearman")
    
    fig = px.imshow(
        corr_matrix, text_auto=".3f", color_continuous_scale="Viridis", zmin=0.5, zmax=1.0,
        title="Spearman Correlation (log₂FC)"
    )
    fig.update_layout(margin=dict(t=50, l=0, r=0, b=0))
    return fig
