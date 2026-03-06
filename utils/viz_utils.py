"""
BioSignal Discovery Engine
Utils: Visualization Functions
================================
Funciones para generación de figuras del pipeline.
Produce: volcano plots, heatmaps, dot plots de pathways y figuras de meta-análisis.
Todas las figuras se guardan como PNG de alta resolución.
"""

from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
from loguru import logger

try:
    import matplotlib
    matplotlib.use("Agg")  # Backend sin display (modo headless/server)
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.colors import Normalize, LinearSegmentedColormap
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    logger.warning("matplotlib no disponible — figuras deshabilitadas")

try:
    import seaborn as sns
    SEABORN_AVAILABLE = True
except ImportError:
    SEABORN_AVAILABLE = False


# Paleta corporativa BioSignal
COLORS = {
    "up": "#d62728",       # Rojo: sobreexpresado
    "down": "#1f77b4",     # Azul: subexpresado
    "ns": "#aaaaaa",       # Gris: no significativo
    "primary": "#1a3a5c",
    "secondary": "#4a7fb5",
    "accent": "#e07b39",
}

DPI = 150
FIGURE_FORMAT = "png"


def _check_matplotlib():
    if not MATPLOTLIB_AVAILABLE:
        logger.error("matplotlib no disponible. Instalar con: pip install matplotlib seaborn")
        return False
    return True


# ------------------------------------------------------------------
# Volcano Plot
# ------------------------------------------------------------------

def plot_volcano(dea_df: pd.DataFrame, output_path: Union[str, Path],
                 gene_col: str = "gene", lfc_col: str = "log2FC",
                 padj_col: str = "padj", lfc_threshold: float = 1.0,
                 padj_threshold: float = 0.05, title: str = "Volcano Plot",
                 label_top_n: int = 15) -> Optional[str]:
    """
    Genera un volcano plot de los resultados DEA.

    Args:
        dea_df:         DataFrame con resultados DEA
        output_path:    Ruta de salida del PNG
        lfc_threshold:  Umbral log2FC para significancia
        padj_threshold: Umbral FDR para significancia
        label_top_n:    Número de genes a etiquetar

    Returns:
        Ruta del archivo generado o None si falla
    """
    if not _check_matplotlib():
        return None

    try:
        df = dea_df.copy().dropna(subset=[lfc_col, padj_col])
        df["-log10padj"] = -np.log10(df[padj_col].clip(lower=1e-300))

        # Clasificar genes
        df["category"] = "NS"
        df.loc[(df[lfc_col] >= lfc_threshold) & (df[padj_col] < padj_threshold), "category"] = "UP"
        df.loc[(df[lfc_col] <= -lfc_threshold) & (df[padj_col] < padj_threshold), "category"] = "DOWN"

        n_up = (df["category"] == "UP").sum()
        n_down = (df["category"] == "DOWN").sum()

        fig, ax = plt.subplots(figsize=(9, 7))

        color_map = {"UP": COLORS["up"], "DOWN": COLORS["down"], "NS": COLORS["ns"]}
        size_map = {"UP": 18, "DOWN": 18, "NS": 5}
        alpha_map = {"UP": 0.85, "DOWN": 0.85, "NS": 0.25}
        zorder_map = {"UP": 3, "DOWN": 3, "NS": 1}

        for cat in ["NS", "DOWN", "UP"]:
            mask = df["category"] == cat
            ax.scatter(
                df.loc[mask, lfc_col],
                df.loc[mask, "-log10padj"],
                c=color_map[cat], s=size_map[cat],
                alpha=alpha_map[cat], linewidths=0,
                zorder=zorder_map[cat],
                label=f"{cat} (n={mask.sum()})"
            )

        # Líneas de umbral
        ax.axhline(-np.log10(padj_threshold), color="black", linestyle="--", linewidth=0.8, alpha=0.5)
        ax.axvline(lfc_threshold, color="black", linestyle="--", linewidth=0.8, alpha=0.5)
        ax.axvline(-lfc_threshold, color="black", linestyle="--", linewidth=0.8, alpha=0.5)

        # Etiquetar top genes
        if gene_col in df.columns and label_top_n > 0:
            sig = df[df["category"] != "NS"].copy()
            sig["score"] = sig["-log10padj"] * sig[lfc_col].abs()
            top = sig.nlargest(label_top_n, "score")

            for _, row in top.iterrows():
                ax.annotate(
                    row[gene_col],
                    xy=(row[lfc_col], row["-log10padj"]),
                    fontsize=6.5, ha="center", va="bottom",
                    xytext=(0, 3), textcoords="offset points",
                    color=color_map[row["category"]],
                )

        ax.set_xlabel("log₂ Fold Change", fontsize=11)
        ax.set_ylabel("-log₁₀(adjusted p-value)", fontsize=11)
        ax.set_title(title, fontsize=13, fontweight="bold", color=COLORS["primary"])
        ax.legend(loc="upper right", fontsize=8, framealpha=0.7)
        ax.text(0.02, 0.97, f"↑ UP: {n_up} | ↓ DOWN: {n_down}", transform=ax.transAxes,
                fontsize=8, va="top", color="grey")

        plt.tight_layout()
        plt.savefig(str(output_path), dpi=DPI, format=FIGURE_FORMAT, bbox_inches="tight")
        plt.close(fig)

        logger.debug(f"[viz_utils] Volcano plot guardado: {output_path}")
        return str(output_path)

    except Exception as e:
        logger.error(f"[viz_utils] Error generando volcano plot: {e}")
        plt.close("all")
        return None


# ------------------------------------------------------------------
# Heatmap de genes top
# ------------------------------------------------------------------

def plot_heatmap(expression_matrix: pd.DataFrame, output_path: Union[str, Path],
                 genes: list = None, sample_meta: pd.DataFrame = None,
                 title: str = "Expression Heatmap", top_n: int = 50) -> Optional[str]:
    """
    Genera un heatmap de expresión de los genes más variables.

    Args:
        expression_matrix: DataFrame genes×muestras (normalizado)
        output_path:       Ruta de salida PNG
        genes:             Subset de genes a mostrar (si None, usa top_n más variables)
        sample_meta:       DataFrame con metadata de muestras (columna 'group')
        top_n:             N genes más variables si genes=None

    Returns:
        Ruta del archivo generado o None
    """
    if not _check_matplotlib() or not SEABORN_AVAILABLE:
        logger.warning("[viz_utils] seaborn requerido para heatmap")
        return None

    try:
        mat = expression_matrix.copy()

        # Seleccionar genes
        if genes is not None:
            common = [g for g in genes if g in mat.index]
            mat = mat.loc[common]
        else:
            variances = mat.var(axis=1)
            mat = mat.loc[variances.nlargest(top_n).index]

        if mat.empty:
            logger.warning("[viz_utils] Matriz vacía para heatmap")
            return None

        # Z-score por gen
        mat_z = mat.subtract(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, 1), axis=0)

        # Colores de columna por grupo
        col_colors = None
        if sample_meta is not None and "group" in sample_meta.columns:
            palette = {"case": COLORS["up"], "control": COLORS["down"]}
            col_colors = sample_meta["group"].map(
                lambda g: COLORS["up"] if "case" in g.lower() or "tumor" in g.lower() or "disease" in g.lower()
                else COLORS["down"]
            )

        fig_height = max(8, min(len(mat) * 0.25, 30))
        cmap = LinearSegmentedColormap.from_list("biosignal", [COLORS["down"], "white", COLORS["up"]])

        g = sns.clustermap(
            mat_z,
            cmap=cmap,
            col_colors=col_colors,
            figsize=(12, fig_height),
            xticklabels=False,
            yticklabels=True,
            dendrogram_ratio=(0.1, 0.15),
            vmin=-2.5, vmax=2.5,
            linewidths=0,
        )
        g.fig.suptitle(title, fontsize=12, fontweight="bold",
                        color=COLORS["primary"], y=1.01)

        plt.savefig(str(output_path), dpi=DPI, format=FIGURE_FORMAT, bbox_inches="tight")
        plt.close("all")

        logger.debug(f"[viz_utils] Heatmap guardado: {output_path}")
        return str(output_path)

    except Exception as e:
        logger.error(f"[viz_utils] Error generando heatmap: {e}")
        plt.close("all")
        return None


# ------------------------------------------------------------------
# Dot Plot de Pathways
# ------------------------------------------------------------------

def plot_pathway_dotplot(enrichment_df: pd.DataFrame, output_path: Union[str, Path],
                          pathway_col: str = "Term", padj_col: str = "Adjusted P-value",
                          score_col: str = "combined_score", gene_ratio_col: str = "gene_ratio",
                          top_n: int = 20, title: str = "Pathway Enrichment") -> Optional[str]:
    """
    Genera un dot plot de enriquecimiento de pathways.

    Eje X: Combined score | Eje Y: Pathways | Tamaño: gene ratio | Color: -log10(padj)
    """
    if not _check_matplotlib():
        return None

    try:
        df = enrichment_df.copy().dropna(subset=[padj_col]).head(top_n)
        if df.empty:
            logger.warning("[viz_utils] Sin datos para dot plot de pathways")
            return None

        df = df.sort_values(score_col if score_col in df.columns else padj_col,
                             ascending=True).tail(top_n)
        df["-log10padj"] = -np.log10(df[padj_col].clip(lower=1e-300))
        df["gene_ratio_val"] = df[gene_ratio_col] if gene_ratio_col in df.columns else 0.1

        fig, ax = plt.subplots(figsize=(9, max(5, len(df) * 0.4)))

        norm = Normalize(vmin=df["-log10padj"].min(), vmax=df["-log10padj"].max())
        cmap = plt.cm.get_cmap("RdYlBu_r")

        for i, (_, row) in enumerate(df.iterrows()):
            color = cmap(norm(row["-log10padj"]))
            size = max(30, row["gene_ratio_val"] * 500)
            ax.scatter(row.get(score_col, 1), i, c=[color], s=size, alpha=0.85, zorder=2)

        ax.set_yticks(range(len(df)))
        pathway_labels = [str(p)[:55] + "..." if len(str(p)) > 55 else str(p)
                          for p in df[pathway_col].tolist()]
        ax.set_yticklabels(pathway_labels, fontsize=8)
        ax.set_xlabel("Combined Score", fontsize=10)
        ax.set_title(title, fontsize=12, fontweight="bold", color=COLORS["primary"])
        ax.grid(axis="x", linestyle="--", alpha=0.4)

        # Colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
        cbar.set_label("-log₁₀(padj)", fontsize=8)

        plt.tight_layout()
        plt.savefig(str(output_path), dpi=DPI, format=FIGURE_FORMAT, bbox_inches="tight")
        plt.close(fig)

        logger.debug(f"[viz_utils] Dot plot guardado: {output_path}")
        return str(output_path)

    except Exception as e:
        logger.error(f"[viz_utils] Error generando dot plot: {e}")
        plt.close("all")
        return None


# ------------------------------------------------------------------
# Forest Plot (meta-análisis)
# ------------------------------------------------------------------

def plot_forest(meta_df: pd.DataFrame, output_path: Union[str, Path],
                gene_col: str = "gene", lfc_col: str = "log2FC",
                top_n: int = 25, title: str = "Forest Plot — Top Consensus Genes") -> Optional[str]:
    """
    Genera un forest plot de los genes consenso del meta-análisis.
    Muestra el log2FC meta con intervalo de confianza 95%.
    """
    if not _check_matplotlib():
        return None

    try:
        df = meta_df.dropna(subset=[lfc_col]).head(top_n).copy()
        if df.empty:
            return None

        # Estimar CI aproximado si no hay SE
        se_col = "lfcSE" if "lfcSE" in df.columns else None
        if se_col:
            df["ci_lower"] = df[lfc_col] - 1.96 * df[se_col]
            df["ci_upper"] = df[lfc_col] + 1.96 * df[se_col]
        else:
            df["ci_lower"] = df[lfc_col] * 0.7
            df["ci_upper"] = df[lfc_col] * 1.3

        df = df.sort_values(lfc_col, ascending=True)
        fig, ax = plt.subplots(figsize=(9, max(6, len(df) * 0.35)))

        colors_list = [COLORS["up"] if v > 0 else COLORS["down"] for v in df[lfc_col]]

        ax.errorbar(
            x=df[lfc_col],
            y=range(len(df)),
            xerr=[df[lfc_col] - df["ci_lower"], df["ci_upper"] - df[lfc_col]],
            fmt="o", color="none", ecolor=colors_list,
            elinewidth=1.5, capsize=3, markersize=6, zorder=2,
        )

        # Scatter sobre los puntos
        ax.scatter(df[lfc_col], range(len(df)), c=colors_list, s=40, zorder=3)

        ax.axvline(0, color="black", linestyle="--", linewidth=0.8, alpha=0.6)
        ax.set_yticks(range(len(df)))
        ax.set_yticklabels(df[gene_col].tolist(), fontsize=8)
        ax.set_xlabel("log₂ Fold Change (meta-analysis)", fontsize=10)
        ax.set_title(title, fontsize=12, fontweight="bold", color=COLORS["primary"])
        ax.grid(axis="x", linestyle=":", alpha=0.4)

        plt.tight_layout()
        plt.savefig(str(output_path), dpi=DPI, format=FIGURE_FORMAT, bbox_inches="tight")
        plt.close(fig)

        logger.debug(f"[viz_utils] Forest plot guardado: {output_path}")
        return str(output_path)

    except Exception as e:
        logger.error(f"[viz_utils] Error generando forest plot: {e}")
        plt.close("all")
        return None