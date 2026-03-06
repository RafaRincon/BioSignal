"""
BioSignal Discovery Engine
Utils: Statistical Functions
==============================
Funciones estadísticas compartidas entre agentes del pipeline.
Incluye corrección múltiple, combinación de p-values, y validación.
"""

from typing import Union
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from loguru import logger


# ------------------------------------------------------------------
# Corrección de p-values múltiples
# ------------------------------------------------------------------

def fdr_bh(pvalues: Union[list, np.ndarray, pd.Series], alpha: float = 0.05) -> np.ndarray:
    """
    Corrección Benjamini-Hochberg (FDR).

    Args:
        pvalues: Array de p-valores
        alpha:   Nivel de significancia

    Returns:
        Array de p-valores ajustados
    """
    pvals = np.array(pvalues, dtype=float)
    pvals_clean = np.where(np.isnan(pvals), 1.0, pvals)
    _, padj, _, _ = multipletests(pvals_clean, alpha=alpha, method="fdr_bh")
    return padj


def bonferroni(pvalues: Union[list, np.ndarray, pd.Series]) -> np.ndarray:
    """Corrección de Bonferroni."""
    pvals = np.array(pvalues, dtype=float)
    n = len(pvals)
    return np.minimum(pvals * n, 1.0)


# ------------------------------------------------------------------
# Combinación de p-values (meta-análisis)
# ------------------------------------------------------------------

def fisher_combined_pvalue(pvalues: Union[list, np.ndarray]) -> float:
    """
    Fisher's combined probability test.
    Combina p-values de estudios independientes.

    H0: Todos los p-values son uniformes (sin efecto)
    Ha: Al menos un estudio tiene efecto real

    Args:
        pvalues: Lista de p-values de K estudios independientes

    Returns:
        p-valor combinado
    """
    pvals = np.array(pvalues, dtype=float)
    pvals = pvals[~np.isnan(pvals)]
    pvals = np.clip(pvals, 1e-300, 1.0)

    if len(pvals) == 0:
        return 1.0

    chi2_stat = -2 * np.sum(np.log(pvals))
    df = 2 * len(pvals)
    combined_pval = 1 - stats.chi2.cdf(chi2_stat, df)
    return float(combined_pval)


def stouffer_combined_pvalue(pvalues: Union[list, np.ndarray],
                              weights: Union[list, np.ndarray] = None) -> float:
    """
    Stouffer's Z-score method para combinar p-values.
    Permite ponderar por tamaño de muestra (más robusto que Fisher).

    Args:
        pvalues: Lista de p-values
        weights: Pesos opcionales (ej. sqrt(n_samples))

    Returns:
        p-valor combinado via Z combinado
    """
    pvals = np.array(pvalues, dtype=float)
    pvals = np.clip(pvals, 1e-300, 1.0 - 1e-10)

    z_scores = stats.norm.ppf(1 - pvals)

    if weights is not None:
        w = np.array(weights, dtype=float)
        w = w / w.sum()
        z_combined = np.sum(w * z_scores) / np.sqrt(np.sum(w ** 2))
    else:
        z_combined = np.sum(z_scores) / np.sqrt(len(z_scores))

    combined_pval = 1 - stats.norm.cdf(z_combined)
    return float(combined_pval)


# ------------------------------------------------------------------
# Estadísticas descriptivas
# ------------------------------------------------------------------

def robust_mean(values: Union[list, np.ndarray], trim: float = 0.1) -> float:
    """
    Media recortada (trimmed mean) para mayor robustez ante outliers.

    Args:
        values: Array de valores
        trim:   Fracción a recortar de cada extremo (0.0 - 0.5)

    Returns:
        Media recortada
    """
    arr = np.array(values, dtype=float)
    arr = arr[~np.isnan(arr)]
    if len(arr) == 0:
        return np.nan
    return float(stats.trim_mean(arr, trim))


def cohens_d(group1: np.ndarray, group2: np.ndarray) -> float:
    """
    Tamaño del efecto de Cohen's d.

    Args:
        group1, group2: Arrays de los dos grupos

    Returns:
        Cohen's d (valor positivo si group1 > group2)
    """
    g1 = np.array(group1, dtype=float)
    g2 = np.array(group2, dtype=float)
    g1 = g1[~np.isnan(g1)]
    g2 = g2[~np.isnan(g2)]

    if len(g1) < 2 or len(g2) < 2:
        return np.nan

    pooled_std = np.sqrt(
        ((len(g1) - 1) * np.var(g1, ddof=1) + (len(g2) - 1) * np.var(g2, ddof=1))
        / (len(g1) + len(g2) - 2)
    )
    if pooled_std == 0:
        return 0.0

    return float((np.mean(g1) - np.mean(g2)) / pooled_std)


# ------------------------------------------------------------------
# Validación
# ------------------------------------------------------------------

def validate_pvalue_vector(pvalues: np.ndarray, label: str = "") -> np.ndarray:
    """
    Valida y limpia un vector de p-valores.
    Reemplaza NaN, valores fuera de [0,1] con 1.0 y registra warnings.
    """
    pvals = np.array(pvalues, dtype=float)
    n_nan = np.isnan(pvals).sum()
    n_out = ((pvals < 0) | (pvals > 1)).sum()

    if n_nan > 0:
        logger.warning(f"[stats_utils] {label}: {n_nan} p-valores NaN → reemplazados con 1.0")
    if n_out > 0:
        logger.warning(f"[stats_utils] {label}: {n_out} p-valores fuera de [0,1] → reemplazados con 1.0")

    pvals = np.where(np.isnan(pvals) | (pvals < 0) | (pvals > 1), 1.0, pvals)
    return pvals


def compute_volcano_significance(log2fc: pd.Series, padj: pd.Series,
                                  lfc_threshold: float = 1.0,
                                  padj_threshold: float = 0.05) -> pd.Series:
    """
    Clasifica genes para plot de volcán.

    Returns:
        pd.Series con categorías: 'UP', 'DOWN', 'NS'
    """
    conditions = [
        (log2fc >= lfc_threshold) & (padj < padj_threshold),
        (log2fc <= -lfc_threshold) & (padj < padj_threshold),
    ]
    choices = ["UP", "DOWN"]
    return pd.Series(np.select(conditions, choices, default="NS"), index=log2fc.index)