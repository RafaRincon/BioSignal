"""
Tests para MetaAnalysisAgent — el agente core del sistema.
"""

import json
import numpy as np
import pandas as pd
import pytest
from pathlib import Path

from agents.meta_analysis import MetaAnalysisAgent


@pytest.fixture
def agent():
    """Instancia de MetaAnalysisAgent con configuración de test."""
    return MetaAnalysisAgent(config={
        "padj_threshold": 0.05,
        "min_dataset_fraction": 0.60,
        "direction_threshold": 0.70,
        "generate_heatmap": False,
    })


@pytest.fixture
def mock_dea_results():
    """Datos de DEA mock para 3 datasets."""
    np.random.seed(42)

    # Genes que deben ser consenso (consistentemente UP)
    consensus_up = ["KRAS", "MUC1", "MMP9", "S100A4", "CEACAM5"]
    # Genes que deben ser consenso (consistentemente DOWN)
    consensus_down = ["CFTR", "SPINK1", "CDH2"]
    # Genes ruido (inconsistentes)
    noise_genes = [f"GENE{i}" for i in range(50)]

    results = {}
    for dataset_id in ["GSE71989", "GSE62452", "GSE15471"]:
        genes = consensus_up + consensus_down + noise_genes
        n = len(genes)

        log2fc = np.random.normal(0, 0.5, n)
        pvalue = np.random.uniform(0.1, 0.9, n)
        padj = pvalue * n  # Rough BH adjustment

        # Genes consenso UP: log2FC positivo, p-value bajo
        for i, g in enumerate(consensus_up):
            idx = genes.index(g)
            log2fc[idx] = np.random.uniform(1.5, 3.0)
            pvalue[idx] = np.random.uniform(0.001, 0.04)
            padj[idx] = pvalue[idx] * 5

        # Genes consenso DOWN: log2FC negativo, p-value bajo
        for g in consensus_down:
            idx = genes.index(g)
            log2fc[idx] = np.random.uniform(-3.0, -1.5)
            pvalue[idx] = np.random.uniform(0.001, 0.04)
            padj[idx] = pvalue[idx] * 5

        results[dataset_id] = pd.DataFrame({
            "gene": genes,
            "log2FC": log2fc,
            "pvalue": np.clip(pvalue, 1e-10, 1.0),
            "padj": np.clip(padj, 1e-10, 1.0),
            "baseMean": np.random.uniform(50, 5000, n),
        })

    return results


class TestFisherCombinedPvalue:
    """Tests para el método de Fisher."""

    def test_combines_significant_pvalues(self, agent):
        """Combinar p-values pequeños debe resultar en p-value aún más pequeño."""
        pvalues = np.array([0.01, 0.02, 0.005])
        result = agent.fisher_combined_pvalue(pvalues)
        assert result < 0.01
        assert 0 < result < 1

    def test_combines_nonsignificant_pvalues(self, agent):
        """Combinar p-values grandes debe resultar en p-value no significativo."""
        pvalues = np.array([0.5, 0.6, 0.7])
        result = agent.fisher_combined_pvalue(pvalues)
        assert result > 0.05

    def test_handles_nan_values(self, agent):
        """NaN en p-values deben ser ignorados."""
        pvalues = np.array([0.01, np.nan, 0.005, np.nan])
        result = agent.fisher_combined_pvalue(pvalues)
        assert not np.isnan(result)
        assert 0 < result < 1

    def test_insufficient_pvalues_returns_nan(self, agent):
        """Menos de 2 p-values válidos debe retornar NaN."""
        pvalues = np.array([0.01, np.nan, np.nan])
        result = agent.fisher_combined_pvalue(pvalues)
        assert np.isnan(result)


class TestBuildGeneMatrix:
    """Tests para la construcción de matrices gene x dataset."""

    def test_builds_correct_dimensions(self, agent, mock_dea_results):
        """La matriz debe tener genes como filas y datasets como columnas."""
        log2fc_matrix, pvalue_matrix = agent.build_gene_matrix(mock_dea_results)

        n_datasets = len(mock_dea_results)
        assert log2fc_matrix.shape[1] == n_datasets
        assert pvalue_matrix.shape[1] == n_datasets

    def test_gene_alignment(self, agent, mock_dea_results):
        """Genes deben estar alineados entre datasets correctamente."""
        log2fc_matrix, _ = agent.build_gene_matrix(mock_dea_results)

        # KRAS debe estar en la matriz
        assert "KRAS" in log2fc_matrix.index

    def test_handles_missing_genes(self, agent):
        """Genes no presentes en un dataset deben aparecer como NaN."""
        results = {
            "GSE001": pd.DataFrame({
                "gene": ["KRAS", "MUC1", "CFTR"],
                "log2FC": [2.0, -1.5, 3.0],
                "pvalue": [0.01, 0.02, 0.03],
                "padj": [0.05, 0.06, 0.07],
            }),
            "GSE002": pd.DataFrame({
                "gene": ["KRAS", "MMP9"],  # Falta MUC1 y CFTR
                "log2FC": [1.8, -2.0],
                "pvalue": [0.01, 0.03],
                "padj": [0.05, 0.08],
            }),
        }
        log2fc_matrix, _ = agent.build_gene_matrix(results)

        # MUC1 debe ser NaN en GSE002
        assert pd.isna(log2fc_matrix.loc["MUC1", "GSE002"])
        # KRAS debe tener valor en ambos datasets
        assert not pd.isna(log2fc_matrix.loc["KRAS", "GSE001"])
        assert not pd.isna(log2fc_matrix.loc["KRAS", "GSE002"])


class TestRunMetaAnalysis:
    """Tests para el pipeline completo de meta-análisis."""

    def test_identifies_consensus_genes(self, agent, mock_dea_results):
        """El meta-análisis debe identificar los genes UP/DOWN consistentes."""
        log2fc_matrix, pvalue_matrix = agent.build_gene_matrix(mock_dea_results)
        significant = agent.run_meta_analysis(log2fc_matrix, pvalue_matrix)

        gene_set = set(significant["gene"].values)
        # Los genes consenso UP deben aparecer
        assert "KRAS" in gene_set or "MUC1" in gene_set

    def test_returns_dataframe(self, agent, mock_dea_results):
        """El resultado debe ser un DataFrame con las columnas esperadas."""
        log2fc_matrix, pvalue_matrix = agent.build_gene_matrix(mock_dea_results)
        result = agent.run_meta_analysis(log2fc_matrix, pvalue_matrix)

        assert isinstance(result, pd.DataFrame)
        required_cols = {"gene", "meta_pvalue", "meta_padj", "mean_log2fc", "direction"}
        assert required_cols.issubset(set(result.columns))

    def test_direction_assignment(self, agent, mock_dea_results):
        """Genes UP deben tener dirección 'UP', genes DOWN 'DOWN'."""
        log2fc_matrix, pvalue_matrix = agent.build_gene_matrix(mock_dea_results)
        significant = agent.run_meta_analysis(log2fc_matrix, pvalue_matrix)

        if "KRAS" in significant["gene"].values:
            kras_row = significant[significant["gene"] == "KRAS"].iloc[0]
            assert kras_row["direction"] == "UP"

        if "CFTR" in significant["gene"].values:
            cftr_row = significant[significant["gene"] == "CFTR"].iloc[0]
            assert cftr_row["direction"] == "DOWN"

    def test_padj_threshold_filter(self, agent, mock_dea_results):
        """Todos los genes en el output deben tener meta_padj < threshold."""
        log2fc_matrix, pvalue_matrix = agent.build_gene_matrix(mock_dea_results)
        significant = agent.run_meta_analysis(log2fc_matrix, pvalue_matrix)

        if len(significant) > 0:
            assert (significant["meta_padj"] < agent.padj_threshold).all()

    def test_dataset_fraction_filter(self, agent, mock_dea_results):
        """Genes deben estar presentes en la fracción mínima de datasets."""
        log2fc_matrix, pvalue_matrix = agent.build_gene_matrix(mock_dea_results)
        significant = agent.run_meta_analysis(log2fc_matrix, pvalue_matrix)

        if len(significant) > 0:
            assert (significant["dataset_fraction"] >= agent.min_dataset_fraction).all()


class TestBuildConsensusGenes:
    """Tests para la construcción de genes consenso."""

    def test_consensus_gene_structure(self, agent, mock_dea_results):
        """Cada gen consenso debe tener la estructura correcta."""
        log2fc_matrix, pvalue_matrix = agent.build_gene_matrix(mock_dea_results)
        significant = agent.run_meta_analysis(log2fc_matrix, pvalue_matrix)
        consensus = agent.build_consensus_genes(significant, log2fc_matrix)

        if consensus:
            first = consensus[0]
            required_keys = {
                "gene", "direction", "meta_padj", "mean_log2fc",
                "n_datasets_present", "dataset_fraction", "per_dataset_log2fc"
            }
            assert required_keys.issubset(set(first.keys()))

    def test_consensus_genes_serializable(self, agent, mock_dea_results):
        """Los genes consenso deben ser serializables a JSON."""
        log2fc_matrix, pvalue_matrix = agent.build_gene_matrix(mock_dea_results)
        significant = agent.run_meta_analysis(log2fc_matrix, pvalue_matrix)
        consensus = agent.build_consensus_genes(significant, log2fc_matrix)

        # Esto no debe lanzar ninguna excepción
        serialized = json.dumps(consensus)
        assert len(serialized) > 0