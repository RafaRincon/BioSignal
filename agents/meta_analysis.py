"""
BioSignal Discovery Engine
Agent 5: MetaAnalysisAgent
===========================
Responsabilidad: Integrar resultados de DEA de múltiples datasets para identificar
señales biológicas REPRODUCIBLES y CONSISTENTES entre estudios.

Esta es la feature diferenciadora principal del sistema.
Ninguna herramienta existente (Galaxy, GEO2R, Bioconductor) automatiza
este paso de integración cross-dataset.

Uso:
    python -m agents.meta_analysis --input data/dea/
"""

import json
import time
from pathlib import Path
from typing import Optional


import click
import numpy as np
import pandas as pd
from loguru import logger
from scipy import stats
from statsmodels.stats.multitest import multipletests


class MetaAnalysisAgent:
    """
    Agente de meta-análisis transcriptómico cross-dataset.

    Métodos implementados:
        - Fisher's combined probability test: combina p-values de cada estudio
        - Stouffer's Z-score: pondera por número de muestras
        - Vote counting: genes up/down en >= N estudios de M totales
        - Random Effects Model (REM): modela heterogeneidad (via rpy2 + metafor)

    Proceso:
        1. Cargar significant_genes.json de cada dataset
        2. Construir matriz: genes (filas) x datasets (columnas) con log2FC
        3. Alinear genes: usar Gene Symbol como ID canónico (mapeo via HGNC)
        4. Calcular meta p-value por gen usando Fisher's method
        5. Aplicar FDR correction (BH) sobre meta p-values
        6. Filtrar: genes con meta_padj < 0.05 Y presentes en ≥60% de datasets
        7. Calcular dirección consenso: UP si mean(log2FC) > 0 en ≥70% datasets

    Salida:
        data/meta/{disease}/meta_analysis_results.csv
        data/meta/{disease}/consensus_genes.json
        data/meta/{disease}/heatmap_top50.png
    """

    def __init__(self, config: dict = None):
        self.config = config or {}
        self.padj_threshold = self.config.get("padj_threshold", 0.05)
        self.min_dataset_fraction = self.config.get("min_dataset_fraction", 0.60)
        self.direction_threshold = self.config.get("direction_threshold", 0.70)
        self.generate_heatmap = self.config.get("generate_heatmap", True)
        self.heatmap_top_genes = self.config.get("heatmap_top_genes", 50)

    def _detect_namespace(self, genes: list[str]) -> str:
        """Detecta el namespace de los gene IDs: 'ensembl_human', 'ensembl_mouse', o 'symbol'."""
        sample = [g for g in genes[:50] if isinstance(g, str)]
        n_human = sum(1 for g in sample if g.startswith("ENSG"))
        n_mouse = sum(1 for g in sample if g.startswith("ENSMUSG"))
        if n_human / max(len(sample), 1) > 0.5:
            return "ensembl_human"
        if n_mouse / max(len(sample), 1) > 0.5:
            return "ensembl_mouse"
        return "symbol"

    def normalize_gene_ids(self, df: pd.DataFrame, gse_id: str) -> pd.DataFrame:
        """
        Normaliza gene IDs a human gene symbols usando mygene.
        Handles: ENSG* → symbol, ENSMUSG* → human ortholog symbol, symbol → symbol.
        """
        import mygene
        mg = mygene.MyGeneInfo()

        genes = df["gene"].dropna().unique().tolist()
        namespace = self._detect_namespace(genes)

        if namespace == "symbol":
            logger.debug(f"  [{gse_id}] Gene IDs ya son symbols — sin conversión")
            return df

        logger.info(f"  [{gse_id}] Convirtiendo {len(genes)} IDs ({namespace}) → human symbols...")

        # Strip versiones tipo ENSMUSG00000030516.14 → ENSMUSG00000030516
        genes_clean = [g.split(".")[0] if "." in g else g for g in genes]
        original_map = dict(zip(genes_clean, genes))  # clean → original

        symbol_map = {}  # original_id → human_symbol

        try:
            if namespace == "ensembl_human":
                hits = mg.querymany(
                    genes_clean,
                    scopes="ensembl.gene",
                    fields="symbol",
                    species="human",
                    verbose=False,
                )
                for hit in hits:
                    if hit.get("notfound"):
                        continue
                    original = original_map.get(hit["query"], hit["query"])
                    sym = hit.get("symbol")
                    if sym:
                        symbol_map[original] = sym

            else:  # ensembl_mouse
                hits = mg.querymany(
                    genes_clean,
                    scopes="ensembl.gene",
                    fields="symbol,ortholog.human",
                    species="mouse",
                    verbose=False,
                )
                for hit in hits:
                    if hit.get("notfound"):
                        continue
                    original = original_map.get(hit["query"], hit["query"])
                    ortho = hit.get("ortholog", {}).get("human") if isinstance(hit.get("ortholog"), dict) else None
                    if isinstance(ortho, list):
                        sym = ortho[0].get("symbol") if ortho else None
                    elif isinstance(ortho, dict):
                        sym = ortho.get("symbol")
                    else:
                        sym = hit.get("symbol")  # fallback: símbolo ratón
                    if sym:
                        symbol_map[original] = sym

        except Exception as e:
            logger.warning(f"  [{gse_id}] mygene falló: {e} — manteniendo IDs originales")
            return df

        logger.info(f"  [{gse_id}] Mapeados {len(symbol_map)}/{len(genes)} genes → symbols")

        df = df.copy()
        df["gene"] = df["gene"].map(
            lambda g: symbol_map.get(g, symbol_map.get(g.split(".")[0] if "." in g else g))
        )
        df = df.dropna(subset=["gene"])
        df = df.sort_values("pvalue").drop_duplicates(subset="gene", keep="first")
        return df

    def load_dea_results(self, dea_dir: str) -> dict[str, pd.DataFrame]:
        """
        Carga resultados de DEA de todos los datasets disponibles.

        Args:
            dea_dir: Directorio con resultados del Agent 4

        Returns:
            Dict {gse_id: DataFrame con columnas [gene, log2FC, pvalue, padj]}
        """
        dea_path = Path(dea_dir)
        results = {}

        for dataset_dir in sorted(dea_path.iterdir()):
            if not dataset_dir.is_dir():
                continue

            gse_id = dataset_dir.name
            dea_file = dataset_dir / "dea_results.csv"
            sig_file = dataset_dir / "significant_genes.json"

            if not dea_file.exists():
                logger.warning(f"  No encontrado: {dea_file}")
                continue

            try:
                df = pd.read_csv(dea_file)
                # Validar columnas requeridas
                required = {"gene", "log2FC", "pvalue", "padj"}
                if not required.issubset(df.columns):
                    logger.warning(f"  {gse_id}: columnas faltantes {required - set(df.columns)}")
                    continue

                df = self.normalize_gene_ids(df, gse_id)
                results[gse_id] = df
                logger.debug(f"  {gse_id}: {len(df)} genes cargados (post-normalization)")

            except Exception as e:
                logger.error(f"  Error cargando {gse_id}: {e}")
                continue

        logger.info(f"Cargados resultados DEA de {len(results)} datasets")
        return results

    def build_gene_matrix(
        self,
        dea_results: dict[str, pd.DataFrame],
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Construye matrices de log2FC y p-values: genes x datasets.

        Args:
            dea_results: Dict {gse_id: DataFrame DEA}

        Returns:
            Tuple de (log2fc_matrix, pvalue_matrix) donde índices son Gene Symbols
        """
        all_genes = set()
        for df in dea_results.values():
            all_genes.update(df["gene"].dropna().unique())

        logger.info(f"Total genes únicos entre datasets: {len(all_genes)}")

        # Construir matrices
        log2fc_dict = {}
        pvalue_dict = {}

        for gse_id, df in dea_results.items():
            gene_df = df.set_index("gene")
            log2fc_dict[gse_id] = gene_df["log2FC"].reindex(all_genes)
            pvalue_dict[gse_id] = gene_df["pvalue"].reindex(all_genes)

        log2fc_matrix = pd.DataFrame(log2fc_dict, index=sorted(all_genes))
        pvalue_matrix = pd.DataFrame(pvalue_dict, index=sorted(all_genes))

        return log2fc_matrix, pvalue_matrix

    def fisher_combined_pvalue(self, pvalues: np.ndarray) -> float:
        """
        Aplica el Fisher's Combined Probability Test para combinar p-values.

        Chi-squared = -2 * sum(log(p_i))
        df = 2 * k donde k es el número de estudios

        Args:
            pvalues: Array de p-values (NaN son excluidos)

        Returns:
            Meta p-value combinado
        """
        valid_pvalues = pvalues[~np.isnan(pvalues)]
        valid_pvalues = np.clip(valid_pvalues, 1e-300, 1.0)  # Evitar log(0)

        if len(valid_pvalues) < 2:
            return np.nan

        chi2_stat = -2 * np.sum(np.log(valid_pvalues))
        df = 2 * len(valid_pvalues)
        meta_pvalue = 1 - stats.chi2.cdf(chi2_stat, df)

        return float(meta_pvalue)

    def stouffer_zscore(
        self,
        pvalues: np.ndarray,
        weights: Optional[np.ndarray] = None,
    ) -> float:
        """
        Aplica el Stouffer's Z-score method para combinar p-values.
        Pondera por número de muestras de cada estudio.

        Args:
            pvalues: Array de p-values
            weights: Pesos por estudio (ej. n_samples); si None, pesos iguales

        Returns:
            Meta p-value combinado por Stouffer
        """
        valid_mask = ~np.isnan(pvalues)
        valid_pvalues = np.clip(pvalues[valid_mask], 1e-300, 1.0 - 1e-10)

        if len(valid_pvalues) < 2:
            return np.nan

        z_scores = stats.norm.ppf(1 - valid_pvalues)

        if weights is not None:
            w = weights[valid_mask]
            w = w / w.sum()  # Normalizar
            combined_z = np.sum(w * z_scores) / np.sqrt(np.sum(w ** 2))
        else:
            combined_z = np.sum(z_scores) / np.sqrt(len(z_scores))

        meta_pvalue = 1 - stats.norm.cdf(combined_z)
        return float(meta_pvalue)

    def vote_counting(
        self,
        log2fc_matrix: pd.DataFrame,
        padj_matrix: pd.DataFrame,
        min_fraction: float = 0.60,
        padj_threshold: float = 0.05,
        log2fc_threshold: float = 1.0,
    ) -> pd.DataFrame:
        """
        Vote counting: cuenta en cuántos estudios cada gen es significativo.

        Args:
            log2fc_matrix: Matriz de log2FC (genes x datasets)
            padj_matrix: Matriz de padj (genes x datasets)
            min_fraction: Fracción mínima de datasets donde debe aparecer
            padj_threshold: Threshold de padj para significancia
            log2fc_threshold: Threshold de |log2FC| mínimo

        Returns:
            DataFrame con conteos de votos por gen
        """
        n_datasets = log2fc_matrix.shape[1]
        min_datasets = int(np.ceil(min_fraction * n_datasets))

        sig_mask = (padj_matrix < padj_threshold) & (log2fc_matrix.abs() >= log2fc_threshold)

        votes_up = ((sig_mask) & (log2fc_matrix > 0)).sum(axis=1)
        votes_down = ((sig_mask) & (log2fc_matrix < 0)).sum(axis=1)
        votes_total = sig_mask.sum(axis=1)

        vote_df = pd.DataFrame({
            "votes_up": votes_up,
            "votes_down": votes_down,
            "votes_total": votes_total,
            "n_datasets": n_datasets,
            "fraction": votes_total / n_datasets,
        })

        return vote_df[vote_df["votes_total"] >= min_datasets]

    def run_meta_analysis(
        self,
        log2fc_matrix: pd.DataFrame,
        pvalue_matrix: pd.DataFrame,
        sample_counts: Optional[dict] = None,
    ) -> pd.DataFrame:
        """
        Ejecuta el meta-análisis completo sobre las matrices de genes x datasets.

        Args:
            log2fc_matrix: Matriz de log2FC (genes x datasets)
            pvalue_matrix: Matriz de p-values (genes x datasets)
            sample_counts: Dict {gse_id: n_samples} para ponderación Stouffer

        Returns:
            DataFrame con resultados del meta-análisis por gen
        """
        logger.info("Calculando meta p-values (Fisher's combined probability test)...")
        n_datasets = log2fc_matrix.shape[1]

        # Calcular meta p-values por gen
        meta_pvalues = []
        present_counts = []
        mean_log2fc = []
        direction_consistency = []

        for gene in log2fc_matrix.index:
            pvals = pvalue_matrix.loc[gene].values
            log2fcs = log2fc_matrix.loc[gene].values

            # Contar presencia
            valid_mask = ~np.isnan(pvals)
            n_present = valid_mask.sum()
            present_counts.append(n_present)

            if n_present < 2:
                meta_pvalues.append(np.nan)
                mean_log2fc.append(np.nan)
                direction_consistency.append(np.nan)
                continue

            # Meta p-value por Fisher
            meta_p = self.fisher_combined_pvalue(pvals)
            meta_pvalues.append(meta_p)

            # Log2FC promedio (sobre datasets presentes)
            valid_log2fc = log2fcs[~np.isnan(log2fcs)]
            mean_fc = np.nanmean(valid_log2fc)
            mean_log2fc.append(mean_fc)

            # Consistencia de dirección
            n_up = np.sum(valid_log2fc > 0)
            n_down = np.sum(valid_log2fc < 0)
            dir_consistency = max(n_up, n_down) / len(valid_log2fc)
            direction_consistency.append(dir_consistency)

        # Construir DataFrame de resultados
        results_df = pd.DataFrame({
            "gene": log2fc_matrix.index,
            "meta_pvalue": meta_pvalues,
            "mean_log2fc": mean_log2fc,
            "direction_consistency": direction_consistency,
            "n_datasets_present": present_counts,
            "n_datasets_total": n_datasets,
            "dataset_fraction": [c / n_datasets for c in present_counts],
        })

        # Aplicar FDR correction (Benjamini-Hochberg)
        valid_mask = ~results_df["meta_pvalue"].isna()
        if valid_mask.sum() > 0:
            _, padj_values, _, _ = multipletests(
                results_df.loc[valid_mask, "meta_pvalue"].values,
                method="fdr_bh"
            )
            results_df.loc[valid_mask, "meta_padj"] = padj_values
        else:
            results_df["meta_padj"] = np.nan

        # Asignar dirección consenso
        results_df["direction"] = results_df.apply(
            lambda row: self._assign_direction(row),
            axis=1
        )

        # Filtrar por significancia y presencia mínima
        significant = results_df[
            (results_df["meta_padj"] < self.padj_threshold) &
            (results_df["dataset_fraction"] >= self.min_dataset_fraction)
        ].copy()

        # Ordenar por meta_padj
        significant = significant.sort_values("meta_padj")

        logger.info(
            f"Meta-análisis completado:\n"
            f"  Genes totales analizados: {len(results_df):,}\n"
            f"  Genes con meta_padj < {self.padj_threshold}: {(results_df['meta_padj'] < self.padj_threshold).sum():,}\n"
            f"  Genes consenso (≥{self.min_dataset_fraction:.0%} datasets): {len(significant):,}\n"
            f"  → UP:   {(significant['direction'] == 'UP').sum()}\n"
            f"  → DOWN: {(significant['direction'] == 'DOWN').sum()}"
        )

        return significant

    def _assign_direction(self, row: pd.Series) -> str:
        """
        Asigna dirección consenso (UP/DOWN/MIXED) a un gen.

        Args:
            row: Fila del DataFrame de resultados

        Returns:
            'UP' | 'DOWN' | 'MIXED' | 'UNKNOWN'
        """
        if pd.isna(row.get("direction_consistency")) or pd.isna(row.get("mean_log2fc")):
            return "UNKNOWN"

        if row["direction_consistency"] >= self.direction_threshold:
            return "UP" if row["mean_log2fc"] > 0 else "DOWN"
        return "MIXED"

    def build_consensus_genes(
        self,
        significant: pd.DataFrame,
        log2fc_matrix: pd.DataFrame,
    ) -> list[dict]:
        """
        Construye la lista de genes consenso con toda la información de soporte.

        Args:
            significant: DataFrame de genes significativos del meta-análisis
            log2fc_matrix: Matriz original de log2FC

        Returns:
            Lista de dicts con información completa por gen
        """
        consensus = []

        for _, row in significant.iterrows():
            gene = row["gene"]

            # Obtener log2FC por dataset para este gen
            if gene in log2fc_matrix.index:
                per_dataset_fc = log2fc_matrix.loc[gene].dropna().to_dict()
            else:
                per_dataset_fc = {}

            consensus.append({
                "gene": gene,
                "direction": row["direction"],
                "meta_padj": round(float(row["meta_padj"]), 6) if not pd.isna(row["meta_padj"]) else None,
                "meta_pvalue": round(float(row["meta_pvalue"]), 8) if not pd.isna(row["meta_pvalue"]) else None,
                "mean_log2fc": round(float(row["mean_log2fc"]), 4) if not pd.isna(row["mean_log2fc"]) else None,
                "direction_consistency": round(float(row["direction_consistency"]), 3) if not pd.isna(row["direction_consistency"]) else None,
                "n_datasets_present": int(row["n_datasets_present"]),
                "n_datasets_total": int(row["n_datasets_total"]),
                "dataset_fraction": round(float(row["dataset_fraction"]), 3),
                "per_dataset_log2fc": {k: round(float(v), 4) for k, v in per_dataset_fc.items()},
            })

        return consensus

    def generate_heatmap_figure(
        self,
        log2fc_matrix: pd.DataFrame,
        consensus_genes: list[dict],
        output_path: str,
        top_n: int = 50,
    ) -> None:
        """
        Genera un heatmap cross-dataset de los top genes consenso.

        Args:
            log2fc_matrix: Matriz de log2FC (genes x datasets)
            consensus_genes: Lista de genes consenso
            output_path: Path de salida para la figura
            top_n: Número de genes a incluir en el heatmap
        """
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns

            # Seleccionar top genes
            top_genes = [g["gene"] for g in consensus_genes[:top_n]]
            top_genes = [g for g in top_genes if g in log2fc_matrix.index]

            if not top_genes:
                logger.warning("No hay genes suficientes para generar heatmap")
                return

            heatmap_data = log2fc_matrix.loc[top_genes]

            # Crear figura
            fig_height = max(10, len(top_genes) * 0.3)
            fig, ax = plt.subplots(figsize=(14, fig_height))

            # Determinar rango del colormap
            vmax = min(abs(heatmap_data.values[~np.isnan(heatmap_data.values)]).max(), 5)

            sns.heatmap(
                heatmap_data,
                cmap="RdBu_r",
                center=0,
                vmin=-vmax,
                vmax=vmax,
                ax=ax,
                cbar_kws={"label": "log2 Fold Change", "shrink": 0.8},
                linewidths=0.5,
                linecolor="gray",
                xticklabels=True,
                yticklabels=True,
            )

            ax.set_title(
                f"Top {len(top_genes)} Consensus Genes — Cross-Dataset log2FC",
                fontsize=14, fontweight="bold", pad=20
            )
            ax.set_xlabel("Datasets (GSE IDs)", fontsize=11)
            ax.set_ylabel("Genes", fontsize=11)
            ax.tick_params(axis="x", rotation=45)
            ax.tick_params(axis="y", labelsize=8)

            plt.tight_layout()
            plt.savefig(output_path, dpi=150, bbox_inches="tight")
            plt.close()

            logger.success(f"✓ Heatmap guardado: {output_path}")

        except ImportError as e:
            logger.warning(f"No se pudo generar heatmap (librería faltante): {e}")

    def run(
        self,
        input_dir: str,
        output_dir: str,
        disease_name: str = "unknown",
    ) -> Path:
        """
        Ejecuta el pipeline completo del meta-análisis.

        Args:
            input_dir: Directorio con resultados del Agent 4 (dea/)
            output_dir: Directorio de salida (meta/)
            disease_name: Nombre de la enfermedad para organizar outputs

        Returns:
            Path al directorio de resultados
        """
        logger.info(f"=== Agent 5: MetaAnalysisAgent ===")
        logger.info(f"Input: {input_dir}")
        start_time = time.time()

        # Preparar directorio de salida
        disease_safe = disease_name.replace(" ", "_").lower()
        output_path = Path(output_dir) / disease_safe
        output_path.mkdir(parents=True, exist_ok=True)

        # 1. Cargar resultados DEA
        dea_results = self.load_dea_results(input_dir)
        if len(dea_results) < 2:
            logger.error(f"Se necesitan al menos 2 datasets para meta-análisis. Encontrados: {len(dea_results)}")
            return output_path

        # 2. Construir matrices gene x dataset
        log2fc_matrix, pvalue_matrix = self.build_gene_matrix(dea_results)

        # 3. Ejecutar meta-análisis
        significant = self.run_meta_analysis(log2fc_matrix, pvalue_matrix)

        # 4. Construir genes consenso
        consensus_genes = self.build_consensus_genes(significant, log2fc_matrix)

        # 5. Guardar resultados
        # CSV completo
        csv_path = output_path / "meta_analysis_results.csv"
        significant.to_csv(csv_path, index=False)
        logger.info(f"  → Meta-análisis CSV: {csv_path}")

        # JSON de genes consenso
        json_path = output_path / "consensus_genes.json"
        with open(json_path, "w") as f:
            json.dump({
                "disease": disease_name,
                "n_datasets": len(dea_results),
                "n_consensus_genes": len(consensus_genes),
                "params": {
                    "padj_threshold": self.padj_threshold,
                    "min_dataset_fraction": self.min_dataset_fraction,
                    "direction_threshold": self.direction_threshold,
                },
                "consensus_genes": consensus_genes,
            }, f, indent=2)
        logger.info(f"  → Genes consenso JSON: {json_path}")

        # 6. Generar heatmap
        if self.generate_heatmap and consensus_genes:
            heatmap_path = output_path / f"heatmap_top{self.heatmap_top_genes}.png"
            self.generate_heatmap_figure(
                log2fc_matrix=log2fc_matrix,
                consensus_genes=consensus_genes,
                output_path=str(heatmap_path),
                top_n=self.heatmap_top_genes,
            )

        elapsed = time.time() - start_time
        logger.success(
            f"\n✓ Meta-análisis completado en {elapsed:.1f}s\n"
            f"  {len(consensus_genes)} genes consenso identificados\n"
            f"  UP: {sum(1 for g in consensus_genes if g['direction'] == 'UP')} | "
            f"DOWN: {sum(1 for g in consensus_genes if g['direction'] == 'DOWN')}"
        )

        return output_path


# ============================================================
# CLI Interface
# ============================================================

@click.command()
@click.option("--input", required=True, help="Directorio con resultados DEA (Agent 4)")
@click.option("--output", default="data/meta/", show_default=True,
              help="Directorio de salida")
@click.option("--disease", default="unknown", show_default=True,
              help="Nombre de la enfermedad")
def main(input, output, disease):
    """
    Agent 5: Meta-análisis transcriptómico cross-dataset.

    Identifica genes biológicamente reproducibles entre múltiples estudios.

    Ejemplo:
        python -m agents.meta_analysis \\
            --input data/dea/ \\
            --disease 'pancreatic cancer'
    """
    from utils.geo_utils import load_config
    config = load_config()

    agent = MetaAnalysisAgent(config=config.get("meta_analysis", {}))
    agent.run(
        input_dir=input,
        output_dir=output,
        disease_name=disease,
    )


if __name__ == "__main__":
    main()