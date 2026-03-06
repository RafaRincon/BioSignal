"""
BioSignal Discovery Engine
Agent 4: DifferentialExpressionAgent
======================================
Responsabilidad: Ejecutar análisis de expresión diferencial (DEA) por dataset.
Método principal: edgeR via rpy2. Fallback automático: limma → t-test Python.

Uso:
    python -m agents.dea --input data/processed/ --output data/dea/
"""

import json
import time
import warnings
from pathlib import Path
from typing import Optional

import click
import numpy as np
import pandas as pd
from loguru import logger

# rpy2 para integración con R
try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri, Formula
    from rpy2.robjects.packages import importr
    from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
    import logging
    rpy2_logger.setLevel(logging.ERROR)  # Silenciar warnings de R
    pandas2ri.activate()
    RPY2_AVAILABLE = True
except ImportError:
    RPY2_AVAILABLE = False
    logger.warning("rpy2 no disponible — DEA usará método Python fallback (limma-like via statsmodels)")


class DifferentialExpressionAgent:
    """
    Agente de análisis de expresión diferencial.

    Método principal: edgeR (rpy2) — robusto para RNA-seq y microarray.
    Fallback automático si edgeR falla:
        edgeR (rpy2) → limma (rpy2) → t-test (Python puro)

    Proceso:
        1. Cargar matrix_normalized.csv y sample_metadata.csv del Agent 3
        2. Ejecutar DEA con edgeR (exactTest o glmQLFTest)
        3. Aplicar filtros: padj < threshold Y |log2FC| > lfc_threshold
        4. Guardar resultados y genes significativos

    Salida por dataset:
        data/dea/{gse_id}/dea_results.csv          # Todos los genes + estadísticas
        data/dea/{gse_id}/significant_genes.json   # Genes que pasan filtros
        data/dea/{gse_id}/volcano_data.csv         # Para plot volcán
        data/dea/{gse_id}/dea_summary.json         # Resumen del análisis
    """

    def __init__(self, config: dict = None):
        self.config = config or {}
        self.padj_threshold = self.config.get("padj_threshold", 0.05)
        self.lfc_threshold = self.config.get("lfc_threshold", 1.0)
        self.prefer_method = self.config.get("method", "edger")  # edger por defecto
        self.min_samples_dea = self.config.get("min_samples_dea", 4)

    # ------------------------------------------------------------------
    # Entrada principal
    # ------------------------------------------------------------------

    def run(self, processed_dir: str, output_dir: str) -> dict:
        """
        Ejecuta DEA en todos los datasets preprocesados que pasaron QC.

        Args:
            processed_dir: Directorio de salida del Agent 3 (data/processed/)
            output_dir:    Directorio destino (data/dea/)

        Returns:
            Resumen {"analyzed": [...], "failed": [...], "skipped": [...]}
        """
        proc_path = Path(processed_dir)
        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        summary = {"analyzed": [], "failed": [], "skipped": []}
        start = time.time()

        dataset_dirs = sorted([d for d in proc_path.iterdir() if d.is_dir()])
        logger.info(f"[DEAAgent] Encontrados {len(dataset_dirs)} datasets procesados")

        for dataset_dir in dataset_dirs:
            gse_id = dataset_dir.name

            # Verificar QC status
            qc_status_file = dataset_dir / "qc_status.txt"
            if qc_status_file.exists() and qc_status_file.read_text().strip() == "FAIL":
                logger.info(f"[{gse_id}] Saltando — QC FAIL")
                summary["skipped"].append(gse_id)
                continue

            try:
                result = self._run_dea_dataset(dataset_dir, out_path / gse_id)
                summary["analyzed"].append({
                    "gse_id": gse_id,
                    "method": result.get("method"),
                    "n_significant": result.get("n_significant", 0),
                    "n_up": result.get("n_up", 0),
                    "n_down": result.get("n_down", 0),
                })
            except Exception as e:
                logger.error(f"[{gse_id}] Error en DEA: {e}")
                summary["failed"].append({"gse_id": gse_id, "error": str(e)})

        elapsed = time.time() - start
        logger.info(
            f"[DEAAgent] Completado en {elapsed:.1f}s | "
            f"Analizados={len(summary['analyzed'])} | "
            f"Saltados={len(summary['skipped'])} | "
            f"Fallidos={len(summary['failed'])}"
        )
        return summary

    # ------------------------------------------------------------------
    # DEA por dataset
    # ------------------------------------------------------------------

    def _run_dea_dataset(self, dataset_dir: Path, out_dir: Path) -> dict:
        """Ejecuta DEA para un dataset individual."""
        out_dir.mkdir(parents=True, exist_ok=True)
        gse_id = dataset_dir.name

        # Cargar datos
        matrix = pd.read_csv(dataset_dir / "matrix_normalized.csv", index_col=0)
        sample_meta = self._load_sample_metadata(dataset_dir)

        if sample_meta is None or "group" not in sample_meta.columns:
            raise ValueError("sample_metadata.csv sin columna 'group'")

        # Alinear muestras
        common_samples = matrix.columns.intersection(sample_meta.index)
        if len(common_samples) < self.min_samples_dea:
            raise ValueError(f"Solo {len(common_samples)} muestras alineadas — mínimo {self.min_samples_dea}")

        matrix = matrix[common_samples]
        sample_meta = sample_meta.loc[common_samples]

        # Cargar QC report para conocer data_type
        qc_report = self._load_qc_report(dataset_dir)
        data_type = qc_report.get("data_type", "unknown")

        # Seleccionar método — edgeR por defecto
        method = self._select_method(data_type)
        logger.info(f"[{gse_id}] Método DEA: {method} | Data type: {data_type}")

        # Ejecutar DEA — cadena: edgeR → limma → t-test
        if method == "edger" and RPY2_AVAILABLE:
            results_df = self._run_edger(matrix, sample_meta, gse_id)
        elif method == "deseq2" and RPY2_AVAILABLE:
            results_df = self._run_deseq2(matrix, sample_meta, gse_id)
        elif method == "limma" and RPY2_AVAILABLE:
            results_df = self._run_limma(matrix, sample_meta, gse_id)
        else:
            logger.warning(f"[{gse_id}] Usando fallback Python (t-test)")
            method = "ttest_python"
            results_df = self._run_ttest_fallback(matrix, sample_meta, gse_id)

        # Filtrar genes significativos
        sig_genes = results_df[
            (results_df["padj"] < self.padj_threshold) &
            (results_df["log2FC"].abs() > self.lfc_threshold)
        ].copy()

        sig_genes["direction"] = np.where(sig_genes["log2FC"] > 0, "UP", "DOWN")

        # Preparar output de genes significativos
        sig_list = sig_genes[["gene", "log2FC", "pvalue", "padj", "direction"]].to_dict(orient="records")

        summary = {
            "gse_id": gse_id,
            "method": method,
            "n_total_genes": len(results_df),
            "n_significant": len(sig_genes),
            "n_up": int((sig_genes["direction"] == "UP").sum()),
            "n_down": int((sig_genes["direction"] == "DOWN").sum()),
            "padj_threshold": self.padj_threshold,
            "lfc_threshold": self.lfc_threshold,
        }

        # Guardar
        results_df.to_csv(out_dir / "dea_results.csv", index=False)
        results_df[["gene", "log2FC", "pvalue", "padj"]].to_csv(out_dir / "volcano_data.csv", index=False)

        with open(out_dir / "significant_genes.json", "w") as f:
            json.dump({"gse_id": gse_id, "genes": sig_list}, f, indent=2, default=str)

        with open(out_dir / "dea_summary.json", "w") as f:
            json.dump(summary, f, indent=2, default=str)

        logger.info(f"[{gse_id}] DEA completo | Significativos: {summary['n_significant']} (UP={summary['n_up']}, DOWN={summary['n_down']})")
        return summary

    # ------------------------------------------------------------------
    # edgeR via rpy2  (método principal)
    # ------------------------------------------------------------------

    def _run_edger(self, matrix: pd.DataFrame, sample_meta: pd.DataFrame, gse_id: str) -> pd.DataFrame:
        """
        Ejecuta edgeR via rpy2.
        Usa exactTest para diseño simple (2 grupos) y glmQLFTest si hay covariables.
        """
        try:
            edgeR = importr("edgeR")
            base = importr("base")

            # edgeR requiere counts crudos (enteros)
            counts_int = matrix.round().astype(int).clip(lower=0)

            groups = sample_meta["group"].tolist()
            ordered_groups = self._order_groups(sample_meta["group"].unique().tolist())
            group_factor = ro.FactorVector(groups, levels=ro.StrVector(ordered_groups))

            # Crear DGEList
            r_counts = pandas2ri.py2rpy(counts_int)
            dge = edgeR.DGEList(counts=r_counts, group=group_factor)

            # Filtrar genes de baja expresión con filterByExpr
            keep = edgeR.filterByExpr(dge)
            dge = base.subset(dge, keep, )

            # Normalización TMM
            dge = edgeR.calcNormFactors(dge, method="TMM")

            # Estimar dispersión
            dge = edgeR.estimateDisp(dge, robust=True)

            # exactTest (2 grupos)
            et = edgeR.exactTest(dge, pair=ro.StrVector(ordered_groups))
            top_table = edgeR.topTags(et, n=matrix.shape[0], sort_by="PValue")
            top_df = pandas2ri.rpy2py(base.as_data_frame(top_table.rx2("table"))).reset_index()

            # Renombrar columnas al estándar interno
            top_df.columns = ["gene"] + list(top_df.columns[1:])
            col_map = {"logFC": "log2FC", "PValue": "pvalue", "FDR": "padj"}
            top_df = top_df.rename(columns=col_map)

            return top_df[["gene", "log2FC", "pvalue", "padj"]].dropna()

        except Exception as e:
            logger.warning(f"[{gse_id}] edgeR falló ({e}), intentando limma")
            return self._run_limma(matrix, sample_meta, gse_id)



    def _run_deseq2(self, matrix: pd.DataFrame, sample_meta: pd.DataFrame, gse_id: str) -> pd.DataFrame:
        """Ejecuta DESeq2 via rpy2. Requiere counts enteros sin normalizar."""
        try:
            deseq2 = importr("DESeq2")
            base = importr("base")

            # DESeq2 requiere counts crudos (enteros)
            counts_int = matrix.round().astype(int)
            counts_int = counts_int.clip(lower=0)

            # Crear colData
            col_data = ro.DataFrame({
                "condition": ro.StrVector(sample_meta["group"].tolist())
            })
            col_data.rownames = ro.StrVector(sample_meta.index.tolist())

            # Crear DESeqDataSet
            r_counts = pandas2ri.py2rpy(counts_int)
            dds = deseq2.DESeqDataSetFromMatrix(
                countData=r_counts,
                colData=col_data,
                design=Formula("~ condition")
            )

            # Correr DESeq
            dds = deseq2.DESeq(dds, quiet=True)

            # Extraer resultados (case vs control)
            groups = sample_meta["group"].unique().tolist()
            contrast = ro.StrVector(["condition"] + self._order_groups(groups))
            res = deseq2.results(dds, contrast=contrast)
            res_df = pandas2ri.rpy2py(base.as_data_frame(res))

            res_df = res_df.reset_index()
            res_df.columns = ["gene", "baseMean", "log2FC", "lfcSE", "stat", "pvalue", "padj"]
            res_df = res_df.dropna(subset=["pvalue", "padj"])

            return res_df[["gene", "log2FC", "pvalue", "padj"]]

        except Exception as e:
            logger.warning(f"[{gse_id}] DESeq2 falló ({e}), intentando limma")
            return self._run_limma(matrix, sample_meta, gse_id)

    # ------------------------------------------------------------------
    # limma via rpy2
    # ------------------------------------------------------------------

    def _run_limma(self, matrix: pd.DataFrame, sample_meta: pd.DataFrame, gse_id: str) -> pd.DataFrame:
        """Ejecuta limma-voom via rpy2. Funciona con datos ya normalizados."""
        try:
            limma = importr("limma")
            base = importr("base")

            groups = sample_meta["group"].tolist()
            ordered_groups = self._order_groups(sample_meta["group"].unique().tolist())
            group_factor = ro.FactorVector(groups, levels=ro.StrVector(ordered_groups))

            design = limma.model_matrix(Formula("~ group_factor"),
                                        data=ro.DataFrame({"group_factor": group_factor}))

            r_matrix = pandas2ri.py2rpy(matrix.astype(float))
            fit = limma.lmFit(r_matrix, design)
            fit = limma.eBayes(fit)

            top_table = limma.topTable(fit, coef=2, number=matrix.shape[0], sort_by="P", adjust_method="BH")
            top_df = pandas2ri.rpy2py(top_table).reset_index()
            top_df.columns = ["gene"] + list(top_df.columns[1:])

            # Renombrar columnas estándar
            col_map = {
                "logFC": "log2FC",
                "P.Value": "pvalue",
                "adj.P.Val": "padj",
            }
            top_df = top_df.rename(columns=col_map)
            return top_df[["gene", "log2FC", "pvalue", "padj"]].dropna()

        except Exception as e:
            logger.warning(f"[{gse_id}] limma falló ({e}), usando fallback Python")
            return self._run_ttest_fallback(matrix, sample_meta, gse_id)

    # ------------------------------------------------------------------
    # Fallback Python: t-test + BH correction
    # ------------------------------------------------------------------

    def _run_ttest_fallback(self, matrix: pd.DataFrame, sample_meta: pd.DataFrame, gse_id: str) -> pd.DataFrame:
        """T-test de Welch + corrección BH. Fallback sin dependencia R."""
        from scipy import stats
        from statsmodels.stats.multitest import multipletests

        groups = sample_meta["group"].unique()
        ordered = self._order_groups(list(groups))
        case_samples = sample_meta[sample_meta["group"] == ordered[0]].index
        ctrl_samples = sample_meta[sample_meta["group"] == ordered[1]].index

        case_mat = matrix[matrix.columns.intersection(case_samples)]
        ctrl_mat = matrix[matrix.columns.intersection(ctrl_samples)]

        results = []
        for gene in matrix.index:
            case_vals = case_mat.loc[gene].dropna().values
            ctrl_vals = ctrl_mat.loc[gene].dropna().values
            if len(case_vals) < 2 or len(ctrl_vals) < 2:
                continue
            t_stat, pval = stats.ttest_ind(case_vals, ctrl_vals, equal_var=False)
            lfc = np.mean(case_vals) - np.mean(ctrl_vals)
            results.append({"gene": gene, "log2FC": lfc, "pvalue": pval})

        res_df = pd.DataFrame(results).dropna()
        if res_df.empty:
            raise ValueError("No se pudieron calcular resultados DEA")

        _, padj, _, _ = multipletests(res_df["pvalue"].fillna(1), method="fdr_bh")
        res_df["padj"] = padj

        return res_df[["gene", "log2FC", "pvalue", "padj"]]

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _select_method(self, data_type: str) -> str:
        """Selecciona el método DEA. edgeR es el default para ambos tipos de dato."""
        if self.prefer_method != "auto":
            return self.prefer_method
        # edgeR funciona bien tanto para RNA-seq como microarray
        return "edger"

    def _order_groups(self, groups: list) -> list:
        """Ordena grupos: case primero, control segundo."""
        case_keywords = ["case", "disease", "tumor", "cancer", "patient", "treated"]
        ctrl_keywords = ["control", "ctrl", "normal", "healthy", "untreated"]

        cases = [g for g in groups if any(k in g.lower() for k in case_keywords)]
        ctrls = [g for g in groups if any(k in g.lower() for k in ctrl_keywords)]

        if cases and ctrls:
            return [cases[0], ctrls[0]]
        return groups[:2]  # fallback: usar el orden tal como viene

    def _load_sample_metadata(self, dataset_dir: Path) -> Optional[pd.DataFrame]:
        for fname in ["sample_metadata.csv", "samples.csv", "phenodata.csv"]:
            fpath = dataset_dir / fname
            if fpath.exists():
                return pd.read_csv(fpath, index_col=0)
        return None

    def _load_qc_report(self, dataset_dir: Path) -> dict:
        qc_path = dataset_dir / "qc_report.json"
        if qc_path.exists():
            with open(qc_path) as f:
                return json.load(f)
        return {}


# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------

@click.command()
@click.option("--input", "processed_dir", required=True, help="Directorio de datos procesados (data/processed/)")
@click.option("--output", "output_dir", default="data/dea/", help="Directorio de salida")
@click.option("--method", default="edger", type=click.Choice(["edger", "deseq2", "limma"]), help="Método DEA")
@click.option("--padj", default=0.05, help="Umbral FDR ajustado")
@click.option("--lfc", default=1.0, help="Umbral log2 Fold Change")
def main(processed_dir, output_dir, method, padj, lfc):
    """Agent 4: Análisis de expresión diferencial (DESeq2/limma/fallback)."""
    config = {"method": method, "padj_threshold": padj, "lfc_threshold": lfc}
    agent = DifferentialExpressionAgent(config=config)
    summary = agent.run(processed_dir, output_dir)
    click.echo(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()