"""
BioSignal Discovery Engine
Agent 3: PreprocessingAgent
============================
Responsabilidad: Normalizar, filtrar y realizar QC de cada dataset descargado.
Soporta RNA-seq (counts) y microarray (intensidades).

Uso:
    python -m agents.preprocess --input data/raw/ --output data/processed/
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


class PreprocessingAgent:
    """
    Agente de preprocesamiento y control de calidad de datasets GEO.

    Pipeline por tipo de dato:
        RNA-seq:
            1. Filtrar genes con expresión baja (< min_counts en < frac_samples)
            2. Normalización TMM-like via log2(CPM + 1)
            3. Detección de outliers (IQR sobre PCA componente 1)
            4. Export: matrix_normalized.csv + qc_report.json

        Microarray:
            1. Verificar normalización previa (quantile-normalized en GEO)
            2. Log2-transform si valores > 100
            3. Filtrar probes de baja varianza (25th percentile)
            4. Colapsar probes múltiples por gen (mean)
            5. Export: matrix_normalized.csv + qc_report.json

    Salida por dataset:
        data/processed/{gse_id}/matrix_normalized.csv
        data/processed/{gse_id}/sample_metadata.csv
        data/processed/{gse_id}/qc_report.json
        data/processed/{gse_id}/qc_status.txt   # PASS | FAIL
    """

    def __init__(self, config: dict = None):
        self.config = config or {}
        self.min_samples_per_group = self.config.get("min_samples_per_group", 3)
        self.min_genes_detected = self.config.get("min_genes_detected", 5000)
        self.max_outlier_fraction = self.config.get("max_outlier_fraction", 0.20)
        self.low_expr_threshold = self.config.get("low_expression_threshold", 10)
        self.low_expr_fraction = self.config.get("low_expression_fraction", 0.70)
        self.variance_percentile = self.config.get("variance_percentile", 25)

    # ------------------------------------------------------------------
    # Entrada principal
    # ------------------------------------------------------------------

    def run(self, raw_dir: str, output_dir: str) -> dict:
        """
        Ejecuta el preprocesamiento de todos los datasets descargados.

        Args:
            raw_dir:    Directorio de salida del Agent 2 (data/raw/)
            output_dir: Directorio destino (data/processed/)

        Returns:
            Resumen {"processed": [...], "failed": [...], "skipped": [...]}
        """
        raw_path = Path(raw_dir)
        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        summary = {"processed": [], "failed": [], "skipped": []}
        start = time.time()

        dataset_dirs = sorted([d for d in raw_path.iterdir() if d.is_dir()])
        logger.info(f"[PreprocessingAgent] Encontrados {len(dataset_dirs)} datasets en {raw_dir}")

        for dataset_dir in dataset_dirs:
            gse_id = dataset_dir.name
            try:
                result = self._process_dataset(dataset_dir, out_path / gse_id)
                if result["status"] == "PASS":
                    summary["processed"].append(gse_id)
                else:
                    summary["skipped"].append(gse_id)
            except Exception as e:
                logger.error(f"[{gse_id}] Error en preprocesamiento: {e}")
                summary["failed"].append({"gse_id": gse_id, "error": str(e)})

        elapsed = time.time() - start
        logger.info(
            f"[PreprocessingAgent] Completado en {elapsed:.1f}s | "
            f"OK={len(summary['processed'])} | "
            f"SKIP={len(summary['skipped'])} | "
            f"FAIL={len(summary['failed'])}"
        )
        return summary

    # ------------------------------------------------------------------
    # Procesamiento por dataset
    # ------------------------------------------------------------------

    def _process_dataset(self, dataset_dir: Path, out_dir: Path) -> dict:
        """Procesa un dataset individual."""
        out_dir.mkdir(parents=True, exist_ok=True)
        gse_id = dataset_dir.name

        # Cargar metadata y detectar tipo de dato
        metadata = self._load_metadata(dataset_dir)
        data_type = metadata.get("data_type", "unknown")

        logger.info(f"[{gse_id}] Tipo: {data_type}")

        # Cargar matriz de expresión
        matrix = self._load_expression_matrix(dataset_dir)
        if matrix is None or matrix.empty:
            return self._fail(out_dir, gse_id, "Matriz de expresión vacía o no encontrada")

        logger.debug(f"[{gse_id}] Dimensiones iniciales: {matrix.shape}")

        # Pipeline según tipo de dato
        if data_type == "RNA-seq":
            matrix, qc = self._process_rnaseq(matrix, gse_id)
        else:
            matrix, qc = self._process_microarray(matrix, gse_id)

        # Validar muestras mínimas por grupo
        sample_meta = self._load_sample_metadata(dataset_dir)
        group_check = self._check_group_balance(sample_meta, gse_id)
        qc.update(group_check)

        # Detección de outliers
        matrix, outliers = self._detect_outliers(matrix, gse_id)
        qc["outliers_removed"] = outliers
        qc["outlier_fraction"] = len(outliers) / max(matrix.shape[1] + len(outliers), 1)

        # Validación de genes mínimos
        qc["genes_final"] = matrix.shape[0]
        qc["samples_final"] = matrix.shape[1]

        status = self._compute_qc_status(qc, gse_id)
        qc["status"] = status

        # Guardar outputs
        matrix.to_csv(out_dir / "matrix_normalized.csv")
        if sample_meta is not None:
            sample_meta.to_csv(out_dir / "sample_metadata.csv")

        with open(out_dir / "qc_report.json", "w") as f:
            json.dump(qc, f, indent=2, default=str)

        (out_dir / "qc_status.txt").write_text(status)
        logger.info(f"[{gse_id}] QC Status: {status} | Genes: {qc['genes_final']} | Muestras: {qc['samples_final']}")

        return {"status": status, "qc": qc}

    # ------------------------------------------------------------------
    # Pipeline RNA-seq
    # ------------------------------------------------------------------

    def _process_rnaseq(self, matrix: pd.DataFrame, gse_id: str) -> tuple[pd.DataFrame, dict]:
        """Normalización y filtrado para datos de RNA-seq (counts)."""
        qc = {"data_type": "RNA-seq"}
        qc["genes_raw"] = matrix.shape[0]
        qc["samples_raw"] = matrix.shape[1]

        # 1. Filtrar genes de baja expresión
        min_samples = max(1, int(matrix.shape[1] * self.low_expr_fraction))
        keep = (matrix >= self.low_expr_threshold).sum(axis=1) >= min_samples
        matrix = matrix.loc[keep]
        qc["genes_after_lowexpr_filter"] = matrix.shape[0]
        logger.debug(f"[{gse_id}] Genes tras filtro baja expresión: {matrix.shape[0]}")

        # 2. Normalización log2(CPM + 1)
        col_sums = matrix.sum(axis=0)
        cpm = matrix.div(col_sums, axis=1) * 1e6
        matrix_norm = np.log2(cpm + 1)
        qc["normalization"] = "log2(CPM+1)"

        return matrix_norm, qc

    # ------------------------------------------------------------------
    # Pipeline Microarray
    # ------------------------------------------------------------------

    def _process_microarray(self, matrix: pd.DataFrame, gse_id: str) -> tuple[pd.DataFrame, dict]:
        """Filtrado y colapso de probes para datos de microarray."""
        qc = {"data_type": "microarray"}
        qc["genes_raw"] = matrix.shape[0]
        qc["samples_raw"] = matrix.shape[1]

        # 1. Log2-transform si no está normalizado
        max_val = matrix.max().max()
        if max_val > 100:
            matrix = np.log2(matrix.clip(lower=1))
            qc["log2_transformed"] = True
            logger.debug(f"[{gse_id}] Aplicado log2 transform (max_val={max_val:.1f})")
        else:
            qc["log2_transformed"] = False

        # 2. Filtrar probes de baja varianza
        variances = matrix.var(axis=1)
        var_threshold = np.percentile(variances, self.variance_percentile)
        matrix = matrix.loc[variances >= var_threshold]
        qc["probes_after_variance_filter"] = matrix.shape[0]

        # 3. Colapsar probes múltiples por gen (tomar la media)
        matrix.index = matrix.index.astype(str).str.split("_").str[0]
        matrix = matrix.groupby(matrix.index).mean()
        qc["genes_after_probe_collapse"] = matrix.shape[0]
        qc["normalization"] = "quantile (GEO) + log2 if needed"

        return matrix, qc

    # ------------------------------------------------------------------
    # Detección de outliers (PCA-based IQR)
    # ------------------------------------------------------------------

    def _detect_outliers(self, matrix: pd.DataFrame, gse_id: str) -> tuple[pd.DataFrame, list]:
        """Detecta muestras outliers usando correlación entre muestras."""
        if matrix.shape[1] < 4:
            return matrix, []

        try:
            corr_matrix = matrix.corr()
            mean_corr = corr_matrix.mean()
            q1 = mean_corr.quantile(0.25)
            q3 = mean_corr.quantile(0.75)
            iqr = q3 - q1
            lower_bound = q1 - 1.5 * iqr

            outlier_samples = mean_corr[mean_corr < lower_bound].index.tolist()

            if len(outlier_samples) / matrix.shape[1] <= self.max_outlier_fraction:
                matrix = matrix.drop(columns=outlier_samples)
                if outlier_samples:
                    logger.info(f"[{gse_id}] Outliers removidos: {outlier_samples}")
                return matrix, outlier_samples
            else:
                logger.warning(f"[{gse_id}] Demasiados outliers ({len(outlier_samples)}), no se eliminan")
                return matrix, []
        except Exception as e:
            logger.warning(f"[{gse_id}] Error en detección de outliers: {e}")
            return matrix, []

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _load_expression_matrix(self, dataset_dir: Path) -> Optional[pd.DataFrame]:
        """Carga la matriz de expresión desde archivos GEO."""
        candidates = [
            ("matrix.tsv", "\t"),
            ("matrix.csv", ","),
            ("expression_matrix.csv", ","),
            ("expression_matrix.tsv", "\t"),
        ]
        for fname, sep in candidates:
            fpath = dataset_dir / fname
            if not fpath.exists():
                continue
            try:
                # Detectar si es GEO Series Matrix (tiene cabeceras con !)
                with open(fpath, "r", encoding="utf-8", errors="replace") as f:
                    first_line = f.readline()

                if first_line.startswith("!"):
                    df = self._parse_geo_series_matrix(fpath, sep)
                else:
                    df = pd.read_csv(fpath, index_col=0, sep=sep)

                if df is not None and not df.empty:
                    df = df.apply(pd.to_numeric, errors="coerce").dropna(how="all")
                    return df

            except Exception as e:
                logger.warning(f"Error cargando {fpath}: {e}")
        return None

    def _parse_geo_series_matrix(self, fpath: Path, sep: str = "\t") -> Optional[pd.DataFrame]:
        """
        Parsea un GEO Series Matrix file con cabeceras que empiezan con '!'.
        Los datos reales están entre:
            !series_matrix_table_begin
            ...datos...
            !series_matrix_table_end
        Si no existe esa sección, el archivo es solo metadata (RNA-seq SRA)
        y retorna None.
        """
        in_table = False
        lines = []

        with open(fpath, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line_stripped = line.strip()
                if line_stripped.lower() == "!series_matrix_table_begin":
                    in_table = True
                    continue
                if line_stripped.lower() == "!series_matrix_table_end":
                    break
                if in_table:
                    lines.append(line_stripped)

        if not lines:
            logger.warning(f"[preprocess] {fpath.name}: sin tabla de expresión (posiblemente RNA-seq SRA sin counts en series matrix)")
            return None

        from io import StringIO
        content = "\n".join(lines)
        df = pd.read_csv(StringIO(content), sep=sep, index_col=0)

        # Limpiar nombre del índice (suele ser "ID_REF")
        df.index.name = "gene"
        logger.debug(f"[preprocess] GEO matrix parseado: {df.shape[0]} genes × {df.shape[1]} muestras")
        return df


    def _load_metadata(self, dataset_dir: Path) -> dict:
        """Carga metadata del dataset e infiere data_type si no está explícito."""
        meta_path = dataset_dir / "metadata.json"
        if not meta_path.exists():
            return {}
        with open(meta_path) as f:
            meta = json.load(f)

        # Inferir data_type desde molecule o type de muestras
        if "data_type" not in meta:
            samples = meta.get("samples", [])
            if samples:
                molecule = samples[0].get("characteristics", {}).get("molecule_ch1", "")
                sample_type = samples[0].get("characteristics", {}).get("type", "")
                if "RNA" in molecule or "SRA" in sample_type:
                    meta["data_type"] = "RNA-seq"
                elif "genomic" in molecule.lower():
                    meta["data_type"] = "microarray"
                else:
                    meta["data_type"] = "RNA-seq"  # default seguro

        return meta

    def _load_sample_metadata(self, dataset_dir: Path) -> Optional[pd.DataFrame]:
        """Carga metadata de muestras. Lee desde metadata.json si no hay CSV separado."""
        # Primero buscar CSV dedicado
        for fname in ["sample_metadata.csv", "samples.csv", "phenodata.csv"]:
            fpath = dataset_dir / fname
            if fpath.exists():
                return pd.read_csv(fpath, index_col=0)

        # Construir desde metadata.json (formato del Agent 2)
        meta_path = dataset_dir / "metadata.json"
        if meta_path.exists():
            with open(meta_path) as f:
                meta = json.load(f)
            samples = meta.get("samples", [])
            if samples:
                rows = []
                for s in samples:
                    rows.append({
                        "sample_id": s.get("gsm_id", ""),
                        "title": s.get("title", ""),
                        "group": s.get("label", "unknown"),  # 'case' | 'control'
                    })
                df = pd.DataFrame(rows).set_index("sample_id")
                return df

        return None

    def _check_group_balance(self, sample_meta: Optional[pd.DataFrame], gse_id: str) -> dict:
        """Verifica que haya suficientes muestras por grupo."""
        if sample_meta is None or "group" not in sample_meta.columns:
            return {"group_check": "skipped_no_metadata"}

        counts = sample_meta["group"].value_counts()
        min_group = counts.min()
        result = {
            "group_counts": counts.to_dict(),
            "min_group_samples": int(min_group),
            "group_balance_ok": bool(min_group >= self.min_samples_per_group),
        }
        if not result["group_balance_ok"]:
            logger.warning(f"[{gse_id}] Grupo con < {self.min_samples_per_group} muestras: {counts.to_dict()}")
        return result

    def _compute_qc_status(self, qc: dict, gse_id: str) -> str:
        """Determina si el dataset pasa o falla el QC."""
        reasons = []

        if qc.get("genes_final", 0) < self.min_genes_detected:
            reasons.append(f"genes_final={qc.get('genes_final')} < {self.min_genes_detected}")

        if not qc.get("group_balance_ok", True):
            reasons.append("grupo desbalanceado")

        if qc.get("outlier_fraction", 0) > self.max_outlier_fraction:
            reasons.append(f"outlier_fraction={qc.get('outlier_fraction'):.2f} > {self.max_outlier_fraction}")

        if reasons:
            logger.warning(f"[{gse_id}] QC FAIL: {'; '.join(reasons)}")
            return "FAIL"
        return "PASS"

    def _fail(self, out_dir: Path, gse_id: str, reason: str) -> dict:
        """Registra un fallo y retorna status FAIL."""
        logger.error(f"[{gse_id}] FAIL: {reason}")
        qc = {"status": "FAIL", "error": reason}
        with open(out_dir / "qc_report.json", "w") as f:
            json.dump(qc, f, indent=2)
        (out_dir / "qc_status.txt").write_text("FAIL")
        return {"status": "FAIL", "qc": qc}


# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------

@click.command()
@click.option("--input", "raw_dir", required=True, help="Directorio con datasets crudos (data/raw/)")
@click.option("--output", "output_dir", default="data/processed/", help="Directorio de salida")
@click.option("--min-genes", default=5000, help="Mínimo genes detectados para PASS")
@click.option("--min-samples", default=3, help="Mínimo muestras por grupo")
def main(raw_dir, output_dir, min_genes, min_samples):
    """Agent 3: Preprocesa y realiza QC de datasets GEO."""
    config = {
        "min_genes_detected": min_genes,
        "min_samples_per_group": min_samples,
    }
    agent = PreprocessingAgent(config=config)
    summary = agent.run(raw_dir, output_dir)
    click.echo(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
