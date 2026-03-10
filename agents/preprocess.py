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
            1. Filtrar genes con expresiÃ³n baja (< min_counts en < frac_samples)
            2. NormalizaciÃ³n TMM-like via log2(CPM + 1)
            3. DetecciÃ³n de outliers (IQR sobre PCA componente 1)
            4. Export: matrix_normalized.csv + qc_report.json

        Microarray:
            1. Verificar normalizaciÃ³n previa (quantile-normalized en GEO)
            2. Log2-transform si valores > 100
            3. Filtrar probes de baja varianza (25th percentile)
            4. Colapsar probes mÃºltiples por gen (mean)
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
        self.low_expr_fraction_large = self.config.get("low_expression_fraction_large", 0.20)
        self.large_dataset_threshold = self.config.get("large_dataset_threshold", 100)
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

        # Cargar matriz de expresiÃ³n
        matrix = self._load_expression_matrix(dataset_dir)
        if matrix is None or matrix.empty:
            return self._fail(out_dir, gse_id, "Matriz de expresiÃ³n vacÃ­a o no encontrada")

        logger.debug(f"[{gse_id}] Dimensiones iniciales: {matrix.shape}")

        # Pipeline segÃºn tipo de dato
        if data_type == "RNA-seq":
            matrix, qc, counts_raw = self._process_rnaseq(matrix, gse_id)
        else:
            matrix, qc, counts_raw = self._process_microarray(matrix, gse_id)

        # Validar muestras mÃ­nimas por grupo
        sample_meta = self._load_sample_metadata(dataset_dir, matrix_cols=matrix.columns.tolist(), out_dir=out_dir)
        group_check = self._check_group_balance(sample_meta, gse_id)
        qc.update(group_check)

        # DetecciÃ³n de outliers
        matrix, outliers = self._detect_outliers(matrix, gse_id)
        qc["outliers_removed"] = outliers
        qc["outlier_fraction"] = len(outliers) / max(matrix.shape[1] + len(outliers), 1)

        # ValidaciÃ³n de genes mÃ­nimos
        qc["genes_final"] = matrix.shape[0]
        qc["samples_final"] = matrix.shape[1]

        status = self._compute_qc_status(qc, gse_id)
        qc["status"] = status

        # Guardar outputs
        matrix.to_csv(out_dir / "matrix_normalized.csv")
        if counts_raw is not None:
            counts_raw.to_csv(out_dir / "matrix_counts.csv")
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

    def _process_rnaseq(self, matrix: pd.DataFrame, gse_id: str) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
        """NormalizaciÃ³n y filtrado para datos de RNA-seq (counts)."""
        qc = {"data_type": "RNA-seq"}
        qc["genes_raw"] = matrix.shape[0]
        qc["samples_raw"] = matrix.shape[1]

        # 1. Filtrar genes de baja expresiÃ³n
        p99_val = matrix.stack().quantile(0.99)
        min_val = matrix.min().min()
        already_normalized = p99_val < 50 and min_val > -5
        if already_normalized:
            variances = matrix.var(axis=1)
            keep = variances >= variances.quantile(0.1)
            matrix = matrix.loc[keep]
            logger.debug(f"[{gse_id}] Datos pre-normalizados, filtro varianza aplicado")
        else:
            # Para datasets grandes (>100 muestras) bajar fraccion minima
            n_samples = matrix.shape[1]
            effective_fraction = self.low_expr_fraction if n_samples <= self.large_dataset_threshold else self.low_expr_fraction_large
            min_samples = max(1, int(n_samples * effective_fraction))
            keep = (matrix >= self.low_expr_threshold).sum(axis=1) >= min_samples
            matrix = matrix.loc[keep]
        qc["genes_after_lowexpr_filter"] = matrix.shape[0]
        logger.debug(f"[{gse_id}] Genes tras filtro baja expresiÃ³n: {matrix.shape[0]}")

        # Guardar counts filtrados (antes de normalizar) para edgeR
        counts_filtered = matrix.copy()

        # 2. NormalizaciÃ³n log2(CPM + 1)
        col_sums = matrix.sum(axis=0)
        cpm = matrix.div(col_sums, axis=1) * 1e6
        matrix_norm = np.log2(cpm + 1)
        qc["normalization"] = "log2(CPM+1)"

        return matrix_norm, qc, counts_filtered

    # ------------------------------------------------------------------
    # Pipeline Microarray
    # ------------------------------------------------------------------

    def _process_microarray(self, matrix: pd.DataFrame, gse_id: str) -> tuple[pd.DataFrame, dict]:
        """Filtrado y colapso de probes para datos de microarray."""
        qc = {"data_type": "microarray"}
        qc["genes_raw"] = matrix.shape[0]
        qc["samples_raw"] = matrix.shape[1]

        # 1. Log2-transform si no estÃ¡ normalizado
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

        # 3. Colapsar probes mÃºltiples por gen (tomar la media)
        matrix.index = matrix.index.astype(str).str.split("_").str[0]
        matrix = matrix.groupby(matrix.index).mean()
        qc["genes_after_probe_collapse"] = matrix.shape[0]
        qc["normalization"] = "quantile (GEO) + log2 if needed"

        return matrix, qc, None

    # ------------------------------------------------------------------
    # DetecciÃ³n de outliers (PCA-based IQR)
    # ------------------------------------------------------------------

    def _detect_outliers(self, matrix: pd.DataFrame, gse_id: str) -> tuple[pd.DataFrame, list]:
        """Detecta muestras outliers usando correlaciÃ³n entre muestras."""
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
            logger.warning(f"[{gse_id}] Error en detecciÃ³n de outliers: {e}")
            return matrix, []

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _load_expression_matrix(self, dataset_dir: Path) -> Optional[pd.DataFrame]:
        """Carga la matriz de expresiÃ³n desde archivos GEO."""
        candidates = [
            ("matrix_counts.tsv", "\t"),  # Counts suplementarios consolidados (SRA)
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
        Los datos reales estÃ¡n entre:
            !series_matrix_table_begin
            ...datos...
            !series_matrix_table_end
        Si no existe esa secciÃ³n, el archivo es solo metadata (RNA-seq SRA)
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
            logger.warning(f"[preprocess] {fpath.name}: sin tabla de expresiÃ³n (posiblemente RNA-seq SRA sin counts en series matrix)")
            return None

        from io import StringIO
        content = "\n".join(lines)
        df = pd.read_csv(StringIO(content), sep=sep, index_col=0)

        # Limpiar nombre del Ã­ndice (suele ser "ID_REF")
        df.index.name = "gene"
        logger.debug(f"[preprocess] GEO matrix parseado: {df.shape[0]} genes Ã— {df.shape[1]} muestras")
        return df


    def _load_metadata(self, dataset_dir: Path) -> dict:
        """Carga metadata del dataset e infiere data_type si no estÃ¡ explÃ­cito."""
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

    def _load_sample_metadata(self, dataset_dir: Path, matrix_cols: list = None, out_dir: Path = None) -> Optional[pd.DataFrame]:
        """
        Carga metadata de muestras. Intenta alinear el Ã­ndice con las columnas
        de la matriz (que pueden ser tÃ­tulos en lugar de GSM IDs).
        """
        search_dirs = [d for d in [out_dir, dataset_dir] if d is not None]
        for fname in ["sample_metadata.csv", "samples.csv", "phenodata.csv"]:
            for search_dir in search_dirs:
                fpath = search_dir / fname
                if fpath.exists():
                    df_cached = pd.read_csv(fpath, index_col=0)
                    if matrix_cols is None or df_cached.index.isin(matrix_cols).any():
                        return df_cached
                    import os; os.remove(fpath)

        meta_path = dataset_dir / "metadata.json"
        if not meta_path.exists():
            return None

        with open(meta_path) as f:
            meta = json.load(f)

        all_samples = meta.get("samples", []) + meta.get("unclassified", [])
        if not all_samples:
            return None

        rows = []
        for s in all_samples:
            rows.append({
                "gsm_id": s.get("gsm_id", ""),
                "title": s.get("title", ""),
                "description": s.get("characteristics", {}).get("description", ""),
                "group": s.get("label", "unclassified"),
            })
        df = pd.DataFrame(rows)

        if matrix_cols:
            # Caso 0: description (ej. WC1, KC1)
            df["desc_clean"] = df["description"].str.strip()
            if df["desc_clean"].isin(matrix_cols).any():
                return df.set_index("desc_clean")[["group"]]
            # Caso 1: columnas son GSM IDs
            if df["gsm_id"].isin(matrix_cols).any():
                return df.set_index("gsm_id")[["group"]]
            # Caso 2: columnas son tÃ­tulos â€” limpiar prefijos tipo "A2712: "
            df["title_clean"] = df["title"].str.split(": ").str[-1].str.strip()
            if df["title_clean"].isin(matrix_cols).any():
                return df.set_index("title_clean")[["group"]]
            # Caso 3: match por substring
            title_map = {}
            for col in matrix_cols:
                for _, row in df.iterrows():
                    tc = str(row["title_clean"])
                    if col in row["title"] or tc in col or col in tc:
                        title_map[col] = row["group"]
                        break
            if title_map:
                return pd.DataFrame({"group": title_map})
            # Caso 4: Capa 2 LLM
            logger.info(f"Intentando Capa 2 (LLM) para {dataset_dir.name}...")
            llm_map = self._llm_map_columns(matrix_cols, all_samples, gse_id=dataset_dir.name)
            llm_groups = set(llm_map.values()) if llm_map else set()
            llm_valid = llm_map and len(llm_groups) >= 2
            if llm_valid:
                result_df = pd.DataFrame({"group": llm_map})
                result_df.to_csv(dataset_dir / "sample_metadata.csv")
                logger.success(f"Capa 2 LLM: {len(llm_map)} columnas mapeadas")
                return result_df
            # Caso 5: Capa 3 SRA RunInfo
            logger.info(f"Capa 2 inválida (grupos={llm_groups}), intentando Capa 3 (SRA)...")
            sra_map = self._sra_map_columns(matrix_cols, all_samples, gse_id=dataset_dir.name)
            if sra_map and len(set(sra_map.values())) >= 2:
                result_df = pd.DataFrame({"group": sra_map})
                result_df.to_csv(dataset_dir / "sample_metadata.csv")
                logger.success(f"Capa 3 SRA: {len(sra_map)} columnas mapeadas")
                return result_df

        return df.set_index("gsm_id")[["group"]]


    def _llm_map_columns(self, matrix_cols: list, all_samples: list, gse_id: str) -> dict:
        """
        Capa 2: Usa Ollama para inferir mapeo columnas -> grupo cuando los casos
        determinísticos fallan.

        Args:
            matrix_cols: Columnas de la matriz (ej. ['GQ01_S46', 'GQ02_S68', ...])
            all_samples: Lista de dicts con gsm_id, title, label de metadata.json
            gse_id: ID del dataset para logging

        Returns:
            dict {columna: 'case'|'control'} o {} si falla
        """
        try:
            import requests as _requests
            import json as _json

            # Preparar contexto para el LLM
            sample_info = []
            for s in all_samples[:30]:  # Limitar contexto
                sample_info.append({
                    "title": s.get("title", ""),
                    "label": s.get("label", "unclassified"),
                })

            cols_preview = matrix_cols[:30]

            prompt = f"""You are a bioinformatics expert. Map matrix column names to sample groups.

Matrix columns (first {len(cols_preview)}): {cols_preview}

Sample metadata:
{_json.dumps(sample_info, indent=2)}

Task: For each matrix column, determine if it belongs to 'case' or 'control' group.
Use the sample titles, GSM IDs, and labels to infer the mapping.
Look for patterns: control=WT/normal/ctrl/vehicle, case=KO/treated/disease/AD/mut.

Respond ONLY with a JSON object mapping column names to groups. Example:
{{"GQ01_S46": "control", "GQ02_S68": "case"}}

Include ALL {len(matrix_cols)} columns. No explanation, just JSON."""

            payload = {
                "model": "gemma3:4b",
                "prompt": prompt,
                "stream": False,
                "options": {"temperature": 0.1}
            }
            resp = _requests.post(
                "http://localhost:11434/api/generate",
                json=payload,
                timeout=60
            )
            if resp.status_code != 200:
                logger.warning(f"[{gse_id}] Ollama error: {resp.status_code}")
                return {}

            text = resp.json().get("response", "").strip()

            # Extraer JSON de la respuesta
            import re as _re
            json_match = _re.search(r'\{[^{}]+\}', text, _re.DOTALL)
            if not json_match:
                logger.warning(f"[{gse_id}] LLM no retornó JSON válido")
                return {}

            mapping = _json.loads(json_match.group())

            # Validar que los valores sean case/control
            valid = {k: v for k, v in mapping.items()
                     if v in ("case", "control") and k in matrix_cols}

            logger.info(f"[{gse_id}] Capa 2 LLM mapeó {len(valid)}/{len(matrix_cols)} columnas")
            return valid

        except Exception as e:
            logger.warning(f"[{gse_id}] Capa 2 LLM falló: {e}")
            return {}


    def _sra_map_columns(self, matrix_cols: list, all_samples: list, gse_id: str) -> dict:
        """
        Capa 3: Busca tabla SRA RunInfo en NCBI para obtener mapeo
        Run/SampleName -> GSM -> label cuando Capa 2 LLM falla.

        Flujo:
            1. esearch SRA para obtener IDs de runs del GSE
            2. efetch RunInfo CSV con SampleName, BioSample, Run
            3. Cruzar SampleName/Run con columnas de la matriz
            4. Mapear GSM -> label desde all_samples

        Returns:
            dict {columna: 'case'|'control'} o {} si falla
        """
        try:
            import requests as _requests
            import csv as _csv
            import io as _io

            base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
            email = self.config.get("ncbi_email", "biosignal@research.com")

            # Paso 1: esearch SRA para el GSE
            logger.info(f"[{gse_id}] Capa 3: consultando SRA RunInfo...")
            # Paso 1a: obtener UID interno de GEO para el GSE
            geo_url = (
                f"{base}/esearch.fcgi?db=gds&term={gse_id}[Accession]"
                f"&retmax=5&retmode=json&email={email}"
            )
            r0 = _requests.get(geo_url, timeout=15)
            geo_ids = r0.json().get("esearchresult", {}).get("idlist", [])
            if not geo_ids:
                logger.warning(f"[{gse_id}] Capa 3: GSE no encontrado en GEO")
                return {}

            # Paso 1b: elink GEO -> SRA para obtener SRA IDs
            elink_url = (
                f"{base}/elink.fcgi?dbfrom=gds&db=sra&id={geo_ids[0]}"
                f"&retmode=json&email={email}"
            )
            r1 = _requests.get(elink_url, timeout=15)
            linksets = r1.json().get("linksets", [])
            ids = []
            for ls in linksets:
                for ld in ls.get("linksetdbs", []):
                    if ld.get("dbto") == "sra":
                        ids.extend(ld.get("links", []))
            if not ids:
                logger.warning(f"[{gse_id}] Capa 3: no se encontraron runs en SRA via elink")
                return {}
            search_url = f"{base}/esearch.fcgi?db=sra&term={gse_id}"  # unused now
            logger.info(f"[{gse_id}] Capa 3: {len(ids)} runs SRA encontrados via elink")

            # Paso 2: efetch RunInfo CSV
            fetch_url = (
                f"{base}/efetch.fcgi?db=sra&id={','.join(ids[:200])}"
                f"&rettype=runinfo&retmode=text&email={email}"
            )
            r2 = _requests.get(fetch_url, timeout=30)
            if r2.status_code != 200 or not r2.text.strip():
                logger.warning(f"[{gse_id}] Capa 3: efetch RunInfo vacío")
                return {}

            # Paso 3: parsear CSV y construir lookup
            reader = _csv.DictReader(_io.StringIO(r2.text))
            rows = list(reader)
            if not rows:
                return {}

            # Construir GSM -> label desde all_samples
            gsm_to_label = {
                s.get("gsm_id", ""): s.get("label", "unclassified")
                for s in all_samples
            }

            # Construir BioSample -> label via GSM
            biosample_to_label = {}
            for row in rows:
                gsm = row.get("SampleName", "")
                bs = row.get("BioSample", "")
                if gsm and bs and gsm in gsm_to_label:
                    biosample_to_label[bs] = gsm_to_label[gsm]

            # Paso 3b: consultar BioSample para obtener nombre de muestra del lab
            # BioSample title suele coincidir con columnas de la matriz
            biosample_to_colname = {}
            unique_biosamples = list({row.get("BioSample","") for row in rows if row.get("BioSample","")})
            for bs_id in unique_biosamples[:100]:
                try:
                    bs_resp = _requests.get(
                        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                        f"?db=biosample&id={bs_id}&retmode=text",
                        timeout=10
                    )
                    # Primera línea contiene el nombre: "1: Males-5xFAD_4"
                    first_line = bs_resp.text.strip().split("\n")[0]
                    if ": " in first_line:
                        col_name = first_line.split(": ", 1)[1].strip()
                        biosample_to_colname[bs_id] = col_name
                except Exception:
                    continue

            # Cruzar col_name -> label
            mapping = {}
            for bs_id, col_name in biosample_to_colname.items():
                label = biosample_to_label.get(bs_id, "")
                if label in ("case", "control"):
                    # Match exacto o parcial con matrix_cols
                    for col in matrix_cols:
                        if col == col_name or col in col_name or col_name in col:
                            mapping[col] = label
                            break

            # Fallback: cruzar campos candidatos directamente
            if not mapping:
                candidate_fields = ["SampleName", "Run", "Sample", "LibraryName", "Experiment"]
                for col in matrix_cols:
                    col_lower = col.lower()
                    for row in rows:
                        matched = False
                        for field in candidate_fields:
                            val = row.get(field, "").strip()
                            if val and (val == col or val.lower() == col_lower
                                        or col in val or val in col):
                                gsm = row.get("SampleName", "") or row.get("Sample", "")
                                label = gsm_to_label.get(gsm, "")
                                if label in ("case", "control"):
                                    mapping[col] = label
                                    matched = True
                                    break
                        if matched:
                            break

            groups = set(mapping.values())
            logger.info(
                f"[{gse_id}] Capa 3 SRA mapeó {len(mapping)}/{len(matrix_cols)} columnas "
                f"| grupos: {groups}"
            )
            return mapping

        except Exception as e:
            logger.warning(f"[{gse_id}] Capa 3 SRA falló: {e}")
            return {}

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
@click.option("--min-genes", default=5000, help="MÃ­nimo genes detectados para PASS")
@click.option("--min-samples", default=3, help="MÃ­nimo muestras por grupo")
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
