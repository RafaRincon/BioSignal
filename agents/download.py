"""
BioSignal Discovery Engine
Agent 2: DatasetDownloadAgent
==============================
Responsabilidad: Descargar los archivos de datos crudos para cada GSE ID
identificado por el DatasetDiscoveryAgent.

Uso:
    python -m agents.download \\
        --input data/discovery/datasets_pancreatic_cancer.json \\
        --output-dir data/raw/ \\
        --parallel 4
"""

import gzip
import hashlib
import json
import shutil
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

import GEOparse
import click
import requests
from loguru import logger

from schemas.dataset import DatasetMetadata, DownloadStatus


class DatasetDownloadAgent:
    NCBI_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
    MATRIX_SUFFIX = "_series_matrix.txt.gz"

    CONTROL_KEYWORDS = [
        "control", "normal", "healthy", "non-tumor", "adjacent",
        "wild type", "wt", "vehicle", "untreated", "mock",
        "_c_", " ckp", "ckp_", "wildtype", "parental", "scramble",
        "ctr", "ctrl", "non-ad", "non_ad", "no infection",
        "empty vector", "gfp", "lacz", "dmso", "sictl", "shctl"
    ]

    CASE_KEYWORDS = [
        "tumor", "cancer", "disease", "patient", "treated",
        "knockdown", "knockout", "kd", "ko", "mutant", "overexpression",
        "cckp", "gckp", "sirna", "shrna", "lnd", "sr4", "inhibitor",
        "alzheimer", " ad ", "ad_", "_ad", "mci", "dementia", "infected",
        "drug", "stimulated", "infected", "transfected", "overexpr"
    ]

    UNCLASSIFIED_THRESHOLD = 0.30

    def __init__(self, config: dict = None, llm_config: dict = None, topic: str = None):
        self.config = config or {}
        self.timeout = self.config.get("timeout_seconds", 300)
        self.max_retries = self.config.get("max_retries", 3)
        self.retry_backoff = self.config.get("retry_backoff", [1, 2, 4])
        self.verify_md5 = self.config.get("verify_md5", True)

        self.topic = topic or ""
        self._classifier = None
        if llm_config:
            try:
                from utils.experiment_classifier import ExperimentDesignClassifier
                self._classifier = ExperimentDesignClassifier(llm_config)
                logger.info(f"[DownloadAgent] ExperimentDesignClassifier activo (backend={llm_config.get('backend')})")
            except Exception as e:
                logger.warning(f"[DownloadAgent] ExperimentDesignClassifier no disponible: {e}")

    def download_dataset(self, gse_id: str, output_dir: str) -> DownloadStatus:
        dataset_dir = Path(output_dir) / gse_id
        dataset_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Descargando {gse_id}...")

        for attempt in range(self.max_retries):
            try:
                result = self._attempt_download(gse_id, dataset_dir)
                if result.success:
                    logger.success(f"✓ {gse_id}: descargado y procesado")
                    return result
            except Exception as e:
                wait_time = self.retry_backoff[min(attempt, len(self.retry_backoff) - 1)]
                logger.warning(
                    f"  Intento {attempt + 1}/{self.max_retries} fallido para {gse_id}: {e}. "
                    f"Esperando {wait_time}s..."
                )
                time.sleep(wait_time)

        logger.error(f"✗ {gse_id}: falló después de {self.max_retries} intentos")
        return DownloadStatus(
            gse_id=gse_id,
            success=False,
            error=f"Failed after {self.max_retries} retries",
            output_dir=str(dataset_dir),
        )

    def _attempt_download(self, gse_id: str, dataset_dir: Path) -> DownloadStatus:
        gse_num = int(gse_id.replace("GSE", ""))
        gse_stub = f"GSE{str(gse_num)[:-3]}nnn"
        matrix_url = (
            f"{self.NCBI_FTP_BASE}/{gse_stub}/{gse_id}/matrix/"
            f"{gse_id}{self.MATRIX_SUFFIX}"
        )

        matrix_gz_path = dataset_dir / f"{gse_id}_series_matrix.txt.gz"
        matrix_path = dataset_dir / "matrix.tsv"

        if not matrix_path.exists():
            self._download_file(matrix_url, matrix_gz_path)
            self._decompress_gz(matrix_gz_path, matrix_path)
            matrix_gz_path.unlink(missing_ok=True)

        gse = GEOparse.get_GEO(geo=gse_id, destdir=str(dataset_dir), silent=True)

        if self._is_sra_rnaseq(gse) and not self._matrix_has_counts(matrix_path):
            logger.info(f"[{gse_id}] RNA-seq SRA detectado — descargando supplementary counts")
            consolidated = self._download_and_consolidate_counts(gse, gse_id, dataset_dir)
            if consolidated is not None and not consolidated.empty:
                logger.success(
                    f"[{gse_id}] Counts consolidados: "
                    f"{consolidated.shape[0]} genes x {consolidated.shape[1]} muestras"
                )
                consolidated.to_csv(matrix_path, sep="\t")
            else:
                logger.warning(f"[{gse_id}] No se pudieron consolidar counts supplementary")

        metadata = self._extract_metadata(gse)
        metadata_path = dataset_dir / "metadata.json"
        with open(metadata_path, "w") as f:
            json.dump(metadata, f, indent=2)

        platform_info = self._extract_platform_info(gse)
        platform_path = dataset_dir / "platform.json"
        with open(platform_path, "w") as f:
            json.dump(platform_info, f, indent=2)

        return DownloadStatus(
            gse_id=gse_id,
            success=True,
            output_dir=str(dataset_dir),
            matrix_path=str(matrix_path),
            metadata_path=str(metadata_path),
            platform_path=str(platform_path),
            n_samples=len(metadata.get("samples", [])),
            platform=platform_info.get("gpl_id", ""),
        )

    def _download_file(self, url: str, dest_path: Path) -> None:
        response = requests.get(url, stream=True, timeout=self.timeout)
        response.raise_for_status()
        with open(dest_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

    def _decompress_gz(self, gz_path: Path, output_path: Path) -> None:
        with gzip.open(gz_path, "rb") as f_in:
            with open(output_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

    def _extract_metadata(self, gse) -> dict:
        samples = []
        unclassified = []

        for gsm_id, gsm in gse.gsms.items():
            characteristics = {}
            for key, values in gsm.metadata.items():
                characteristics[key] = values[0] if values else ""

            sample_text = " ".join([
                characteristics.get("title", ""),
                characteristics.get("source_name_ch1", ""),
                characteristics.get("characteristics_ch1", ""),
                characteristics.get("description", ""),
            ]).lower()

            label = self._classify_sample(sample_text)

            sample_info = {
                "gsm_id": gsm_id,
                "title": characteristics.get("title", ""),
                "label": label,
                "characteristics": characteristics,
            }

            if label == "unclassified":
                unclassified.append(sample_info)
            else:
                samples.append(sample_info)

        gse_title = gse.metadata.get("title", [""])[0]
        n_case = sum(1 for s in samples if s["label"] == "case")
        n_total = len(samples) + len(unclassified)
        unclassified_fraction = len(unclassified) / max(n_total, 1)

        classifier_result = None
        if self._classifier and (n_case == 0 or unclassified_fraction > self.UNCLASSIFIED_THRESHOLD):
            logger.info(
                f"[{gse.name}] Activando ExperimentDesignClassifier "
                f"(n_case={n_case}, unclassified={len(unclassified)}/{n_total})"
            )
            all_samples_for_classifier = [
                {"gsm_id": s["gsm_id"], "title": s["title"]}
                for s in (samples + unclassified)
            ]
            classifier_result = self._classifier.classify(
                gse_id=gse.name,
                gse_title=gse_title,
                topic=self.topic,
                samples=all_samples_for_classifier,
            )

            if classifier_result and classifier_result["valid"]:
                classification_map = classifier_result["classification"]
                reclassified = []
                still_unclassified = []

                for s in (samples + unclassified):
                    new_label = classification_map.get(s["gsm_id"], "unclassified")
                    if new_label == "exclude":
                        new_label = "unclassified"
                    s_updated = {**s, "label": new_label}
                    if new_label in ("case", "control"):
                        reclassified.append(s_updated)
                    else:
                        still_unclassified.append(s_updated)

                samples = reclassified
                unclassified = still_unclassified
                logger.success(
                    f"[{gse.name}] Reclasificado: "
                    f"case={sum(1 for s in samples if s['label']=='case')} "
                    f"control={sum(1 for s in samples if s['label']=='control')} "
                    f"unclassified={len(unclassified)} | "
                    f"axis='{classifier_result.get('relevant_axis')}'"
                )
            else:
                reason = classifier_result.get("validation_message", "desconocido") if classifier_result else "sin respuesta"
                logger.warning(f"[{gse.name}] ExperimentDesignClassifier inválido: {reason}")

        return {
            "gse_id": gse.name,
            "title": gse_title,
            "samples": samples,
            "unclassified": unclassified,
            "n_case": sum(1 for s in samples if s["label"] == "case"),
            "n_control": sum(1 for s in samples if s["label"] == "control"),
            "n_unclassified": len(unclassified),
            "classifier_used": classifier_result is not None,
            "classifier_result": {
                "axes": classifier_result.get("axes", []),
                "relevant_axis": classifier_result.get("relevant_axis", ""),
                "reasoning": classifier_result.get("reasoning", ""),
                "valid": classifier_result.get("valid", False),
            } if classifier_result else None,
        }

    def _classify_sample(self, sample_text: str) -> str:
        case_keywords_sorted = sorted(self.CASE_KEYWORDS, key=len, reverse=True)
        control_keywords_sorted = sorted(self.CONTROL_KEYWORDS, key=len, reverse=True)

        case_matches = [kw for kw in case_keywords_sorted if kw in sample_text]
        control_matches = [kw for kw in control_keywords_sorted if kw in sample_text]

        if case_matches and control_matches:
            if len(case_matches[0]) >= len(control_matches[0]):
                return "case"
            return "control"
        elif case_matches:
            return "case"
        elif control_matches:
            return "control"
        return "unclassified"

    def _is_sra_rnaseq(self, gse) -> bool:
        for gsm_id, gsm in gse.gsms.items():
            sample_type = gsm.metadata.get("type", [""])[0]
            library_strategy = gsm.metadata.get("library_strategy", [""])[0]
            if sample_type == "SRA" or "RNA-Seq" in library_strategy:
                return True
        return False

    def _matrix_has_counts(self, matrix_path: Path) -> bool:
        if not matrix_path.exists():
            return False
        in_table = False
        with open(matrix_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                stripped = line.strip().lower()
                if stripped == "!series_matrix_table_begin":
                    in_table = True
                    continue
                if stripped == "!series_matrix_table_end":
                    break
                if in_table and stripped and not stripped.startswith("!"):
                    next_line = f.readline().strip()
                    if next_line and not next_line.startswith("!"):
                        return True
                    return False
        return False

    def _download_and_consolidate_counts(self, gse, gse_id: str, dataset_dir: Path):
        import tarfile
        import pandas as pd

        suppl_dir = dataset_dir / "suppl"
        suppl_dir.mkdir(exist_ok=True)

        gse_num = int(gse_id.replace("GSE", ""))
        gse_stub = f"GSE{str(gse_num)[:-3]}nnn"
        raw_tar_url = (
            f"{self.NCBI_FTP_BASE}/{gse_stub}/{gse_id}/suppl/{gse_id}_RAW.tar"
        )

        tar_path = suppl_dir / f"{gse_id}_RAW.tar"
        try:
            logger.info(f"[{gse_id}] Descargando {gse_id}_RAW.tar...")
            self._download_file(raw_tar_url, tar_path)

            with tarfile.open(tar_path) as tar:
                count_files = [
                    m for m in tar.getmembers()
                    if any(m.name.endswith(ext) for ext in [".txt.gz", ".tsv.gz", ".count.gz", ".counts.gz"])
                ]
                tar.extractall(path=suppl_dir, members=count_files)

            tar_path.unlink(missing_ok=True)

        except Exception as e:
            logger.warning(f"[{gse_id}] No se pudo descargar RAW.tar: {e}")
            gse_suppl = gse.metadata.get("supplementary_file", [])
            if gse_suppl:
                logger.info(f"[{gse_id}] Intentando supplementary a nivel GSE: {len(gse_suppl)} archivos")
                for url in gse_suppl:
                    if not url or url == "NONE":
                        continue
                    fname = url.split("/")[-1]
                    dest = suppl_dir / fname
                    if not dest.exists():
                        try:
                            https_url = url.replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")
                            self._download_file(https_url, dest)
                            logger.debug(f"[{gse_id}] Descargado: {fname}")
                        except Exception as e2:
                            logger.warning(f"[{gse_id}] Error descargando {fname}: {e2}")
            else:
                self._download_gsm_supplementary(gse, suppl_dir)

        return self._consolidate_count_files(suppl_dir, gse, gse_id)

    def _download_gsm_supplementary(self, gse, suppl_dir: Path):
        for gsm_id, gsm in gse.gsms.items():
            suppl_files = gsm.metadata.get("supplementary_file_1", [])
            for url in suppl_files:
                if not url or url == "NONE":
                    continue
                fname = url.split("/")[-1]
                dest = suppl_dir / fname
                if dest.exists():
                    continue
                try:
                    self._download_file(url, dest)
                    logger.debug(f"Descargado: {fname}")
                except Exception as e:
                    logger.warning(f"Error descargando {fname}: {e}")

    def _consolidate_count_files(self, suppl_dir: Path, gse, gse_id: str):
        import pandas as pd

        all_files = (
            list(suppl_dir.glob("*.gz")) +
            list(suppl_dir.glob("*.txt")) +
            list(suppl_dir.glob("*.tsv")) +
            list(suppl_dir.glob("*.csv"))
        )
        all_files = [
            f for f in all_files
            if not any(x in f.name for x in ["barcodes", "features", "matrix.mtx"])
        ]

        if not all_files:
            logger.warning(f"[{gse_id}] Sin archivos de counts en {suppl_dir}")
            return None

        frames = []
        for fpath in sorted(all_files):
            try:
                compression = "gzip" if fpath.suffix == ".gz" else None
                sep = "," if fpath.name.endswith(".csv.gz") or fpath.name.endswith(".csv") else "\t"

                df = pd.read_csv(fpath, sep=sep, index_col=0,
                                  compression=compression, comment="#")

                numeric_cols = df.select_dtypes(include="number").columns.tolist()
                if len(numeric_cols) == 0:
                    df = pd.read_csv(fpath, sep=sep, index_col=0,
                                      compression=compression, comment="#", header=None)
                    numeric_cols = df.select_dtypes(include="number").columns.tolist()

                if len(numeric_cols) == 0:
                    logger.warning(f"Sin columnas numéricas en {fpath.name}")
                    continue

                df = df[numeric_cols]

                if len(numeric_cols) == 1:
                    gsm_id = self._extract_gsm_from_filename(fpath.stem)
                    col_name = gsm_id if gsm_id else fpath.stem
                    df.columns = [col_name]

                df = df[~df.index.astype(str).str.startswith("__")]
                frames.append(df)
                logger.debug(f"[{gse_id}] Leído {fpath.name}: {df.shape}")

            except Exception as e:
                logger.warning(f"Error leyendo {fpath.name}: {e}")

        if not frames:
            return None

        # FIX: filtrar frames incompatibles antes de concatenar.
        # Algunos datasets mezclan archivos de counts con archivos de anotación/metadata
        # (ej. GEOgeneinfo, GEOsamples) con n_rows muy distintos — causan
        # "Reindexing only valid with uniquely valued Index objects" en pd.concat.
        if len(frames) == 1:
            matrix = frames[0]
        else:
            frames_sorted = sorted(frames, key=lambda f: len(f), reverse=True)
            ref_frame = frames_sorted[0]
            ref_nrows = len(ref_frame)

            compatible = [ref_frame]
            for f in frames_sorted[1:]:
                if len(f) < ref_nrows * 0.5:
                    logger.debug(
                        f"[{gse_id}] Descartando frame ({len(f)} filas) "
                        f"incompatible con referencia ({ref_nrows} filas)"
                    )
                    continue
                if not f.index.is_unique:
                    logger.warning(
                        f"[{gse_id}] Frame con índice duplicado descartado "
                        f"({len(f)} filas, {len(f.columns)} cols)"
                    )
                    continue
                compatible.append(f)

            if len(compatible) == 1:
                matrix = compatible[0]
            else:
                if not compatible[0].index.is_unique:
                    compatible[0] = compatible[0][~compatible[0].index.duplicated(keep="first")]
                    logger.warning(f"[{gse_id}] Duplicados en índice de referencia — conservando primera ocurrencia")
                matrix = pd.concat(compatible, axis=1)

        matrix.index.name = "gene"
        matrix = matrix.select_dtypes(include="number")
        return matrix if not matrix.empty else None

    def _extract_gsm_from_filename(self, stem: str) -> str:
        import re
        match = re.search(r"(GSM\d+)", stem)
        return match.group(1) if match else ""

    def _extract_platform_info(self, gse) -> dict:
        platform_info = {
            "gpl_id": "",
            "title": "",
            "organism": "",
            "technology": "",
        }

        if gse.gpls:
            gpl_id = list(gse.gpls.keys())[0]
            gpl = gse.gpls[gpl_id]
            platform_info.update({
                "gpl_id": gpl_id,
                "title": gpl.metadata.get("title", [""])[0],
                "organism": gpl.metadata.get("organism", [""])[0],
                "technology": gpl.metadata.get("technology", [""])[0],
            })

        return platform_info

    def run_parallel(self, input_json: str, output_dir: str, parallel: int = 4) -> dict:
        logger.info(f"=== Agent 2: DatasetDownloadAgent ===")

        with open(input_json) as f:
            discovery_output = json.load(f)

        datasets = discovery_output.get("datasets", [])
        gse_ids = [d["gse_id"] for d in datasets if d.get("gse_id")]

        logger.info(f"Descargando {len(gse_ids)} datasets con {parallel} workers...")
        start_time = time.time()

        results = {"success": [], "failed": [], "review": []}

        with ThreadPoolExecutor(max_workers=parallel) as executor:
            future_to_gse = {
                executor.submit(self.download_dataset, gse_id, output_dir): gse_id
                for gse_id in gse_ids
            }

            for future in as_completed(future_to_gse):
                gse_id = future_to_gse[future]
                try:
                    status = future.result()
                    if status.success:
                        if status.n_unclassified and status.n_unclassified > 0:
                            results["review"].append(gse_id)
                            logger.warning(f"  → {gse_id}: {status.n_unclassified} muestras sin clasificar → data/review/")
                        else:
                            results["success"].append(gse_id)
                    else:
                        results["failed"].append(gse_id)
                except Exception as e:
                    logger.error(f"Error procesando {gse_id}: {e}")
                    results["failed"].append(gse_id)

        elapsed = time.time() - start_time
        logger.success(
            f"\n✓ Descarga completada en {elapsed:.0f}s:\n"
            f"  Exitosos:      {len(results['success'])}\n"
            f"  Para revisión: {len(results['review'])}\n"
            f"  Fallidos:      {len(results['failed'])}"
        )

        if results["failed"]:
            logger.warning(f"Datasets fallidos: {', '.join(results['failed'])}")

        return results


# ============================================================
# CLI Interface
# ============================================================

@click.command()
@click.option("--input", required=True, help="Path al JSON de descubrimiento (Agent 1)")
@click.option("--output-dir", default="data/raw/", show_default=True,
              help="Directorio de salida")
@click.option("--parallel", default=4, show_default=True,
              help="Número de descargas paralelas")
@click.option("--topic", default="", help="Tema de búsqueda para clasificación semántica (ej. 'alzheimer')")
@click.option("--config", "config_path", default="config/settings.yaml",
              help="Path a settings.yaml")
def main(input, output_dir, parallel, topic, config_path):
    """
    Agent 2: Descarga datasets de GEO/NCBI.

    Ejemplo:
        python -m agents.download \\
            --input data/discovery/datasets_alzheimer.json \\
            --topic "alzheimer"
    """
    from utils.geo_utils import load_config
    config = load_config(config_path)

    agent = DatasetDownloadAgent(
        config=config.get("download", {}),
        llm_config=config.get("llm", {}),
        topic=topic,
    )
    agent.run_parallel(
        input_json=input,
        output_dir=output_dir,
        parallel=parallel,
    )


if __name__ == "__main__":
    main()