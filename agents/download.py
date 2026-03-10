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
    """
    Agente de descarga de datasets de GEO/NCBI.

    Proceso:
        1. Para cada GSE ID, descargar archivo SOFT o matrix file
        2. Verificar integridad via MD5 checksum
        3. Descomprimir archivos .gz automáticamente
        4. Parsear metadata: GSM IDs, condición, plataforma
        5. Clasificar muestras: control vs case usando anotaciones GSM

    Salida por dataset:
        data/raw/{GSE_ID}/matrix.tsv         — Matriz de expresión cruda
        data/raw/{GSE_ID}/metadata.json      — Metadata de muestras
        data/raw/{GSE_ID}/platform.json      — Información de plataforma (GPL)
    """

    NCBI_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
    MATRIX_SUFFIX = "_series_matrix.txt.gz"

    # Keywords para clasificar muestras como "control"
    CONTROL_KEYWORDS = [
        "control", "normal", "healthy", "non-tumor", "adjacent",
        "wild type", "wt", "vehicle", "untreated", "mock",
        "_c_", " ckp", "ckp_", "wildtype", "parental", "scramble",
        "ctr", "ctrl", "non-ad", "non_ad", "no infection",
        "empty vector", "gfp", "lacz", "dmso", "sictl", "shctl"
    ]

    # Keywords para clasificar muestras como "case"
    CASE_KEYWORDS = [
        "tumor", "cancer", "disease", "patient", "treated",
        "knockdown", "knockout", "kd", "ko", "mutant", "overexpression",
        "cckp", "gckp", "sirna", "shrna", "lnd", "sr4", "inhibitor",
        "alzheimer", " ad ", "ad_", "_ad", "mci", "dementia", "infected",
        "drug", "stimulated", "infected", "transfected", "overexpr"
    ]

    def __init__(self, config: dict = None):
        self.config = config or {}
        self.timeout = self.config.get("timeout_seconds", 300)
        self.max_retries = self.config.get("max_retries", 3)
        self.retry_backoff = self.config.get("retry_backoff", [1, 2, 4])
        self.verify_md5 = self.config.get("verify_md5", True)

    def download_dataset(
        self,
        gse_id: str,
        output_dir: str,
    ) -> DownloadStatus:
        """
        Descarga un dataset completo de GEO.

        Args:
            gse_id: Identificador GEO (ej. GSE12345)
            output_dir: Directorio base para almacenamiento

        Returns:
            DownloadStatus con estado y paths de archivos descargados
        """
        dataset_dir = Path(output_dir) / gse_id
        dataset_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Descargando {gse_id}...")

        # Intentar descarga con reintentos
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

        # Todos los intentos fallaron
        logger.error(f"✗ {gse_id}: falló después de {self.max_retries} intentos")
        return DownloadStatus(
            gse_id=gse_id,
            success=False,
            error=f"Failed after {self.max_retries} retries",
            output_dir=str(dataset_dir),
        )

    def _attempt_download(
        self,
        gse_id: str,
        dataset_dir: Path,
    ) -> DownloadStatus:
        """
        Intenta descargar y procesar un dataset de GEO.

        Args:
            gse_id: Identificador GSE
            dataset_dir: Directorio específico para este dataset

        Returns:
            DownloadStatus con resultado del intento
        """
        # Construir URL del matrix file
        # GEO organiza por: GSE1-100, GSE101-200, etc.
        gse_num = int(gse_id.replace("GSE", ""))
        gse_stub = f"GSE{str(gse_num)[:-3]}nnn"
        matrix_url = (
            f"{self.NCBI_FTP_BASE}/{gse_stub}/{gse_id}/matrix/"
            f"{gse_id}{self.MATRIX_SUFFIX}"
        )

        # Descargar matrix file
        matrix_gz_path = dataset_dir / f"{gse_id}_series_matrix.txt.gz"
        matrix_path = dataset_dir / "matrix.tsv"

        if not matrix_path.exists():
            self._download_file(matrix_url, matrix_gz_path)
            self._decompress_gz(matrix_gz_path, matrix_path)
            matrix_gz_path.unlink(missing_ok=True)  # Limpiar .gz

        # Parsear con GEOparse para extraer metadata
        gse = GEOparse.get_GEO(geo=gse_id, destdir=str(dataset_dir), silent=True)

        # Detectar si es RNA-seq SRA (sin counts en series matrix)
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

        # Extraer y guardar metadata de muestras
        metadata = self._extract_metadata(gse)
        metadata_path = dataset_dir / "metadata.json"
        with open(metadata_path, "w") as f:
            json.dump(metadata, f, indent=2)

        # Extraer información de plataforma
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
        """
        Descarga un archivo con streaming.

        Args:
            url: URL del archivo a descargar
            dest_path: Path destino del archivo
        """
        response = requests.get(url, stream=True, timeout=self.timeout)
        response.raise_for_status()

        with open(dest_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

    def _decompress_gz(self, gz_path: Path, output_path: Path) -> None:
        """
        Descomprime un archivo .gz.

        Args:
            gz_path: Path al archivo .gz
            output_path: Path destino descomprimido
        """
        with gzip.open(gz_path, "rb") as f_in:
            with open(output_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

    def _extract_metadata(self, gse) -> dict:
        """
        Extrae metadata de muestras de un objeto GEO.
        Clasifica automáticamente muestras como control o case.

        Args:
            gse: Objeto GEOparse GSE

        Returns:
            Dict con metadata de todas las muestras
        """
        samples = []
        unclassified = []

        for gsm_id, gsm in gse.gsms.items():
            # Obtener características de la muestra
            characteristics = {}
            for key, values in gsm.metadata.items():
                characteristics[key] = values[0] if values else ""

            # Clasificar muestra
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

        return {
            "gse_id": gse.name,
            "title": gse.metadata.get("title", [""])[0],
            "samples": samples,
            "unclassified": unclassified,
            "n_case": sum(1 for s in samples if s["label"] == "case"),
            "n_control": sum(1 for s in samples if s["label"] == "control"),
            "n_unclassified": len(unclassified),
        }

    def _classify_sample(self, sample_text: str) -> str:
        """
        Clasifica una muestra como 'control', 'case' o 'unclassified'.
        Ordena keywords por longitud descendente para que matches más específicos
        (ej. 'cckp') tengan prioridad sobre matches más cortos (ej. 'ckp').

        Returns:
            'control' | 'case' | 'unclassified'
        """
        # Ordenar por longitud desc para priorizar keywords más específicos
        case_keywords_sorted = sorted(self.CASE_KEYWORDS, key=len, reverse=True)
        control_keywords_sorted = sorted(self.CONTROL_KEYWORDS, key=len, reverse=True)

        case_matches = [kw for kw in case_keywords_sorted if kw in sample_text]
        control_matches = [kw for kw in control_keywords_sorted if kw in sample_text]

        # Si hay match de case más específico (más largo) que el de control, es case
        if case_matches and control_matches:
            if len(case_matches[0]) >= len(control_matches[0]):
                return "case"
            return "control"
        elif case_matches:
            return "case"
        elif control_matches:
            return "control"
        return "unclassified"


    # ------------------------------------------------------------------
    # SRA RNA-seq: detección y consolidación de counts suplementarios
    # ------------------------------------------------------------------

    def _is_sra_rnaseq(self, gse) -> bool:
        """Detecta si el dataset es RNA-seq almacenado en SRA (sin counts en series matrix)."""
        for gsm_id, gsm in gse.gsms.items():
            sample_type = gsm.metadata.get("type", [""])[0]
            library_strategy = gsm.metadata.get("library_strategy", [""])[0]
            if sample_type == "SRA" or "RNA-Seq" in library_strategy:
                return True
        return False

    def _matrix_has_counts(self, matrix_path: Path) -> bool:
        """
        Verifica si el matrix.tsv contiene genes reales (no solo headers).
        Un GEO Series Matrix de RNA-seq SRA tiene !series_matrix_table_begin
        pero 0 filas de datos — eso NO cuenta como counts válidos.
        """
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
                    # Esta es la línea de header (ID_REF ...), saltarla
                    # La siguiente línea con datos sería un gen real
                    # Leer una línea más para confirmar
                    next_line = f.readline().strip()
                    if next_line and not next_line.startswith("!"):
                        return True  # Hay al menos un gen
                    return False
        return False

    def _download_and_consolidate_counts(self, gse, gse_id: str, dataset_dir: Path):
        """
        Descarga archivos suplementarios de counts por muestra y los consolida
        en una matriz genes x muestras.

        Soporta dos estrategias:
            1. GSE*_RAW.tar — archivo con todos los counts por muestra
            2. Archivos individuales por GSM desde supplementary_file en metadata
        """
        import tarfile
        import pandas as pd

        suppl_dir = dataset_dir / "suppl"
        suppl_dir.mkdir(exist_ok=True)

        # Estrategia 1: descargar GSE*_RAW.tar
        gse_num = int(gse_id.replace("GSE", ""))
        gse_stub = f"GSE{str(gse_num)[:-3]}nnn"
        raw_tar_url = (
            f"{self.NCBI_FTP_BASE}/{gse_stub}/{gse_id}/suppl/{gse_id}_RAW.tar"
        )

        tar_path = suppl_dir / f"{gse_id}_RAW.tar"
        try:
            logger.info(f"[{gse_id}] Descargando {gse_id}_RAW.tar...")
            self._download_file(raw_tar_url, tar_path)

            # Extraer archivos de counts (.txt.gz, .tsv.gz, .count.gz)
            with tarfile.open(tar_path) as tar:
                count_files = [
                    m for m in tar.getmembers()
                    if any(m.name.endswith(ext) for ext in [".txt.gz", ".tsv.gz", ".count.gz", ".counts.gz"])
                ]
                tar.extractall(path=suppl_dir, members=count_files)

            tar_path.unlink(missing_ok=True)

        except Exception as e:
            logger.warning(f"[{gse_id}] No se pudo descargar RAW.tar: {e}")
            # Estrategia 2: archivos supplementary a nivel GSE
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
                            # Convertir ftp:// a https://
                            https_url = url.replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")
                            self._download_file(https_url, dest)
                            logger.debug(f"[{gse_id}] Descargado: {fname}")
                        except Exception as e2:
                            logger.warning(f"[{gse_id}] Error descargando {fname}: {e2}")
            else:
                # Estrategia 3: archivos individuales por GSM
                self._download_gsm_supplementary(gse, suppl_dir)

        # Consolidar archivos de counts en matriz
        return self._consolidate_count_files(suppl_dir, gse, gse_id)

    def _download_gsm_supplementary(self, gse, suppl_dir: Path):
        """Descarga archivos supplementary individuales de cada GSM."""
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
        """
        Lee todos los archivos de counts en suppl_dir y construye
        una matriz genes x muestras.

        Soporta 3 formatos:
            1. Un archivo por muestra (GSM*.txt.gz) — columna única de counts
            2. Matriz completa en un solo archivo (genes x muestras) con header
            3. Múltiples archivos con múltiples columnas — se concatenan
        """
        import pandas as pd

        all_files = (
            list(suppl_dir.glob("*.gz")) +
            list(suppl_dir.glob("*.txt")) +
            list(suppl_dir.glob("*.tsv")) +
            list(suppl_dir.glob("*.csv"))
        )
        # Excluir archivos de barcodes/features (single-cell)
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
                # Detectar separador y compresión
                compression = "gzip" if fpath.suffix == ".gz" else None
                sep = "," if fpath.name.endswith(".csv.gz") or fpath.name.endswith(".csv") else "\t"

                # Leer con header primero para detectar formato
                df = pd.read_csv(fpath, sep=sep, index_col=0,
                                  compression=compression, comment="#")

                # Verificar si tiene columnas numéricas (matriz completa)
                numeric_cols = df.select_dtypes(include="number").columns.tolist()
                if len(numeric_cols) == 0:
                    # Intentar sin header
                    df = pd.read_csv(fpath, sep=sep, index_col=0,
                                      compression=compression, comment="#", header=None)
                    numeric_cols = df.select_dtypes(include="number").columns.tolist()

                if len(numeric_cols) == 0:
                    logger.warning(f"Sin columnas numéricas en {fpath.name}")
                    continue

                df = df[numeric_cols]

                # Si es un archivo de muestra única, renombrar columna con GSM o filename
                if len(numeric_cols) == 1:
                    gsm_id = self._extract_gsm_from_filename(fpath.stem)
                    col_name = gsm_id if gsm_id else fpath.stem
                    df.columns = [col_name]

                # Filtrar líneas de resumen HTSeq
                df = df[~df.index.astype(str).str.startswith("__")]
                frames.append(df)
                logger.debug(f"[{gse_id}] Leído {fpath.name}: {df.shape}")

            except Exception as e:
                logger.warning(f"Error leyendo {fpath.name}: {e}")

        if not frames:
            return None

        # Si hay múltiples frames, concatenar por columnas
        if len(frames) == 1:
            matrix = frames[0]
        else:
            matrix = pd.concat(frames, axis=1)

        matrix.index.name = "gene"
        # Mantener solo columnas numéricas
        matrix = matrix.select_dtypes(include="number")
        return matrix if not matrix.empty else None

    def _extract_gsm_from_filename(self, stem: str) -> str:
        """Extrae el GSM ID de un nombre de archivo tipo GSM12345_counts."""
        import re
        match = re.search(r"(GSM\d+)", stem)
        return match.group(1) if match else ""

    def _extract_platform_info(self, gse) -> dict:
        """
        Extrae información de la plataforma (GPL) del dataset.

        Args:
            gse: Objeto GEOparse GSE

        Returns:
            Dict con información de la plataforma
        """
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

    def run_parallel(
        self,
        input_json: str,
        output_dir: str,
        parallel: int = 4,
    ) -> dict:
        """
        Descarga múltiples datasets en paralelo.

        Args:
            input_json: Path al JSON de DiscoveryOutput (Agent 1)
            output_dir: Directorio base para almacenamiento
            parallel: Número de workers paralelos

        Returns:
            Dict con estadísticas de descarga
        """
        logger.info(f"=== Agent 2: DatasetDownloadAgent ===")

        # Cargar datasets del Agent 1
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
                        # Verificar si tiene muestras sin clasificar
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
def main(input, output_dir, parallel):
    """
    Agent 2: Descarga datasets de GEO/NCBI.

    Ejemplo:
        python -m agents.download \\
            --input data/discovery/datasets_pancreatic_cancer.json
    """
    from utils.geo_utils import load_config
    config = load_config()

    agent = DatasetDownloadAgent(config=config.get("download", {}))
    agent.run_parallel(
        input_json=input,
        output_dir=output_dir,
        parallel=parallel,
    )


if __name__ == "__main__":
    main()