"""
BioSignal Discovery Engine
Agent 1: DatasetDiscoveryAgent
==============================
Responsabilidad: Dado un término de enfermedad, identificar los datasets de expresión
génica más relevantes disponibles públicamente en GEO/NCBI.

Fuente de datos: NCBI E-utilities (https://eutils.ncbi.nlm.nih.gov)
No requiere API key para <3 req/seg; para mayor throughput registrar API key gratuita.

Uso:
    python -m agents.discovery \\
        --disease 'pancreatic cancer' \\
        --max-datasets 20 \\
        --min-samples 10 \\
        --output data/discovery/
"""

import json
import requests
import time
import click
from pathlib import Path
from typing import Optional
from loguru import logger
from pydantic import ValidationError

from utils.geo_utils import NCBIClient
from schemas.dataset import DatasetMetadata, DiscoveryOutput
from utils.geo_utils import load_disease_aliases


class DatasetDiscoveryAgent:
    """
    Agente de descubrimiento de datasets en GEO/NCBI.

    Proceso:
        1. Consultar NCBI E-utilities API (esearch + esummary) con query:
           '{disease_name} expression profiling'
        2. Filtrar por tipo: expression profiling by high throughput sequencing / array
        3. Filtrar por número mínimo de muestras (sample_count >= min_samples)
        4. Rankear por: número de citas, número de muestras, año de publicación
        5. Retornar lista de GSE IDs con metadata

    Salida:
        data/discovery/datasets_{disease_name}.json
        Schema: [{gse_id, title, organism, sample_count, platform, year, pmid}]
    """

    def __init__(self, config: dict = None):
        self.config = config or {}
        self.ncbi = NCBIClient(
            email=self.config.get("ncbi_email", ""),
            api_key=self.config.get("ncbi_api_key", ""),
        )
        self.disease_aliases = load_disease_aliases()

    def resolve_disease_term(self, disease_name: str) -> str:
        """
        Resuelve el nombre de enfermedad al término canónico de GEO.

        Args:
            disease_name: Nombre coloquial de la enfermedad

        Returns:
            Término de búsqueda optimizado para GEO/NCBI
        """
        term = disease_name.lower().strip()
        # Buscar en aliases
        for category in self.disease_aliases.values():
            if isinstance(category, dict) and term in category:
                resolved = category[term]
                logger.info(f"Enfermedad '{disease_name}' → '{resolved}'")
                return resolved
        # Si no hay alias, construir query genérico
        return f"{disease_name} expression profiling"

    def search_datasets(
        self,
        disease_name: str,
        max_datasets: int = 20,
        min_samples: int = 10,
        data_types: list = None,
    ) -> list[DatasetMetadata]:
        """
        Busca y filtra datasets de GEO para una enfermedad dada.

        Args:
            disease_name: Nombre de la enfermedad a buscar
            max_datasets: Número máximo de datasets a recuperar
            min_samples: Mínimo de muestras por dataset
            data_types: Lista de tipos de datos ['RNA-seq', 'microarray']

        Returns:
            Lista de DatasetMetadata validados y rankeados
        """
        if data_types is None:
            data_types = ["RNA-seq", "microarray"]

        query_term = self.resolve_disease_term(disease_name)
        logger.info(f"Buscando datasets para: '{query_term}'")
        logger.info(f"Parámetros: max={max_datasets}, min_samples={min_samples}, tipos={data_types}")

        # Paso 1: esearch — obtener GSE IDs
        gse_ids = self.ncbi.esearch(
            db="gds",
            term=query_term,
            retmax=max_datasets * 3,  # Buscar 3x para filtrar después
        )
        logger.info(f"Encontrados {len(gse_ids)} IDs crudos en GEO")

        if not gse_ids:
            logger.warning(f"No se encontraron datasets para '{query_term}'")
            return []

        # Paso 2: esummary — obtener metadata de cada dataset
        raw_metadata = self.ncbi.esummary(db="gds", ids=gse_ids)

        # Paso 3: Filtrar y validar
        datasets = []
        for meta in raw_metadata:
            try:
                dataset = self._parse_and_filter(
                    meta,
                    min_samples=min_samples,
                    data_types=data_types,
                )
                if dataset:
                    datasets.append(dataset)
            except (ValidationError, KeyError) as e:
                logger.debug(f"Dataset descartado: {e}")
                continue

        logger.info(f"Datasets válidos después de filtros: {len(datasets)}")

        # Paso 3b: Capa 1 scoring de descargabilidad
        logger.info("Verificando descargabilidad de datasets...")
        for d in datasets:
            d.downloadability_score = self._check_downloadability(d.gse_id)
            logger.debug(f"  [{d.gse_id}] Downloadability score: {d.downloadability_score}")

        # Paso 4: Rankear y limitar
        ranked = self._rank_datasets(datasets)
        final = ranked[:max_datasets]

        logger.success(
            f"✓ Descubrimiento completado: {len(final)} datasets seleccionados "
            f"de {len(datasets)} candidatos"
        )
        return final

    def _parse_and_filter(
        self,
        raw_meta: dict,
        min_samples: int,
        data_types: list,
    ) -> Optional[DatasetMetadata]:
        """
        Parsea metadata cruda de NCBI y aplica filtros de calidad.

        Args:
            raw_meta: Metadata cruda del dataset desde esummary
            min_samples: Mínimo de muestras requeridas
            data_types: Tipos de datos permitidos

        Returns:
            DatasetMetadata validado o None si no pasa filtros
        """
        # Extraer tipo de dato
        gds_type = raw_meta.get("gdstype", "").lower()

        # Filtrar por tipo de datos
        is_rnaseq = "high throughput sequencing" in gds_type or "rna-seq" in gds_type
        is_array = "array" in gds_type or "microarray" in gds_type

        type_ok = (
            ("RNA-seq" in data_types and is_rnaseq) or
            ("microarray" in data_types and is_array)
        )
        if not type_ok:
            return None

        # Extraer número de muestras
        sample_count = int(raw_meta.get("n_samples", 0) or 0)
        if sample_count < min_samples:
            return None

        # Determinar tipo de datos
        detected_type = "RNA-seq" if is_rnaseq else "microarray"

        # Construir DatasetMetadata
        dataset = DatasetMetadata(
            gse_id=raw_meta.get("accession", ""),
            title=raw_meta.get("title", "")[:200],
            organism=raw_meta.get("taxon", "Homo sapiens"),
            sample_count=sample_count,
            platform=raw_meta.get("GPL", ""),
            year=int(str(raw_meta.get("pdat", "2000"))[:4]),
            pmid=raw_meta.get("pubmedids", [None])[0] if raw_meta.get("pubmedids") else None,
            data_type=detected_type,
            summary=raw_meta.get("summary", "")[:500],
        )
        return dataset


    def _check_downloadability(self, gse_id: str) -> int:
        """
        Capa 1: Scoring de descargabilidad consultando supplementary files de GEO.
        Scores: 3=counts consolidados, 2=RAW.tar, 1=solo series matrix, 0=sin archivos
        """
        try:
            gse_num = int(gse_id.replace("GSE", ""))
            prefix = str(gse_num)[:-3] if gse_num >= 1000 else "0"
            ftp_url = (
                f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{prefix}nnn"
                f"/{gse_id}/suppl/"
            )
            resp = requests.get(ftp_url, timeout=10)
            if resp.status_code != 200:
                return 1
            body = resp.text.lower()
            count_keywords = [
                "count", "counts", "rawcounts", "raw_count", "genecounts",
                "featurecounts", "starcounts", "count_matrix", "readcounts",
            ]
            has_counts = any(kw in body for kw in count_keywords)
            has_tabular = ".csv" in body or ".tsv" in body or ".txt" in body
            if has_counts and has_tabular:
                return 3
            if "_raw.tar" in body:
                return 2
            return 1
        except Exception as e:
            logger.debug(f"[{gse_id}] _check_downloadability error: {e}")
            return 1

    def _rank_datasets(self, datasets: list[DatasetMetadata]) -> list[DatasetMetadata]:
        """
        Rankea datasets por: citations (40%), sample_count (40%), year (20%).

        Args:
            datasets: Lista de DatasetMetadata a rankear

        Returns:
            Lista ordenada de mayor a menor score
        """
        if not datasets:
            return []

        # Normalizar valores para ranking
        max_samples = max(d.sample_count for d in datasets) or 1
        max_year = max(d.year for d in datasets) or 2000
        min_year = min(d.year for d in datasets) or 2000
        year_range = max(max_year - min_year, 1)

        for dataset in datasets:
            sample_score = dataset.sample_count / max_samples
            year_score = (dataset.year - min_year) / year_range
            # PMIDs como proxy de citaciones: datasets con PMID tienen +0.5
            citation_score = 1.0 if dataset.pmid else 0.5

            dl_score = getattr(dataset, 'downloadability_score', 1) / 3.0
            dataset.rank_score = (
                0.35 * citation_score +
                0.35 * sample_score +
                0.15 * year_score +
                0.15 * dl_score
            )

        return sorted(datasets, key=lambda d: d.rank_score, reverse=True)

    def save_results(
        self,
        datasets: list[DatasetMetadata],
        disease_name: str,
        output_dir: str,
    ) -> Path:
        """
        Guarda los resultados del descubrimiento en JSON.

        Args:
            datasets: Lista de datasets descubiertos
            disease_name: Nombre de la enfermedad
            output_dir: Directorio de salida

        Returns:
            Path al archivo JSON generado
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Sanitizar nombre para el archivo
        safe_name = disease_name.replace(" ", "_").replace("/", "_").lower()
        output_path = output_dir / f"datasets_{safe_name}.json"

        # Construir output completo
        output = DiscoveryOutput(
            disease_name=disease_name,
            total_found=len(datasets),
            datasets=datasets,
        )

        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(output.model_dump(), f, indent=2, default=str)

        logger.success(f"✓ Resultados guardados: {output_path}")
        logger.info(f"  → {len(datasets)} datasets | "
                    f"{sum(1 for d in datasets if d.data_type == 'RNA-seq')} RNA-seq | "
                    f"{sum(1 for d in datasets if d.data_type == 'microarray')} microarray")
        return output_path

    def run(
        self,
        disease_name: str,
        max_datasets: int = 20,
        min_samples: int = 10,
        data_types: list = None,
        output_dir: str = "data/discovery/",
    ) -> Path:
        """
        Ejecuta el pipeline completo del agente.

        Args:
            disease_name: Nombre de la enfermedad
            max_datasets: Máximo datasets a recuperar
            min_samples: Mínimo muestras por dataset
            data_types: Tipos de datos permitidos
            output_dir: Directorio de salida

        Returns:
            Path al archivo JSON de resultados
        """
        logger.info(f"=== Agent 1: DatasetDiscoveryAgent ===")
        logger.info(f"Enfermedad: '{disease_name}'")
        start_time = time.time()

        datasets = self.search_datasets(
            disease_name=disease_name,
            max_datasets=max_datasets,
            min_samples=min_samples,
            data_types=data_types,
        )

        output_path = self.save_results(
            datasets=datasets,
            disease_name=disease_name,
            output_dir=output_dir,
        )

        elapsed = time.time() - start_time
        logger.info(f"Tiempo total Agent 1: {elapsed:.1f}s")
        return output_path


# ============================================================
# CLI Interface
# ============================================================

@click.command()
@click.option("--disease", required=True, help="Nombre de la enfermedad a analizar")
@click.option("--max-datasets", default=20, show_default=True, help="Máximo datasets a recuperar")
@click.option("--min-samples", default=10, show_default=True, help="Mínimo muestras por dataset")
@click.option("--data-types", default="RNA-seq,microarray", show_default=True,
              help="Tipos de datos separados por coma")
@click.option("--output", default="data/discovery/", show_default=True,
              help="Directorio de salida")
def main(disease, max_datasets, min_samples, data_types, output):
    """
    Agent 1: Descubrimiento de datasets GEO para una enfermedad dada.

    Ejemplo:
        python -m agents.discovery --disease 'pancreatic cancer'
    """
    from utils.geo_utils import load_config
    config = load_config()

    agent = DatasetDiscoveryAgent(config=config.get("discovery", {}))
    agent.run(
        disease_name=disease,
        max_datasets=max_datasets,
        min_samples=min_samples,
        data_types=data_types.split(","),
        output_dir=output,
    )


if __name__ == "__main__":
    main()
