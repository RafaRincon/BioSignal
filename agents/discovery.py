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

from dotenv import load_dotenv
from pathlib import Path as _Path
load_dotenv(_Path(__file__).parent.parent / ".env")
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
import os

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
        # LLM config para clasificación semántica de archivos
        llm_cfg = self.config.get("llm", {})
        self.ollama_url = llm_cfg.get("ollama_url", "http://localhost:11434")
        self.llm_model = llm_cfg.get("model", "gemma3:4b")
        self.llm_backend = llm_cfg.get("backend", "ollama")
        self.llm_client = None
        if self.llm_backend == "openai":
            import openai as _openai
            api_key = os.getenv("OPENAI_API_KEY", llm_cfg.get("api_key", ""))
            base_url = llm_cfg.get("openai_base_url", "")
            api_version = llm_cfg.get("openai_api_version", "")
            if api_version and base_url:
                self.llm_client = _openai.AzureOpenAI(api_key=api_key, azure_endpoint=base_url, api_version=api_version)
            elif base_url:
                self.llm_client = _openai.OpenAI(api_key=api_key, base_url=base_url)
            else:
                self.llm_client = _openai.OpenAI(api_key=api_key)

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
        return f'"{disease_name}"[Title]'

    def search_datasets(
        self,
        disease_name: str,
        max_datasets: int = 20,
        min_samples: int = 10,
        max_samples: int = 1000,
        data_types: list = None,
        organism: str = None,
        year_start: int = None,
        year_end: int = None,
    ) -> list[DatasetMetadata]:

        if data_types is None:
            data_types = ["RNA-seq"]

        DATA_TYPE_FILTERS = {
            "RNA-seq": '"Expression profiling by high throughput sequencing"[DataSet Type]',
            "microarray": '"Expression profiling by array"[DataSet Type]',
        }

        query_term = self.resolve_disease_term(disease_name)

        # Filtro nativo de tipo de dato
        if data_types and len(data_types) == 1 and data_types[0] in DATA_TYPE_FILTERS:
            query_term = f'{query_term} AND {DATA_TYPE_FILTERS[data_types[0]]}'
        elif data_types and len(data_types) > 1:
            type_filters = [DATA_TYPE_FILTERS[dt] for dt in data_types if dt in DATA_TYPE_FILTERS]
            if type_filters:
                query_term = f'{query_term} AND ({" OR ".join(type_filters)})'

        # Filtro de organismo
        if organism:
            query_term = f'{query_term} AND "{organism}"[Organism]'

        # Filtro de rango de muestras
        query_term = f'{query_term} AND {min_samples}:{max_samples}[Number of Samples]'

        # Filtro de años
        if year_start or year_end:
            y_start = year_start or 2000
            y_end = year_end or 2026
            query_term = f'{query_term} AND {y_start}:{y_end}[PDAT]'

        # Solo Series GSE
        query_term = f'{query_term} AND gse[Entry Type]'

        logger.info(f"Query NCBI: {query_term}")
        logger.info(f"Parámetros: max={max_datasets}, min_samples={min_samples}, tipos={data_types}")

        # B?squeda iterativa hasta completar max_datasets procesables (score >= 3)
        processable = []
        seen_ids = set()
        max_attempts = 3
        attempt = 0

        while len(processable) < max_datasets and attempt < max_attempts:
            attempt += 1
            retmax = max_datasets * (attempt + 1)

            gse_ids = self.ncbi.esearch(
                db="gds",
                term=query_term,
                retmax=retmax,
            )
            logger.info(f"Intento {attempt}: {len(gse_ids)} IDs crudos en GEO")

            if not gse_ids:
                logger.warning(f"No se encontraron datasets para '{query_term}'")
                break

            # Solo procesar IDs nuevos
            new_ids = [i for i in gse_ids if i not in seen_ids]
            seen_ids.update(new_ids)

            if not new_ids:
                logger.info("No hay IDs nuevos disponibles ? agotados los resultados NCBI")
                break

            raw_metadata = self.ncbi.esummary(db="gds", ids=new_ids)

            # Evaluar cada ID nuevo individualmente
            for meta in raw_metadata:
                if len(processable) >= max_datasets:
                    break
                try:
                    dataset = self._parse_and_filter(
                        meta,
                        min_samples=min_samples,
                        data_types=data_types,
                    )
                    if not dataset:
                        continue
                    dataset.downloadability_score, dataset.supplementary_files = self._check_downloadability(dataset.gse_id, summary=getattr(dataset, "summary", ""))
                    dataset._dl_checked = True
                    logger.debug(f"  [{dataset.gse_id}] Score: {dataset.downloadability_score} | Archivos: {dataset.supplementary_files}")
                    if dataset.downloadability_score >= 3:
                        processable.append(dataset)
                except (ValidationError, KeyError) as e:
                    logger.debug(f"Dataset descartado: {e}")
                    continue

            if len(processable) >= max_datasets:
                break

            remaining = max_datasets - len(processable)
            logger.info(f"Intento {attempt}: {len(processable)} procesables, faltan {remaining} ? ampliando b?squeda")

        if len(processable) < max_datasets:
            logger.warning(
                f"Solo se encontraron {len(processable)} datasets procesables "
                f"de {max_datasets} solicitados ? RAW.tar excluidos autom?ticamente"
            )

        # Rankear y limitar
        ranked = self._rank_datasets(processable)
        final = ranked[:max_datasets]

        logger.success(
            f"✓ Descubrimiento completado: {len(final)} datasets seleccionados"
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
    
    def _classify_files_with_llm(self, files: list[str]) -> tuple[int, str]:
        """
        Usa LLM local para clasificar semánticamente archivos suplementarios GEO.
        El LLM clasifica en categorías de texto; el código asigna el score.
        Retorna (score, reason). Fallback a (-1, "") si falla.
        """
        CATEGORY_TO_SCORE = {
            "RAW_COUNTS": 4,
            "NORMALIZED": 3,
            "OTHER": 2,
            "RAW_ARCHIVE": 1,
        }
        try:
            prompt = (
                "You are analyzing GEO supplementary filenames for a bioinformatics pipeline.\n\n"
                f"Files: {files}\n\n"
                "Classify the BEST available file into exactly one category:\n\n"
                "RAW_COUNTS: unnormalized integer count matrix ready for DESeq2. "
                "Examples: GSE254877_raw_counts.txt.gz, GSE312093_gene_count.txt.gz, GSE309948_COUNTS_093022.txt.gz. "
                "Key signals: words like 'count', 'counts', 'raw_count' in the filename WITHOUT 'normalized'.\n\n"
                "NORMALIZED: expression matrix that has been normalized. "
                "Examples: GSE254877_normalized_counts.csv.gz, GSE293164_normalized_gene_counts.csv.gz, GSE254877_fpkm.txt.gz. "
                "Key signals: words like 'normalized', 'normalised', 'fpkm', 'tpm', 'cpm', 'rpkm' anywhere in the filename.\n\n"
                "RAW_ARCHIVE: a .tar archive of raw sequencing or microarray files (FASTQ, BAM, CEL). "
                "Examples: GSE255433_RAW.tar, GSE307642_RAW.tar. "
                "Key signals: filename ends in _RAW.tar or RAW.tar. "
                "WARNING: do NOT confuse with RAW_COUNTS — a .tar file is never a count matrix.\n\n"
                "OTHER: any other file that does not fit the above categories.\n\n"
                "RULES:\n"
                "- If multiple files are present, return the category of the BEST one (RAW_COUNTS > NORMALIZED > OTHER > RAW_ARCHIVE).\n"
                "- If any filename contains 'normalized', 'fpkm', 'tpm', 'cpm', or 'rpkm', category MUST be NORMALIZED unless a RAW_COUNTS file also exists.\n"
                "- If any filename ends in '.gct.gz' or '.gct', category MUST be OTHER ? GCT format is not supported by the pipeline.\n"
                "- If the only expression file ends in _RAW.tar, category MUST be RAW_ARCHIVE.\n\n"
                "Respond ONLY with valid JSON: "
                "{\"category\": \"<RAW_COUNTS|NORMALIZED|RAW_ARCHIVE|OTHER>\", \"reason\": \"<one sentence>\"}"
            )
            raw = self._call_llm(prompt)
            raw = raw.replace("`json", "").replace("`", "").strip()
            result = json.loads(raw)
            category = result.get("category", "").upper()
            reason = result.get("reason", "")
            if category not in CATEGORY_TO_SCORE:
                return -1, ""
            return CATEGORY_TO_SCORE[category], reason
        except Exception as e:
            logger.debug(f"LLM file classification failed: {e}")
            return -1, ""

    def _call_llm(self, prompt: str) -> str:
        """Llama al LLM configurado y retorna el texto de respuesta."""
        if self.llm_backend == "openai":
            base_url = self.config.get("llm", {}).get("openai_base_url", "")
            api_version = self.config.get("llm", {}).get("openai_api_version", "")
            import os as _os; api_key = _os.getenv("OPENAI_API_KEY", self.config.get("llm", {}).get("api_key", ""))
            url = f"{base_url}/openai/deployments/{self.llm_model}/chat/completions?api-version={api_version}"
            r = requests.post(
                url,
                headers={"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"},
                json={"messages": [{"role": "user", "content": prompt}], "max_completion_tokens": 500, "response_format": {"type": "json_object"}},
                timeout=30,
            )
            r.raise_for_status()
            return r.json()["choices"][0]["message"]["content"].strip()
        else:
            r = requests.post(
                f"{self.ollama_url}/api/generate",
                json={"model": self.llm_model, "prompt": prompt, "stream": False, "format": "json"},
                timeout=30,
            )
            r.raise_for_status()
            return r.json().get("response", "").strip()

    def _check_case_control(self, gse_id: str, summary: str = "") -> bool:
        if not summary.strip():
            return True
        prompt = (
            "You are evaluating a GEO dataset for a transcriptomic meta-analysis pipeline.\n\n"
            f"Dataset: {gse_id}\n"
            f"Description: {summary.strip()[:600]}\n\n"
            "Does this dataset have at least TWO distinct biological groups that can be compared "
            "in a differential expression analysis (e.g. treated vs untreated, disease vs healthy, "
            "knockout vs wildtype, condition A vs condition B)?\n\n"
            "STRICT RULES:\n"
            "- Answer TRUE if the description explicitly mentions two or more distinct comparable groups.\n"
            "- Answer FALSE if all samples come from a single condition with no comparison group mentioned.\n"
            "- Answer FALSE if the description is too vague to determine any group structure.\n\n"
            "Respond ONLY with valid JSON: "
            '{"has_case_control": true, "reason": "<one sentence>"}'
        )
        try:
            raw = self._call_llm(prompt)
            result = json.loads(raw)
            has_cc = result.get("has_case_control", True)
            reason = result.get("reason", "")
            logger.debug(f"  [{gse_id}] caso/control={has_cc} ? {reason}")
            return has_cc
        except Exception as e:
            logger.debug(f"[{gse_id}] _check_case_control error: {e}")
            return True

    def _check_downloadability(self, gse_id: str, summary: str = "") -> tuple[int, list[str]]:
        """
        Verifica archivos suplementarios en FTP de NCBI y los clasifica con LLM.

        Scores:
            4 = count matrix crudo (listo para DESeq2/edgeR)
            3 = count matrix normalizado (FPKM, TPM, CPM — usable con limma)
            2 = series matrix u otros tabulares
            1 = solo RAW.tar (no usable sin infraestructura adicional)
            0 = sin archivos suplementarios

        Returns:
            (score, lista de nombres de archivos encontrados)
        """
        # Verificar diseno caso/control antes de evaluar archivos
        if not self._check_case_control(gse_id, summary=summary):
            logger.debug(f"  [{gse_id}] Sin diseno caso/control ? excluido")
            return 0, []

        files_found = []
        try:
            gse_num = int(gse_id.replace("GSE", ""))
            prefix = str(gse_num)[:-3] if gse_num >= 1000 else "0"
            ftp_url = (
                f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{prefix}nnn"
                f"/{gse_id}/suppl/"
            )
            resp = requests.get(ftp_url, timeout=10)
            if resp.status_code != 200:
                return 2, []

            # Extraer nombres de archivos del HTML del FTP
            exts = ('.gz', '.tar', '.csv', '.txt', '.tsv', '.xlsx', '.h5', '.loom')
            for part in resp.text.split('href='):
                chunk = part.split('"')
                if len(chunk) > 1:
                    fname = chunk[1]
                    if any(fname.lower().endswith(e) for e in exts) and fname not in ('..', ''):
                        files_found.append(fname)

            if not files_found:
                return 2, []

            # Pre-selección: si hay count crudo Y normalizado, quedarse solo con el count
            files_for_llm = files_found
            body = " ".join(files_found).lower()
            has_raw_count = any(k in body for k in ["gene_count", "rawcount", "raw_count",
                                                     "featurecount", "starcounts", "counts_"])
            has_normalized = any(k in body for k in ["normalized", "fpkm", "tpm", "cpm", "rpkm"])
            if has_raw_count and has_normalized:
                files_for_llm = [f for f in files_found if not any(
                    k in f.lower() for k in ["fpkm", "tpm", "cpm", "rpkm", "normalized"]
                )]

            # Clasificación semántica con LLM
            # Excluir RAW.tar antes de pasar al LLM ? nunca es el mejor archivo
            files_for_llm = [f for f in files_for_llm if not f.endswith("_RAW.tar") and f != "filelist.txt"]
            if not files_for_llm:
                return 1, files_found  # solo habia RAW.tar
            score, reason = self._classify_files_with_llm(files_for_llm)
            if score >= 0:
                logger.debug(f"  [{gse_id}] LLM score: {score} — {reason}")
                return score, files_found

            # Fallback: reglas deterministas si LLM falla
            logger.debug(f"  [{gse_id}] LLM falló — usando reglas deterministas")
            body = " ".join(files_found).lower()
            if any(k in body for k in ["normalized", "normalised", "fpkm", "tpm", "cpm", "rpkm"]):
                return 3, files_found
            if any(k in body for k in ["rawcount", "raw_count", "gene_count", "featurecount", "counts_"]):
                return 4, files_found
            if "_raw.tar" in body:
                return 1, files_found
            return 2, files_found

        except Exception as e:
            logger.debug(f"[{gse_id}] _check_downloadability error: {e}")
            return 2, []

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
        max_samples: int = 1000,
        data_types: list = None,
        output_dir: str = "data/discovery/",
        organism: str = None,
        year_start: int = None,
        year_end: int = None,
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
            max_samples=max_samples,
            data_types=data_types,
            organism=organism,
            year_start=year_start,
            year_end=year_end,
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
@click.option("--data-types", default="RNA-seq", show_default=True,
              help="Tipos de datos: RNA-seq | microarray | RNA-seq,microarray")
@click.option("--output", default="data/discovery/", show_default=True,
              help="Directorio de salida")
@click.option("--organism", default=None, show_default=True,
              help="Filtro por organismo: 'Homo sapiens', 'Mus musculus', etc.")
@click.option("--max-samples", default=1000, show_default=True,
              help="Máximo muestras por dataset")
@click.option("--year-start", default=None, type=int, show_default=True,
              help="Año inicio para filtro de publicación, ej: 2015")
@click.option("--year-end", default=None, type=int, show_default=True,
              help="Año fin para filtro de publicación, ej: 2026")
def main(disease, max_datasets, min_samples, data_types, output, organism, max_samples, year_start, year_end):
    """
    Agent 1: Descubrimiento de datasets GEO para una enfermedad dada.

    Ejemplo:
        python -m agents.discovery --disease 'pancreatic cancer'
    """
    from utils.geo_utils import load_config
    config = load_config()

    agent = DatasetDiscoveryAgent(config={**config.get("discovery", {}), "llm": config.get("llm", {})})
    agent.run(
        disease_name=disease,
        max_datasets=max_datasets,
        min_samples=min_samples,
        max_samples=max_samples,
        data_types=data_types.split(","),
        output_dir=output,
        organism=organism,
        year_start=year_start,
        year_end=year_end,
    )


if __name__ == "__main__":
    main()
