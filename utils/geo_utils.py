"""
BioSignal Discovery Engine
Utilities: GEO/NCBI Access Functions
=====================================
Funciones de utilidad para acceder a GEO y NCBI E-utilities.
"""

import json
import os
import time
from pathlib import Path
from typing import Optional

import requests
import yaml
from loguru import logger


def load_config(config_path: str = "config/settings.yaml") -> dict:
    """
    Carga la configuración global del proyecto.

    Args:
        config_path: Path al archivo settings.yaml

    Returns:
        Dict de configuración
    """
    config_file = Path(config_path)
    if not config_file.exists():
        logger.warning(f"Config no encontrada en {config_path}. Usando defaults.")
        return {}
    with open(config_file, encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def load_disease_aliases(
    aliases_path: str = "config/disease_aliases.json",
) -> dict:
    """
    Carga el mapeo de nombres de enfermedades a términos GEO.

    Args:
        aliases_path: Path al archivo de aliases

    Returns:
        Dict de aliases por categoría
    """
    aliases_file = Path(aliases_path)
    if not aliases_file.exists():
        return {}
    with open(aliases_file) as f:
        return json.load(f)


class NCBIClient:
    """
    Cliente para la API de NCBI E-utilities.

    Endpoints utilizados:
        - esearch: Búsqueda de IDs en bases de datos NCBI
        - esummary: Obtener metadata de IDs
        - efetch: Descargar datos completos

    Límites de uso:
        - Sin API key: máximo 3 requests/seg
        - Con API key: máximo 10 requests/seg
    """

    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    def __init__(self, email: str = "", api_key: str = ""):
        self.email = email or os.getenv("NCBI_EMAIL", "biosignal@example.com")
        self.api_key = api_key or os.getenv("NCBI_API_KEY", "")
        self._rate_limit_delay = 0.11 if self.api_key else 0.34  # 9/seg vs 3/seg

    def _build_params(self, **kwargs) -> dict:
        """Construye parámetros base para todas las requests."""
        params = {
            "email": self.email,
            "tool": "biosignal-discovery",
            "retmode": "json",
        }
        if self.api_key:
            params["api_key"] = self.api_key
        params.update(kwargs)
        return params

    def _get(self, endpoint: str, params: dict, retries: int = 3) -> dict:
        """
        Realiza una request GET con rate limiting y reintentos.

        Args:
            endpoint: Nombre del endpoint (esearch, esummary, etc.)
            params: Parámetros de la request
            retries: Número de reintentos en caso de error

        Returns:
            JSON response como dict
        """
        url = f"{self.BASE_URL}{endpoint}.fcgi"
        time.sleep(self._rate_limit_delay)  # Rate limiting

        for attempt in range(retries):
            try:
                response = requests.get(url, params=params, timeout=30)
                response.raise_for_status()
                return response.json()
            except requests.exceptions.RequestException as e:
                if attempt < retries - 1:
                    wait = 2 ** attempt
                    logger.debug(f"NCBI request failed ({e}). Retrying in {wait}s...")
                    time.sleep(wait)
                else:
                    raise

        return {}

    def esearch(
        self,
        db: str,
        term: str,
        retmax: int = 100,
        sort: str = "relevance",
    ) -> list[str]:
        """
        Busca en una base de datos NCBI y retorna lista de IDs.

        Args:
            db: Base de datos (gds, pubmed, gene, etc.)
            term: Término de búsqueda
            retmax: Máximo número de resultados
            sort: Orden de resultados

        Returns:
            Lista de IDs como strings
        """
        params = self._build_params(
            db=db,
            term=term,
            retmax=retmax,
            sort=sort,
            usehistory="y",
        )

        try:
            result = self._get("esearch", params)
            ids = result.get("esearchresult", {}).get("idlist", [])
            logger.debug(f"esearch '{term}': {len(ids)} IDs encontrados")
            return ids
        except Exception as e:
            logger.error(f"esearch falló para '{term}': {e}")
            return []

    def esummary(
        self,
        db: str,
        ids: list[str],
        batch_size: int = 100,
    ) -> list[dict]:
        """
        Obtiene metadata/summary para una lista de IDs.

        Args:
            db: Base de datos NCBI
            ids: Lista de IDs a consultar
            batch_size: Tamaño de lote para queries masivas

        Returns:
            Lista de dicts con metadata de cada ID
        """
        all_results = []

        # Procesar en lotes para evitar URLs muy largas
        for i in range(0, len(ids), batch_size):
            batch = ids[i:i + batch_size]
            params = self._build_params(
                db=db,
                id=",".join(batch),
            )

            try:
                result = self._get("esummary", params)
                doc_set = result.get("result", {})

                for uid in doc_set.get("uids", []):
                    if uid in doc_set:
                        all_results.append(doc_set[uid])

                logger.debug(f"esummary batch {i//batch_size + 1}: {len(batch)} IDs")

            except Exception as e:
                logger.error(f"esummary falló para batch {i//batch_size + 1}: {e}")
                continue

        return all_results

    def efetch(
        self,
        db: str,
        ids: list[str],
        rettype: str = "text",
    ) -> str:
        """
        Descarga datos completos para una lista de IDs.

        Args:
            db: Base de datos NCBI
            ids: Lista de IDs
            rettype: Tipo de retorno (text, fasta, etc.)

        Returns:
            Contenido como string
        """
        params = self._build_params(
            db=db,
            id=",".join(ids),
            rettype=rettype,
        )
        params.pop("retmode")  # efetch no usa retmode=json

        response = requests.get(
            f"{self.BASE_URL}efetch.fcgi",
            params=params,
            timeout=60,
        )
        response.raise_for_status()
        return response.text