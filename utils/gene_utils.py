"""
BioSignal Discovery Engine
Utils: Gene ID Mapping and Normalization
==========================================
Funciones para mapeo, normalización y validación de identificadores de genes.
Soporta conversión entre Ensembl, Entrez, Probe IDs y símbolos HGNC.
"""

import re
import time
from typing import Optional
from pathlib import Path

import pandas as pd
import requests
from loguru import logger

# ------------------------------------------------------------------
# Normalización de símbolos génicos
# ------------------------------------------------------------------

def normalize_gene_symbol(symbol: str) -> str:
    """
    Normaliza un símbolo génico a formato HGNC estándar.
    Elimina versiones, sufijos de plataforma y caracteres especiales.

    Examples:
        "TP53_at"   → "TP53"
        "GAPDH.1"   → "GAPDH"
        "hsa-mir-21"→ "MIR21"
    """
    if not symbol or not isinstance(symbol, str):
        return ""

    s = symbol.strip()

    # Eliminar sufijos de probes de microarray: _at, _s_at, _x_at
    s = re.sub(r'_[sx]?_?at$', '', s, flags=re.IGNORECASE)

    # Eliminar versiones con punto o guión bajo: GENE.1, GENE_1
    s = re.sub(r'[._]\d+$', '', s)

    # Eliminar prefijos de microRNA Entrez
    s = re.sub(r'^hsa-mir-', 'MIR', s, flags=re.IGNORECASE)
    s = re.sub(r'^hsa-miR-', 'MIR', s)

    # Convertir a mayúsculas (convención HGNC)
    s = s.upper()

    # Eliminar caracteres no válidos
    s = re.sub(r'[^A-Z0-9\-]', '', s)

    return s.strip()


def normalize_gene_list(genes: list) -> list:
    """Normaliza una lista de símbolos génicos."""
    normalized = [normalize_gene_symbol(g) for g in genes]
    return [g for g in normalized if g]  # eliminar vacíos


# ------------------------------------------------------------------
# Filtros de genes de referencia (housekeeping, pseudogenes, etc.)
# ------------------------------------------------------------------

# Genes housekeeping comúnmente excluidos del análisis diferencial
HOUSEKEEPING_GENES = frozenset([
    "GAPDH", "ACTB", "B2M", "HPRT1", "RPS18", "RPLP0",
    "PPIA", "TFRC", "TBP", "GUSB", "HMBS", "SDHA",
])

# Prefijos de pseudogenes y genes no codificantes a filtrar opcionalmente
PSEUDO_PREFIXES = ("LOC", "LINC", "SNORD", "SNRNA", "MIR")


def filter_housekeeping(gene_list: list, remove_housekeeping: bool = False,
                         remove_pseudogenes: bool = False) -> list:
    """
    Filtra genes de referencia y pseudogenes de una lista.

    Args:
        gene_list:          Lista de símbolos HGNC
        remove_housekeeping: Si True, elimina genes housekeeping clásicos
        remove_pseudogenes:  Si True, elimina LOC*, LINC*, SNORD*, etc.

    Returns:
        Lista filtrada
    """
    result = gene_list.copy()

    if remove_housekeeping:
        before = len(result)
        result = [g for g in result if g not in HOUSEKEEPING_GENES]
        logger.debug(f"[gene_utils] Eliminados {before - len(result)} genes housekeeping")

    if remove_pseudogenes:
        before = len(result)
        result = [g for g in result if not any(g.startswith(p) for p in PSEUDO_PREFIXES)]
        logger.debug(f"[gene_utils] Eliminados {before - len(result)} pseudogenes/ncRNA")

    return result


# ------------------------------------------------------------------
# Mapeo de IDs via MyGene.info API
# ------------------------------------------------------------------

def map_entrez_to_symbol(entrez_ids: list, species: str = "human") -> dict:
    """
    Convierte Entrez Gene IDs a símbolos HGNC via MyGene.info.

    Args:
        entrez_ids: Lista de Entrez IDs (como strings o ints)
        species:    Organismo (default: 'human')

    Returns:
        Dict {entrez_id: symbol}
    """
    if not entrez_ids:
        return {}

    mapping = {}
    chunk_size = 1000
    base_url = "https://mygene.info/v3/gene"

    for i in range(0, len(entrez_ids), chunk_size):
        chunk = entrez_ids[i:i + chunk_size]
        ids_str = ",".join(str(e) for e in chunk)

        try:
            r = requests.post(
                base_url,
                data={"ids": ids_str, "fields": "symbol,alias", "species": species},
                timeout=20,
            )
            r.raise_for_status()
            results = r.json()

            for entry in results:
                if "notfound" in entry:
                    continue
                entrez = str(entry.get("_id", ""))
                symbol = entry.get("symbol", "")
                if entrez and symbol:
                    mapping[entrez] = symbol

            time.sleep(0.2)

        except Exception as e:
            logger.warning(f"[gene_utils] Error MyGene.info chunk {i}: {e}")

    logger.info(f"[gene_utils] Mapeados {len(mapping)}/{len(entrez_ids)} Entrez IDs a símbolos")
    return mapping


def map_ensembl_to_symbol(ensembl_ids: list, species: str = "human") -> dict:
    """
    Convierte Ensembl Gene IDs (ENSG...) a símbolos HGNC via MyGene.info.

    Args:
        ensembl_ids: Lista de IDs tipo 'ENSG00000141510'
        species:     Organismo

    Returns:
        Dict {ensembl_id: symbol}
    """
    if not ensembl_ids:
        return {}

    mapping = {}
    chunk_size = 500
    base_url = "https://mygene.info/v3/query"

    for i in range(0, len(ensembl_ids), chunk_size):
        chunk = ensembl_ids[i:i + chunk_size]

        try:
            r = requests.post(
                base_url,
                data={
                    "q": ",".join(chunk),
                    "scopes": "ensembl.gene",
                    "fields": "symbol,ensembl.gene",
                    "species": species,
                    "size": chunk_size,
                },
                timeout=20,
            )
            r.raise_for_status()
            results = r.json()

            for entry in results.get("hits", []):
                if entry.get("notfound"):
                    continue
                ensembl = entry.get("query", "")
                symbol = entry.get("symbol", "")
                if ensembl and symbol:
                    mapping[ensembl] = symbol

            time.sleep(0.2)

        except Exception as e:
            logger.warning(f"[gene_utils] Error Ensembl mapeo chunk {i}: {e}")

    logger.info(f"[gene_utils] Mapeados {len(mapping)}/{len(ensembl_ids)} Ensembl IDs a símbolos")
    return mapping


# ------------------------------------------------------------------
# Detección automática del tipo de ID
# ------------------------------------------------------------------

def detect_gene_id_type(ids: list) -> str:
    """
    Detecta el tipo de identificador génico en una lista.

    Returns:
        'ensembl' | 'entrez' | 'symbol' | 'probe' | 'unknown'
    """
    sample = [str(x) for x in ids[:50] if x]

    ensembl_count = sum(1 for s in sample if s.startswith("ENSG") or s.startswith("ENSMUSG"))
    entrez_count = sum(1 for s in sample if s.isdigit())
    probe_count = sum(1 for s in sample if re.search(r'_at$|_x_at$|_s_at$', s, re.I))

    if ensembl_count / len(sample) > 0.5:
        return "ensembl"
    if entrez_count / len(sample) > 0.5:
        return "entrez"
    if probe_count / len(sample) > 0.3:
        return "probe"
    return "symbol"


def auto_map_to_symbol(df: pd.DataFrame, gene_col: str = None, species: str = "human") -> pd.DataFrame:
    """
    Detecta automáticamente el tipo de ID génico y convierte a símbolos HGNC.

    Args:
        df:       DataFrame con columna de genes/IDs
        gene_col: Nombre de la columna de genes (si None, usa el índice)
        species:  Organismo

    Returns:
        DataFrame con columna 'gene_symbol' añadida
    """
    if gene_col is None:
        ids = df.index.tolist()
    else:
        ids = df[gene_col].tolist()

    id_type = detect_gene_id_type(ids)
    logger.info(f"[gene_utils] Tipo de ID detectado: {id_type}")

    if id_type == "symbol":
        # Normalizar directamente
        symbols = [normalize_gene_symbol(str(i)) for i in ids]
        mapping = {str(original): normalized for original, normalized in zip(ids, symbols)}

    elif id_type == "ensembl":
        mapping = map_ensembl_to_symbol(ids, species)

    elif id_type == "entrez":
        mapping = map_entrez_to_symbol(ids, species)

    elif id_type == "probe":
        # Intentar extraer símbolo del probe name (microarray Affymetrix)
        mapping = {}
        for probe in ids:
            clean = normalize_gene_symbol(str(probe))
            if clean:
                mapping[str(probe)] = clean

    else:
        logger.warning("[gene_utils] Tipo de ID desconocido — sin mapeo disponible")
        mapping = {}

    # Aplicar mapeo
    if gene_col is None:
        df = df.copy()
        df["gene_symbol"] = [mapping.get(str(i), str(i)) for i in ids]
    else:
        df = df.copy()
        df["gene_symbol"] = df[gene_col].map(lambda x: mapping.get(str(x), str(x)))

    mapped = df["gene_symbol"].notna().sum()
    logger.info(f"[gene_utils] Genes mapeados: {mapped}/{len(df)}")
    return df


# ------------------------------------------------------------------
# Colapso de genes duplicados
# ------------------------------------------------------------------

def collapse_duplicate_genes(df: pd.DataFrame, gene_col: str = "gene_symbol",
                               value_cols: list = None, method: str = "mean") -> pd.DataFrame:
    """
    Colapsa filas con el mismo símbolo génico.

    Args:
        df:         DataFrame con posibles duplicados
        gene_col:   Columna que contiene los símbolos
        value_cols: Columnas numéricas a colapsar (si None, todas las numéricas)
        method:     'mean' | 'max' | 'first' (método de colapso)

    Returns:
        DataFrame sin duplicados
    """
    if gene_col not in df.columns:
        return df

    if value_cols is None:
        value_cols = df.select_dtypes(include="number").columns.tolist()

    n_before = len(df)
    n_dupes = df.duplicated(subset=[gene_col]).sum()

    if n_dupes == 0:
        return df

    if method == "mean":
        df_out = df.groupby(gene_col)[value_cols].mean().reset_index()
    elif method == "max":
        df_out = df.groupby(gene_col)[value_cols].max().reset_index()
    elif method == "first":
        df_out = df.drop_duplicates(subset=[gene_col], keep="first")
    else:
        raise ValueError(f"Método desconocido: {method}")

    logger.debug(f"[gene_utils] Colapsados {n_dupes} duplicados: {n_before} → {len(df_out)} genes únicos")
    return df_out