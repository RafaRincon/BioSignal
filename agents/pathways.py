"""
BioSignal Discovery Engine
Agent 6: PathwayEnrichmentAgent
=================================
Responsabilidad: Identificar pathways biológicos y procesos celulares
significativamente alterados en los genes consenso del meta-análisis.

Usa la API de Enrichr (gseapy) para consultar KEGG, GO, Reactome y MSigDB Hallmarks.

Uso:
    python -m agents.pathways --input data/meta/ --output data/pathways/
"""

import json
import time
from pathlib import Path
from typing import Optional

import click
import numpy as np
import pandas as pd
from loguru import logger

try:
    import gseapy as gp
    GSEAPY_AVAILABLE = True
except ImportError:
    GSEAPY_AVAILABLE = False
    logger.warning("gseapy no disponible — enriquecimiento de pathways deshabilitado")


# Bases de datos de pathways a consultar
PATHWAY_DATABASES = {
    "kegg": "KEGG_2021_Human",
    "go_bp": "GO_Biological_Process_2023",
    "reactome": "Reactome_2022",
    "hallmarks": "MSigDB_Hallmark_2020",
}


class PathwayEnrichmentAgent:
    """
    Agente de enriquecimiento de pathways biológicos.

    Proceso:
        1. Cargar consensus_genes.json del Agent 5 (meta-análisis)
        2. Separar genes UP y DOWN regulados
        3. Ejecutar Over-Representation Analysis (ORA) via Enrichr para cada lista
        4. Filtrar: padj < 0.05 y gene_ratio >= 0.05
        5. Rankear por combined_score = -log10(pval) × z-score
        6. Exportar resumen top pathways

    Salida:
        data/pathways/{disease}/kegg_enrichment.csv
        data/pathways/{disease}/go_enrichment.csv
        data/pathways/{disease}/reactome_enrichment.csv
        data/pathways/{disease}/hallmarks_enrichment.csv
        data/pathways/{disease}/top_pathways_summary.json
        data/pathways/{disease}/enrichment_status.json
    """

    def __init__(self, config: dict = None):
        self.config = config or {}
        self.padj_threshold = self.config.get("padj_threshold", 0.05)
        self.min_gene_ratio = self.config.get("min_gene_ratio", 0.05)
        self.top_n_pathways = self.config.get("top_n_pathways", 20)
        self.enrichr_timeout = self.config.get("enrichr_timeout", 30)
        self.organism = self.config.get("organism", "human")

    # ------------------------------------------------------------------
    # Entrada principal
    # ------------------------------------------------------------------

    def run(self, meta_dir: str, output_dir: str) -> dict:
        """
        Ejecuta enriquecimiento de pathways sobre resultados del meta-análisis.

        Args:
            meta_dir:   Directorio de salida del Agent 5 (data/meta/)
            output_dir: Directorio destino (data/pathways/)

        Returns:
            Resumen por enfermedad analizada
        """
        meta_path = Path(meta_dir)
        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        summary = {"enriched": [], "failed": [], "skipped": []}
        start = time.time()

        if not GSEAPY_AVAILABLE:
            logger.error("gseapy no instalado. Instalar con: pip install gseapy")
            return {"error": "gseapy_not_available"}

        disease_dirs = sorted([d for d in meta_path.iterdir() if d.is_dir()])
        if self.config.get("disease_filter"):
            disease_dirs = [d for d in disease_dirs if d.name == self.config["disease_filter"]]
        logger.info(f"[PathwayAgent] Encontrados {len(disease_dirs)} resultados de meta-análisis")

        for disease_dir in disease_dirs:
            disease_name = disease_dir.name
            try:
                result = self._enrich_disease(disease_dir, out_path / disease_name)
                summary["enriched"].append({
                    "disease": disease_name,
                    "n_pathways_significant": result.get("n_significant_total", 0),
                })
            except Exception as e:
                logger.error(f"[{disease_name}] Error en pathway enrichment: {e}")
                summary["failed"].append({"disease": disease_name, "error": str(e)})

        elapsed = time.time() - start
        logger.info(f"[PathwayAgent] Completado en {elapsed:.1f}s")
        return summary

    # ------------------------------------------------------------------
    # Enriquecimiento por enfermedad
    # ------------------------------------------------------------------

    def _enrich_disease(self, disease_dir: Path, out_dir: Path) -> dict:
        """Ejecuta el enriquecimiento para una enfermedad."""
        out_dir.mkdir(parents=True, exist_ok=True)
        disease_name = disease_dir.name

        # Cargar genes consenso
        consensus_path = disease_dir / "consensus_genes.json"
        if not consensus_path.exists():
            raise FileNotFoundError(f"consensus_genes.json no encontrado en {disease_dir}")

        with open(consensus_path) as f:
            consensus = json.load(f)

        genes_up = [g["gene"] for g in consensus.get("consensus_genes", []) if g.get("direction") == "UP"]
        genes_down = [g["gene"] for g in consensus.get("consensus_genes", []) if g.get("direction") == "DOWN"]
        all_genes = genes_up + genes_down

        logger.info(f"[{disease_name}] Genes UP={len(genes_up)}, DOWN={len(genes_down)}")

        if len(all_genes) < 5:
            logger.warning(f"[{disease_name}] Muy pocos genes ({len(all_genes)}) para enriquecimiento significativo")

        all_results = {}
        n_significant_total = 0

        for db_key, db_name in PATHWAY_DATABASES.items():
            try:
                logger.info(f"[{disease_name}] Consultando {db_name}...")
                enr_df = self._run_enrichr(all_genes, db_name, disease_name)

                if enr_df is not None and not enr_df.empty:
                    # Filtrar significativos
                    sig_df = enr_df[enr_df["Adjusted P-value"] < self.padj_threshold].copy()

                    # Calcular gene_ratio
                    sig_df["gene_ratio"] = sig_df["Overlap"].apply(self._parse_overlap_ratio)
                    sig_df = sig_df[sig_df["gene_ratio"] >= self.min_gene_ratio]

                    # Calcular combined_score
                    sig_df["combined_score"] = sig_df.apply(
                        lambda r: self._compute_combined_score(r), axis=1
                    )
                    sig_df = sig_df.sort_values("combined_score", ascending=False)

                    # Guardar CSV
                    out_csv = out_dir / f"{db_key}_enrichment.csv"
                    sig_df.to_csv(out_csv, index=False)
                    n_significant_total += len(sig_df)
                    all_results[db_key] = sig_df.head(self.top_n_pathways).to_dict(orient="records")
                    logger.info(f"[{disease_name}] {db_name}: {len(sig_df)} pathways significativos")
                else:
                    logger.warning(f"[{disease_name}] Sin resultados de {db_name}")

            except Exception as e:
                logger.warning(f"[{disease_name}] Error en {db_name}: {e}")

            time.sleep(0.5)  # Rate limiting educado con la API de Enrichr

        # Construir resumen top pathways
        top_summary = self._build_top_summary(all_results, disease_name)
        top_summary["disease"] = disease_name
        top_summary["genes_analyzed"] = len(all_genes)
        top_summary["n_significant_total"] = n_significant_total

        with open(out_dir / "top_pathways_summary.json", "w") as f:
            json.dump(top_summary, f, indent=2, default=str)

        status = {"disease": disease_name, "n_significant_total": n_significant_total, "databases_queried": list(PATHWAY_DATABASES.keys())}
        with open(out_dir / "enrichment_status.json", "w") as f:
            json.dump(status, f, indent=2)

        logger.info(f"[{disease_name}] Enriquecimiento completo | Total significativos: {n_significant_total}")
        return status

    # ------------------------------------------------------------------
    # Enrichr via gseapy
    # ------------------------------------------------------------------

    def _run_enrichr(self, gene_list: list, db_name: str, disease_name: str) -> Optional[pd.DataFrame]:
        """Llama a la API de Enrichr con reintentos."""
        max_retries = 3
        for attempt in range(max_retries):
            try:
                enr = gp.enrichr(
                    gene_list=gene_list,
                    gene_sets=db_name,
                    organism=self.organism,
                    outdir=None,
                    no_plot=True,
                    verbose=False,
                )
                if hasattr(enr, "results") and enr.results is not None:
                    return enr.results
                return None
            except Exception as e:
                if attempt < max_retries - 1:
                    wait = 2 ** attempt
                    logger.warning(f"[{disease_name}] Enrichr intento {attempt+1} falló: {e}. Reintentando en {wait}s")
                    time.sleep(wait)
                else:
                    raise

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _parse_overlap_ratio(self, overlap_str: str) -> float:
        """Convierte '5/100' a 0.05."""
        try:
            parts = str(overlap_str).split("/")
            if len(parts) == 2:
                return int(parts[0]) / int(parts[1])
        except Exception:
            pass
        return 0.0

    def _compute_combined_score(self, row) -> float:
        """Combined score = -log10(pval) × |z-score| (similar a Enrichr interno)."""
        try:
            pval = max(row.get("P-value", 1e-300), 1e-300)
            z = abs(row.get("Z-score", 0) or 0)
            return -np.log10(pval) * z
        except Exception:
            return 0.0

    def _build_top_summary(self, all_results: dict, disease_name: str) -> dict:
        """Construye resumen consolidado de top pathways por base de datos."""
        summary = {"top_pathways_by_db": {}}

        for db_key, records in all_results.items():
            if not records:
                continue
            top_paths = []
            for r in records[:10]:  # top 10 por DB
                top_paths.append({
                    "pathway": r.get("Term", ""),
                    "padj": round(r.get("Adjusted P-value", 1), 6),
                    "combined_score": round(r.get("combined_score", 0), 3),
                    "genes": r.get("Genes", "").split(";")[:10] if r.get("Genes") else [],
                })
            summary["top_pathways_by_db"][db_key] = top_paths

        return summary


# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------

@click.command()
@click.option("--input", "meta_dir", required=True, help="Directorio con resultados del meta-análisis (data/meta/)")
@click.option("--output", "output_dir", default="data/pathways/", help="Directorio de salida")
@click.option("--padj", default=0.05, help="Umbral FDR para filtrado")
@click.option("--top-n", default=20, help="Top N pathways a reportar")
@click.option("--disease", default=None, help="Filtrar por enfermedad especifica (ej. alzheimer)")
def main(meta_dir, output_dir, padj, top_n, disease):
    """Agent 6: Enriquecimiento de pathways (KEGG, GO, Reactome, Hallmarks)."""
    config = {"padj_threshold": padj, "top_n_pathways": top_n, "disease_filter": disease}
    agent = PathwayEnrichmentAgent(config=config)
    summary = agent.run(meta_dir, output_dir)
    click.echo(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
