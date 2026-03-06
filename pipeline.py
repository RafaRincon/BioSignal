"""
BioSignal Discovery Engine
Pipeline Orchestrator
======================
Orquesta la ejecución secuencial de los 8 agentes del sistema.

Uso:
    from biosignal.pipeline import BioSignalPipeline

    pipeline = BioSignalPipeline()
    pipeline.run(
        disease="pancreatic cancer",
        max_datasets=15,
        output_dir="data/",
    )
"""

import json
import time
from pathlib import Path

from loguru import logger

from agents.discovery import DatasetDiscoveryAgent
from agents.download import DatasetDownloadAgent
from agents.meta_analysis import MetaAnalysisAgent
from utils.geo_utils import load_config


class BioSignalPipeline:
    """
    Orquestador del pipeline completo de BioSignal Discovery Engine.

    Ejecuta los 8 agentes en secuencia:
        Agent 1: DatasetDiscoveryAgent   — Identifica datasets GEO
        Agent 2: DatasetDownloadAgent    — Descarga datos crudos
        Agent 3: PreprocessingAgent      — Normaliza y QC
        Agent 4: DifferentialExpressionAgent — DEA por dataset
        Agent 5: MetaAnalysisAgent       — Meta-análisis cross-dataset
        Agent 6: PathwayEnrichmentAgent  — Enriquecimiento de pathways
        Agent 7: InsightGenerationAgent  — Hipótesis con LLM
        Agent 8: ReportGenerationAgent   — Reporte final

    Principios de diseño:
        - Modularidad: cada agente puede re-ejecutarse independientemente
        - Idempotencia: re-ejecución sin efectos secundarios
        - Fail-safe: si un dataset falla, el pipeline continúa
        - Reproducibilidad: seeds fijos, versiones pinadas
    """

    def __init__(self, config_path: str = "config/settings.yaml"):
        self.config = load_config(config_path)
        self.pipeline_config = self.config.get("pipeline", {})

        # Configurar logging
        log_level = self.pipeline_config.get("log_level", "INFO")
        logger.remove()
        logger.add(
            lambda msg: print(msg, end=""),
            level=log_level,
            format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | {message}",
            colorize=True,
        )
        logger.add(
            "logs/pipeline_{time}.log",
            level="DEBUG",
            rotation="100 MB",
            retention="30 days",
        )

    def run(
        self,
        disease: str,
        max_datasets: int = 15,
        min_samples: int = 10,
        parallel: int = 4,
        output_dir: str = "data/",
        report_format: str = "pdf,json,csv",
        start_from_agent: int = 1,
        stop_after_agent: int = 8,
    ) -> dict:
        """
        Ejecuta el pipeline completo para una enfermedad dada.

        Args:
            disease: Nombre de la enfermedad a analizar
            max_datasets: Máximo datasets a descubrir/analizar
            min_samples: Mínimo de muestras por dataset
            parallel: Número de workers para descarga paralela
            output_dir: Directorio base para todos los outputs
            report_format: Formatos de reporte separados por coma
            start_from_agent: Número de agente desde donde iniciar (1-8)
            stop_after_agent: Número de agente donde detener (1-8)

        Returns:
            Dict con resultados y paths de cada agente
        """
        logger.info("=" * 60)
        logger.info("🧬 BioSignal Discovery Engine v1.0.0")
        logger.info("=" * 60)
        logger.info(f"Enfermedad: '{disease}'")
        logger.info(f"Parámetros: max_datasets={max_datasets}, min_samples={min_samples}")
        logger.info(f"Agentes a ejecutar: {start_from_agent} → {stop_after_agent}")

        pipeline_start = time.time()
        results = {"disease": disease, "agents": {}}

        disease_safe = disease.replace(" ", "_").lower()
        base_dir = Path(output_dir)

        # ============================================================
        # Agent 1: DatasetDiscoveryAgent
        # ============================================================
        if start_from_agent <= 1 <= stop_after_agent:
            logger.info("\n" + "─" * 40)
            t0 = time.time()
            discovery_agent = DatasetDiscoveryAgent(
                config=self.config.get("discovery", {})
            )
            discovery_output = discovery_agent.run(
                disease_name=disease,
                max_datasets=max_datasets,
                min_samples=min_samples,
                output_dir=str(base_dir / "discovery"),
            )
            results["agents"]["1_discovery"] = {
                "output": str(discovery_output),
                "elapsed_s": round(time.time() - t0, 1),
                "status": "success",
            }

        # ============================================================
        # Agent 2: DatasetDownloadAgent
        # ============================================================
        if start_from_agent <= 2 <= stop_after_agent:
            logger.info("\n" + "─" * 40)
            t0 = time.time()

            discovery_json = base_dir / "discovery" / f"datasets_{disease_safe}.json"
            if not discovery_json.exists():
                logger.error(f"No encontrado output del Agent 1: {discovery_json}")
                results["agents"]["2_download"] = {"status": "skipped_missing_input"}
            else:
                download_agent = DatasetDownloadAgent(
                    config=self.config.get("download", {})
                )
                download_results = download_agent.run_parallel(
                    input_json=str(discovery_json),
                    output_dir=str(base_dir / "raw"),
                    parallel=parallel,
                )
                results["agents"]["2_download"] = {
                    **download_results,
                    "elapsed_s": round(time.time() - t0, 1),
                    "status": "success",
                }

        # ============================================================
        # Agent 3-4: Preprocessing + DEA (placeholder - importar cuando estén listos)
        # ============================================================
        if start_from_agent <= 3 <= stop_after_agent:
            logger.info("\n" + "─" * 40)
            logger.info("[Agent 3] PreprocessingAgent — pendiente de implementación completa")
            # TODO: from biosignal.agents.preprocess import PreprocessingAgent
            results["agents"]["3_preprocessing"] = {"status": "pending"}

        if start_from_agent <= 4 <= stop_after_agent:
            logger.info("\n" + "─" * 40)
            logger.info("[Agent 4] DifferentialExpressionAgent — pendiente de implementación completa")
            # TODO: from biosignal.agents.dea import DifferentialExpressionAgent
            results["agents"]["4_dea"] = {"status": "pending"}

        # ============================================================
        # Agent 5: MetaAnalysisAgent
        # ============================================================
        if start_from_agent <= 5 <= stop_after_agent:
            logger.info("\n" + "─" * 40)
            dea_dir = base_dir / "dea"
            if not dea_dir.exists() or not any(dea_dir.iterdir()):
                logger.warning(f"No hay resultados DEA en {dea_dir}. Saltando Agent 5.")
                results["agents"]["5_meta_analysis"] = {"status": "skipped_no_dea"}
            else:
                t0 = time.time()
                meta_agent = MetaAnalysisAgent(
                    config=self.config.get("meta_analysis", {})
                )
                meta_output = meta_agent.run(
                    input_dir=str(dea_dir),
                    output_dir=str(base_dir / "meta"),
                    disease_name=disease,
                )
                results["agents"]["5_meta_analysis"] = {
                    "output": str(meta_output),
                    "elapsed_s": round(time.time() - t0, 1),
                    "status": "success",
                }

        # ============================================================
        # Agent 6-8: Pathways, Insights, Report (placeholders)
        # ============================================================
        for agent_num, agent_name in [
            (6, "PathwayEnrichmentAgent"),
            (7, "InsightGenerationAgent"),
            (8, "ReportGenerationAgent"),
        ]:
            if start_from_agent <= agent_num <= stop_after_agent:
                logger.info("\n" + "─" * 40)
                logger.info(f"[Agent {agent_num}] {agent_name} — pendiente de implementación completa")
                results["agents"][f"{agent_num}_{agent_name.lower()}"] = {"status": "pending"}

        # ============================================================
        # Resumen final
        # ============================================================
        total_elapsed = time.time() - pipeline_start
        results["total_elapsed_s"] = round(total_elapsed, 1)

        logger.info("\n" + "=" * 60)
        logger.success("✓ Pipeline completado")
        logger.info(f"  Enfermedad: {disease}")
        logger.info(f"  Tiempo total: {total_elapsed/60:.1f} minutos")
        logger.info(f"  Outputs en: {output_dir}")

        # Guardar resumen del run
        summary_path = base_dir / f"pipeline_summary_{disease_safe}.json"
        with open(summary_path, "w") as f:
            json.dump(results, f, indent=2, default=str)
        logger.info(f"  Resumen guardado: {summary_path}")
        logger.info("=" * 60)

        return results