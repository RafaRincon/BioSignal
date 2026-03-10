"""
BioSignal Discovery Engine
Audit Trail — AuditWriter
==========================
Módulo compartido por todos los agentes para registrar decisiones,
criterios y evidencia de forma estandarizada por run.

Cada agente instancia un AuditWriter al inicio de su ejecución y llama:
    writer.add_decision(...)  — por cada decisión tomada
    writer.add_warning(...)   — por cada condición borderline
    writer.add_error(...)     — por cada fallo no-crítico
    writer.save()             — al finalizar el agente

Salida:
    data/runs/{run_id}/audit/{agent_id}.json

Uso típico:
    from utils.audit import AuditWriter

    writer = AuditWriter(
        agent_id="03_preprocess",
        agent_name="PreprocessingAgent",
        run_id=run_id,
        run_dir=run_dir,
        config_snapshot=config,
        data_type="transcriptomics",
    )
    writer.add_decision(
        subject="GSE226507",
        decision="ACCEPT",
        reason="Todos los criterios de calidad superados",
        evidence={
            "samples_per_group": {"observed": 8, "threshold": 3, "pass": True},
        },
    )
    writer.save()
"""

import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from loguru import logger


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _utcnow() -> str:
    """ISO-8601 timestamp en UTC."""
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _config_hash(config: dict) -> str:
    """SHA-256 corto del snapshot de configuración (primeros 16 chars)."""
    raw = json.dumps(config, sort_keys=True, default=str).encode()
    return "sha256:" + hashlib.sha256(raw).hexdigest()[:16]


# ---------------------------------------------------------------------------
# AuditWriter
# ---------------------------------------------------------------------------

class AuditWriter:
    """
    Escritor de trazabilidad para un agente en un run específico.

    Genera un archivo JSON en:
        {run_dir}/audit/{agent_id}.json

    El archivo sigue el esquema AuditEntry documentado en
    docs/AUDIT_TRAIL_ARCHITECTURE.md
    """

    def __init__(
        self,
        agent_id: str,
        agent_name: str,
        run_id: str,
        run_dir: Path,
        config_snapshot: dict,
        data_type: str = "transcriptomics",
    ):
        self.agent_id = agent_id
        self.agent_name = agent_name
        self.run_id = run_id
        self.run_dir = Path(run_dir)
        self.data_type = data_type

        self._started_at = _utcnow()
        self._start_ts = time.time()

        # Guardar snapshot de configuración relevante para este agente
        self._config_snapshot = config_snapshot or {}
        self._config_hash = _config_hash(self._config_snapshot)

        self._decisions: list[dict] = []
        self._warnings: list[str] = []
        self._errors: list[str] = []

        # Crear directorio de auditoría
        self._audit_dir = self.run_dir / "audit"
        self._audit_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # API pública
    # ------------------------------------------------------------------

    def add_decision(
        self,
        subject: str,
        decision: str,
        reason: str,
        evidence: Optional[dict[str, Any]] = None,
    ) -> None:
        """
        Registra una decisión tomada sobre un sujeto (dataset, gen, pathway, etc.).

        Args:
            subject:  Identificador del sujeto (ej. "GSE226507", "TREM2")
            decision: Resultado de la decisión (ej. "ACCEPT", "REJECT", "WARN",
                      "edgeR", "PASS", "FAIL")
            reason:   Explicación en lenguaje natural
            evidence: Dict de métricas con estructura:
                      { "metric_name": {"observed": val, "threshold": val, "pass": bool} }
        """
        entry = {
            "subject": subject,
            "decision": decision,
            "reason": reason,
        }
        if evidence:
            entry["evidence"] = evidence
        self._decisions.append(entry)

    def add_warning(self, message: str, subject: Optional[str] = None) -> None:
        """Registra una condición borderline que no impide continuar."""
        prefix = f"[{subject}] " if subject else ""
        self._warnings.append(f"{prefix}{message}")
        logger.warning(f"[{self.agent_id}] AUDIT WARNING: {prefix}{message}")

    def add_error(self, message: str, subject: Optional[str] = None) -> None:
        """Registra un error no-crítico (el agente continuó)."""
        prefix = f"[{subject}] " if subject else ""
        self._errors.append(f"{prefix}{message}")

    def save(self) -> Path:
        """
        Serializa el AuditEntry completo a disco.

        Returns:
            Path al archivo JSON generado.
        """
        elapsed = round(time.time() - self._start_ts, 2)
        entry = {
            "agent_id": self.agent_id,
            "agent_name": self.agent_name,
            "run_id": self.run_id,
            "data_type": self.data_type,
            "started_at": self._started_at,
            "completed_at": _utcnow(),
            "duration_seconds": elapsed,
            "config_snapshot": self._config_snapshot,
            "config_hash": self._config_hash,
            "decisions": self._decisions,
            "warnings": self._warnings,
            "errors": self._errors,
        }

        out_path = self._audit_dir / f"{self.agent_id}.json"
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(entry, f, indent=2, default=str, ensure_ascii=False)

        logger.info(
            f"[{self.agent_id}] Audit guardado → {out_path} "
            f"| decisions={len(self._decisions)} "
            f"| warnings={len(self._warnings)}"
        )
        return out_path

    # ------------------------------------------------------------------
    # Helpers para construir evidencia estándar
    # ------------------------------------------------------------------

    @staticmethod
    def metric(observed: Any, threshold: Any, pass_: bool) -> dict:
        """
        Construye un item de evidencia estándar.

        Uso:
            evidence={
                "samples_per_group": AuditWriter.metric(8, 3, True),
                "balance_ratio": AuditWriter.metric(1.0, 5.0, True),
            }
        """
        return {
            "observed": observed,
            "threshold": threshold,
            "pass": pass_,
        }

    @staticmethod
    def metric_range(observed: Any, min_val: Any, max_val: Any, pass_: bool) -> dict:
        """
        Evidencia para métricas con rango aceptable (min, max).

        Uso:
            "n_sig_genes": AuditWriter.metric_range(1200, 10, 8000, True)
        """
        return {
            "observed": observed,
            "threshold_min": min_val,
            "threshold_max": max_val,
            "pass": pass_,
        }


# ---------------------------------------------------------------------------
# RunManifest — escritor del resumen ejecutivo del run
# ---------------------------------------------------------------------------

class RunManifest:
    """
    Genera y actualiza el run_manifest.json en la raíz del directorio del run.

    El manifest se va actualizando conforme avanzan los agentes:
        manifest = RunManifest(run_id, run_dir, analysis_context, config)
        manifest.update_agent("03_preprocess", {"datasets_passed_qc": 3})
        manifest.finalize(results_summary)
    """

    def __init__(
        self,
        run_id: str,
        run_dir: Path,
        analysis_context: dict,
        pipeline_config: dict,
        pipeline_version: str = "1.1.0",
    ):
        self.run_id = run_id
        self.run_dir = Path(run_dir)
        self.run_dir.mkdir(parents=True, exist_ok=True)

        self._manifest = {
            "run_id": run_id,
            "pipeline_version": pipeline_version,
            "analysis_context": analysis_context,
            "started_at": _utcnow(),
            "completed_at": None,
            "config_hash": _config_hash(pipeline_config),
            "agents": {},
            "datasets": {},
            "results_summary": {},
            "warnings_total": 0,
            "audit_files": [],
        }
        self._save()

    def update_agent(self, agent_id: str, data: dict) -> None:
        """Actualiza el estado de un agente en el manifest."""
        self._manifest["agents"][agent_id] = {
            "completed_at": _utcnow(),
            **data,
        }
        self._save()

    def update_datasets(self, datasets: dict) -> None:
        """
        Actualiza conteos de datasets.

        Args:
            datasets: {
                "discovered": 12,
                "downloaded": 5,
                "passed_qc": 3,
                "used_in_meta": 3,
            }
        """
        self._manifest["datasets"].update(datasets)
        self._save()

    def finalize(self, results_summary: dict, warnings_total: int = 0) -> Path:
        """
        Cierra el manifest con resultados finales.

        Args:
            results_summary: {
                "consensus_genes": 4031,
                "pathways": 143,
                "top_targets": {"therapeutic": [...], "biomarker_diagnostic": [...]}
            }
            warnings_total: suma de warnings de todos los agentes

        Returns:
            Path al manifest guardado.
        """
        self._manifest["completed_at"] = _utcnow()
        self._manifest["results_summary"] = results_summary
        self._manifest["warnings_total"] = warnings_total

        # Listar archivos de auditoría presentes
        audit_dir = self.run_dir / "audit"
        if audit_dir.exists():
            self._manifest["audit_files"] = sorted(
                [f.name for f in audit_dir.glob("*.json")]
            )

        return self._save()

    def _save(self) -> Path:
        out_path = self.run_dir / "run_manifest.json"
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(self._manifest, f, indent=2, default=str, ensure_ascii=False)
        return out_path