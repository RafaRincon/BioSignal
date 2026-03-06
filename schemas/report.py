"""
BioSignal Discovery Engine
Schemas: Final Report Validation
===================================
Pydantic schemas para validación y serialización del reporte final
generado por el Agent 8.
"""

from datetime import datetime
from typing import Literal, Optional
from pydantic import BaseModel, Field, field_validator, model_validator


# ------------------------------------------------------------------
# Sub-schemas para secciones del reporte
# ------------------------------------------------------------------

class TherapeuticTarget(BaseModel):
    """Target terapéutico identificado por el Agent 7."""

    gene: str = Field(..., description="Símbolo HGNC del gen target")
    confidence: Literal["High", "Medium", "Low"] = Field(
        ..., description="Nivel de confianza basado en evidencia"
    )
    rationale: str = Field(default="", description="Justificación biológica")
    known_drugs: list[str] = Field(default_factory=list, description="Fármacos conocidos que modulan este target")
    druggability_score: Literal["high", "medium", "low", "unknown"] = Field(
        default="unknown", description="Score de drogabilidad"
    )
    chembl_id: Optional[str] = Field(default=None, description="ID de ChEMBL del target")

    @field_validator("gene")
    @classmethod
    def validate_gene(cls, v: str) -> str:
        return v.strip().upper()


class KeyMechanism(BaseModel):
    """Mecanismo biológico clave identificado."""

    mechanism: str = Field(..., description="Nombre del mecanismo")
    evidence_strength: Literal["High", "Medium", "Low"] = Field(...)
    genes_involved: list[str] = Field(default_factory=list)
    pathways_involved: list[str] = Field(default_factory=list)
    description: str = Field(default="", description="Descripción mecanicista")


class NovelHypothesis(BaseModel):
    """Hipótesis novedosa generada por el LLM."""

    hypothesis: str = Field(..., description="Enunciado de la hipótesis")
    genes_involved: list[str] = Field(default_factory=list)
    testable_prediction: str = Field(default="", description="Predicción experimental validable")
    confidence: Literal["High", "Medium", "Low", "Speculative"] = Field(default="Speculative")


class BiomarkerCandidate(BaseModel):
    """Gen candidato a biomarcador."""

    gene: str = Field(..., description="Símbolo HGNC")
    biomarker_type: Literal["diagnostic", "prognostic", "predictive"] = Field(...)
    rationale: str = Field(default="")

    @field_validator("gene")
    @classmethod
    def validate_gene(cls, v: str) -> str:
        return v.strip().upper()


class ConsensusGene(BaseModel):
    """Gen consenso del meta-análisis (Agent 5)."""

    gene: str = Field(..., description="Símbolo HGNC")
    log2FC: float = Field(..., description="Log2FC meta-analítico")
    meta_pvalue: float = Field(..., ge=0.0, le=1.0)
    meta_padj: float = Field(..., ge=0.0, le=1.0)
    direction: Literal["UP", "DOWN"] = Field(...)
    n_datasets: int = Field(..., ge=1, description="Número de datasets donde es significativo")
    consistency_score: Optional[float] = Field(
        default=None, ge=0.0, le=1.0,
        description="Fracción de datasets donde mantiene la misma dirección"
    )

    @field_validator("gene")
    @classmethod
    def validate_gene(cls, v: str) -> str:
        return v.strip().upper()


class PipelineMetrics(BaseModel):
    """Métricas de ejecución del pipeline completo."""

    datasets_analyzed: int = Field(default=0, ge=0)
    datasets_passed_qc: int = Field(default=0, ge=0)
    datasets_failed: int = Field(default=0, ge=0)
    total_samples: int = Field(default=0, ge=0)
    total_genes_tested: int = Field(default=0, ge=0)
    execution_time_seconds: Optional[float] = Field(default=None)
    pipeline_version: str = Field(default="1.0.0")


# ------------------------------------------------------------------
# Schema principal del reporte
# ------------------------------------------------------------------

class BioSignalReport(BaseModel):
    """
    Schema completo del reporte final de BioSignal Discovery Engine.
    Consolida outputs de los 8 agentes del pipeline.
    """

    # Metadata
    disease: str = Field(..., description="Enfermedad o condición analizada")
    analysis_date: str = Field(default_factory=lambda: datetime.utcnow().strftime("%Y-%m-%d"))
    generated_at: datetime = Field(default_factory=datetime.utcnow)
    pipeline_version: str = Field(default="1.0.0")
    report_id: Optional[str] = Field(default=None, description="Identificador único del reporte")

    # Sección 1: Resumen ejecutivo
    executive_summary: str = Field(default="", description="Resumen ejecutivo generado por LLM")

    # Sección 2: Genes consenso (Agent 5)
    consensus_genes: list[ConsensusGene] = Field(
        default_factory=list,
        description="Genes diferenciados reproducibles entre datasets"
    )
    total_consensus_genes: int = Field(default=0)
    consensus_genes_up: list[str] = Field(default_factory=list)
    consensus_genes_down: list[str] = Field(default_factory=list)

    # Sección 3: Pathways (Agent 6)
    top_pathways_kegg: list[dict] = Field(default_factory=list)
    top_pathways_go: list[dict] = Field(default_factory=list)
    top_pathways_reactome: list[dict] = Field(default_factory=list)
    top_pathways_hallmarks: list[dict] = Field(default_factory=list)

    # Sección 4: Insights del LLM (Agent 7)
    key_mechanisms: list[KeyMechanism] = Field(default_factory=list)
    therapeutic_targets: list[TherapeuticTarget] = Field(default_factory=list)
    novel_hypotheses: list[NovelHypothesis] = Field(default_factory=list)
    biomarker_candidates: list[BiomarkerCandidate] = Field(default_factory=list)

    # Sección 5: Calidad y métricas
    pipeline_metrics: PipelineMetrics = Field(default_factory=PipelineMetrics)
    limitations: list[str] = Field(default_factory=list)
    recommended_next_steps: list[str] = Field(default_factory=list)

    # Archivos generados
    output_files: dict[str, str] = Field(
        default_factory=dict,
        description="Mapa formato → ruta del archivo generado"
    )

    @model_validator(mode="after")
    def populate_derived_fields(self) -> "BioSignalReport":
        """Calcula campos derivados automáticamente."""
        self.total_consensus_genes = len(self.consensus_genes)
        self.consensus_genes_up = [g.gene for g in self.consensus_genes if g.direction == "UP"]
        self.consensus_genes_down = [g.gene for g in self.consensus_genes if g.direction == "DOWN"]
        return self

    @field_validator("disease")
    @classmethod
    def validate_disease(cls, v: str) -> str:
        v = v.strip()
        if not v:
            raise ValueError("El campo 'disease' no puede estar vacío")
        return v

    def top_targets_by_confidence(self, confidence: str = "High") -> list[TherapeuticTarget]:
        """Filtra targets por nivel de confianza."""
        return [t for t in self.therapeutic_targets if t.confidence == confidence]

    def to_summary_dict(self) -> dict:
        """Retorna un diccionario resumido para logging y CLI output."""
        return {
            "disease": self.disease,
            "analysis_date": self.analysis_date,
            "n_datasets": self.pipeline_metrics.datasets_analyzed,
            "n_samples": self.pipeline_metrics.total_samples,
            "n_consensus_genes": self.total_consensus_genes,
            "n_genes_up": len(self.consensus_genes_up),
            "n_genes_down": len(self.consensus_genes_down),
            "n_therapeutic_targets": len(self.therapeutic_targets),
            "n_high_confidence_targets": len(self.top_targets_by_confidence("High")),
            "n_novel_hypotheses": len(self.novel_hypotheses),
            "output_files": self.output_files,
        }

    model_config = {"str_strip_whitespace": True}