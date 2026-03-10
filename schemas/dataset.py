"""
BioSignal Discovery Engine
Schemas: Dataset Validation
============================
Pydantic schemas para validación y serialización de datasets.
"""

from datetime import datetime
from typing import Optional
from pydantic import BaseModel, Field, field_validator


class DatasetMetadata(BaseModel):
    """Schema para metadata de un dataset GEO."""

    gse_id: str = Field(..., description="Identificador GEO, ej. GSE12345")
    title: str = Field(default="", description="Título del estudio")
    organism: str = Field(default="Homo sapiens", description="Organismo")
    sample_count: int = Field(..., ge=0, description="Número total de muestras")
    platform: str = Field(default="", description="Plataforma GPL")
    year: int = Field(default=2000, ge=1990, le=2030, description="Año de publicación")
    pmid: Optional[str] = Field(default=None, description="PubMed ID asociado")
    data_type: str = Field(default="unknown", description="'RNA-seq' | 'microarray'")
    summary: str = Field(default="", description="Resumen del estudio")
    rank_score: float = Field(default=0.0, description="Score de ranking (0-1)")
    downloadability_score: int = Field(default=1, description="Score de descargabilidad (0-3)")

    @field_validator("gse_id")
    @classmethod
    def validate_gse_id(cls, v: str) -> str:
        if not v.startswith("GSE"):
            raise ValueError(f"GSE ID inválido: {v}. Debe comenzar con 'GSE'")
        return v.strip()

    @field_validator("data_type")
    @classmethod
    def validate_data_type(cls, v: str) -> str:
        valid = {"RNA-seq", "microarray", "unknown"}
        if v not in valid:
            return "unknown"
        return v

    model_config = {"str_strip_whitespace": True}


class DiscoveryOutput(BaseModel):
    """Schema para el output completo del DatasetDiscoveryAgent."""

    disease_name: str
    timestamp: datetime = Field(default_factory=datetime.utcnow)
    total_found: int
    datasets: list[DatasetMetadata]


class DownloadStatus(BaseModel):
    """Schema para el estado de descarga de un dataset."""

    gse_id: str
    success: bool
    output_dir: str
    error: Optional[str] = None
    matrix_path: Optional[str] = None
    metadata_path: Optional[str] = None
    platform_path: Optional[str] = None
    n_samples: Optional[int] = None
    n_case: Optional[int] = None
    n_control: Optional[int] = None
    n_unclassified: Optional[int] = None
    platform: Optional[str] = None
