"""
BioSignal Discovery Engine
Schemas: DEA Result Validation
================================
Pydantic schemas para validación y serialización de resultados
del análisis de expresión diferencial (Agent 4).
"""

from datetime import datetime
from typing import Literal, Optional
from pydantic import BaseModel, Field, field_validator, model_validator


class DEAGeneResult(BaseModel):
    """Schema para el resultado DEA de un gen individual."""

    gene: str = Field(..., description="Símbolo HGNC del gen")
    log2FC: float = Field(..., description="Log2 Fold Change (case/control)")
    pvalue: float = Field(..., ge=0.0, le=1.0, description="P-valor sin ajustar")
    padj: float = Field(..., ge=0.0, le=1.0, description="P-valor ajustado (FDR)")
    direction: Literal["UP", "DOWN", "NS"] = Field(
        default="NS", description="Dirección de cambio: UP | DOWN | NS (no significativo)"
    )
    base_mean: Optional[float] = Field(default=None, description="Expresión media (DESeq2)")
    stat: Optional[float] = Field(default=None, description="Estadístico del test")

    @field_validator("gene")
    @classmethod
    def validate_gene_symbol(cls, v: str) -> str:
        v = v.strip().upper()
        if not v:
            raise ValueError("Símbolo génico vacío")
        return v

    @field_validator("log2FC")
    @classmethod
    def validate_lfc(cls, v: float) -> float:
        if abs(v) > 20:
            # log2FC extremo — posible error numérico, recortar
            return max(-20.0, min(20.0, v))
        return v

    @model_validator(mode="after")
    def set_direction(self) -> "DEAGeneResult":
        """Establece dirección automáticamente si no viene definida."""
        if self.direction == "NS":
            if self.log2FC > 0:
                self.direction = "UP"
            elif self.log2FC < 0:
                self.direction = "DOWN"
        return self

    model_config = {"str_strip_whitespace": True}


class DEASummary(BaseModel):
    """Schema para el resumen del análisis DEA de un dataset."""

    gse_id: str = Field(..., description="ID del dataset GEO")
    method: Literal["deseq2", "limma", "ttest_python"] = Field(
        ..., description="Método estadístico utilizado"
    )
    timestamp: datetime = Field(default_factory=datetime.utcnow)

    n_total_genes: int = Field(..., ge=0, description="Total de genes analizados")
    n_significant: int = Field(..., ge=0, description="Genes significativos (padj < threshold AND |LFC| > threshold)")
    n_up: int = Field(..., ge=0, description="Genes sobreexpresados")
    n_down: int = Field(..., ge=0, description="Genes subexpresados")

    padj_threshold: float = Field(default=0.05, ge=0.0, le=1.0)
    lfc_threshold: float = Field(default=1.0, ge=0.0)

    n_samples_case: Optional[int] = Field(default=None, ge=0)
    n_samples_control: Optional[int] = Field(default=None, ge=0)
    data_type: Optional[str] = Field(default=None, description="'RNA-seq' | 'microarray'")

    @field_validator("n_significant")
    @classmethod
    def check_significant_vs_total(cls, v: int, info) -> int:
        total = info.data.get("n_total_genes", 0)
        if total > 0 and v > total:
            raise ValueError(f"n_significant ({v}) > n_total_genes ({total})")
        return v

    @model_validator(mode="after")
    def check_up_down_sum(self) -> "DEASummary":
        if self.n_up + self.n_down > self.n_significant:
            raise ValueError(
                f"n_up ({self.n_up}) + n_down ({self.n_down}) > n_significant ({self.n_significant})"
            )
        return self


class SignificantGenesOutput(BaseModel):
    """Schema para el archivo significant_genes.json generado por el Agent 4."""

    gse_id: str
    genes: list[DEAGeneResult] = Field(default_factory=list)
    n_genes: int = Field(default=0)

    @model_validator(mode="after")
    def set_n_genes(self) -> "SignificantGenesOutput":
        self.n_genes = len(self.genes)
        return self

    def genes_up(self) -> list[DEAGeneResult]:
        """Retorna genes sobreexpresados."""
        return [g for g in self.genes if g.direction == "UP"]

    def genes_down(self) -> list[DEAGeneResult]:
        """Retorna genes subexpresados."""
        return [g for g in self.genes if g.direction == "DOWN"]

    def to_dataframe(self):
        """Convierte a pandas DataFrame."""
        import pandas as pd
        return pd.DataFrame([g.model_dump() for g in self.genes])