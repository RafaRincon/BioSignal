"""
BioSignal Discovery Engine
Agent 8: ReportGenerationAgent
================================
Responsabilidad: Consolidar todos los resultados del pipeline y generar
el reporte final en múltiples formatos: PDF, JSON y CSV.

Uso:
    python -m agents.report --disease "pancreatic cancer" --output data/reports/
"""

import json
import time
from datetime import datetime
from pathlib import Path
from typing import Optional

import click
import pandas as pd
from loguru import logger

try:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4, letter
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import inch, cm
    from reportlab.platypus import (
        SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
        HRFlowable, PageBreak, KeepTogether,
    )
    from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
    REPORTLAB_AVAILABLE = True
except ImportError:
    REPORTLAB_AVAILABLE = False
    logger.warning("reportlab no disponible. Instalar con: pip install reportlab")


class ReportGenerationAgent:
    """
    Agente de generación de reportes finales.

    Consolida datos de todos los agentes anteriores y genera:
        - PDF: reporte ejecutivo con tablas de genes, pathways y targets
        - JSON: datos estructurados completos del análisis
        - CSV: genes significativos listos para análisis posterior

    Estructura del reporte PDF:
        1. Portada con metadata del análisis
        2. Resumen ejecutivo (del Agent 7)
        3. Tabla de genes consenso (Agent 5)
        4. Pathways enriquecidos (Agent 6)
        5. Targets terapéuticos (Agent 7)
        6. Hipótesis y próximos pasos (Agent 7)
        7. Métricas del pipeline

    Salida:
        data/reports/{disease}/report_{disease}_{date}.pdf
        data/reports/{disease}/report_{disease}_{date}.json
        data/reports/{disease}/consensus_genes_{disease}.csv
        data/reports/{disease}/report_summary.json
    """

    VERSION = "1.0.0"

    def __init__(self, config: dict = None):
        self.config = config or {}
        self.formats = self.config.get("formats", ["json"])
        self.top_genes_table = self.config.get("top_genes_table", 30)
        self.top_pathways_table = self.config.get("top_pathways_table", 15)
        self.top_targets_table = self.config.get("top_targets_table", 10)

    # ------------------------------------------------------------------
    # Entrada principal
    # ------------------------------------------------------------------

    def run(self, disease_name: str, data_dir: str, output_dir: str) -> dict:
        """
        Genera el reporte final consolidando todos los resultados del pipeline.

        Args:
            disease_name: Nombre de la enfermedad analizada
            data_dir:     Directorio raíz de datos (data/)
            output_dir:   Directorio destino para reportes (data/reports/)

        Returns:
            Rutas de archivos generados
        """
        data_path = Path(data_dir)
        out_path = Path(output_dir) / self._sanitize_name(disease_name)
        out_path.mkdir(parents=True, exist_ok=True)

        start = time.time()
        date_str = datetime.utcnow().strftime("%Y%m%d")
        safe_name = self._sanitize_name(disease_name)

        logger.info(f"[ReportAgent] Generando reporte para: {disease_name}")

        # 1. Recopilar todos los datos del pipeline
        report_data = self._collect_pipeline_data(disease_name, data_path)

        generated_files = {}

        # 2. Generar JSON
        if "json" in self.formats:
            json_path = out_path / f"report_{safe_name}_{date_str}.json"
            self._generate_json(report_data, json_path)
            generated_files["json"] = str(json_path)
            logger.info(f"[ReportAgent] JSON generado: {json_path}")

        # 3. Generar CSV de genes
        if "csv" in self.formats:
            csv_path = out_path / f"consensus_genes_{safe_name}.csv"
            self._generate_csv(report_data, csv_path)
            generated_files["csv"] = str(csv_path)
            logger.info(f"[ReportAgent] CSV generado: {csv_path}")

        # 4. Generar PDF
        if "pdf" in self.formats:
            pdf_path = out_path / f"report_{safe_name}_{date_str}.pdf"
            if REPORTLAB_AVAILABLE:
                self._generate_pdf(report_data, pdf_path, disease_name)
                generated_files["pdf"] = str(pdf_path)
                logger.info(f"[ReportAgent] PDF generado: {pdf_path}")
            else:
                logger.warning("[ReportAgent] reportlab no disponible — PDF omitido")

        elapsed = time.time() - start

        # 5. Guardar resumen del reporte
        summary = {
            "disease": disease_name,
            "generated_at": datetime.utcnow().isoformat() + "Z",
            "pipeline_version": self.VERSION,
            "files_generated": generated_files,
            "elapsed_seconds": round(elapsed, 2),
            "n_genes": len(report_data.get("consensus_genes", [])),
            "n_datasets": report_data.get("pipeline_metrics", {}).get("datasets_analyzed", 0),
        }

        with open(out_path / "report_summary.json", "w") as f:
            json.dump(summary, f, indent=2)

        logger.info(f"[ReportAgent] Completado en {elapsed:.1f}s | Archivos: {list(generated_files.keys())}")
        return summary

    # ------------------------------------------------------------------
    # Recopilación de datos
    # ------------------------------------------------------------------

    def _collect_pipeline_data(self, disease_name: str, data_path: Path) -> dict:
        """Consolida datos de todos los agentes del pipeline."""
        safe = self._sanitize_name(disease_name)

        # Agent 5: Meta-análisis
        consensus_genes = self._load_json_safe(data_path / "meta" / safe / "consensus_genes.json")
        meta_results = self._load_csv_safe(data_path / "meta" / safe / "meta_analysis_results.csv")

        # Agent 6: Pathways
        top_pathways = self._load_json_safe(data_path / "pathways" / safe / "top_pathways_summary.json")

        # Agent 7: Insights
        insights = self._load_json_safe(data_path / "insights" / safe / "biological_insights.json")
        therapeutic_targets = self._load_json_safe(data_path / "insights" / safe / "therapeutic_targets.json")

        # Métricas del pipeline
        pipeline_metrics = self._compute_pipeline_metrics(data_path, safe)

        return {
            "disease": disease_name,
            "analysis_date": datetime.utcnow().strftime("%Y-%m-%d"),
            "pipeline_version": self.VERSION,
            "consensus_genes": consensus_genes.get("genes", []) if consensus_genes else [],
            "meta_results_df": meta_results,
            "top_pathways": top_pathways if top_pathways else {},
            "executive_summary": insights.get("executive_summary", "") if insights else "",
            "key_mechanisms": insights.get("key_mechanisms", []) if insights else [],
            "therapeutic_targets": therapeutic_targets.get("targets", []) if therapeutic_targets else [],
            "novel_hypotheses": insights.get("novel_hypotheses", []) if insights else [],
            "biomarker_candidates": insights.get("biomarker_candidates", []) if insights else [],
            "recommended_next_steps": insights.get("recommended_next_steps", []) if insights else [],
            "limitations": insights.get("limitations", []) if insights else [],
            "pipeline_metrics": pipeline_metrics,
        }

    def _compute_pipeline_metrics(self, data_path: Path, safe_name: str) -> dict:
        """Calcula métricas del pipeline a partir de los archivos de salida."""
        metrics = {}

        # Contar datasets analizados en DEA
        dea_path = data_path / "dea"
        if dea_path.exists():
            dea_dirs = [d for d in dea_path.iterdir() if d.is_dir()]
            metrics["datasets_analyzed"] = len(dea_dirs)
            # Sumar muestras totales
            total_sig = 0
            for d in dea_dirs:
                summary_f = d / "dea_summary.json"
                if summary_f.exists():
                    s = self._load_json_safe(summary_f) or {}
                    total_sig += s.get("n_significant", 0)
            metrics["total_significant_genes_sum"] = total_sig

        # Total muestras desde discovery
        disc_path = data_path / "discovery"
        if disc_path.exists():
            for f in disc_path.glob("*.json"):
                disc_data = self._load_json_safe(f) or {}
                datasets = disc_data.get("datasets", [])
                metrics["total_samples"] = sum(d.get("sample_count", 0) for d in datasets)
                break

        return metrics

    # ------------------------------------------------------------------
    # Generación JSON
    # ------------------------------------------------------------------

    def _generate_json(self, report_data: dict, output_path: Path):
        """Genera el JSON estructurado completo del reporte."""
        export = {
            "disease": report_data["disease"],
            "analysis_date": report_data["analysis_date"],
            "pipeline_version": report_data["pipeline_version"],
            "datasets_analyzed": report_data["pipeline_metrics"].get("datasets_analyzed", 0),
            "total_samples": report_data["pipeline_metrics"].get("total_samples", 0),
            "executive_summary": report_data["executive_summary"],
            "consensus_genes": report_data["consensus_genes"][:self.top_genes_table],
            "top_pathways": report_data["top_pathways"],
            "therapeutic_targets": report_data["therapeutic_targets"][:self.top_targets_table],
            "key_mechanisms": report_data["key_mechanisms"],
            "novel_hypotheses": report_data["novel_hypotheses"],
            "biomarker_candidates": report_data["biomarker_candidates"],
            "limitations": report_data["limitations"],
            "recommended_next_steps": report_data["recommended_next_steps"],
            "pipeline_metrics": report_data["pipeline_metrics"],
        }

        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(export, f, indent=2, ensure_ascii=False, default=str)

    # ------------------------------------------------------------------
    # Generación CSV
    # ------------------------------------------------------------------

    def _generate_csv(self, report_data: dict, output_path: Path):
        """Genera CSV de genes consenso."""
        genes = report_data.get("consensus_genes", [])
        if not genes:
            logger.warning("[ReportAgent] Sin genes consenso para exportar a CSV")
            return

        df = pd.DataFrame(genes)

        # Asegurar columnas estándar
        for col in ["gene", "log2FC", "meta_pvalue", "meta_padj", "direction", "n_datasets"]:
            if col not in df.columns:
                df[col] = ""

        df = df.sort_values("meta_padj") if "meta_padj" in df.columns else df
        df.to_csv(output_path, index=False, float_format="%.6f")

    # ------------------------------------------------------------------
    # Generación PDF
    # ------------------------------------------------------------------

    def _generate_pdf(self, report_data: dict, output_path: Path, disease_name: str):
        """Genera el reporte PDF completo usando reportlab."""
        doc = SimpleDocTemplate(
            str(output_path),
            pagesize=A4,
            rightMargin=2 * cm,
            leftMargin=2 * cm,
            topMargin=2.5 * cm,
            bottomMargin=2 * cm,
        )

        styles = getSampleStyleSheet()
        story = []

        # --- Estilos personalizados ---
        title_style = ParagraphStyle("BioTitle", parent=styles["Title"],
                                      fontSize=22, spaceAfter=6, textColor=colors.HexColor("#1a3a5c"))
        subtitle_style = ParagraphStyle("BioSubtitle", parent=styles["Normal"],
                                         fontSize=12, spaceAfter=20, textColor=colors.HexColor("#4a7fb5"))
        h1_style = ParagraphStyle("BioH1", parent=styles["Heading1"],
                                   fontSize=14, spaceBefore=14, spaceAfter=6,
                                   textColor=colors.HexColor("#1a3a5c"))
        h2_style = ParagraphStyle("BioH2", parent=styles["Heading2"],
                                   fontSize=11, spaceBefore=10, spaceAfter=4,
                                   textColor=colors.HexColor("#2c5f8a"))
        body_style = ParagraphStyle("BioBody", parent=styles["Normal"],
                                     fontSize=9, spaceAfter=6, leading=14, alignment=TA_JUSTIFY)
        caption_style = ParagraphStyle("BioCaption", parent=styles["Normal"],
                                        fontSize=8, textColor=colors.grey, spaceAfter=4)

        # ---- PORTADA ----
        story.append(Spacer(1, 1.5 * cm))
        story.append(Paragraph("BioSignal Discovery Engine", title_style))
        story.append(Paragraph("Biological Signal Discovery Report", subtitle_style))
        story.append(HRFlowable(width="100%", thickness=2, color=colors.HexColor("#1a3a5c")))
        story.append(Spacer(1, 0.5 * cm))

        meta_data = [
            ["Disease / Condition:", disease_name],
            ["Analysis Date:", report_data["analysis_date"]],
            ["Pipeline Version:", report_data["pipeline_version"]],
            ["Datasets Analyzed:", str(report_data["pipeline_metrics"].get("datasets_analyzed", "N/A"))],
            ["Total Samples:", str(report_data["pipeline_metrics"].get("total_samples", "N/A"))],
            ["Genes Identified:", str(len(report_data["consensus_genes"]))],
        ]
        meta_table = Table(meta_data, colWidths=[5 * cm, 12 * cm])
        meta_table.setStyle(TableStyle([
            ("FONTSIZE", (0, 0), (-1, -1), 9),
            ("FONTNAME", (0, 0), (0, -1), "Helvetica-Bold"),
            ("TEXTCOLOR", (0, 0), (0, -1), colors.HexColor("#1a3a5c")),
            ("ROWBACKGROUNDS", (0, 0), (-1, -1), [colors.HexColor("#f0f5fa"), colors.white]),
            ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#ccddee")),
            ("LEFTPADDING", (0, 0), (-1, -1), 8),
            ("RIGHTPADDING", (0, 0), (-1, -1), 8),
            ("TOPPADDING", (0, 0), (-1, -1), 5),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ]))
        story.append(meta_table)
        story.append(PageBreak())

        # ---- SECCIÓN 1: RESUMEN EJECUTIVO ----
        story.append(Paragraph("1. Executive Summary", h1_style))
        story.append(HRFlowable(width="100%", thickness=0.5, color=colors.HexColor("#4a7fb5")))
        story.append(Spacer(1, 0.3 * cm))

        exec_summary = report_data.get("executive_summary") or "Executive summary not available."
        story.append(Paragraph(exec_summary, body_style))
        story.append(Spacer(1, 0.5 * cm))

        # ---- SECCIÓN 2: GENES CONSENSO ----
        story.append(Paragraph("2. Consensus Differentially Expressed Genes", h1_style))
        story.append(HRFlowable(width="100%", thickness=0.5, color=colors.HexColor("#4a7fb5")))
        story.append(Spacer(1, 0.3 * cm))

        genes = report_data.get("consensus_genes", [])[:self.top_genes_table]
        if genes:
            gene_table_data = [["Gene", "log2FC", "Meta padj", "Direction", "# Datasets"]]
            for g in genes:
                gene_table_data.append([
                    g.get("gene", ""),
                    f"{g.get('log2FC', 0):.3f}",
                    f"{g.get('meta_padj', g.get('padj', '?')):.2e}",
                    g.get("direction", "?"),
                    str(g.get("n_datasets", "?")),
                ])

            gene_table = Table(gene_table_data, colWidths=[3.5 * cm, 2.5 * cm, 3 * cm, 2.5 * cm, 3 * cm])
            gene_table.setStyle(self._base_table_style())
            story.append(gene_table)
            story.append(Paragraph(f"Table 1. Top {len(genes)} consensus genes reproducibly dysregulated across datasets.", caption_style))
        else:
            story.append(Paragraph("No consensus genes data available.", body_style))

        story.append(Spacer(1, 0.5 * cm))

        # ---- SECCIÓN 3: PATHWAYS ----
        story.append(Paragraph("3. Enriched Biological Pathways", h1_style))
        story.append(HRFlowable(width="100%", thickness=0.5, color=colors.HexColor("#4a7fb5")))
        story.append(Spacer(1, 0.3 * cm))

        pathways = report_data.get("top_pathways", {})
        top_by_db = pathways.get("top_pathways_by_db", {})

        for db_key, db_paths in top_by_db.items():
            if not db_paths:
                continue
            story.append(Paragraph(f"3.{list(top_by_db.keys()).index(db_key)+1} {db_key.upper()}", h2_style))
            pw_data = [["Pathway", "padj", "Combined Score"]]
            for p in db_paths[:8]:
                pw_data.append([
                    p.get("pathway", "")[:60],
                    f"{p.get('padj', 1):.2e}",
                    f"{p.get('combined_score', 0):.2f}",
                ])
            pw_table = Table(pw_data, colWidths=[10 * cm, 2.5 * cm, 3 * cm])
            pw_table.setStyle(self._base_table_style())
            story.append(pw_table)
            story.append(Spacer(1, 0.3 * cm))

        story.append(PageBreak())

        # ---- SECCIÓN 4: TARGETS TERAPÉUTICOS ----
        story.append(Paragraph("4. Therapeutic Target Candidates", h1_style))
        story.append(HRFlowable(width="100%", thickness=0.5, color=colors.HexColor("#4a7fb5")))
        story.append(Spacer(1, 0.3 * cm))

        targets = report_data.get("therapeutic_targets", [])[:self.top_targets_table]
        if targets:
            tgt_data = [["Gene", "Confidence", "Known Drugs", "Rationale"]]
            for t in targets:
                drugs = ", ".join(t.get("known_drugs", [])) or "None reported"
                rationale = (t.get("rationale", "") or "")[:80]
                tgt_data.append([
                    t.get("gene", ""),
                    t.get("confidence", "?"),
                    drugs[:30],
                    rationale,
                ])
            tgt_table = Table(tgt_data, colWidths=[2.5 * cm, 2.5 * cm, 4 * cm, 7 * cm])
            tgt_table.setStyle(self._base_table_style())
            story.append(tgt_table)
        else:
            story.append(Paragraph("Therapeutic targets not available.", body_style))

        story.append(Spacer(1, 0.5 * cm))

        # ---- SECCIÓN 5: HIPÓTESIS ----
        story.append(Paragraph("5. Novel Biological Hypotheses", h1_style))
        story.append(HRFlowable(width="100%", thickness=0.5, color=colors.HexColor("#4a7fb5")))
        story.append(Spacer(1, 0.3 * cm))

        for i, hyp in enumerate(report_data.get("novel_hypotheses", []), 1):
            story.append(Paragraph(f"<b>Hypothesis {i}:</b> {hyp.get('hypothesis', '')}", body_style))
            if hyp.get("testable_prediction"):
                story.append(Paragraph(f"<i>Validation:</i> {hyp.get('testable_prediction', '')}", body_style))
            story.append(Spacer(1, 0.2 * cm))

        story.append(Spacer(1, 0.3 * cm))

        # ---- SECCIÓN 6: PRÓXIMOS PASOS ----
        story.append(Paragraph("6. Recommended Next Steps", h1_style))
        story.append(HRFlowable(width="100%", thickness=0.5, color=colors.HexColor("#4a7fb5")))
        story.append(Spacer(1, 0.3 * cm))

        for step in report_data.get("recommended_next_steps", []):
            story.append(Paragraph(f"• {step}", body_style))

        story.append(Spacer(1, 0.3 * cm))
        story.append(Paragraph("<b>Limitations:</b>", body_style))
        for lim in report_data.get("limitations", []):
            story.append(Paragraph(f"• {lim}", body_style))

        # Footer via onFirstPage / onLaterPages
        def add_footer(canvas, doc):
            canvas.saveState()
            canvas.setFont("Helvetica", 7)
            canvas.setFillColor(colors.grey)
            canvas.drawString(2 * cm, 1 * cm,
                              f"BioSignal Discovery Engine v{self.VERSION} | {disease_name} | {report_data['analysis_date']} | Confidential")
            canvas.drawRightString(A4[0] - 2 * cm, 1 * cm, f"Page {doc.page}")
            canvas.restoreState()

        doc.build(story, onFirstPage=add_footer, onLaterPages=add_footer)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _base_table_style(self) -> TableStyle:
        """Estilo base para tablas del reporte."""
        return TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1a3a5c")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE", (0, 0), (-1, -1), 8),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.HexColor("#f0f5fa"), colors.white]),
            ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#ccddee")),
            ("ALIGN", (1, 0), (-1, -1), "CENTER"),
            ("ALIGN", (0, 0), (0, -1), "LEFT"),
            ("LEFTPADDING", (0, 0), (-1, -1), 6),
            ("RIGHTPADDING", (0, 0), (-1, -1), 6),
            ("TOPPADDING", (0, 0), (-1, -1), 4),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ])

    def _sanitize_name(self, name: str) -> str:
        """Convierte nombre de enfermedad a nombre de archivo válido."""
        return name.lower().replace(" ", "_").replace("/", "_").replace("'", "")

    def _load_json_safe(self, path: Path) -> Optional[dict]:
        if path and path.exists():
            try:
                with open(path) as f:
                    return json.load(f)
            except Exception as e:
                logger.warning(f"Error cargando {path}: {e}")
        return None

    def _load_csv_safe(self, path: Path) -> Optional[pd.DataFrame]:
        if path and path.exists():
            try:
                return pd.read_csv(path)
            except Exception as e:
                logger.warning(f"Error cargando {path}: {e}")
        return None


# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------

@click.command()
@click.option("--disease", required=True, help="Nombre de la enfermedad analizada")
@click.option("--data-dir", default="data/", help="Directorio raíz de datos")
@click.option("--output", "output_dir", default="data/reports/", help="Directorio de salida")
@click.option("--formats", default="json", help="Formatos de salida separados por coma (json, pdf, csv)")
def main(disease, data_dir, output_dir, formats):
    """Agent 8: Genera el reporte final del análisis en PDF, JSON y CSV."""
    config = {"formats": [f.strip() for f in formats.split(",")]}
    agent = ReportGenerationAgent(config=config)
    summary = agent.run(disease, data_dir, output_dir)
    click.echo(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()