# Changelog — BioSignal Discovery Engine

Todos los cambios notables del proyecto serán documentados en este archivo.

Formato basado en [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [1.0.0] — 2026-03-01 — MVP Local

### Agregado
- **Agent 1 (DatasetDiscoveryAgent):** Búsqueda automática en NCBI/GEO via E-utilities
- **Agent 2 (DatasetDownloadAgent):** Descarga paralela con reintentos y clasificación de muestras
- **Agent 3 (PreprocessingAgent):** Normalización TMM/quantile, QC, detección de outliers, ComBat
- **Agent 4 (DifferentialExpressionAgent):** DESeq2/limma/edgeR via pydeseq2 y rpy2
- **Agent 5 (MetaAnalysisAgent):** Meta-análisis Fisher + Stouffer + vote counting — feature diferenciadora
- **Agent 6 (PathwayEnrichmentAgent):** Enrichment en KEGG, GO, Reactome, MSigDB via gseapy
- **Agent 7 (InsightGenerationAgent):** LLM (Claude/GPT-4/llama3) para hipótesis terapéuticas
- **Agent 8 (ReportGenerationAgent):** PDF + JSON + CSV + Markdown via reportlab
- **CLI completa:** `biosignal run`, `biosignal status`, `biosignal doctor`, `biosignal list-diseases`
- **Schemas Pydantic v2** para validación de todos los inputs/outputs
- **Suite de tests** con pytest para los módulos core
- **Documentación completa:** README, INSTALLATION, AGENTS, API, TROUBLESHOOTING, ROADMAP
- **Aliases de enfermedades:** 40+ enfermedades con términos canónicos GEO predefinidos
- **Configuración YAML** con todos los parámetros ajustables

### Soporte
- Python 3.11+, R 4.3+
- Ubuntu 22.04 / macOS 13+
- LLM backends: Ollama (local), Anthropic Claude, OpenAI GPT-4

---

## Próximas Versiones

### [1.1.0] — Mes 3 — Validación Científica (Planificado)
- Benchmarking contra papers publicados conocidos
- Reproducción de genes reportados para cáncer pancreático, Alzheimer, AML
- Métricas de precisión/recall vs literatura
- Panel de validación con 10+ enfermedades

### [2.0.0] — Mes 4-5 — API REST (Planificado)
- FastAPI backend con autenticación JWT
- Rate limiting y request queuing
- Soporte para 10 requests concurrentes
- Endpoints: `/analyze`, `/status/{job_id}`, `/results/{job_id}`
- Docker Compose para deployment

### [3.0.0] — Mes 6-8 — UI Web (Planificado)
- Interface React con visualizaciones interactivas
- Volcano plots, heatmaps y pathway maps interactivos
- Dashboard de jobs y resultados históricos
- Onboarding para usuarios sin bioinformática

### [4.0.0] — Mes 9-12 — Multi-omics (Planificado)
- Integración proteómica (UniProt, PhosphoSitePlus)
- Integración GWAS (GWAS Catalog)
- Datos clínicos (GEO clinical data)
- Knowledge graph multi-omics con Neo4j

### [5.0.0] — Mes 13-15 — SaaS Beta (Planificado)
- Producto comercial completo
- Pricing tier: Academic / Biotech / Enterprise
- SOC2 compliance para datos propietarios
- 5 primeros clientes pagos