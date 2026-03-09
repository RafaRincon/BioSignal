# Changelog — BioSignal Discovery Engine

Todos los cambios notables del proyecto serán documentados en este archivo.

Formato basado en [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [1.0.1] — 2026-03-08 — MVP Debugging Session

### Corregido
- ✅ **Agent 2 (download.py):** Lógica de clasificación de muestras — keywords más específicos (`cckp`, `gckp`, `sr4`, `lnd`) con prioridad por longitud para evitar falsos matches (ej. `cckp` sobreescribe `ckp`)
- ✅ **Agent 3 (preprocess.py):** `_load_sample_metadata` ahora acepta `matrix_cols` para alinear índice del metadata con columnas de la matriz — resuelve mismatch GSM IDs vs títulos de muestra
- ✅ **Agent 3 (preprocess.py):** Ahora guarda `matrix_counts.csv` (counts crudos pre-normalización) además de `matrix_normalized.csv` — necesario para edgeR
- ✅ **Agent 3 (preprocess.py):** Incluye muestras `unclassified` del `metadata.json` al construir `sample_metadata.csv`
- ✅ **Agent 4 (dea.py):** Alineación de matrix columns con sample_meta index antes de llamar a R (resuelve `Two subscripts required` por muestras extra en metadata)
- ✅ **Agent 4 (dea.py):** `_run_edger` ahora usa counts crudos (`matrix_counts.csv`) y pasa `matrix_norm` a fallbacks limma/ttest
- ✅ **Agent 4 (dea.py):** Firma de `_run_edger` actualizada a `(matrix, matrix_norm, sample_meta, gse_id)`
- ✅ **Agent 4 (dea.py):** Subsetting de DGEList via `ro.r("dge[keep,,keep.lib.sizes=FALSE]")` — reemplaza `base.subset()` y `dge.rx()` que fallan en Windows
- ✅ **Agent 5 (meta_analysis.py):** Normalización de gene IDs cross-dataset via `mygene` — ENSG→symbol, ENSMUSG→human ortholog, symbol→symbol — resuelve 0 overlap entre datasets
- ✅ **Agent 5 (meta_analysis.py):** Usa `mygene.MyGeneInfo()` en lugar de requests POST (endpoint `/v3/querymany` no accesible vía requests en este entorno)
- ✅ **Agent 6 (pathways.py):** Clave `consensus_genes` en lugar de `genes` al leer `consensus_genes.json`
- ✅ **Agent 6 (pathways.py):** `organism` configurable via `settings.yaml` en lugar de hardcodeado como `"Human"`
- ✅ **Agent 6 (pathways.py):** CLI debe recibir `data/meta/` (directorio padre) no `data/meta/pancreatic_cancer/`

### Agregado
- **settings.yaml:** Campo `organism: "human"` bajo sección `pathways`

### Pipeline — Primera ejecución end-to-end completada (2026-03-08)
- 3 datasets procesados (GSE282985, GSE291806, GSE317494) | 129 muestras totales
- edgeR exitoso: 1,283 + 16 + 1,093 genes DEA significativos
- 6,552 genes consenso identificados (meta_padj < 0.05, ≥2/3 datasets)
- 593 pathways significativos (KEGG:71, GO:197, Reactome:295, MSigDB:30)
- 2 targets terapéuticos generados por Ollama gemma3:4b
- Reporte JSON generado con executive summary automático

### Corregido (Agent 7 y 8)
- ✅ **Agent 7 (insights.py):** `open(path)` → `open(path, encoding="utf-8")` para settings.yaml
- ✅ **Agent 7 (insights.py):** `data.get("genes", [])` → `data.get("consensus_genes", [])`
- ✅ **Agent 7 (insights.py):** `open(...biological_insights.json, "w")` → `open(..., "w", encoding="utf-8")` — evita cp1252 en Windows
- ✅ **Agent 6 (pathways.py):** CLI requiere `data/meta/` (directorio padre), no subdirectorio de enfermedad
- ✅ **Agent 8 (report.py):** `consensus_genes.get("genes", [])` → `consensus_genes.get("consensus_genes", [])` — n_genes ahora correcto (6,552)
- ✅ **Agent 8 (report.py):** `open(path)` → `open(path, encoding="utf-8")` en `_load_json_safe`

---

## DEUDA TÉCNICA Y MEJORAS PENDIENTES

### 🔴 Alta Prioridad (antes de validación científica)

**Agent 5 — Meta-análisis demasiado permisivo**
- 6,552 genes consenso de 29,891 (22%) es excesivo para un meta-análisis robusto
- Causa: con 3 datasets, "≥60% datasets" = presente en ≥2, umbral muy bajo
- Fix: ajustar `min_dataset_fraction` dinámicamente según N datasets, o usar umbral absoluto mínimo de 3 datasets
- Fix alternativo: revisar GSE291806 (diseño ratón, solo 16 genes DEA) — puede estar inflando overlap por ruido

**Agent 6 — Organism debería inferirse automáticamente**
- Actualmente hardcodeado en `settings.yaml` como `"human"`
- Correcto: Agent 1 detecta `organism_ch1` de GEO metadata → guarda en `metadata.json` → Agent 2 lo propaga → Agent 5 lo incluye en `consensus_genes.json` → Agent 6 lo lee automáticamente
- Impacto: pipeline multi-organismo fallará silenciosamente si el usuario no cambia settings.yaml manualmente

**Agent 4 — edgeR/limma aún caen a t-test en algunos casos**
- limma falla con `module 'limma' has no attribute 'model_matrix'` — el método `model_matrix` debe llamarse via `ro.r("model.matrix(...)")`
- edgeR funciona pero GSE291806 produce solo 16-19 genes significativos — revisar diseño experimental del dataset

**GSE291806 — Clasificación de muestras cuestionable**
- Dataset de ratón (ENSMUSG) con experimento de knockout de genes de cáncer pancreático
- Solo 16 genes DEA significativos — posiblemente el contraste case/control no es biológicamente apropiado
- Considerar excluirlo del meta-análisis o revisar manualmente la clasificación

### 🟡 Media Prioridad (v1.1.0)

**Agent 2 — Clasificación de muestras basada en keywords es frágil**
- Funciona para datasets estándar pero falla con nomenclaturas crípticas de líneas celulares
- Mejora: usar LLM (Agent 7) para clasificar muestras ambiguas basándose en el abstract del paper asociado
- Alternativa: interfaz de revisión manual para muestras `unclassified`

**Agent 3 — Solo guarda matrix_counts.csv para RNA-seq**
- `_process_microarray` retorna `None` como tercer elemento — consistente pero podría confundir
- Mejorar: guardar siempre ambos archivos con nombres descriptivos (`matrix_raw.csv`, `matrix_normalized.csv`)

**Agent 5 — Conversión de IDs de ratón usa símbolo de ratón como fallback**
- Si no hay ortholog humano disponible en MyGeneInfo, se usa el símbolo de ratón directamente
- Esto puede introducir genes incorrectos en el análisis humano
- Fix: drop genes sin ortholog humano confirmado cuando el dataset es de ratón

**Pipeline — Windows-specific issues**
- `"sh" no se reconoce` y `cffi mode "ANY"` son warnings de rpy2 en Windows — no bloquean pero ensucian el output
- Considerar suprimir estos warnings con filtro de logging al inicio

### 🟢 Baja Prioridad (v2.0.0+)

- Agent 1: agregar campo `organism` al `metadata.json` desde `organism_ch1` de GEO
- Agent 5: propagar `organism` al `consensus_genes.json`
- Agent 6: leer `organism` de `consensus_genes.json` en lugar de `settings.yaml`
- Cache de conversiones MyGeneInfo para evitar llamadas repetidas en re-runs
- Soporte multi-organismo en el mismo análisis (ej. human + mouse datasets juntos con ortholog mapping)

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