# Changelog — BioSignal Discovery Engine

Todos los cambios notables del proyecto serán documentados en este archivo.

Formato basado en [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [1.0.3] — 2026-03-10 — Aislamiento por Enfermedad y Calidad de Datos

### Corregido
- ✅ **Agent 6 (pathways.py):** Flag `--disease` añadido al CLI — ahora solo enriquece la enfermedad especificada, no todas las carpetas en `data/meta/`
- ✅ **Agent 7 (insights.py):** Flag `--disease` añadido al CLI — evita procesar análisis de otras enfermedades
- ✅ **Agent 8 (report.py):** Ya tenía `--disease` — confirmado correcto ✅
- ✅ **Agent 3 (preprocess.py):** `effective_fraction` para filtro de baja expresión ahora configurable via `settings.yaml` — elimina lógica hardcodeada `max(0.1, 20.0/n_samples)`

### Agregado
- **settings.yaml:** Parámetros `low_expression_fraction_large` (default 0.20) y `large_dataset_threshold` (default 100) bajo sección `preprocessing`

### Comportamiento actual
- Datasets ≤100 muestras → fracción 0.70 (estricto)
- Datasets >100 muestras → fracción 0.20 (permisivo, configurable)
- Ambos valores ajustables desde `settings.yaml` sin tocar código

---

## [1.0.2] — 2026-03-09 — Pipeline Alzheimer + Fixes Críticos

### Corregido
- ✅ **Agent 1 (discovery.py):** Nuevo método `_check_downloadability(gse_id)` — consulta FTP NCBI antes de rankear; score 3=counts consolidados, 2=RAW.tar, 1=series matrix, 0=sin archivos; pesa 15% en ranking final
- ✅ **schemas/dataset.py:** Campo `downloadability_score: int = Field(default=1)` añadido
- ✅ **Agent 3 (preprocess.py) — Capa 2:** Nuevo método `_llm_map_columns()` — usa Ollama gemma3:4b para mapear columnas ambiguas a grupos caso/control cuando capas determinísticas fallan
- ✅ **Agent 3 (preprocess.py) — Capa 3:** Nuevo método `_sra_map_columns()` — consulta SRA RunInfo + BioSample para recuperar metadata cuando LLM falla
- ✅ **Agent 3 (preprocess.py):** Fix cache — busca en `out_dir` primero, luego `dataset_dir`; si índice no coincide → borra y continúa
- ✅ **Agent 3 (preprocess.py):** `_process_microarray` ahora retorna 3 valores (`matrix, qc, None`)
- ✅ **Agent 3 (preprocess.py):** Detección de datos pre-normalizados usa percentil p99 en lugar de max
- ✅ **Agent 4 (dea.py):** `_select_method` asigna limma (no edgeR) para datos de microarray
- ✅ **Agent 5 (meta_analysis.py):** Fix tipos mixtos int/str — `sorted(str(g) for g in all_genes)`
- ✅ **Agent 2 (download.py):** Keywords expandidos — CONTROL: `ctr, ctrl, non-ad`; CASE: `alzheimer, ad, mci, dementia`

### Pipeline — Alzheimer ejecutado end-to-end (2026-03-09)
- 2 datasets con señal real (GSE226507 humano, GSE294900 murino Sox9-KO)
- 4,031 genes consenso | 143 pathways significativos
- Genes clave: TREM2, CX3CR1, CSF1, CTSB, SPARC, COL6A1/COL6A2

---

## [1.0.1] — 2026-03-08 — MVP Debugging Session

### Corregido
- ✅ **Agent 2 (download.py):** Clasificación de muestras — keywords con prioridad por longitud
- ✅ **Agent 3 (preprocess.py):** `_load_sample_metadata` acepta `matrix_cols` para alinear índices
- ✅ **Agent 3 (preprocess.py):** Guarda `matrix_counts.csv` para edgeR
- ✅ **Agent 4 (dea.py):** Alineación matrix/sample_meta antes de llamar a R
- ✅ **Agent 4 (dea.py):** Subsetting DGEList via `ro.r("dge[keep,,keep.lib.sizes=FALSE]")`
- ✅ **Agent 5 (meta_analysis.py):** Normalización cross-dataset via `mygene` (ENSG→symbol, ENSMUSG→human ortholog)
- ✅ **Agent 6 (pathways.py):** Clave `consensus_genes` en JSON; `organism` configurable
- ✅ **Agent 7/8:** Fixes encoding UTF-8, claves JSON correctas

### Pipeline — Primera ejecución end-to-end (2026-03-08)
- GSE282985, GSE291806, GSE317494 | 129 muestras
- 6,552 genes consenso | 593 pathways | Targets: HIPK2, MARCKS

---

## DEUDA TÉCNICA Y MEJORAS PENDIENTES

### 🔴 Alta Prioridad

**[CALIDAD DE DATOS] No existe estándar para definir qué datasets son "buenos"**

El pipeline actualmente decide internamente con reglas fijas si un dataset es aceptable. Esto genera variabilidad entre runs: datasets distintos pueden entrar/salir del análisis sin que el usuario lo sepa ni lo controle.

**Solución propuesta (v1.1.0):** Scorecard de calidad por dataset, visible al usuario antes de continuar el pipeline. El usuario decide qué incluir — el sistema recomienda pero no impone.

**Parámetros clave de calidad a exponer:**

| Nivel | Parámetro | Umbral sugerido | Por qué importa |
|-------|-----------|-----------------|-----------------|
| Dataset | Muestras por grupo | Mínimo 3, ideal ≥6 | DEA estadísticamente confiable |
| Dataset | Balance caso/control | Ratio ≤5:1 | Evita sesgo estadístico |
| Dataset | Tipo de tejido | Mismo entre datasets | Meta-análisis válido |
| Dataset | Organismo | Humano preferido | Interpretabilidad clínica |
| Expresión | % genes detectables | >30% del genoma | Proxy de profundidad de secuenciación |
| Expresión | Correlación entre réplicas | >0.85 | Reproducibilidad del experimento |
| Expresión | Varianza intra-grupo | Baja relativa | Muestras homogéneas, sin mezcla de condiciones |
| DEA | Genes significativos | Entre 10 y 8,000 | Señal biológica real (muy pocos o muchos = problema) |
| DEA | Balance UP/DOWN | Sin extremos >90% | Descarta artefactos técnicos |
| DEA | Distribución p-values | Uniforme con pico en 0 | DEA estadísticamente correcto |

**Agent 5 — Meta-análisis poco robusto con <3 datasets**
- Fix: requerir mínimo absoluto de 3 datasets; ajustar `min_dataset_fraction` dinámicamente

**Agent 6 — 0 pathways con <50 genes consenso**
- Enrichr requiere ~50-100 genes mínimo; advertir al usuario cuando lista es insuficiente

**Agent 7 — JSON malformado de Ollama**
- gemma3:4b genera JSON inválido en respuestas largas — fix: retry con temperatura reducida

### 🟡 Media Prioridad

- **Agent 2:** Clasificación de muestras frágil con nomenclaturas crípticas — interfaz revisión manual
- **Agent 6:** `organism` debería inferirse automáticamente desde GEO metadata
- **Agent 4:** limma falla en Windows — `model.matrix` debe llamarse via `ro.r()`
- **report.py:** `datetime.utcnow()` deprecated — reemplazar por `datetime.now(datetime.UTC)`

### 🟢 Baja Prioridad

- Capa 3 SRA: batching de consultas BioSample (actualmente individual = lento)
- Cache de conversiones MyGeneInfo para re-runs
- Directorio `unknown` huérfano en `data/meta/` procesado innecesariamente por Agent 6

---

## [1.0.0] — 2026-03-01 — MVP Local

### Agregado
- 8 agentes del pipeline (Discovery → Report)
- CLI completa, schemas Pydantic v2, suite de tests
- Soporte: Python 3.11+, R 4.3+, Windows/macOS/Linux
- LLM backends: Ollama, Anthropic Claude, OpenAI GPT-4

---

## Próximas Versiones

### [1.1.0] — Scorecard de Calidad de Datasets
- Scorecard por dataset con los 10 parámetros identificados
- Revisión manual antes de DEA — el usuario decide qué incluir
- Recomendación automática con justificación (aceptar / advertir / rechazar)
- Flag `--interactive` en CLI

### [2.0.0] — API REST
- FastAPI + JWT, Docker Compose, endpoints `/analyze`, `/status`, `/results`

### [3.0.0] — UI Web
- React + visualizaciones interactivas (volcano plots, heatmaps, pathway maps)

### [4.0.0] — Multi-omics
- Proteómica, GWAS, datos clínicos, knowledge graph Neo4j

### [5.0.0] — SaaS Beta
- Pricing: Academic / Biotech / Enterprise | SOC2 compliance