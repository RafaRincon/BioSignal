# BioSignal Discovery Engine — Arquitectura de Trazabilidad y Calidad de Datos
**Versión:** 1.1.0-design  
**Fecha:** 2026-03-09  
**Estado:** Diseño aprobado — pendiente implementación por fases

---

## Problema

El pipeline actual toma decisiones de calidad de forma silenciosa y dispersa. No existe un momento donde el usuario pueda ver qué datasets entraron al análisis, por qué se descartaron otros, o cómo el sistema llegó a un gen target. Dos runs con los mismos parámetros pueden producir resultados distintos sin explicación visible.

**Principio rector:** Cada decisión del pipeline debe responder tres preguntas:
- **¿Qué** decidió el sistema?
- **¿Por qué** lo decidió (criterio + valor observado)?
- **¿Con qué datos** lo decidió?

---

## Estructura de Archivos por Run

Cada ejecución del pipeline genera una carpeta de auditoría completa:

```
data/runs/{run_id}/
├── run_manifest.json          ← resumen ejecutivo del run completo
└── audit/
    ├── 01_discovery.json
    ├── 02_download.json
    ├── 03_preprocess.json
    ├── 04_dea.json
    ├── 05_meta.json
    ├── 06_pathways.json
    ├── 07_insights.json
    └── 08_report.json
```

El audit trail es **siempre activo** — no requiere flag. Cada run genera su propio directorio identificado por `{topic}_{timestamp}`.

---

## Esquema Base — `AuditEntry`

Todos los archivos de auditoría comparten esta estructura raíz:

```json
{
  "agent_id": "03_preprocess",
  "agent_name": "PreprocessingAgent",
  "run_id": "alzheimer_20260309_143022",
  "timestamp": "2026-03-09T14:32:11Z",
  "duration_seconds": 47.3,
  "data_type": "transcriptomics",
  "config_snapshot": { },
  "decisions": [ ],
  "warnings": [ ],
  "errors": [ ]
}
```

El campo `data_type` aparece en cada `AuditEntry` para soportar futuros tipos de dato (proteómica, epigenómica, multi-omics) sin cambiar el esquema.

### Estructura de una `decision`

Cada decisión registra la evidencia observada vs el umbral configurado:

```json
{
  "subject": "GSE226507",
  "decision": "ACCEPT",
  "reason": "Todos los criterios de calidad superados",
  "evidence": {
    "samples_per_group": { "observed": 8, "threshold": 3, "pass": true },
    "balance_ratio":     { "observed": 1.0, "threshold": 5.0, "pass": true },
    "genes_detected":    { "observed": 15420, "threshold": 5000, "pass": true },
    "outlier_fraction":  { "observed": 0.08, "threshold": 0.20, "pass": true }
  }
}
```

---

## Detalle por Agente

### Agent 1 — Discovery (`01_discovery.json`)

```
decisions:
  - Por cada dataset candidato:
      score total, score por dimensión (relevancia, downloadability, fecha),
      razón de inclusión o descarte del top-N
  - Query GEO usada, n resultados crudos vs n seleccionados

warnings:
  - Datasets con downloadability_score = 0
  - Datasets sin fecha de publicación
```

### Agent 2 — Download (`02_download.json`)

```
decisions:
  - Por cada muestra:
      label asignado (case / control / unknown),
      keyword que hizo match,
      fuente del label (title / description / sra)
  - Por cada dataset:
      n_case, n_control, n_unknown, método de descarga

warnings:
  - Lista explícita de muestras clasificadas como unknown
  - Datasets con >20% unknowns
```

### Agent 3 — Preprocess (`03_preprocess.json`) ← más crítico

```
decisions:
  - Scorecard completa por dataset (10 parámetros vs umbrales):
      samples_per_group, balance_ratio, tissue_type, organism,
      pct_genes_detected, replicate_correlation, intragroup_variance,
      n_sig_genes, up_down_balance, pvalue_distribution
  - PASS / WARN / FAIL con valor observado por métrica
  - Método de normalización aplicado y criterio de selección
  - n genes filtrados y razón

warnings:
  - Métricas borderline (pasan pero <110% del umbral mínimo)
  - Datasets con correlación entre réplicas 0.85–0.90
```

**Parámetros de la scorecard:**

| Parámetro | Umbral sugerido | Por qué importa |
|-----------|-----------------|-----------------|
| Muestras por grupo | Mínimo 3, ideal ≥6 | DEA estadísticamente confiable |
| Balance caso/control | Ratio ≤5:1 | Evita sesgo estadístico |
| Tipo de tejido | Mismo entre datasets | Meta-análisis válido |
| Organismo | Humano preferido | Interpretabilidad clínica |
| % genes detectables | >30% del genoma | Proxy de profundidad de secuenciación |
| Correlación entre réplicas | >0.85 | Reproducibilidad del experimento |
| Varianza intra-grupo | Baja relativa | Muestras homogéneas |
| Genes significativos | Entre 10 y 8,000 | Señal biológica real |
| Balance UP/DOWN | Sin extremos >90% | Descarta artefactos técnicos |
| Distribución p-values | Uniforme con pico en 0 | DEA estadísticamente correcto |

### Agent 4 — DEA (`04_dea.json`)

```
decisions:
  - Método elegido (edgeR / limma / t-test) y criterio de selección
  - Si hubo fallback: qué falló y por qué se activó el siguiente
  - n genes significativos, ratio UP/DOWN
  - Distribución p-values: % en [0–0.05], [0.05–0.5], [0.5–1.0]

warnings:
  - Ratio UP/DOWN > 80/20 (posible artefacto técnico)
  - n genes sig < 10 o > 8,000
  - Fallback activado (edgeR falló → limma)
```

### Agent 5 — Meta-analysis (`05_meta.json`)

```
decisions:
  - Por cada gen consenso:
      n datasets que lo reportaron,
      efecto promedio ponderado,
      dirección consistente o inconsistente entre datasets
  - min_dataset_fraction usado y justificación
  - n datasets incluidos vs descartados

warnings:
  - Genes con dirección inconsistente (ej. UP en 2 datasets, DOWN en 1)
  - n datasets < 3 (meta-análisis estadísticamente débil)
```

### Agent 6 — Pathways (`06_pathways.json`)

```
decisions:
  - Por cada pathway significativo:
      FDR, n genes input que lo sustentan,
      % de cobertura del pathway total,
      base de datos fuente (KEGG / Reactome / Enrichr)

warnings:
  - n genes input < 50 (Enrichr poco confiable)
  - Pathways muy genéricos con cobertura >40%
```

### Agent 7 — Insights (`07_insights.json`)

```
decisions:
  - Prompt exacto enviado al LLM (reproducible)
  - Modelo + versión + temperatura usados
  - n tokens input / output
  - target_classes evaluadas según catálogo

warnings:
  - JSON inválido detectado + n retries realizados
  - Respuesta truncada por límite de tokens
```

---

## `run_manifest.json` — Vista Ejecutiva

```json
{
  "run_id": "alzheimer_20260309_143022",
  "pipeline_version": "1.1.0",
  "analysis_context": {
    "domain": "disease",
    "topic": "alzheimer",
    "organism": "human",
    "data_type": "transcriptomics",
    "objective": "therapeutic_target_discovery",
    "target_classes": ["therapeutic", "biomarker_diagnostic"]
  },
  "started_at": "2026-03-09T14:30:00Z",
  "completed_at": "2026-03-09T15:12:44Z",
  "config_hash": "sha256:a3f9...",
  "datasets": {
    "discovered": 12,
    "downloaded": 5,
    "passed_qc": 3,
    "used_in_meta": 3
  },
  "results_summary": {
    "consensus_genes": 4031,
    "pathways": 143,
    "top_targets": {
      "therapeutic": ["TREM2", "CX3CR1", "CTSB"],
      "biomarker_diagnostic": ["SPARC", "COL6A1"]
    }
  },
  "warnings_total": 4,
  "audit_files": [
    "01_discovery.json",
    "02_download.json",
    "03_preprocess.json",
    "04_dea.json",
    "05_meta.json",
    "06_pathways.json",
    "07_insights.json",
    "08_report.json"
  ]
}
```

El campo `config_hash` es determinístico: dos runs con el mismo hash sobre los mismos datos producen resultados idénticos. Es la garantía de reproducibilidad.

---

## Versatilidad — Extensión sin Romper el Esquema

### `analysis_context` — Intención declarada

El campo `analysis_context` en el manifest declara *para qué* se ejecuta el pipeline. Permite que el mismo sistema analice transcriptómica de enfermedad, resistencia a antibióticos, aging, o respuesta ambiental sin cambiar código.

Campos clave:
- `domain`: `disease` | `microbiology` | `aging` | `environmental` | `pharmacology` | ...
- `objective`: `therapeutic_target_discovery` | `biomarker_discovery` | `mechanism_elucidation` | ...
- `data_type`: `transcriptomics` | `proteomics` | `epigenomics` | `multi_omics`
- `target_classes`: lista de clases activas del catálogo

### `target_classes` — Catálogo extensible en `settings.yaml`

Los tipos de targets son configurables sin tocar código:

```yaml
target_classes:
  therapeutic:
    description: "Gen druggable con evidencia funcional"
    evidence_required: [differential_expression, pathway_enrichment]

  biomarker_diagnostic:
    description: "Gen con expresión consistente entre datasets"
    evidence_required: [consensus_fraction > 0.8]

  transcription_factor:
    description: "Regulador maestro upstream"
    evidence_required: [tf_database_hit, high_connectivity]

  biomarker_prognostic:
    description: "Gen asociado a evolución clínica"
    evidence_required: [survival_correlation]

  crispr_target:
    description: "Gen esencial con baja toxicidad off-target"
    evidence_required: [essentiality_score]
```

El Agent 7 (insights) y el reporte final leen este catálogo — el reporte muestra secciones solo para las clases activas en el run.

### `data_type` como dimensión de primer nivel

| `data_type` | Fuente de datos | Agent 4 usa |
|-------------|-----------------|-------------|
| `transcriptomics` | GEO RNA-seq / microarray | edgeR / limma |
| `proteomics` | PRIDE / PhosphoSitePlus | limma / t-test |
| `epigenomics` | GEO ATAC-seq / ChIP-seq | DESeq2 especializado |
| `multi_omics` | combinado | por definir |

El audit de Agent 4 registra `data_type` — cuando cambie el método estadístico, la trazabilidad explica el criterio.

---

## Límites Deliberados (No Over-Engineering)

Lo que **no** se abstrae en esta versión:

- Los agentes no se convierten en clases genéricas intercambiables
- No se implementa un sistema de plugins
- Proteómica y epigenómica no se soportan en código — solo el esquema no las impide
- `settings.yaml` es el único punto de extensión

---

## Fases de Implementación

| Versión | Alcance | Valor entregado |
|---------|---------|-----------------|
| **1.1.0** | Esquema base `AuditEntry` + Agent 3 scorecard completa | Abre la caja más opaca del pipeline |
| **1.2.0** | Agent 1 + Agent 2 | Trazabilidad desde el origen de los datos |
| **1.3.0** | Agent 4 + Agent 5 | Confianza en DEA y consenso de genes |
| **1.4.0** | Agent 6 + Agent 7 + `run_manifest.json` completo | Trail completo, base para UI v3.0.0 |

---

## Relación con Roadmap General

```
v1.1.0 – v1.4.0  →  Audit Trail (este documento)
v2.0.0           →  API REST — los JSON de auditoría se sirven como endpoints
v3.0.0           →  UI Web — los JSON se convierten en vistas interactivas
v4.0.0           →  Multi-omics — data_type y target_classes ya preparados
v5.0.0           →  SaaS — audit trail como evidencia auditable para clientes enterprise
```

---

*BioSignal Discovery Engine — Documento de arquitectura interno*  
*No distribuir sin autorización*