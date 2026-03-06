# 🤖 Documentación de Agentes — BioSignal Discovery Engine

Cada agente es un módulo Python independiente con interfaz estandarizada,
idempotente y re-ejecutable de forma aislada.

---

## Agent 1 — DatasetDiscoveryAgent

**Archivo:** `biosignal/agents/discovery.py`  
**Tiempo esperado:** 30-60 segundos  
**Dependencias externas:** NCBI E-utilities API

### Responsabilidad
Dado un término de enfermedad, identificar los datasets de expresión génica más relevantes en GEO/NCBI.

### Parámetros de Entrada

| Parámetro | Tipo | Default | Descripción |
|---|---|---|---|
| `disease_name` | str | — | Nombre canónico de la enfermedad |
| `max_datasets` | int | 20 | Máximo datasets a recuperar |
| `min_samples` | int | 10 | Mínimo muestras por dataset |
| `data_types` | list | `["RNA-seq", "microarray"]` | Tipos de datos aceptados |

### Proceso
1. Consultar NCBI E-utilities (`esearch` + `esummary`) con query: `'{disease} expression profiling'`
2. Filtrar por tipo: `expression profiling by high throughput sequencing` / `array`
3. Filtrar por `sample_count >= min_samples`
4. Rankear por: citas (40%) + muestras (40%) + año (20%)
5. Retornar lista rankeada de GSE IDs con metadata

### Salida
```
data/discovery/datasets_{disease_name}.json
```
```json
{
  "disease_name": "pancreatic cancer",
  "total_found": 18,
  "datasets": [
    {
      "gse_id": "GSE71989",
      "title": "Gene expression profiling of pancreatic cancer",
      "organism": "Homo sapiens",
      "sample_count": 45,
      "platform": "GPL570",
      "year": 2021,
      "pmid": "34521234",
      "data_type": "microarray",
      "rank_score": 0.87
    }
  ]
}
```

### CLI
```bash
python -m biosignal.agents.discovery \
  --disease 'pancreatic cancer' \
  --max-datasets 20 \
  --min-samples 10 \
  --output data/discovery/
```

---

## Agent 2 — DatasetDownloadAgent

**Archivo:** `biosignal/agents/download.py`  
**Tiempo esperado:** 10-45 minutos (según tamaño datasets)  
**Dependencias externas:** NCBI FTP, GEOparse

### Responsabilidad
Descargar archivos de datos crudos de GEO para cada GSE ID identificado.

### Parámetros de Entrada

| Parámetro | Tipo | Default | Descripción |
|---|---|---|---|
| `input_json` | str | — | Path al JSON del Agent 1 |
| `output_dir` | str | `data/raw/` | Directorio de almacenamiento |
| `parallel` | int | 4 | Workers paralelos de descarga |

### Proceso
1. Para cada GSE ID: descargar `{GSE_ID}_series_matrix.txt.gz`
2. Verificar integridad MD5
3. Descomprimir `.gz` automáticamente
4. Parsear metadata (GSM IDs, condición, plataforma)
5. Clasificar muestras como `control` vs `case` por keywords
6. Datasets con muestras sin clasificar → `data/review/` para revisión manual

### Clasificación Automática de Muestras

| Label | Keywords detectados |
|---|---|
| `control` | control, normal, healthy, non-tumor, adjacent, wild type |
| `case` | tumor, cancer, disease, patient, treated, knockdown, knockout |
| `unclassified` | keywords ambiguos o ausentes |

### Manejo de Errores

| Error | Acción |
|---|---|
| Timeout | Reintentar hasta 3 veces (backoff: 1s, 2s, 4s) |
| Dataset corrupto | Marcar como `failed`, continuar con siguientes |
| Sin clasificación clara | Mover a `data/review/` |

### Salida por Dataset
```
data/raw/{GSE_ID}/
├── matrix.tsv          # Matriz de expresión cruda (genes × muestras)
├── metadata.json       # Metadata de muestras con clasificación
└── platform.json       # Información de plataforma GPL
```

---

## Agent 3 — PreprocessingAgent

**Archivo:** `biosignal/agents/preprocess.py`  
**Tiempo esperado:** 2-10 min/dataset

### Responsabilidad
Normalizar, filtrar y realizar control de calidad de cada dataset independientemente.

### Sub-pipeline de Preprocesamiento

| Paso | Acción | Método |
|---|---|---|
| 1 | Carga de matriz cruda | Detección automática de formato (counts, TPM, RPKM, log2) |
| 2 | Filtración genes baja expresión | Remover genes con <10 counts en >70% de muestras |
| 3 | Normalización | TMM para RNA-seq; quantile para microarray |
| 4 | Detección outliers | PCA + Mahalanobis distance; z-score >3 marcados |
| 5 | Anotación genes | Probe IDs → Gene Symbols via bioDBnet |
| 6 | Detección batch effects | PVE analysis; si PVE batch >30%, aplicar ComBat |
| 7 | Export normalizado | log2-transformado, listo para DEA |

### Parámetros de Calidad Mínimos

| Parámetro | Threshold MVP | Acción si falla |
|---|---|---|
| Muestras por grupo (min) | ≥3 caso, ≥3 control | Excluir dataset del pipeline |
| Genes detectados | ≥5,000 | Warning + continuar |
| Muestras outlier | ≤20% del total | Remover outliers + warning |
| Batch PVE | >30% aplica ComBat | Corrección automática |

### Salida
```
data/processed/{GSE_ID}/
├── expr_normalized.parquet    # Matriz normalizada (genes × muestras)
├── qc_report.json             # Métricas QC
└── sample_labels.csv          # Clasificación final de muestras
```

---

## Agent 4 — DifferentialExpressionAgent

**Archivo:** `biosignal/agents/dea.py`  
**Tiempo esperado:** 1-5 min/dataset

### Responsabilidad
Ejecutar análisis de expresión diferencial (caso vs. control) en cada dataset procesado.

### Lógica de Selección de Método

| Tipo de dato | Método | Implementación |
|---|---|---|
| RNA-seq (counts) | DESeq2 | `pydeseq2` (nativo Python) o `rpy2` |
| Microarray (log2) | limma (voom) | `rpy2` |
| Ambiguo / desconocido | Wilcoxon rank-sum | `scipy.stats` |

### Pipeline DESeq2 (RNA-seq)
```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.stat_models import DeseqStats

dds = DeseqDataSet(counts=counts_df, metadata=meta_df, design_factors='condition')
dds.deseq2()
stat_res = DeseqStats(dds, contrast=['condition', 'case', 'control'])
stat_res.summary()
results = stat_res.results_df  # baseMean, log2FC, lfcSE, stat, pvalue, padj
```

### Criterios de Significancia
- `padj < 0.05` (FDR Benjamini-Hochberg)
- `|log2FoldChange| >= 1.0` (cambio mínimo 2x)
- `baseMean >= 10` (expresión mínima detectable)

### Salida
```
data/dea/{GSE_ID}/
├── dea_results.csv         # Tabla completa DEA (todos los genes)
├── significant_genes.json  # Solo genes significativos con FC y padj
└── volcano_plot.png        # Figura volcano plot
```

---

## Agent 5 — MetaAnalysisAgent ⭐

**Archivo:** `biosignal/agents/meta_analysis.py`  
**Tiempo esperado:** 1-3 minutos  
**Feature diferenciadora principal del sistema**

### Responsabilidad
Integrar resultados DEA de múltiples datasets para identificar señales biológicas **reproducibles y consistentes** entre estudios.

> 📌 Ninguna herramienta existente (Galaxy, GEO2R, Bioconductor) automatiza este paso de integración cross-dataset.

### Métodos Implementados

| Método | Descripción | Uso |
|---|---|---|
| **Fisher's combined probability** | Combina p-values independientes de cada estudio | Default |
| **Stouffer's Z-score** | Pondera por número de muestras | Con `--method stouffer` |
| **Vote counting** | Genes up/down en ≥N estudios de M totales | Complementario |
| **Random Effects Model** | Modela heterogeneidad entre estudios | Opcional, via rpy2+metafor |

### Pipeline de Integración
1. Cargar `significant_genes.json` de cada dataset
2. Construir matriz: genes (filas) × datasets (columnas) con log2FC
3. Alinear genes por Gene Symbol canónico
4. Calcular meta p-value por gen (Fisher's method)
5. Aplicar FDR correction (Benjamini-Hochberg) sobre meta p-values
6. Filtrar: genes con `meta_padj < 0.05` **Y** presentes en ≥60% de datasets
7. Calcular dirección consenso: UP si mean(log2FC) > 0 en ≥70% de datasets

### Salida
```
data/meta/{disease}/
├── meta_analysis_results.csv   # Tabla completa de meta-análisis
├── consensus_genes.json        # Genes top con score de consistencia
└── heatmap_top50.png           # Heatmap cross-dataset
```

---

## Agent 6 — PathwayEnrichmentAgent

**Archivo:** `biosignal/agents/pathways.py`  
**Tiempo esperado:** 2-5 minutos

### Responsabilidad
Identificar pathways biológicos significativamente alterados en los genes consenso.

### Bases de Datos

| Base de Datos | Tipo | Herramienta |
|---|---|---|
| KEGG | Pathway ORA | `gseapy.enrichr` — `KEGG_2021_Human` |
| Gene Ontology (GO) | BP, MF, CC enrichment | `gseapy.enrichr` — `GO_Biological_Process_2023` |
| Reactome | Pathway ORA | `gseapy.enrichr` — `Reactome_2022` |
| MSigDB Hallmarks | Cancer hallmarks | `gseapy.enrichr` — `MSigDB_Hallmark_2020` |

### Proceso
1. Separar consensus_genes en upregulated / downregulated
2. ORA (Over-Representation Analysis) para cada lista
3. Filtrar: `adjusted_pval < 0.05` y `gene_ratio >= 0.05`
4. Rankear por: `combined_score = log(pval) × z-score`
5. Generar dot plot de top 20 pathways

---

## Agent 7 — InsightGenerationAgent

**Archivo:** `biosignal/agents/insights.py`  
**Tiempo esperado:** 3-10 minutos (LLM local), 1-3 min (API)

### Responsabilidad
Usar un LLM para interpretar resultados biológicos y generar hipótesis terapéuticas.

### LLM Backends

| Backend | Configuración | Velocidad |
|---|---|---|
| **Ollama llama3:8b** | `LLM_BACKEND=local` | ~10 min/query (CPU) |
| **Ollama llama3:70b** | GPU RTX 4090 requerida | ~2 min/query |
| **Claude 3.5 Sonnet** | `ANTHROPIC_API_KEY=...` | ~30 seg/query |
| **GPT-4o** | `OPENAI_API_KEY=...` | ~30 seg/query |

### Herramientas del Agente (Tool Calling)

| Tool | Función |
|---|---|
| `pubmed_search` | Busca papers sobre genes/pathways identificados |
| `uniprot_query` | Obtiene función de proteína e interacciones |
| `chembl_search` | Identifica compuestos que modulan los targets |
| `pathway_context` | Recupera descripción biológica de pathways |

### Output Generado
1. Top 5 targets terapéuticos con rationale biológico
2. Mecanismos clave dysregulados en la enfermedad
3. Fármacos/compuestos conocidos contra estos genes (ChEMBL)
4. 3 hipótesis experimentales testables
5. Nivel de confianza por hallazgo (High/Medium/Low)

---

## Agent 8 — ReportGenerationAgent

**Archivo:** `biosignal/agents/report.py`  
**Tiempo esperado:** 1-3 minutos

### Responsabilidad
Consolidar todos los outputs en un reporte final estructurado y autoexplicativo.

### Formatos de Salida

| Formato | Descripción |
|---|---|
| **PDF** | Reporte completo con figuras, tablas y texto (via ReportLab) |
| **JSON** | Todos los resultados estructurados en un archivo indexado |
| **CSV** | Tablas de genes y pathways para Excel/R/Python |
| **Markdown** | Reporte legible en texto plano |

### Estructura del Reporte PDF

1. **Resumen ejecutivo** (1 pág): disease, datasets analizados, top findings
2. **Metodología**: datasets, QC summary, métodos estadísticos
3. **Resultados — Genes**: tabla top 50 con FC, padj, consistencia
4. **Resultados — Pathways**: dot plot KEGG + GO + Reactome
5. **Visualizaciones**: heatmap, PCA, volcano plots
6. **Targets terapéuticos**: gen, función, compuestos conocidos, evidencia
7. **Hipótesis IA**: con nivel de confianza y referencias
8. **Anexo**: datasets, parámetros, versiones de software