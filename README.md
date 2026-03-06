# 🧬 BioSignal Discovery Engine
### *Plataforma de Descubrimiento de Señales Biológicas basada en Transcriptómica Multi-Dataset e IA*

> **Versión:** 1.0 — MVP Local  
> **Estado:** Prototipo funcional sin UI  
> **Clasificación:** Confidencial — Uso Interno  
> **Fecha:** Marzo 2026

---

## 📋 Descripción

BioSignal Discovery Engine automatiza el análisis transcriptómico multi-dataset para identificar genes, pathways y targets terapéuticos consistentes entre múltiples estudios públicos. El sistema opera completamente en local, construido sobre un stack Python modular con agentes de IA orquestados.

**Propuesta de valor única:**  
El usuario escribe `'pancreatic cancer'` y el sistema automáticamente descarga, normaliza, analiza e integra múltiples datasets de GEO, entregando genes consistentemente dysregulados, pathways clave y targets terapéuticos candidatos — **en horas, no semanas**.

---

## 🚨 El Problema que Resuelve

| Problema | Impacto |
|----------|---------|
| 80% del tiempo de investigación se pierde limpiando e integrando datos | Pérdida masiva de productividad |
| Los investigadores biológicos no dominan bioinformática computacional | Barrera de acceso al análisis |
| Las herramientas existentes (Galaxy, GEO2R, Bioconductor) producen datos pero no insights | Falta de interpretación accionable |
| No existe plataforma que analice múltiples datasets en forma cruzada automáticamente | Gap de mercado crítico |

---

## 🏗️ Arquitectura del Sistema

```
INPUT: disease_name (string)
       │
[Agent 1] DatasetDiscoveryAgent   → Identifica datasets GEO relevantes
       │
[Agent 2] DatasetDownloadAgent    → Descarga SOFT/MATRIX files de NCBI
       │
[Agent 3] PreprocessingAgent      → Normaliza, filtra y QC por dataset
       │
[Agent 4] DifferentialExpressionAgent → DEA por dataset (DESeq2/limma/edgeR)
       │
[Agent 5] MetaAnalysisAgent       → Integra DEA cross-dataset (Fisher, REM)
       │
[Agent 6] PathwayEnrichmentAgent  → KEGG, GO, Reactome enrichment
       │
[Agent 7] InsightGenerationAgent  → LLM interpreta señales biológicas
       │
[Agent 8] ReportGenerationAgent   → PDF + JSON + CSV
       │
OUTPUT: biological_insight_report
```

---

## 📁 Estructura del Repositorio

```
biosignal-discovery/
├── README.md                       # Este archivo
├── CHANGELOG.md                    # Historial de cambios
├── requirements.txt                # Dependencias Python pinadas
├── requirements_r.txt              # Paquetes R requeridos
├── environment.yml                 # Entorno conda completo
├── install_r_packages.R            # Script instalación R
├── .env.example                    # Variables de entorno ejemplo
│
├── config/
│   ├── settings.yaml               # Parámetros globales del pipeline
│   └── disease_aliases.json        # Mapeo nombres enfermedades → términos GEO
│
├── biosignal/
│   ├── __init__.py
│   ├── cli.py                      # Interfaz de línea de comandos principal
│   ├── pipeline.py                 # Orquestador del pipeline completo
│   ├── agents/
│   │   ├── __init__.py
│   │   ├── discovery.py            # Agent 1: DatasetDiscoveryAgent
│   │   ├── download.py             # Agent 2: DatasetDownloadAgent
│   │   ├── preprocess.py           # Agent 3: PreprocessingAgent
│   │   ├── dea.py                  # Agent 4: DifferentialExpressionAgent
│   │   ├── meta_analysis.py        # Agent 5: MetaAnalysisAgent
│   │   ├── pathways.py             # Agent 6: PathwayEnrichmentAgent
│   │   ├── insights.py             # Agent 7: InsightGenerationAgent
│   │   └── report.py               # Agent 8: ReportGenerationAgent
│   ├── utils/
│   │   ├── __init__.py
│   │   ├── geo_utils.py            # Funciones de acceso a GEO/NCBI
│   │   ├── stats_utils.py          # Funciones estadísticas comunes
│   │   ├── gene_utils.py           # Mapeo y normalización de gene IDs
│   │   └── viz_utils.py            # Generación de figuras
│   └── schemas/
│       ├── __init__.py
│       ├── dataset.py              # Pydantic schema para dataset
│       ├── dea_result.py           # Pydantic schema para resultados DEA
│       └── report.py               # Pydantic schema para reporte final
│
├── data/
│   ├── discovery/                  # Outputs Agent 1 (datasets JSON)
│   ├── raw/                        # Outputs Agent 2 (archivos GEO)
│   ├── processed/                  # Outputs Agent 3 (matrices normalizadas)
│   ├── dea/                        # Outputs Agent 4 (resultados DEA)
│   ├── meta/                       # Outputs Agent 5 (meta-análisis)
│   ├── pathways/                   # Outputs Agent 6 (enrichment)
│   ├── insights/                   # Outputs Agent 7 (insights LLM)
│   └── reports/                    # Outputs Agent 8 (reportes finales)
│
├── tests/
│   ├── __init__.py
│   ├── test_discovery.py
│   ├── test_preprocess.py
│   ├── test_meta_analysis.py
│   └── conftest.py
│
├── notebooks/
│   └── exploration.ipynb           # Exploración y debugging interactivo
│
└── docs/
    ├── ARCHITECTURE.md             # Arquitectura detallada
    ├── AGENTS.md                   # Documentación de cada agente
    ├── API.md                      # Referencia de la CLI
    ├── INSTALLATION.md             # Guía de instalación paso a paso
    ├── TROUBLESHOOTING.md          # Problemas frecuentes y soluciones
    └── ROADMAP.md                  # Hoja de ruta del proyecto
```

---

## ⚡ Inicio Rápido

### 1. Instalación

```bash
# Clonar repositorio
git clone https://github.com/tu-org/biosignal-discovery.git
cd biosignal-discovery

# Crear entorno conda
conda env create -f environment.yml
conda activate biosignal

# Instalar paquetes R
Rscript install_r_packages.R

# Configurar variables de entorno
cp .env.example .env
# Editar .env con tus API keys
```

### 2. Ejecución del Pipeline Completo

```bash
python -m biosignal.cli run \
  --disease 'pancreatic cancer' \
  --max-datasets 15 \
  --min-samples 10 \
  --parallel 4 \
  --output-dir data/ \
  --report-format pdf,json,csv
```

### 3. Ejecución por Agente Individual

```bash
python -m biosignal.agents.discovery --disease 'alzheimer disease'
python -m biosignal.agents.download  --input data/discovery/
python -m biosignal.agents.preprocess --input data/raw/
python -m biosignal.agents.dea       --input data/processed/
python -m biosignal.agents.meta_analysis --input data/dea/
python -m biosignal.agents.pathways  --input data/meta/
python -m biosignal.agents.insights  --input data/pathways/
python -m biosignal.agents.report    --input data/insights/
```

---

## 💻 Requerimientos de Hardware

| Componente | Mínimo (MVP) | Recomendado (Producción) |
|------------|--------------|--------------------------|
| CPU | 4 cores / 8 threads | 16+ cores |
| RAM | 16 GB | 64 GB+ |
| Almacenamiento | 100 GB SSD | 1 TB NVMe SSD |
| GPU (LLM local) | Opcional (CPU mode) | NVIDIA RTX 4090 / A100 |
| Red | 10 Mbps+ | 100 Mbps+ |
| OS | Ubuntu 22.04 / macOS 13+ | Ubuntu 22.04 LTS |

> 💡 Para LLM local con GPU RTX 4090 (24GB VRAM): llama3:70b en 4-bit quantization.  
> Sin GPU: usar llama3:8b en CPU (~10 min/query) o API de Claude/OpenAI.

---

## 📊 SLA del MVP Local

| Métrica | Target | Condición |
|---------|--------|-----------|
| Disponibilidad del pipeline | 99% (local) | No depende de servidores externos |
| Tasa de éxito por dataset | ≥ 80% | Datasets con QC aprobado |
| Reproducibilidad | 100% | Mismo input → mismo output |
| Tiempo total (10 datasets) | ≤ 6 horas | Hardware mínimo |
| Genes consenso identificados | ≥ 50 genes | Para enfermedades con ≥5 datasets en GEO |
| Pathways significativos | ≥ 10 pathways | padj < 0.05, ≥2 de 4 bases de datos |

---

## 🗺️ Hoja de Ruta

| Fase | Timeline | Entregables | Criterio de Éxito |
|------|----------|-------------|-------------------|
| **MVP Local** | Mes 1-2 | Pipeline funcional, 8 agentes, reporte PDF | Análisis completo de cáncer pancreático con 10+ datasets |
| Validación | Mes 3 | Comparación con papers publicados | Reproducir genes de 3 papers conocidos |
| **API REST** | Mes 4-5 | FastAPI backend, autenticación | 10 requests concurrentes sin degradación |
| UI Web | Mes 6-8 | Interface usuario, visualizaciones | Usuario sin bioinformática produce reporte en <5 clicks |
| **Multi-omics** | Mes 9-12 | Proteómica, GWAS, datos clínicos | Knowledge graph multi-omics funcional |
| SaaS Beta | Mes 13-15 | Producto comercial | 5 clientes pagos (academia o biotech) |

---

## 📚 Documentación

- [Arquitectura Detallada](docs/ARCHITECTURE.md)
- [Guía de Instalación](docs/INSTALLATION.md)
- [Documentación de Agentes](docs/AGENTS.md)
- [Referencia de la CLI](docs/API.md)
- [Solución de Problemas](docs/TROUBLESHOOTING.md)

---

## 🔗 Fuentes y Bases de Datos

| Recurso | URL | Uso |
|---------|-----|-----|
| NCBI GEO | ncbi.nlm.nih.gov/geo | Fuente principal de datasets |
| NCBI E-utilities | eutils.ncbi.nlm.nih.gov | API de búsqueda programática |
| UniProt | uniprot.org | Función de proteínas e interacciones |
| ChEMBL | ebi.ac.uk/chembl | Compuestos bioactivos y targets |
| KEGG | genome.jp/kegg | Pathways biológicos |
| Gene Ontology | geneontology.org | Términos funcionales de genes |
| Reactome | reactome.org | Pathways moleculares curados |
| MSigDB | gsea-msigdb.org | Colecciones de gene sets |
| Enrichr | maayanlab.cloud/Enrichr | API de pathway enrichment |
| PubMed | pubmed.ncbi.nlm.nih.gov | Literatura científica |
| Ollama | ollama.com | Servidor LLM local |

---

## 📄 Licencia

Confidencial — Uso Interno. Ver LICENSE.txt para términos completos.

---

*BioSignal Discovery Engine v1.0 | Documento Confidencial | Marzo 2026*