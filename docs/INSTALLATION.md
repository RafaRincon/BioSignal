# 📦 Guía de Instalación — BioSignal Discovery Engine

## Requisitos Previos

| Componente | Versión Mínima | Verificar |
|---|---|---|
| Python | 3.11+ | `python --version` |
| R | 4.3+ | `R --version` |
| Conda/Mamba | cualquier | `conda --version` |
| Git | 2.x+ | `git --version` |
| RAM | 16 GB | — |
| Disco | 100 GB libres | `df -h` |

---

## Paso 1 — Clonar el Repositorio

```bash
git clone https://github.com/tu-org/biosignal-discovery.git
cd biosignal-discovery
```

---

## Paso 2 — Crear el Entorno Conda

```bash
# Crear entorno desde environment.yml
conda env create -f environment.yml

# Activar entorno
conda activate biosignal

# Verificar activación
which python  # Debe apuntar al entorno conda
```

> 💡 **Tip:** Si prefieres pip puro sin conda:
> ```bash
> python -m venv .venv
> source .venv/bin/activate  # Linux/macOS
> # .venv\Scripts\activate   # Windows
> pip install -r requirements.txt
> ```

---

## Paso 3 — Instalar Paquetes R

```bash
# Desde el directorio del proyecto, con el entorno conda activado
Rscript install_r_packages.R
```

El script instalará automáticamente:
- **Bioconductor:** DESeq2, limma, edgeR, GEOquery, sva, clusterProfiler
- **CRAN:** tidyverse, metafor, data.table
- **Annotation:** org.Hs.eg.db, biomaRt, AnnotationDbi

> ⚠️ La primera instalación puede tardar **30-60 minutos** por la compilación de paquetes R.

---

## Paso 4 — Configurar Variables de Entorno

```bash
# Copiar archivo de ejemplo
cp .env.example .env

# Editar con tu editor preferido
nano .env  # o vim .env / code .env
```

Configura al menos:

```env
# Requerido para E-utilities (evita bloqueos por rate limit)
NCBI_EMAIL=tu_email@ejemplo.com

# Opcional: obtener en https://www.ncbi.nlm.nih.gov/account/
NCBI_API_KEY=tu_api_key

# LLM Backend: local | anthropic | openai
LLM_BACKEND=local
```

---

## Paso 5 — Configurar LLM Local (Opcional pero Recomendado)

Para usar el Agent 7 (InsightGenerationAgent) en modo local:

```bash
# Instalar Ollama
curl -fsSL https://ollama.com/install.sh | sh

# Descargar modelo (elige según tu RAM disponible)
ollama pull llama3:8b    # 8GB VRAM/RAM — modo por defecto
ollama pull llama3:70b   # 40GB VRAM — requiere GPU RTX 4090

# Verificar que Ollama está corriendo
ollama list
```

> 💡 **Sin GPU:** llama3:8b funciona en CPU pero tarda ~10 min/query.  
> **Con GPU RTX 4090:** llama3:70b en 4-bit quantization (~2 min/query).

---

## Paso 6 — Verificar Instalación

```bash
# Verificar todas las dependencias
python -m biosignal.cli doctor

# Listar enfermedades con aliases predefinidos
python -m biosignal.cli list-diseases
```

Deberías ver todos los checks en verde ✓.

---

## Paso 7 — Primera Ejecución (Test)

```bash
# Prueba rápida: solo Agent 1 (descubrimiento, ~1 min)
python -m biosignal.cli run \
  --disease 'pancreatic cancer' \
  --max-datasets 5 \
  --start-agent 1 \
  --stop-agent 1

# Ver archivos generados
ls data/discovery/
cat data/discovery/datasets_pancreatic_cancer.json | head -50
```

---

## Paso 8 — Pipeline Completo

```bash
python -m biosignal.cli run \
  --disease 'pancreatic cancer' \
  --max-datasets 15 \
  --min-samples 10 \
  --parallel 4 \
  --output-dir data/ \
  --report-format pdf,json,csv
```

---

## Solución de Problemas de Instalación

### Error: `rpy2` no puede encontrar R

```bash
# Verificar que R está en el PATH
which R
R --version

# Configurar R_HOME si es necesario
export R_HOME=$(R RHOME)
echo "export R_HOME=$(R RHOME)" >> ~/.bashrc
```

### Error: Paquetes R de Bioconductor no instalan

```bash
# Actualizar BiocManager
R -e 'BiocManager::install(version = "3.18")'

# Instalar paquete específico manualmente
R -e 'BiocManager::install("DESeq2")'
```

### Error: `conda env create` falla

```bash
# Probar con mamba (más rápido)
conda install mamba -n base -c conda-forge
mamba env create -f environment.yml
```

### Error de SSL en descargas NCBI

```bash
# En macOS, instalar certificados de Python
/Applications/Python*/Install\ Certificates.command

# En Linux, actualizar certificados
sudo apt-get update && sudo apt-get install ca-certificates
```

---

## Verificación de Hardware para LLM Local

```bash
# Verificar GPU disponible
nvidia-smi

# Verificar VRAM disponible
nvidia-smi --query-gpu=memory.free,memory.total --format=csv
```

| GPU | VRAM | Modelo recomendado |
|-----|------|-------------------|
| RTX 4090 | 24 GB | llama3:70b (4-bit) |
| RTX 3080/3090 | 10-24 GB | llama3:8b o llama3:13b |
| Sin GPU | CPU | llama3:8b (lento) |
| — | API | claude-sonnet-4 / gpt-4o |