#!/usr/bin/env Rscript
# ============================================================
# BioSignal Discovery Engine
# Script de instalación de paquetes R
# Uso: Rscript install_r_packages.R
# ============================================================

cat("=== BioSignal Discovery Engine: Instalación de Paquetes R ===\n")
cat("R version:", R.version.string, "\n\n")

# --- Verificar versión de R ---
if (as.numeric(R.version$major) < 4) {
  stop("ERROR: Se requiere R 4.0 o superior. Version actual: ", R.version$major)
}

# --- Instalar BiocManager si no existe ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Instalando BiocManager...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install(version = "3.22", ask = FALSE)

# --- Paquetes CRAN ---
cran_packages <- c(
  "tidyverse",
  "jsonlite",
  "data.table",
  "metafor"       # Random Effects Model para meta-análisis
)

cat("\n[1/4] Instalando paquetes CRAN...\n")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Instalando %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  } else {
    cat(sprintf("  ✓ %s ya instalado\n", pkg))
  }
}

# --- Paquetes Bioconductor: Core ---
bioc_core <- c(
  "BiocGenerics",
  "S4Vectors",
  "IRanges",
  "GenomicRanges",
  "SummarizedExperiment",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "biomaRt"
)

cat("\n[2/4] Instalando paquetes Bioconductor core...\n")
for (pkg in bioc_core) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Instalando %s...\n", pkg))
    BiocManager::install(pkg, ask = FALSE, quiet = TRUE)
  } else {
    cat(sprintf("  ✓ %s ya instalado\n", pkg))
  }
}

# --- Paquetes Bioconductor: Análisis ---
bioc_analysis <- c(
  "DESeq2",          # v1.42+ - RNA-seq differential expression
  "limma",           # v3.58+ - Microarray and RNA-seq (voom)
  "edgeR",           # v4.0+  - RNA-seq differential expression
  "GEOquery",        # v2.70+ - GEO data access
  "sva",             # v3.50+ - ComBat batch correction
  "clusterProfiler", # v4.10+ - GO/KEGG enrichment
  "ReactomePA",      # Reactome pathway analysis
  "DOSE",            # Disease ontology analysis
  "enrichplot"       # Visualization for enrichment results
)

cat("\n[3/4] Instalando paquetes Bioconductor de análisis...\n")
for (pkg in bioc_analysis) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Instalando %s...\n", pkg))
    BiocManager::install(pkg, ask = FALSE, quiet = TRUE)
  } else {
    cat(sprintf("  ✓ %s ya instalado\n", pkg))
  }
}

# --- Verificación final ---
cat("\n[4/4] Verificando instalaciones...\n")
all_packages <- c(cran_packages, bioc_core, bioc_analysis)
failed <- c()

for (pkg in all_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    failed <- c(failed, pkg)
    cat(sprintf("  ✗ FALLO: %s\n", pkg))
  } else {
    version <- tryCatch(
      as.character(packageVersion(pkg)),
      error = function(e) "?"
    )
    cat(sprintf("  ✓ %s (%s)\n", pkg, version))
  }
}

# --- Reporte final ---
cat("\n=== Resumen de Instalación ===\n")
cat(sprintf("Total paquetes: %d\n", length(all_packages)))
cat(sprintf("Instalados correctamente: %d\n", length(all_packages) - length(failed)))

if (length(failed) > 0) {
  cat(sprintf("Fallos: %d\n", length(failed)))
  cat("Paquetes con error:\n")
  for (pkg in failed) cat(sprintf("  - %s\n", pkg))
  cat("\nIntenta instalar manualmente los paquetes fallidos con:\n")
  cat('BiocManager::install(c("', paste(failed, collapse = '", "'), '"))\n', sep = "")
} else {
  cat("✓ Todos los paquetes instalados correctamente!\n")
}

cat("\nInstalación completada.\n")