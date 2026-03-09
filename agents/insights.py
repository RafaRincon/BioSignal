"""
BioSignal Discovery Engine
Agent 7: InsightGenerationAgent
=================================
Responsabilidad: Usar LLM configurable para interpretar resultados biológicos
y generar hipótesis terapéuticas en lenguaje natural.

Backend LLM seleccionable via settings.yaml (llm.backend):
    - "anthropic"  → Claude API (requiere ANTHROPIC_API_KEY)
    - "openai"     → GPT-4o (requiere OPENAI_API_KEY)
    - "ollama"     → Ollama local, ej. llama3:70b (requiere Ollama corriendo)
    - "basic"      → Fallback sin LLM (sin API key)

Integra evidencia de PubMed, UniProt y ChEMBL para contextualizar los hallazgos.

Uso:
    python -m agents.insights --input data/pathways/ --meta data/meta/ --output data/insights/
"""

import json
import os
import time
from pathlib import Path
from typing import Optional

import click
import requests
import yaml
from loguru import logger

# Intentar importar SDKs opcionales
try:
    import anthropic
    ANTHROPIC_AVAILABLE = True
except ImportError:
    ANTHROPIC_AVAILABLE = False

try:
    import openai
    OPENAI_AVAILABLE = True
except ImportError:
    OPENAI_AVAILABLE = False


# Prompt del sistema para el agente biomédicco
SYSTEM_PROMPT = """You are a world-class biomedical research expert specializing in transcriptomics,
systems biology, and drug discovery. Your role is to interpret multi-dataset transcriptomic
analysis results and generate actionable biological insights.

Given a disease name, differentially expressed genes, and enriched pathways, you must:
1. Identify the most biologically relevant and reproducible signals
2. Explain the mechanistic basis of the key dysregulated pathways
3. Propose therapeutic hypotheses based on the gene expression patterns
4. Identify druggable targets among the consensus genes
5. Highlight any unexpected or particularly interesting findings

Be precise, cite biological mechanisms, and distinguish between well-established findings
and novel hypotheses. Use standard gene nomenclature (HGNC symbols).

Always structure your response in clearly labeled sections.
Output ONLY valid JSON as specified — no preamble, no markdown backticks."""


class InsightGenerationAgent:
    """
    Agente de generación de insights biológicos usando Claude API.

    Herramientas externas integradas:
        - PubMed E-utilities: busca publicaciones relevantes
        - UniProt REST API: función de proteínas e interacciones
        - ChEMBL API: compuestos bioactivos que modulan los targets

    Proceso:
        1. Cargar consensus_genes.json (Agent 5) y top_pathways_summary.json (Agent 6)
        2. Consultar PubMed, UniProt y ChEMBL para los top genes
        3. Construir contexto enriquecido
        4. Llamar a Claude API con el contexto completo
        5. Parsear y guardar respuesta estructurada

    Salida:
        data/insights/{disease}/biological_insights.json
        data/insights/{disease}/therapeutic_targets.json
        data/insights/{disease}/insight_summary.txt
    """

    # Backends soportados
    SUPPORTED_BACKENDS = ("anthropic", "openai", "ollama", "basic")

    def __init__(self, config: dict = None, settings_path: str = "settings.yaml"):
        # Cargar settings.yaml y sobreescribir con config explícita
        self.settings = self._load_settings(settings_path)
        llm_cfg = self.settings.get("llm", {})

        self.config = config or {}
        self.top_genes_for_context = self.config.get("top_genes_for_context", llm_cfg.get("top_genes_for_context", 20))
        self.max_tokens = self.config.get("max_tokens", llm_cfg.get("max_tokens", 4096))
        self.ncbi_email = self.config.get("ncbi_email", os.getenv("NCBI_EMAIL", "research@biosignal.ai"))

        # Backend LLM: settings.yaml > config explícita > fallback "basic"
        self.backend = self.config.get("backend") or llm_cfg.get("backend", "basic")
        if self.backend not in self.SUPPORTED_BACKENDS:
            logger.warning(f"[InsightAgent] Backend '{self.backend}' desconocido — usando 'basic'")
            self.backend = "basic"

        self.client = self._init_client(llm_cfg)
        logger.info(f"[InsightAgent] Backend LLM: {self.backend}")

    def _load_settings(self, settings_path: str) -> dict:
        """Carga settings.yaml si existe."""
        path = Path(settings_path)
        if not path.exists():
            # Buscar en rutas comunes del proyecto
            for candidate in ["config/settings.yaml", "../settings.yaml"]:
                if Path(candidate).exists():
                    path = Path(candidate)
                    break
            else:
                return {}
        try:
            with open(path, encoding="utf-8") as f:
                return yaml.safe_load(f) or {}
        except Exception as e:
            logger.warning(f"[InsightAgent] No se pudo cargar {path}: {e}")
            return {}

    def _init_client(self, llm_cfg: dict):
        """Inicializa el cliente LLM según el backend configurado."""
        if self.backend == "anthropic":
            api_key = os.getenv("ANTHROPIC_API_KEY", llm_cfg.get("api_key", ""))
            if not api_key:
                logger.warning("[InsightAgent] ANTHROPIC_API_KEY no configurada — fallback a 'basic'")
                self.backend = "basic"
                return None
            if not ANTHROPIC_AVAILABLE:
                logger.warning("[InsightAgent] SDK anthropic no instalado — fallback a 'basic'")
                self.backend = "basic"
                return None
            self.model = llm_cfg.get("model", "claude-sonnet-4-20250514")
            return anthropic.Anthropic(api_key=api_key)

        elif self.backend == "openai":
            api_key = os.getenv("OPENAI_API_KEY", llm_cfg.get("api_key", ""))
            if not api_key or not OPENAI_AVAILABLE:
                logger.warning("[InsightAgent] OpenAI no disponible — fallback a 'basic'")
                self.backend = "basic"
                return None
            self.model = llm_cfg.get("model", "gpt-4o")
            return openai.OpenAI(api_key=api_key)

        elif self.backend == "ollama":
            self.ollama_url = llm_cfg.get("ollama_url", "http://localhost:11434")
            self.model = llm_cfg.get("model", "llama3:70b")
            return "ollama"  # Placeholder — usa requests directamente

        else:  # "basic"
            return None

    # ------------------------------------------------------------------
    # Entrada principal
    # ------------------------------------------------------------------

    def run(self, pathways_dir: str, meta_dir: str, output_dir: str) -> dict:
        """
        Genera insights biológicos para cada enfermedad analizada.

        Args:
            pathways_dir: Directorio del Agent 6 (data/pathways/)
            meta_dir:     Directorio del Agent 5 (data/meta/)
            output_dir:   Directorio destino (data/insights/)

        Returns:
            Resumen de insights generados
        """
        pathways_path = Path(pathways_dir)
        meta_path = Path(meta_dir)
        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        summary = {"generated": [], "failed": []}
        start = time.time()

        disease_dirs = sorted([d for d in pathways_path.iterdir() if d.is_dir()])
        logger.info(f"[InsightAgent] Procesando {len(disease_dirs)} enfermedades")

        for disease_dir in disease_dirs:
            disease_name = disease_dir.name
            try:
                result = self._generate_insights(
                    disease_name=disease_name,
                    pathways_dir=disease_dir,
                    meta_dir=meta_path / disease_name,
                    out_dir=out_path / disease_name,
                )
                summary["generated"].append({"disease": disease_name, "targets": result.get("n_targets", 0)})
            except Exception as e:
                logger.error(f"[{disease_name}] Error generando insights: {e}")
                summary["failed"].append({"disease": disease_name, "error": str(e)})

        elapsed = time.time() - start
        logger.info(f"[InsightAgent] Completado en {elapsed:.1f}s")
        return summary

    # ------------------------------------------------------------------
    # Generación de insights por enfermedad
    # ------------------------------------------------------------------

    def _generate_insights(self, disease_name: str, pathways_dir: Path,
                           meta_dir: Path, out_dir: Path) -> dict:
        """Genera insights para una enfermedad."""
        out_dir.mkdir(parents=True, exist_ok=True)

        # 1. Cargar datos del pipeline
        consensus_genes = self._load_consensus_genes(meta_dir)
        top_pathways = self._load_top_pathways(pathways_dir)

        if not consensus_genes:
            raise ValueError("No se encontraron genes consenso del meta-análisis")

        # 2. Seleccionar top genes para contexto
        top_genes = consensus_genes[:self.top_genes_for_context]
        gene_symbols = [g["gene"] for g in top_genes]

        logger.info(f"[{disease_name}] Enriqueciendo contexto para {len(gene_symbols)} genes top")

        # 3. Recopilar contexto externo
        pubmed_context = self._query_pubmed(disease_name, gene_symbols[:5])
        uniprot_context = self._query_uniprot(gene_symbols[:10])
        chembl_context = self._query_chembl(gene_symbols[:10])

        # 4. Construir prompt enriquecido
        user_prompt = self._build_prompt(
            disease_name=disease_name,
            consensus_genes=top_genes,
            top_pathways=top_pathways,
            pubmed_context=pubmed_context,
            uniprot_context=uniprot_context,
            chembl_context=chembl_context,
        )

        # 5. Llamar al LLM configurado
        if self.client is None and self.backend == "basic":
            logger.warning(f"[{disease_name}] Modo basic — generando insights sin LLM")
            insights = self._generate_basic_insights(disease_name, top_genes, top_pathways)
        else:
            insights = self._call_llm(user_prompt, disease_name)

        # 6. Post-procesar y guardar
        insights["disease"] = disease_name
        insights["genes_analyzed"] = len(consensus_genes)
        insights["generated_at"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

        therapeutic_targets = insights.get("therapeutic_targets", [])

        with open(out_dir / "biological_insights.json", "w", encoding="utf-8") as f:
            json.dump(insights, f, indent=2, ensure_ascii=False)

        with open(out_dir / "therapeutic_targets.json", "w") as f:
            json.dump({"disease": disease_name, "targets": therapeutic_targets}, f, indent=2)

        # Resumen legible en texto plano
        summary_txt = self._format_summary_text(insights)
        (out_dir / "insight_summary.txt").write_text(summary_txt, encoding="utf-8")

        logger.info(f"[{disease_name}] Insights generados | Targets: {len(therapeutic_targets)}")
        return {"n_targets": len(therapeutic_targets)}

    # ------------------------------------------------------------------
    # Claude API call
    # ------------------------------------------------------------------

    def _call_llm(self, prompt: str, disease_name: str) -> dict:
        """Despacha la llamada al backend LLM configurado."""
        if self.backend == "anthropic":
            return self._call_anthropic(prompt, disease_name)
        elif self.backend == "openai":
            return self._call_openai(prompt, disease_name)
        elif self.backend == "ollama":
            return self._call_ollama(prompt, disease_name)
        else:
            return self._generate_basic_insights(disease_name, [], {})

    def _call_anthropic(self, prompt: str, disease_name: str) -> dict:
        """Llama a Anthropic Claude API."""
        try:
            message = self.client.messages.create(
                model=self.model,
                max_tokens=self.max_tokens,
                system=SYSTEM_PROMPT,
                messages=[{"role": "user", "content": prompt}],
            )
            raw_text = message.content[0].text.strip()
            raw_text = raw_text.replace("```json", "").replace("```", "").strip()
            return json.loads(raw_text)
        except json.JSONDecodeError as e:
            logger.error(f"[{disease_name}] JSON parse error (Anthropic): {e}")
            return {"error": "json_parse_error"}
        except Exception as e:
            logger.error(f"[{disease_name}] Error Anthropic API: {e}")
            raise

    def _call_openai(self, prompt: str, disease_name: str) -> dict:
        """Llama a OpenAI API (GPT-4o)."""
        try:
            response = self.client.chat.completions.create(
                model=self.model,
                max_tokens=self.max_tokens,
                messages=[
                    {"role": "system", "content": SYSTEM_PROMPT},
                    {"role": "user", "content": prompt},
                ],
                response_format={"type": "json_object"},
            )
            raw_text = response.choices[0].message.content.strip()
            return json.loads(raw_text)
        except Exception as e:
            logger.error(f"[{disease_name}] Error OpenAI API: {e}")
            raise

    def _call_ollama(self, prompt: str, disease_name: str) -> dict:
        """Llama a Ollama local via HTTP."""
        try:
            payload = {
                "model": self.model,
                "prompt": f"{SYSTEM_PROMPT}\n\n{prompt}",
                "stream": False,
                "format": "json",
            }
            r = requests.post(
                f"{self.ollama_url}/api/generate",
                json=payload,
                timeout=300,  # Modelos locales pueden ser lentos
            )
            r.raise_for_status()
            raw_text = r.json().get("response", "").strip()
            raw_text = raw_text.replace("```json", "").replace("```", "").strip()
            return json.loads(raw_text)
        except Exception as e:
            logger.error(f"[{disease_name}] Error Ollama ({self.ollama_url}): {e}")
            raise

    # ------------------------------------------------------------------
    # Construcción del prompt
    # ------------------------------------------------------------------

    def _build_prompt(self, disease_name: str, consensus_genes: list,
                      top_pathways: dict, pubmed_context: str,
                      uniprot_context: str, chembl_context: str) -> str:
        """Construye el prompt enriquecido para Claude."""

        genes_text = json.dumps(consensus_genes[:15], indent=2)
        pathways_text = json.dumps(top_pathways, indent=2)

        return f"""Analyze the following multi-dataset transcriptomic findings for {disease_name}:

## CONSENSUS DIFFERENTIALLY EXPRESSED GENES (top {min(15, len(consensus_genes))} genes, reproducible across datasets):
{genes_text}

## TOP ENRICHED BIOLOGICAL PATHWAYS:
{pathways_text}

## PUBMED CONTEXT (recent publications):
{pubmed_context}

## UNIPROT PROTEIN FUNCTIONS:
{uniprot_context}

## CHEMBL DRUGGABILITY DATA:
{chembl_context}

Please provide a comprehensive biological interpretation. Respond ONLY with a JSON object
with this exact structure (no markdown, no backticks):
{{
  "executive_summary": "2-3 sentence summary of the main biological findings",
  "key_mechanisms": [
    {{
      "mechanism": "name of biological mechanism",
      "evidence_strength": "High|Medium|Low",
      "genes_involved": ["GENE1", "GENE2"],
      "pathways_involved": ["pathway1"],
      "description": "mechanistic explanation"
    }}
  ],
  "therapeutic_targets": [
    {{
      "gene": "GENE_SYMBOL",
      "confidence": "High|Medium|Low",
      "rationale": "why this is a good target",
      "known_drugs": ["drug1", "drug2"],
      "druggability_score": "high|medium|low|unknown"
    }}
  ],
  "novel_hypotheses": [
    {{
      "hypothesis": "description of the novel finding or hypothesis",
      "genes_involved": ["GENE1"],
      "testable_prediction": "how this could be experimentally validated"
    }}
  ],
  "biomarker_candidates": [
    {{
      "gene": "GENE_SYMBOL",
      "biomarker_type": "diagnostic|prognostic|predictive",
      "rationale": "why this gene is a promising biomarker"
    }}
  ],
  "limitations": ["limitation1", "limitation2"],
  "recommended_next_steps": ["step1", "step2"]
}}"""

    # ------------------------------------------------------------------
    # Contexto externo: PubMed
    # ------------------------------------------------------------------

    def _query_pubmed(self, disease: str, genes: list, max_results: int = 5) -> str:
        """Busca publicaciones recientes en PubMed sobre la enfermedad y genes top."""
        try:
            query = f"{disease}[Title/Abstract] AND ({' OR '.join(genes[:3])}[Title/Abstract])"
            base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

            # Search
            search_url = f"{base_url}esearch.fcgi"
            params = {
                "db": "pubmed", "term": query, "retmax": max_results,
                "sort": "relevance", "email": self.ncbi_email,
                "usehistory": "y", "retmode": "json"
            }
            r = requests.get(search_url, params=params, timeout=15)
            data = r.json()
            ids = data.get("esearchresult", {}).get("idlist", [])

            if not ids:
                return "No se encontraron publicaciones recientes relevantes."

            # Fetch summaries
            summary_url = f"{base_url}esummary.fcgi"
            params = {"db": "pubmed", "id": ",".join(ids), "retmode": "json", "email": self.ncbi_email}
            r = requests.get(summary_url, params=params, timeout=15)
            summaries = r.json().get("result", {})

            lines = []
            for pmid in ids:
                if pmid in summaries:
                    s = summaries[pmid]
                    title = s.get("title", "")
                    year = s.get("pubdate", "")[:4]
                    lines.append(f"- [{year}] {title} (PMID:{pmid})")

            return "\n".join(lines) if lines else "Sin publicaciones encontradas."

        except Exception as e:
            logger.warning(f"Error consultando PubMed: {e}")
            return "PubMed no disponible."

    # ------------------------------------------------------------------
    # Contexto externo: UniProt
    # ------------------------------------------------------------------

    def _query_uniprot(self, genes: list) -> str:
        """Obtiene funciones de proteínas desde UniProt REST API."""
        results = []
        base_url = "https://rest.uniprot.org/uniprotkb/search"

        for gene in genes[:8]:  # Limitar para no saturar la API
            try:
                params = {
                    "query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
                    "fields": "gene_names,protein_name,function",
                    "format": "json",
                    "size": 1,
                }
                r = requests.get(base_url, params=params, timeout=10)
                data = r.json()

                if data.get("results"):
                    entry = data["results"][0]
                    protein_name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
                    function = ""
                    for comment in entry.get("comments", []):
                        if comment.get("commentType") == "FUNCTION":
                            texts = comment.get("texts", [])
                            if texts:
                                function = texts[0].get("value", "")[:200]
                                break
                    if protein_name:
                        results.append(f"- {gene}: {protein_name}. {function}")

                time.sleep(0.2)  # Rate limiting

            except Exception as e:
                logger.debug(f"UniProt error para {gene}: {e}")

        return "\n".join(results) if results else "Datos UniProt no disponibles."

    # ------------------------------------------------------------------
    # Contexto externo: ChEMBL
    # ------------------------------------------------------------------

    def _query_chembl(self, genes: list) -> str:
        """Busca compuestos bioactivos que modulan los genes target en ChEMBL."""
        results = []
        base_url = "https://www.ebi.ac.uk/chembl/api/data"

        for gene in genes[:5]:
            try:
                # Buscar target por nombre de gen
                r = requests.get(
                    f"{base_url}/target/search.json",
                    params={"q": gene, "limit": 1},
                    timeout=10
                )
                data = r.json()
                targets = data.get("targets", [])

                if targets:
                    target = targets[0]
                    chembl_id = target.get("target_chembl_id", "")
                    target_name = target.get("pref_name", gene)

                    # Buscar moléculas contra este target
                    r2 = requests.get(
                        f"{base_url}/activity.json",
                        params={"target_chembl_id": chembl_id, "limit": 3, "pchembl_value__gte": 6},
                        timeout=10
                    )
                    activities = r2.json().get("activities", [])
                    drugs = list({a.get("molecule_pref_name", "") for a in activities if a.get("molecule_pref_name")})[:3]

                    if drugs:
                        results.append(f"- {gene} ({target_name}): compuestos activos → {', '.join(drugs)}")
                    else:
                        results.append(f"- {gene}: target identificado en ChEMBL, sin compuestos de alta afinidad reportados")

                time.sleep(0.3)

            except Exception as e:
                logger.debug(f"ChEMBL error para {gene}: {e}")

        return "\n".join(results) if results else "Datos ChEMBL no disponibles."

    # ------------------------------------------------------------------
    # Fallback sin API
    # ------------------------------------------------------------------

    def _generate_basic_insights(self, disease_name: str, genes: list, pathways: dict) -> dict:
        """Genera insights básicos sin LLM cuando la API no está disponible."""
        up_genes = [g["gene"] for g in genes if g.get("direction") == "UP"][:5]
        down_genes = [g["gene"] for g in genes if g.get("direction") == "DOWN"][:5]

        return {
            "executive_summary": (
                f"Análisis transcriptómico de {disease_name} identificó "
                f"{len(genes)} genes consenso reproducibles entre datasets. "
                f"Genes sobreexpresados incluyen: {', '.join(up_genes)}. "
                f"Genes subexpresados: {', '.join(down_genes)}."
            ),
            "key_mechanisms": [],
            "therapeutic_targets": [
                {"gene": g["gene"], "confidence": "Medium", "rationale": "Gen diferenciado reproducible", "known_drugs": [], "druggability_score": "unknown"}
                for g in genes[:5]
            ],
            "novel_hypotheses": [],
            "biomarker_candidates": [],
            "limitations": ["Insights generados sin LLM — requiere ANTHROPIC_API_KEY"],
            "recommended_next_steps": ["Configurar ANTHROPIC_API_KEY para análisis completo"],
            "generated_by": "fallback_basic",
        }

    # ------------------------------------------------------------------
    # Helpers de carga
    # ------------------------------------------------------------------

    def _load_consensus_genes(self, meta_dir: Path) -> list:
        consensus_path = meta_dir / "consensus_genes.json"
        if not consensus_path.exists():
            return []
        with open(consensus_path, encoding="utf-8") as f:
            data = json.load(f)
        return data.get("consensus_genes", [])

    def _load_top_pathways(self, pathways_dir: Path) -> dict:
        summary_path = pathways_dir / "top_pathways_summary.json"
        if not summary_path.exists():
            return {}
        with open(summary_path) as f:
            return json.load(f)

    def _format_summary_text(self, insights: dict) -> str:
        """Formatea los insights como texto legible."""
        lines = [
            f"BioSignal Discovery Engine — Insights Biológicos",
            f"Enfermedad: {insights.get('disease', 'N/A')}",
            f"Generado: {insights.get('generated_at', '')}",
            "=" * 60,
            "",
            "RESUMEN EJECUTIVO:",
            insights.get("executive_summary", ""),
            "",
            "TARGETS TERAPÉUTICOS:",
        ]
        for t in insights.get("therapeutic_targets", []):
            drugs = ", ".join(t.get("known_drugs", [])) or "Ninguno reportado"
            lines.append(f"  - {t['gene']} (Confianza: {t.get('confidence','?')}) | Fármacos: {drugs}")

        lines += ["", "PRÓXIMOS PASOS RECOMENDADOS:"]
        for step in insights.get("recommended_next_steps", []):
            lines.append(f"  • {step}")

        return "\n".join(lines)


# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------

@click.command()
@click.option("--input", "pathways_dir", required=True, help="Directorio pathways (data/pathways/)")
@click.option("--meta", "meta_dir", required=True, help="Directorio meta-análisis (data/meta/)")
@click.option("--output", "output_dir", default="data/insights/", help="Directorio de salida")
@click.option("--top-genes", default=20, help="Número de genes top para contexto LLM")
def main(pathways_dir, meta_dir, output_dir, top_genes):
    """Agent 7: Generación de insights biológicos con Claude API."""
    config = {"top_genes_for_context": top_genes}
    agent = InsightGenerationAgent(config=config)
    summary = agent.run(pathways_dir, meta_dir, output_dir)
    click.echo(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
