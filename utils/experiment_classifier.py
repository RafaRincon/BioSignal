"""
BioSignal Discovery Engine
ExperimentDesignClassifier
===========================
Clasifica muestras de un experimento GEO en 'case' / 'control' / 'exclude'
usando razonamiento semántico sobre el diseño experimental.

Se activa como fallback en Agent 2 cuando la clasificación por keywords
falla (n_case=0 o n_unclassified > umbral configurado).

Es agnóstico al dominio: funciona para enfermedades, tratamientos farmacológicos,
infecciones, condiciones ambientales, estudios de desarrollo, etc.

Backends soportados: ollama | anthropic | openai
Configuración vía bloque `llm:` en settings.yaml (mismo que Agent 7).

Uso interno (Agent 2):
    from utils.experiment_classifier import ExperimentDesignClassifier

    classifier = ExperimentDesignClassifier(llm_config)
    result = classifier.classify(
        gse_id="GSE307339",
        gse_title="Inflammation-enhanced synapse...",
        topic="alzheimer",
        samples=[{"gsm_id": "GSM001", "title": "APP-Veh1"}, ...]
    )
    # result = {
    #   "classification": {"GSM001": "case", "GSM002": "control"},
    #   "axes": ["genotype", "treatment"],
    #   "relevant_axis": "genotype",
    #   "reasoning": "...",
    #   "source": "llm_classifier",
    #   "valid": True
    # }
"""

import json
import re
from typing import Optional

from loguru import logger

# ---------------------------------------------------------------------------
# Prompt del sistema — agnóstico al dominio
# ---------------------------------------------------------------------------

_SYSTEM_PROMPT = """You are a computational biologist expert in experimental design.
Your task is to classify biological samples into case/control groups
based solely on the experimental context — regardless of biological domain.
Output ONLY valid JSON. No explanation outside the JSON."""

_USER_PROMPT_TEMPLATE = """## Study information
Title: "{gse_title}"
Research topic: "{topic}"

## Samples
{sample_list}

## Instructions

Step 1 — Identify experimental axes:
  Examine ALL sample names and extract the variables that distinguish them.
  Examples: genotype, treatment, tissue, timepoint, dose, cell_line, condition, strain.

Step 2 — Select the relevant axis:
  Determine which axis best represents the contrast of interest given the research topic.
  The relevant axis separates the biological condition under study from its baseline
  reference — regardless of domain (disease, drug, environment, infection, etc).

Step 3 — Classify each sample:
  - "case"    = condition of interest (disease, treatment, perturbation, mutant,
                exposed, infected, knockout, overexpressed, stimulated, etc.)
  - "control" = baseline reference (healthy, untreated, wild-type, vehicle,
                unexposed, naive, scramble, empty vector, etc.)
  - "exclude" = cannot be cleanly assigned (intermediate timepoint, mixed condition,
                technical replicate, secondary comparison unrelated to the topic)

## Critical rules
- Base classification ONLY on experimental design, not on topic assumptions.
- If multiple axes are equally relevant, choose the one with the clearest contrast.
- Samples belonging to a secondary comparison unrelated to the main topic → "exclude".
- Never force classification when ambiguous — use "exclude" and state why in reasoning.
- Every sample in the list must appear in the classification output.

## Response — valid JSON only, no markdown, no backticks:
{{
  "axes": ["axis1", "axis2"],
  "relevant_axis": "axis_name",
  "reasoning": "one concise sentence explaining the classification logic",
  "classification": {{
    "GSM000001": "case",
    "GSM000002": "control",
    "GSM000003": "exclude"
  }}
}}"""


# ---------------------------------------------------------------------------
# ExperimentDesignClassifier
# ---------------------------------------------------------------------------

class ExperimentDesignClassifier:
    """
    Clasifica muestras de un experimento GEO usando razonamiento LLM
    sobre el diseño experimental completo.

    Se activa únicamente cuando la clasificación por keywords de Agent 2 falla.
    """

    MAX_RETRIES = 2
    # Umbral mínimo de muestras por grupo para considerar la clasificación válida
    MIN_SAMPLES_PER_GROUP = 2

    def __init__(self, llm_config: dict):
        """
        Args:
            llm_config: Bloque `llm:` de settings.yaml
                {
                  "backend": "ollama",
                  "model": "gemma3:4b",
                  "ollama_url": "http://localhost:11434",
                  ...
                }
        """
        self.backend = llm_config.get("backend", "ollama")
        self.model = llm_config.get("model", "gemma3:4b")
        self.ollama_url = llm_config.get("ollama_url", "http://localhost:11434")
        self.max_tokens = llm_config.get("max_tokens", 1024)
        self._llm_config = llm_config

        # Inicializar cliente según backend
        self._client = self._init_client()

    # ------------------------------------------------------------------
    # API pública
    # ------------------------------------------------------------------

    def classify(
        self,
        gse_id: str,
        gse_title: str,
        topic: str,
        samples: list[dict],
    ) -> dict:
        """
        Clasifica muestras en case/control/exclude.

        Args:
            gse_id:    Identificador GEO (para logging)
            gse_title: Título del estudio
            topic:     Tema de la búsqueda (ej. "alzheimer", "antibiotic resistance")
            samples:   Lista de dicts con al menos {"gsm_id": str, "title": str}

        Returns:
            {
                "classification": {"GSM001": "case", ...},
                "axes": [...],
                "relevant_axis": "...",
                "reasoning": "...",
                "source": "llm_classifier",
                "valid": bool,
                "validation_message": str,
                "n_case": int,
                "n_control": int,
                "n_exclude": int,
            }
        """
        if not samples:
            return self._empty_result("No hay muestras para clasificar")

        if self.backend == "basic" or self._client is None:
            return self._empty_result("Backend LLM no disponible — clasificación omitida")

        prompt = self._build_prompt(gse_title, topic, samples)

        llm_result = None
        for attempt in range(1, self.MAX_RETRIES + 1):
            try:
                raw = self._call_llm(prompt, gse_id, attempt)
                llm_result = self._parse_response(raw, gse_id)
                if llm_result:
                    break
            except Exception as e:
                logger.warning(f"[{gse_id}] ExperimentClassifier intento {attempt}: {e}")

        if not llm_result:
            return self._empty_result("LLM no retornó clasificación válida tras reintentos")

        # Asegurar que todos los GSM tengan clasificación
        all_gsm = {s["gsm_id"] for s in samples}
        classified = llm_result.get("classification", {})
        missing = all_gsm - set(classified.keys())
        for gsm in missing:
            classified[gsm] = "exclude"
            logger.debug(f"[{gse_id}] {gsm} sin clasificar por LLM → 'exclude'")

        llm_result["classification"] = classified

        # Contar grupos
        n_case    = sum(1 for v in classified.values() if v == "case")
        n_control = sum(1 for v in classified.values() if v == "control")
        n_exclude = sum(1 for v in classified.values() if v == "exclude")

        # Validar
        valid, msg = self._validate(n_case, n_control)

        logger.info(
            f"[{gse_id}] ExperimentClassifier: "
            f"axis='{llm_result.get('relevant_axis')}' | "
            f"case={n_case} control={n_control} exclude={n_exclude} | "
            f"valid={valid}"
        )

        return {
            "classification": classified,
            "axes": llm_result.get("axes", []),
            "relevant_axis": llm_result.get("relevant_axis", ""),
            "reasoning": llm_result.get("reasoning", ""),
            "source": "llm_classifier",
            "valid": valid,
            "validation_message": msg,
            "n_case": n_case,
            "n_control": n_control,
            "n_exclude": n_exclude,
        }

    # ------------------------------------------------------------------
    # Construcción del prompt
    # ------------------------------------------------------------------

    def _build_prompt(self, gse_title: str, topic: str, samples: list[dict]) -> str:
        """Construye el prompt con solo título de muestra — no metadata completa."""
        sample_lines = "\n".join(
            f"  {s['gsm_id']}: {s['title']}"
            for s in samples
        )
        return _USER_PROMPT_TEMPLATE.format(
            gse_title=gse_title,
            topic=topic,
            sample_list=sample_lines,
        )

    # ------------------------------------------------------------------
    # Llamadas LLM por backend
    # ------------------------------------------------------------------

    def _call_llm(self, prompt: str, gse_id: str, attempt: int) -> str:
        """Despacha al backend correcto y retorna texto crudo."""
        # Reducir temperatura en reintentos para respuestas más conservadoras
        temperature = 0.1 if attempt == 1 else 0.0

        if self.backend == "ollama":
            return self._call_ollama(prompt, temperature)
        elif self.backend == "anthropic":
            return self._call_anthropic(prompt, temperature)
        elif self.backend == "openai":
            return self._call_openai(prompt, temperature)
        else:
            raise ValueError(f"Backend no soportado: {self.backend}")

    def _call_ollama(self, prompt: str, temperature: float) -> str:
        import requests as _req
        payload = {
            "model": self.model,
            "prompt": f"{_SYSTEM_PROMPT}\n\n{prompt}",
            "stream": False,
            "format": "json",
            "options": {"temperature": temperature},
        }
        r = _req.post(
            f"{self.ollama_url}/api/generate",
            json=payload,
            timeout=240,
        )
        r.raise_for_status()
        return r.json().get("response", "")

    def _call_anthropic(self, prompt: str, temperature: float) -> str:
        message = self._client.messages.create(
            model=self.model,
            max_tokens=self.max_tokens,
            temperature=temperature,
            system=_SYSTEM_PROMPT,
            messages=[{"role": "user", "content": prompt}],
        )
        return message.content[0].text

    def _call_openai(self, prompt: str, temperature: float) -> str:
        response = self._client.chat.completions.create(
            model=self.model,
            max_tokens=self.max_tokens,
            temperature=temperature,
            messages=[
                {"role": "system", "content": _SYSTEM_PROMPT},
                {"role": "user", "content": prompt},
            ],
            response_format={"type": "json_object"},
        )
        return response.choices[0].message.content

    # ------------------------------------------------------------------
    # Parsing y validación
    # ------------------------------------------------------------------

    def _parse_response(self, raw: str, gse_id: str) -> Optional[dict]:
        """Extrae y valida JSON de la respuesta cruda del LLM."""
        if not raw:
            return None

        # Limpiar markdown fences si los hay
        text = raw.replace("```json", "").replace("```", "").strip()

        # Intentar parse directo
        try:
            data = json.loads(text)
        except json.JSONDecodeError:
            # Intentar extraer el primer JSON objeto completo
            match = re.search(r'\{.*\}', text, re.DOTALL)
            if not match:
                logger.warning(f"[{gse_id}] ExperimentClassifier: sin JSON en respuesta")
                return None
            try:
                data = json.loads(match.group())
            except json.JSONDecodeError as e:
                logger.warning(f"[{gse_id}] ExperimentClassifier: JSON inválido: {e}")
                return None

        # Validar campos mínimos
        if "classification" not in data:
            logger.warning(f"[{gse_id}] ExperimentClassifier: falta campo 'classification'")
            return None

        # Normalizar valores a case/control/exclude
        valid_labels = {"case", "control", "exclude"}
        normalized = {}
        for gsm, label in data["classification"].items():
            label_lower = str(label).lower().strip()
            if label_lower not in valid_labels:
                label_lower = "exclude"
            normalized[gsm] = label_lower
        data["classification"] = normalized

        return data

    def _validate(self, n_case: int, n_control: int) -> tuple[bool, str]:
        """Verifica que haya suficientes muestras en ambos grupos."""
        if n_case < self.MIN_SAMPLES_PER_GROUP:
            return False, f"Insuficientes case={n_case} (mínimo={self.MIN_SAMPLES_PER_GROUP})"
        if n_control < self.MIN_SAMPLES_PER_GROUP:
            return False, f"Insuficientes control={n_control} (mínimo={self.MIN_SAMPLES_PER_GROUP})"
        return True, "OK"

    # ------------------------------------------------------------------
    # Inicialización de clientes
    # ------------------------------------------------------------------

    def _init_client(self):
        """Inicializa cliente LLM según backend."""
        if self.backend == "anthropic":
            try:
                import anthropic as _anthropic
                import os
                api_key = os.getenv("ANTHROPIC_API_KEY", "")
                if not api_key:
                    logger.warning("[ExperimentClassifier] ANTHROPIC_API_KEY no configurada → fallback basic")
                    self.backend = "basic"
                    return None
                return _anthropic.Anthropic(api_key=api_key)
            except ImportError:
                logger.warning("[ExperimentClassifier] SDK anthropic no instalado → fallback basic")
                self.backend = "basic"
                return None

        elif self.backend == "openai":
            try:
                import openai as _openai
                import os
                api_key = os.getenv("OPENAI_API_KEY", "")
                if not api_key:
                    logger.warning("[ExperimentClassifier] OPENAI_API_KEY no configurada → fallback basic")
                    self.backend = "basic"
                    return None
                return _openai.OpenAI(api_key=api_key)
            except ImportError:
                logger.warning("[ExperimentClassifier] SDK openai no instalado → fallback basic")
                self.backend = "basic"
                return None

        elif self.backend == "ollama":
            return "ollama"  # usa requests directamente

        else:
            self.backend = "basic"
            return None

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _empty_result(self, reason: str) -> dict:
        """Retorna resultado vacío con razón."""
        logger.warning(f"[ExperimentClassifier] {reason}")
        return {
            "classification": {},
            "axes": [],
            "relevant_axis": "",
            "reasoning": reason,
            "source": "llm_classifier",
            "valid": False,
            "validation_message": reason,
            "n_case": 0,
            "n_control": 0,
            "n_exclude": 0,
        }