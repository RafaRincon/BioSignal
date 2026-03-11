/**
 * services/api.js
 * ─────────────────────────────────────────────────────────────
 * Capa de acceso al backend Express.
 * Todos los componentes consumen datos ÚNICAMENTE desde aquí.
 * Para cambiar endpoints o mocks, edita solo este archivo.
 */

const BASE = '/api'

async function request(path, options = {}) {
  const res = await fetch(`${BASE}${path}`, {
    headers: { 'Content-Type': 'application/json' },
    ...options,
  })
  if (!res.ok) {
    const err = await res.json().catch(() => ({ message: res.statusText }))
    throw new Error(err.message ?? `HTTP ${res.status}`)
  }
  return res.json()
}

// ── Reports ────────────────────────────────────────────────────────────────

/**
 * Devuelve el reporte completo de una enfermedad.
 * Backend lee: data/reports/{disease}/report_{disease}_{date}.json
 * @param {string} disease  ej: "alzheimer"
 * @returns {Promise<object>}
 */
export function fetchReport(disease) {
  return request(`/reports/${disease}`)
}

/**
 * Lista todos los reportes disponibles.
 * @returns {Promise<Array<{disease, date, path}>>}
 */
export function fetchReportList() {
  return request('/reports')
}

// ── Pipeline ───────────────────────────────────────────────────────────────

/**
 * Lanza un nuevo análisis.
 * @param {object} params  { disease, organism, min_samples, max_datasets, dea_method, data_types }
 * @returns {Promise<{ run_id: string }>}
 */
export function runPipeline(params) {
  return request('/pipeline/run', { method: 'POST', body: JSON.stringify(params) })
}

/**
 * Consulta el estado de un run activo.
 * @param {string} runId
 * @returns {Promise<{ status, current_agent, progress, log }>}
 */
export function fetchPipelineStatus(runId) {
  return request(`/pipeline/status/${runId}`)
}

// ── Datasets ───────────────────────────────────────────────────────────────

/**
 * Busca datasets en GEO/NCBI para una enfermedad con filtros.
 * @param {object} params  { disease, organism, min_samples, max_datasets, data_types, year_min, year_max }
 * @returns {Promise<Array<DatasetMetadata>>}
 */
export function searchDatasets(params) {
  return request('/datasets/search', { method: 'POST', body: JSON.stringify(params) })
}
