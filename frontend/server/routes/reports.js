/**
 * server/routes/reports.js
 * ─────────────────────────────────────────────────────────────
 * GET  /api/reports          → lista reportes disponibles
 * GET  /api/reports/:disease → devuelve el JSON más reciente
 *
 * Estructura esperada en disco:
 *   data/reports/{disease}/report_{disease}_{date}.json
 *   ej: data/reports/alzheimer/report_alzheimer_20260310.json
 */

import express from 'express'
import fs      from 'fs'
import path    from 'path'

const router  = express.Router()
const DATA_DIR = path.resolve('..', 'data', 'reports')

/** Devuelve el archivo más reciente de una carpeta */
function latestReport(dir) {
  const files = fs.readdirSync(dir)
    .filter(f => f.endsWith('.json'))
    .sort()                      // ISO date → lexicographic sort funciona
    .reverse()
  return files[0] ?? null
}

// GET /api/reports
router.get('/', (_req, res) => {
  if (!fs.existsSync(DATA_DIR)) return res.json([])

  const diseases = fs.readdirSync(DATA_DIR, { withFileTypes: true })
    .filter(d => d.isDirectory())
    .map(d => {
      const dir    = path.join(DATA_DIR, d.name)
      const latest = latestReport(dir)
      return latest
        ? { disease: d.name, file: latest, path: path.join(dir, latest) }
        : null
    })
    .filter(Boolean)

  res.json(diseases)
})

// GET /api/reports/:disease
router.get('/:disease', (req, res) => {
  const { disease } = req.params
  const dir = path.join(DATA_DIR, disease)

  if (!fs.existsSync(dir)) {
    return res.status(404).json({ message: `No report found for disease: ${disease}` })
  }

  const latest = latestReport(dir)
  if (!latest) {
    return res.status(404).json({ message: `Report directory exists but is empty: ${disease}` })
  }

  try {
    const raw    = fs.readFileSync(path.join(dir, latest), 'utf-8')
    const report = JSON.parse(raw)
    res.json(report)
  } catch (err) {
    res.status(500).json({ message: `Failed to parse report: ${err.message}` })
  }
})

export default router
