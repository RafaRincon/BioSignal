/**
 * server/index.js
 * ─────────────────────────────────────────────────────────────
 * Express backend. Corre en puerto 3001.
 * El frontend Vite corre en 5173 con proxy → 3001.
 *
 * Arrancar:   node server/index.js
 * Dev mode:   nodemon server/index.js
 */

import express  from 'express'
import cors     from 'cors'
import path     from 'path'
import { fileURLToPath } from 'url'
import reportsRouter from './routes/reports.js'

const __dirname = path.dirname(fileURLToPath(import.meta.url))
const app       = express()
const PORT      = process.env.PORT ?? 3001

app.use(cors())
app.use(express.json())

// ── API routes ────────────────────────────────────────────────
app.use('/api/reports',  reportsRouter)

// ── Health check ──────────────────────────────────────────────
app.get('/api/health', (_req, res) => res.json({ status: 'ok' }))

// ── En producción: sirve el build de Vite ────────────────────
if (process.env.NODE_ENV === 'production') {
  const distDir = path.join(__dirname, '..', 'dist')
  app.use(express.static(distDir))
  app.get('*', (_req, res) => res.sendFile(path.join(distDir, 'index.html')))
}

app.listen(PORT, () => {
  console.log(`🧬 BioSignal backend running on http://localhost:${PORT}`)
})
