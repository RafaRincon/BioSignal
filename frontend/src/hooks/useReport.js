/**
 * hooks/useReport.js
 * ─────────────────────────────────────────────────────────────
 * Hook que encapsula fetch + loading + error del reporte.
 * Los componentes de UI no conocen fetch ni la API directamente.
 */

import { useState, useEffect, useCallback } from 'react'
import { fetchReport } from '../services/api'

export function useReport(disease) {
  const [data,    setData]    = useState(null)
  const [loading, setLoading] = useState(true)
  const [error,   setError]   = useState(null)

  const load = useCallback(async () => {
    if (!disease) return
    setLoading(true)
    setError(null)
    try {
      const report = await fetchReport(disease)
      setData(report)
    } catch (err) {
      setError(err.message)
    } finally {
      setLoading(false)
    }
  }, [disease])

  useEffect(() => { load() }, [load])

  return { data, loading, error, reload: load }
}
