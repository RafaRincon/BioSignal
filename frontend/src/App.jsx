/**
 * App.jsx
 * ─────────────────────────────────────────────────────────────
 * Root de la aplicación. Gestiona navegación entre tabs
 * y orquesta el hook de datos.
 */

import { useState } from 'react'
import Sidebar    from './components/Sidebar'
import OverviewTab from './components/overview/OverviewTab'
import { useReport } from './hooks/useReport'

// Placeholder para tabs aún no implementadas
function PlaceholderTab({ name }) {
  return (
    <div style={{
      display: 'flex', alignItems: 'center', justifyContent: 'center',
      minHeight: 300, flexDirection: 'column', gap: 12,
      color: 'var(--ink-4)', fontFamily: 'var(--font-mono)', fontSize: 13,
    }}>
      <span style={{ fontSize: 32 }}>🧬</span>
      <span>{name} — coming soon</span>
    </div>
  )
}

export default function App() {
  const [activeTab, setActiveTab] = useState('Overview')

  // TODO: En producción, obtener disease del contexto/URL params
  const disease = 'alzheimer'
  const { data: report, loading, error } = useReport(disease)

  const runInfo = report
    ? { id: `${report.disease}_v1`, date: report.analysis_date, version: report.pipeline_version }
    : null

  return (
    <div className="app">
      <Sidebar active={activeTab} onSelect={setActiveTab} runInfo={runInfo} />

      <main className="main">
        {loading && (
          <div className="loading-state">
            <div className="spinner" />
            <span>Loading report…</span>
          </div>
        )}

        {error && (
          <div className="error-state">
            ⚠ Failed to load report: {error}
          </div>
        )}

        {!loading && !error && report && (
          <>
            {activeTab === 'Overview'   && <OverviewTab report={report} />}
            {activeTab === 'Datasets'   && <PlaceholderTab name="Datasets" />}
            {activeTab === 'Genes'      && <PlaceholderTab name="Genes" />}
            {activeTab === 'Pathways'   && <PlaceholderTab name="Pathways" />}
            {activeTab === 'Targets'    && <PlaceholderTab name="Targets" />}
            {activeTab === 'Biomarkers' && <PlaceholderTab name="Biomarkers" />}
            {activeTab === 'Hypotheses' && <PlaceholderTab name="Hypotheses" />}
          </>
        )}
      </main>
    </div>
  )
}
