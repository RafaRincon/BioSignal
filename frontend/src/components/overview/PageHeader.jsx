/**
 * components/overview/PageHeader.jsx
 * Breadcrumb, título, meta-chips y botones de acción.
 */
import { fmtNum, capitalize } from '../../utils/format'

export default function PageHeader({ report }) {
  const {
    disease,
    analysis_date,
    datasets_analyzed,
    total_samples,
    pipeline_metrics,
    top_pathways,
  } = report

  const sigGenes    = pipeline_metrics?.total_significant_genes_sum ?? 0
  const sigPathways = top_pathways?.n_significant_total ?? 0
  const nTargets    = report.therapeutic_targets?.length ?? 0

  const handleExport = (format) => {
    const blob = new Blob([JSON.stringify(report, null, 2)], { type: 'application/json' })
    const url  = URL.createObjectURL(blob)
    const a    = document.createElement('a')
    a.href     = url
    a.download = `report_${disease}_${analysis_date}.${format}`
    a.click()
    URL.revokeObjectURL(url)
  }

  return (
    <div className="page-header">
      <div className="breadcrumb">
        <span>Analyses</span>
        {' / '}
        <span style={{ color: 'var(--ink-2)' }}>{capitalize(disease)}'s Disease</span>
      </div>

      <div className="page-eyebrow">Multi-dataset transcriptomics · meta-analysis</div>

      <h1 className="page-title">
        {capitalize(disease)}'s Disease <em>Discovery Report</em>
      </h1>

      <div className="page-meta">
        <div className="meta-chip"><strong>{datasets_analyzed}</strong> datasets analyzed</div>
        <div className="meta-sep" />
        <div className="meta-chip"><strong>{fmtNum(total_samples)}</strong> total samples</div>
        <div className="meta-sep" />
        <div className="meta-chip"><strong>{fmtNum(sigGenes)}</strong> consensus genes</div>
        {sigPathways > 0 && <>
          <div className="meta-sep" />
          <div className="meta-chip"><strong>{fmtNum(sigPathways)}</strong> significant pathways</div>
        </>}
        <div className="meta-sep" />
        <div className="meta-chip"><strong>{nTargets}</strong> therapeutic targets</div>
      </div>

      <div className="header-actions">
        <button className="btn" onClick={() => handleExport('json')}>↓ Export JSON</button>
        <button className="btn" onClick={() => handleExport('csv')}>↓ Export CSV</button>
        <button className="btn btn-primary">⊕ New Analysis</button>
      </div>
    </div>
  )
}
