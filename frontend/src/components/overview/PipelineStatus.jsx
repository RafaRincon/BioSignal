/**
 * components/overview/PipelineStatus.jsx
 * Track visual de los 8 agentes del pipeline.
 */

const AGENTS = ['Discovery', 'Download', 'Preprocess', 'DEA', 'Meta-Analysis', 'Pathways', 'Insights', 'Report']

export default function PipelineStatus({ report }) {
  const version = report.pipeline_version ?? '1.0.0'
  const date    = report.analysis_date    ?? ''

  return (
    <div className="card" style={{ marginBottom: 24 }}>
      <div className="card-header">
        <div>
          <div className="card-title">Pipeline</div>
          <div className="card-sub">8 agents · all completed successfully · {date} · v{version}</div>
        </div>
        <span className="badge badge-green">✓ Complete</span>
      </div>

      <div className="pipeline-track">
        {AGENTS.map((name, i) => (
          <div key={name} style={{ display: 'flex', alignItems: 'center', flex: 1 }}>
            <div className="step-wrap done">
              <div className="step-node">
                <div className="step-circle">{String(i + 1).padStart(2, '0')}</div>
                <div className="step-name">{name}</div>
              </div>
            </div>
            {i < AGENTS.length - 1 && <div className="step-connector done" />}
          </div>
        ))}
      </div>
    </div>
  )
}
