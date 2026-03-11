/**
 * components/overview/BiomarkerCandidates.jsx
 */

const TYPE_BADGE = {
  diagnostic:  'badge-blue',
  prognostic:  'badge-violet',
  predictive:  'badge-amber',
}

const CONFIDENCE_COLOR = {
  High:   'var(--emerald)',
  Medium: 'var(--amber)',
  Low:    'var(--ink-4)',
}

export default function BiomarkerCandidates({ biomarkers = [] }) {
  if (!biomarkers.length) return null

  return (
    <div className="card">
      <div className="card-header">
        <div>
          <div className="card-title">Biomarker Candidates</div>
          <div className="card-sub">diagnostic · prognostic · predictive</div>
        </div>
        <span className="badge badge-blue">{biomarkers.length}</span>
      </div>

      <div className="card-body">
        {biomarkers.map((b, i) => (
          <div key={i} className="bm-item">
            <div className="bm-gene">{b.gene}</div>
            <div className="bm-body">
              <div className="bm-type-row">
                <span className={`badge ${TYPE_BADGE[b.biomarker_type] ?? 'badge-gray'}`}>
                  {b.biomarker_type}
                </span>
                <span className="bm-confidence" style={{ color: CONFIDENCE_COLOR[b.confidence] }}>
                  {b.confidence} confidence
                </span>
              </div>
              {b.rationale && (
                <p className="bm-rationale">{b.rationale}</p>
              )}
            </div>
          </div>
        ))}
      </div>
    </div>
  )
}
