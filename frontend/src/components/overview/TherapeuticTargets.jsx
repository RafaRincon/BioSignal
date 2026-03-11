/**
 * components/overview/TherapeuticTargets.jsx
 * Cards de targets con barra de color lateral según confianza.
 */

const CONFIDENCE_COLOR = {
  High:   'var(--emerald)',
  Medium: 'var(--amber)',
  Low:    'var(--rose)',
}

export default function TherapeuticTargets({ targets = [] }) {
  return (
    <div className="card">
      <div className="card-header">
        <div>
          <div className="card-title">Therapeutic Targets</div>
          <div className="card-sub">ChEMBL · druggability scoring</div>
        </div>
        <span className="badge badge-violet">{targets.length} identified</span>
      </div>

      <div className="card-body">
        {targets.map((t, i) => {
          const color = CONFIDENCE_COLOR[t.confidence] ?? 'var(--ink-4)'
          return (
            <div
              key={i}
              className="target-item"
              style={{ '--target-color': color }}
            >
              <div className="target-top">
                <span className="target-gene">{t.gene}</span>
                <div style={{ display: 'flex', gap: 6, alignItems: 'center' }}>
                  <span className="badge" style={{
                    background: `color-mix(in srgb, ${color} 10%, white)`,
                    color, border: `1px solid color-mix(in srgb, ${color} 25%, white)`,
                  }}>
                    {t.confidence}
                  </span>
                  <span className="badge badge-gray">
                    Druggability: {t.druggability_score}
                  </span>
                </div>
              </div>

              <p className="target-rationale">{t.rationale}</p>

              <div className="target-tags">
                {t.known_drugs?.length > 0
                  ? t.known_drugs.map(d => <span key={d} className="tag">{d}</span>)
                  : <span className="tag" style={{ color: 'var(--ink-4)' }}>No known drugs</span>
                }
              </div>
            </div>
          )
        })}
      </div>
    </div>
  )
}
