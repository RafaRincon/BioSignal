/**
 * components/overview/NovelHypotheses.jsx
 * Bloque de hipótesis generadas por LLM con predicción testeable.
 */

export default function NovelHypotheses({ hypotheses = [] }) {
  if (!hypotheses.length) return null

  return (
    <div className="card">
      <div className="card-header">
        <div>
          <div className="card-title">Novel Hypotheses</div>
          <div className="card-sub">LLM-generated · require experimental validation</div>
        </div>
        <span className="badge badge-amber">{hypotheses.length} generated</span>
      </div>

      <div className="card-body" style={{ display: 'flex', flexDirection: 'column', gap: 14 }}>
        {hypotheses.map((h, i) => (
          <div key={i} className="hypothesis-block">
            <div className="hypothesis-eyebrow">Hypothesis {i + 1}</div>
            <p className="hypothesis-text">{h.hypothesis}</p>

            {h.testable_prediction && (
              <div className="hypothesis-predict">
                <div className="predict-label">Testable Prediction</div>
                {h.testable_prediction}
              </div>
            )}

            {h.genes_involved?.length > 0 && (
              <div className="gene-chips" style={{ marginTop: 10 }}>
                {h.genes_involved.map(g => (
                  <span key={g} className="gene-chip">{g}</span>
                ))}
              </div>
            )}
          </div>
        ))}
      </div>
    </div>
  )
}
