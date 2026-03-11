/**
 * components/overview/KeyMechanisms.jsx
 * Lista de mecanismos clave con chips de genes.
 */
import { evidencePips } from '../../utils/format'

function EvidencePips({ strength }) {
  const filled = evidencePips(strength)
  const color  = { 3: 'var(--emerald)', 2: 'var(--amber)', 1: 'var(--rose)' }[filled] ?? 'var(--ink-4)'
  return (
    <div style={{ display: 'flex', gap: 3, alignItems: 'center' }}>
      {[1, 2, 3].map(i => (
        <div key={i} style={{
          width: 8, height: 8, borderRadius: 2,
          background: i <= filled ? color : 'var(--bg-subtle)',
          border: `1px solid ${i <= filled ? color : 'var(--border)'}`,
        }} />
      ))}
      <span style={{ fontFamily: 'var(--font-mono)', fontSize: 10, color: 'var(--ink-4)', marginLeft: 4 }}>
        {strength}
      </span>
    </div>
  )
}

export default function KeyMechanisms({ mechanisms = [] }) {
  return (
    <div className="card">
      <div className="card-header">
        <div>
          <div className="card-title">Key Mechanisms</div>
          <div className="card-sub">Evidence-ranked biological processes</div>
        </div>
        <span className="badge badge-blue">{mechanisms.length}</span>
      </div>

      <div className="card-body">
        {mechanisms.map((m, i) => (
          <div key={i} className="mech-item">
            <span className="mech-num">{String(i + 1).padStart(2, '0')}</span>
            <div>
              <div style={{ display: 'flex', alignItems: 'center', gap: 10, marginBottom: 6 }}>
                <div className="mech-name">{m.mechanism}</div>
                <EvidencePips strength={m.evidence_strength} />
              </div>
              <p className="mech-desc">{m.description}</p>
              <div className="gene-chips">
                {m.genes_involved?.slice(0, 8).map(g => (
                  <span key={g} className={`gene-chip ${i % 3 === 1 ? 'violet' : i % 3 === 2 ? 'emerald' : ''}`}>
                    {g}
                  </span>
                ))}
                {(m.genes_involved?.length ?? 0) > 8 && (
                  <span className="gene-chip more">+{m.genes_involved.length - 8}</span>
                )}
              </div>
            </div>
          </div>
        ))}
      </div>
    </div>
  )
}
