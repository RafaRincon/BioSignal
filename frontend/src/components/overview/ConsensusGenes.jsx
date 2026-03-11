/**
 * components/overview/ConsensusGenes.jsx
 * Tabla de genes con barras de log2FC y dots de consistencia.
 */
import { dirClass, pct, fmtNum } from '../../utils/format'

function MiniBar({ value }) {
  const max   = 2
  const width = pct(Math.abs(value), max)
  const isUp  = value >= 0
  return (
    <div style={{ display: 'flex', alignItems: 'center', gap: 6 }}>
      <div style={{ width: 60, height: 5, background: 'var(--bg-subtle)', borderRadius: 3, overflow: 'hidden' }}>
        <div style={{
          width: `${width}%`, height: '100%', borderRadius: 3,
          background: isUp
            ? 'linear-gradient(90deg,#10b981,#34d399)'
            : 'linear-gradient(90deg,#f43f5e,#fb7185)',
        }} />
      </div>
      <span style={{
        fontFamily: 'var(--font-mono)', fontSize: 10,
        color: isUp ? 'var(--emerald)' : 'var(--rose)',
        minWidth: 46,
      }}>
        {value > 0 ? '+' : ''}{value.toFixed(3)}
      </span>
    </div>
  )
}

function ConsistencyDots({ value }) {
  const filled = Math.round(value * 4)
  return (
    <div style={{ display: 'flex', gap: 3 }}>
      {[0, 1, 2, 3].map(i => (
        <div key={i} style={{
          width: 7, height: 7, borderRadius: '50%',
          background: i < filled ? 'var(--blue)' : 'var(--bg-subtle)',
          border: '1px solid var(--border)',
        }} />
      ))}
    </div>
  )
}

export default function ConsensusGenes({ genes = [] }) {
  const up    = genes.filter(g => g.direction === 'UP').length
  const down  = genes.filter(g => g.direction === 'DOWN').length
  const mixed = genes.filter(g => g.direction === 'MIXED').length

  return (
    <div className="card">
      <div className="card-header">
        <div>
          <div className="card-title">Consensus Genes</div>
          <div className="card-sub">
            <span style={{ color: 'var(--emerald)', fontWeight: 600 }}>{up} UP</span>
            {' · '}
            <span style={{ color: 'var(--rose)', fontWeight: 600 }}>{down} DOWN</span>
            {' · '}
            <span style={{ color: 'var(--violet)', fontWeight: 600 }}>{mixed} MIXED</span>
          </div>
        </div>
        <span className="badge badge-blue">{genes.length} total</span>
      </div>

      <table>
        <thead>
          <tr>
            <th>Gene</th>
            <th>Direction</th>
            <th>log₂FC</th>
            <th>Consistency</th>
            <th>Datasets</th>
          </tr>
        </thead>
        <tbody>
          {genes.map((g) => (
            <tr key={g.gene}>
              <td>
                <span style={{ fontFamily: 'var(--font-mono)', fontSize: 12, fontWeight: 600, color: 'var(--ink)' }}>
                  {g.gene}
                </span>
              </td>
              <td>
                <span className={`badge ${dirClass(g.direction)}`}>{g.direction}</span>
              </td>
              <td><MiniBar value={g.mean_log2fc} /></td>
              <td><ConsistencyDots value={g.direction_consistency} /></td>
              <td>
                <span className="num">{g.n_datasets_present}/{g.n_datasets_total}</span>
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  )
}
