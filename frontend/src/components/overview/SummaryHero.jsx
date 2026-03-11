/**
 * components/overview/SummaryHero.jsx
 * Bloque destacado con el executive summary y stats de genes.
 */
import { fmtNum } from '../../utils/format'

/** Resalta keywords clave en el summary con <mark> */
function highlightSummary(text = '') {
  const keywords = [
    'oxidative phosphorylation', 'DNA repair', 'Myc target network',
    'unfolded protein response', 'protein translation', 'therapeutic targets',
  ]
  let result = text
  keywords.forEach((kw) => {
    result = result.replace(
      new RegExp(`(${kw})`, 'gi'),
      '<mark>$1</mark>'
    )
  })
  return result
}

export default function SummaryHero({ report }) {
  const { executive_summary, consensus_genes = [] } = report

  const up    = consensus_genes.filter(g => g.direction === 'UP').length
  const down  = consensus_genes.filter(g => g.direction === 'DOWN').length

  return (
    <div className="summary-hero">
      <div>
        <div style={{
          fontFamily: 'var(--font-mono)', fontSize: 10, letterSpacing: '0.1em',
          textTransform: 'uppercase', color: 'var(--ink-4)', marginBottom: 10,
        }}>
          Analysis Summary
        </div>
        <p
          className="summary-text"
          dangerouslySetInnerHTML={{ __html: highlightSummary(executive_summary) }}
        />
      </div>

      <div className="summary-stats">
        <div className="stat-row">
          <div className="stat-label">Consensus UP</div>
          <div className="stat-value" style={{ color: 'var(--emerald)' }}>{fmtNum(up)}</div>
          <div className="stat-sub">upregulated genes</div>
        </div>
        <div className="stat-row">
          <div className="stat-label">Consensus DOWN</div>
          <div className="stat-value" style={{ color: 'var(--rose)' }}>{fmtNum(down)}</div>
          <div className="stat-sub">downregulated genes</div>
        </div>
      </div>
    </div>
  )
}
