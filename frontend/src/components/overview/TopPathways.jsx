/**
 * components/overview/TopPathways.jsx
 * Lista de pathways top por base de datos con barras visuales.
 */
import { useState } from 'react'
import { cleanPathwayName, dbColor } from '../../utils/format'

const DB_TABS = [
  { key: 'kegg',      label: 'KEGG'     },
  { key: 'go_bp',     label: 'GO:BP'    },
  { key: 'reactome',  label: 'Reactome' },
  { key: 'hallmarks', label: 'Hallmarks'},
]

export default function TopPathways({ topPathways = {} }) {
  const byDb = topPathways.top_pathways_by_db ?? {}
  const [activeDb, setActiveDb] = useState('kegg')

  const pathways = byDb[activeDb] ?? []
  const maxGenes = Math.max(...pathways.map(p => p.genes?.length ?? 0), 1)
  const colors   = dbColor(activeDb)

  return (
    <div className="card">
      <div className="card-header">
        <div>
          <div className="card-title">Top Enriched Pathways</div>
          <div className="card-sub">{topPathways.genes_analyzed?.toLocaleString()} genes analyzed · {topPathways.n_significant_total?.toLocaleString()} significant</div>
        </div>
      </div>

      {/* DB Tabs */}
      <div style={{ display: 'flex', borderBottom: '1px solid var(--border)', padding: '0 20px' }}>
        {DB_TABS.map(({ key, label }) => {
          const c = dbColor(key)
          return (
            <button
              key={key}
              onClick={() => setActiveDb(key)}
              style={{
                background: 'none', border: 'none', cursor: 'pointer',
                padding: '8px 12px', fontSize: 12, fontWeight: 500,
                color: activeDb === key ? c.text : 'var(--ink-4)',
                borderBottom: `2px solid ${activeDb === key ? c.text : 'transparent'}`,
                transition: 'all 0.15s',
                fontFamily: 'var(--font-mono)',
              }}
            >
              {label}
            </button>
          )
        })}
      </div>

      <div className="card-body">
        {pathways.slice(0, 8).map((p, i) => {
          const genesCount = p.genes?.length ?? 0
          const width = (genesCount / maxGenes) * 100
          return (
            <div key={i} className="pathway-item">
              <span className="pathway-rank">{i + 1}</span>
              <span className="pathway-name">{cleanPathwayName(p.pathway)}</span>
              <div className="pathway-bar-wrap">
                <div className="pathway-bar">
                  <div
                    className="pathway-bar-fill"
                    style={{ width: `${width}%`, background: colors.text }}
                  />
                </div>
                <span style={{ fontFamily: 'var(--font-mono)', fontSize: 10, color: 'var(--ink-4)', minWidth: 22, textAlign: 'right' }}>
                  {genesCount}
                </span>
              </div>
            </div>
          )
        })}
      </div>
    </div>
  )
}
