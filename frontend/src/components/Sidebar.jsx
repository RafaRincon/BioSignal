/**
 * components/Sidebar.jsx
 */

const NAV_ANALYSIS = [
  { label: 'Overview',  icon: <><rect x="1" y="1" width="6" height="6" rx="1.5"/><rect x="9" y="1" width="6" height="6" rx="1.5"/><rect x="1" y="9" width="6" height="6" rx="1.5"/><rect x="9" y="9" width="6" height="6" rx="1.5"/></> },
  { label: 'Datasets',  icon: <path d="M2 4h12M2 8h12M2 12h8"/> },
  { label: 'Genes',     icon: <><circle cx="8" cy="8" r="6"/><path d="M8 5v3l2 2"/></> },
  { label: 'Pathways',  icon: <path d="M2 12 L5 7 L8 9 L11 4 L14 6"/> },
]

const NAV_DISCOVERY = [
  { label: 'Targets',    icon: <><circle cx="8" cy="6" r="3"/><path d="M3 14c0-2.76 2.24-5 5-5s5 2.24 5 5"/></> },
  { label: 'Biomarkers', icon: <><path d="M8 2v4M8 10v4M2 8h4M10 8h4"/><circle cx="8" cy="8" r="2.5"/></> },
  { label: 'Hypotheses', icon: <><path d="M2 14L8 2l6 12"/><path d="M5 9h6"/></> },
]

const NAV_TOOLS = [
  { label: 'New Analysis', icon: <><circle cx="8" cy="8" r="6"/><path d="M8 5v3l2 2"/></> },
  { label: 'Export Report',icon: <path d="M8 2v8M5 7l3 3 3-3M2 13h12"/> },
]

function NavGroup({ label, items, active, onSelect }) {
  return (
    <>
      <div className="nav-section-label">{label}</div>
      {items.map(({ label: l, icon }) => (
        <div
          key={l}
          className={`nav-item ${active === l ? 'active' : ''}`}
          onClick={() => onSelect(l)}
        >
          <svg className="nav-icon" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.5">
            {icon}
          </svg>
          {l}
        </div>
      ))}
    </>
  )
}

export default function Sidebar({ active, onSelect, runInfo }) {
  return (
    <aside className="sidebar">
      <div className="logo">
        <div className="logo-row">
          <div className="logo-mark">
            <svg viewBox="0 0 16 16" fill="none">
              <circle cx="8" cy="4"  r="2.2" fill="white" opacity="0.9"/>
              <circle cx="4" cy="11" r="1.8" fill="white" opacity="0.7"/>
              <circle cx="12"cy="11" r="1.8" fill="white" opacity="0.7"/>
              <line x1="8" y1="6.2"  x2="4"  y2="9.2"  stroke="white" strokeWidth="1"   opacity="0.6"/>
              <line x1="8" y1="6.2"  x2="12" y2="9.2"  stroke="white" strokeWidth="1"   opacity="0.6"/>
              <line x1="5.8"y1="11"  x2="10.2"y2="11"  stroke="white" strokeWidth="0.8" opacity="0.4"/>
            </svg>
          </div>
          <span className="logo-name">BioSignal</span>
        </div>
        <div className="logo-tagline">Discovery Engine</div>
      </div>

      <nav className="nav">
        <NavGroup label="Analysis"  items={NAV_ANALYSIS}  active={active} onSelect={onSelect} />
        <NavGroup label="Discovery" items={NAV_DISCOVERY} active={active} onSelect={onSelect} />
        <NavGroup label="Tools"     items={NAV_TOOLS}     active={active} onSelect={onSelect} />
      </nav>

      {runInfo && (
        <div className="sidebar-bottom">
          <div className="run-pill">
            <div className="run-pill-top">
              <div className="run-status-dot" />
              <span className="run-id">{runInfo.id}</span>
            </div>
            <div className="run-meta">{runInfo.date} · pipeline v{runInfo.version}</div>
          </div>
        </div>
      )}
    </aside>
  )
}
