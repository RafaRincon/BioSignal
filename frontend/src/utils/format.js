/**
 * utils/format.js
 * Funciones puras de formato. Sin dependencias de React.
 */

/** Capitaliza primera letra */
export const capitalize = (s = '') => s.charAt(0).toUpperCase() + s.slice(1)

/** Formatea número con comas */
export const fmtNum = (n) => (n ?? 0).toLocaleString()

/** log2FC → color CSS según dirección */
export const fc2color = (fc) =>
  fc > 0.1 ? 'var(--emerald)' : fc < -0.1 ? 'var(--rose)' : 'var(--ink-4)'

/** Direction string → clase CSS */
export const dirClass = (d) =>
  ({ UP: 'badge-green', DOWN: 'badge-rose', MIXED: 'badge-violet' }[d] ?? 'badge-gray')

/** Evidence strength → número de pips */
export const evidencePips = (s) => ({ High: 3, Medium: 2, Low: 1 }[s] ?? 0)

/** Trunca texto largo */
export const truncate = (s = '', max = 60) =>
  s.length > max ? s.slice(0, max) + '…' : s

/** Porcentaje sobre un máximo */
export const pct = (val, max) => (max > 0 ? Math.min((val / max) * 100, 100) : 0)

/** Parsea nombre de pathway — quita IDs como (GO:0019646) o R-HSA-xxx */
export const cleanPathwayName = (name = '') =>
  name.replace(/\s*\(GO:\d+\)/g, '').replace(/\s*R-HSA-\d+/g, '').trim()

/** Color por base de datos de pathways */
export const dbColor = (db = '') =>
  ({
    kegg:      { bg: 'var(--blue-bg)',    text: 'var(--blue)',    border: '#c7d9fb' },
    go_bp:     { bg: 'var(--violet-bg)',  text: 'var(--violet)',  border: '#ddd6fe' },
    reactome:  { bg: 'var(--emerald-bg)', text: 'var(--emerald)', border: '#a7f3d0' },
    hallmarks: { bg: 'var(--amber-bg)',   text: 'var(--amber)',   border: '#fde68a' },
  })[db.toLowerCase().replace(/[^a-z]/g, '_')] ?? { bg: 'var(--bg-subtle)', text: 'var(--ink-3)', border: 'var(--border)' }
