/**
 * components/overview/OverviewTab.jsx
 * ─────────────────────────────────────────────────────────────
 * Orquestador de la pestaña Overview.
 * Solo importa sub-componentes y pasa los datos del reporte.
 * NO contiene lógica de fetch (eso vive en useReport.js).
 */

import PageHeader          from './PageHeader'
import PipelineStatus      from './PipelineStatus'
import SummaryHero         from './SummaryHero'
import ConsensusGenes      from './ConsensusGenes'
import TopPathways         from './TopPathways'
import TherapeuticTargets  from './TherapeuticTargets'
import KeyMechanisms       from './KeyMechanisms'
import NovelHypotheses     from './NovelHypotheses'
import BiomarkerCandidates from './BiomarkerCandidates'

export default function OverviewTab({ report }) {
  return (
    <>
      <PageHeader report={report} />

      <PipelineStatus report={report} />

      <SummaryHero report={report} />

      {/* Row 1: Genes + Targets */}
      <div className="grid-main">
        <ConsensusGenes genes={report.consensus_genes ?? []} />

        <div style={{ display: 'flex', flexDirection: 'column', gap: 20 }}>
          <TherapeuticTargets targets={report.therapeutic_targets ?? []} />
          <BiomarkerCandidates biomarkers={report.biomarker_candidates ?? []} />
        </div>
      </div>

      {/* Row 2: Pathways + Mechanisms */}
      <div className="grid-2">
        <TopPathways topPathways={report.top_pathways ?? {}} />
        <KeyMechanisms mechanisms={report.key_mechanisms ?? []} />
      </div>

      {/* Row 3: Hypotheses */}
      <NovelHypotheses hypotheses={report.novel_hypotheses ?? []} />
    </>
  )
}
