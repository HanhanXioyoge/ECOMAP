// scripts/web/views/analysis.js — Analysis track (Phase 6 documentation).
//
// Phase 6 replaces the prior three-card task UI with hand-authored doc
// sections (spec §8.4). The only interactive surface is the Workflow "Run
// in MATLAB →" button, which dispatches a `matlab:run` CustomEvent on the
// root element. The downstream task-card wiring lives on a feature branch
// and is intentionally absent here.

import {
  mountTrackShell,
  codeBlock,
  docHeader,
  docSection,
} from './shared.js';

// -----------------------------------------------------------------------------
// View
// -----------------------------------------------------------------------------

export function mountAnalysis(rootEl) {
  if (!rootEl) return () => {};
  const cleanups = [];

  const content = document.createElement('div');

  content.appendChild(docHeader(
    'Analysis',
    'Probe flux ranges, essential genes and protein allocation.',
  ));

  // ----- Overview ----------------------------------------------------------
  content.appendChild(docSection({
    anchor: 'overview',
    title: 'Overview',
    paragraphs: [
      'Analysis runs read-only queries against a calibrated ecGEM: ecFVA computes flux ranges at a given biomass fraction; single-gene knockout scans identify essential genes; protein-usage reports break down which enzymes consume the proteome.',
      'Every task is idempotent — running the same Analysis job twice on the same ecGEM returns the same numbers — which makes it safe to drive from notebooks.',
    ],
  }));

  // ----- Getting started ----------------------------------------------------
  content.appendChild(docSection({
    anchor: 'getting-started',
    title: 'Getting started',
    paragraphs: [
      'Paste the modelId returned by Reconstruction or Calibration. For FVA, set the biomass fraction (default 0.9). For knockouts, paste a comma-separated gene list or leave blank to scan all genes.',
      'No MATLAB install required — the bridge proxies to a worker pool and streams progress via SSE. Results render inline as soon as the job completes.',
    ],
  }));

  // ----- Workflow ----------------------------------------------------------
  const runBtn = (anchor, track) => () => {
    rootEl.dispatchEvent(new CustomEvent('matlab:run', {
      detail: { track, anchor },
      bubbles: true,
    }));
  };

  content.appendChild(docSection({
    anchor: 'workflow',
    title: 'Workflow',
    paragraphs: [
      'Pick a task — FVA, knockout, or protein usage — point at a calibrated ecGEM, run. Results render as either an FVA bar chart, a per-gene growth table, or a per-reaction allocation table.',
      'Three independent jobs run in parallel under the hood; the UI shows three live status dots, one per task.',
    ],
    code: {
      lang: 'matlab',
      body: [
        'ranges = analysis.fva(ecGEM, ...',
        '    \'fraction\', 0.9);',
        'ko      = analysis.ko(ecGEM, ...',
        '    \'genes\', {{\'gapA\',\'pfkA\',\'aceE\'}});',
        'usage   = analysis.protein(ecGEM);',
      ].join('\n'),
    },
    runButton: {
      label: 'Run in MATLAB →',
      track: 'analysis',
      onClick: runBtn('workflow', 'analysis'),
    },
  }));

  // ----- Parameters ---------------------------------------------------------
  content.appendChild(docSection({
    anchor: 'parameters',
    title: 'Parameters',
    paragraphs: [
      'fraction — biomass fraction for FVA (default 0.9, range 0–1). Lower values surface more reaction variability but cost more solver time.',
      'geneList — comma-separated genes for knockout scan; empty list scans every gene. conditions — exchange-rxns / media override for protein-usage (e.g. "glucose=10,aerobic").',
    ],
  }));

  // ----- API ----------------------------------------------------------------
  content.appendChild(docSection({
    anchor: 'api',
    title: 'API',
    paragraphs: [
      'POST /api/jobs with algo=fva / ko / prot. All three return a jobId; the Analysis track watches it and renders the corresponding chart or table inline.',
    ],
    code: {
      lang: 'matlab',
      body: [
        'jid = analysis.submit(ecGEM, struct( ...',
        '    \'task\', \'fva\', ...',
        '    \'fraction\', 0.9));',
        'ranges = analysis.await(jid);',
      ].join('\n'),
    },
  }));

  // ----- Examples -----------------------------------------------------------
  content.appendChild(docSection({
    anchor: 'examples',
    title: 'Examples',
    paragraphs: [
      'Run ecFVA at 0.5 / 0.9 / 0.99 biomass fractions and compare the growth-coupled reaction ranges; the script ships as scripts/examples/ecfva_compare.m.',
      'A worked E. coli essential-gene screen with the calibrated yeast model is in scripts/examples/ecoli_ko.m; expected F-statistics match the published paper.',
    ],
  }));

  // ----- Extras -------------------------------------------------------------
  content.appendChild(docSection({
    anchor: 'fva',
    title: 'FVA',
    paragraphs: [
      'Flux Variability Analysis sweeps each reaction independently between its min and max feasible flux while keeping biomass ≥ fraction × v_biomass_max. The result exposes which reactions are growth-coupled and which are fully variable.',
      'ecFVA extends the standard formulation with the enzyme-pool constraint so the variability respects proteome limits.',
    ],
  }));

  content.appendChild(docSection({
    anchor: 'sensitivity',
    title: 'Sensitivity',
    paragraphs: [
      'Local sensitivity scans each parameter (kcat, f, σ) ± 10% and reports the resulting biomass change. Output is a ranked list of which parameters matter most for the current ecGEM.',
      'Global sensitivity (Sobol indices) is available via the Calibration track\'s Bayesian runner.',
    ],
  }));

  // ----- Mount shell --------------------------------------------------------
  const shellCleanup = mountTrackShell(rootEl, {
    track: 'analysis',
    title: 'Analysis',
    sections: [
      { id: 'overview',    label: 'Overview' },
      { id: 'getting-started', label: 'Getting started' },
      { id: 'workflow',    label: 'Workflow' },
      { id: 'parameters',  label: 'Parameters' },
      { id: 'api',         label: 'API' },
      { id: 'examples',    label: 'Examples' },
      { id: 'fva',         label: 'FVA',         isExtra: true },
      { id: 'sensitivity', label: 'Sensitivity', isExtra: true },
    ],
    content,
  });

  return function teardown() {
    cleanups.forEach((fn) => { try { fn(); } catch (_) {} });
    if (shellCleanup) { try { shellCleanup(); } catch (_) {} }
    if (rootEl) rootEl.innerHTML = '';
  };
}
