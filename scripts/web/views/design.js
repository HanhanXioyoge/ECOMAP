// scripts/web/views/design.js — Design track (Phase 6 documentation).
//
// Phase 6 replaces the prior four-step wizard with hand-authored doc sections
// (spec §8.4). The only interactive surface is the Workflow "Run in MATLAB →"
// button, which dispatches a `matlab:run` CustomEvent on the root element.
// The downstream wizard wiring lives on a feature branch and is intentionally
// absent here.

import {
  mountTrackShell,
  codeBlock,
  docHeader,
  docSection,
} from './shared.js';

function okoPlusPanel(rootEl) {
  const panel = document.createElement('section');
  panel.id = 'oko-plus';
  panel.className = 'doc-section oko-plus-workflow';
  panel.innerHTML = `
    <h2>OKO+</h2>
    <p>Build cross-species UniKP/CatPred kcat intervals, inspect coverage, then run the OKO+ MILP.</p>
    <div class="form-grid">
      <label>ecModel ID<input class="input" data-oko="model" placeholder="uploaded modelId" /></label>
      <label>Biomass reaction<input class="input" data-oko="biomass" placeholder="BIOMASS" /></label>
      <label>Target reaction<input class="input" data-oko="target" placeholder="EX_product_e" /></label>
      <label>Predictor<select class="input" data-oko="predictor"><option>UniKP</option><option>CatPred</option></select></label>
      <label>Max homologs<input class="input" data-oko="max" type="number" min="2" max="500" value="100" /></label>
    </div>
    <div class="button-row">
      <button class="button button-primary" data-oko-action="build">1. Build intervals</button>
      <button class="button" data-oko-action="run" disabled>2. Run OKO+</button>
    </div>
    <p data-oko="status">Ready.</p>
    <pre class="code-block" data-oko="qc" hidden></pre>`;

  let intervalJobId = '';
  const value = (name) => panel.querySelector(`[data-oko="${name}"]`)?.value.trim() || '';
  const status = panel.querySelector('[data-oko="status"]');
  const qc = panel.querySelector('[data-oko="qc"]');
  const runButton = panel.querySelector('[data-oko-action="run"]');

  async function submit(algo, params) {
    const response = await fetch('/api/jobs', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        modelId: value('model'), biomass: value('biomass'), target: value('target'),
        algo, params,
      }),
    });
    if (!response.ok) throw new Error(await response.text());
    return (await response.json()).jobId;
  }

  function awaitJob(jid) {
    return new Promise((resolve, reject) => {
      const stream = new EventSource(`/api/jobs/${jid}/events`);
      stream.onmessage = async (event) => {
        const message = JSON.parse(event.data);
        if (message.type === 'log') status.textContent = message.line;
        if (message.type === 'error') { stream.close(); reject(new Error(message.message)); }
        if (message.type === 'done') {
          stream.close();
          const response = await fetch(`/api/jobs/${jid}/result`);
          if (!response.ok) reject(new Error(await response.text()));
          else resolve(await response.json());
        }
      };
      stream.onerror = () => { stream.close(); reject(new Error('Job event stream disconnected.')); };
    });
  }

  panel.querySelector('[data-oko-action="build"]').addEventListener('click', async () => {
    try {
      if (!value('model')) throw new Error('ecModel ID is required.');
      status.textContent = 'Submitting interval build…';
      intervalJobId = await submit('okoplus-build', {
        predictors: [value('predictor')], max_homologs: Number(value('max') || 100),
      });
      const result = await awaitJob(intervalJobId);
      status.textContent = `Intervals ready (job ${intervalJobId}).`;
      qc.hidden = false;
      qc.textContent = JSON.stringify({ runDir: result.run_dir, qc: result.qc }, null, 2);
      runButton.disabled = false;
    } catch (error) { status.textContent = `Error: ${error.message}`; }
  });

  runButton.addEventListener('click', async () => {
    try {
      status.textContent = 'Submitting OKO+ MILP…';
      const jid = await submit('okoplus', {
        intervalJobId, predictor: value('predictor'), profile: 'auto',
      });
      const result = await awaitJob(jid);
      status.textContent = `OKO+ completed (job ${jid}).`;
      qc.hidden = false;
      qc.textContent = JSON.stringify(result, null, 2);
    } catch (error) { status.textContent = `Error: ${error.message}`; }
  });

  return panel;
}

// -----------------------------------------------------------------------------
// View
// -----------------------------------------------------------------------------

export function mountDesign(rootEl) {
  if (!rootEl) return () => {};
  const cleanups = [];

  const content = document.createElement('div');

  content.appendChild(docHeader(
    'Design',
    'Find intervention targets with FSEOF, OptKnock and OptForce.',
  ));

  content.appendChild(okoPlusPanel(rootEl));

  // ----- Overview ----------------------------------------------------------
  content.appendChild(docSection({
    anchor: 'overview',
    title: 'Overview',
    paragraphs: [
      'Design turns a calibrated ecGEM into a ranked list of strain interventions. Three algorithms ship: FSEOF (Flux Scanning based on Enforced Objective Flux), OptKnock (bilevel knock-out), and OptForce (mandatory up-/down-regulation).',
      'All three produce a candidate table with reaction-level interventions and a growth-curve overlay that compares reference vs mutant.',
    ],
  }));

  // ----- Getting started ----------------------------------------------------
  content.appendChild(docSection({
    anchor: 'getting-started',
    title: 'Getting started',
    paragraphs: [
      'Paste the modelId from Reconstruction or Calibration and pick a target reaction (the exchange reaction of the metabolite you want to overproduce). Pick an algorithm based on the intervention type you want.',
      'OptKnock deletes reactions; OptForce perturbs their bounds; FSEOF ranks candidates without committing to deletions.',
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
      'Specify the target exchange reaction, the minimum biomass yield, and (for OptKnock) the number of deletions. The bridge enqueues a bilevel job; OptKnock can take a few minutes on yeast-scale models.',
      'Results stream back as a candidate table (reaction, fold-change, biomass) plus a side-by-side growth curve.',
    ],
    code: {
      lang: 'matlab',
      body: [
        'targets = optKnock(ecGEM, ...',
        '    \'targetRxn\', \'EX_etoh_e\', ...',
        '    \'numKnockouts\', 5, ...',
        '    \'solver\', \'gurobi\');',
      ].join('\n'),
    },
    runButton: {
      label: 'Run in MATLAB →',
      track: 'design',
      onClick: runBtn('workflow', 'design'),
    },
  }));

  // ----- Parameters ---------------------------------------------------------
  content.appendChild(docSection({
    anchor: 'parameters',
    title: 'Parameters',
    paragraphs: [
      'target — exchange reaction id (e.g. "EX_etoh_e"). biomass — minimum biomass yield (default 0.05). numKnockouts — OptKnock deletions (default 5, max 10 for tractability).',
      'iterations / coefficient — FSEOF only. coefficient is the enforced-objective multiplier (default 0.9). solver — gurobi | cplex | glpk (auto-detected).',
    ],
  }));

  // ----- API ----------------------------------------------------------------
  content.appendChild(docSection({
    anchor: 'api',
    title: 'API',
    paragraphs: [
      'POST /api/jobs with algo=optknock / optforce / fseof / ecfseof. The job result is a candidate list (reaction, fold-change, biomass) and a growth-curve payload.',
    ],
    code: {
      lang: 'matlab',
      body: [
        'jid = design.submit(ecGEM, struct( ...',
        '    \'algo\', \'optknock\', ...',
        '    \'targetRxn\', \'EX_etoh_e\', ...',
        '    \'numKnockouts\', 5));',
        'candidates = design.await(jid);',
      ].join('\n'),
    },
  }));

  // ----- Examples -----------------------------------------------------------
  content.appendChild(docSection({
    anchor: 'examples',
    title: 'Examples',
    paragraphs: [
      'A worked ethanol-overproduction case on the yeast consensus model ships as scripts/examples/ethanol_optknock.m; expected top-3 interventions match the published strain.',
      'A lysine-overproduction FSEOF scan on the E. coli core model is in scripts/examples/lysine_fseof.m with a comparison bar chart.',
    ],
  }));

  // ----- Extras -------------------------------------------------------------
  content.appendChild(docSection({
    anchor: 'optknock',
    title: 'OptKnock',
    paragraphs: [
      'OptKnock solves a bilevel MILP: the outer problem maximises target flux, the inner problem minimises the number of reaction deletions subject to biomass ≥ minBiomass. Candidates are returned ranked by target flux at the optimum.',
      'The runner uses gurobi when available and falls back to cplex / glpk otherwise. Number of deletions is bounded by 10 for tractability.',
    ],
  }));

  content.appendChild(docSection({
    anchor: 'optforce',
    title: 'OptForce',
    paragraphs: [
      'OptForce decomposes the problem: identify MUST sets (reactions that must change), then solve two LPs to find MUST-UP and MUST-DOWN sets that drive target flux.',
      'The output is a compact regulatory strategy rather than a set of deletions, which makes OptForce the natural choice when genetic interventions are expensive.',
    ],
  }));

  content.appendChild(docSection({
    anchor: 'fseof',
    title: 'FSEOF',
    paragraphs: [
      'Flux Scanning based on Enforced Objective Flux pushes biomass from 100% down to 50% in 10 steps and records each reaction\'s flux trajectory. Reactions whose flux rises monotonically are over-expression (OE) candidates; reactions whose flux falls are knock-down (KD) candidates.',
      'FSEOF is fast (no MILP), scales to genome-scale models, and is the recommended starting point before committing to OptKnock.',
    ],
  }));

  // ----- Mount shell --------------------------------------------------------
  const shellCleanup = mountTrackShell(rootEl, {
    track: 'design',
    title: 'Design',
    sections: [
      { id: 'oko-plus', label: 'OKO+' },
      { id: 'overview',    label: 'Overview' },
      { id: 'getting-started', label: 'Getting started' },
      { id: 'workflow',    label: 'Workflow' },
      { id: 'parameters',  label: 'Parameters' },
      { id: 'api',         label: 'API' },
      { id: 'examples',    label: 'Examples' },
      { id: 'optknock',    label: 'OptKnock', isExtra: true },
      { id: 'optforce',    label: 'OptForce', isExtra: true },
      { id: 'fseof',       label: 'FSEOF',    isExtra: true },
    ],
    content,
  });

  return function teardown() {
    cleanups.forEach((fn) => { try { fn(); } catch (_) {} });
    if (shellCleanup) { try { shellCleanup(); } catch (_) {} }
    if (rootEl) rootEl.innerHTML = '';
  };
}
