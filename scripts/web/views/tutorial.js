import { deepLearningKcatHtml } from './shared.js';

const SECTIONS = [
  { id: 'start', title: 'Start' },
  { id: 'reconstruction', title: 'Reconstruction' },
  { id: 'calibration', title: 'Calibration' },
  { id: 'data', title: 'Data matrix' },
  { id: 'next', title: 'Analysis and Design' },
];

function sidebarLinks(prefix) {
  return SECTIONS.map((section, idx) => `
    <li><a href="#${prefix}-${section.id}" data-section="${section.id}" ${idx === 0 ? 'class="is-active"' : ''}>${section.title}</a></li>
  `).join('');
}

function bindSidebar(main, prefix) {
  const links = main.querySelectorAll('.page__sidebar-list a');
  links.forEach((link) => {
    link.addEventListener('click', (event) => {
      event.preventDefault();
      const target = main.querySelector(`#${prefix}-${link.dataset.section}`);
      if (target) target.scrollIntoView({ behavior: 'smooth', block: 'start' });
      links.forEach((item) => item.classList.toggle('is-active', item === link));
    });
  });
  if (typeof IntersectionObserver === 'undefined') return;
  const sections = SECTIONS.map((section) => main.querySelector(`#${prefix}-${section.id}`));
  const observer = new IntersectionObserver((entries) => {
    entries.forEach((entry) => {
      if (!entry.isIntersecting) return;
      const id = entry.target.id.replace(`${prefix}-`, '');
      links.forEach((item) => item.classList.toggle('is-active', item.dataset.section === id));
    });
  }, { rootMargin: '-30% 0px -60% 0px', threshold: 0 });
  sections.forEach((section) => section && observer.observe(section));
}

export function mountTutorial(main) {
  main.innerHTML = `
    <section class="page page--tutorial">
      <aside class="page__sidebar">
        <h2 class="page__sidebar-title">Tutorial</h2>
        <ul class="page__sidebar-list" role="list">${sidebarLinks('tutorial')}</ul>
      </aside>
      <main class="page__main">
        <article id="tutorial-start" class="page__section">
          <h1>Create a project first</h1>
          <p>Start the web UI with <code>ecomapWeb('start')</code>, open Projects, then create a project before entering the four modules. Upload a GEM or ecGEM, confirm the model type, and fill the shared project parameters once.</p>
          <p>Each module can still be run on its own inside the project. The module loads the project parameter set as defaults, then saves any page edits as a run-scoped snapshot instead of forcing the shared parameter manager to reinitialize.</p>
        </article>

        <article id="tutorial-reconstruction" class="page__section">
          <h1>Reconstruction tutorial</h1>
          <ol class="workflow-steps">
            <li><strong>Create a GEM project.</strong><p>Use a COBRA-compatible <code>.mat</code>, SBML, JSON, or YAML file and select GEM as the model type. The upload creates the model id used by Reconstruction.</p></li>
            <li><strong>Configure parameters.</strong><p>Choose organism preset, proteome fraction <code>f</code>, <code>sigma</code>, and one of the basic, isozyme, or integrated S-matrix topologies.</p></li>
            <li><strong>Convert and annotate.</strong><p>Generate the ecGEM, then attach protein-complex, EC-number, and metabolite information.</p></li>
            <li><strong>Predict and merge kcat.</strong><p>Run ${deepLearningKcatHtml()} models, compare predicted and database values, and merge database, prediction, and median sources.</p></li>
            <li><strong>Check growth.</strong><p>Set carbon source and biomass reaction, then export the reconstructed ecGEM for Calibration.</p></li>
          </ol>
          <figure class="figure-placeholder">Figure placeholder: Reconstruction tutorial screenshots and paper reconstruction diagram.</figure>
        </article>

        <article id="tutorial-calibration" class="page__section">
          <h1>Calibration tutorial</h1>
          <ol class="workflow-steps">
            <li><strong>Open an ecGEM-ready project.</strong><p>Calibration can start from Reconstruction output or from a project created directly with an uploaded ecGEM.</p></li>
            <li><strong>Select methods.</strong><p>The UI always runs selected methods in the order Sluice, KcatRepo, Sensitivity, Bayesian, PRESTO, and GAUKS.</p></li>
            <li><strong>Add method data.</strong><p>Sensitivity requires glucose uptake and matched growth-rate data. PRESTO requires proteomics, growth, and total protein files. GAUKS requires unconstrained maximum growth data.</p></li>
            <li><strong>Set parameters.</strong><p>Configure exchange reactions, target growth, Sensitivity factor, Bayesian iterations, and the GAUKS biomass reaction.</p></li>
            <li><strong>Run and inspect.</strong><p>The run log records every selected method and returns calibrated parameter groups when MATLAB completes.</p></li>
          </ol>
        </article>

        <article id="tutorial-data" class="page__section">
          <h1>Data matrix</h1>
          <table class="data-table">
            <thead><tr><th>Method</th><th>Required data</th><th>Important parameters</th></tr></thead>
            <tbody>
              <tr><td>Sluice</td><td>exchange reaction list</td><td>exchange reactions</td></tr>
              <tr><td>KcatRepo</td><td>ecGEM with initial kcat values</td><td>repository group name</td></tr>
              <tr><td>Sensitivity</td><td>glucose uptake data and matched growth-rate data</td><td>glucose exchange reaction, target growth, factor, multi-condition flag</td></tr>
              <tr><td>Bayesian</td><td>scenario data</td><td>max iterations, processes, samples per generation, reject ratio</td></tr>
              <tr><td>PRESTO</td><td>proteomics, growth, total protein</td><td>data file handles</td></tr>
              <tr><td>GAUKS</td><td>unconstrained maximum growth data</td><td>GEM ID, biomass reaction</td></tr>
            </tbody>
          </table>
        </article>

        <article id="tutorial-next" class="page__section">
          <h1>Analysis and Design</h1>
          <p>Analysis and Design keep their current web behavior in this round. Use Analysis to inspect ecFVA, knockout, and protein usage outputs; use Design for existing FSEOF, OptKnock, and OptForce workflows. Their deeper redesign will follow after Reconstruction and Calibration stabilize.</p>
        </article>
      </main>
    </section>
  `;
  bindSidebar(main, 'tutorial');
}
