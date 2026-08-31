import { deepLearningKcatHtml } from './shared.js';

const SECTIONS = [
  { id: 'overview', title: 'Overview' },
  { id: 'architecture', title: 'ecGEM architecture' },
  { id: 'kcat', title: 'kcat assignment' },
  { id: 'calibration', title: 'Calibration' },
  { id: 'modules', title: 'Modules' },
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

export function mountDocs(main) {
  main.innerHTML = `
    <section class="page page--docs">
      <aside class="page__sidebar">
        <h2 class="page__sidebar-title">Docs</h2>
        <ul class="page__sidebar-list" role="list">${sidebarLinks('docs')}</ul>
      </aside>
      <main class="page__main">
        <article id="docs-overview" class="page__section">
          <h1>ECOMAP documentation</h1>
          <p>ECOMAP (Constrained Optimization of Metabolic Models & Analysis Pipeline) supports automated reconstruction, calibration, analysis, and rational design of enzyme-constrained genome-scale metabolic models. The current web direction keeps the paper's scientific pipeline, then wraps every run in an independent project with one uploaded GEM or ecGEM and one shared parameter set.</p>
          <p>The central artifact is an ecGEM: a genome-scale metabolic model extended with enzyme usage, molecular weight, turnover number, and a finite proteome budget.</p>
          <figure class="figure-placeholder">Figure placeholder: paper workflow overview.</figure>
        </article>

        <article id="docs-architecture" class="page__section">
          <h1>ecGEM architecture</h1>
          <p>Reconstruction supports basic, isozyme, and integrated S-matrix layouts. The basic form adds enzyme capacity at reaction level, the isozyme form separates alternative enzymes, and the integrated form preserves enzyme choice while carrying a unified protein pool constraint.</p>
          <p>This lets users compare model structures before committing to calibration. In the web app, a project created from a GEM starts in Reconstruction; a project created from an ecGEM can go directly to Calibration, Analysis, or Design.</p>
          <figure class="figure-placeholder">Figure placeholder: three S-matrix architectures from the paper.</figure>
        </article>

        <article id="docs-kcat" class="page__section">
          <h1>kcat assignment</h1>
          <p>Turnover numbers can come from curated databases, user measurements, and ${deepLearningKcatHtml()} prediction models such as DLKcat, UniKP, and CatPred. Reconstruction compares these sources and merges them into the ecGEM before growth evaluation.</p>
          <p>When multiple sources exist, the web workflow keeps the methods independent: prediction, comparison, and merge can be run separately so users can inspect intermediate assumptions.</p>
          <figure class="figure-placeholder">Figure placeholder: hierarchical kcat assignment and prediction comparison.</figure>
        </article>

        <article id="docs-calibration" class="page__section">
          <h1>Calibration</h1>
          <p>Calibration refines the ecGEM with method-specific evidence. Sluice prepares the model structure, KcatRepo manages parameter groups, Sensitivity uses glucose uptake and matched growth-rate data, Bayesian calibration fits scenario data, PRESTO uses proteomics/growth/total-protein data, and GAUKS requires unconstrained maximum growth data.</p>
          <p>Selected methods run in a fixed order: Sluice, KcatRepo, Sensitivity, Bayesian, PRESTO, and GAUKS. This keeps dependencies predictable while still allowing each module to run separately inside the current project.</p>
          <figure class="figure-placeholder">Figure placeholder: calibration methods and data dependencies.</figure>
        </article>

        <article id="docs-modules" class="page__section">
          <h1>Modules</h1>
          <ol class="module-list">
            <li><strong>Reconstruction</strong> converts GEM files into enzyme-constrained ecGEMs.</li>
            <li><strong>Calibration</strong> tunes enzyme parameters against growth, uptake, proteomics, and unconstrained growth evidence.</li>
            <li><strong>Analysis</strong> remains available in its current form for ecFVA, knockout, and protein-usage inspection.</li>
            <li><strong>Design</strong> remains available in its current form for FSEOF, OptKnock, and OptForce-oriented strain design.</li>
          </ol>
          <p>The Modules page is now a catalog. Operational work begins from Projects, where the Current Project menu exposes Overview, Reconstruction, Calibration, Analysis, Design, Parameters, Files, and Runs.</p>
        </article>
      </main>
    </section>
  `;
  bindSidebar(main, 'docs');
}
