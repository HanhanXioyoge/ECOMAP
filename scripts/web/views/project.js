import { mountAnalysis } from './analysis.js';
import { mountCalib } from './calib.js';
import { mountDesign } from './design.js';
import { mountRecon } from './recon.js';
import { projectDisplayOrganism, projectModuleStatus, rememberCurrentProject } from './projectState.js';

const MODULES = [
  { id: 'reconstruction', title: 'Reconstruction', mount: mountRecon, icon: '<polygon points="12,2 22,22 2,22" fill="currentColor"/>' },
  { id: 'calibration', title: 'Calibration', mount: mountCalib, icon: '<circle cx="12" cy="12" r="10" fill="currentColor"/>' },
  { id: 'analysis', title: 'Analysis', mount: mountAnalysis, icon: '<rect x="3" y="3" width="8" height="8" fill="currentColor"/><rect x="13" y="3" width="8" height="8" fill="currentColor"/><rect x="3" y="13" width="8" height="8" fill="currentColor"/><rect x="13" y="13" width="8" height="8" fill="currentColor"/>' },
  { id: 'design', title: 'Design', mount: mountDesign, icon: '<polygon points="2,2 22,12 2,22" fill="currentColor"/>' },
];

function escapeHtml(value) {
  return String(value ?? '').replace(/[&<>"']/g, (ch) => ({
    '&': '&amp;',
    '<': '&lt;',
    '>': '&gt;',
    '"': '&quot;',
    "'": '&#39;',
  }[ch]));
}

function renderOverview(root, project) {
  root.innerHTML = `
    <section class="page page--project">
      <header class="page__header">
        <p class="section-eyebrow">CURRENT PROJECT</p>
        <h1>${escapeHtml(project.name)}</h1>
        <p>${escapeHtml(project.stage || 'Ready')}</p>
      </header>

      <section class="project-summary">
        <div class="project-summary__item"><span>Model type</span><strong>${escapeHtml(project.modelType)}</strong></div>
        <div class="project-summary__item"><span>Model ID</span><strong>${escapeHtml(project.modelId)}</strong></div>
        <div class="project-summary__item"><span>Organism</span><strong>${escapeHtml(projectDisplayOrganism(project))}</strong></div>
      </section>

      <section class="modules-grid project-modules">
        ${MODULES.map((item) => `
          <a class="module-card" href="#project/${encodeURIComponent(project.projectId)}/${item.id}">
            <svg class="module-card__icon module__icon" viewBox="0 0 24 24" aria-hidden="true">${item.icon}</svg>
            <h2 class="module-card__title">${item.title}</h2>
            <p class="module-card__desc">${escapeHtml(projectModuleStatus(project, item.id))}</p>
            <span class="module-card__cta">Open module -></span>
          </a>
        `).join('')}
      </section>

      <section class="workbench-card project-params">
        <h2>Project parameters</h2>
        <pre class="param-preview">${escapeHtml(JSON.stringify(project.params || {}, null, 2))}</pre>
      </section>
    </section>
  `;
}

function renderProjectDetail(root, project, view) {
  const titleMap = {
    parameters: 'Project parameters',
    files: 'Project files',
    runs: 'Project runs',
  };
  const files = Object.entries(project.files || {});
  const runs = project.runs || [];
  const body = view === 'parameters'
    ? `<pre class="param-preview">${escapeHtml(JSON.stringify(project.params || {}, null, 2))}</pre>`
    : view === 'files'
      ? (files.length
        ? `<dl class="project-card__meta">${files.map(([key, value]) => `<div><dt>${escapeHtml(key)}</dt><dd>${escapeHtml(value)}</dd></div>`).join('')}</dl>`
        : '<p class="recent__empty">No project files recorded yet.</p>')
      : (runs.length
        ? `<pre class="param-preview">${escapeHtml(JSON.stringify(runs, null, 2))}</pre>`
        : '<p class="recent__empty">No runs recorded yet.</p>');
  root.innerHTML = `
    <section class="page page--project">
      <header class="page__header">
        <p class="section-eyebrow">CURRENT PROJECT</p>
        <h1>${escapeHtml(titleMap[view] || 'Project detail')}</h1>
        <p>${escapeHtml(project.name)}</p>
      </header>
      <section class="workbench-card project-params">
        ${body}
      </section>
      <a class="btn btn-secondary" href="#project/${encodeURIComponent(project.projectId)}">Back to overview</a>
    </section>
  `;
}

function mountModule(root, project, view) {
  const mod = MODULES.find((item) => item.id === view);
  if (!mod) {
    renderOverview(root, project);
    return;
  }
  const teardown = mod.mount(root, {
    project,
    projectParams: project.params || {},
    modelId: project.modelType === 'GEM' ? project.modelId : '',
    gemModelId: project.modelType === 'GEM' ? project.modelId : project.gemModelId || '',
    ecModelId: project.modelType === 'ecGEM' ? project.modelId : project.ecModelId || '',
  });
  const banner = document.createElement('div');
  banner.className = 'project-context-banner';
  banner.innerHTML = `
    <a href="#project/${encodeURIComponent(project.projectId)}">Project overview</a>
    <span>${escapeHtml(project.name)}</span>
    <strong>${escapeHtml(project.stage || 'Ready')}</strong>
  `;
  root.prepend(banner);
  return teardown;
}

export function mountProject(main, projectId, view = 'overview') {
  if (!projectId) {
    main.innerHTML = '<section class="page"><p>No project selected.</p><a class="btn btn-secondary" href="#projects">Open projects</a></section>';
    return;
  }
  main.innerHTML = '<section class="page"><p class="recent__empty">Loading project...</p></section>';
  fetch(`/api/projects/${encodeURIComponent(projectId)}`, { headers: { Accept: 'application/json' } })
    .then((res) => (res.ok ? res.json() : null))
    .then((project) => {
      if (!project) {
        main.innerHTML = '<section class="page"><p>Project not found.</p><a class="btn btn-secondary" href="#projects">Back to projects</a></section>';
        return;
      }
      rememberCurrentProject(project.projectId);
      if (view === 'overview') renderOverview(main, project);
      else if (['parameters', 'files', 'runs'].includes(view)) renderProjectDetail(main, project, view);
      else mountModule(main, project, view);
    })
    .catch(() => {
      main.innerHTML = '<section class="page"><p>Project not found.</p><a class="btn btn-secondary" href="#projects">Back to projects</a></section>';
    });
}
