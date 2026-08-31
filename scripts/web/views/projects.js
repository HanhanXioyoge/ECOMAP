import {
  buildProjectPayload,
  createProjectDraft,
  projectDisplayOrganism,
  rememberCurrentProject,
  validateProjectDraft,
} from './projectState.js';

function escapeHtml(value) {
  return String(value ?? '').replace(/[&<>"']/g, (ch) => ({
    '&': '&amp;',
    '<': '&lt;',
    '>': '&gt;',
    '"': '&quot;',
    "'": '&#39;',
  }[ch]));
}

function projectRow(project) {
  return `
    <article class="project-card">
      <div>
        <h2>${escapeHtml(project.name)}</h2>
        <p>${escapeHtml(project.stage || 'Ready')}</p>
      </div>
      <dl class="project-card__meta">
        <div><dt>Model</dt><dd>${escapeHtml(project.modelType || '')}</dd></div>
        <div><dt>Organism</dt><dd>${escapeHtml(projectDisplayOrganism(project))}</dd></div>
        <div><dt>Status</dt><dd>${escapeHtml(project.status || 'ready')}</dd></div>
      </dl>
      <a class="btn btn-secondary" href="#project/${encodeURIComponent(project.projectId)}">Open project</a>
    </article>
  `;
}

async function loadProjects(root) {
  const list = root.querySelector('[data-role="projects-list"]');
  if (!list) return;
  list.innerHTML = '<p class="recent__empty">Loading projects...</p>';
  try {
    const res = await fetch('/api/projects?limit=20', { headers: { Accept: 'application/json' } });
    const data = res.ok ? await res.json() : { projects: [] };
    const projects = data.projects || [];
    list.innerHTML = projects.length
      ? projects.map(projectRow).join('')
      : '<p class="recent__empty">No projects yet.</p>';
  } catch (_) {
    list.innerHTML = '<p class="recent__empty">No projects yet.</p>';
  }
}

function readDraft(root, draft) {
  draft.name = root.querySelector('[data-role="project-name"]')?.value.trim() || '';
  draft.modelType = root.querySelector('[data-role="project-model-type"]')?.value || '';
  return draft;
}

export function mountProjects(main) {
  const draft = createProjectDraft();
  main.innerHTML = `
    <section class="page page--projects">
      <header class="page__header">
        <h1>Projects</h1>
        <p>Create a local project workspace first, then configure parameters and run modules inside that project.</p>
      </header>

      <section class="project-create">
        <h2>New project</h2>
        <div class="project-form">
          <label><span class="label">Project name</span><input class="input" data-role="project-name" placeholder="Ecoli demo" /></label>
          <label><span class="label">Input model type</span>
            <select class="input" data-role="project-model-type">
              <option value="">Select type</option>
              <option value="GEM">GEM</option>
              <option value="ecGEM">ecGEM</option>
            </select>
          </label>
          <p class="help-text" data-role="project-create-status">The model file is uploaded inside the selected module.</p>
          <button type="button" class="btn btn-primary" data-role="create-project">Create project</button>
        </div>
      </section>

      <section class="projects-list-shell">
        <h2>All projects</h2>
        <div class="projects-list" data-role="projects-list"></div>
      </section>
    </section>
  `;

  main.querySelector('[data-role="create-project"]')?.addEventListener('click', async () => {
    readDraft(main, draft);
    const status = main.querySelector('[data-role="project-create-status"]');
    const errors = validateProjectDraft(draft);
    if (errors.length) {
      if (status) status.textContent = `Missing: ${errors.join(', ')}`;
      return;
    }
    const res = await fetch('/api/projects', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(buildProjectPayload(draft)),
    });
    const created = await res.json().catch(() => ({}));
    if (res.ok && created.projectId) {
      rememberCurrentProject(created.projectId);
      location.hash = `#project/${created.projectId}`;
    } else if (status) {
      status.textContent = created.detail || created.error || 'Project creation failed.';
      loadProjects(main);
    }
  });

  loadProjects(main);
}
