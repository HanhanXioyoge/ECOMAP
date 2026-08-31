// Hash router: dispatches window.location.hash to a mount function.
// Each mount*() is responsible for clearing <main> before re-rendering.

import { mountHome }        from './views/home.js';
import { mountDocs }        from './views/docs.js';
import { mountTutorial }    from './views/tutorial.js';
import { mountModules }     from './views/modules.js';
import { mountProjects }    from './views/projects.js';
import { mountProject }     from './views/project.js';
import { mountRecon }       from './views/recon.js';
import { mountCalib }       from './views/calib.js';
import { mountAnalysis }    from './views/analysis.js';
import { mountDesign }      from './views/design.js';
import { getCurrentProjectId } from './views/projectState.js';

const routes = {
  '':                mountHome,
  '#docs':           mountDocs,
  '#tutorial':       mountTutorial,
  '#modules':        mountModules,
  '#projects':       mountProjects,
  '#reconstruction': mountProjects,
  '#calibration':    mountProjects,
  '#analysis':       mountProjects,
  '#design':         mountProjects,
};

function navigate() {
  const hash = window.location.hash || '';
  const main = document.getElementById('main');
  if (main) main.innerHTML = '';
  try {
    const projectRoute = parseProjectRoute(hash);
    if (projectRoute) {
      mountProject(main, projectRoute.projectId, projectRoute.view);
    } else {
      const mount = routes[hash] || routes[''];
      mount(main);
    }
  } catch (err) {
    console.error('[router] mount failed for', hash, err);
  }
  updateCurrentProjectNav(hash);
  // Update topnav active state.
  document.querySelectorAll('.topnav-item').forEach((btn) => {
    const track = btn.getAttribute('data-track') || '';
    const wantsActive =
      (track === '' && (hash === '' || hash === '#')) ||
      (track === 'projects' && (hash === '#projects' || hash.startsWith('#project/'))) ||
      (track !== '' && track !== 'projects' && hash === '#' + track);
    btn.classList.toggle('is-active', wantsActive);
  });
}

function parseProjectRoute(hash) {
  const m = String(hash || '').match(/^#project\/([^/]+)(?:\/([^/]+))?$/);
  if (!m) {
    const current = String(hash || '').match(/^#project-current(?:\/([^/]+))?$/);
    if (!current) return null;
    const projectId = getCurrentProjectId();
    return projectId ? { projectId, view: current[1] || 'overview' } : { projectId: '', view: 'overview' };
  }
  return {
    projectId: decodeURIComponent(m[1]),
    view: decodeURIComponent(m[2] || 'overview'),
  };
}

function updateCurrentProjectNav(hash) {
  const currentId = parseProjectRoute(hash)?.projectId || getCurrentProjectId();
  document.querySelectorAll('#nav-current-project-wrap').forEach((wrap) => {
    const link = wrap.querySelector('#nav-current-project');
    wrap.classList.toggle('has-current-project', Boolean(currentId));
    if (!link) return;
    if (!currentId) {
      link.setAttribute('href', '#projects');
      wrap.querySelectorAll('[data-project-view]').forEach((item) => {
        item.setAttribute('href', '#projects');
      });
      return;
    }
    link.setAttribute('href', `#project/${encodeURIComponent(currentId)}`);
    wrap.querySelectorAll('[data-project-view]').forEach((item) => {
      const view = item.getAttribute('data-project-view') || 'overview';
      item.setAttribute('href', view === 'overview'
        ? `#project/${encodeURIComponent(currentId)}`
        : `#project/${encodeURIComponent(currentId)}/${view}`);
    });
  });
}

export function initRouter() {
  window.addEventListener('hashchange', navigate);
  window.addEventListener('popstate',  navigate);
  navigate();
}
