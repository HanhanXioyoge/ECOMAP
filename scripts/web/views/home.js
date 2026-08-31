// mountHome: clones the #home-view-template into <main> and wires up
// recent-projects fetch / language translation.

import { t, getLang } from './shared.js';
import { projectDisplayOrganism } from './projectState.js';

let _switchTrack = null;
export function registerSwitchTrack(fn) { _switchTrack = fn; }

export function mountHome() {
  const main = document.getElementById('main');
  if (!main) return () => {};
  const tpl = document.getElementById('home-view-template');
  if (!tpl) return () => {};

  main.innerHTML = '';
  main.appendChild(tpl.content.cloneNode(true));
  main.dataset.track = '';
  main.classList.add('home-active');

  for (const el of main.querySelectorAll('[data-i18n]')) {
    el.textContent = t(el.getAttribute('data-i18n'), el.textContent);
  }

  fetch('/api/projects?limit=3', { headers: { Accept: 'application/json' } })
    .then((r) => (r.ok ? r.json() : { projects: [] }))
    .then((data) => renderRecentProjects(data.projects || []))
    .catch(() => renderRecentProjects([]));

  return () => {};
}

function _statusClass(status) {
  switch (status) {
    case 'done':
    case 'ready':
      return 'badge--done';
    case 'running':
      return 'badge--running';
    case 'error':
      return 'badge--error';
    case 'queued':
    default:
      return 'badge--queued';
  }
}

function _statusLabel(status) {
  return t(`home_recent_status_${status}`, status || 'ready');
}

function _timeAgo(iso) {
  if (!iso) return '-';
  const then = new Date(iso).getTime();
  if (Number.isNaN(then)) return '-';
  const sec = Math.max(1, Math.floor((Date.now() - then) / 1000));
  const en = getLang() === 'en';
  if (sec < 60) return en ? 'just now' : '刚刚';
  if (sec < 3600) return en ? `${Math.floor(sec / 60)} min ago` : `${Math.floor(sec / 60)} 分钟前`;
  if (sec < 86400) return en ? `${Math.floor(sec / 3600)} h ago` : `${Math.floor(sec / 3600)} 小时前`;
  const days = Math.floor(sec / 86400);
  return en ? `${days} day${days === 1 ? '' : 's'} ago` : `${days} 天前`;
}

function escapeHtml(s) {
  return String(s ?? '').replace(/[&<>"']/g, (c) =>
    ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c])
  );
}

function renderRecentProjects(projects) {
  const list = document.getElementById('recent-projects-list');
  if (!list) return;
  if (!projects.length) {
    list.innerHTML = `<p class="recent__empty" data-i18n="home_recent_empty">${t('home_recent_empty', 'No projects yet.')}</p>`;
    return;
  }
  list.innerHTML = projects.map((project) => {
    const id = project.projectId || project.project_id || '';
    const name = project.name || id;
    const model = [project.modelType || 'Model', project.modelId || ''].filter(Boolean).join(' ');
    return `
      <a class="recent__row" data-project-id="${escapeHtml(id)}" href="#project/${encodeURIComponent(id)}">
        <span class="recent__row-id">${escapeHtml(name)}</span>
        <span class="recent__row-title">${escapeHtml(model)} <small>${escapeHtml(projectDisplayOrganism(project))}</small></span>
        <span class="recent__row-time">${escapeHtml(_timeAgo(project.updatedAt || project.createdAt))}</span>
        <span class="recent__row-algo">${escapeHtml(project.stage || 'Ready')}</span>
        <span><span class="badge ${_statusClass(project.status)}">${escapeHtml(_statusLabel(project.status))}</span></span>
      </a>
    `;
  }).join('');
}
