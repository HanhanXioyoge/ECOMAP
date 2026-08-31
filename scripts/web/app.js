// =============================================================================
// ECOMAP — top-level app shell
//
// Wires the four-track navigation, language toggle, and per-track view
// dispatcher. Per-track views live in ./views/{recon,calib,analysis,design}.js
// and each exports a mountXxx(rootEl) function.
//
// Exposes bindNav / bindLang / switchTrack / applyLang as named exports so
// tests can drive the wiring without booting the full bootstrap pipeline.
// =============================================================================

import { setLang, getLang, t, loadI18n } from './views/shared.js';
import { initRouter }    from './router.js';
import { mountHome, registerSwitchTrack } from './views/home.js';
import { mountProjects } from './views/projects.js';
import { mountRecon }    from './views/recon.js';
import { mountCalib }    from './views/calib.js';
import { mountAnalysis } from './views/analysis.js';
import { mountDesign }   from './views/design.js';

const VIEWS = {
  recon:    mountRecon,
  calib:    mountCalib,
  analysis: mountAnalysis,
  design:   mountDesign,
};

const TRACKS = ['recon', 'calib', 'analysis', 'design'];

function goHash(hash) {
  if (location.hash === hash) {
    window.dispatchEvent(new Event('hashchange'));
  } else {
    location.hash = hash;
  }
}

function replaceHash(hash) {
  if (history?.pushState) {
    try {
      history.pushState(null, '', `${location.pathname}${location.search}${hash}`);
      return;
    } catch (_) {}
  }
  location.hash = hash;
}

// -----------------------------------------------------------------------------
// Bootstrap — runs once on module load (skipped under test when window missing)
// -----------------------------------------------------------------------------

async function bootstrap() {
  // Load the i18n dictionary from the backend. Use the absolute /api/i18n
  // route rather than a relative path into the static tree: the static mount
  // root IS scripts/web/, so a relative '../scripts/web/...' resolves to
  // /scripts/web/... which does not exist and 404s.
  try {
    const [zh, en] = await Promise.all([
      fetch('/api/i18n/zh').then(r => r.ok ? r.json() : null).catch(() => null),
      fetch('/api/i18n/en').then(r => r.ok ? r.json() : null).catch(() => null),
    ]);
    const dict = mergeDicts(zh, en);
    if (Object.keys(dict).length) {
      loadI18n(dict);
    } else {
      // Never fail silently here: an empty dictionary means the language
      // toggle is inert even though the page renders fine from the
      // fallback strings baked into index.html.
      console.warn('[ecomap] i18n dictionary empty - language switching disabled');
    }
  } catch (_) {
    // keep fallback strings in HTML
  }

  // Wire switchTrack into the home view so card clicks work.
  registerSwitchTrack(switchTrack);
  bindNav();
  bindLang();
  applyLang();
  renderLangLabel();
  initRouter();
}

function mergeDicts(zh, en) {
  const dict = {};
  const keys = new Set([...Object.keys(zh ?? {}), ...Object.keys(en ?? {})]);
  for (const k of keys) {
    dict[k] = { zh: zh?.[k] ?? '', en: en?.[k] ?? '' };
  }
  return dict;
}

function currentTrackFromUrl() {
  const m = location.hash.match(/^#track=(recon|calib|analysis|design)$/);
  return m ? m[1] : null;
}

// -----------------------------------------------------------------------------
// Wiring (exported for tests)
// -----------------------------------------------------------------------------

// Track the current view's teardown so switching pages cancels any
// pending fetches / timers / event listeners the previous view set up.
let currentTeardown = null;

export function bindNav() {
  for (const btn of document.querySelectorAll('.topnav-item[data-track]')) {
    const track = btn.dataset.track;
    btn.addEventListener('click', (event) => {
      if (!track) {
        event.preventDefault();
        showHome();
      } else if (VIEWS[track]) {
        event.preventDefault();
        switchTrack(track);
      } else if (track === 'projects') {
        event.preventDefault();
        showProjects();
      }
    });
  }
  // Logo button (the ECOMAP brand) returns to home.
  const brand = document.getElementById('brand-home');
  if (brand) {
    brand.addEventListener('click', () => showHome());
  }
  // Listen for showHome dispatched from shared.js (e.g. track "← Home" link).
  // The handler mounts home via the internal mountHome() and restores the
  // home tab's active state on the topnav.
  window.addEventListener('ecomap:show-home', () => {
    showHome();
    const home = document.getElementById('nav-home');
    if (home) {
      document.querySelectorAll('.topnav-item').forEach((i) => i.classList.remove('is-active'));
      home.classList.add('is-active');
    }
  });
  // Mobile drawer toggle (Phase 8 §8.3). Run once; bindTopnav is idempotent
  // when #topnav-burger / #topnav-drawer are already wired up.
  bindTopnav();
  bindRouteLinks();
}

function bindRouteLinks() {
  if (document.documentElement.dataset.ecomapRouteLinksBound === '1') return;
  document.addEventListener('click', (event) => {
    const link = event.target.closest('a[href="#projects"]');
    if (!link) return;
    event.preventDefault();
    showProjects();
  });
  document.documentElement.dataset.ecomapRouteLinksBound = '1';
}

/** Wire the mobile burger button + drawer overlay (Phase 8 §8.3). */
export function bindTopnav() {
  const burger = document.getElementById('topnav-burger');
  const drawer = document.getElementById('topnav-drawer');
  if (!burger || !drawer) return () => {};
  if (burger.dataset.bound === '1') return () => {};

  // Populate drawer once by cloning the visible topnav items.
  if (!drawer.dataset.ready) {
    const src = document.querySelector('.topnav-items');
    if (src) {
      Array.from(src.children).forEach((child) => {
        const clone = child.cloneNode(true);
        drawer.appendChild(clone);
      });
    }
    drawer.dataset.ready = '1';
  }

  const open  = () => { drawer.classList.add('is-open'); burger.setAttribute('aria-expanded', 'true'); };
  const close = () => { drawer.classList.remove('is-open'); burger.setAttribute('aria-expanded', 'false'); };
  const onBurger = () => (drawer.classList.contains('is-open') ? close() : open());
  const onLinkClick = (e) => {
    const btn = e.target.closest('.topnav-item[data-track]');
    if (btn) {
      if (btn.dataset.track === '') {
        e.preventDefault();
        showHome();
      } else if (VIEWS[btn.dataset.track]) {
        e.preventDefault();
        switchTrack(btn.dataset.track);
      } else if (btn.dataset.track === 'projects') {
        e.preventDefault();
        showProjects();
      }
    }
    if (e.target.closest('a, button')) close();
  };
  burger.addEventListener('click', onBurger);
  drawer.addEventListener('click', onLinkClick);
  burger.dataset.bound = '1';
  return () => {
    burger.removeEventListener('click', onBurger);
    drawer.removeEventListener('click', onLinkClick);
    burger.dataset.bound = '';
  };
}

/** Show the home/landing page. Tears down whatever track view is mounted. */
export function showHome() {
  if (typeof currentTeardown === 'function') {
    currentTeardown();
    currentTeardown = null;
  }
  const main = document.getElementById('main');
  if (main) {
    main.dataset.track = '';
    // The #home-view markup lives statically in index.html; do NOT clear it.
    // Just unmount any track view that may have replaced it via switchTrack.
  }
  // Deactivate every track item; activate the "主页" item if present.
  for (const btn of document.querySelectorAll('.topnav-item')) {
    btn.classList.toggle('is-active', btn.dataset.track === '');
  }
  currentTeardown = mountHome() ?? null;
  if (history && history.replaceState) {
    try { history.replaceState(null, '', location.pathname + location.search); } catch (_) {}
  }
}

export function showProjects() {
  if (typeof currentTeardown === 'function') {
    currentTeardown();
    currentTeardown = null;
  }
  const main = document.getElementById('main');
  if (main) {
    main.innerHTML = '';
    main.dataset.track = 'projects';
    currentTeardown = mountProjects(main) ?? null;
  }
  for (const btn of document.querySelectorAll('.topnav-item')) {
    btn.classList.toggle('is-active', btn.dataset.track === 'projects');
  }
  replaceHash('#projects');
}

export function bindLang() {
  // New segmented pill [EN | 中]. Two buttons share state via getLang().
  const en = document.getElementById('lang-en');
  const zh = document.getElementById('lang-zh');
  if (!en || !zh) return;
  const onPick = (target) => {
    setLang(target);
    applyLang();
    renderLangLabel();
  };
  en.addEventListener('click', () => onPick('en'));
  zh.addEventListener('click', () => onPick('zh'));
}

export function applyLang() {
  document.documentElement.lang = getLang() === 'zh' ? 'zh-CN' : 'en';
  for (const el of document.querySelectorAll('[data-i18n]')) {
    const key = el.getAttribute('data-i18n');
    const fallback = el.textContent;
    el.textContent = t(key, fallback);
  }
  document.title = 'ECOMAP';
}

export function renderLangLabel() {
  // Update the segmented pill's active state.
  const en = document.getElementById('lang-en');
  const zh = document.getElementById('lang-zh');
  if (!en || !zh) return;
  const isEN = getLang() === 'en';
  en.classList.toggle('is-active', isEN);
  en.setAttribute('aria-pressed', isEN ? 'true' : 'false');
  zh.classList.toggle('is-active', !isEN);
  zh.setAttribute('aria-pressed', !isEN ? 'true' : 'false');
}

export function switchTrack(track) {
  if (!VIEWS[track]) return;
  for (const btn of document.querySelectorAll('.topnav-item')) {
    btn.classList.toggle('is-active', btn.dataset.track === track);
  }
  const main = document.getElementById('main');
  if (main) main.dataset.track = track;
  if (typeof currentTeardown === 'function') currentTeardown();
  if (main) main.innerHTML = '';
  currentTeardown = (main ? VIEWS[track](main) : null) ?? null;
  if (history?.replaceState) {
    try { history.replaceState(null, '', `#track=${track}`); } catch (_) {}
  }
}

export function getTracks() {
  return TRACKS.slice();
}

// -----------------------------------------------------------------------------
// Module bootstrap — only fire when there's a real DOM to drive (skips in tests)
// -----------------------------------------------------------------------------

if (typeof window !== 'undefined' && typeof document !== 'undefined') {
  bootstrap().catch((e) => {
    const main = document.getElementById('main');
    if (main) main.innerHTML = `<div class="alert">${e.message}</div>`;
  });
}
