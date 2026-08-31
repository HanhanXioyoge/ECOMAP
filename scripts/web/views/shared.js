// =============================================================================
// ECOMAP / MDP — shared view helpers
//
// ES module consumed by every track view. Provides:
//   - i18n dictionary + DOM refresher (B13 will replace with proper loader)
//   - SSE job watcher
//   - data table renderer
//   - canvas2D chart drawers (FVA, RMSE history, scatter regression, growth)
//   - file uploader
//
// No external dependencies. All charts are hand-rolled against <canvas> 2D.
// =============================================================================

// -----------------------------------------------------------------------------
// i18n — placeholder until B13 wires the real bilingual loader.
// -----------------------------------------------------------------------------

const LANG_KEY = 'ecomap_lang';

function readStoredLang() {
  try {
    const v = localStorage.getItem(LANG_KEY);
    return (v === 'en' || v === 'zh') ? v : null;
  } catch (_) { return null; }
}

let _i18n = { lang: readStoredLang() || 'en', dict: {} };

export function setLang(lang) {
  _i18n.lang = (lang === 'en') ? 'en' : 'zh';
  try { localStorage.setItem(LANG_KEY, _i18n.lang); } catch (_) {}
  refreshI18nDom();
}

export function getLang() {
  return _i18n.lang;
}

export function t(key, fallback) {
  const v = (_i18n.dict && _i18n.dict[key]) ? _i18n.dict[key][_i18n.lang] : undefined;
  if (v == null || v === '') {
    return (fallback != null) ? fallback : key;
  }
  return v;
}

export function loadI18n(dict) {
  _i18n.dict = dict || {};
  refreshI18nDom();
}

// -----------------------------------------------------------------------------
// Shared module run parameters
// -----------------------------------------------------------------------------

export function createParameterManager(opts = {}) {
  const defaults = { ...(opts.defaults || {}) };
  const projectParams = { ...(opts.projectParams || {}) };
  let runParams = { ...(opts.runParams || {}) };

  const merged = () => ({ ...defaults, ...projectParams, ...runParams });

  return {
    get(key) {
      return merged()[key];
    },
    set(key, value) {
      runParams[key] = value;
      return this.snapshot();
    },
    merge(values = {}) {
      runParams = { ...runParams, ...values };
      return this.snapshot();
    },
    resetRun(values = {}) {
      runParams = { ...(values || {}) };
      return this.snapshot();
    },
    snapshot() {
      return merged();
    },
  };
}

export function deepLearningKcatHtml() {
  return 'Deep Learning <em>k</em><sub>cat</sub>';
}

// -----------------------------------------------------------------------------
// App-shell navigation helpers
// -----------------------------------------------------------------------------
//
// showHome() is dispatched from places where shared.js has a UI hook into
// the app shell (e.g. the track-page "← Home" back link) without creating a
// circular import on app.js. The actual home-mount logic lives in app.js
// and listens for the `ecomap:show-home` window event.
export function showHome() {
  if (typeof window === 'undefined') return;
  window.dispatchEvent(new CustomEvent('ecomap:show-home'));
}

function refreshI18nDom() {
  if (typeof document === 'undefined') return;
  const nodes = document.querySelectorAll('[data-i18n]');
  for (const el of nodes) {
    const key = el.getAttribute('data-i18n');
    el.textContent = t(key, el.textContent);
  }
}

// -----------------------------------------------------------------------------
// Jobs / SSE
// -----------------------------------------------------------------------------

export function watchJob(jid, onEvent) {
  const src = new EventSource(`/api/jobs/${encodeURIComponent(jid)}/events`);
  src.onmessage = (msg) => {
    try {
      const data = JSON.parse(msg.data);
      onEvent(data);
    } catch (err) {
      // best-effort: surface a structured error so the caller can decide
      onEvent({ event: 'sse-parse-error', message: String(err && err.message || err) });
    }
  };
  src.onerror = () => {
    onEvent({ event: 'sse-error', message: 'connection lost' });
  };
  return () => src.close();
}

// -----------------------------------------------------------------------------
// Rendering — data table
// -----------------------------------------------------------------------------

export function renderTable(rootEl, headers, rows) {
  if (!rootEl) return;
  // Clear previous content.
  rootEl.innerHTML = '';

  const table = document.createElement('table');
  table.className = 'data-table';

  // Header row
  const thead = document.createElement('thead');
  const headRow = document.createElement('tr');
  for (const h of headers) {
    const th = document.createElement('th');
    th.textContent = String(h);
    headRow.appendChild(th);
  }
  thead.appendChild(headRow);
  table.appendChild(thead);

  // Body
  const tbody = document.createElement('tbody');
  if (!rows || rows.length === 0) {
    const tr = document.createElement('tr');
    const td = document.createElement('td');
    td.colSpan = headers.length;
    td.className = 'text-mute';
    td.textContent = 'No data';
    tr.appendChild(td);
    tbody.appendChild(tr);
  } else {
    for (const row of rows) {
      const tr = document.createElement('tr');
      for (const h of headers) {
        const td = document.createElement('td');
        const cell = row ? row[h] : '';
        td.textContent = (cell == null) ? '' : String(cell);
        tr.appendChild(td);
      }
      tbody.appendChild(tr);
    }
  }
  table.appendChild(tbody);

  rootEl.appendChild(table);
}

// -----------------------------------------------------------------------------
// Canvas2D drawing helpers — small, dependency-free charts.
// -----------------------------------------------------------------------------

const TRACK_COLORS = {
  recon:    '#2563EB',
  calib:    '#7C3AED',
  analysis: '#0891B2',
  design:   '#EA580C',
  success:  '#16A34A',
  warn:     '#D97706',
  error:    '#DC2626',
  mute:     '#64748B',
  grid:     'rgba(15, 23, 42, 0.08)',
};

function pickCtx(canvasEl) {
  if (!canvasEl) throw new Error('canvas element required');
  const ctx = canvasEl.getContext('2d');
  if (!ctx) throw new Error('2d context unavailable');
  return ctx;
}

function sizeCanvas(canvasEl) {
  // Honor the rendered CSS size while keeping a sensible device-pixel backing.
  const rect = canvasEl.getBoundingClientRect();
  const dpr = (typeof window !== 'undefined' && window.devicePixelRatio) || 1;
  const cssW = Math.max(1, Math.round(rect.width));
  const cssH = Math.max(1, Math.round(rect.height));
  if (canvasEl.width !== cssW * dpr || canvasEl.height !== cssH * dpr) {
    canvasEl.width = cssW * dpr;
    canvasEl.height = cssH * dpr;
  }
  return { ctx: canvasEl.getContext('2d'), w: cssW, h: cssH, dpr };
}

function clearAndPrepare(ctx, w, h, dpr) {
  ctx.save();
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
}

function drawAxes(ctx, w, h, padding) {
  ctx.strokeStyle = TRACK_COLORS.grid;
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(padding.left, padding.top);
  ctx.lineTo(padding.left, h - padding.bottom);
  ctx.lineTo(w - padding.right, h - padding.bottom);
  ctx.stroke();
}

function drawAxisLabel(ctx, text, x, y, align = 'center') {
  ctx.fillStyle = TRACK_COLORS.mute;
  ctx.font = '11px ui-monospace, Consolas, monospace';
  ctx.textAlign = align;
  ctx.textBaseline = 'middle';
  ctx.fillText(text, x, y);
}

function niceTicks(min, max, count) {
  if (!isFinite(min) || !isFinite(max) || min === max) {
    return [min || 0];
  }
  const span = max - min;
  const step0 = span / Math.max(1, count);
  const mag = Math.pow(10, Math.floor(Math.log10(step0)));
  const norm = step0 / mag;
  let step;
  if (norm < 1.5) step = 1 * mag;
  else if (norm < 3) step = 2 * mag;
  else if (norm < 7) step = 5 * mag;
  else step = 10 * mag;
  const start = Math.ceil(min / step) * step;
  const ticks = [];
  for (let v = start; v <= max + step * 0.5; v += step) {
    ticks.push(Number(v.toFixed(10)));
  }
  return ticks;
}

// --- FVA bar ----------------------------------------------------------------
//
// ranges: [{ min, max, mean }, ...]
// labels: string[]
export function drawFvaBar(canvasEl, ranges, labels) {
  const { ctx, w, h, dpr } = sizeCanvas(canvasEl);
  clearAndPrepare(ctx, w, h, dpr);

  const padding = { top: 16, right: 16, bottom: 16, left: 120 };
  const n = (ranges || []).length;
  if (n === 0) {
    drawAxisLabel(ctx, 'No FVA data', w / 2, h / 2);
    ctx.restore();
    return;
  }
  const rowH = (h - padding.top - padding.bottom) / n;
  const barAreaW = w - padding.left - padding.right;

  let globalMin = Infinity, globalMax = -Infinity;
  for (const r of ranges) {
    if (r.min < globalMin) globalMin = r.min;
    if (r.max > globalMax) globalMax = r.max;
  }
  if (!isFinite(globalMin) || !isFinite(globalMax) || globalMin === globalMax) {
    globalMin = (globalMin || 0) - 1;
    globalMax = (globalMax || 0) + 1;
  }
  const pad = (globalMax - globalMin) * 0.05;
  const x0 = globalMin - pad;
  const x1 = globalMax + pad;
  const xToPx = (v) => padding.left + ((v - x0) / (x1 - x0)) * barAreaW;

  // gridlines + ticks
  const ticks = niceTicks(globalMin, globalMax, 5);
  ctx.strokeStyle = TRACK_COLORS.grid;
  ctx.lineWidth = 1;
  ctx.font = '10px ui-monospace, Consolas, monospace';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  ctx.fillStyle = TRACK_COLORS.mute;
  for (const t of ticks) {
    const x = xToPx(t);
    ctx.beginPath();
    ctx.moveTo(x, padding.top);
    ctx.lineTo(x, h - padding.bottom);
    ctx.stroke();
    ctx.fillText(Number(t.toFixed(3)).toString(), x, h - padding.bottom + 2);
  }

  // bars
  for (let i = 0; i < n; i++) {
    const r = ranges[i];
    const cy = padding.top + rowH * (i + 0.5);
    const barH = Math.max(8, rowH * 0.5);

    // label
    ctx.fillStyle = TRACK_COLORS.mute;
    ctx.font = '11px Inter, sans-serif';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    ctx.fillText(String(labels[i] || ''), padding.left - 8, cy);

    // translucent range band
    ctx.fillStyle = TRACK_COLORS.analysis + '33'; // ~20% opacity
    const xMin = xToPx(r.min);
    const xMax = xToPx(r.max);
    ctx.fillRect(Math.min(xMin, xMax), cy - barH / 2, Math.abs(xMax - xMin), barH);

    // mean line
    ctx.strokeStyle = TRACK_COLORS.analysis;
    ctx.lineWidth = 2;
    ctx.beginPath();
    const xMean = xToPx(r.mean);
    ctx.moveTo(xMean, cy - barH / 2 - 2);
    ctx.lineTo(xMean, cy + barH / 2 + 2);
    ctx.stroke();
  }

  ctx.restore();
}

// --- RMSE history line chart -----------------------------------------------
//
// points: [{ iteration, rmse }, ...]   (iteration is x, rmse is y)
export function drawRmseHistory(canvasEl, points) {
  const { ctx, w, h, dpr } = sizeCanvas(canvasEl);
  clearAndPrepare(ctx, w, h, dpr);

  const padding = { top: 16, right: 16, bottom: 32, left: 48 };
  drawAxes(ctx, w, h, padding);

  const pts = (points || []).filter(p => isFinite(p.iteration) && isFinite(p.rmse));
  if (pts.length === 0) {
    drawAxisLabel(ctx, 'No RMSE history', w / 2, h / 2);
    drawAxisLabel(ctx, 'Iteration', (padding.left + w - padding.right) / 2, h - 8);
    drawAxisLabel(ctx, 'RMSE', 12, padding.top + 4, 'left');
    ctx.restore();
    return;
  }

  let xMin = pts[0].iteration, xMax = pts[0].iteration;
  let yMin = pts[0].rmse, yMax = pts[0].rmse;
  for (const p of pts) {
    if (p.iteration < xMin) xMin = p.iteration;
    if (p.iteration > xMax) xMax = p.iteration;
    if (p.rmse < yMin) yMin = p.rmse;
    if (p.rmse > yMax) yMax = p.rmse;
  }
  if (xMin === xMax) xMax = xMin + 1;
  const yPad = Math.max((yMax - yMin) * 0.1, 1e-6);
  yMin -= yPad; yMax += yPad;
  const xPad = (xMax - xMin) * 0.02;
  xMin -= xPad; xMax += xPad;

  const xToPx = (v) => padding.left + ((v - xMin) / (xMax - xMin)) * (w - padding.left - padding.right);
  const yToPx = (v) => padding.top + (1 - (v - yMin) / (yMax - yMin)) * (h - padding.top - padding.bottom);

  // grid
  ctx.strokeStyle = TRACK_COLORS.grid;
  ctx.lineWidth = 1;
  const yTicks = niceTicks(yMin, yMax, 4);
  for (const t of yTicks) {
    const y = yToPx(t);
    ctx.beginPath();
    ctx.moveTo(padding.left, y);
    ctx.lineTo(w - padding.right, y);
    ctx.stroke();
  }

  // line
  ctx.strokeStyle = TRACK_COLORS.calib;
  ctx.lineWidth = 2;
  ctx.beginPath();
  for (let i = 0; i < pts.length; i++) {
    const x = xToPx(pts[i].iteration);
    const y = yToPx(pts[i].rmse);
    if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  }
  ctx.stroke();

  // points
  ctx.fillStyle = TRACK_COLORS.calib;
  for (const p of pts) {
    ctx.beginPath();
    ctx.arc(xToPx(p.iteration), yToPx(p.rmse), 3, 0, Math.PI * 2);
    ctx.fill();
  }

  // axis labels
  drawAxisLabel(ctx, 'Iteration', (padding.left + w - padding.right) / 2, h - 12);
  ctx.save();
  ctx.translate(12, (padding.top + h - padding.bottom) / 2);
  ctx.rotate(-Math.PI / 2);
  drawAxisLabel(ctx, 'RMSE', 0, 0);
  ctx.restore();

  ctx.restore();
}

// --- Scatter + regression line ---------------------------------------------
//
// points: [{ predicted, measured }, ...]   (log10 kcat values)
export function drawScatterRegression(canvasEl, points) {
  const { ctx, w, h, dpr } = sizeCanvas(canvasEl);
  clearAndPrepare(ctx, w, h, dpr);

  const padding = { top: 16, right: 16, bottom: 36, left: 48 };
  drawAxes(ctx, w, h, padding);

  const pts = (points || []).filter(p => isFinite(p.predicted) && isFinite(p.measured));
  if (pts.length === 0) {
    drawAxisLabel(ctx, 'No scatter data', w / 2, h / 2);
    ctx.restore();
    return;
  }

  let mn = Math.min(...pts.map(p => Math.min(p.predicted, p.measured)));
  let mx = Math.max(...pts.map(p => Math.max(p.predicted, p.measured)));
  const pad = Math.max((mx - mn) * 0.05, 1e-6);
  mn -= pad; mx += pad;

  const xToPx = (v) => padding.left + ((v - mn) / (mx - mn)) * (w - padding.left - padding.right);
  const yToPx = (v) => padding.top + (1 - (v - mn) / (mx - mn)) * (h - padding.top - padding.bottom);

  // grid
  ctx.strokeStyle = TRACK_COLORS.grid;
  ctx.lineWidth = 1;
  const ticks = niceTicks(mn, mx, 5);
  for (const t of ticks) {
    const x = xToPx(t), y = yToPx(t);
    ctx.beginPath();
    ctx.moveTo(padding.left, y);
    ctx.lineTo(w - padding.right, y);
    ctx.moveTo(x, padding.top);
    ctx.lineTo(x, h - padding.bottom);
    ctx.stroke();
  }

  // identity y = x reference
  ctx.strokeStyle = TRACK_COLORS.mute;
  ctx.lineWidth = 1;
  ctx.setLineDash([4, 4]);
  ctx.beginPath();
  ctx.moveTo(xToPx(mn), yToPx(mn));
  ctx.lineTo(xToPx(mx), yToPx(mx));
  ctx.stroke();
  ctx.setLineDash([]);

  // least-squares fit y = m*x + b
  let m = 0, b = 0;
  if (pts.length >= 2) {
    const n = pts.length;
    let sx = 0, sy = 0, sxy = 0, sxx = 0;
    for (const p of pts) {
      sx += p.predicted; sy += p.measured;
      sxy += p.predicted * p.measured;
      sxx += p.predicted * p.predicted;
    }
    const denom = (n * sxx - sx * sx) || 1e-12;
    m = (n * sxy - sx * sy) / denom;
    b = (sy - m * sx) / n;
  }

  // regression line clipped to axis range
  const fitY1 = m * mn + b;
  const fitY2 = m * mx + b;
  ctx.strokeStyle = TRACK_COLORS.recon;
  ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.moveTo(xToPx(mn), yToPx(fitY1));
  ctx.lineTo(xToPx(mx), yToPx(fitY2));
  ctx.stroke();

  // scatter dots
  ctx.fillStyle = TRACK_COLORS.recon + 'AA';
  for (const p of pts) {
    ctx.beginPath();
    ctx.arc(xToPx(p.predicted), yToPx(p.measured), 3, 0, Math.PI * 2);
    ctx.fill();
  }

  // axis labels
  drawAxisLabel(ctx, 'Predicted kcat (log10)', (padding.left + w - padding.right) / 2, h - 14);
  ctx.save();
  ctx.translate(12, (padding.top + h - padding.bottom) / 2);
  ctx.rotate(-Math.PI / 2);
  drawAxisLabel(ctx, 'Measured kcat (log10)', 0, 0);
  ctx.restore();

  ctx.restore();
}

// --- Growth-curve overlay --------------------------------------------------
//
// seriesMap: { ref: [{t, v}, ...], mut1: [...], ... }
// opts:      { xLabel = 'Time (h)', yLabel = 'OD600' }
export function drawGrowthCompare(canvasEl, seriesMap, opts) {
  const { ctx, w, h, dpr } = sizeCanvas(canvasEl);
  clearAndPrepare(ctx, w, h, dpr);

  opts = opts || {};
  const xLabel = opts.xLabel || 'Time (h)';
  const yLabel = opts.yLabel || 'OD600';

  const padding = { top: 16, right: 110, bottom: 32, left: 48 };
  drawAxes(ctx, w, h, padding);

  const names = Object.keys(seriesMap || {});
  const allPoints = [];
  for (const n of names) {
    for (const p of (seriesMap[n] || [])) {
      if (isFinite(p.t) && isFinite(p.v)) allPoints.push(p);
    }
  }
  if (allPoints.length === 0) {
    drawAxisLabel(ctx, 'No growth data', w / 2, h / 2);
    ctx.restore();
    return;
  }

  let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity;
  for (const p of allPoints) {
    if (p.t < xMin) xMin = p.t;
    if (p.t > xMax) xMax = p.t;
    if (p.v < yMin) yMin = p.v;
    if (p.v > yMax) yMax = p.v;
  }
  if (xMin === xMax) xMax = xMin + 1;
  const yPad = Math.max((yMax - yMin) * 0.08, 1e-6);
  yMin -= yPad; yMax += yPad;

  const xToPx = (v) => padding.left + ((v - xMin) / (xMax - xMin)) * (w - padding.left - padding.right);
  const yToPx = (v) => padding.top + (1 - (v - yMin) / (yMax - yMin)) * (h - padding.top - padding.bottom);

  // grid
  ctx.strokeStyle = TRACK_COLORS.grid;
  ctx.lineWidth = 1;
  const yTicks = niceTicks(yMin, yMax, 4);
  for (const t of yTicks) {
    const y = yToPx(t);
    ctx.beginPath();
    ctx.moveTo(padding.left, y);
    ctx.lineTo(w - padding.right, y);
    ctx.stroke();
  }

  // series palette — prefer track accents in declared order
  const palette = [
    TRACK_COLORS.recon,
    TRACK_COLORS.calib,
    TRACK_COLORS.analysis,
    TRACK_COLORS.design,
    TRACK_COLORS.success,
    TRACK_COLORS.warn,
    TRACK_COLORS.error,
  ];

  // draw each series
  names.forEach((name, idx) => {
    const series = (seriesMap[name] || []).filter(p => isFinite(p.t) && isFinite(p.v));
    if (series.length === 0) return;
    const color = palette[idx % palette.length];
    ctx.strokeStyle = color;
    ctx.lineWidth = 2;
    ctx.beginPath();
    for (let i = 0; i < series.length; i++) {
      const x = xToPx(series[i].t);
      const y = yToPx(series[i].v);
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
  });

  // legend top-right
  ctx.font = '11px Inter, sans-serif';
  ctx.textBaseline = 'middle';
  ctx.textAlign = 'left';
  const lx = w - padding.right + 8;
  let ly = padding.top + 4;
  names.forEach((name, idx) => {
    const color = palette[idx % palette.length];
    ctx.strokeStyle = color;
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(lx, ly);
    ctx.lineTo(lx + 18, ly);
    ctx.stroke();
    ctx.fillStyle = TRACK_COLORS.mute;
    ctx.fillText(String(name), lx + 22, ly);
    ly += 14;
  });

  // axis labels
  drawAxisLabel(ctx, xLabel, (padding.left + w - padding.right) / 2, h - 12);
  ctx.save();
  ctx.translate(12, (padding.top + h - padding.bottom) / 2);
  ctx.rotate(-Math.PI / 2);
  drawAxisLabel(ctx, yLabel, 0, 0);
  ctx.restore();

  ctx.restore();
}

// -----------------------------------------------------------------------------
// File upload
// -----------------------------------------------------------------------------

export async function uploadFile(file) {
  if (!file) throw new Error('file is required');
  const fd = new FormData();
  fd.append('file', file);
  const res = await fetch('/api/uploads', { method: 'POST', body: fd });
  let body;
  try { body = await res.json(); } catch (_) { body = null; }
  if (!res.ok) {
    const err = (body && body.error) ? body : { error_code: 'upload_failed', error_message: `HTTP ${res.status}` };
    return err;
  }
  return body || {};
}

// -----------------------------------------------------------------------------
// Track shell (spec §8) — three-column layout shared by every track view.
// -----------------------------------------------------------------------------
//
// Wrap a track view's content in the standard three-column layout.
//   rootEl : HTMLElement (the <main> slot)
//   opts   : { track: 'recon'|'calib'|'analysis'|'design',
//              title: string,             // aside heading
//              sections: [{ id, label, isExtra?: boolean }],
//              content: HTMLElement|DocumentFragment }
// Returns a cleanup function.
export function mountTrackShell(rootEl, opts) {
  if (!rootEl) return () => {};
  opts = opts || {};
  rootEl.innerHTML = '';
  rootEl.dataset.track = opts.track || '';

  const shell = document.createElement('div');
  shell.className = 'track-shell';

  // -- Aside --
  const aside = document.createElement('aside');
  aside.className = 'track-aside';
  aside.innerHTML = `
    <h2 class="track-aside__title">${escapeHtml(opts.title || '')}</h2>
    <nav class="track-aside__list" role="navigation" aria-label="${escapeHtml(opts.title || 'Track')} sections"></nav>
  `;
  const asideList = aside.querySelector('.track-aside__list');
  (opts.sections || []).forEach((s) => {
    const a = document.createElement('a');
    a.className = 'track-aside__link';
    a.href = `#${s.id}`;
    a.textContent = s.label;
    a.dataset.section = s.id;
    if (s.isExtra) a.classList.add('is-extra');
    asideList.appendChild(a);
  });
  const back = document.createElement('a');
  back.className = 'track-aside__link track-aside__back';
  back.href = '#home';
  back.textContent = '← Home';
  aside.appendChild(back);

  // -- Main --
  const main = document.createElement('main');
  main.className = 'track-main';
  if (opts.content) main.appendChild(opts.content);

  // -- TOC --
  const toc = document.createElement('aside');
  toc.className = 'track-toc';
  toc.innerHTML = `
    <div class="track-toc__title">ON THIS PAGE</div>
    <nav class="track-toc__list" role="navigation" aria-label="On this page"></nav>
  `;
  const tocList = toc.querySelector('.track-toc__list');
  const h2s = Array.from(main.querySelectorAll('h2'));
  h2s.forEach((h2) => {
    if (!h2.id) h2.id = slugify(h2.textContent);
    const link = document.createElement('a');
    link.className = 'track-toc__link';
    link.href = `#${h2.id}`;
    link.textContent = h2.textContent;
    link.dataset.section = h2.id;
    tocList.appendChild(link);
  });

  shell.appendChild(aside);
  shell.appendChild(main);
  shell.appendChild(toc);
  rootEl.appendChild(shell);

  // -- Intersection observer: highlight current section --
  const headings = main.querySelectorAll('[data-anchor], h2[id]');
  const links = new Map();
  asideList.querySelectorAll('a').forEach((a) => links.set(a.dataset.section, a));
  tocList.querySelectorAll('a').forEach((a) => links.set(a.dataset.section, a));
  const observer = (typeof IntersectionObserver === 'undefined')
    ? { observe() {}, disconnect() {} }
    : new IntersectionObserver((entries) => {
      entries.forEach((entry) => {
        if (entry.isIntersecting) {
          const id = entry.target.id;
          links.forEach((a) => a.classList?.toggle('is-current', a.dataset.section === id));
        }
      });
    }, { rootMargin: '-30% 0px -60% 0px', threshold: 0 });
  headings.forEach((h) => observer.observe(h));

  // -- Smooth scroll --
  const onClick = (e) => {
    const a = e.target.closest('a[href="#home"]');
    if (a) {
      e.preventDefault();
      showHome();
      return;
    }
    const a2 = e.target.closest('a[href^="#"]');
    if (!a2) return;
    const id = a2.getAttribute('href').slice(1);
    const target = document.getElementById(id);
    if (target) {
      e.preventDefault();
      target.scrollIntoView({ behavior: 'smooth', block: 'start' });
    }
  };
  shell.addEventListener('click', onClick);
  return () => {
    observer.disconnect();
    shell.removeEventListener('click', onClick);
  };
}

// --- small helpers used by mountTrackShell ---
function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, (c) =>
    ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c])
  );
}
function slugify(s) {
  return String(s).toLowerCase().trim().replace(/[^a-z0-9]+/g, '-').replace(/(^-|-$)/g, '');
}

// -----------------------------------------------------------------------------
// Code block (spec §9) — black panel with lang label + copy button.
// -----------------------------------------------------------------------------
//
// lang : string label rendered top-left (e.g. 'matlab')
// code : code text body (no escaping needed — passed through escapeHtml)
// Returns an HTMLPreElement that the caller appends into a section.
export function codeBlock(lang, code) {
  const pre = document.createElement('pre');
  pre.className = 'code-block';
  pre.innerHTML = `
    <span class="code-block__lang">${escapeHtml(lang || '')}</span>
    <button type="button" class="code-block__copy" aria-label="Copy code">⧉</button>
    <code>${escapeHtml(code || '')}</code>
  `;
  const copyBtn = pre.querySelector('.code-block__copy');
  copyBtn.addEventListener('click', async () => {
    try {
      await navigator.clipboard.writeText(String(code || ''));
      const prev = copyBtn.textContent;
      copyBtn.textContent = '✓';
      setTimeout(() => { copyBtn.textContent = prev; }, 1200);
    } catch (_) { /* clipboard rejected */ }
  });
  return pre;
}

// -----------------------------------------------------------------------------
// Doc section builder (Phase 6 helper) — small factory for the
// <section data-anchor="..."> blocks used by every track page.
// -----------------------------------------------------------------------------
//
// opts : { anchor, title, paragraphs?: string[], code?: { lang, body },
//          runButton?: { label, onClick } }
// Returns a HTMLElement ready to append into the main content column.
export function docSection(opts) {
  const sec = document.createElement('section');
  sec.className = 'doc-section';
  sec.setAttribute('data-anchor', opts.anchor);

  const h2 = document.createElement('h2');
  h2.id = opts.anchor;
  h2.textContent = opts.title || '';
  sec.appendChild(h2);

  for (const p of (opts.paragraphs || [])) {
    const pe = document.createElement('p');
    pe.textContent = p;
    sec.appendChild(pe);
  }

  if (opts.code) {
    sec.appendChild(codeBlock(opts.code.lang, opts.code.body));
  }

  if (opts.runButton) {
    const btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'btn btn-primary doc-section__run';
    btn.setAttribute('data-track', opts.runButton.track || '');
    btn.textContent = opts.runButton.label || 'Run in MATLAB →';
    btn.addEventListener('click', opts.runButton.onClick || (() => {}));
    sec.appendChild(btn);
  }

  return sec;
}

// -----------------------------------------------------------------------------
// Doc page H1 header (Phase 6 helper) — used at the top of every track body.
// -----------------------------------------------------------------------------
export function docHeader(title, tagline) {
  const head = document.createElement('header');
  head.className = 'doc-header';
  const h1 = document.createElement('h1');
  h1.textContent = title;
  head.appendChild(h1);
  if (tagline) {
    const tag = document.createElement('p');
    tag.className = 'doc-header__tag';
    tag.textContent = tagline;
    head.appendChild(tag);
  }
  return head;
}
