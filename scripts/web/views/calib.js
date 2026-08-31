import {
  createParameterManager,
  uploadFile,
} from './shared.js';

export const CALIB_METHOD_ORDER = [
  'sluice',
  'kcatRepo',
  'sensitivity',
  'bayesian',
  'presto',
  'gauks',
];

const CALIB_METHODS = {
  sluice: {
    title: 'Sluice',
    desc: 'Apply the Sluice structure before downstream parameter calibration.',
  },
  kcatRepo: {
    title: 'KcatRepo',
    desc: 'Create a shared kcat repository with the initial group for later method outputs.',
  },
  sensitivity: {
    title: 'Sensitivity',
    desc: 'Tune around glucose uptake and matched growth-rate measurements.',
  },
  bayesian: {
    title: 'Bayesian',
    desc: 'Run ABC-style Bayesian calibration on selected scenarios.',
  },
  presto: {
    title: 'PRESTO',
    desc: 'Fit proteomics, growth, and total protein constraints.',
  },
  gauks: {
    title: 'GAUKS',
    desc: 'Use GAUKS with unconstrained maximum growth data.',
  },
};

const CALIB_DEFAULTS = {
  exRxns: 'EX_glc__D_e',
  glucoseExchange: 'EX_glc__D_e',
  targetGrowth: 0.3,
  factor: 1.1,
  multiCondition: false,
  maxIter: 200,
  proc: 4,
  numPerGen: 100,
  rejectNum: 0.05,
  runGauksAfterBayesian: false,
  biomassReaction: 'biomass',
};

const METHOD_REQUIREMENTS = {
  sluice: ['exchangeRxns'],
  kcatRepo: ['ecModelId'],
  sensitivity: ['glucoseUptake', 'growthRate', 'glucoseExchange', 'targetGrowth', 'factor', 'multiCondition'],
  bayesian: ['scenarioData', 'maxIter', 'proc', 'numPerGen', 'rejectNum'],
  presto: ['proteomics', 'growth', 'totalProtein'],
  gauks: ['gemModelId', 'biomassReaction', 'unconstrainedMaxGrowth'],
};

export function getCalibRunMethods(selected) {
  const set = selected instanceof Set ? selected : new Set(selected || []);
  return CALIB_METHOD_ORDER.filter((method) => set.has(method));
}

export function methodRequiresData(method) {
  return [...(METHOD_REQUIREMENTS[method] || [])];
}

export function createCalibState(opts = {}) {
  return {
    project: opts.project || null,
    ecModelId: opts.ecModelId || '',
    gemModelId: opts.gemModelId || '',
    selectedMethods: new Set(opts.selectedMethods || ['sluice', 'kcatRepo', 'sensitivity']),
    params: createParameterManager({
      defaults: CALIB_DEFAULTS,
      projectParams: opts.projectParams || {},
      runParams: opts.runParams || {},
    }),
    data: {
      sensitivity: { glucoseUptake: '', growthRate: '' },
      bayesian: { scenarioData: '' },
      presto: { proteomics: '', growth: '', totalProtein: '' },
      gauks: { unconstrainedMaxGrowth: '' },
      ...(opts.data || {}),
    },
    log: [],
  };
}

export function buildCalibPayload(state) {
  return {
    projectId: state.project?.projectId || state.project?.project_id || '',
    ecModelId: state.ecModelId,
    gemModelId: state.gemModelId,
    methods: getCalibRunMethods(state.selectedMethods),
    params: state.params.snapshot(),
    data: JSON.parse(JSON.stringify(state.data || {})),
  };
}

function setText(root, selector, text) {
  const el = root.querySelector(selector);
  if (el) el.textContent = text;
}

function appendLog(root, state, line) {
  state.log.push(line);
  const log = root.querySelector('[data-role="calib-log"]');
  if (log) log.textContent = state.log.join('\n');
}

function syncCalibForm(root, state) {
  const form = root.querySelector('[data-role="calib-params"]');
  if (!form) return;
  const field = (name) => form.elements.namedItem(name);
  state.ecModelId = root.querySelector('[data-role="ec-model-id"]')?.value.trim() || state.ecModelId;
  state.gemModelId = field('gemModelId')?.value.trim() || '';
  state.params.merge({
    exRxns: field('exRxns')?.value || '',
    glucoseExchange: field('glucoseExchange')?.value || 'EX_glc__D_e',
    targetGrowth: Number(field('targetGrowth')?.value || 0),
    factor: Number(field('factor')?.value || 1),
    multiCondition: Boolean(field('multiCondition')?.checked),
    maxIter: Number(field('maxIter')?.value || 200),
    proc: Number(field('proc')?.value || 4),
    numPerGen: Number(field('numPerGen')?.value || 100),
    rejectNum: Number(field('rejectNum')?.value || 0.05),
    runGauksAfterBayesian: Boolean(field('runGauksAfterBayesian')?.checked),
    biomassReaction: field('biomassReaction')?.value || 'biomass',
  });
  setText(root, '[data-role="calib-queue"]', getCalibRunMethods(state.selectedMethods).join(' -> ') || 'Select at least one method.');
  setText(root, '[data-role="calib-param-preview"]', JSON.stringify(state.params.snapshot(), null, 2));
}

function hydrateCalibForm(root, state) {
  const params = state.params.snapshot();
  const form = root.querySelector('[data-role="calib-params"]');
  if (!form) return;
  const setValue = (name, value) => {
    const field = form.elements.namedItem(name);
    if (field && value != null) field.value = value;
  };
  setValue('exRxns', params.exRxns);
  setValue('glucoseExchange', params.glucoseExchange);
  setValue('targetGrowth', params.targetGrowth);
  setValue('factor', params.factor);
  setValue('maxIter', params.maxIter);
  setValue('proc', params.proc);
  setValue('numPerGen', params.numPerGen);
  setValue('rejectNum', params.rejectNum);
  setValue('gemModelId', state.gemModelId);
  setValue('biomassReaction', params.biomassReaction);
  const multi = form.elements.namedItem('multiCondition');
  if (multi) multi.checked = Boolean(params.multiCondition);
  const gauks = form.elements.namedItem('runGauksAfterBayesian');
  if (gauks) gauks.checked = Boolean(params.runGauksAfterBayesian);
  const ecInput = root.querySelector('[data-role="ec-model-id"]');
  if (ecInput && state.ecModelId) ecInput.value = state.ecModelId;
  setText(root, '[data-role="ec-model-status"]', state.ecModelId ? `Using project ecGEM ${state.ecModelId}` : 'No ecGEM loaded.');
}

async function handleDataUpload(root, state, input) {
  const file = input.files?.[0];
  if (!file) return;
  const method = input.dataset.method;
  const field = input.dataset.field;
  setText(root, `[data-status-for="${method}.${field}"]`, `Uploading ${file.name}...`);
  try {
    const up = await uploadFile(file);
    if (up.error_code) throw new Error(up.error_message || up.error_code);
    if (!state.data[method]) state.data[method] = {};
    state.data[method][field] = up.path || up.uploadId || up.name || file.name;
    setText(root, `[data-status-for="${method}.${field}"]`, up.name || file.name);
  } catch (err) {
    setText(root, `[data-status-for="${method}.${field}"]`, err.message || String(err));
  }
}

async function runCalibration(root, state) {
  syncCalibForm(root, state);
  const payload = buildCalibPayload(state);
  if (!payload.ecModelId) {
    appendLog(root, state, 'calibration: err_model_format');
    return;
  }
  appendLog(root, state, `calibration: ${payload.methods.join(' -> ')} queued`);
  try {
    const res = await fetch('/api/calib/run', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    const body = await res.json().catch(() => ({}));
    if (!res.ok || body.ok === false) throw new Error(body.error || body.detail || `HTTP ${res.status}`);
    appendLog(root, state, 'calibration: done');
  } catch (err) {
    appendLog(root, state, `calibration: ${err.message || err}`);
  }
}

function dataField(method, field, label, accept = '.csv,.tsv,.txt,.mat') {
  return `
    <label class="data-upload">
      <span>${label}</span>
      <input type="file" data-role="method-data" data-method="${method}" data-field="${field}" accept="${accept}" />
      <small data-status-for="${method}.${field}">No file</small>
    </label>
  `;
}

export function mountCalib(rootEl, opts = {}) {
  if (!rootEl) return () => {};
  const state = createCalibState(opts);
  rootEl.dataset.track = 'calib';
  rootEl.innerHTML = `
    <section class="track-page">
      <header class="track-header" data-track="calib">
        <h2>Calibration</h2>
        <p>Run selected calibration methods independently, with shared run parameters and method-specific data checks.</p>
      </header>

      <div class="module-workbench module-workbench--calib">
        <aside class="workbench-panel">
          <h3>Model and parameters</h3>
          <label class="label">ecGEM ID</label>
          <input class="input" data-role="ec-model-id" placeholder="ecGEM ID from upload or Reconstruction" />
          <label class="label">Upload ecGEM</label>
          <input type="file" data-role="ec-model-upload" accept=".mat,.xml,.json,.yml,.yaml" />
          <p class="help-text" data-role="ec-model-status">No ecGEM loaded.</p>

          <form data-role="calib-params">
            <label class="label">Sluice exchange reactions</label>
            <input class="input" name="exRxns" value="EX_glc__D_e" />
            <label class="label">Glucose exchange reaction</label>
            <input class="input" name="glucoseExchange" value="EX_glc__D_e" />
            <div class="workbench-grid-2">
              <label><span class="label">Target growth</span><input class="input" name="targetGrowth" type="number" step="0.01" value="0.3" /></label>
              <label><span class="label">Sensitivity factor</span><input class="input" name="factor" type="number" step="0.01" value="1.1" /></label>
            </div>
            <label class="check-line"><input type="checkbox" name="multiCondition" /> Multi-condition Sensitivity</label>
            <div class="workbench-grid-2">
              <label><span class="label">Bayesian max iter</span><input class="input" name="maxIter" type="number" value="200" /></label>
              <label><span class="label">Processes</span><input class="input" name="proc" type="number" value="4" /></label>
            </div>
            <div class="workbench-grid-2">
              <label><span class="label">Samples per generation</span><input class="input" name="numPerGen" type="number" value="100" /></label>
              <label><span class="label">Reject ratio</span><input class="input" name="rejectNum" type="number" step="0.01" value="0.05" /></label>
            </div>
            <label class="check-line"><input type="checkbox" name="runGauksAfterBayesian" /> Run GAUKS after Bayesian</label>
            <label class="label">GEM ID for GAUKS</label>
            <input class="input" name="gemModelId" placeholder="unconstrained GEM ID" />
            <label class="label">Biomass reaction for GAUKS</label>
            <input class="input" name="biomassReaction" value="biomass" />
          </form>
          <pre class="param-preview" data-role="calib-param-preview"></pre>
        </aside>

        <main class="workbench-main">
          <section class="calib-method-grid" aria-label="Calibration methods">
            ${CALIB_METHOD_ORDER.map((method) => `
              <article class="calib-method-card ${state.selectedMethods.has(method) ? 'is-selected' : ''}" data-method="${method}">
                <label>
                  <input type="checkbox" data-role="method-toggle" value="${method}" ${state.selectedMethods.has(method) ? 'checked' : ''} />
                  <strong>${CALIB_METHODS[method].title}</strong>
                </label>
                <p>${CALIB_METHODS[method].desc}</p>
                <small>Requires: ${methodRequiresData(method).join(', ')}</small>
              </article>
            `).join('')}
          </section>

          <section class="workbench-card">
            <h3>Method-specific data</h3>
            <p>Sensitivity requires glucose uptake and matched growth-rate data, plus the parameters above. GAUKS requires unconstrained maximum growth data.</p>
            <div class="data-grid">
              ${dataField('sensitivity', 'glucoseUptake', 'Sensitivity glucose uptake')}
              ${dataField('sensitivity', 'growthRate', 'Sensitivity growth-rate')}
              ${dataField('bayesian', 'scenarioData', 'Bayesian scenario data')}
              ${dataField('presto', 'proteomics', 'PRESTO proteomics')}
              ${dataField('presto', 'growth', 'PRESTO growth')}
              ${dataField('presto', 'totalProtein', 'PRESTO total protein')}
              ${dataField('gauks', 'unconstrainedMaxGrowth', 'GAUKS unconstrained maximum growth')}
            </div>
          </section>

          <section class="workbench-card">
            <h3>Run queue</h3>
            <p data-role="calib-queue"></p>
            <button type="button" class="btn btn-primary" data-track="calib" data-role="run-calib">Run calibration</button>
            <pre class="run-log" data-role="calib-log">Ready.</pre>
          </section>
        </main>
      </div>
    </section>
  `;

  hydrateCalibForm(rootEl, state);
  syncCalibForm(rootEl, state);

  const paramsForm = rootEl.querySelector('[data-role="calib-params"]');
  paramsForm?.addEventListener('input', () => syncCalibForm(rootEl, state));
  paramsForm?.addEventListener('change', () => syncCalibForm(rootEl, state));

  rootEl.querySelector('[data-role="ec-model-id"]')?.addEventListener('input', (event) => {
    state.ecModelId = event.target.value.trim();
    syncCalibForm(rootEl, state);
  });

  rootEl.querySelector('[data-role="ec-model-upload"]')?.addEventListener('change', async (event) => {
    const file = event.target.files?.[0];
    if (!file) return;
    setText(rootEl, '[data-role="ec-model-status"]', `Uploading ${file.name}...`);
    try {
      const up = await uploadFile(file);
      if (up.error_code) throw new Error(up.error_message || up.error_code);
      state.ecModelId = up.modelId || up.uploadId || '';
      const input = rootEl.querySelector('[data-role="ec-model-id"]');
      if (input) input.value = state.ecModelId;
      setText(rootEl, '[data-role="ec-model-status"]', `Loaded ${up.name || file.name}`);
    } catch (err) {
      setText(rootEl, '[data-role="ec-model-status"]', err.message || String(err));
    }
  });

  rootEl.querySelectorAll('[data-role="method-toggle"]').forEach((input) => {
    input.addEventListener('change', () => {
      if (input.checked) state.selectedMethods.add(input.value);
      else state.selectedMethods.delete(input.value);
      input.closest('.calib-method-card')?.classList.toggle('is-selected', input.checked);
      syncCalibForm(rootEl, state);
    });
  });

  rootEl.querySelectorAll('[data-role="method-data"]').forEach((input) => {
    input.addEventListener('change', () => handleDataUpload(rootEl, state, input));
  });

  rootEl.querySelector('[data-role="run-calib"]')?.addEventListener('click', () => runCalibration(rootEl, state));

  return function teardown() {
    rootEl.innerHTML = '';
  };
}
