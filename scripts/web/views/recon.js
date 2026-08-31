import {
  createParameterManager,
  deepLearningKcatHtml,
  uploadFile,
} from './shared.js';

const RECON_DEFAULTS = {
  organism: 'custom',
  topology: 'integrated',
  topologies: ['integrated'],
  sigma: 0.5,
  f: 0.46,
  runComplexAnnotation: true,
  runEcAnnotation: true,
  runMetaNetXIntegration: true,
  annotationStages: ['ec', 'metabolite'],
  metaboliteSources: ['A', 'B', 'C'],
  deepLearningModels: ['DLKcat', 'UniKP', 'CatPred'],
  kcatReferenceTopology: '',
  kcatPredictionModel: 'CatPred',
  customKcatRxnNameType: '',
  useCustomKcatFile: false,
  useLoggedMedian: true,
  carbonSource: 'glucose',
  biomassReaction: 'biomass',
  InitialModel: '',
  modeltype: 'Tradition',
  Ptot: 0.56,
  org_name: '',
  c_source: '',
  bioRxn: '',
  taxonomicID: '',
  'uniprot.type': 'proteome',
  'uniprot.ID': '',
  'uniprot.geneIDfield': 'gene_oln',
  'uniprot.reviewed': false,
  'PRESTO.runParallel': true,
  'PRESTO.ncpu': 8,
  'PRESTO.nIter': 50,
  'PRESTO.epsilon': 1e5,
  'PRESTO.lambda': 1e-5,
  'PRESTO.theta': 0.6,
};

const NUMERIC_PARAMS = new Set([
  'sigma', 'Ptot', 'f', 'figResolutionDPI', 'medianThreshold',
  'PRESTO.ncpu', 'PRESTO.nIter', 'PRESTO.epsilon', 'PRESTO.lambda', 'PRESTO.theta',
]);
const INTEGER_PARAMS = new Set(['figResolutionDPI', 'PRESTO.ncpu', 'PRESTO.nIter']);
const BOOLEAN_PARAMS = new Set(['useCustomKcatFile', 'useLoggedMedian', 'runComplexAnnotation', 'runEcAnnotation', 'runMetaNetXIntegration', 'uniprot.reviewed', 'PRESTO.runParallel']);
const MULTI_PARAMS = new Set(['annotationStages', 'deepLearningModels', 'topologies', 'metaboliteSources']);

const PARAM_LABELS = {
  projectName: 'Project name',
  InitialModel: 'Initial model',
  modeltype: 'MATLAB model type',
  org_name: 'Organism',
  taxonomicID: 'Taxonomic ID',
  c_source: 'Carbon source reaction',
  bioRxn: 'Biomass reaction',
  runComplexAnnotation: 'Use complex information',
  runEcAnnotation: 'Annotate EC numbers',
  runMetaNetXIntegration: 'Integrate MetaNetX identifiers',
  kcatReferenceTopology: 'kcat reference topology',
  kcatPredictionModel: 'kcat prediction source',
  customKcatRxnNameType: 'Custom kcat reaction naming',
  useLoggedMedian: 'Use logged median fallback',
  biomassReaction: 'Biomass reaction',
  carbonSource: 'Carbon source',
  'uniprot.type': 'UniProt lookup type',
  'uniprot.ID': 'UniProt ID',
  'uniprot.geneIDfield': 'Gene ID field',
  'uniprot.reviewed': 'Reviewed UniProt entries only',
  'PRESTO.runParallel': 'Run PRESTO in parallel',
  'PRESTO.ncpu': 'PRESTO CPU cores',
  'PRESTO.nIter': 'PRESTO iterations',
  'PRESTO.epsilon': 'PRESTO epsilon',
  'PRESTO.lambda': 'PRESTO lambda',
  'PRESTO.theta': 'PRESTO theta',
};

export const RECON_STEPS = [
  { id: 'init', titleKey: 'Project parameters', action: 'init', nextValid: () => true },
  { id: 'load', titleKey: 'Upload model', action: 'load', nextValid: (s) => Boolean(s.modelId) },
  { id: 'convert', titleKey: 'S-matrix topology', action: 'convert', nextValid: (s) => Boolean(s.modelId) },
  { id: 'annotate', titleKey: 'Annotation', action: 'annotate', nextValid: (s) => hasEcModel(s) },
  { id: 'deepLearningKcat', titleKey: 'Deep Learning kcat', action: 'deepLearningKcat', nextValid: (s) => hasEcModel(s) },
  { id: 'compare', titleKey: 'Compare kcat', action: 'compare', nextValid: (s) => hasEcModel(s) },
  { id: 'merge', titleKey: 'Merge kcat', action: 'merge', nextValid: (s) => hasEcModel(s) },
  { id: 'growth', titleKey: 'Growth check', action: 'growth', nextValid: (s) => hasEcModel(s) },
];

const PROJECT_STEPS = [
  {
    id: 'initialize',
    legacyAction: 'init',
    title: 'Initialize project',
    functionName: 'InitializeECOMAPproject',
    description: 'Configure the shared parameter manager for this project. Values are written back to the project ParameterManagement.m file.',
    params: [
      'InitialModel', 'modeltype', 'sigma', 'Ptot', 'f',
      'org_name', 'uniprot.type', 'uniprot.ID', 'uniprot.geneIDfield', 'uniprot.reviewed',
      'taxonomicID', 'c_source', 'bioRxn',
      'PRESTO.runParallel', 'PRESTO.ncpu', 'PRESTO.nIter', 'PRESTO.epsilon', 'PRESTO.lambda', 'PRESTO.theta',
    ],
  },
  {
    id: 'loadGem',
    legacyAction: 'load',
    title: 'Load GEM',
    functionName: 'loadModel',
    description: 'Load the initial GEM from the project models directory.',
    params: ['InitialModel', 'modeltype'],
  },
  {
    id: 'convert',
    legacyAction: 'convert',
    title: 'Convert to ecGEMs',
    functionName: 'convertecModel',
    description: 'Build basic, isozyme, and integrated ecGEM variants in the tutorial order.',
    params: ['topologies', 'sigma', 'Ptot', 'f'],
  },
  {
    id: 'annotate',
    legacyAction: 'annotate',
    title: 'Apply annotations',
    functionName: 'getComplexdata / getECnumber / getMetinfo / addMetMetaNetXID',
    description: 'Apply complex, EC number, metabolite, and MetaNetX annotations.',
    params: ['runComplexAnnotation', 'runEcAnnotation', 'metaboliteSources', 'runMetaNetXIntegration'],
  },
  {
    id: 'predictKcat',
    legacyAction: 'deepLearningKcat',
    title: 'Predict Deep Learning kcat',
    functionName: 'writeInputFile / ExecutePrediction',
    description: 'Generate prediction input tables from the selected reference topology, then execute the selected Deep Learning kcat predictors.',
    params: ['kcatReferenceTopology', 'deepLearningModels'],
  },
  {
    id: 'compareKcat',
    legacyAction: 'compare',
    title: 'Compare kcat predictions',
    functionName: 'BuildKcatMatches / AnalyzeKcatMatches',
    description: 'Compare predicted kcat values against database values and save benchmark figures.',
    params: ['kcatReferenceTopology', 'deepLearningModels', 'figFormat', 'figResolutionDPI'],
  },
  {
    id: 'mergeKcat',
    legacyAction: 'merge',
    title: 'Integrate and assign kcat values',
    functionName: 'getPrediction / completeKcatMatch / mergeKcats / selectKcatValue',
    description: 'Use the same reference topology as prediction, remap only when multiple topology naming spaces are present, then assign final kcat values.',
    params: ['kcatReferenceTopology', 'kcatPredictionModel', 'customKcatRxnNameType', 'useCustomKcatFile', 'medianThreshold', 'useLoggedMedian'],
  },
  {
    id: 'growthSave',
    legacyAction: 'growth',
    title: 'Growth validation and save models',
    functionName: 'solveLP / saveModel',
    description: 'Validate growth under protein constraints and save final ecGEM files.',
    params: ['carbonSource', 'c_source', 'bioRxn', 'biomassReaction'],
  },
];

function hasEcModel(state) {
  return Boolean(Object.keys(state.ecModelIds || {}).length || state.ecModelId);
}

function normalizeTopologyParams(params = {}) {
  const next = { ...(params || {}) };
  if (!Array.isArray(next.topologies) && next.topology) {
    next.topologies = next.topology === 'all' ? ['basic', 'isozyme', 'integrated'] : [next.topology];
  }
  if (Array.isArray(next.topologies) && !next.topology) {
    next.topology = legacyTopologyValue(next.topologies);
  }
  if (next.annotationOptions && typeof next.annotationOptions === 'object') {
    if (typeof next.runComplexAnnotation !== 'boolean' && typeof next.annotationOptions.runComplexAnnotation === 'boolean') {
      next.runComplexAnnotation = next.annotationOptions.runComplexAnnotation;
    }
    if (typeof next.runEcAnnotation !== 'boolean' && typeof next.annotationOptions.runEcAnnotation === 'boolean') {
      next.runEcAnnotation = next.annotationOptions.runEcAnnotation;
    }
    if (typeof next.runMetaNetXIntegration !== 'boolean' && typeof next.annotationOptions.runMetaNetXIntegration === 'boolean') {
      next.runMetaNetXIntegration = next.annotationOptions.runMetaNetXIntegration;
    }
    if (!Array.isArray(next.metaboliteSources) && Array.isArray(next.annotationOptions.metaboliteSources)) {
      next.metaboliteSources = next.annotationOptions.metaboliteSources;
    }
  }
  if (Array.isArray(next.annotationStages)) {
    next.annotationStages = selectedAnnotationStages(next);
  }
  return next;
}

export function createReconState(opts = {}) {
  return {
    step: 0,
    project: opts.project || null,
    modelId: opts.modelId || '',
    ecModelId: opts.ecModelId || '',
    ecModelIds: { ...(opts.ecModelIds || {}) },
    params: createParameterManager({
      defaults: RECON_DEFAULTS,
      projectParams: normalizeTopologyParams(opts.projectParams || {}),
      runParams: normalizeTopologyParams(opts.runParams || {}),
    }),
    uploads: {},
    log: [],
    projectState: opts.projectState || null,
    projectFiles: opts.projectFiles || [],
    activeProjectStep: opts.activeProjectStep || PROJECT_STEPS[0].id,
  };
}

export function canAdvance(state, idx) {
  const step = RECON_STEPS[idx];
  if (!step) return false;
  if (step.id === 'convert' && !selectedTopologies(state.params.snapshot()).length) return false;
  return step.nextValid(state);
}

export function validateStep(state, idx) {
  if (canAdvance(state, idx)) return null;
  const id = RECON_STEPS[idx]?.id;
  if (id === 'load') return 'err_model_format';
  if (id === 'convert') return 'err_param_invalid';
  if (id === 'annotate' || id === 'deepLearningKcat' || id === 'compare' || id === 'merge' || id === 'growth') {
    return 'err_param_invalid';
  }
  return 'err_init_fail';
}

export function buildReconPayload(state, action) {
  const params = state.params.snapshot();
  const kcatReferenceTopology = selectedKcatReferenceTopology(params, state);
  if (kcatReferenceTopology) {
    params.kcatReferenceTopology = kcatReferenceTopology;
    if (!params.customKcatRxnNameType) params.customKcatRxnNameType = kcatReferenceTopology;
  }
  params.annotationOptions = {
    runComplexAnnotation: Boolean(params.runComplexAnnotation),
    runEcAnnotation: Boolean(params.runEcAnnotation),
    metaboliteSources: selectedMetaboliteSources(params),
    runMetaNetXIntegration: Boolean(params.runMetaNetXIntegration),
  };
  params.annotationStages = annotationStagesFromOptions(params);
  const referenceEcModelId = ecModelIdForTopology(state, kcatReferenceTopology);
  const ecModelId = ['deepLearningKcat', 'compare'].includes(action)
    ? referenceEcModelId || state.ecModelId || firstEcModelId(state)
    : state.ecModelId || firstEcModelId(state);
  return {
    action,
    projectId: state.project?.projectId || state.project?.project_id || '',
    modelId: state.modelId,
    ecModelId,
    ecModelIds: { ...(state.ecModelIds || {}) },
    params,
  };
}

function firstEcModelId(state) {
  const ids = Object.values(state.ecModelIds || {}).filter(Boolean);
  return ids[0] || '';
}

function ecModelIdForTopology(state, topology) {
  if (!topology) return '';
  return (state.ecModelIds || {})[topology] || '';
}

function setText(root, selector, text) {
  const el = root.querySelector(selector);
  if (el) el.textContent = text;
}

function appendLog(root, state, line) {
  state.log.push(line);
  renderRunEvents(root, state);
}

function runEventTone(line) {
  const text = String(line);
  if (text.includes(': done') || text.includes('saved') || text.includes('found')) return 'done';
  if (text.includes(': queued') || text.includes('saving')) return 'running';
  if (text.includes('err_') || text.includes('not found') || text.includes('HTTP')) return 'error';
  return 'info';
}

function runEventRows(log) {
  if (!log.length) {
    return [{ line: 'Ready to run the selected step.', label: 'Ready' }];
  }
  let eventNumber = 0;
  const activeEvents = new Map();
  return log.map((line) => {
    const text = String(line);
    const action = text.split(':')[0] || 'event';
    const isStart = text.includes(': queued') || text.includes('saving parameter manager');
    const isResult = text.includes(': done')
      || text.includes('saved')
      || text.includes('err_')
      || text.includes('not found')
      || text.includes('HTTP')
      || (activeEvents.has(action) && !isStart);
    if (isStart || !activeEvents.has(action)) {
      eventNumber += 1;
      activeEvents.set(action, eventNumber);
    }
    const currentEvent = activeEvents.get(action);
    if (isResult) {
      activeEvents.delete(action);
      return { line, label: `Event ${currentEvent}_result` };
    }
    return { line, label: `Event ${currentEvent}` };
  });
}

function renderRunEvents(root, state) {
  const list = root.querySelector('[data-role="run-events"]');
  if (!list) return;
  const rows = runEventRows(state.log);
  list.innerHTML = rows.map(({ line, label }) => {
    const tone = runEventTone(line);
    return `<li class="run-event run-event--${tone}">
      <span class="run-event__dot" aria-hidden="true"></span>
      <div>
        <strong>${escapeHtml(label)}</strong>
        <span>${escapeHtml(line)}</span>
      </div>
    </li>`;
  }).join('');
}

function escapeHtml(value) {
  return String(value ?? '').replace(/[&<>"']/g, (ch) => ({
    '&': '&amp;',
    '<': '&lt;',
    '>': '&gt;',
    '"': '&quot;',
    "'": '&#39;',
  }[ch]));
}

function statusLabel(status) {
  if (status === 'completed') return 'Completed';
  if (status === 'ready') return 'Ready';
  if (status === 'locked') return 'Locked';
  return 'Needs review';
}

function activeStepFromProjectState(projectState) {
  const steps = projectState?.steps || [];
  return steps.find((step) => step.status === 'ready')?.id
    || steps.find((step) => step.status !== 'completed')?.id
    || steps[steps.length - 1]?.id
    || PROJECT_STEPS[0].id;
}

function projectStepState(state, stepId) {
  return (state.projectState?.steps || []).find((step) => step.id === stepId)
    || { id: stepId, status: 'ready', outputs: [] };
}

function projectFilesForStep(state, stepId) {
  const step = projectStepState(state, stepId);
  const outputs = new Set(step.outputs || []);
  return (state.projectFiles || []).filter((file) => outputs.has(file.relativePath));
}

function previewUrl(projectId, relativePath) {
  return `/api/projects/${encodeURIComponent(projectId)}/file/${encodeURIComponent(relativePath)}`;
}

function nestedParam(params, name) {
  const parts = String(name).split('.');
  let cursor = params;
  for (const part of parts) {
    if (!cursor || typeof cursor !== 'object' || !(part in cursor)) return undefined;
    cursor = cursor[part];
  }
  return cursor;
}

function paramValue(params, name) {
  if (Object.prototype.hasOwnProperty.call(params, name)) {
    const flat = params[name];
    if (flat !== '' && flat != null) return flat;
  }
  const nested = nestedParam(params, name);
  if (nested !== undefined) return nested;
  return params[name];
}

function selectedTopologies(params) {
  if (Array.isArray(params?.topologies)) return params.topologies.filter(Boolean);
  if (params?.topology === 'all') return ['basic', 'isozyme', 'integrated'];
  if (params?.topology) return [params.topology];
  return [];
}

function availableKcatTopologies(params, state = null) {
  const fromModels = Object.keys(state?.ecModelIds || {}).filter((key) => state.ecModelIds[key]);
  const values = fromModels.length ? fromModels : selectedTopologies(params);
  return ['integrated', 'isozyme', 'basic'].filter((topology) => values.includes(topology));
}

function selectedKcatReferenceTopology(params, state = null) {
  const available = availableKcatTopologies(params, state);
  const requested = params?.kcatReferenceTopology || params?.customKcatRxnNameType;
  if (requested && available.includes(requested)) return requested;
  return available[0] || '';
}

function selectedMetaboliteSources(params) {
  if (Array.isArray(params?.metaboliteSources)) return params.metaboliteSources.filter(Boolean);
  const options = params?.annotationOptions || {};
  if (Array.isArray(options.metaboliteSources)) return options.metaboliteSources.filter(Boolean);
  return ['A', 'B', 'C'];
}

function annotationStagesFromOptions(params) {
  const stages = [];
  if (params?.runEcAnnotation !== false) stages.push('ec');
  stages.push('metabolite');
  return stages;
}

function selectedAnnotationStages(params) {
  const stages = Array.isArray(params?.annotationStages) ? params.annotationStages : [];
  return stages.map((stage) => (stage === 'met' ? 'metabolite' : stage)).filter((stage) => stage !== 'complex');
}

function legacyTopologyValue(topologies) {
  const values = Array.isArray(topologies) ? topologies.filter(Boolean) : [];
  if (values.length === 1) return values[0];
  if (values.length === 3 && ['basic', 'isozyme', 'integrated'].every((id) => values.includes(id))) return 'all';
  return values[0] || '';
}

function renderParamControl(name, params, state = null) {
  const value = paramValue(params, name);
  if (name === 'topologies') {
    const selected = new Set(selectedTopologies(params));
    return `<fieldset class="option-stack recon-param-wide"><legend>Model structures</legend>
      ${[
        ['basic', 'basic'],
        ['isozyme', 'isozyme'],
        ['integrated', 'integrated'],
      ].map(([id, label]) => `<label><input type="checkbox" name="topologies" value="${id}"${selected.has(id) ? ' checked' : ''} /> ${label}</label>`).join('')}
    </fieldset>`;
  }
  if (name === 'modeltype') {
    const choices = ['Tradition', 'ECOMAP', 'sMOMENT', 'ECMpy', 'GECKO'];
    return `<label><span class="label">Model type</span><select class="input" name="modeltype">
      ${choices.map((choice) => `<option value="${choice}"${value === choice ? ' selected' : ''}>${choice}</option>`).join('')}
    </select></label>`;
  }
  if (name === 'uniprot.type') {
    const choices = ['proteome', 'taxonomy'];
    return `<label><span class="label">${PARAM_LABELS[name]}</span><select class="input" name="${escapeHtml(name)}">
      ${choices.map((choice) => `<option value="${choice}"${value === choice ? ' selected' : ''}>${choice}</option>`).join('')}
    </select></label>`;
  }
  if (name === 'metaboliteSources') {
    const selected = new Set(selectedMetaboliteSources(params));
    return `<fieldset class="option-stack recon-param-wide"><legend>Metabolite sources</legend>
      <p class="option-stack__hint">Leave all unchecked to use the local metInfo.tsv cache only.</p>
      ${[
        ['A', 'A - chem_prop.tsv'],
        ['B', 'B - ChEBI API'],
        ['C', 'C - PubChem'],
      ].map(([id, label]) => `<label><input type="checkbox" name="metaboliteSources" value="${id}"${selected.has(id) ? ' checked' : ''} /> ${label}</label>`).join('')}
    </fieldset>`;
  }
  if (name === 'deepLearningModels') {
    const selected = new Set(value || []);
    return `<fieldset class="option-stack recon-param-wide"><legend>${deepLearningKcatHtml()} models</legend>
      ${['DLKcat', 'UniKP', 'CatPred'].map((id) => `<label><input type="checkbox" name="deepLearningModels" value="${id}"${selected.has(id) ? ' checked' : ''} /> ${id}</label>`).join('')}
    </fieldset>`;
  }
  if (name === 'kcatReferenceTopology') {
    const choices = availableKcatTopologies(params, state);
    const selected = selectedKcatReferenceTopology(params, state);
    if (choices.length <= 1) {
      return `<div class="readonly-field recon-param-wide">
        <span class="label">${escapeHtml(PARAM_LABELS[name])}</span>
        <strong>${escapeHtml(selected || 'No ecGEM topology available')}</strong>
        <small>${escapeHtml(selected ? 'Single-topology runs use this model naming directly.' : 'Run Convert to ecGEMs first.')}</small>
      </div>`;
    }
    return `<label><span class="label">${escapeHtml(PARAM_LABELS[name])}</span><select class="input" name="kcatReferenceTopology">
      ${choices.map((choice) => `<option value="${choice}"${choice === selected ? ' selected' : ''}>${choice}</option>`).join('')}
    </select></label>`;
  }
  if (name === 'kcatPredictionModel') {
    const models = Array.isArray(params.deepLearningModels) && params.deepLearningModels.length
      ? params.deepLearningModels
      : ['DLKcat', 'UniKP', 'CatPred'];
    const selected = models.includes(value) ? value : models[models.length - 1];
    return `<label><span class="label">${escapeHtml(PARAM_LABELS[name])}</span><select class="input" name="kcatPredictionModel">
      ${models.map((choice) => `<option value="${choice}"${choice === selected ? ' selected' : ''}>${choice}</option>`).join('')}
    </select></label>`;
  }
  if (name === 'customKcatRxnNameType') {
    const choices = availableKcatTopologies(params, state);
    const selected = choices.includes(value) ? value : selectedKcatReferenceTopology(params, state);
    if (choices.length <= 1) {
      return `<div class="readonly-field recon-param-wide">
        <span class="label">${escapeHtml(PARAM_LABELS[name])}</span>
        <strong>${escapeHtml(selected || 'No ecGEM topology available')}</strong>
        <small>Custom kcat rows should use the same reaction naming as this topology.</small>
      </div>`;
    }
    return `<label><span class="label">${escapeHtml(PARAM_LABELS[name])}</span><select class="input" name="customKcatRxnNameType">
      ${choices.map((choice) => `<option value="${choice}"${choice === selected ? ' selected' : ''}>${choice}</option>`).join('')}
    </select></label>`;
  }
  if (BOOLEAN_PARAMS.has(name)) {
    return `<label class="check-line recon-param-wide"><input type="checkbox" name="${escapeHtml(name)}"${value ? ' checked' : ''} /> ${escapeHtml(PARAM_LABELS[name] || name)}</label>`;
  }
  if (NUMERIC_PARAMS.has(name)) {
    const fallback = name === 'figResolutionDPI' ? 300 : name === 'medianThreshold' ? 5 : value ?? '';
    const step = INTEGER_PARAMS.has(name) ? '1' : 'any';
    return `<label><span class="label">${escapeHtml(PARAM_LABELS[name] || name)}</span><input class="input" name="${escapeHtml(name)}" type="number" step="${step}" value="${escapeHtml(fallback)}" /></label>`;
  }
  if (name === 'figFormat') {
    const choices = ['png', 'svg', 'pdf'];
    return `<label><span class="label">Figure format</span><select class="input" name="figFormat">
      ${choices.map((choice) => `<option value="${choice}"${value === choice ? ' selected' : ''}>${choice}</option>`).join('')}
    </select></label>`;
  }
  return `<label><span class="label">${escapeHtml(PARAM_LABELS[name] || name)}</span><input class="input" name="${escapeHtml(name)}" value="${escapeHtml(value ?? '')}" /></label>`;
}

function renderStepOutputs(state, stepId) {
  const files = projectFilesForStep(state, stepId);
  if (!files.length) return '<p class="recent__empty">No files have been recorded for this step yet.</p>';
  return `<div class="recon-output-list">${files.map((file) => `
    <a class="recon-output-item" href="${previewUrl(state.project?.projectId || '', file.relativePath)}" target="_blank" rel="noreferrer">
      <span>${escapeHtml(file.name)}</span>
      <small>${escapeHtml(file.relativePath)}</small>
    </a>
  `).join('')}</div>`;
}

function renderProjectFiles(state) {
  const files = state.projectFiles || [];
  if (!files.length) return '<p class="recent__empty">No project files have been detected yet.</p>';
  return files.slice(0, 18).map((file) => `
    <article class="project-file-chip">
      ${file.previewable ? `<img src="${previewUrl(state.project?.projectId || '', file.relativePath)}" alt="${escapeHtml(file.name)}" />` : ''}
      <div>
        <strong>${escapeHtml(file.name)}</strong>
        <span>${escapeHtml(file.kind || 'file')} · ${escapeHtml(file.relativePath)}</span>
      </div>
    </article>
  `).join('');
}

function renderStepper(root, state) {
  const projectParams = state.params.snapshot();
  const activeSpec = PROJECT_STEPS.find((step) => step.id === state.activeProjectStep) || PROJECT_STEPS[0];
  const activeState = projectStepState(state, activeSpec.id);
  const projectId = state.project?.projectId || '';
  const stepButtons = PROJECT_STEPS.map((step, idx) => {
    const stepState = projectStepState(state, step.id);
    const active = step.id === activeSpec.id ? ' is-active' : '';
    return `<button type="button" class="recon-step-btn${active}" data-project-step="${step.id}" data-action="${step.legacyAction}">
      <span>${idx + 1}</span>
      <strong>${escapeHtml(step.title)}</strong>
      <small class="step-status step-status--${escapeHtml(stepState.status)}">${statusLabel(stepState.status)}</small>
    </button>`;
  }).join('');
  root.querySelector('[data-role="recon-stepper"]').innerHTML = stepButtons;
  root.querySelector('[data-role="active-step-panel"]').innerHTML = `
    <div class="active-step__head">
      <div>
        <p class="section-eyebrow">${escapeHtml(activeSpec.functionName)}</p>
        <h3>${escapeHtml(activeSpec.title)}</h3>
      </div>
      <span class="step-status step-status--${escapeHtml(activeState.status)}">${statusLabel(activeState.status)}</span>
    </div>
    <p>${escapeHtml(activeSpec.description)}</p>
    <form data-role="recon-params" class="recon-step-params">
      ${activeSpec.params.map((name) => renderParamControl(name, projectParams, state)).join('')}
    </form>
    <div class="active-step__actions">
      ${activeSpec.id === 'loadGem' ? '<button type="button" class="btn btn-secondary" data-role="find-project-model">Find in models</button>' : ''}
      <button type="button" class="btn btn-primary" data-role="run-active-step" data-action="${activeSpec.legacyAction}"${activeState.status === 'locked' ? ' disabled' : ''}>Run step</button>
      <button type="button" class="btn btn-secondary" data-role="next-project-step"${activeState.status !== 'completed' ? ' disabled' : ''}>Next</button>
    </div>
    <section class="step-output-card">
      <h4>Outputs from this step</h4>
      <p class="help-text">Files produced or used by the selected workflow step.</p>
      ${renderStepOutputs(state, activeSpec.id)}
    </section>
  `;
  root.querySelector('[data-role="project-file-list"]').innerHTML = renderProjectFiles(state);
  root.querySelector('[data-role="project-summary-line"]').textContent = projectId
    ? `${projectId} · ${state.projectFiles.length} files detected`
    : 'Independent run';
  setText(root, '[data-role="param-preview"]', JSON.stringify(state.params.snapshot(), null, 2));
  renderRunEvents(root, state);
  bindDynamicStepperEvents(root, state);
}

function scrollActiveStepTop(root) {
  const panel = root.querySelector('[data-role="active-step-panel"]');
  if (panel?.scrollIntoView) {
    panel.scrollIntoView({ block: 'start', behavior: 'smooth' });
  }
}

function bindDynamicStepperEvents(root, state) {
  root.querySelectorAll('[data-project-step]').forEach((btn) => {
    btn.addEventListener('click', () => {
      state.activeProjectStep = btn.dataset.projectStep;
      renderStepper(root, state);
      scrollActiveStepTop(root);
    });
  });
  root.querySelector('[data-role="run-active-step"]')?.addEventListener('click', (event) => {
    runReconAction(root, state, event.currentTarget.dataset.action);
  });
  root.querySelector('[data-role="find-project-model"]')?.addEventListener('click', () => {
    syncParamsFromForm(root, state);
    const initial = paramValue(state.params.snapshot(), 'InitialModel');
    const target = (state.projectFiles || []).find((file) => file.relativePath === `models/${initial}` || file.name === initial);
    if (target) {
      state.modelId = target.relativePath;
      const input = root.querySelector('[data-role="model-id"]');
      if (input) input.value = state.modelId;
      setText(root, '[data-role="upload-status"]', `Using ${target.relativePath}`);
      appendLog(root, state, `load: found ${target.relativePath}`);
    } else {
      appendLog(root, state, `load: ${initial || 'InitialModel'} not found in models/`);
    }
  });
  root.querySelector('[data-role="next-project-step"]')?.addEventListener('click', () => {
    const idx = PROJECT_STEPS.findIndex((step) => step.id === state.activeProjectStep);
    state.activeProjectStep = PROJECT_STEPS[Math.min(idx + 1, PROJECT_STEPS.length - 1)]?.id || state.activeProjectStep;
    renderStepper(root, state);
    scrollActiveStepTop(root);
  });
  root.querySelector('[data-role="recon-params"]')?.addEventListener('input', () => syncParamsFromForm(root, state));
  root.querySelector('[data-role="recon-params"]')?.addEventListener('change', () => syncParamsFromForm(root, state));
}

async function loadProjectReconState(root, state, opts = {}) {
  const projectId = state.project?.projectId;
  const previousActiveStep = state.activeProjectStep;
  if (!projectId) {
    state.projectState = { steps: PROJECT_STEPS.map((step) => ({ id: step.id, title: step.title, status: 'ready', outputs: [] })) };
    state.projectFiles = [];
    state.activeProjectStep = opts.preserveActive ? previousActiveStep : PROJECT_STEPS[0].id;
    renderStepper(root, state);
    return;
  }
  try {
    const res = await fetch(`/api/projects/${encodeURIComponent(projectId)}/reconstruction/state`, { headers: { Accept: 'application/json' } });
    const body = res.ok ? await res.json() : {};
    state.projectState = body;
    state.projectFiles = body.files || [];
    state.activeProjectStep = opts.preserveActive ? previousActiveStep : activeStepFromProjectState(body);
  } catch (_) {
    state.projectState = { steps: PROJECT_STEPS.map((step) => ({ id: step.id, title: step.title, status: 'ready', outputs: [] })) };
    state.projectFiles = [];
    state.activeProjectStep = opts.preserveActive ? previousActiveStep : PROJECT_STEPS[0].id;
  }
  renderStepper(root, state);
}

function syncParamsFromForm(root, state) {
  const form = root.querySelector('[data-role="recon-params"]');
  if (!form) return;
  const groups = new Map();
  Array.from(form.elements || []).forEach((control) => {
    if (!control.name) return;
    if (!groups.has(control.name)) groups.set(control.name, []);
    groups.get(control.name).push(control);
  });
  const next = {};
  groups.forEach((controls, name) => {
    if (MULTI_PARAMS.has(name)) {
      next[name] = controls.filter((control) => control.checked).map((control) => control.value);
      return;
    }
    if (BOOLEAN_PARAMS.has(name)) {
      next[name] = Boolean(controls[0]?.checked);
      return;
    }
    const raw = controls[0]?.value ?? '';
    if (NUMERIC_PARAMS.has(name)) {
      const parsed = INTEGER_PARAMS.has(name) ? parseInt(raw, 10) : Number(raw);
      next[name] = Number.isFinite(parsed) ? parsed : '';
      return;
    }
    next[name] = raw;
  });
  if (Object.prototype.hasOwnProperty.call(next, 'topologies')) {
    next.topology = legacyTopologyValue(next.topologies);
    const reference = selectedKcatReferenceTopology({ ...state.params.snapshot(), ...next });
    next.kcatReferenceTopology = reference;
    next.customKcatRxnNameType = reference;
  }
  if (Object.prototype.hasOwnProperty.call(next, 'runComplexAnnotation')
      || Object.prototype.hasOwnProperty.call(next, 'runEcAnnotation')
      || Object.prototype.hasOwnProperty.call(next, 'runMetaNetXIntegration')
      || Object.prototype.hasOwnProperty.call(next, 'metaboliteSources')) {
    next.annotationOptions = {
      runComplexAnnotation: Boolean(next.runComplexAnnotation ?? state.params.get('runComplexAnnotation')),
      runEcAnnotation: Boolean(next.runEcAnnotation ?? state.params.get('runEcAnnotation')),
      metaboliteSources: selectedMetaboliteSources({ ...state.params.snapshot(), ...next }),
      runMetaNetXIntegration: Boolean(next.runMetaNetXIntegration ?? state.params.get('runMetaNetXIntegration')),
    };
    next.annotationStages = annotationStagesFromOptions({ ...state.params.snapshot(), ...next });
  }
  state.params.merge(next);
  setText(root, '[data-role="param-preview"]', JSON.stringify(state.params.snapshot(), null, 2));
}

function hydrateReconForm(root, state) {
  const params = state.params.snapshot();
  const form = root.querySelector('[data-role="recon-params"]');
  if (!form) return;
  const setValue = (name, value) => {
    const field = form.elements.namedItem(name);
    if (field && value != null) field.value = value;
  };
  setValue('organism', params.organism);
  setValue('sigma', params.sigma);
  setValue('f', params.f);
  setValue('kcatReferenceTopology', selectedKcatReferenceTopology(params, state));
  setValue('kcatPredictionModel', params.kcatPredictionModel);
  setValue('customKcatRxnNameType', params.customKcatRxnNameType || selectedKcatReferenceTopology(params, state));
  setValue('carbonSource', params.carbonSource);
  setValue('biomassReaction', params.biomassReaction);
  const syncChecks = (name, selected) => {
    const values = new Set(selected || []);
    form.querySelectorAll(`input[name="${name}"]`).forEach((input) => {
      input.checked = values.has(input.value);
    });
  };
  syncChecks('topologies', selectedTopologies(params));
  syncChecks('annotationStages', params.annotationStages);
  syncChecks('metaboliteSources', selectedMetaboliteSources(params));
  syncChecks('deepLearningModels', params.deepLearningModels);
  const complex = form.elements.namedItem('runComplexAnnotation');
  if (complex) complex.checked = Boolean(params.runComplexAnnotation);
  const ec = form.elements.namedItem('runEcAnnotation');
  if (ec) ec.checked = Boolean(params.runEcAnnotation);
  const metanetx = form.elements.namedItem('runMetaNetXIntegration');
  if (metanetx) metanetx.checked = Boolean(params.runMetaNetXIntegration);
  const custom = form.elements.namedItem('useCustomKcatFile');
  if (custom) custom.checked = Boolean(params.useCustomKcatFile);
  const modelInput = root.querySelector('[data-role="model-id"]');
  if (modelInput && state.modelId) modelInput.value = state.modelId;
  setText(root, '[data-role="upload-status"]', state.modelId ? `Using project model ${state.modelId}` : 'No model loaded.');
}

async function runReconAction(root, state, action) {
  syncParamsFromForm(root, state);
  const stepIdx = RECON_STEPS.findIndex((s) => s.action === action);
  const error = validateStep(state, stepIdx);
  if (error) {
    appendLog(root, state, `${action}: ${error}`);
    return;
  }
  if (action === 'init' && state.project?.projectId) {
    appendLog(root, state, 'init: saving parameter manager');
    try {
      const res = await fetch(`/api/projects/${encodeURIComponent(state.project.projectId)}/params`, {
        method: 'PUT',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ params: state.params.snapshot(), loadManager: true }),
      });
      const body = await res.json().catch(() => ({}));
      if (!res.ok) throw new Error(body.detail || body.error || `HTTP ${res.status}`);
      state.project = body;
      state.modelId = body.modelId || state.modelId;
      appendLog(root, state, 'init: parameter manager saved');
      await loadProjectReconState(root, state, { preserveActive: true });
    } catch (err) {
      appendLog(root, state, `init: ${err.message || err}`);
    }
    return;
  }
  appendLog(root, state, `${action}: queued`);
  try {
    const res = await fetch('/api/recon/run', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(buildReconPayload(state, action)),
    });
    const body = await res.json().catch(() => ({}));
    if (!res.ok || body.ok === false) throw new Error(body.error || body.detail || `HTTP ${res.status}`);
    const result = body.result || body.results?.[0] || {};
    if (result.modelId) state.modelId = result.modelId;
    if (result.ecModelId) state.ecModelId = result.ecModelId;
    if (result.ecModelIds) state.ecModelIds = { ...state.ecModelIds, ...result.ecModelIds };
    appendLog(root, state, `${action}: done`);
  } catch (err) {
    appendLog(root, state, `${action}: ${err.message || err}`);
  }
}

async function uploadModelForState(state, file) {
  const projectId = state.project?.projectId;
  if (!projectId) return uploadFile(file);
  const fd = new FormData();
  fd.append('file', file);
  const res = await fetch(`/api/projects/${encodeURIComponent(projectId)}/models`, { method: 'POST', body: fd });
  let body;
  try { body = await res.json(); } catch (_) { body = null; }
  if (!res.ok) {
    return body || { error_code: 'upload_failed', error_message: `HTTP ${res.status}` };
  }
  return body || {};
}

export function mountRecon(rootEl, opts = {}) {
  if (!rootEl) return () => {};
  const state = createReconState(opts);
  rootEl.dataset.track = 'recon';
  rootEl.innerHTML = `
    <section class="track-page">
      <header class="track-header" data-track="recon">
        <h2>Reconstruction</h2>
        <p>Build project-scoped ecGEMs through the ordered tutorial workflow, including ${deepLearningKcatHtml()} prediction, output review, and validation files. Later steps consume an ecGEM ID; conversion can produce ecGEM IDs for all topologies.</p>
      </header>

      <div class="module-workbench module-workbench--recon recon-project-workbench">
        <aside class="workbench-panel recon-step-nav">
          <div class="recon-panel-head">
            <h3>Workflow</h3>
            <p data-role="project-summary-line">Loading project outputs...</p>
          </div>
          <div class="recon-step-strip" data-role="recon-stepper" aria-label="Reconstruction steps"></div>
        </aside>

        <main class="workbench-main">
          <section class="workbench-card recon-active-step recon-step-body" data-role="active-step-panel">
            <p class="recent__empty">Loading reconstruction workflow...</p>
          </section>
          <section class="workbench-card run-details-card">
            <div class="run-details-card__head">
              <h3>Execution details</h3>
              <p class="help-text">Recent actions and validation messages for the selected workflow.</p>
            </div>
            <ol class="run-events" data-role="run-events" aria-live="polite"></ol>
          </section>
        </main>

        <aside class="workbench-panel recon-output-rail">
          <h3>Project outputs</h3>
          <p class="help-text">All detected files for this project. Use this as a quick overview across steps.</p>
          <div class="project-file-list" data-role="project-file-list"></div>

          <div class="model-input-compact">
            <h4>Manual model input</h4>
            <input type="file" data-role="model-upload" accept=".mat,.xml,.json,.yml,.yaml" />
            <input class="input" data-role="model-id" placeholder="modelId" />
            <p class="help-text" data-role="upload-status">No model loaded.</p>
          </div>
          <pre class="param-preview" data-role="param-preview"></pre>
        </aside>
      </div>
    </section>
  `;

  rootEl.querySelector('[data-role="model-id"]')?.addEventListener('input', (event) => {
    state.modelId = event.target.value.trim();
    setText(rootEl, '[data-role="upload-status"]', state.modelId ? `Using modelId ${state.modelId}` : 'No model loaded.');
  });

  rootEl.querySelector('[data-role="model-upload"]')?.addEventListener('change', async (event) => {
    const file = event.target.files?.[0];
    if (!file) return;
    setText(rootEl, '[data-role="upload-status"]', `Uploading ${file.name}...`);
    try {
      const up = await uploadModelForState(state, file);
      if (up.error_code) throw new Error(up.error_message || up.error_code);
      state.modelId = up.projectRelativePath || up.modelId || up.uploadId || '';
      state.uploads.model = up;
      const modelInput = rootEl.querySelector('[data-role="model-id"]');
      if (modelInput) modelInput.value = state.modelId;
      setText(rootEl, '[data-role="upload-status"]', `Loaded ${up.name || file.name}`);
      await loadProjectReconState(rootEl, state);
    } catch (err) {
      setText(rootEl, '[data-role="upload-status"]', err.message || String(err));
    }
  });

  loadProjectReconState(rootEl, state);

  return function teardown() {
    rootEl.innerHTML = '';
  };
}
