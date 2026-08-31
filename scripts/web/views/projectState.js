export function stageForModelType(modelType) {
  return modelType === 'ecGEM' ? 'Ready for Calibration' : 'Ready for Reconstruction';
}

export function createProjectDraft(opts = {}) {
  return {
    name: opts.name || '',
    modelType: opts.modelType || '',
  };
}

export function validateProjectDraft(draft) {
  const errors = [];
  if (!draft?.name?.trim()) errors.push('name');
  if (!['GEM', 'ecGEM'].includes(draft?.modelType)) errors.push('modelType');
  return errors;
}

export function buildProjectPayload(draft) {
  return {
    name: draft.name.trim(),
    modelType: draft.modelType,
    stage: stageForModelType(draft.modelType),
  };
}

export function projectModuleStatus(project, moduleId) {
  const modelType = project?.modelType || 'GEM';
  const hasEcGem = modelType === 'ecGEM' || Boolean(project?.ecModelId || project?.ecGEMId);
  if (moduleId === 'reconstruction') {
    return modelType === 'GEM' ? 'Ready' : 'Optional';
  }
  if (moduleId === 'calibration') {
    return hasEcGem ? 'Ready' : 'Needs reconstructed ecGEM';
  }
  if (moduleId === 'analysis' || moduleId === 'design') {
    return hasEcGem ? 'Ready' : 'Needs calibrated ecGEM';
  }
  return 'Ready';
}

export function projectDisplayOrganism(project) {
  const params = project?.params || {};
  return params.org_name || project?.org_name || project?.organism || params.organism || '';
}

export function rememberCurrentProject(projectId) {
  try { localStorage.setItem('ecomap_current_project', projectId || ''); } catch (_) {}
}

export function getCurrentProjectId() {
  try { return localStorage.getItem('ecomap_current_project') || ''; } catch (_) { return ''; }
}
