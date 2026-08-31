function specs = algorithmRegistry()
% algorithmRegistry  Declarative catalog of all algorithms, organised by track.
%   Each entry drives the UI (algorithm cards, params) and execution.
%   Add a new bridge by adding an entry here (its run function already exists).
%
% Tracks:
%   recon     - reconstruction (T6-T14)
%   calib     - calibration (T15-T20)
%   analysis  - analysis tasks (T21-T23a)
%   design    - design algorithms (T23b, T23c, algFseof, algOptknock, algOptforce)
%
% Entry fields:
%   id           char  - stable id used by the frontend (and inside MATLAB code)
%   nameKey      char  - i18n key for the algorithm's display name
%   descKey      char  - i18n key for the algorithm's short description
%   track        char  - one of 'recon' | 'calib' | 'analysis' | 'design'
%   supports     cell  - model types this entry applies to: {'GEM'}, {'ecGEM'}, etc.
%   solverNeed   char  - 'LP' | 'MILP' | 'QP' | 'NA'
%   runnable     bool  - true if the bridge can run end-to-end from the UI
%   params       struct array (see schema below)
%   bridgeName   char  - literal MATLAB function name (used by auto-orchestration)
%
% Param schema fields:
%     key       char  - MATLAB-side name (forwarded to the algorithm)
%     labelKey  char  - i18n key for the parameter's short label
%     helpKey   char  - i18n key for the parameter's longer description
%     type      char  - 'int' | 'double' (UI enforces numeric input)
%     default   num   - default value
%     min,max   num   - inclusive bounds
%
% NOTE: `supports` MUST be wrapped in double braces ({{...}}) when passed to
% struct(). A single-brace cell would be expanded into a 1xN struct array of
% separate fields, which silently breaks the registry at runtime. The `e`
% helper below takes care of this wrap.

    p = @(key,labelKey,helpKey,type,def,mn,mx) struct( ...
        'key',key,'labelKey',labelKey,'helpKey',helpKey,'type',type, ...
        'default',def,'min',mn,'max',mx);

    % 0x0 struct array so it concatenates cleanly with non-empty params.
    emptyParams = struct('key','','labelKey','','helpKey','','type','', ...
        'default',0,'min',0,'max',0);
    emptyParams(1) = [];

    fseofParams = [ ...
        p('Iterations',  'param_iterations',     'param_iterations_help',     'int',    10,   2,    100), ...
        p('Coefficient', 'param_coefficient',    'param_coefficient_help',    'double', 0.9,  0,    0.999) ];
    ecFvaParams = [ ...
        p('fraction', 'param_fva_fraction', 'param_fva_fraction_help', 'double', 0.9, 0, 1) ];
    optknockParams = [ ...
        p('maxCandidates','param_max_candidates','param_max_candidates_help','int',   200, 1, 5000), ...
        p('numDel',       'param_num_del',       'param_num_del_help',       'int',     5, 1, 20), ...
        p('minGrowth',    'param_min_growth',    'param_min_growth_help',    'double',0.1, 0, 1) ];
    optforceParams = [ ...
        p('k',            'param_optforce_k',     'param_optforce_k_help',     'int',   2,   1, 10), ...
        p('nsets',        'param_optforce_nsets', 'param_optforce_nsets_help', 'int',   1,   1, 5), ...
        p('maxCandidates','param_max_candidates','param_max_candidates_help','int', 500,   1, 5000) ];
    okoParams = [ ...
        p('targetProductFold',          'param_target_fold',         'param_target_fold_help',          'double', 2.0,   1.0,  10.0), ...
        p('growthFraction',             'param_growth_fraction',     'param_growth_fraction_help',      'double', 0.99,  0.0,  1.0),  ...
        p('enzymeReferenceGrowthFraction','param_ref_growth',        'param_ref_growth_help',           'double', 0.98,  0.0,  1.0),  ...
        p('referenceGrowthFraction',    'param_lp1_growth',          'param_lp1_growth_help',           'double', 0.99,  0.0,  1.0),  ...
        p('referenceProductFraction',   'param_lp1_product',         'param_lp1_product_help',          'double', 0.99,  0.0,  1.0),  ...
        p('enzymeTolerance',            'param_enzyme_tol',          'param_enzyme_tol_help',           'double', 0.1,   0.0,  1.0),  ...
        p('significance',               'param_significance',        'param_significance_help',         'double', 1e-8,  1e-12,1e-3), ...
        p('foldLimit',                  'param_fold_limit',          'param_fold_limit_help',           'double', 10.0,  1.01, 100.0),...
        p('warmupProductFold',          'param_warmup_product',      'param_warmup_product_help',       'double', 2.0,   1.0,  10.0), ...
        p('timeLimit',                  'param_time_limit',          'param_time_limit_help',           'double', 900.0, 10.0, 7200.0) ];
    okoPlusParams = [okoParams; ...
        p('abundanceWeight',            'param_abundance_weight',    'param_abundance_weight_help',     'double', 10.0,  0.0,  100.0) ];
    sensitivityParams = [ ...
        p('targetGrowth','param_target_growth','param_target_growth_help','double',0.4, 0, 5), ...
        p('factor',      'param_factor',       'param_factor_help',       'double',1.1, 0.01, 100) ];

    % Helper: build one registry entry.
    % CRITICAL: `supports` is wrapped in {{ }} (cell-of-cell) so struct() stores
    % it as a single CELL field. A bare cell would expand into a 1xN struct array.
    e = @(id,nameKey,descKey,track,supports,solverNeed,runnable,params,bridgeName) ...
        struct('id',id,'nameKey',nameKey,'descKey',descKey, ...
               'track',track,'supports',{supports},'solverNeed',solverNeed, ...
               'runnable',runnable,'params',params,'bridgeName',bridgeName);

    % Build specs incrementally with a cell-array accumulator (end+1 pattern).
    % We deliberately do NOT use an inline `[ ... ]` literal here because:
    %   (a) inline `[ ... ]` does not allow `;`-separated rows with comments
    %       between them (each `%` row break ends the current row, which causes
    %       vertcat to fail with non-rectangular dimensions);
    %   (b) the accumulator pattern lets us comment freely between rows.
    local_specs = {};

    % ---- recon track ----
    local_specs{end+1} = e('mdpInitProject',      'algo_init_project_name',      'algo_init_project_desc',      'recon',   {'GEM','ecGEM'}, 'NA',   false, emptyParams, 'mdpInitProject');
    local_specs{end+1} = e('mdpLoadModel',        'algo_load_model_name',        'algo_load_model_desc',        'recon',   {'GEM'},         'NA',   true,  emptyParams, 'mdpLoadModel');
    local_specs{end+1} = e('mdpConvertecModel',   'algo_convertec_name',         'algo_convertec_desc',         'recon',   {'GEM'},         'NA',   true,  emptyParams, 'mdpConvertecModel');
    local_specs{end+1} = e('mdpAnnotate',         'algo_annotate_name',          'algo_annotate_desc',          'recon',   {'ecGEM'},       'NA',   true,  emptyParams, 'mdpAnnotate');
    local_specs{end+1} = e('mdpDlPredict',        'algo_dl_predict_name',        'algo_dl_predict_desc',        'recon',   {'ecGEM'},       'NA',   true,  emptyParams, 'mdpDlPredict');
    local_specs{end+1} = e('mdpKcatCompare',      'algo_kcat_compare_name',      'algo_kcat_compare_desc',      'recon',   {'ecGEM'},       'NA',   true,  emptyParams, 'mdpKcatCompare');
    local_specs{end+1} = e('mdpKcatMerge',        'algo_kcat_merge_name',        'algo_kcat_merge_desc',        'recon',   {'ecGEM'},       'NA',   true,  emptyParams, 'mdpKcatMerge');
    local_specs{end+1} = e('mdpGrowthPredict',    'algo_growth_predict_name',    'algo_growth_predict_desc',    'recon',   {'ecGEM'},       'LP',   true,  emptyParams, 'mdpGrowthPredict');

    % ---- calib track ----
    local_specs{end+1} = e('mdpApplySluice',       'algo_apply_sluice_name',       'algo_apply_sluice_desc',       'calib', {'ecGEM'}, 'NA',   true, emptyParams,      'mdpApplySluice');
    local_specs{end+1} = e('mdpKcatRepoInit',      'algo_kcat_repo_init_name',     'algo_kcat_repo_init_desc',     'calib', {'ecGEM'}, 'NA',   true, emptyParams,      'mdpKcatRepoInit');
    local_specs{end+1} = e('mdpSensitivityTuning', 'algo_sensitivity_tuning_name', 'algo_sensitivity_tuning_desc', 'calib', {'ecGEM'}, 'NA',   true, sensitivityParams, 'mdpSensitivityTuning');
    local_specs{end+1} = e('mdpGauks',             'algo_gauks_name',              'algo_gauks_desc',              'calib', {'ecGEM'}, 'LP',   true, emptyParams,      'mdpGauks');
    local_specs{end+1} = e('mdpBayesian',          'algo_bayesian_name',           'algo_bayesian_desc',           'calib', {'ecGEM'}, 'LP',   true, emptyParams,      'mdpBayesian');
    local_specs{end+1} = e('mdpPresto',            'algo_presto_name',             'algo_presto_desc',             'calib', {'ecGEM'}, 'QP',   true, emptyParams,      'mdpPresto');

    % ---- analysis track ----
    local_specs{end+1} = e('mdpEcFva',            'algo_ecfva_name',            'algo_ecfva_desc',             'analysis', {'ecGEM'}, 'LP',  true, ecFvaParams,    'mdpEcFva');
    local_specs{end+1} = e('mdpKnockout',         'algo_knockout_name',         'algo_knockout_desc',          'analysis', {'ecGEM'}, 'LP',  true, emptyParams,    'mdpKnockout');
    local_specs{end+1} = e('mdpProteinAnalysis',  'algo_protein_analysis_name', 'algo_protein_analysis_desc',  'analysis', {'ecGEM'}, 'NA',  true, emptyParams,    'mdpProteinAnalysis');

    % ---- design track (bridges) ----
    local_specs{end+1} = e('mdpRunFseof',         'algo_fseof_name',            'algo_fseof_desc',             'design',   {'GEM','ecGEM'}, 'LP',   true,  fseofParams,    'mdpRunFseof');
    local_specs{end+1} = e('mdpRunOptknock',      'algo_optknock_name',         'algo_optknock_desc',          'design',   {'GEM'},         'MILP', true,  optknockParams, 'mdpRunOptknock');
    local_specs{end+1} = e('mdpRunOptforce',      'algo_optforce_name',         'algo_optforce_desc',          'design',   {'GEM'},         'MILP', true,  optforceParams, 'mdpRunOptforce');

    % ---- design track (OKO / OKO+) ----
    local_specs{end+1} = e('mdpRunOko',           'algo_oko_name',              'algo_oko_desc',               'design',   {'ecGEM'},       'MILP', true,  okoParams,      'mdpRunOko');
    local_specs{end+1} = e('mdpRunOkoPlus',       'algo_oko_plus_name',         'algo_oko_plus_desc',          'design',   {'ecGEM'},       'MILP', true,  okoPlusParams,  'mdpRunOkoPlus');

    % ---- design track (pure algorithms living in scripts/StrainDesign/algorithms/) ----
    local_specs{end+1} = e('fseof',               'algo_fseof_name',            'algo_fseof_desc',             'design',   {'GEM','ecGEM'}, 'LP',   false, fseofParams,    'algFseof');
    local_specs{end+1} = e('optknock',            'algo_optknock_name',         'algo_optknock_desc',          'design',   {'GEM'},         'MILP', false, optknockParams, 'algOptknock');
    local_specs{end+1} = e('optforce',            'algo_optforce_name',         'algo_optforce_desc',          'design',   {'GEM'},         'MILP', false, optforceParams, 'algOptforce');
    local_specs{end+1} = e('oko',                 'algo_oko_name',              'algo_oko_desc',               'design',   {'ecGEM'},       'MILP', false, okoParams,      'algOko');
    local_specs{end+1} = e('okoplus',             'algo_oko_plus_name',         'algo_oko_plus_desc',          'design',   {'ecGEM'},       'MILP', false, okoPlusParams,  'algOkoPlus');

    % Convert the cell array of entries into a vertical struct array.
    specs = vertcat(local_specs{:});
end
