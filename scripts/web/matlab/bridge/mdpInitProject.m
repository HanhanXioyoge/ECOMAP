function result = mdpInitProject(project_name, project_path)
%MDPINITPROJECT Initialise a new ECOMAP project.
%   result = mdpInitProject(project_name, project_path) returns the standard
%   bridge envelope (see CONTRACT.md) carrying:
%       .project_id            (char)  --  UUID-style handle for the new project
%       .param_template_path   (char)  --  path to the auto-generated
%                                       <project>ParameterManagement.m
%       .default_params        (struct)  --  {c_source, bioRxn, sigma, Ptot, f}
    addpath_once(fullfile(project_path, 'scripts', 'Reconstruction'));
    addpath_once(fullfile(project_path, 'scripts', 'ParameterManagement'));
    bridge_log('mdpInitProject', 'Initializing project %s at %s', project_name, project_path);
    try
        InitializeECOMAPproject(project_name, project_path);
    catch err
        result = make_err('err_init_fail', err.message);
        return;
    end
    safe_project_name = matlab.lang.makeValidName(project_name);
    template_path = fullfile(project_path, 'projects', project_name, [safe_project_name 'ParameterManagement.m']);
    if ~exist(template_path, 'file')
        result = make_err('err_init_fail', 'ParameterManagement template not generated');
        return;
    end
    defaults = struct( ...
        'c_source', 'EX_glc__D_e', ...
        'bioRxn',   'BIOMASS_Ec_iML1515_core_75p37M', ...
        'sigma',    0.5, ...
        'Ptot',     0.3, ...
        'f',        0.5);
    payload = struct('project_id', char(java.util.UUID.randomUUID.toString), ...
                     'param_template_path', template_path, ...
                     'default_params', defaults);
    result = make_ok(payload);
end
