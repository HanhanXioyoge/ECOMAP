%{
[ecModel, noSMILES, noInChIKey] = getMetinfo(ecModel, 2);
[ecModel, foundComplex, proposedComplex] = applyComplexdata(ecModel, complexInfo);
ecModel = convertecModel(model, 'integrated');
DeepLearningModel = {'DLKcat', 'UniKP', 'CatPred'};
OUT = dbKcatMatch(DeepLearningModel);
%}
%{
model = loadModel(modelname, 'Tradition');
ecModel_integrated    = convertecModel(model, 'integrated');
DLKcat();
UniKP();
CatPred();
%}

%{
DeepLearningModel = {'DLKcat', 'UniKP', 'CatPred'};

project_name = 'eciML1515';
project_path = findECOMAProot;
InitializeECOMAPproject(project_name, project_path);
ParameterManagerLocation = fullfile(findECOMAProot, project_name, [project_name 'ParameterManagement.m']); 
ParameterManager.getParams(ParameterManagerLocation);
ExecutePrediction(DeepLearningModel);

project_name = 'ecYeast';
project_path = findECOMAProot;
InitializeECOMAPproject(project_name, project_path);
ParameterManagerLocation = fullfile(findECOMAProot, project_name, [project_name 'ParameterManagement.m']); 
ParameterManager.getParams(ParameterManagerLocation);
ExecutePrediction(DeepLearningModel);

project_name = 'eciCW773';
project_path = findECOMAProot;
InitializeECOMAPproject(project_name, project_path);
ParameterManagerLocation = fullfile(findECOMAProot, project_name, [project_name 'ParameterManagement.m']); 
ParameterManager.getParams(ParameterManagerLocation);
ExecutePrediction(DeepLearningModel);

project_name = 'ecHuman';
project_path = findECOMAProot;
InitializeECOMAPproject(project_name, project_path);
ParameterManagerLocation = fullfile(findECOMAProot, project_name, [project_name 'ParameterManagement.m']); 
ParameterManager.getParams(ParameterManagerLocation);
ExecutePrediction(DeepLearningModel);
%}

%{
out_eciML1515 = dbKcatMatch(DeepLearningModel);
out_ecYeast = dbKcatMatch(DeepLearningModel);
out_ecHuman = dbKcatMatch(DeepLearningModel);
out_eciCW773 = dbKcatMatch(DeepLearningModel);
%}
%{
UseproposedComplex=true;
[ecModel, foundComplex, proposedComplex] = applyComplexdata(ecModel, complexInfo, 'select');
%}

% OUTc = BenchmarkComplexImpact(out_eciML1515_3, complex_name);
% OUT = AnalyzeKcatMatches(MATCH, DeepLearningModel,true, true);
MATCH = BuildKcatMatches(DeepLearningModel, complex_name);
% [ecModel, foundComplex_test, proposedComplex_test] = applyComplexdata(ecModel, complexInfo, "select");