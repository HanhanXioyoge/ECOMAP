function OUT = AnalyzeKcatMatches(MATCH, DeepLearningModel, saveAnalysisMat, saveFigs, figFormat, figResolutionDPI, parameters)
% AnalyzeKcatMatches
% ------------------
% Purpose
%   Given MATCH tables produced by BuildKcatMatches (one table per model),
%   compute per-model statistics on:
%     (A) all matched pairs (overall)
%     (B) the intersection ("common set") of matched samples shared by ALL
%         selected models, where a sample is defined by (ProteinID + metabolite identity).
%   Additionally, generate per-model scatter plots on the COMMON set with the
%   number of points annotated in the figure, and optionally save the OUT
%   struct and figures to parameters.analysisDir.
%
% Inputs
%   MATCH            : struct with fields {DLKcat, UniKP, CatPred} (a subset is fine).
%                      Each field is a table with columns from BuildKcatMatches, e.g.:
%                      ReactionName, ProteinID, ec, MetaNetXID, InChIKey, Substrate,
%                      n_group, predicted_kcat_log10, exp_kcat_log10, isComplex,
%                      predicted_kcat, exp_kcat.
%   DeepLearningModel: string|char|cellstr — subset of {'DLKcat','UniKP','CatPred'}.
%   parameters       : struct (if empty, ParameterManager.getParams() will be used).
%                      If saving, requires parameters.analysisDir.
%   saveAnalysisMat  : logical (default false) — save OUT as AnalyzeKcatMatches.mat
%   saveFigs         : logical (default false) — save per-model figures (COMMON set only)
%   figFormat        : char (default 'png') — {'png','tif','jpg','pdf','svg'}
%   figResolutionDPI : integer (default 300)
%
% Output
%   OUT : struct with fields
%     .models               : model list used
%     .overall.(tag)        : stats on ALL matched pairs of that model (r, p, n, MAE_log10)
%     .common.(tag)         : stats on COMMON-set pairs only (r, p, n, MAE_log10)
%     .summary_overall      : table summarizing overall stats per model
%     .summary_common       : table summarizing common-set stats per model
%     .commonKeysCount      : number of common samples across ALL models
%     .commonTables.(tag)   : the subset table per model restricted to the common set
%
% Notes
%   - "Common set" samples are identified by a stable key:
%       sampleKey = ProteinID + "§" + key, where key prioritizes
%         InChIKey > MetaNetXID > normalized(Substrate).
%   - Figures show ONLY the COMMON set (to enable fair comparison across models).
%   - If the common set is empty, figures indicate "No common matches".
%
% Example
%   OUT = AnalyzeKcatMatches(MATCH, {'DLKcat','UniKP','CatPred'}, params, true, true, 'png', 300);

    % ------------------------- Defaults & args -------------------------
    if nargin < 6 || isempty(figResolutionDPI), figResolutionDPI = 300; end
    if nargin < 5 || isempty(figFormat),        figFormat        = 'png'; end
    if nargin < 4 || isempty(saveFigs),         saveFigs         = false; end
    if nargin < 3 || isempty(saveAnalysisMat),  saveAnalysisMat  = false; end

    if nargin < 7 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters), error('ParameterManager is not set.'); end
    end

    if isstring(DeepLearningModel) || ischar(DeepLearningModel)
        DeepLearningModel = cellstr(DeepLearningModel);
    end
    if nargin < 2 || isempty(DeepLearningModel)
        % Use only models present in MATCH
        possible = {'DLKcat','UniKP','CatPred'};
        present  = possible(isfield(MATCH, possible));
        if isempty(present)
            error('AnalyzeKcatMatches:NoModels','MATCH has no recognized model fields.');
        end
        DeepLearningModel = present;
    end
    DeepLearningModel = unique(strtrim(DeepLearningModel(:)'));
    validTags = {'DLKcat','UniKP','CatPred'};
    if any(~ismember(DeepLearningModel, validTags))
        error('DeepLearningModel must be a subset of {DLKcat, UniKP, CatPred}.');
    end
    % Ensure MATCH has the requested models
    for k = 1:numel(DeepLearningModel)
        tag = DeepLearningModel{k};
        if ~isfield(MATCH, tag)
            error('AnalyzeKcatMatches:MissingField','MATCH.%s not found.', tag);
        end
    end

    % Validate output directory if saving
    analysisDir = getfield_def(parameters, 'analysisDir', "");
    if (saveAnalysisMat || saveFigs)
        if analysisDir == ""
            error('AnalyzeKcatMatches:MissinganalysisDir', ...
                 'parameters.analysisDir is required when saving mat or figures.');
        end
        analysisDir = fullfile(analysisDir, 'AnalyzeKcatMatches');
        if ~exist(analysisDir, 'dir'), mkdir(analysisDir); end
    end

    % ------------------------- Build per-model sample keys -------------------------
    % For each model table, construct a stable "sampleKey" (unique within that table):
    %   sampleKey = ProteinID + "§" + (InChIKey | MetaNetXID | normalized(Substrate))
    modelKeys   = containers.Map();   % tag -> string array of keys
    keyedTables = struct();           % tag -> table with added sampleKey
    for k = 1:numel(DeepLearningModel)
        tag = DeepLearningModel{k};
        T   = MATCH.(tag);

        if isempty(T)
            % Keep consistent shape (no rows)
            T.sampleKey = strings(0,1);
            modelKeys(tag) = T.sampleKey;
            keyedTables.(tag) = T;
            continue;
        end

        % Ensure string columns exist (robustness in case upstream changed)
        reqStr = {'ProteinID','InChIKey','MetaNetXID','Substrate'};
        for c = 1:numel(reqStr)
            cn = reqStr{c};
            if ~ismember(cn, T.Properties.VariableNames)
                T.(cn) = strings(height(T),1);
            else
                if ~isstring(T.(cn)), T.(cn) = string(T.(cn)); end
                T.(cn)(ismissing(T.(cn))) = "";
            end
        end

        % Build normalized substrate (for fallback)
        SubN = normalize_substrate_name(T.Substrate);

        % Choose key priority
        useIK   = T.InChIKey ~= "";
        useMNX  = ~useIK & (T.MetaNetXID ~= "");
        useName = ~useIK & ~useMNX & (SubN ~= "");

        key = strings(height(T),1);
        key(useIK)   = "IK:"  + T.InChIKey(useIK);
        key(useMNX)  = "MNX:" + T.MetaNetXID(useMNX);
        key(useName) = "NM:"  + SubN(useName);

        % Final sampleKey
        T.sampleKey = T.ProteinID + "§" + key;

        % Store
        modelKeys(tag)   = T.sampleKey;
        keyedTables.(tag)= T;
    end

    % ------------------------- Compute the common set (intersection) -------------------------
    if isempty(DeepLearningModel)
        error('AnalyzeKcatMatches:EmptyModelList','No models provided.');
    end
    % Start with keys of the first model
    tag1 = DeepLearningModel{1};
    commonKeys = unique(modelKeys(tag1), 'stable');
    for k = 2:numel(DeepLearningModel)
        tagk = DeepLearningModel{k};
        keys = unique(modelKeys(tagk), 'stable');
        commonKeys = intersect(commonKeys, keys, 'stable');
        if isempty(commonKeys), break; end
    end

    % ------------------------- Per-model statistics -------------------------
    OUT = struct();
    OUT.models          = DeepLearningModel;
    OUT.overall         = struct();
    OUT.common          = struct();
    OUT.commonKeysCount = numel(commonKeys);
    OUT.commonTables    = struct();

    % overall and common summary tables (collect rows as cell arrays)
    overall_rows = {};
    common_rows  = {};

    for k = 1:numel(DeepLearningModel)
        tag = DeepLearningModel{k};
        T   = keyedTables.(tag);

        % ---- Overall stats (all matched rows) ----
        [r_all, p_all, n_all, mae_all] = compute_stats(T.predicted_kcat_log10, T.exp_kcat_log10);
        OUT.overall.(tag) = struct('corr_r', r_all, 'corr_p', p_all, 'nPoints', n_all, 'MAE_log10', mae_all);
        overall_rows(end+1, :) = {tag, r_all, p_all, n_all, mae_all}; %#ok<AGROW>

        % ---- Common set subset & stats ----
        if ~isempty(T)
            mask = ismember(T.sampleKey, commonKeys);
            Tc   = T(mask, :);
        else
            Tc   = T;
        end
        OUT.commonTables.(tag) = Tc;
        OUT.overallTables.(tag) = T;

        [r_c, p_c, n_c, mae_c] = compute_stats(Tc.predicted_kcat_log10, Tc.exp_kcat_log10);
        OUT.common.(tag) = struct('corr_r', r_c, 'corr_p', p_c, 'nPoints', n_c, 'MAE_log10', mae_c);
        common_rows(end+1, :) = {tag, r_c, p_c, n_c, mae_c}; %#ok<AGROW>
    end

    OUT.summary_overall = cell2table(overall_rows, ...
        'VariableNames', {'Model','Correlation_r','p_value','nPoints','MAE_log10'});

    OUT.summary_common  = cell2table(common_rows, ...
        'VariableNames', {'Model','Correlation_r','p_value','nPoints','MAE_log10'});

    % ------------------------- Plot per model on COMMON set -------------------------
    % One figure per model, using only the common set. Annotate N.
    for k = 1:numel(DeepLearningModel)
        tag = DeepLearningModel{k};
        Tc  = OUT.commonTables.(tag);

        if ~isempty(Tc)
            x = Tc.predicted_kcat_log10;
            y = Tc.exp_kcat_log10;
            valid = isfinite(x) & isfinite(y);
            x = x(valid); y = y(valid);
        else
            x = []; y = [];
        end

        figH = figure('Color','w'); 
        hold on;

        if ~isempty(x)
            scatter(x, y, ...
                15, ...
                [0.20 0.45 0.80], ...
                'filled', ...
                'MarkerFaceAlpha', 0.6, ...
                'MarkerEdgeColor', [0.20 0.45 0.80], ...
                'MarkerEdgeAlpha', 0.6);

            minVal = min([x; y]); maxVal = max([x; y]);
            if ~isfinite(minVal) || ~isfinite(maxVal), minVal = -1; maxVal = 1; end
            if minVal == maxVal
                span = max(1e-3, abs(maxVal));
                minVal = minVal - 0.5*span; maxVal = maxVal + 0.5*span;
            end
            commonLim = [minVal maxVal];
            plot(commonLim, commonLim, 'k--', 'LineWidth', 1.2);

            xlim(commonLim); ylim(commonLim);
            xt = xticks; yticks(xt);

            % Annotate R, R^2 and N inside axes
            stats = OUT.common.(tag);
            R2   = stats.corr_r.^2;
            pval = stats.corr_p;
            N    = stats.nPoints;
            MAE  = stats.MAE_log10;
            
            dx = (commonLim(2) - commonLim(1));
            tx = commonLim(1) + 0.02*dx;
            ty = commonLim(2) - 0.02*dx;
            
            txt = sprintf('R^2 = %.3f\np_value = %.2g\nN = %d\nMAE = %.3f', R2, pval, N, MAE);
            text(tx, ty, txt, ...
                'VerticalAlignment','top', ...
                'FontSize',11, ...
                'FontWeight','bold', ...
                'BackgroundColor','w', ...
                'Margin',4, ...
                'EdgeColor',[0.5 0.5 0.5], ...
                'LineWidth',0.5, ...
                'Interpreter','none'); 
        else
            % No data for this model in the common set
            commonLim = [0 1];
            plot(commonLim, commonLim, 'k--', 'LineWidth', 1.2);
            xlim(commonLim); ylim(commonLim);
            text(mean(commonLim), mean(commonLim), 'No common matches', ...
                'HorizontalAlignment','center','FontWeight','bold');
        end

        axis square; box on;
        xlabel('Predicted k_{cat} (log_{10} s^{-1})', 'FontWeight','bold');
        ylabel('Experimental k_{cat} (log_{10} s^{-1})', 'FontWeight','bold');
        title(sprintf('%s vs Experimental k_{cat} (Common set)', tag), 'FontWeight','bold');

        ax = gca;
        ax.LineWidth   = 1.2;
        ax.FontSize    = 12;
        ax.FontName    = 'Arial';
        ax.TickDir     = 'out';
        ax.TickLength  = [0.015 0.015];
        ax.XMinorTick  = 'off';
        ax.YMinorTick  = 'off';
        set(gca, 'LooseInset', get(gca,'TightInset'));
        hold off;

        if saveFigs
            outFigName = fullfile(analysisDir, sprintf('Benchmark_%s.%s', tag, figFormat));
            try
                if any(strcmpi(figFormat, {'png','jpg','tif'}))
                    exportgraphics(figH, outFigName, 'Resolution', figResolutionDPI);
                else
                    exportgraphics(figH, outFigName);
                end
            catch
                % Fallback for older MATLAB
                switch lower(figFormat)
                    case {'png','tif','tiff','jpg','jpeg','bmp'}
                        print(figH, outFigName, ['-d' lower(figFormat)], ['-r' num2str(figResolutionDPI)]);
                    case 'pdf'
                        print(figH, outFigName, '-dpdf');
                    otherwise
                        warning('AnalyzeKcatMatches:FigureSave', ...
                                'Unknown figFormat "%s"; saving as PNG fallback.', figFormat);
                        print(figH, [outFigName '.png'], '-dpng', ['-r' num2str(figResolutionDPI)]);
                end
            end
        end
    end

    % ------------------------- Save OUT (optional) -------------------------
    if saveAnalysisMat
        outMatPath = fullfile(analysisDir, 'AnalyzeKcatMatches.mat');
        try
            OUT_saved = OUT;
            save(outMatPath, 'OUT_saved', '-v7.3');
        catch ME
            warning('AnalyzeKcatMatches:SaveMatFailed', 'Failed to save AnalyzeKcatMatches.mat: %s', ME.message);
        end
    end
end


% ============================ Helpers ============================
function sN = normalize_substrate_name(s)
% Normalize metabolite names for robust equality:
% - lowercase, trim
% - strip trailing bracketed notes: "xxx [c]" or "(cytosol)"
% - unify separators to single space
% - collapse multiple spaces
    if ~isstring(s), s = string(s); end
    sN = lower(strtrim(s));
    sN = regexprep(sN, '\s*(\[[^\]]*\]|\([^\)]*\))\s*$', '', 'once');
    sN = regexprep(sN, '[_\-\,;]+', ' ');
    sN = regexprep(sN, '\s+', ' ');
    sN = strtrim(sN);
end

function val = getfield_def(S, fname, def)
% Safe S.(fname) with default.
    if ~isstruct(S), val = def; return; end
    if ~isfield(S, fname) || isempty(S.(fname)), val = def; else, val = S.(fname); end
end
