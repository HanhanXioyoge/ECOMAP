classdef ReconstructionApp < handle
    % ReconstructionApp
    % GUI Application for ECOMAP Enzyme-Constraint Model Reconstruction
    % New workflow based on eciML1515_reconstruction_tutorial.mlx

    properties
        UIFigure
        StepList
        ContentPanels
        StatusText
        Model
        Parameters
        KcatList
        StepCompleted(9, 1) logical = false(9, 1)
        NavButtons
    end

    properties (Access = private)
        StepNames = {...
            '1. Initialize Project', ...
            '2. Configure Parameters', ...
            '3. Import Model', ...
            '4. Model Conversion', ...
            '5. Fill ecModel Info', ...
            '6. DL kcat Prediction', ...
            '7. Integrate kcat', ...
            '8. Growth Rate Solution', ...
            '9. Save Model'};
        CurrentStep = 1;
    end

    methods

        function app = ReconstructionApp()
            app.Model = [];
            app.Parameters = struct();
            app.KcatList = [];

            app.createFigure();
            app.updateStep(1);
        end
    end

    methods (Access = private)

        function createFigure(app)
            app.UIFigure = figure('Name', 'ECOMAP Reconstruction Tool', ...
                'Position', [100 100 1280 800], ...
                'Resize', 'off', ...
                'Color', [0.95 0.95 0.95], ...
                'WindowStyle', 'normal', ...
                'MenuBar', 'none', ...
                'ToolBar', 'none', ...
                'NumberTitle', 'off', ...
                'CloseRequestFcn', @app.onClose);

            % ============ SIDEBAR ============
            sidebar = uipanel(app.UIFigure);
            sidebar.Title = '';
            sidebar.BorderWidth = 0;
            sidebar.Units = 'pixels';
            sidebar.Position = [0 0 200 800];
            sidebar.BackgroundColor = [0.20 0.25 0.35];

            % Logo area
            logoPanel = uipanel(sidebar);
            logoPanel.Title = '';
            logoPanel.BorderWidth = 0;
            logoPanel.Units = 'pixels';
            logoPanel.Position = [0 740 200 60];
            logoPanel.BackgroundColor = [0.15 0.18 0.28];

            t = uicontrol(logoPanel, 'Style', 'text');
            t.String = 'ECOMAP';
            t.FontSize = 22;
            t.FontWeight = 'bold';
            t.ForegroundColor = [0.35 0.75 0.85];
            t.Units = 'pixels';
            t.Position = [15 18 170 32];
            t.HorizontalAlignment = 'left';

            t2 = uicontrol(logoPanel, 'Style', 'text');
            t2.String = 'Enzyme-Constraint Model';
            t2.FontSize = 8;
            t2.ForegroundColor = [0.65 0.65 0.65];
            t2.Units = 'pixels';
            t2.Position = [15 3 170 15];
            t2.HorizontalAlignment = 'left';

            % Steps header
            t = uicontrol(sidebar, 'Style', 'text');
            t.String = 'WORKFLOW';
            t.FontSize = 9;
            t.FontWeight = 'bold';
            t.ForegroundColor = [0.45 0.55 0.65];
            t.Units = 'pixels';
            t.Position = [15 710 170 15];

            % Step list
            app.StepList = uicontrol(sidebar, 'Style', 'listbox');
            app.StepList.String = app.StepNames;
            app.StepList.FontSize = 11;
            app.StepList.Units = 'pixels';
            app.StepList.Position = [10 50 180 635];
            app.StepList.BackgroundColor = [0.25 0.30 0.40];
            app.StepList.ForegroundColor = [0.90 0.90 0.90];
            app.StepList.Callback = @app.onStepChanged;
            app.StepList.Value = 1;

            % Bottom version info
            t = uicontrol(sidebar, 'Style', 'text');
            t.String = 'v1.0';
            t.FontSize = 8;
            t.ForegroundColor = [0.35 0.35 0.35];
            t.Units = 'pixels';
            t.Position = [85 8 50 15];
            t.HorizontalAlignment = 'center';

            % ============ MAIN CONTENT ============
            % Background panel - unified light gray color
            mainBg = uipanel(app.UIFigure);
            mainBg.Title = '';
            mainBg.BorderWidth = 0;
            mainBg.Units = 'pixels';
            mainBg.Position = [200 25 1065 760];
            mainBg.BackgroundColor = [0.92 0.94 0.98];

            % Content panels - same color as background
            app.ContentPanels = cell(9, 1);
            for i = 1:9
                app.ContentPanels{i} = uipanel(mainBg);
                app.ContentPanels{i}.Title = '';
                app.ContentPanels{i}.BorderWidth = 0;
                app.ContentPanels{i}.Units = 'pixels';
                app.ContentPanels{i}.Position = [0 0 1065 760];
                app.ContentPanels{i}.BackgroundColor = [0.92 0.94 0.98];
                app.ContentPanels{i}.Visible = 'off';
            end

            % ============ STATUS BAR ============
            statusBar = uipanel(app.UIFigure);
            statusBar.Title = '';
            statusBar.BorderWidth = 0;
            statusBar.Units = 'pixels';
            statusBar.Position = [200 0 1065 25];
            statusBar.BackgroundColor = [0.92 0.92 0.92];

            app.StatusText = uicontrol(statusBar, 'Style', 'text');
            app.StatusText.String = 'Ready';
            app.StatusText.FontSize = 10;
            app.StatusText.ForegroundColor = [0.30 0.30 0.30];
            app.StatusText.Units = 'pixels';
            app.StatusText.Position = [15 4 900 18];
            app.StatusText.HorizontalAlignment = 'left';

            % Initialize NavButtons struct
            app.NavButtons = struct();
        end

        function updateStep(app, stepNum)
            app.CurrentStep = stepNum;

            % Update step list visual selection
            app.StepList.Value = stepNum;

            % Update step names to show completion status
            app.updateStepListDisplay();

            for i = 1:9
                app.ContentPanels{i}.Visible = 'off';
            end
            app.ContentPanels{stepNum}.Visible = 'on';

            switch stepNum
                case 1, app.buildStep1_Initialize();
                case 2, app.buildStep2_Configure();
                case 3, app.buildStep3_ImportModel();
                case 4, app.buildStep4_Conversion();
                case 5, app.buildStep5_FillInfo();
                case 6, app.buildStep6_DLPrediction();
                case 7, app.buildStep7_IntegrateKcat();
                case 8, app.buildStep8_GrowthRate();
                case 9, app.buildStep9_Save();
            end

            % Update navigation buttons
            app.updateNavigationButtons();
            drawnow;
        end

        function updateNavigationButtons(app)
            % Update Previous button state
            if isfield(app.NavButtons, 'prevBtn') && isvalid(app.NavButtons.prevBtn)
                if app.CurrentStep == 1
                    app.NavButtons.prevBtn.Enable = 'off';
                else
                    app.NavButtons.prevBtn.Enable = 'on';
                end
            end

            % Update Next button state
            if isfield(app.NavButtons, 'nextBtn') && isvalid(app.NavButtons.nextBtn)
                if app.CurrentStep == 9
                    app.NavButtons.nextBtn.Enable = 'off';
                else
                    % Next is enabled only if current step is completed
                    app.NavButtons.nextBtn.Enable = iif(app.StepCompleted(app.CurrentStep), 'on', 'off');
                end
            end
        end

        % Helper function to build consistent step header
        function buildStepHeader(app, parent, stepNum, title, description)
            % Header bar - use axes with rectangle for solid background
            headerPanel = uipanel(parent);
            headerPanel.Title = '';
            headerPanel.BorderWidth = 0;
            headerPanel.Units = 'pixels';
            headerPanel.Position = [0 685 1065 70];
            headerPanel.BackgroundColor = [0.20 0.28 0.45];

            % Use axes with rectangle to ensure solid background
            ax = axes('Parent', headerPanel, 'Units', 'pixels', 'Position', [0 0 1065 70], 'Color', [0.20 0.28 0.45], 'XTick', [], 'YTick', []);
            rectangle('Parent', ax, 'Position', [0 0 1065 70], 'FaceColor', [0.20 0.28 0.45], 'EdgeColor', 'none');

            t = uicontrol(headerPanel, 'Style', 'text');
            t.String = sprintf('STEP %d', stepNum);
            t.FontSize = 10;
            t.ForegroundColor = [0.55 0.65 0.75];
            t.BackgroundColor = [0.20 0.28 0.45];
            t.Units = 'pixels';
            t.Position = [25 48 60 15];

            t = uicontrol(headerPanel, 'Style', 'text');
            t.String = title;
            t.FontSize = 18;
            t.FontWeight = 'bold';
            t.ForegroundColor = [1 1 1];
            t.BackgroundColor = [0.20 0.28 0.45];
            t.Units = 'pixels';
            t.Position = [25 18 700 32];

            t = uicontrol(headerPanel, 'Style', 'text');
            t.String = description;
            t.FontSize = 9;
            t.ForegroundColor = [0.60 0.60 0.60];
            t.BackgroundColor = [0.20 0.28 0.45];
            t.Units = 'pixels';
            t.Position = [25 5 800 12];
        end

        function updateStepListDisplay(app)
            % Update step list to show completion status
            stepNames = app.StepNames;
            for i = 1:9
                if app.StepCompleted(i)
                    % Add checkmark for completed steps
                    if ~contains(stepNames{i}, '✓')
                        stepNames{i} = ['✓ ', stepNames{i}];
                    end
                else
                    % Remove checkmark if step is no longer completed
                    if contains(stepNames{i}, '✓')
                        stepNames{i} = strrep(stepNames{i}, '✓ ', '');
                    end
                end
            end
            app.StepList.String = stepNames;

            % Update listbox styling based on current step and completion
            if ~isempty(app.StepList) && isvalid(app.StepList)
                % Set background color for the entire list to dark theme
                app.StepList.BackgroundColor = [0.22 0.28 0.38];
            end
        end

        function setStatus(app, text)
            if ~isempty(app.StatusText) && isvalid(app.StatusText)
                app.StatusText.String = text;
            end
        end

        % -----------------------------------------------------------------
        % STEP 1: Initialize Project
        % -----------------------------------------------------------------
        function buildStep1_Initialize(app)
            p = app.ContentPanels{1};

            % Use consistent header (header at y=700, height 65)
            app.buildStepHeader(p, 1, 'Initialize Project', 'Create project folder structure');

            % Layout: content from y=30 to y=680 (below header at y=700)
            % Row 1: y=380-520 (140px high)
            % Row 2: y=180-340 (140px high)

            % Row 1: Project Settings Card (y=380, height=140)
            row1Y = 500;
            row1H = 140;
            settingsPanel = uipanel('Parent', p);
            settingsPanel.Title = 'Project Settings';
            settingsPanel.FontSize = 10;
            settingsPanel.Units = 'pixels';
            settingsPanel.Position = [30 row1Y 700 row1H];
            settingsPanel.BackgroundColor = [0.92 0.94 0.98];

            % Project name
            t = uicontrol('Parent', settingsPanel, 'Style', 'text');
            t.String = 'Project Name:';
            t.FontSize = 10;
            t.HorizontalAlignment = 'left';
            t.Units = 'pixels';
            t.Position = [20 90 100 22];

            projEdit = uicontrol('Parent', settingsPanel, 'Style', 'edit');
            projEdit.Units = 'pixels';
            projEdit.Position = [130 88 300 26];
            projEdit.BackgroundColor = [1 1 1];
            projEdit.String = 'ecModelGEM';
            projEdit.Tag = 'projName';

            % Project path
            t = uicontrol('Parent', settingsPanel, 'Style', 'text');
            t.String = 'Project Path:';
            t.FontSize = 10;
            t.HorizontalAlignment = 'left';
            t.Units = 'pixels';
            t.Position = [20 50 100 22];

            pathEdit = uicontrol('Parent', settingsPanel, 'Style', 'edit');
            pathEdit.Units = 'pixels';
            pathEdit.Position = [130 48 300 26];
            pathEdit.BackgroundColor = [1 1 1];
            pathEdit.String = pwd;
            pathEdit.Tag = 'projPath';

            btn = uicontrol('Parent', settingsPanel, 'Style', 'pushbutton');
            btn.String = 'Browse';
            btn.Units = 'pixels';
            btn.Position = [440 46 80 28];
            btn.BackgroundColor = [0.3 0.3 0.4];
            btn.ForegroundColor = [1 1 1];
            btn.Callback = @(~,~) app.browseProjectPath(pathEdit);

            % Create button
            btn = uicontrol('Parent', settingsPanel, 'Style', 'pushbutton');
            btn.String = 'Create Project';
            btn.Units = 'pixels';
            btn.Position = [20 8 150 35];
            btn.BackgroundColor = [0.2 0.5 0.4];
            btn.ForegroundColor = [1 1 1];
            btn.FontSize = 11;
            btn.Callback = {@app.createProjectCallback, projEdit, pathEdit};

            % Row 2: Folder Structure info card (y=170, height=200)
            row2Y = 270;
            row2H = 200;
            infoPanel = uipanel('Parent', p);
            infoPanel.Title = 'Folder Structure';
            infoPanel.FontSize = 10;
            infoPanel.Units = 'pixels';
            infoPanel.Position = [30 row2Y 700 row2H];
            infoPanel.BackgroundColor = [0.92 0.94 0.98];

            t = uicontrol('Parent', infoPanel, 'Style', 'text');
            t.Position = [15 10 670 160];
            t.FontSize = 9;
            t.FontName = 'Courier New';
            t.HorizontalAlignment = 'left';
            t.String = {'ProjectRoot/', ...
                        '  data/', ...
                        '    - uniprot.tsv, kcatData/, growth_rates.tsv', ...
                        '  models/', ...
                        '    - *.xml, *.json, *.yml', ...
                        '  analysis/', ...
                        '    - MultiCondition_Summary.tsv'};

            % Navigation buttons - consistent position at bottom
            app.NavButtons.prevBtn = uicontrol('Parent', p, 'Style', 'pushbutton');
            app.NavButtons.prevBtn.String = 'Previous';
            app.NavButtons.prevBtn.Units = 'pixels';
            app.NavButtons.prevBtn.Position = [800 30 100 35];
            app.NavButtons.prevBtn.Enable = 'off';
            app.NavButtons.prevBtn.BackgroundColor = [0.5 0.5 0.55];

            app.NavButtons.nextBtn = uicontrol('Parent', p, 'Style', 'pushbutton');
            app.NavButtons.nextBtn.String = 'Next';
            app.NavButtons.nextBtn.Units = 'pixels';
            app.NavButtons.nextBtn.Position = [910 30 100 35];
            app.NavButtons.nextBtn.BackgroundColor = [0.2 0.4 0.6];
            app.NavButtons.nextBtn.ForegroundColor = [1 1 1];
            app.NavButtons.nextBtn.Enable = iif(app.StepCompleted(1), 'on', 'off');
            app.NavButtons.nextBtn.Callback = @(~,~) app.updateStep(2);
        end

        function browseProjectPath(app, pathEdit)
            folder = uigetdir(pwd, 'Select Project Folder');
            if folder ~= 0
                pathEdit.String = folder;
            end
        end

        function createProjectCallback(app, hObj, event, projEdit, pathEdit)
            name = strtrim(projEdit.String);
            path = strtrim(pathEdit.String);

            if isempty(name) || isempty(path)
                msgbox('Please enter project name and path.', 'Error', 'warn');
                return;
            end

            % Validate that path exists
            if ~exist(path, 'dir')
                msgbox('The specified project path does not exist. Please select a valid folder.', 'Error', 'warn');
                return;
            end

            % Check if path already ends with name (avoid nested path)
            % Normalize paths for comparison
            pathNormalized = strrep(path, '\', '/');
            nameNormalized = strrep(name, '\', '/');
            if endsWith(pathNormalized, ['/' nameNormalized])
                % Path already contains name, use path directly as fullPath
                fullPath = path;
                parentPath = strrep(path, [filesep name], '');
            else
                % Normal case: combine path and name
                fullPath = fullfile(path, name);
                parentPath = path;
            end

            if exist(fullPath, 'dir')
                msgbox('Project folder already exists.', 'Warning', 'warn');
                % Still allow proceeding - mark as completed
            else
                % Create folders
                mkdir(fullPath, 'data');
                mkdir(fullPath, 'models');
                mkdir(fullPath, 'analysis');
                mkdir(fullPath, 'data/kcatData');

                % Create .keep files
                fclose(fopen(fullfile(fullPath, 'data', '.keep'), 'w'));
                fclose(fopen(fullfile(fullPath, 'models', '.keep'), 'w'));
                fclose(fopen(fullfile(fullPath, 'analysis', '.keep'), 'w'));

                app.setStatus(['Created: ', fullPath]);
            end

            % Verify fullPath was constructed correctly
            if isempty(strfind(fullPath, '\')) || isempty(strfind(fullPath, name))
                error('Project path construction failed. Please check the project name and path.');
            end

            % Store in parameters - store parent path and project name separately
            app.Parameters.projectPath = parentPath;  % Parent directory
            app.Parameters.projectName = name;        % Project name
            app.Parameters.path = fullPath;          % Full project path
            app.Parameters.modelDir = fullfile(fullPath, 'models');
            app.Parameters.dataDir = fullfile(fullPath, 'data');
            app.Parameters.outputDir = fullfile(fullPath, 'analysis');

            % Mark step as completed and enable Next button
            app.StepCompleted(1) = true;
            app.updateStepListDisplay();
            app.NavButtons.nextBtn.Enable = 'on';

            msgbox({['Project ready:'], fullPath, '', 'Click "Next" to continue.'}, 'Success', 'help');
        end

        % -----------------------------------------------------------------
        % STEP 2: Configure Parameters
        % -----------------------------------------------------------------
        function buildStep2_Configure(app)
            p = app.ContentPanels{2};

            % Validate Step 1 was completed
            if ~app.StepCompleted(1)
                msgbox('Please complete Step 1 (Initialize Project) first.', 'Error', 'warn');
                app.updateStep(1);
                return;
            end

            % Use consistent header (header at y=700, height 65)
            app.buildStepHeader(p, 2, 'Configure Parameters', 'Set enzyme constraint and organism parameters');

            % Layout: content from y=30 to y=680 (below header at y=700)
            % Three panels side by side at y=400, height=180

            panelY = 450;
            panelH = 180;
            panelGap = 20;

            % ===== Basic Parameters Card =====
            leftPanel = uipanel('Parent', p);
            leftPanel.Title = 'Basic Parameters';
            leftPanel.FontSize = 10;
            leftPanel.Units = 'pixels';
            leftPanel.Position = [30 panelY+30 320 panelH];
            leftPanel.BackgroundColor = [0.92 0.94 0.98];

            y = 130;
            fields = {'sigma', 'Ptot', 'f'};
            defaults = {0.5, 0.55, 0.55};
            descriptions = {'Avg enzyme saturation', 'Total protein (g/gDCW)', 'Enzyme fraction'};

            for i = 1:3
                t = uicontrol('Parent', leftPanel, 'Style', 'text');
                t.String = fields{i};
                t.FontSize = 9;
                t.HorizontalAlignment = 'left';
                t.Position = [15 y 60 18];

                edit = uicontrol('Parent', leftPanel, 'Style', 'edit');
                edit.Position = [80 y 60 22];
                edit.BackgroundColor = [1 1 1];
                edit.String = num2str(defaults{i});
                edit.Tag = fields{i};

                t = uicontrol('Parent', leftPanel, 'Style', 'text');
                t.String = descriptions{i};
                t.FontSize = 8;
                t.ForegroundColor = [0.5 0.5 0.5];
                t.HorizontalAlignment = 'left';
                t.Position = [150 y-2 160 18];
                y = y - 35;
            end

            % ===== Organism Information Card =====
            centerPanel = uipanel('Parent', p);
            centerPanel.Title = 'Organism Information';
            centerPanel.FontSize = 10;
            centerPanel.Units = 'pixels';
            centerPanel.Position = [370 panelY 320 panelH+30];
            centerPanel.BackgroundColor = [0.92 0.94 0.98];

            y = 155;
            orgFields = {'org_name', 'uniprot_type', 'uniprot_ID', 'geneIDfield', 'taxonomicID'};
            orgDefaults = {'Escherichia coli', 'proteome', 'UP000000625', 'gene_oln', '83333'};
            orgLabels = {'Organism', 'UniProt Type', 'UniProt ID', 'Gene ID', 'Taxonomic ID'};

            for i = 1:5
                t = uicontrol('Parent', centerPanel, 'Style', 'text');
                t.String = orgLabels{i};
                t.FontSize = 9;
                t.HorizontalAlignment = 'left';
                t.Units = 'pixels';
                t.Position = [15 y 90 18];

                edit = uicontrol('Parent', centerPanel, 'Style', 'edit');
                edit.Units = 'pixels';
                edit.Position = [110 y 190 22];
                edit.BackgroundColor = [1 1 1];
                edit.String = orgDefaults{i};
                edit.Tag = orgFields{i};
                y = y - 28;
            end

            % ===== Reaction IDs & Model Type Card =====
            rightPanel = uipanel('Parent', p);
            rightPanel.Title = 'Reaction & Model Settings';
            rightPanel.FontSize = 10;
            rightPanel.Units = 'pixels';
            rightPanel.Position = [710 panelY+30 320 panelH];
            rightPanel.BackgroundColor = [0.92 0.94 0.98];

            y = 130;
            rxnFields = {'c_source', 'bioRxn'};
            rxnDefaults = {'EX_glc__D_e', 'biomass'};
            rxnLabels = {'Carbon Source', 'Biomass Reaction'};

            for i = 1:2
                t = uicontrol('Parent', rightPanel, 'Style', 'text');
                t.String = rxnLabels{i};
                t.FontSize = 9;
                t.HorizontalAlignment = 'left';
                t.Position = [15 y 100 18];

                edit = uicontrol('Parent', rightPanel, 'Style', 'edit');
                edit.Position = [120 y 180 22];
                edit.BackgroundColor = [1 1 1];
                edit.String = rxnDefaults{i};
                edit.Tag = rxnFields{i};
                y = y - 30;
            end

            % Model type
            t = uicontrol('Parent', rightPanel, 'Style', 'text');
            t.String = 'Model Type:';
            t.FontSize = 9;
            t.HorizontalAlignment = 'left';
            t.Position = [15 y 80 18];

            modelTypePopup = uicontrol('Parent', rightPanel, 'Style', 'popupmenu');
            modelTypePopup.String = 'ECOMAP|Tradition|sMOMENT|ECMpy|GECKO';
            modelTypePopup.Value = 2;
            modelTypePopup.Position = [100 y 100 22];
            modelTypePopup.BackgroundColor = [1 1 1];
            modelTypePopup.Tag = 'modeltype';

            % ecModel type
            y = y - 25;
            t = uicontrol('Parent', rightPanel, 'Style', 'text');
            t.String = 'ecModel:';
            t.FontSize = 9;
            t.HorizontalAlignment = 'left';
            t.Position = [15 y 80 18];

            ecTypePopup = uicontrol('Parent', rightPanel, 'Style', 'popupmenu');
            ecTypePopup.String = 'basic|isozyme|integrated';
            ecTypePopup.Value = 3;
            ecTypePopup.Position = [100 y 100 22];
            ecTypePopup.BackgroundColor = [1 1 1];
            ecTypePopup.Tag = 'ecModelType';

            % Save params button
            btn = uicontrol('Parent', p, 'Style', 'pushbutton');
            btn.String = 'Save Parameters';
            btn.Units = 'pixels';
            btn.Position = [30 30 150 40];
            btn.BackgroundColor = [0.2 0.5 0.4];
            btn.ForegroundColor = [1 1 1];
            btn.FontSize = 11;
            btn.Callback = @app.saveParametersCallback;

            % Navigation buttons - consistent bottom position
            app.NavButtons.prevBtn = uicontrol('Parent', p, 'Style', 'pushbutton');
            app.NavButtons.prevBtn.String = 'Previous';
            app.NavButtons.prevBtn.Units = 'pixels';
            app.NavButtons.prevBtn.Position = [800 30 100 35];
            app.NavButtons.prevBtn.BackgroundColor = [0.5 0.5 0.55];
            app.NavButtons.prevBtn.Callback = @(~,~) app.updateStep(1);

            app.NavButtons.nextBtn = uicontrol('Parent', p, 'Style', 'pushbutton');
            app.NavButtons.nextBtn.String = 'Next';
            app.NavButtons.nextBtn.Units = 'pixels';
            app.NavButtons.nextBtn.Position = [910 30 100 35];
            app.NavButtons.nextBtn.BackgroundColor = [0.2 0.4 0.6];
            app.NavButtons.nextBtn.ForegroundColor = [1 1 1];
            app.NavButtons.nextBtn.Enable = iif(app.StepCompleted(2), 'on', 'off');
            app.NavButtons.nextBtn.Callback = @(~,~) app.updateStep(3);
        end

        function saveParametersCallback(app, hObj, event)
            % Validate Step 1 was completed
            if ~app.StepCompleted(1)
                msgbox('Please complete Step 1 (Initialize Project) first.', 'Error', 'warn');
                return;
            end

            % Find all edit fields and popupmenus in step 2
            p = app.ContentPanels{2};

            % Helper function to get string value from uicontrol
            getStrVal = @(h) get(h, 'String');

            % Basic params - use (1) to handle array case
            sigmaEdit = findobj(p, 'Tag', 'sigma');
            PtotEdit = findobj(p, 'Tag', 'Ptot');
            fEdit = findobj(p, 'Tag', 'f');

            % Organism info
            orgNameEdit = findobj(p, 'Tag', 'org_name');
            uniprotTypeEdit = findobj(p, 'Tag', 'uniprot_type');
            uniprotIDEdit = findobj(p, 'Tag', 'uniprot_ID');
            geneIDFieldEdit = findobj(p, 'Tag', 'geneIDfield');
            taxIDEdit = findobj(p, 'Tag', 'taxonomicID');

            % Reaction IDs
            cSourceEdit = findobj(p, 'Tag', 'c_source');
            bioRxnEdit = findobj(p, 'Tag', 'bioRxn');

            % Model types
            modelTypePopup = findobj(p, 'Tag', 'modeltype');
            ecTypePopup = findobj(p, 'Tag', 'ecModelType');

            modelTypes = {'ECOMAP', 'Tradition', 'sMOMENT', 'ECMpy', 'GECKO'};
            ecTypes = {'basic', 'isozyme', 'integrated'};

            % Save to parameters - handle array case for edit fields
            sigmaStr = sigmaEdit.String;
            if iscell(sigmaStr), sigmaStr = sigmaStr{1}; end
            app.Parameters.sigma = str2double(sigmaStr);

            PtotStr = PtotEdit.String;
            if iscell(PtotStr), PtotStr = PtotStr{1}; end
            app.Parameters.Ptot = str2double(PtotStr);

            fStr = fEdit.String;
            if iscell(fStr), fStr = fStr{1}; end
            app.Parameters.f = str2double(fStr);

            % String fields - handle cell array case
            orgNameStr = orgNameEdit.String;
            if iscell(orgNameStr), orgNameStr = orgNameStr{1}; end
            app.Parameters.org_name = orgNameStr;

            uniprotTypeStr = uniprotTypeEdit.String;
            if iscell(uniprotTypeStr), uniprotTypeStr = uniprotTypeStr{1}; end
            app.Parameters.uniprot.type = uniprotTypeStr;

            uniprotIDStr = uniprotIDEdit.String;
            if iscell(uniprotIDStr), uniprotIDStr = uniprotIDStr{1}; end
            app.Parameters.uniprot.ID = uniprotIDStr;

            geneIDFieldStr = geneIDFieldEdit.String;
            if iscell(geneIDFieldStr), geneIDFieldStr = geneIDFieldStr{1}; end
            app.Parameters.uniprot.geneIDfield = geneIDFieldStr;

            taxIDStr = taxIDEdit.String;
            if iscell(taxIDStr), taxIDStr = taxIDStr{1}; end
            app.Parameters.taxonomicID = taxIDStr;

            cSourceStr = cSourceEdit.String;
            if iscell(cSourceStr), cSourceStr = cSourceStr{1}; end
            app.Parameters.c_source = cSourceStr;

            bioRxnStr = bioRxnEdit.String;
            if iscell(bioRxnStr), bioRxnStr = bioRxnStr{1}; end
            app.Parameters.bioRxn = bioRxnStr;

            app.Parameters.modeltype = modelTypes{modelTypePopup.Value};
            app.Parameters.ecModelType = ecTypes{ecTypePopup.Value};

            % Ensure uniprot structure exists
            if ~isfield(app.Parameters, 'uniprot')
                app.Parameters.uniprot = struct();
            end

            % Set default for reviewed if not already set
            if ~isfield(app.Parameters.uniprot, 'reviewed')
                app.Parameters.uniprot.reviewed = true;
            end

            % Initialize PRESTO structure with defaults
            app.Parameters.PRESTO = struct();
            app.Parameters.PRESTO.runParallel = false;
            app.Parameters.PRESTO.ncpu = 4;
            app.Parameters.PRESTO.nIter = 100;
            app.Parameters.PRESTO.epsilon = 0.1;
            app.Parameters.PRESTO.lambda = 0.5;
            app.Parameters.PRESTO.theta = 1.0;

            % Generate ParameterManagement.m file
            try
                % Ensure projectPath is set correctly
                if ~isstruct(app.Parameters)
                    error('Parameters not initialized. Please complete Step 1 first.');
                end

                hasProjectPath = isfield(app.Parameters, 'projectPath') && ~isempty(app.Parameters.projectPath);
                hasProjectName = isfield(app.Parameters, 'projectName') && ~isempty(app.Parameters.projectName);

                if ~hasProjectPath && hasProjectName
                    % Try to extract from path if available
                    if isfield(app.Parameters, 'path') && ~isempty(app.Parameters.path)
                        [app.Parameters.projectPath, ~] = fileparts(app.Parameters.path);
                        hasProjectPath = true;
                    end
                end

                if hasProjectPath && hasProjectName
                    % Ensure path is correctly constructed from projectPath and projectName
                    app.Parameters.path = fullfile(app.Parameters.projectPath, app.Parameters.projectName);
                    app.generateParameterManagementFile();

                    % Initialize ParameterManager with the generated file
                    paramFile = fullfile(app.Parameters.path, [app.Parameters.projectName 'ParameterManagement.m']);
                    if isfile(paramFile)
                        try
                            ParameterManager.getParams(paramFile);
                            app.setStatus('ParameterManager initialized');
                        catch ME2
                            warning('Failed to initialize ParameterManager: %s', ME2.message);
                        end
                    end
                else
                    projPathStr = '';
                    projNameStr = '';
                    if isfield(app.Parameters, 'projectPath'), projPathStr = app.Parameters.projectPath; end
                    if isfield(app.Parameters, 'projectName'), projNameStr = app.Parameters.projectName; end
                    error('Please complete Step 1 (Initialize Project) first. projectPath=%s, projectName=%s', ...
                        projPathStr, projNameStr);
                end
            catch ME
                warning('Failed to generate ParameterManagement.m: %s', ME.message);
            end

            % Mark step as completed and enable Next button
            app.StepCompleted(2) = true;
            app.updateStepListDisplay();
            app.NavButtons.nextBtn.Enable = 'on';

            app.setStatus('Parameters saved');
            msgbox('Parameters saved. Click "Next" to continue.', 'Success', 'help');
        end

        function generateParameterManagementFile(app)
            % Generate ParameterManagement.m file from template
            projectName = app.Parameters.projectName;
            fullPath = app.Parameters.path;
            % projectPath = app.Parameters.projectPath;% This is the complete path (parent + projectName)

            % Extract parent directory and project name
            [parentPath, ~] = fileparts(fullPath);

            % Read Template.m contents
            templatePath = fullfile(fileparts(fileparts(fileparts(which('launchReconstructionGUI.m')))), ...
                'scripts', 'ParameterManagement', 'Template.m');

            if ~exist(templatePath, 'file')
                error('Template.m not found at: %s', templatePath);
            end

            fid = fopen(templatePath, 'r');
            f = fread(fid, '*char')';
            fclose(fid);

            % Replace key values - KEY_PATH is parent dir, KEY_NAME is project name
            f = strrep(f, 'KEY_Template', [projectName 'ParameterManagement']);
            f = strrep(f, 'KEY_PATH', parentPath);
            f = strrep(f, 'KEY_NAME', projectName);

            % Fill in Project information
            InitialModel = '';
            if isfield(app.Parameters, 'InitialModel'), InitialModel = app.Parameters.InitialModel; end
            f = strrep(f, "obj.params.InitialModel = '';", sprintf("obj.params.InitialModel = '%s';", InitialModel));

            modeltype = '';
            if isfield(app.Parameters, 'modeltype'), modeltype = app.Parameters.modeltype; end
            f = strrep(f, "obj.params.modeltype  = '';", sprintf("obj.params.modeltype  = '%s';", modeltype));

            % Fill in Enzyme constraint parameters (sigma, Ptot, f)
            sigma = 0.5;
            if isfield(app.Parameters, 'sigma'), sigma = app.Parameters.sigma; end
            f = strrep(f, 'obj.params.sigma   = ;', sprintf('obj.params.sigma   = %g;', sigma));

            Ptot = 0.55;
            if isfield(app.Parameters, 'Ptot'), Ptot = app.Parameters.Ptot; end
            f = strrep(f, 'obj.params.Ptot    = ;', sprintf('obj.params.Ptot    = %g;', Ptot));

            f_val = 0.55;
            if isfield(app.Parameters, 'f'), f_val = app.Parameters.f; end
            f = strrep(f, 'obj.params.f       = ;', sprintf('obj.params.f       = %g;', f_val));

            % Fill in Organism metadata
            org_name = '';
            if isfield(app.Parameters, 'org_name'), org_name = app.Parameters.org_name; end
            f = strrep(f, "obj.params.org_name        = '';", sprintf("obj.params.org_name        = '%s';", org_name));

            uniprot_type = 'proteome';
            uniprot_ID = '';
            geneIDfield = 'gene_oln';
            reviewed = true;
            if isfield(app.Parameters, 'uniprot')
                if isfield(app.Parameters.uniprot, 'type'), uniprot_type = app.Parameters.uniprot.type; end
                if isfield(app.Parameters.uniprot, 'ID'), uniprot_ID = app.Parameters.uniprot.ID; end
                if isfield(app.Parameters.uniprot, 'geneIDfield'), geneIDfield = app.Parameters.uniprot.geneIDfield; end
                if isfield(app.Parameters.uniprot, 'reviewed'), reviewed = app.Parameters.uniprot.reviewed; end
            end
            f = strrep(f, "obj.params.uniprot.type        = '';", sprintf("obj.params.uniprot.type        = '%s';", uniprot_type));
            f = strrep(f, "obj.params.uniprot.ID          = '';", sprintf("obj.params.uniprot.ID          = '%s';", uniprot_ID));
            f = strrep(f, "obj.params.uniprot.geneIDfield = '';", sprintf("obj.params.uniprot.geneIDfield = '%s';", geneIDfield));
            f = strrep(f, 'obj.params.uniprot.reviewed    = true;', sprintf('obj.params.uniprot.reviewed    = %s;', iif(reviewed, 'true', 'false')));

            taxonomicID = '';
            if isfield(app.Parameters, 'taxonomicID'), taxonomicID = app.Parameters.taxonomicID; end
            f = strrep(f, "obj.params.taxonomicID         = '';", sprintf("obj.params.taxonomicID         = '%s';", taxonomicID));

            % Fill in Core reaction IDs
            c_source = 'EX_glc__D_e';
            if isfield(app.Parameters, 'c_source'), c_source = app.Parameters.c_source; end
            f = strrep(f, "obj.params.c_source = '';", sprintf("obj.params.c_source = '%s';", c_source));

            bioRxn = 'biomass';
            if isfield(app.Parameters, 'bioRxn'), bioRxn = app.Parameters.bioRxn; end
            f = strrep(f, "obj.params.bioRxn   = '';", sprintf("obj.params.bioRxn   = '%s';", bioRxn));

            % Fill in PRESTO options
            PRESTO = struct();
            if isfield(app.Parameters, 'PRESTO'), PRESTO = app.Parameters.PRESTO; end

            runParallel = false;
            if isfield(PRESTO, 'runParallel'), runParallel = PRESTO.runParallel; end
            f = strrep(f, 'obj.params.PRESTO.runParallel = ;', sprintf('obj.params.PRESTO.runParallel = %s;', iif(runParallel, 'true', 'false')));

            ncpu = 4;
            if isfield(PRESTO, 'ncpu'), ncpu = PRESTO.ncpu; end
            f = strrep(f, 'obj.params.PRESTO.ncpu = ;', sprintf('obj.params.PRESTO.ncpu = %d;', ncpu));

            nIter = 100;
            if isfield(PRESTO, 'nIter'), nIter = PRESTO.nIter; end
            f = strrep(f, 'obj.params.PRESTO.nIter = ;', sprintf('obj.params.PRESTO.nIter = %d;', nIter));

            epsilon = 0.1;
            if isfield(PRESTO, 'epsilon'), epsilon = PRESTO.epsilon; end
            f = strrep(f, 'obj.params.PRESTO.epsilon = ;', sprintf('obj.params.PRESTO.epsilon = %g;', epsilon));

            lambda = 0.5;
            if isfield(PRESTO, 'lambda'), lambda = PRESTO.lambda; end
            f = strrep(f, 'obj.params.PRESTO.lambda = ;', sprintf('obj.params.PRESTO.lambda = %g;', lambda));

            theta = 1.0;
            if isfield(PRESTO, 'theta'), theta = PRESTO.theta; end
            f = strrep(f, 'obj.params.PRESTO.theta = ;', sprintf('obj.params.PRESTO.theta = %g;', theta));

            % Save the class file
            filename = fullfile(fullPath, [projectName 'ParameterManagement.m']);
            fid = fopen(filename, 'w');
            fwrite(fid, f);
            fclose(fid);

            app.setStatus(['Generated: ', filename]);
        end

        % -----------------------------------------------------------------
        % STEP 3: Import Model
        % -----------------------------------------------------------------
        function buildStep3_ImportModel(app)
            p = app.ContentPanels{3};

            % Use consistent header (header at y=700, height 65)
            app.buildStepHeader(p, 3, 'Import Model', 'Select model file and configure import settings');

            % Layout: content from y=30 to y=680 (below header at y=700)
            % Row 1: y=520-600 (File Format panel)
            % Row 2: y=430-480 (Model file selection)
            % Row 3: y=370-400 (Checkbox)
            % Row 4: y=30-340 (Info panel)

            % Row 1: File format card
            row1Y = 600;
            fmtPanel = uipanel('Parent', p);
            fmtPanel.Title = 'File Format';
            fmtPanel.FontSize = 10;
            fmtPanel.Units = 'pixels';
            fmtPanel.Position = [30 row1Y 700 55];
            fmtPanel.BackgroundColor = [0.92 0.94 0.98];

            fmtGrp = uibuttongroup('Parent', fmtPanel);
            fmtGrp.Units = 'pixels';
            fmtGrp.Position = [15 12 670 32];
            fmtGrp.BackgroundColor = [0.92 0.94 0.98];

            r1 = uicontrol('Parent', fmtGrp, 'Style', 'radiobutton');
            r1.String = 'SBML (.xml)';
            r1.Units = 'pixels';
            r1.Position = [0 5 100 22];
            r1.Tag = 'xml';
            r1.Value = 1;

            r2 = uicontrol('Parent', fmtGrp, 'Style', 'radiobutton');
            r2.String = 'JSON (.json)';
            r2.Units = 'pixels';
            r2.Position = [110 5 100 22];
            r2.Tag = 'json';

            r3 = uicontrol('Parent', fmtGrp, 'Style', 'radiobutton');
            r3.String = 'YAML (.yml)';
            r3.Units = 'pixels';
            r3.Position = [220 5 80 22];
            r3.Tag = 'yml';

            fmtGrp.SelectedObject = r1;

            % Row 2: Model file selection
            row2Y = 550;
            t = uicontrol('Parent', p, 'Style', 'text');
            t.String = 'Model File:';
            t.FontSize = 10;
            t.HorizontalAlignment = 'left';
            t.Units = 'pixels';
            t.Position = [30 row2Y 80 22];

            modelEdit = uicontrol('Parent', p, 'Style', 'edit');
            modelEdit.Units = 'pixels';
            modelEdit.Position = [120 row2Y 480 26];
            modelEdit.BackgroundColor = [1 1 1];
            modelEdit.Tag = 'modelFile';

            btn = uicontrol('Parent', p, 'Style', 'pushbutton');
            btn.String = 'Browse';
            btn.Units = 'pixels';
            btn.Position = [610 row2Y 80 26];
            btn.BackgroundColor = [0.3 0.3 0.4];
            btn.ForegroundColor = [1 1 1];
            btn.Callback = {@app.browseModelCallback, modelEdit, fmtGrp};

            % Row 3: Checkbox
            row3Y = 520;
            checkInitial = uicontrol('Parent', p, 'Style', 'checkbox');
            checkInitial.String = 'Use as Initial Model (for ParameterManagement)';
            checkInitial.FontSize = 10;
            checkInitial.Units = 'pixels';
            checkInitial.Position = [30 row3Y 400 22];
            checkInitial.Value = 1;
            checkInitial.Tag = 'useAsInitial';

            % Row 4: Info panel
            infoPanelY = 40;
            infoPanelH = 450;
            infoPanel = uipanel('Parent', p);
            infoPanel.Title = 'Model Information';
            infoPanel.FontSize = 10;
            infoPanel.Units = 'pixels';
            infoPanel.Position = [30 infoPanelY 700 infoPanelH];
            infoPanel.BackgroundColor = [0.92 0.94 0.98];

            infoBox = uicontrol('Parent', infoPanel, 'Style', 'edit');
            infoBox.Units = 'pixels';
            infoBox.Position = [10 0 680 infoPanelH-20];
            infoBox.FontSize = 9;
            infoBox.FontName = 'Courier New';
            infoBox.BackgroundColor = [0.92 0.94 0.98];
            infoBox.HorizontalAlignment = 'left';
            infoBox.Max = 10;
            infoBox.String = {'Model information will appear here after loading...', '', 'Tip: For eciML1515, use iML1515.xml from BiGG database'};
            infoBox.Tag = 'infoBox';

            % Load button
            btn = uicontrol('Parent', p, 'Style', 'pushbutton');
            btn.String = 'Load Model';
            btn.Units = 'pixels';
            btn.Position = [30 30 130 40];
            btn.BackgroundColor = [0.2 0.5 0.4];
            btn.ForegroundColor = [1 1 1];
            btn.FontSize = 11;
            btn.Callback = {@app.loadModelCallback, modelEdit, infoBox, checkInitial};

            % Navigation buttons
            app.NavButtons.prevBtn = uicontrol('Parent', p, 'Style', 'pushbutton');
            app.NavButtons.prevBtn.String = 'Previous';
            app.NavButtons.prevBtn.Units = 'pixels';
            app.NavButtons.prevBtn.Position = [800 30 100 35];
            app.NavButtons.prevBtn.BackgroundColor = [0.5 0.5 0.55];
            app.NavButtons.prevBtn.Callback = @(~,~) app.updateStep(2);

            app.NavButtons.nextBtn = uicontrol('Parent', p, 'Style', 'pushbutton');
            app.NavButtons.nextBtn.String = 'Next';
            app.NavButtons.nextBtn.Units = 'pixels';
            app.NavButtons.nextBtn.Position = [910 30 100 35];
            app.NavButtons.nextBtn.BackgroundColor = [0.2 0.4 0.6];
            app.NavButtons.nextBtn.ForegroundColor = [1 1 1];
            app.NavButtons.nextBtn.Enable = iif(app.StepCompleted(3), 'on', 'off');
            app.NavButtons.nextBtn.Callback = @(~,~) app.updateStep(4);
        end

        function browseModelCallback(app, hObj, event, modelEdit, fmtGrp)
            if isempty(fmtGrp.SelectedObject)
                tag = 'xml';
            else
                tag = fmtGrp.SelectedObject.Tag;
            end

            switch tag
                case 'xml'
                    filters = {'*.xml', 'SBML Files (*.xml)'};
                case 'json'
                    filters = {'*.json', 'JSON Files (*.json)'};
                case 'yml'
                    filters = {'*.yml;*.yaml', 'YAML Files (*.yml, *.yaml)'};
            end

            % Use model's directory if available, otherwise use pwd
            if isfield(app.Parameters, 'modelDir') && isfolder(app.Parameters.modelDir)
                startDir = app.Parameters.modelDir;
            else
                startDir = pwd;
            end

            [file, path] = uigetfile(filters, 'Select Model File', startDir);
            if file ~= 0
                modelEdit.String = fullfile(path, file);
            end
        end

        function loadModelCallback(app, hObj, event, modelEdit, infoBox, checkInitial)
            filename = modelEdit.String;

            if isempty(filename)
                msgbox('Please select a model file.', 'Error', 'warn');
                return;
            end

            app.setStatus(['Loading: ', filename]);
            drawnow;

            try
                if isempty(app.Parameters)
                    app.Parameters = struct();
                end

                modelType = 'TRADITION';
                if isfield(app.Parameters, 'modeltype')
                    modelType = app.Parameters.modeltype;
                end

                model = loadModel(filename, modelType, fileparts(filename), app.Parameters);
                app.Model = model;

                if checkInitial.Value
                    app.Parameters.InitialModel = filename;
                end

                % Copy model file to project's models directory
                if isfield(app.Parameters, 'modelDir') && isfolder(app.Parameters.modelDir)
                    [~, name, ext] = fileparts(filename);
                    destFile = fullfile(app.Parameters.modelDir, [name ext]);
                    if ~strcmp(filename, destFile)
                        copyfile(filename, destFile);
                        app.setStatus(['Model copied to: ', destFile]);
                    end
                end

                info = sprintf('Model loaded successfully!\n\n');
                info = [info, sprintf('ID: %s\n', model.id)];
                info = [info, sprintf('Type: %s\n', model.type)];
                info = [info, sprintf('Reactions: %d\n', numel(model.rxns))];
                info = [info, sprintf('Metabolites: %d\n', numel(model.mets))];
                info = [info, sprintf('Genes: %d', numel(model.genes))];

                infoBox.String = splitlines(info);
                app.setStatus('Model loaded');

                % Mark step as completed and enable Next button
                app.StepCompleted(3) = true;
                app.updateStepListDisplay();
                app.NavButtons.nextBtn.Enable = 'on';

                msgbox({info, '', 'Click "Next" to continue.'}, 'Success', 'help');

            catch ME
                app.setStatus('Load failed');
                msgbox(ME.message, 'Error', 'error');
            end
        end
        % -----------------------------------------------------------------
        function buildStep4_Conversion(app)
            p = app.ContentPanels{4};

            % Use consistent header
            app.buildStepHeader(p, 4, 'Model Conversion', 'Convert model to enzyme-constraint format');

            contentY = 650;

            % ecModel type selection - left aligned
            t = uicontrol('Parent', p, 'Style', 'text');
            t.String = 'ecModel Type:';
            t.FontSize = 10;
            t.HorizontalAlignment = 'left';
            t.Units = 'pixels';
            t.Position = [30 contentY 100 25];

            ecTypePopup = uicontrol('Parent', p, 'Style', 'popupmenu');
            ecTypePopup.String = 'basic|isozyme|integrated';
            ecTypePopup.Units = 'pixels';
            ecTypePopup.Position = [140 contentY 150 28];
            ecTypePopup.BackgroundColor = [1 1 1];
            ecTypePopup.Tag = 'ecTypePopup';
            if isfield(app.Parameters, 'ecModelType')
                switch app.Parameters.ecModelType
                    case 'basic', ecTypePopup.Value = 1;
                    case 'isozyme', ecTypePopup.Value = 2;
                    case 'integrated', ecTypePopup.Value = 3;
                end
            else
                ecTypePopup.Value = 3;
            end

            % Parameters display - left aligned card
            contentY = contentY - 30;
            sigma = 0.5; Ptot = 0.55; f = 0.55;
            if isfield(app.Parameters, 'sigma'), sigma = app.Parameters.sigma; end
            if isfield(app.Parameters, 'Ptot'), Ptot = app.Parameters.Ptot; end
            if isfield(app.Parameters, 'f'), f = app.Parameters.f; end

            paramPanel = uipanel('Parent', p);
            paramPanel.Title = 'Conversion Parameters';
            paramPanel.FontSize = 10;
            paramPanel.Units = 'pixels';
            paramPanel.Position = [30 contentY-120 500 150];
            paramPanel.BackgroundColor = [0.92 0.94 0.98];

            t = uicontrol('Parent', paramPanel, 'Style', 'text');
            t.Position = [15 15 470 120];
            t.FontSize = 10;
            t.HorizontalAlignment = 'left';
            t.Units = 'pixels';
            t.String = {sprintf('sigma = %.3f (enzyme saturation)', sigma), '', ...
                       sprintf('Ptot = %.3f (total protein, g/gDCW)', Ptot), '', ...
                       sprintf('f = %.3f (enzyme fraction)', f), '', ...
                       sprintf('Protein Budget = %.4f g/gDCW', sigma * Ptot * f * 1000)};

            % Run button - left aligned
            btn = uicontrol('Parent', p, 'Style', 'pushbutton');
            btn.String = 'Run Conversion';
            btn.Units = 'pixels';
            btn.Position = [30 30 150 40];
            btn.BackgroundColor = [0.2 0.5 0.4];
            btn.ForegroundColor = [1 1 1];
            btn.FontSize = 11;
            btn.Callback = @app.runConversionCallback;

            % Progress label
            progressLabel = uicontrol('Parent', p, 'Style', 'text');
            progressLabel.String = '';
            progressLabel.FontSize = 10;
            progressLabel.Units = 'pixels';
            progressLabel.Position = [195 38 300 25];
            progressLabel.Tag = 'progressLabel';

            % Results panel - left aligned
            resultsPanel = uipanel('Parent', p);
            resultsPanel.Title = 'Conversion Results';
            resultsPanel.FontSize = 10;
            resultsPanel.Units = 'pixels';
            resultsPanel.Position = [30 210 700 280];
            resultsPanel.BackgroundColor = [0.92 0.94 0.98];

            resultsText = uicontrol('Parent', resultsPanel, 'Style', 'text');
            resultsText.Position = [10 0 680 260];
            resultsText.FontSize = 9;
            resultsText.FontName = 'Courier New';
            resultsText.HorizontalAlignment = 'left';
            resultsText.Units = 'pixels';
            resultsText.String = {'Run conversion to see results'};
            resultsText.Tag = 'resultsText';

            % Navigation buttons - consistent bottom position
            app.NavButtons.prevBtn = uicontrol('Parent', p, 'Style', 'pushbutton');
            app.NavButtons.prevBtn.String = 'Previous';
            app.NavButtons.prevBtn.Units = 'pixels';
            app.NavButtons.prevBtn.Position = [800 30 100 35];
            app.NavButtons.prevBtn.BackgroundColor = [0.5 0.5 0.55];
            app.NavButtons.prevBtn.Callback = @(~,~) app.updateStep(3);

            app.NavButtons.nextBtn = uicontrol('Parent', p, 'Style', 'pushbutton');
            app.NavButtons.nextBtn.String = 'Next';
            app.NavButtons.nextBtn.Units = 'pixels';
            app.NavButtons.nextBtn.Position = [910 30 100 35];
            app.NavButtons.nextBtn.BackgroundColor = [0.2 0.4 0.6];
            app.NavButtons.nextBtn.ForegroundColor = [1 1 1];
            app.NavButtons.nextBtn.Enable = iif(app.StepCompleted(4), 'on', 'off');
            app.NavButtons.nextBtn.Callback = @(~,~) app.updateStep(5);
        end

        function runConversionCallback(app, hObj, event)
            if isempty(app.Model)
                msgbox('Please load a model first.', 'Error', 'warn');
                return;
            end

            progressLabel = findobj(app.ContentPanels{4}, 'Tag', 'progressLabel');
            resultsText = findobj(app.ContentPanels{4}, 'Tag', 'resultsText');
            ecTypePopup = findobj(app.ContentPanels{4}, 'Tag', 'ecTypePopup');

            if ~isempty(progressLabel), progressLabel.String = 'Converting...'; end
            drawnow;

            try
                app.setStatus('Converting model...');

                % Get ecModel type from popupmenu selection
                ecTypes = {'basic', 'isozyme', 'integrated'};
                ecModelType = ecTypes{ecTypePopup.Value};

                ecModel = convertecModel(app.Model, ecModelType, app.Parameters);
                app.Model = ecModel;

                if ~isempty(progressLabel), progressLabel.String = 'Conversion complete!'; end

                % Show results
                m = app.Model;
                txt = {sprintf('Model ID: %s', m.id), ...
                       sprintf('Model Type: %s', m.type), '', ...
                       sprintf('Original Reactions: %d', numel(m.rxns)), ...
                       sprintf('Original Metabolites: %d', numel(m.mets)), ...
                       sprintf('Original Genes: %d', numel(m.genes)), ''};
                if isfield(m, 'enzymeConstraints')
                    ec = m.enzymeConstraints;
                    txt{end+1} = sprintf('EC Reactions: %d', numel(ec.rxns));
                    txt{end+1} = sprintf('EC Genes: %d', numel(ec.genes));
                    if isfield(ec, 'sigma')
                        txt{end+1} = sprintf('Sigma: %.3f', ec.sigma);
                    end
                end

                if ~isempty(resultsText), resultsText.String = txt; end
                app.setStatus('Model converted');

                % Mark step as completed and enable Next button
                app.StepCompleted(4) = true;
                app.updateStepListDisplay();
                app.NavButtons.nextBtn.Enable = 'on';

                msgbox('Model converted successfully! Click "Next" to continue.', 'Success', 'help');

            catch ME
                if ~isempty(progressLabel), progressLabel.String = 'Failed'; end
                app.setStatus('Conversion failed');
                msgbox(ME.message, 'Error', 'error');
            end
        end

        % -----------------------------------------------------------------
        % STEP 5: Fill ecModel Info
        % -----------------------------------------------------------------
        function buildStep5_FillInfo(app)
            p = app.ContentPanels{5};

            % Use consistent header
            app.buildStepHeader(p, 5, 'Fill ecModel Information', 'Fill complex info, EC numbers, and SMILES/InChIKey data');

            % Layout: Left side has 3 panels stacked, Right side has results panel below them
            panelWidth = 500;
            panelHeight = 100;
            panelGap = 15;
            startY = 620;

            % === Section 1: Complex Enzyme Information ===
            p1Y = startY - 50;
            complexPanel = uipanel('Parent', p);
            complexPanel.Title = '1. Complex Enzyme Information';
            complexPanel.FontSize = 10;
            complexPanel.Units = 'pixels';
            complexPanel.Position = [30 p1Y panelWidth panelHeight];
            complexPanel.BackgroundColor = [0.92 0.94 0.98];

            t = uicontrol('Parent', complexPanel, 'Style', 'text');
            t.String = 'Apply complex enzyme data (from ComplexPortal)';
            t.FontSize = 9;
            t.Units = 'pixels';
            t.Position = [15 55 470 20];

            btnComplex = uicontrol('Parent', complexPanel, 'Style', 'pushbutton');
            btnComplex.String = 'Run';
            btnComplex.Units = 'pixels';
            btnComplex.Position = [15 15 80 30];
            btnComplex.BackgroundColor = [0.3 0.5 0.7];
            btnComplex.ForegroundColor = [1 1 1];
            btnComplex.Callback = @app.runComplexCallback;

            t = uicontrol('Parent', complexPanel, 'Style', 'text');
            t.String = 'Status:';
            t.FontSize = 9;
            t.Units = 'pixels';
            t.Position = [110 18 50 18];

            complexStatus = uicontrol('Parent', complexPanel, 'Style', 'text');
            complexStatus.String = 'Not run';
            complexStatus.FontSize = 9;
            complexStatus.ForegroundColor = [0.5 0.5 0.5];
            complexStatus.Units = 'pixels';
            complexStatus.Position = [160 18 320 18];
            complexStatus.Tag = 'complexStatus';

            % === Section 2: EC Number Data ===
            p2Y = p1Y - panelHeight - panelGap;
            ecPanel = uipanel('Parent', p);
            ecPanel.Title = '2. EC Number Data';
            ecPanel.FontSize = 10;
            ecPanel.Units = 'pixels';
            ecPanel.Position = [30 p2Y panelWidth 125];
            ecPanel.BackgroundColor = [0.92 0.94 0.98];

            t = uicontrol('Parent', ecPanel, 'Style', 'text');
            t.String = 'Fill EC numbers from UniProt database';
            t.FontSize = 9;
            t.Units = 'pixels';
            t.Position = [15 75 470 20];

            checkEC = uicontrol('Parent', ecPanel, 'Style', 'checkbox');
            checkEC.String = 'Fill EC numbers';
            checkEC.FontSize = 9;
            checkEC.Units = 'pixels';
            checkEC.Position = [15 45 150 18];
            checkEC.Value = 1;
            checkEC.Tag = 'checkEC';

            btnEC = uicontrol('Parent', ecPanel, 'Style', 'pushbutton');
            btnEC.String = 'Run';
            btnEC.Units = 'pixels';
            btnEC.Position = [15 5 80 30];
            btnEC.BackgroundColor = [0.3 0.5 0.7];
            btnEC.ForegroundColor = [1 1 1];
            btnEC.Callback = {@app.runECnumberCallback, checkEC};

            t = uicontrol('Parent', ecPanel, 'Style', 'text');
            t.String = 'Status:';
            t.FontSize = 9;
            t.Units = 'pixels';
            t.Position = [110 8 50 18];

            ecStatus = uicontrol('Parent', ecPanel, 'Style', 'text');
            ecStatus.String = 'Not run';
            ecStatus.FontSize = 9;
            ecStatus.ForegroundColor = [0.5 0.5 0.5];
            ecStatus.Units = 'pixels';
            ecStatus.Position = [160 8 320 18];
            ecStatus.Tag = 'ecStatus';

            % === Section 3: SMILES/InChIKey Data ===
            p3Y = p2Y - panelHeight - panelGap;
            smilesPanel = uipanel('Parent', p);
            smilesPanel.Title = '3. SMILES/InChIKey Data';
            smilesPanel.FontSize = 10;
            smilesPanel.Units = 'pixels';
            smilesPanel.Position = [30 p3Y panelWidth 125];
            smilesPanel.BackgroundColor = [0.92 0.94 0.98];

            t = uicontrol('Parent', smilesPanel, 'Style', 'text');
            t.String = 'Fill metabolite SMILES/InChIKey (PubChem/MetaNetX)';
            t.FontSize = 9;
            t.Units = 'pixels';
            t.Position = [15 75 470 20];

            checkSmiles = uicontrol('Parent', smilesPanel, 'Style', 'checkbox');
            checkSmiles.String = 'Fill SMILES/InChIKey';
            checkSmiles.FontSize = 9;
            checkSmiles.Units = 'pixels';
            checkSmiles.Position = [15 45 150 18];
            checkSmiles.Value = 1;
            checkSmiles.Tag = 'checkSmiles';

            btnSmiles = uicontrol('Parent', smilesPanel, 'Style', 'pushbutton');
            btnSmiles.String = 'Run';
            btnSmiles.Units = 'pixels';
            btnSmiles.Position = [15 5 80 30];
            btnSmiles.BackgroundColor = [0.3 0.5 0.7];
            btnSmiles.ForegroundColor = [1 1 1];
            btnSmiles.Callback = {@app.runSmilesCallback, checkSmiles};

            t = uicontrol('Parent', smilesPanel, 'Style', 'text');
            t.String = 'Status:';
            t.FontSize = 9;
            t.Units = 'pixels';
            t.Position = [110 8 50 18];

            smilesStatus = uicontrol('Parent', smilesPanel, 'Style', 'text');
            smilesStatus.String = 'Not run';
            smilesStatus.FontSize = 9;
            smilesStatus.ForegroundColor = [0.5 0.5 0.5];
            smilesStatus.Units = 'pixels';
            smilesStatus.Position = [160 8 320 18];
            smilesStatus.Tag = 'smilesStatus';

            % Results summary - Below the three panels on the left
            resultsY = p3Y - panelHeight - panelGap;
            resultsPanel = uipanel('Parent', p);
            resultsPanel.Title = 'Results Summary';
            resultsPanel.FontSize = 10;
            resultsPanel.Units = 'pixels';
            resultsPanel.Position = [560 340 500 330];
            resultsPanel.BackgroundColor = [0.92 0.94 0.98];

            resultsText = uicontrol('Parent', resultsPanel, 'Style', 'text');
            resultsText.Units = 'pixels';
            resultsText.Position = [10 0 500 310];
            resultsText.FontSize = 9;
            resultsText.FontName = 'Courier New';
            resultsText.HorizontalAlignment = 'left';
            resultsText.String = {'Step 5: Fill ecModel Information', '', ...
                '1. Complex Enzyme Information: From ComplexPortal/RHAMN database', ...
                '2. EC Number Data: From UniProt database', ...
                '3. SMILES/InChIKey Data: From PubChem/MetaNetX', '', ...
                'Click each section button to run that step.'};
            resultsText.Tag = 'step5Results';

            % Navigation buttons - consistent bottom position
            app.NavButtons.prevBtn = uicontrol('Parent', p, 'Style', 'pushbutton');
            app.NavButtons.prevBtn.String = 'Previous';
            app.NavButtons.prevBtn.Units = 'pixels';
            app.NavButtons.prevBtn.Position = [800 30 100 35];
            app.NavButtons.prevBtn.BackgroundColor = [0.5 0.5 0.55];
            app.NavButtons.prevBtn.Callback = @(~,~) app.updateStep(4);

            app.NavButtons.nextBtn = uicontrol('Parent', p, 'Style', 'pushbutton');
            app.NavButtons.nextBtn.String = 'Next';
            app.NavButtons.nextBtn.Units = 'pixels';
            app.NavButtons.nextBtn.Position = [910 30 100 35];
            app.NavButtons.nextBtn.BackgroundColor = [0.2 0.4 0.6];
            app.NavButtons.nextBtn.ForegroundColor = [1 1 1];
            app.NavButtons.nextBtn.Enable = iif(app.StepCompleted(5), 'on', 'off');
            app.NavButtons.nextBtn.Callback = @(~,~) app.updateStep(6);
        end

        function runComplexCallback(app, hObj, event, checkComplex)
            if isempty(app.Model)
                msgbox('Please convert a model first.', 'Error', 'warn');
                return;
            end

            complexStatus = findobj(app.ContentPanels{5}, 'Style', 'text', 'Tag', 'complexStatus');
            resultsText = findobj(app.ContentPanels{5}, 'Style', 'text', 'Tag', 'step5Results');
            if numel(complexStatus) > 1, complexStatus = complexStatus(1); end
            if numel(resultsText) > 1, resultsText = resultsText(1); end

            try
                complexStatus.String = 'Running...';
                complexStatus.ForegroundColor = [0.8 0.5 0.2];
                drawnow;

                % Get parameters - use ParameterManager if app.Parameters is incomplete
                if ~isempty(app.Parameters) && isfield(app.Parameters, 'dataDir') && ~isempty(app.Parameters.dataDir)
                    params = app.Parameters;
                else
                    % Try to get from ParameterManager
                    try
                        params = ParameterManager.getParams();
                    catch ME2
                        error('Please complete Step 2 (Configure Parameters) first to initialize parameters. Error: %s', ME2.message);
                    end
                end

                % Verify dataDir exists
                if ~isfield(params, 'dataDir') || isempty(params.dataDir)
                    error('dataDir is not set. Please check that Step 2 was completed correctly.');
                end

                % Debug: show actual paths being used
                app.setStatus(sprintf('dataDir: %s', params.dataDir));

                % Ensure data directory exists
                if ~exist(params.dataDir, 'dir')
                    error('Data directory does not exist: %s. Please complete Step 1 (Initialize Project) first.', params.dataDir);
                end

                % Get complex data
                % complexInfo = getComplexdata(params.taxonomicID, params);
                complexInfo = getComplexdata();

                % Apply complex data directly with third parameter = false (no proposed)
                [model, foundComplex, ~, ~] = applyComplexdata(app.Model, complexInfo, false);
                app.Model = model;

                nFound = height(foundComplex);

                complexStatus.String = sprintf('Done! Found: %d', nFound);
                complexStatus.ForegroundColor = [0.2 0.5 0.3];

                txt = {sprintf('Complex Info:'), '', ...
                    sprintf('  Found complexes: %d', nFound), ''};
                if ~isempty(resultsText)
                    oldTxt = resultsText.String;
                    resultsText.String = [oldTxt; {''}; txt(:)];
                end

                % Mark complex section as done
                complexStatus.Tag = 'complexDone';

                % Check if all sections done
                app.checkStep5Completion();

            catch ME
                complexStatus.String = 'Failed';
                complexStatus.ForegroundColor = [0.8 0.2 0.2];
                msgbox(ME.message, 'Complex Info Error', 'error');
            end
        end

        function runECnumberCallback(app, hObj, event, checkEC)
            if isempty(app.Model)
                msgbox('Please convert a model first.', 'Error', 'warn');
                return;
            end

            ecStatus = findobj(app.ContentPanels{5}, 'Style', 'text', 'Tag', 'ecStatus');
            resultsText = findobj(app.ContentPanels{5}, 'Style', 'text', 'Tag', 'step5Results');
            if numel(ecStatus) > 1, ecStatus = ecStatus(1); end
            if numel(resultsText) > 1, resultsText = resultsText(1); end

            try
                ecStatus.String = 'Running...';
                ecStatus.ForegroundColor = [0.8 0.5 0.2];
                drawnow;

                % Get parameters - use ParameterManager if app.Parameters is incomplete
                if ~isempty(app.Parameters) && isfield(app.Parameters, 'dataDir') && ~isempty(app.Parameters.dataDir)
                    params = app.Parameters;
                else
                    % Try to get from ParameterManager
                    try
                        params = ParameterManager.getParams();
                    catch ME2
                        error('Please complete Step 2 (Configure Parameters) first to initialize parameters. Error: %s', ME2.message);
                    end
                end

                % Verify dataDir exists
                if ~isfield(params, 'dataDir') || isempty(params.dataDir)
                    error('dataDir is not set. Please check that Step 2 was completed correctly.');
                end

                % Fill EC numbers
                model = getECnumber(app.Model, params);
                app.Model = model;

                ecStatus.String = 'Done!';
                ecStatus.ForegroundColor = [0.2 0.5 0.3];

                txt = {sprintf('EC Numbers:'), '', ...
                    sprintf('  EC numbers filled from UniProt'), ''};
                if ~isempty(resultsText)
                    oldTxt = resultsText.String;
                    resultsText.String = [oldTxt; {''}; txt(:)];
                end

                % Check if all sections done
                app.checkStep5Completion();

            catch ME
                ecStatus.String = 'Failed';
                ecStatus.ForegroundColor = [0.8 0.2 0.2];
                msgbox(ME.message, 'EC Number Error', 'error');
            end
        end

        function runSmilesCallback(app, hObj, event, checkSmiles)
            if isempty(app.Model)
                msgbox('Please convert a model first.', 'Error', 'warn');
                return;
            end

            smilesStatus = findobj(app.ContentPanels{5}, 'Style', 'text', 'Tag', 'smilesStatus');
            resultsText = findobj(app.ContentPanels{5}, 'Style', 'text', 'Tag', 'step5Results');
            if numel(smilesStatus) > 1, smilesStatus = smilesStatus(1); end
            if numel(resultsText) > 1, resultsText = resultsText(1); end

            try
                smilesStatus.String = 'Running...';
                smilesStatus.ForegroundColor = [0.8 0.5 0.2];
                drawnow;

                % Get parameters - use ParameterManager if app.Parameters is incomplete
                if ~isempty(app.Parameters) && isfield(app.Parameters, 'dataDir') && ~isempty(app.Parameters.dataDir)
                    params = app.Parameters;
                else
                    % Try to get from ParameterManager
                    try
                        params = ParameterManager.getParams();
                    catch ME2
                        error('Please complete Step 2 (Configure Parameters) first to initialize parameters. Error: %s', ME2.message);
                    end
                end

                % Verify dataDir exists
                if ~isfield(params, 'dataDir') || isempty(params.dataDir)
                    error('dataDir is not set. Please check that Step 2 was completed correctly.');
                end

                % Fill SMILES/InChIKey
                [model, noSMILES, noInChIKey] = getMetinfo(app.Model, 1, 'AB', params);

                % Add MetaNetX IDs
                [model, noMNX] = addMetMetaNetXID(model, 0);

                app.Model = model;

                nSMILES = sum(~cellfun(@isempty, model.metSmiles));
                nInChIKey = sum(~cellfun(@isempty, model.metInChIKey));
                nMNX = sum(~cellfun(@isempty, model.metMetaNetXID));

                smilesStatus.String = sprintf('Done! SMILES: %d, InChIKey: %d, MNX: %d', nSMILES, nInChIKey, nMNX);
                smilesStatus.ForegroundColor = [0.2 0.5 0.3];

                txt = {sprintf('SMILES/InChIKey/MetaNetX:'), '', ...
                    sprintf('  Metabolites with SMILES: %d', nSMILES), ...
                    sprintf('  Metabolites with InChIKey: %d', nInChIKey), ...
                    sprintf('  Metabolites with MNX ID: %d', nMNX), ...
                    sprintf('  Missing SMILES: %d', numel(noSMILES)), ...
                    sprintf('  Missing InChIKey: %d', numel(noInChIKey)), ...
                    sprintf('  Missing MNX ID: %d', numel(noMNX)), ''};
                if ~isempty(resultsText)
                    oldTxt = resultsText.String;
                    resultsText.String = [oldTxt; {''}; txt(:)];
                end

                % Check if all sections done
                app.checkStep5Completion();

            catch ME
                smilesStatus.String = 'Failed';
                smilesStatus.ForegroundColor = [0.8 0.2 0.2];
                msgbox(ME.message, 'SMILES/InChIKey Error', 'error');
            end
        end

        function checkStep5Completion(app)
            % Enable next button when all 3 sections are done
            complexStatus = findobj(app.ContentPanels{5}, 'Style', 'text', 'Tag', 'complexStatus');
            ecStatus = findobj(app.ContentPanels{5}, 'Style', 'text', 'Tag', 'ecStatus');
            smilesStatus = findobj(app.ContentPanels{5}, 'Style', 'text', 'Tag', 'smilesStatus');
            if numel(complexStatus) > 1, complexStatus = complexStatus(1); end
            if numel(ecStatus) > 1, ecStatus = ecStatus(1); end
            if numel(smilesStatus) > 1, smilesStatus = smilesStatus(1); end

            allDone = true;
            if ~isempty(complexStatus) && ~strcmp(complexStatus.String(1:3), 'Don')
                allDone = false;
            end
            if ~isempty(ecStatus) && ~strcmp(ecStatus.String(1:3), 'Don')
                allDone = false;
            end
            if ~isempty(smilesStatus) && ~strcmp(smilesStatus.String(1:3), 'Don')
                allDone = false;
            end

            if allDone
                app.StepCompleted(5) = true;
                app.updateStepListDisplay();
                app.NavButtons.nextBtn.Enable = 'on';
                msgbox('Step 5 completed! Click "Next" to continue.', 'Success', 'help');
            end
        end

        % -----------------------------------------------------------------
        % STEP 6: DL kcat Prediction
        % -----------------------------------------------------------------
        function buildStep6_DLPrediction(app)
            p = app.ContentPanels{6};

            % Use helper to build consistent header
            app.buildStepHeader(p, 6, 'Deep Learning kcat Prediction', ...
                'Generate DL input file, run Docker-based prediction, and select kcat values');

            % Layout: 3 panels in a row, status below each, results at bottom
            panelWidth = 320;
            panelHeight = 180;
            panelY = 490;
            panelGap = 20;

            % === Section 1: Generate DL Input File ===
            inputPanel = uipanel('Parent', p);
            inputPanel.Title = '1. Generate Input';
            inputPanel.FontSize = 10;
            inputPanel.Units = 'pixels';
            inputPanel.Position = [30 panelY panelWidth panelHeight];
            inputPanel.BackgroundColor = [0.92 0.94 0.98];

            t = uicontrol('Parent', inputPanel, 'Style', 'text');
            t.String = 'Select model:';
            t.FontSize = 9;
            t.Units = 'pixels';
            t.Position = [15 135 80 18];

            toolPopup1 = uicontrol('Parent', inputPanel, 'Style', 'popupmenu');
            toolPopup1.String = 'CatPred|DLKcat|UniKP';
            toolPopup1.Value = 1;
            toolPopup1.Units = 'pixels';
            toolPopup1.Position = [15 108 150 24];
            toolPopup1.Tag = 'dlToolPopup1';
            toolPopup1.BackgroundColor = [1 1 1];

            btnGenInput = uicontrol('Parent', inputPanel, 'Style', 'pushbutton');
            btnGenInput.String = 'Generate Input';
            btnGenInput.Units = 'pixels';
            btnGenInput.Position = [15 60 290 40];
            btnGenInput.BackgroundColor = [0.3 0.5 0.7];
            btnGenInput.ForegroundColor = [1 1 1];
            btnGenInput.FontSize = 10;
            btnGenInput.Callback = @(~,~) app.generateDLInput(toolPopup1);

            inputStatus = uicontrol('Parent', inputPanel, 'Style', 'text');
            inputStatus.String = 'Status: Not run';
            inputStatus.FontSize = 9;
            inputStatus.ForegroundColor = [0.5 0.5 0.5];
            inputStatus.Units = 'pixels';
            inputStatus.Position = [15 25 290 18];
            inputStatus.Tag = 'inputStatus';

            % === Section 2: Run DL Prediction ===
            predictPanel = uipanel('Parent', p);
            predictPanel.Title = '2. Run Prediction';
            predictPanel.FontSize = 10;
            predictPanel.Units = 'pixels';
            predictPanel.Position = [30 + panelWidth + panelGap panelY panelWidth panelHeight];
            predictPanel.BackgroundColor = [0.92 0.94 0.98];

            t = uicontrol('Parent', predictPanel, 'Style', 'text');
            t.String = 'Requires Docker running';
            t.FontSize = 9;
            t.ForegroundColor = [0.6 0.4 0.4];
            t.Units = 'pixels';
            t.Position = [15 135 290 18];

            t = uicontrol('Parent', predictPanel, 'Style', 'text');
            t.String = 'Select model:';
            t.FontSize = 9;
            t.Units = 'pixels';
            t.Position = [15 108 80 18];

            toolPopup2 = uicontrol('Parent', predictPanel, 'Style', 'popupmenu');
            toolPopup2.String = 'CatPred|DLKcat|UniKP';
            toolPopup2.Value = 1;
            toolPopup2.Units = 'pixels';
            toolPopup2.Position = [15 82 150 24];
            toolPopup2.Tag = 'dlToolPopup2';
            toolPopup2.BackgroundColor = [1 1 1];

            btnRun = uicontrol('Parent', predictPanel, 'Style', 'pushbutton');
            btnRun.String = 'Run Prediction';
            btnRun.Units = 'pixels';
            btnRun.Position = [15 35 290 40];
            btnRun.BackgroundColor = [0.2 0.5 0.4];
            btnRun.ForegroundColor = [1 1 1];
            btnRun.FontSize = 10;
            btnRun.Callback = @(~,~) app.runDLPrediction(toolPopup2);

            % === Section 3: Select kcat Values ===
            selectPanel = uipanel('Parent', p);
            selectPanel.Title = '3. Select kcat';
            selectPanel.FontSize = 10;
            selectPanel.Units = 'pixels';
            selectPanel.Position = [30 + (panelWidth + panelGap)*2 panelY panelWidth panelHeight];
            selectPanel.BackgroundColor = [0.92 0.94 0.98];

            t = uicontrol('Parent', selectPanel, 'Style', 'text');
            t.String = 'Load kcat prediction results';
            t.FontSize = 9;
            t.ForegroundColor = [0.4 0.4 0.4];
            t.Units = 'pixels';
            t.Position = [15 135 290 18];

            btnSelect = uicontrol('Parent', selectPanel, 'Style', 'pushbutton');
            btnSelect.String = 'Select kcat File';
            btnSelect.Units = 'pixels';
            btnSelect.Position = [15 35 290 40];
            btnSelect.BackgroundColor = [0.5 0.5 0.3];
            btnSelect.ForegroundColor = [1 1 1];
            btnSelect.FontSize = 10;
            btnSelect.Callback = @app.selectKcatCallback;

            selectStatus = uicontrol('Parent', selectPanel, 'Style', 'text');
            selectStatus.String = 'Status: No file';
            selectStatus.FontSize = 9;
            selectStatus.ForegroundColor = [0.5 0.5 0.5];
            selectStatus.Units = 'pixels';
            selectStatus.Position = [15 10 290 18];
            selectStatus.Tag = 'selectStatus';

            % Docker check button - below section 2
            dockerStatus = uicontrol('Parent', p, 'Style', 'text');
            dockerStatus.String = 'Docker: Not checked';
            dockerStatus.FontSize = 10;
            dockerStatus.ForegroundColor = [0.5 0.5 0.5];
            dockerStatus.Units = 'pixels';
            dockerStatus.Position = [130 + panelWidth + panelGap - 200 438 400 40];
            dockerStatus.Tag = 'dockerStatus6';

            btnCheck = uicontrol('Parent', p, 'Style', 'pushbutton');
            btnCheck.String = 'Check Docker';
            btnCheck.Units = 'pixels';
            btnCheck.Position = [30 + panelWidth + panelGap + 315 448 100 24];
            btnCheck.BackgroundColor = [0.5 0.5 0.5];
            btnCheck.ForegroundColor = [1 1 1];
            btnCheck.Callback = @(~,~) app.checkDockerStatus(dockerStatus);

            % Results Summary Panel - below the three sections
            resultsY = panelY - panelHeight - 30 - 100;
            resultsPanel = uipanel('Parent', p);
            resultsPanel.Title = 'Results & Information';
            resultsPanel.FontSize = 10;
            resultsPanel.Units = 'pixels';
            resultsPanel.Position = [30 resultsY 1000 250];
            resultsPanel.BackgroundColor = [0.92 0.94 0.98];

            resultsText = uicontrol('Parent', resultsPanel, 'Style', 'text');
            resultsText.Position = [15 0 970 250];
            resultsText.FontSize = 9;
            resultsText.FontName = 'Courier New';
            resultsText.HorizontalAlignment = 'left';
            resultsText.Units = 'pixels';
            resultsText.String = {'Instructions:', '', ...
                '1. Generate DL Input File: Select CatPred/DLKcat/UniKP, click Generate Input', ...
                '2. Run DL Prediction: Ensure Docker is running, select same model, click Run', ...
                '3. Select kcat Values: Load kcat prediction results for Step 7', '', ...
                'Note: Each model can be run independently.'};
            resultsText.Tag = 'step6Results';

            % Navigation buttons - consistent bottom position
            app.NavButtons.prevBtn = uicontrol('Parent', p, 'Style', 'pushbutton');
            app.NavButtons.prevBtn.String = 'Previous';
            app.NavButtons.prevBtn.Units = 'pixels';
            app.NavButtons.prevBtn.Position = [800 30 100 35];
            app.NavButtons.prevBtn.BackgroundColor = [0.5 0.5 0.55];
            app.NavButtons.prevBtn.Callback = @(~,~) app.updateStep(5);

            app.NavButtons.nextBtn = uicontrol('Parent', p, 'Style', 'pushbutton');
            app.NavButtons.nextBtn.String = 'Next';
            app.NavButtons.nextBtn.Units = 'pixels';
            app.NavButtons.nextBtn.Position = [910 30 100 35];
            app.NavButtons.nextBtn.BackgroundColor = [0.2 0.4 0.6];
            app.NavButtons.nextBtn.ForegroundColor = [1 1 1];
            app.NavButtons.nextBtn.Enable = iif(app.StepCompleted(6), 'on', 'off');
            app.NavButtons.nextBtn.Callback = @(~,~) app.updateStep(7);
        end

        function generateDLInput(app, toolPopup)
            if isempty(app.Model)
                msgbox('Please convert a model first.', 'Error', 'warn');
                return;
            end

            inputStatus = findobj(app.ContentPanels{6}, 'Style', 'text', 'Tag', 'inputStatus');
            if numel(inputStatus) > 1, inputStatus = inputStatus(1); end

            tools = {'CatPred', 'DLKcat', 'UniKP'};
            tool = tools{toolPopup.Value};

            try
                % Get parameters
                if ~isempty(app.Parameters) && isfield(app.Parameters, 'dataDir') && ~isempty(app.Parameters.dataDir)
                    params = app.Parameters;
                else
                    params = ParameterManager.getParams();
                end

                % Check if input file already exists
                inputFile = fullfile(params.dataDir, 'kcatData', [tool '_input.csv']);
                if exist(inputFile, 'file')
                    inputStatus.String = sprintf('Exists: %s_input.csv', tool);
                    inputStatus.ForegroundColor = [0.2 0.5 0.3];
                    app.setStatus([tool ' input file already exists']);
                    msgbox([tool ' input file already exists. No need to regenerate.'], 'Info', 'info');
                    return;
                end

                inputStatus.String = 'Generating...';
                inputStatus.ForegroundColor = [0.8 0.5 0.2];
                drawnow;

                % Generate input file
                Table = writeInputFile(app.Model, tool, params);

                inputStatus.String = sprintf('Done! Rows: %d', height(Table));
                inputStatus.ForegroundColor = [0.2 0.5 0.3];

                app.setStatus([tool ' input file generated']);

            catch ME
                inputStatus.String = 'Failed';
                inputStatus.ForegroundColor = [0.8 0.2 0.2];
                msgbox(ME.message, 'Generate Input Error', 'error');
            end
        end

        function runDLPrediction(app, toolPopup)
            result = DockerChecker.checkDocker();
            if ~result.running
                msgbox('Docker is not running. Please start Docker Desktop manually.', 'Warning', 'warn');
                return;
            end

            resultsText = findobj(app.ContentPanels{6}, 'Tag', 'step6Results');
            if numel(resultsText) > 1, resultsText = resultsText(1); end

            tools = {'CatPred', 'DLKcat', 'UniKP'};
            tool = tools{toolPopup.Value};

            try
                if ~isfield(app.Parameters, 'dataDir') || isempty(app.Parameters.dataDir)
                    params = ParameterManager.getParams();
                else
                    params = app.Parameters;
                end

                % Check if selected model's prediction file already exists
                kcatDir = fullfile(params.dataDir, 'kcatData');
                predFile = fullfile(kcatDir, [tool '_kcat.csv']);

                if exist(predFile, 'file')
                    resultsText.String = {[tool ' prediction already exists'], '', 'No need to run again.'};
                    app.setStatus([tool ' prediction file already exists']);
                    return;
                end

                resultsText.String = {['Running ' tool '...'], '', 'Please wait...'};
                drawnow;

                % Run the selected tool
                switch tool
                    case 'CatPred'
                        CatPred(kcatDir, params);
                    case 'DLKcat'
                        DLKcat(kcatDir, params);
                    case 'UniKP'
                        UniKP(kcatDir, params);
                end

                resultsText.String = {[tool ' prediction completed!'], '', 'Go to Step 7 to integrate kcat values.'};
                app.setStatus([tool ' kcat prediction completed']);

                % Mark step as completed
                app.StepCompleted(6) = true;
                app.updateStepListDisplay();
                app.NavButtons.nextBtn.Enable = 'on';

            catch ME
                resultsText.String = {[tool ' prediction failed:'], ME.message};
                msgbox(ME.message, 'Prediction Error', 'error');
            end
        end

        function selectKcatCallback(app, hObj, event)
            [file, path] = uigetfile('*.csv;*.xlsx;*.mat', 'Select kcat Prediction Results');
            if file == 0, return; end

            selectStatus = findobj(app.ContentPanels{6}, 'Style', 'text', 'Tag', 'selectStatus');
            if numel(selectStatus) > 1, selectStatus = selectStatus(1); end

            try
                fullPath = fullfile(path, file);
                if endsWith(file, '.csv')
                    kcatData = readtable(fullPath);
                elseif endsWith(file, '.xlsx')
                    kcatData = readtable(fullPath);
                else
                    kcatData = load(fullPath);
                end

                app.KcatList = kcatData;
                selectStatus.String = sprintf('Loaded: %s', file);
                selectStatus.ForegroundColor = [0.2 0.5 0.3];

                msgbox('kcat data loaded successfully!', 'Success', 'help');

            catch ME
                selectStatus.String = 'Failed';
                selectStatus.ForegroundColor = [0.8 0.2 0.2];
                msgbox(ME.message, 'Load Error', 'error');
            end
        end

        function checkDockerStatus(app, statusLabel)
            result = DockerChecker.checkDocker();
            statusLabel.String = ['Docker: ' DockerChecker.getStatusMessage(result)];
            if result.installed && result.running
                statusLabel.ForegroundColor = [0.2 0.5 0.3];
            else
                statusLabel.ForegroundColor = [0.8 0.2 0.2];
            end
        end

        % -----------------------------------------------------------------
        % STEP 7: Integrate kcat
        % -----------------------------------------------------------------
        function buildStep7_IntegrateKcat(app)
            p = app.ContentPanels{7};

            % Use consistent header
            app.buildStepHeader(p, 7, 'Integrate kcat Values', 'Select and integrate kcat values');

            contentY = 620;

            % Kcat Data selection card - left aligned
            kcatPanel = uipanel('Parent', p);
            kcatPanel.Title = 'Kcat Data';
            kcatPanel.FontSize = 10;
            kcatPanel.Units = 'pixels';
            kcatPanel.Position = [30 contentY-180 700 180];
            kcatPanel.BackgroundColor = [0.92 0.94 0.98];

            t = uicontrol('Parent', kcatPanel, 'Style', 'text');
            t.String = 'Kcat Data:';
            t.FontSize = 10;
            t.Units = 'pixels';
            t.Position = [15 135 80 25];

            kcatEdit = uicontrol('Parent', kcatPanel, 'Style', 'edit');
            kcatEdit.Units = 'pixels';
            kcatEdit.Position = [100 135 450 28];
            kcatEdit.BackgroundColor = [1 1 1];
            kcatEdit.Tag = 'kcatFile';

            btnBrowse = uicontrol('Parent', kcatPanel, 'Style', 'pushbutton');
            btnBrowse.String = 'Browse';
            btnBrowse.Units = 'pixels';
            btnBrowse.Position = [560 135 80 28];
            btnBrowse.BackgroundColor = [0.3 0.3 0.4];
            btnBrowse.ForegroundColor = [1 1 1];
            btnBrowse.Callback = @(~,~) app.browseKcatFile(kcatEdit);

            btnLoad = uicontrol('Parent', kcatPanel, 'Style', 'pushbutton');
            btnLoad.String = 'Load';
            btnLoad.Units = 'pixels';
            btnLoad.Position = [650 135 35 28];
            btnLoad.BackgroundColor = [0.3 0.5 0.7];
            btnLoad.ForegroundColor = [1 1 1];
            btnLoad.Callback = {@app.loadKcatCallback, kcatEdit};

            % Selection criteria
            t = uicontrol('Parent', kcatPanel, 'Style', 'text');
            t.String = 'Selection:';
            t.FontSize = 10;
            t.Units = 'pixels';
            t.Position = [15 85 70 25];

            criteriaPopup = uicontrol('Parent', kcatPanel, 'Style', 'popupmenu');
            criteriaPopup.String = 'max (recommended)|min|median|mean';
            criteriaPopup.Value = 1;
            criteriaPopup.Units = 'pixels';
            criteriaPopup.Position = [90 85 120 26];
            criteriaPopup.BackgroundColor = [1 1 1];
            criteriaPopup.Tag = 'criteriaPopup';

            t = uicontrol('Parent', kcatPanel, 'Style', 'text');
            t.String = 'Overwrite:';
            t.FontSize = 10;
            t.Units = 'pixels';
            t.Position = [230 85 80 25];

            overwritePopup = uicontrol('Parent', kcatPanel, 'Style', 'popupmenu');
            overwritePopup.String = 'true (always)|false (if zero)|ifHigher';
            overwritePopup.Value = 1;
            overwritePopup.Units = 'pixels';
            overwritePopup.Position = [315 85 150 26];
            overwritePopup.BackgroundColor = [1 1 1];
            overwritePopup.Tag = 'overwritePopup';

            % Run button
            btn = uicontrol('Parent', kcatPanel, 'Style', 'pushbutton');
            btn.String = 'Integrate kcat';
            btn.Units = 'pixels';
            btn.Position = [15 25 150 40];
            btn.BackgroundColor = [0.2 0.5 0.4];
            btn.ForegroundColor = [1 1 1];
            btn.FontSize = 11;
            btn.Callback = {@app.integrateKcatCallback, criteriaPopup, overwritePopup};

            % Results panel - left aligned
            resultsPanel = uipanel('Parent', p);
            resultsPanel.Title = 'Integration Results';
            resultsPanel.FontSize = 10;
            resultsPanel.Units = 'pixels';
            resultsPanel.Position = [30 30 700 290];
            resultsPanel.BackgroundColor = [0.92 0.94 0.98];

            resultsText = uicontrol('Parent', resultsPanel, 'Style', 'text');
            resultsText.Position = [10 10 680 265];
            resultsText.FontSize = 9;
            resultsText.FontName = 'Courier New';
            resultsText.HorizontalAlignment = 'left';
            resultsText.Units = 'pixels';
            resultsText.String = {'Load kcat data and run integration'};
            resultsText.Tag = 'integrateResults';

            % Navigation buttons - consistent bottom position
            app.NavButtons.prevBtn = uicontrol('Parent', p, 'Style', 'pushbutton');
            app.NavButtons.prevBtn.String = 'Previous';
            app.NavButtons.prevBtn.Units = 'pixels';
            app.NavButtons.prevBtn.Position = [800 30 100 35];
            app.NavButtons.prevBtn.BackgroundColor = [0.5 0.5 0.55];
            app.NavButtons.prevBtn.Callback = @(~,~) app.updateStep(6);

            app.NavButtons.nextBtn = uicontrol('Parent', p, 'Style', 'pushbutton');
            app.NavButtons.nextBtn.String = 'Next';
            app.NavButtons.nextBtn.Units = 'pixels';
            app.NavButtons.nextBtn.Position = [910 30 100 35];
            app.NavButtons.nextBtn.BackgroundColor = [0.2 0.4 0.6];
            app.NavButtons.nextBtn.ForegroundColor = [1 1 1];
            app.NavButtons.nextBtn.Enable = iif(app.StepCompleted(7), 'on', 'off');
            app.NavButtons.nextBtn.Callback = @(~,~) app.updateStep(8);
        end

        function browseKcatFile(app, kcatEdit)
            [file, path] = uigetfile({'*.csv', 'CSV Files'}, 'Select Kcat Data');
            if file ~= 0
                kcatEdit.String = fullfile(path, file);
            end
        end

        function loadKcatCallback(app, hObj, event, kcatEdit)
            filename = kcatEdit.String;
            if isempty(filename)
                msgbox('Please select a kcat file.', 'Error', 'warn');
                return;
            end

            try
                data = readtable(filename);
                app.KcatList = table2struct(data);
                app.setStatus(sprintf('Loaded %d kcat values', height(data)));
                msgbox(sprintf('Loaded %d kcat values', height(data)), 'Success', 'help');
            catch ME
                msgbox(ME.message, 'Error', 'error');
            end
        end

        function integrateKcatCallback(app, hObj, event, criteriaPopup, overwritePopup)
            if isempty(app.Model)
                msgbox('Please load and convert a model first.', 'Error', 'warn');
                return;
            end

            if isempty(app.KcatList)
                msgbox('Please load kcat data first.', 'Error', 'warn');
                return;
            end

            resultsText = findobj(app.ContentPanels{7}, 'Tag', 'integrateResults');

            criteria = {'max', 'min', 'median', 'mean'};
            overwrite = {'true', 'false', 'ifhigher'};

            try
                nBefore = sum(app.Model.enzymeConstraints.kcat > 0);

                model = selectKcatValue(app.Model, app.KcatList, ...
                    criteria{criteriaPopup.Value}, overwrite{overwritePopup.Value});

                app.Model = model;

                nAfter = sum(model.enzymeConstraints.kcat > 0);

                txt = {sprintf('Integration completed!'), '', ...
                       sprintf('Reactions with kcat: %d -> %d', nBefore, nAfter), ''};

                % Source breakdown
                if isfield(model.enzymeConstraints, 'source')
                    [C,~,ic] = unique(model.enzymeConstraints.source);
                    counts = accumarray(ic, 1);
                    txt{end+1} = 'Source breakdown:';
                    for i = 1:numel(C)
                        txt{end+1} = sprintf('  %s: %d', C{i}, counts(i));
                    end
                end

                if ~isempty(resultsText), resultsText.String = txt; end
                app.setStatus('kcat integration completed');

                % Mark step as completed and enable Next button
                app.StepCompleted(7) = true;
                app.updateStepListDisplay();
                app.NavButtons.nextBtn.Enable = 'on';
                msgbox('kcat integration completed! Click "Next" to continue.', 'Success', 'help');

            catch ME
                app.setStatus('Integration failed');
                msgbox(ME.message, 'Error', 'error');
            end
        end

        % -----------------------------------------------------------------
        % STEP 8: Growth Rate Solution
        % -----------------------------------------------------------------
        function buildStep8_GrowthRate(app)
            p = app.ContentPanels{8};

            % Use consistent header
            app.buildStepHeader(p, 8, 'Growth Rate Solution', 'Solve for growth rate using integrated ecModel');

            contentY = 620;

            % Solver options card - left aligned
            optionsPanel = uipanel('Parent', p);
            optionsPanel.Title = 'Solver Options';
            optionsPanel.FontSize = 10;
            optionsPanel.Units = 'pixels';
            optionsPanel.Position = [30 contentY-150 500 150];
            optionsPanel.BackgroundColor = [0.92 0.94 0.98];

            t = uicontrol('Parent', optionsPanel, 'Style', 'text');
            t.String = 'Objective:';
            t.FontSize = 10;
            t.Units = 'pixels';
            t.Position = [15 105 70 25];

            objPopup = uicontrol('Parent', optionsPanel, 'Style', 'popupmenu');
            objPopup.String = 'Biomass (max)|Minimize glucose uptake';
            objPopup.Value = 1;
            objPopup.Units = 'pixels';
            objPopup.Position = [90 103 180 26];
            objPopup.BackgroundColor = [1 1 1];
            objPopup.Tag = 'objPopup';

            t = uicontrol('Parent', optionsPanel, 'Style', 'text');
            t.String = 'Solver:';
            t.FontSize = 10;
            t.Units = 'pixels';
            t.Position = [15 60 70 25];

            solverPopup = uicontrol('Parent', optionsPanel, 'Style', 'popupmenu');
            solverPopup.String = 'gurobi|cplex|glpk';
            solverPopup.Value = 1;
            solverPopup.Units = 'pixels';
            solverPopup.Position = [90 58 100 26];
            solverPopup.BackgroundColor = [1 1 1];
            solverPopup.Tag = 'solverPopup';

            % Run button
            btn = uicontrol('Parent', optionsPanel, 'Style', 'pushbutton');
            btn.String = 'Solve Growth Rate';
            btn.Units = 'pixels';
            btn.Position = [15 10 150 40];
            btn.BackgroundColor = [0.2 0.5 0.4];
            btn.ForegroundColor = [1 1 1];
            btn.FontSize = 11;
            btn.Callback = @app.solveGrowthRateCallback;

            % Results panel - left aligned
            resultsPanel = uipanel('Parent', p);
            resultsPanel.Title = 'Growth Rate Results';
            resultsPanel.FontSize = 10;
            resultsPanel.Units = 'pixels';
            resultsPanel.Position = [30 30 700 290];
            resultsPanel.BackgroundColor = [0.92 0.94 0.98];

            resultsText = uicontrol('Parent', resultsPanel, 'Style', 'text');
            resultsText.Position = [10 10 680 265];
            resultsText.FontSize = 10;
            resultsText.HorizontalAlignment = 'left';
            resultsText.Units = 'pixels';
            resultsText.String = {'Run solver to get growth rate prediction'};
            resultsText.Tag = 'growthResults';

            % Navigation buttons - consistent bottom position
            app.NavButtons.prevBtn = uicontrol('Parent', p, 'Style', 'pushbutton');
            app.NavButtons.prevBtn.String = 'Previous';
            app.NavButtons.prevBtn.Units = 'pixels';
            app.NavButtons.prevBtn.Position = [800 30 100 35];
            app.NavButtons.prevBtn.BackgroundColor = [0.5 0.5 0.55];
            app.NavButtons.prevBtn.Callback = @(~,~) app.updateStep(7);

            app.NavButtons.nextBtn = uicontrol('Parent', p, 'Style', 'pushbutton');
            app.NavButtons.nextBtn.String = 'Next';
            app.NavButtons.nextBtn.Units = 'pixels';
            app.NavButtons.nextBtn.Position = [910 30 100 35];
            app.NavButtons.nextBtn.BackgroundColor = [0.2 0.4 0.6];
            app.NavButtons.nextBtn.ForegroundColor = [1 1 1];
            app.NavButtons.nextBtn.Enable = iif(app.StepCompleted(8), 'on', 'off');
            app.NavButtons.nextBtn.Callback = @(~,~) app.updateStep(9);
        end

        function solveGrowthRateCallback(app, hObj, event)
            if isempty(app.Model)
                msgbox('Please complete previous steps first.', 'Error', 'warn');
                return;
            end

            resultsText = findobj(app.ContentPanels{8}, 'Tag', 'growthResults');

            try
                if ~isempty(resultsText)
                    resultsText.String = {'Solving for growth rate...', '', 'This requires COBRA Toolbox and a solver (Gurobi/CPLEX).'};
                end
                drawnow;

                % Placeholder - actual implementation would call FBA/simulation
                % Using model.enzymeConstraints for constraint

                txt = {'Growth Rate Solution (Preliminary)', '', ...
                       'Note: Full implementation requires COBRA Toolbox', '', ...
                       'Model constraints applied:', ...
                       sprintf('  - Total protein: %.3f g/gDCW', app.Parameters.Ptot), ...
                       sprintf('  - Enzyme fraction: %.3f', app.Parameters.f), ...
                       sprintf('  - Sigma: %.3f', app.Parameters.sigma), '', ...
                       'The model is ready for FBA simulation.', '', ...
                       'Use verifyModel() to validate before simulation.'};

                if ~isempty(resultsText), resultsText.String = txt; end
                app.setStatus('Growth rate solution ready');

                % Mark step as completed and enable Next button
                app.StepCompleted(8) = true;
                app.updateStepListDisplay();
                app.NavButtons.nextBtn.Enable = 'on';
                msgbox('Growth rate solution completed (preliminary). Click "Next" to continue.', 'Success', 'help');

            catch ME
                if ~isempty(resultsText)
                    resultsText.String = {'Solution failed', '', ME.message};
                end
                msgbox(ME.message, 'Error', 'error');
            end
        end

        % -----------------------------------------------------------------
        % STEP 9: Save Model
        % -----------------------------------------------------------------
        function buildStep9_Save(app)
            p = app.ContentPanels{9};

            % Use consistent header
            app.buildStepHeader(p, 9, 'Save Model', 'Export the reconstructed ecModel');

            contentY = 620;

            % Export settings card - left aligned
            exportPanel = uipanel('Parent', p);
            exportPanel.Title = 'Export Settings';
            exportPanel.FontSize = 10;
            exportPanel.Units = 'pixels';
            exportPanel.Position = [30 contentY-250 500 250];
            exportPanel.BackgroundColor = [0.92 0.94 0.98];

            % Format selection
            t = uicontrol('Parent', exportPanel, 'Style', 'text');
            t.String = 'Format:';
            t.FontSize = 10;
            t.Units = 'pixels';
            t.Position = [15 200 60 25];

            fmtGrp = uibuttongroup('Parent', exportPanel);
            fmtGrp.Units = 'pixels';
            fmtGrp.Position = [80 200 380 30];
            fmtGrp.BackgroundColor = [0.92 0.94 0.98];

            r1 = uicontrol('Parent', fmtGrp, 'Style', 'radiobutton');
            r1.String = 'MATLAB (.mat)';
            r1.Units = 'pixels';
            r1.Position = [0 5 120 22];
            r1.Tag = 'mat';
            r1.Value = 1;

            r2 = uicontrol('Parent', fmtGrp, 'Style', 'radiobutton');
            r2.String = 'SBML (.xml)';
            r2.Units = 'pixels';
            r2.Position = [130 5 100 22];
            r2.Tag = 'xml';

            r3 = uicontrol('Parent', fmtGrp, 'Style', 'radiobutton');
            r3.String = 'JSON (.json)';
            r3.Units = 'pixels';
            r3.Position = [240 5 100 22];
            r3.Tag = 'json';

            % Output folder
            t = uicontrol('Parent', exportPanel, 'Style', 'text');
            t.String = 'Output:';
            t.FontSize = 10;
            t.Units = 'pixels';
            t.Position = [15 160 70 25];

            outputEdit = uicontrol('Parent', exportPanel, 'Style', 'edit');
            outputEdit.Units = 'pixels';
            outputEdit.Position = [90 160 250 28];
            outputEdit.BackgroundColor = [1 1 1];
            if isfield(app.Parameters, 'modelDir')
                outputEdit.String = app.Parameters.modelDir;
            else
                outputEdit.String = pwd;
            end
            outputEdit.Tag = 'outputFolder';

            btn = uicontrol('Parent', exportPanel, 'Style', 'pushbutton');
            btn.String = 'Browse';
            btn.Units = 'pixels';
            btn.Position = [350 160 70 28];
            btn.BackgroundColor = [0.3 0.3 0.4];
            btn.ForegroundColor = [1 1 1];
            btn.Callback = @(~,~) app.browseOutputPath(outputEdit);

            % Filename
            t = uicontrol('Parent', exportPanel, 'Style', 'text');
            t.String = 'Filename:';
            t.FontSize = 10;
            t.Units = 'pixels';
            t.Position = [15 120 70 25];

            modelId = 'ecModel';
            if ~isempty(app.Model) && isfield(app.Model, 'id')
                modelId = app.Model.id;
            end
            fnameEdit = uicontrol('Parent', exportPanel, 'Style', 'edit');
            fnameEdit.Units = 'pixels';
            fnameEdit.Position = [90 120 200 28];
            fnameEdit.BackgroundColor = [1 1 1];
            fnameEdit.String = [modelId, '_ec.mat'];
            fnameEdit.Tag = 'saveFilename';

            % Save button
            btn = uicontrol('Parent', exportPanel, 'Style', 'pushbutton');
            btn.String = 'Save Model';
            btn.Units = 'pixels';
            btn.Position = [15 45 150 45];
            btn.BackgroundColor = [0.2 0.5 0.4];
            btn.ForegroundColor = [1 1 1];
            btn.FontSize = 11;
            btn.Callback = {@app.saveModelCallback, fnameEdit, outputEdit, fmtGrp};

            % Summary panel - left aligned
            summaryPanel = uipanel('Parent', p);
            summaryPanel.Title = 'Model Summary';
            summaryPanel.FontSize = 10;
            summaryPanel.Units = 'pixels';
            summaryPanel.Position = [30 30 700 270];
            summaryPanel.BackgroundColor = [0.92 0.94 0.98];

            summaryText = uicontrol('Parent', summaryPanel, 'Style', 'text');
            summaryText.Position = [10 10 680 245];
            summaryText.FontSize = 9;
            summaryText.FontName = 'Courier New';
            summaryText.HorizontalAlignment = 'left';
            summaryText.Units = 'pixels';

            if ~isempty(app.Model)
                m = app.Model;
                txt = {sprintf('Model ID: %s', m.id), ...
                       sprintf('Model Type: %s', m.type), '', ...
                       sprintf('Reactions: %d', numel(m.rxns)), ...
                       sprintf('Metabolites: %d', numel(m.mets)), ...
                       sprintf('Genes: %d', numel(m.genes)), ''};
                if isfield(m, 'enzymeConstraints')
                    ec = m.enzymeConstraints;
                    txt{end+1} = sprintf('EC Reactions: %d', numel(ec.rxns));
                    if isfield(ec, 'kcat')
                        txt{end+1} = sprintf('Reactions with kcat: %d', sum(ec.kcat > 0));
                    end
                end
                txt{end+1} = '';
                txt{end+1} = sprintf('Project: %s', app.Parameters.projectName);
                summaryText.String = txt;
            else
                summaryText.String = {'No model to save'};
            end

            % Navigation - Step 9 is final, only Previous button
            app.NavButtons.prevBtn = uicontrol('Parent', p, 'Style', 'pushbutton');
            app.NavButtons.prevBtn.String = 'Previous';
            app.NavButtons.prevBtn.Units = 'pixels';
            app.NavButtons.prevBtn.Position = [800 30 100 35];
            app.NavButtons.prevBtn.BackgroundColor = [0.5 0.5 0.55];
            app.NavButtons.prevBtn.Callback = @(~,~) app.updateStep(8);
        end

        function browseOutputPath(app, outputEdit)
            folder = uigetdir(pwd, 'Select Output Directory');
            if folder ~= 0
                outputEdit.String = folder;
            end
        end

        function saveModelCallback(app, hObj, event, fnameEdit, outputEdit, fmtGrp)
            filename = fnameEdit.String;
            outputFolder = outputEdit.String;
            tag = get(fmtGrp.SelectedObject, 'Tag');

            if isempty(app.Model)
                msgbox('No model to save.', 'Error', 'warn');
                return;
            end

            [~, name, ext] = fileparts(filename);
            if isempty(ext)
                switch tag
                    case 'mat', ext = '.mat';
                    case 'xml', ext = '.xml';
                    case 'json', ext = '.json';
                end
            end

            outputPath = fullfile(outputFolder, [name, ext]);

            try
                switch tag
                    case 'mat'
                        modelToSave = app.Model;
                        save(outputPath, 'modelToSave');
                    otherwise
                        msgbox('Export format not yet implemented.', 'Info', 'help');
                        return;
                end
                app.setStatus(['Saved: ', outputPath]);

                % Mark final step as completed
                app.StepCompleted(9) = true;
                app.updateStepListDisplay();
                msgbox({['Model saved to:'], outputPath, '', 'Reconstruction workflow complete!'}, 'Success', 'help');
            catch ME
                msgbox(ME.message, 'Save Error', 'error');
            end
        end

        % -----------------------------------------------------------------
        % Event Handlers
        % -----------------------------------------------------------------
        function onStepChanged(app, hObj, event)
            stepNum = get(hObj, 'Value');
            if stepNum ~= app.CurrentStep
                % DEBUG MODE: Allow free navigation to any step
                app.updateStep(stepNum);
            end
        end

        function onClose(app, hObj, event)
            selection = questdlg('Close?', 'Confirm', 'Yes', 'No', 'No');
            if strcmp(selection, 'Yes')
                delete(app.UIFigure);
            end
        end
    end
end

% Helper function
function out = iif(cond, trueVal, falseVal)
    if cond
        out = trueVal;
    else
        out = falseVal;
    end
end
