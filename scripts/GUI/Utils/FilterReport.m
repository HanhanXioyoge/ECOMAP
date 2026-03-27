classdef FilterReport
    % FilterReport
    % Utility class for visualizing enzyme constraint reaction filter results
    % Academic/minimalist style

    methods (Static)

        function createReportPanel(report, parent)
            % createReportPanel
            %   Creates a panel displaying the filter report summary
            %
            % Inputs:
            %   report - struct from buildEnzConstrRxnSet()

            if ~isfield(report, 'summary')
                return;
            end

            % Create container panel
            mainPanel = uix.VBox('Parent', parent, ...
                'Padding', 10, ...
                'Spacing', 10);

            % Title
            uicontrol('Parent', mainPanel, ...
                'Style', 'text', ...
                'String', 'Reaction Filter Report', ...
                'FontSize', 11, ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left');

            % Summary stats
            summaryGrid = uix.Grid('Parent', mainPanel, ...
                'Padding', 5, ...
                'Spacing', [5 5]);

            s = report.summary;
            ModelStats.addStatRow(summaryGrid, 'Total Reactions', num2str(s.totalRxns));
            ModelStats.addStatRow(summaryGrid, 'Kept', num2str(s.totalKeep));
            ModelStats.addStatRow(summaryGrid, 'Filtered', num2str(s.totalDrop));

            summaryGrid.ColumnWidth = {120, 60};

            % Percentage bar
            if s.totalRxns > 0
                keptRatio = s.totalKeep / s.totalRxns;
                FilterReport.createRetentionBar(mainPanel, keptRatio);
            end

            % Rule breakdown table
            if isfield(report, 'rules')
                FilterReport.createRuleTable(mainPanel, report.rules);
            end

            % Set sizes
            mainPanel.Heights = [20, 80, 40, -1];
        end

        function createRetentionBar(parent, ratio)
            % createRetentionBar
            %   Creates a simple retention rate bar

            ax = axes('Parent', parent, ...
                'Units', 'normalized', ...
                'Position', [0.05 0.05 0.9 0.15], ...
                'Color', [0.95 0.95 0.95], ...
                'XTick', [], ...
                'YTick', [], ...
                'Box', 'on');

            hold(ax, 'on');

            % Background (total)
            bar(ax, 1, 1, 'FaceColor', [0.85 0.85 0.85], 'BarWidth', 0.6);

            % Kept portion
            if ratio > 0
                bar(ax, 1, ratio, 'FaceColor', [0.2 0.5 0.3], 'BarWidth', 0.6);
            end

            % Label
            text(0.5, 0.5, sprintf('Retention: %.1f%%', ratio * 100), ...
                'Parent', ax, ...
                'FontSize', 9, ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Color', 'white');

            hold(ax, 'off');
            axis(ax, 'off');
        end

        function createRuleTable(parent, rulesTable)
            % createRuleTable
            %   Creates a table showing filter rule statistics

            % Container
            tablePanel = uix.VBox('Parent', parent, ...
                'Padding', 5);

            uicontrol('Parent', tablePanel, ...
                'Style', 'text', ...
                'String', 'Filter Rules', ...
                'FontSize', 9, ...
                'FontWeight', 'bold');

            % Convert table for uitable
            if istable(rulesTable)
                colnames = rulesTable.Properties.VariableNames;
                data = table2cell(rulesTable);
            else
                colnames = {'Rule', 'Description', 'Enabled', 'Dropped'};
                data = {};
            end

            hTable = uitable('Parent', tablePanel, ...
                'Data', data, ...
                'ColumnName', colnames, ...
                'FontSize', 8, ...
                'RowName', [], ...
                'ColumnWidth', 'auto', ...
                'BackgroundColor', [1 1 1], ...
                'ForegroundColor', [0.2 0.2 0.2]);

            tablePanel.Heights = [15, -1];
        end

        function hBar = createRuleDropBar(report, parent)
            % createRuleDropBar
            %   Creates a horizontal bar chart showing drops by rule
            %
            % Output:
            %   hBar - axes handle

            hBar = axes('Parent', parent, ...
                'Units', 'normalized', ...
                'Position', [0.15 0.2 0.75 0.6]);

            if ~isfield(report, 'summary') || ~isfield(report.summary, 'byRule')
                return;
            end

            rules = report.summary.byRule;
            fields = fieldnames(rules);

            ruleLabels = FilterReport.getRuleShortNames(fields);
            ruleCounts = cell2mat(struct2cell(rules));

            % Sort by count descending
            [sortedCounts, sortIdx] = sort(ruleCounts, 'descend');
            sortedLabels = ruleLabels(sortIdx);

            % Remove zero entries
            nonZero = sortedCounts > 0;
            sortedCounts = sortedCounts(nonZero);
            sortedLabels = sortedLabels(nonZero);

            if isempty(sortedCounts)
                return;
            end

            yvals = 1:numel(sortedCounts);

            % Color gradient: red for high drop, orange for medium
            colors = [linspace(0.8, 0.4, numel(sortedCounts))', ...
                      linspace(0.2, 0.4, numel(sortedCounts))', ...
                      linspace(0.2, 0.2, numel(sortedCounts))'];

            barh(yvals, sortedCounts, 0.6, 'FaceColor', 'flat', 'CData', colors);

            set(hBar, 'YTick', yvals, ...
                'YTickLabel', sortedLabels, ...
                'FontSize', 8, ...
                'TickDir', 'out', ...
                'Box', 'off');

            xlabel(hBar, 'Reactions Dropped', 'FontSize', 8);

            % Value labels
            for i = 1:numel(sortedCounts)
                text(sortedCounts(i) + 0.5, yvals(i), ...
                    num2str(sortedCounts(i)), ...
                    'Parent', hBar, ...
                    'FontSize', 8, ...
                    'VerticalAlignment', 'middle');
            end

            title(hBar, 'Reactions Dropped by Filter Rule', 'FontSize', 9, 'FontWeight', 'bold');
        end

        function labels = getRuleShortNames(ruleKeys)
            % getRuleShortNames
            %   Maps rule keys to short display names

            shortNames = containers.Map(...
                {'R0_gene', 'R2_exch', 'R3_pure', 'R4_proton', 'R5_stoich', 'R6_name'}, ...
                {'No Gene', 'Exchange', 'Pure Trans.', 'Proton Coupling', 'Stoich. ID', 'Name Filter'});

            labels = cell(numel(ruleKeys), 1);
            for i = 1:numel(ruleKeys)
                if shortNames.isKey(ruleKeys{i})
                    labels{i} = shortNames(ruleKeys{i});
                else
                    labels{i} = ruleKeys{i};
                end
            end
        end
    end
end
