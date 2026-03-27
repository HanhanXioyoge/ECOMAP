classdef ModelStats
    % ModelStats
    % Utility class for visualizing model statistics in reconstruction workflow
    % Academic/minimalist style

    methods (Static)

        function createStatsPanel(model, parent)
            % createStatsPanel
            %   Creates a panel showing model statistics
            %
            % Inputs:
            %   model  - COBRA/ecModel structure
            %   parent - uipanel handle to contain the stats

            if ~isfield(model, 'rxns'), return; end

            % Calculate stats
            nRxns = numel(model.rxns);
            nMets = numel(model.mets);
            nGenes = numel(model.genes);

            % Reaction stats
            if isfield(model, 'rev')
                nReversible = sum(model.rev);
                nIrreversible = nRxns - nReversible;
            else
                nReversible = 0;
                nIrreversible = nRxns;
            end

            % Gene-associated reactions
            if isfield(model, 'rxnGeneMat')
                nGeneRxns = sum(any(model.rxnGeneMat ~= 0, 2));
            else
                nGeneRxns = 0;
            end

            % Create grid layout
            grid = uix.Grid('Parent', parent, ...
                'Padding', [15 10 15 10], ...
                'Spacing', [10 8]);

            % Row 1: Title
            uicontrol('Parent', grid, ...
                'Style', 'text', ...
                'String', 'Model Overview', ...
                'FontSize', 12, ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left');

            % Empty spacer
            uicontrol('Parent', grid, 'Style', 'text', 'String', '');

            % Stats rows
            ModelStats.addStatRow(grid, 'Reactions', num2str(nRxns));
            ModelStats.addStatRow(grid, 'Metabolites', num2str(nMets));
            ModelStats.addStatRow(grid, 'Genes', num2str(nGenes));
            ModelStats.addStatRow(grid, 'Reversible', num2str(nReversible));
            ModelStats.addStatRow(grid, 'Gene-associated', num2str(nGeneRxns));

            % Enzyme constraints info
            if isfield(model, 'enzymeConstraints')
                ec = model.enzymeConstraints;
                if isfield(ec, 'rxns')
                    ModelStats.addStatRow(grid, 'EC Reactions', num2str(numel(ec.rxns)));
                end
                if isfield(ec, 'genes')
                    ModelStats.addStatRow(grid, 'EC Genes', num2str(numel(ec.genes)));
                end
                if isfield(ec, 'sigma')
                    ModelStats.addStatRow(grid, 'Sigma', sprintf('%.3f', ec.sigma));
                end
                if isfield(ec, 'Ptot')
                    ModelStats.addStatRow(grid, 'Ptot', sprintf('%.3f', ec.Ptot));
                end
                if isfield(ec, 'f')
                    ModelStats.addStatRow(grid, 'f', sprintf('%.3f', ec.f));
                end
            end

            % Set column widths
            grid.ColumnWidth = {120, 80};
            grid.RowHeight = {20, 15, 20, 20, 20, 20, 20};
        end

        function addStatRow(grid, label, value)
            % addStatRow
            %   Helper to add a label-value row to a grid

            uicontrol('Parent', grid, ...
                'Style', 'text', ...
                'String', label, ...
                'FontSize', 10, ...
                'HorizontalAlignment', 'left');

            uicontrol('Parent', grid, ...
                'Style', 'text', ...
                'String', value, ...
                'FontSize', 10, ...
                'FontWeight', 'normal', ...
                'HorizontalAlignment', 'right');
        end

        function hBar = createRxnTypeBar(model, parent)
            % createRxnTypeBar
            %   Creates a horizontal bar chart showing reaction types
            %
            % Output:
            %   hBar - axes handle

            hBar = axes('Parent', parent, ...
                'Units', 'normalized', ...
                'Position', [0.1 0.1 0.8 0.4]);

            categories = {};
            counts = {};

            % Count exchange reactions
            if isfield(model, 'rxns')
                exchIdx = contains(lower(model.rxns), 'exchange');
                if any(exchIdx)
                    categories{end+1} = 'Exchange';
                    counts{end+1} = sum(exchIdx);
                end

                % Count transport
                transportMask = false(numel(model.rxns), 1);
                keywords = {'transport', 'transporter', ' pump ', ' permease'};
                for i = 1:numel(keywords)
                    transportMask = transportMask | contains(lower(model.rxnNames), keywords{i});
                end
                if any(transportMask)
                    categories{end+1} = 'Transport';
                    counts{end+1} = sum(transportMask);
                end

                % Count metabolic
                metabolicMask = ~exchIdx & ~transportMask;
                if any(metabolicMask)
                    categories{end+1} = 'Metabolic';
                    counts{end+1} = sum(metabolicMask);
                end
            end

            if isempty(categories)
                return;
            end

            % Create horizontal bar
            yvals = 1:numel(categories);
            barh(yvals, cell2mat(counts), 0.6, 'FaceColor', [0.3 0.5 0.7]);

            set(hBar, 'YTick', yvals, ...
                'YTickLabel', categories, ...
                'FontSize', 9, ...
                'TickDir', 'out', ...
                'Box', 'off');

            xlabel(hBar, 'Count', 'FontSize', 9);

            % Add value labels
            for i = 1:numel(counts)
                text(cell2mat(counts(i)), yvals(i), ...
                    ['  ', num2str(cell2mat(counts(i)))], ...
                    'Parent', hBar, ...
                    'FontSize', 9, ...
                    'VerticalAlignment', 'middle');
            end
        end

        function hPie = createGeneAssociationPie(model, parent)
            % createGeneAssociationPie
            %   Creates a pie chart showing gene-associated vs orphan reactions
            %
            % Output:
            %   hPie - axes handle

            hPie = axes('Parent', parent, ...
                'Units', 'normalized', ...
                'Position', [0.1 0.1 0.8 0.4]);

            if ~isfield(model, 'rxnGeneMat')
                return;
            end

            geneAssociated = sum(any(model.rxnGeneMat ~= 0, 2));
            orphan = numel(model.rxns) - geneAssociated;

            if geneAssociated == 0 && orphan == 0
                return;
            end

            slices = [geneAssociated, orphan];
            labels = {sprintf('Gene-associated\n%d', geneAssociated), ...
                      sprintf('Orphan\n%d', orphan)};

            pie(hPie, slices, {'#3498db', '#95a5a6'});

            title(hPie, 'Reaction Gene Association', 'FontSize', 10, 'FontWeight', 'bold');
        end
    end
end
