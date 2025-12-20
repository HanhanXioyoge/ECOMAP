function model = updateProtPool(model, updateFfactor, Ptot, f, sigma, parameters)
% updateProtPool
%   This function updates the protein pool in an ECOMAP ecModel based on the 
%   fraction of protein mass allocated to enzymes (f-factor). It adjusts the 
%   model's protein pool exchange reaction using proteomics data, model parameters, 
%   and the current protein fraction.
%
% Inputs:
%   model           - ECOMAP ecModel (including enzymeConstraints structure)
%   updateFfactor   - Boolean flag to specify whether to update the f-factor 
%                     based on proteomics data (optional, default: false)
%   Ptot            - Total protein content in the cell (optional, default: model.enzymeConstraints.Ptot)
%   f               - Initial f-factor (fraction of protein mass for enzymes) (optional, default: model.enzymeConstraints.f)
%   sigma           - Scaling factor for the protein pool (optional, default: model.enzymeConstraints.sigma)
%   parameters      - Struct containing additional parameters (must include 'dataDir') (optional)
%
% Output:
%   model           - Updated ecModel with modified protein pool exchange (prot_pool_exchange)
%
% Usage: 
%   model = updateProtPool(model, updateFfactor, Ptot, f, sigma, parameters);

    % ------------------------- Parameter Defaults -------------------------
    % Ensure parameters are set if not provided
    if nargin < 6 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters)
            error('ParameterManager is not set.');
        end
    end
    
    % Check if the data directory is provided
    if ~isfield(parameters, 'dataDir') || isempty(parameters.dataDir)
        error('parameters.dataDir is required and must point to the folder containing CSV files.');
    end
    dataDir = parameters.dataDir;
    
    % Set default values for optional inputs
    if nargin < 2 || isempty(updateFfactor)
        updateFfactor = false;  % Default: do not update f-factor
    end

    if nargin < 5 || isempty(sigma)
        sigma = model.enzymeConstraints.sigma;  % Default scaling factor
    end
    if nargin < 4 || isempty(f)
        f = model.enzymeConstraints.f;  % Default enzyme protein fraction
    end
    if nargin < 3 || isempty(Ptot)
        Ptot = model.enzymeConstraints.Ptot;  % Default total protein content
    end

    % -------------------- Update f-factor Using Proteomics Data --------------------
    % Update the f-factor if specified by the user
    if updateFfactor
        if exist(fullfile(dataDir, 'paxDB.tsv'), 'file')
            protData = fullfile(dataDir, 'paxDB.tsv');
            
            % Load UniProt data if necessary
            uniprot_Path = fullfile(parameters.dataDir, 'uniprot.tsv');
            if ~isfile(uniprot_Path)
                DownloadUniProtData(parameters.uniprot, parameters.dataDir);  % Download UniProt data if not available
            end
            dbStruct = ParseUniProtData(uniprot_Path);  % Parse UniProt data
            
            % Parse the PaxDB file if available
            if ischar(protData) && endsWith(protData, 'paxDB.tsv')
                % Read the PaxDB file and skip header lines
                fID = fopen(fullfile(protData), 'r');
                fileContent = textscan(fID, '%s', 'delimiter', '\n');
                headerLines = find(startsWith(fileContent{1}, '#'), 1, 'last');
                fclose(fID);
                
                % Read actual data: gene names and abundance levels
                fID = fopen(fullfile(protData), 'r');
                fileContent = textscan(fID, '%s %s %f', 'delimiter', '\t', 'HeaderLines', headerLines);
                genes = fileContent{2};
                genes = regexprep(genes, '^\d+\.','');  % Clean up gene IDs
                level = fileContent{3};
                fclose(fID);
                
                % Match genes to UniProt IDs
                a = cellfun(@(gene) any(contains(dbStruct.genes, gene)), genes);
                b = cellfun(@(gene) find(contains(dbStruct.genes, gene), 1, 'first'), genes, 'UniformOutput', false);
                b_values = cell2mat(b(a));
                uniprot = dbStruct.ID(b_values);
                level(~a) = [];  % Remove levels for unmatched genes
                clear protData;
                
                % Store matched UniProt data in protData structure
                protData.uniprotIDs = uniprot;
                protData.level = level;
                
                % Retrieve molecular weight (MW) and abundance data
                [~, idx] = ismember(protData.uniprotIDs, dbStruct.ID);
                protData.MW = dbStruct.MW(idx);
                protData.abundances = protData.level .* protData.MW;
            end
            
            % Calculate average abundances and total protein content
            avgAbundances = mean(protData.abundances, 2);
            totalProt = sum(avgAbundances, 'omitnan');
            
            % Calculate total enzyme mass from model
            enzymes = model.enzymeConstraints.enzymes;
            enzymesInModel = ismember(protData.uniprotIDs, enzymes);
            totalEnz = sum(avgAbundances(enzymesInModel), 'omitnan');
            
            % Update the f-factor based on proteomics data
            f_old = f;
            f = totalEnz / totalProt;
            model.enzymeConstraints.f = f;  % Update f-factor in the model
            fprintf('The enzyme protein fraction has been updated from %.2f to %.2f\n', f_old, f);

        else
            % If PaxDB data is not found, return the original f-value and issue a warning
            warning('No proteomics data is available. The original f-value of %.2f is retained.\n', f);
        end
    end

    % -------------------- Update Protein Pool Exchange Reaction --------------------
    % Update the lower bound (lb) for the prot_pool_exchange reaction based on f-factor
    prot_pool_old = model.lb(strcmp(model.rxns, 'prot_pool_exchange'));
    prot_pool_new = -(Ptot * f * sigma * 1000);  % Adjust the lower bound for the reaction
    model.lb(strcmp(model.rxns, 'prot_pool_exchange')) = prot_pool_new;  % Update protein pool exchange rate
    fprintf('The protein pool constraint lower bound has been updated from %.2f to %.2f\n', prot_pool_old, prot_pool_new);
end
