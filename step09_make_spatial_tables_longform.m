%RUN EXAMPLE
% outDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_spatial";

function step09_make_spatial_tables_longform(outDir)

if nargin < 1 || isempty(outDir)
    outDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_spatial";
end

bandOrder = ["delta","theta","alpha","beta","low_gamma","high_gamma"];

Tall = table();

for b = bandOrder
    betaPath  = fullfile(outDir, "BETA_dev_" + b + ".csv");
    pPath     = fullfile(outDir, "P_dev_"    + b + ".csv");
    nsubPath  = fullfile(outDir, "Nsubj_"    + b + ".csv");
    npairPath = fullfile(outDir, "Npairs_"   + b + ".csv");
    nobsPath  = fullfile(outDir, "Nobs_"     + b + ".csv");

    if ~exist(betaPath,"file")
        fprintf("[SKIP] Missing %s\n", betaPath);
        continue;
    end

    % --- Read matrices ---
    Tbeta   = readtable(betaPath, "ReadRowNames", true);
    regions = string(Tbeta.Properties.VariableNames);
    rows    = string(Tbeta.Properties.RowNames);
    BETA    = table2array(Tbeta);

    P = nan(size(BETA));
    if exist(pPath,"file")
        Tp = readtable(pPath, "ReadRowNames", true);
        P  = table2array(Tp);
    end

    Nsubj = nan(size(BETA));
    if exist(nsubPath,"file")
        Tn = readtable(nsubPath, "ReadRowNames", true);
        Nsubj = table2array(Tn);
    end

    Npairs = nan(size(BETA));
    if exist(npairPath,"file")
        Tn = readtable(npairPath, "ReadRowNames", true);
        Npairs = table2array(Tn);
    end

    Nobs = nan(size(BETA));
    if exist(nobsPath,"file")
        Tn = readtable(nobsPath, "ReadRowNames", true);
        Nobs = table2array(Tn);
    end

    % --- Build long-form rows via cell array (safe) ---
    nR = numel(rows);
    nC = numel(regions);
    nTot = nR * nC;

    C = cell(nTot, 8); % band, regA, regB, beta, p, nSubj, nPairs, nObs
    k = 1;

    for i = 1:nR
        for j = 1:nC
            C{k,1} = string(b);
            C{k,2} = rows(i);
            C{k,3} = regions(j);
            C{k,4} = BETA(i,j);
            C{k,5} = P(i,j);
            C{k,6} = Nsubj(i,j);
            C{k,7} = Npairs(i,j);
            C{k,8} = Nobs(i,j);
            k = k + 1;
        end
    end

    Tlong = cell2table(C, "VariableNames", ...
        ["band","regA","regB","beta_dev","p_dev","nSubj","nPairs","nObs"]);

    % ensure numeric columns are numeric
    Tlong.beta_dev = double(Tlong.beta_dev);
    Tlong.p_dev    = double(Tlong.p_dev);
    Tlong.nSubj    = double(Tlong.nSubj);
    Tlong.nPairs   = double(Tlong.nPairs);
    Tlong.nObs     = double(Tlong.nObs);

    % keep only modeled cells (beta not NaN)
    Tlong = Tlong(~isnan(Tlong.beta_dev), :);

    Tall = [Tall; Tlong]; %#ok<AGROW>
end

if isempty(Tall)
    error("No long-form rows created. Check that BETA_dev_*.csv exist and contain numeric values.");
end

% Simple significance flag (raw p)
Tall.sig_p005 = Tall.p_dev < 0.05;

% Save
csvOut = fullfile(outDir, "SpatialCoupling_longform_ALLbands.csv");
matOut = fullfile(outDir, "SpatialCoupling_longform_ALLbands.mat");

writetable(Tall, csvOut);
save(matOut, "Tall", "-v7.3");

fprintf("[OK] Saved long-form spatial table to:\n%s\n", csvOut);

end

