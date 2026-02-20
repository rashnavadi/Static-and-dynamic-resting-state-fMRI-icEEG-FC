outDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_spatial";
T = readtable(fullfile(outDir,"SpatialCoupling_longform_ALLbands.csv"));

bandOrder = ["delta","theta","alpha","beta","low_gamma","high_gamma"];

% make sure band is string/categorical comparable
T.band = string(T.band);

Summary = table('Size',[numel(bandOrder) 4], ...
    'VariableTypes', {'string','string','string','double'}, ...
    'VariableNames', {'Band','Significant_region_pairs','Beta_range','Mean_beta'});

for i = 1:numel(bandOrder)
    b = bandOrder(i);
    Tb = T(T.band == b, :);

    % significant pairs using the precomputed flag
    sigMask = Tb.sig_p005 == 1;
    sigVals = Tb.beta_dev(sigMask);

    nSig = sum(sigMask);
    nTot = height(Tb);   % should be 64 if 8x8 region pairs

    if isempty(sigVals)
        brange = "NA";
        mBeta = NaN;
    else
        brange = sprintf("%.3f \x2192 %.3f", min(sigVals), max(sigVals));
        mBeta  = mean(sigVals);
    end

    Summary.Band(i) = b;
    Summary.Significant_region_pairs(i) = sprintf("%d / %d", nSig, nTot);
    Summary.Beta_range(i) = brange;
    Summary.Mean_beta(i) = mBeta;
end

disp(Summary)
