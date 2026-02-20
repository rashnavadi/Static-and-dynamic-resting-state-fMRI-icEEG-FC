%% Fig_Run1a_ICE017_plot_ordered_fMRI_heatmap.m
clear; clc;

subjID  = "ICE017";
runKey  = "Run1a";

rootDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/FC_rebuild";
subjDir = fullfile(rootDir, subjID);
fmriFile = fullfile(subjDir, "FMRI_FC_" + subjID + ".mat");
load(fmriFile, "OUT");

names = string(OUT.(runKey).seedNames);
% Z     = OUT.(runKey).FC_static_z;
Z = OUT.(runKey).FC_static;  % r


% --- build ordering ---
[ord, grp, grpLabels] = order_ICE017_seednames(names);

names_ord = names(ord);
Z_ord     = Z(ord, ord);

% --- plot + save ---
outDir = fullfile(rootDir, "RESULTS_REPRESENTATIVE", "ordered_heatmaps_fMRI", subjID);
if ~exist(outDir,"dir"), mkdir(outDir); end
outPng = fullfile(outDir, sprintf("%s_fMRI_static_z_ORDERED_%s_N%d.png", subjID, runKey, numel(names)));

plot_heatmap_z_noDiag_withNames_andGroupLines(Z_ord, names_ord, grp(ord), grpLabels, ...
    sprintf("%s fMRI static FC_z (ordered, no diag) | %s | N=%d", subjID, runKey, size(Z_ord,1)), ...
    outPng);

fprintf("[OK] Saved: %s\n", outPng);

%% ===================== helpers =====================

function [ord, grp, grpLabels] = order_ICE017_seednames(names)
    % For ICE017 Run1a names like: DRMOF1-DRMOF2
    % We sort by region group rank, then by contact index (left side).
    n = numel(names);

    region = strings(n,1);
    idxNum = zeros(n,1);

    for i = 1:n
        nm = upper(names(i));
        left = split(nm, "-");
        left = left(1);                     % e.g., DRMOF1
        tok = regexp(left, '^(DR)([A-Z]+)(\d+)$', 'tokens', 'once');
        if isempty(tok)
            region(i) = "UNK";
            idxNum(i) = 0;
        else
            region(i) = string(tok{2});     % e.g., MOF, LOF, FOP, AIN, PIN, PT
            idxNum(i) = str2double(tok{3}); % e.g., 1,2,3...
        end
    end

    % Region rank for this subject (most interpretable blocks)
    % OF: MOF then LOF
    % then OP (FOP)
    % then Insula (AIN then PIN)
    % then PT
    regRank = zeros(n,1);
    for i = 1:n
        switch region(i)
            case "MOF", regRank(i) = 1;
            case "LOF", regRank(i) = 2;
            case "FOP", regRank(i) = 3;
            case "AIN", regRank(i) = 4;
            case "PIN", regRank(i) = 5;
            case "PT",  regRank(i) = 6;
            otherwise,  regRank(i) = 99;
        end
    end

    % Group id (for drawing separator lines)
    grp = regRank;  % since regRank is already group-coded

    % Final ordering: group then index
    T = table(regRank, idxNum, (1:n)', region);
    T = sortrows(T, {'regRank','idxNum','region'});
    ord = T.Var3;

    % For legend labels in separators
    grpLabels = containers.Map('KeyType','double','ValueType','char');
    grpLabels(1) = 'MOF';
    grpLabels(2) = 'LOF';
    grpLabels(3) = 'FOP';
    grpLabels(4) = 'AIN';
    grpLabels(5) = 'PIN';
    grpLabels(6) = 'PT';
end

function plot_heatmap_z_noDiag_withNames_andGroupLines(Z, names, grp, grpLabels, ttl, outPng)
    names = string(names);
    n = size(Z,1);

    % remove diagonal (don’t let it dominate)
    Z(1:n+1:end) = NaN;

    % use robust clim based on off-diagonal only (keeps contrast)
    v = Z(~isnan(Z));
    clim = prctile(v, [2 98]);
    if clim(1) == clim(2)
        clim = [clim(1)-0.1, clim(2)+0.1];
    end

    fig = figure("Visible","off");
    h = imagesc(Z);
    axis image;

    cb = colorbar;
    cb.Location = 'eastoutside';
    cb.TickDirection = 'out';

    % Make space for colorbar labels
    ax = gca;
    ax.Position(3) = ax.Position(3) * 0.88;  % shrink width → more space on right

    title(ttl, "Interpreter","none");

    % show ALL labels (N=24 is fine)
    set(gca, "XTick", 1:n, "XTickLabel", names, ...
             "YTick", 1:n, "YTickLabel", names, ...
             "TickLabelInterpreter","none", ...
             "FontSize", 7);
    xtickangle(90);

    set(h, "AlphaData", ~isnan(Z));
    %     caxis(clim);
    caxis([-.5 .5]);


    % draw separator lines between groups
    hold on;
    boundaries = find(diff(grp) ~= 0);
    for b = boundaries'
        x = b + 0.5;
        plot([0.5, n+0.5], [x, x], 'k-', 'LineWidth', 1.2);
        plot([x, x], [0.5, n+0.5], 'k-', 'LineWidth', 1.2);
    end
    hold off;

    exportgraphics(fig, outPng, "Resolution", 300);
    close(fig);
end
