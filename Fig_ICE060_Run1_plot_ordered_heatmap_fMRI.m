%% Fig_Run1a_ICE060_plot_ordered_fMRI_heatmap.m
clear; clc;

subjID  = "ICE057";
runKey  = "Run1";

rootDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/FC_rebuild";
subjDir = fullfile(rootDir, subjID);
fmriFile = fullfile(subjDir, "FMRI_FC_" + subjID + ".mat");
load(fmriFile, "OUT");

names = string(OUT.(runKey).seedNames);
% Z     = OUT.(runKey).FC_static_z;
Z = OUT.(runKey).FC_static;  % r


% --- build ordering ---
[ord, grp, grpLabels] = order_ICE060_seednames(names);

names_ord = names(ord);
Z_ord     = Z(ord, ord);

% --- plot + save ---
outDir = fullfile(rootDir, "RESULTS_REPRESENTATIVE", "ordered_heatmaps_fMRI", subjID);
if ~exist(outDir,"dir"), mkdir(outDir); end
outPng = fullfile(outDir, sprintf("%s_fMRI_static_r_ORDERED_%s_N%d.png", subjID, runKey, numel(names)));

plot_heatmap_z_noDiag_withNames_andGroupLines(Z_ord, names_ord, grp(ord), grpLabels, ...
    sprintf("%s fMRI static FC_r (ordered, no diag) | %s | N=%d", subjID, runKey, size(Z_ord,1)), ...
    outPng);

fprintf("[OK] Saved: %s\n", outPng);

%% ===================== helpers =====================

function [ord, grp, grpLabels] = order_ICE060_seednames(names)
% Anatomy-first ordering for ICE060 (bipolar labels like DLA1-DLA2, DRAIN1-DRAIN2, etc.)
% Returns:
%   ord  : permutation indices
%   grp  : group ID per channel in the ORIGINAL order (so use grp(ord) after ordering)
%   grpLabels : containers.Map(groupID -> label)

names = string(names(:));
n = numel(names);

% ---- parse LEFT contact of bipolar name ----
lat   = strings(n,1);   % "L" or "R"
code  = strings(n,1);   % "AIN","PIN","A","AH","PH", etc.
idxNum = nan(n,1);

for i = 1:n
    nm   = upper(names(i));
    left = extractBefore(nm, "-");   % e.g., "DRAH6"

    tok = regexp(left, '^(DR|DL)([A-Z]+)(\d+)$', 'tokens', 'once');
    if isempty(tok)
        lat(i) = "UNK";
        code(i) = "UNK";
        idxNum(i) = inf;
    else
        lat(i)  = string(tok{1});            % "DR" or "DL"
        lat(i)  = replace(lat(i),"DL","L");
        lat(i)  = replace(lat(i),"DR","R");  % now "L" or "R"

        code(i) = string(tok{2});            % "AIN","PIN","A","AH","PH", ...
        idxNum(i) = str2double(tok{3});
    end
end


% ---- map prefix -> anatomy group label ----
% Edit these rules ONCE based on what ICE060 actually has.
% ---- map prefix -> anatomy group label (MORE GRANULAR) ----
anat = strings(n,1);

for i = 1:n
    if lat(i) == "UNK"
        anat(i) = "UNK"; 
        continue;
    end

    c = code(i);

    % Insula
    if c == "AIN"
        anat(i) = lat(i) + "_AIN";
    elseif c == "PIN"
        anat(i) = lat(i) + "_PIN";
    elseif contains(c,"IN")
        anat(i) = lat(i) + "_IN";

    % Temporal
    elseif c == "A"
        anat(i) = lat(i) + "_AMYG";     % DLA*, DRA*
    elseif c == "AH"
        anat(i) = lat(i) + "_HIPP";     % DLAH*, DRAH*
    elseif c == "PH"
        anat(i) = lat(i) + "_PHG";      % DLPH*, DRPH*
    else
        anat(i) = "UNK";
    end
end


% ---- define ordering of anatomy groups (edit if needed) ----
% Typical: insula -> left MTL -> right MTL
grpOrder = ["L_AIN","L_PIN","L_IN", ...
            "L_AMYG","L_HIPP","L_PHG", ...
            "R_AIN","R_PIN","R_IN", ...
            "R_AMYG","R_HIPP","R_PHG", ...
            "UNK"];


% rank each channel by group order
rank = nan(n,1);
for k = 1:numel(grpOrder)
    rank(anat == grpOrder(k)) = k;
end
rank(isnan(rank)) = numel(grpOrder);

% ---- sort ----
T = table(rank, idxNum, (1:n)', anat);
T = sortrows(T, {'rank','idxNum','anat'});
ord = T.Var3;

% ---- group IDs (contiguous 1..G in the ORDERED list) ----
anat_ord = anat(ord);
[uniqGroups,~,grp_ord] = unique(anat_ord, 'stable');

% convert grp_ord (ordered) back to "original index" vector grp
grp = nan(n,1);
grp(ord) = grp_ord;

% labels map (groupID -> groupName)
grpLabels = containers.Map('KeyType','double','ValueType','char');
for g = 1:numel(uniqGroups)
    grpLabels(g) = char(uniqGroups(g));
end
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
