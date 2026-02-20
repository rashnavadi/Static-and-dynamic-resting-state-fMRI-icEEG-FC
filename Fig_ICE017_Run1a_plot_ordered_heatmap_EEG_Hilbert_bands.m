%% Fig_Run1a_ICE017_plot_ordered_heatmap_EEG_Hilbert_bands.m
clear; clc;

subjID  = "ICE017";
runKey  = "Run1a";


rootDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/FC_rebuild";
subjDir = fullfile(rootDir, subjID);

eegFile = fullfile(subjDir, "EEG_FC_Hilbert_" + subjID + ".mat");
load(eegFile, "OUT");

% ----- bands to plot -----
bandOrder = ["delta","theta","alpha","beta","low_gamma","high_gamma"];

% ----- get names (channel labels) -----
names = get_run_channel_names(OUT.(runKey));
names = string(names);

% ----- build ordering (same as fMRI) -----
[ord, grp, grpLabels] = order_ICE017_seednames(names);

names_ord = names(ord);

% ----- output dir -----
outDir = fullfile(rootDir, "RESULTS_REPRESENTATIVE", "ordered_heatmaps_EEG_Hilbert", subjID);
if ~exist(outDir,"dir"), mkdir(outDir); end

% ============================== 
% ===== ALSO make one combined panel (2x3) =====
panelPng = fullfile(outDir, sprintf("%s_EEG_Hilbert_static_r_ORDERED_%s_ALLBANDS_PANEL_N%d.png", ...
    subjID, runKey, numel(names)));
panelPdf = erase(panelPng, ".png") + ".pdf";

fs = 6;              % tick label font size (panel needs small text)
clim = [-0.5 0.5];        % keep consistent with your per-band plots

figP = figure("Color","w","Position",[50 50 1700 900],"Visible","off");
tl = tiledlayout(figP, 2, 3, "Padding","compact", "TileSpacing","compact");
% tl = tiledlayout(figP, 2, 3, "Padding","tight", "TileSpacing","tight");
title(tl, sprintf("%s EEG Hilbert static FC_r (ordered, no diag) | %s | ALL BANDS | N=%d", ...
    subjID, runKey, numel(names)), "Interpreter","none", "FontWeight","bold");

axList = gobjects(numel(bandOrder),1);
hmList = gobjects(numel(bandOrder),1);
% ============================== 

% ----- plot each band -----
for b = 1:numel(bandOrder)
    band = bandOrder(b);

    [R, ok] = get_eeg_band_matrix(OUT.(runKey), band);
    if ~ok
        fprintf("[MISS] %s %s: band %s not found in EEG OUT.\n", subjID, runKey, band);
        continue;
    end

    % reorder
    R_ord = R(ord, ord);
    % ============================== 
    % ---- also draw into the combined panel ----
    ax = nexttile(tl);
    axList(b) = ax;

    % remove diagonal for display
    Rp = R_ord;
    n = size(Rp,1);
    Rp(1:n+1:end) = NaN;

    hm = imagesc(ax, Rp);
    hmList(b) = hm;

    axis(ax,"image");
    colormap(ax, parula);
    caxis(ax, clim);

    % alpha mask for NaNs
    set(hm, "AlphaData", ~isnan(Rp));

    % band title
    ttl = "$" + band_to_greek_title(band) + "$";
    title(ax, char(ttl), "Interpreter","latex", "FontWeight","bold", "FontSize", 30);


    % ---- tick formatting ----
    set(ax, ...
        "XTick", 1:n, ...
        "YTick", 1:n, ...
        "TickLabelInterpreter","none", ...
        "FontSize", fs, ...
        "TickDir","out", ...
        "LineWidth", 0.8);

    xtickangle(ax, 90);

    % Determine row (2 rows x 3 columns)
    if b <= 3
        row = 1;
    else
        row = 2;
    end

    col = mod(b-1,3) + 1;

    % -------- X TICKS --------
    if row == 1
        % FIRST ROW → show labels on TOP only
        ax.XAxisLocation = 'top';
        ax.XTickLabel = names_ord;
        ax.XTickLabelRotation = 90;
    else
        % SECOND ROW → show labels on bottom only
        ax.XAxisLocation = 'bottom';
        ax.XTickLabel = names_ord;
        ax.XTickLabelRotation = 90;
    end

    % -------- Y TICKS --------
    % Show Y labels only on first column of each row
    if col == 1
        ax.YTickLabel = names_ord;
    else
        ax.YTickLabel = repmat({''}, n, 1);
    end


    % group boundary lines (same logic as your helper)
    hold(ax,"on");
    boundaries = find(diff(grp(ord)) ~= 0);
    for bb = boundaries'
        x = bb + 0.5;
        plot(ax, [0.5, n+0.5], [x, x], 'k-', 'LineWidth', 1.0);
        plot(ax, [x, x], [0.5, n+0.5], 'k-', 'LineWidth', 1.0);
    end
    hold(ax,"off");
    % ==============================

    outPng = fullfile(outDir, sprintf("%s_EEG_Hilbert_static_r_ORDERED_%s_%s_N%d.png", ...
        subjID, runKey, band, numel(names)));

    ttl = sprintf("%s EEG Hilbert static FC_r (ordered, no diag) | %s | %s | N=%d", ...
        subjID, runKey, band, size(R_ord,1));

    plot_heatmap_r_noDiag_withNames_andGroupLines(R_ord, names_ord, grp(ord), grpLabels, ttl, outPng);

    fprintf("[OK] Saved: %s\n", outPng);
end

% ==============================
% ---- one shared colorbar for the whole panel ----
cb = colorbar(axList(1));
cb.Layout.Tile = "east";                 % puts it on the right of the tiledlayout
cb.TickDirection = "out";
cb.FontSize = 14;
ylabel(cb, "FC (r)", "FontWeight","bold");

% export panel
exportgraphics(figP, panelPng, "Resolution", 300);
try
    exportgraphics(figP, panelPdf, "ContentType","vector");
catch
    exportgraphics(figP, panelPdf, "ContentType","image", "Resolution", 300);
end
close(figP);

fprintf("[OK] Saved combined panel:\n%s\n", panelPng);
 % ==============================


%% ===================== helpers =====================

function names = get_run_channel_names(runStruct)
    % Try common name fields used in your pipeline
    cand = ["seedNames","chanNames","channelNames","bipolarNames","bipNames","labels"];
    f = string(fieldnames(runStruct));
    names = string.empty;

    for c = cand
        if any(f == c)
            names = string(runStruct.(c));
            names = names(:);
            return;
        end
    end

    % fallback
    if isfield(runStruct,"FC_static") && isnumeric(runStruct.FC_static)
        n = size(runStruct.FC_static,1);
    else
        n = 0;
    end
    names = "Ch" + string(1:max(n,1));
end

function [R, ok] = get_eeg_band_matrix(runStruct, band)
    ok = false;
    R  = [];

    band = lower(string(band));

    % Map your canonical band names to possible variants in saved struct
    variants = band_variants(band);

    % 1) direct fields like FC_static_delta
    fn = string(fieldnames(runStruct));
    for v = variants
        cand = "FC_static_" + v;
        if any(lower(fn) == lower(cand))
            R = runStruct.(fn(lower(fn) == lower(cand)));
            ok = true; return;
        end
    end

    % 2) containers that might hold band subfields
    containers = ["FC_static","FC_static_r","FC","R","R_static","FC_r"];
    for c = containers
        if isfield(runStruct, c)
            X = runStruct.(c);

            % If the container itself is the matrix (single band)
            if isnumeric(X) && ismatrix(X) && size(X,1) == size(X,2)
                R = X; ok = true; return;
            end

            % If the container is a struct, look for band fields inside (case-insensitive)
            if isstruct(X)
                xfn = string(fieldnames(X));
                for v = variants
                    hit = find(lower(xfn) == lower(v), 1);
                    if ~isempty(hit)
                        R = X.(xfn(hit));
                        ok = true; return;
                    end
                end
            end
        end
    end

    % 3) sometimes bands are top-level fields: runStruct.alpha, runStruct.beta, ...
    for v = variants
        hit = find(lower(fn) == lower(v), 1);
        if ~isempty(hit)
            X = runStruct.(fn(hit));
            % could be matrix or struct with FC field
            if isnumeric(X) && ismatrix(X) && size(X,1) == size(X,2)
                R = X; ok = true; return;
            elseif isstruct(X)
                if isfield(X,"FC_static")
                    R = X.FC_static; ok = true; return;
                elseif isfield(X,"FC")
                    R = X.FC; ok = true; return;
                end
            end
        end
    end
end

function vars = band_variants(band)
    % returns lowercase strings
    switch band
        case "delta"
            vars = ["delta","d"];
        case "theta"
            vars = ["theta","t"];
        case "alpha"
            vars = ["alpha","a"];
        case "beta"
            vars = ["beta","b"];
        case "low_gamma"
            vars = ["low_gamma","lowgamma","lowGamma","gamma_low","lgamma","lg", "gammalow","gammaLow","lowg","low_g"];
        case "high_gamma"
            vars = ["high_gamma","highgamma","highGamma","gamma_high","hgamma","hg", "gammahigh","gammaHigh","highg","high_g"];
        otherwise
            vars = band;
    end
end


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

    grp = regRank;

    T = table(regRank, idxNum, (1:n)', region);
    T = sortrows(T, {'regRank','idxNum','region'});
    ord = T.Var3;

    grpLabels = containers.Map('KeyType','double','ValueType','char');
    grpLabels(1) = 'MOF';
    grpLabels(2) = 'LOF';
    grpLabels(3) = 'FOP';
    grpLabels(4) = 'AIN';
    grpLabels(5) = 'PIN';
    grpLabels(6) = 'PT';
end

function plot_heatmap_r_noDiag_withNames_andGroupLines(R, names, grp, grpLabels, ttl, outPng)
    names = string(names);
    n = size(R,1);

    % remove diagonal
    R(1:n+1:end) = NaN;

    fig = figure("Visible","off");
    h = imagesc(R);
    axis image;

    cb = colorbar;
    cb.Location = 'eastoutside';
    cb.TickDirection = 'out';

    ax = gca;
    ax.Position(3) = ax.Position(3) * 0.88;

    title(ttl, "Interpreter","none");

    set(gca, "XTick", 1:n, "XTickLabel", names, ...
             "YTick", 1:n, "YTickLabel", names, ...
             "TickLabelInterpreter","none", ...
             "FontSize", 10);
    xtickangle(90);

    set(h, "AlphaData", ~isnan(R));

    % fixed r range like your fMRI r plots
    caxis([-1 1]);

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

function ttl = band_to_greek_title(band)
band = string(lower(band));
switch band
    case "delta"
        ttl = "\delta";
    case "theta"
        ttl = "\theta";
    case "alpha"
        ttl = "\alpha";
    case "beta"
        ttl = "\beta";
    case "low_gamma"
        ttl = "Low \gamma";
    case "high_gamma"
        ttl = "High \gamma";
    otherwise
        ttl = band;
end
end



