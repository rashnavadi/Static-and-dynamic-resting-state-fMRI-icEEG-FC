
%RUN EXAMPLE
% outDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_spatial";
% figDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_spatial/figures";
function step09_make_spatial_figures_fromCSVs_v2_ROIsreordered(outDir, figDir)
% ---------- Figure font sizes (paper-ready) ----------
fs_roi   = 14;   % ROI tick labels
fs_title = 16;   % per-band titles
fs_super = 20;   % main figure title
fs_cb    = 16;   % colorbar label + ticks


if nargin < 1 || isempty(outDir)
    outDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_spatial";
end
if nargin < 2 || isempty(figDir)
    figDir = fullfile(outDir, "figs_coupling_only");
end
if ~exist(figDir,"dir"); mkdir(figDir); end

bandOrder = ["delta","theta","alpha","beta","low_gamma","high_gamma"];
pThresh = 0.05;

bandLabelMap = containers.Map( ...
    ["delta","theta","alpha","beta","low_gamma","high_gamma"], ...
    ["\delta","\theta","\alpha","\beta","low \gamma","high \gamma"] );


betaMasked_all = cell(numel(bandOrder),1);
regions_ref = [];

% -----------------------------
% 1) Load + mask all bands
% -----------------------------
for bi = 1:numel(bandOrder)
    band = string(bandOrder(bi));
    betaPath = fullfile(outDir, "BETA_dev_" + band + ".csv");
    pPath    = fullfile(outDir, "P_dev_"    + band + ".csv");

    if ~exist(betaPath,"file") || ~exist(pPath,"file")
        fprintf("[SKIP] Missing beta/p files for %s\n", band);
        continue;
    end

    Tbeta = readtable(betaPath, "ReadRowNames", true);
    Tp    = readtable(pPath,    "ReadRowNames", true);

    regions = string(Tbeta.Properties.VariableNames);
    betaMat = table2array(Tbeta);
    pMat    = table2array(Tp);

    % ---------- REORDER REGIONS (Left block then Right block) ----------
    regions_ordered = ["LF","LMT","LLT","LI","RF","RMT","RLT","RI"];
    [tf, idx] = ismember(regions_ordered, regions);
    assert(all(tf), "Some desired regions not found in CSV labels: %s", strjoin(regions_ordered(~tf), ", "));

    betaMat = betaMat(idx, idx);
    pMat    = pMat(idx, idx);
    regions = regions_ordered;   % overwrite for consistency
    % ---------------------------------------------------------------

    if isempty(regions_ref)
        regions_ref = regions;   % now regions_ref is the ORDERED one
    end

    betaMasked = betaMat;
    betaMasked(pMat > pThresh) = NaN;

    betaMasked_all{bi} = betaMasked;
end

% sanity
if all(cellfun(@isempty, betaMasked_all))
    error("No band matrices loaded. Check outDir and filenames BETA_dev_*.csv / P_dev_*.csv");
end

% -----------------------------
% 2) Shared color limits (ROBUST, zero-centered) — more contrast
% -----------------------------
allVals = [];
for bi = 1:numel(betaMasked_all)
    M = betaMasked_all{bi};
    if isempty(M); continue; end
    allVals = [allVals; M(~isnan(M))]; %#ok<AGROW>
end

if isempty(allVals)
    clim = [-0.2 0.2];   % fallback
else
    absVals = abs(allVals);

    % Pick a percentile that gives more contrast (try 85–95)
    lim = prctile(absVals, 65);

    if lim == 0
        lim = 1e-6;
    end
    clim = [-lim lim];
end

% -----------------------------
% 3) One single combined figure
% -----------------------------
fAll = figure("Color","w", "Position",[50 50 1700 950]);
tl = tiledlayout(fAll, 2, 3, "Padding","compact", "TileSpacing","compact");

ax = gobjects(numel(bandOrder),1);

for bi = 1:numel(bandOrder)
    band = string(bandOrder(bi));
    M = betaMasked_all{bi};

    ax(bi) = nexttile(tl);

    if isempty(M)
        axis off;
        bandLabel = bandLabelMap(band);
        title(ax(bi), bandLabel, ...
            "Interpreter","tex", ...
            "FontWeight","bold", ...
            "FontSize", fs_title);
        continue;
    end

    imagesc(ax(bi), M);
    axis(ax(bi), "image");
    colormap(ax(bi), parula);
    caxis(ax(bi), clim);
    set(ax(bi), "Box","on");
    
    set(ax(bi), ...
        "XTick", 1:numel(regions_ref), ...
        "XTickLabel", regions_ref, ...
        "YTick", 1:numel(regions_ref), ...
        "YTickLabel", regions_ref, ...
        "TickDir","out", ...
        "FontSize", fs_roi, ...
        "LineWidth", 1.2);

    xtickangle(ax(bi), 45);


    xtickangle(ax(bi), 45);

    bandLabel = bandLabelMap(band);
    title(ax(bi), bandLabel, ...
        "Interpreter","tex", ...
        "FontWeight","bold", ...
        "FontSize", fs_title);


    % grid
    hold(ax(bi), "on");
    nR = numel(regions_ref);
    for k = 0.5:1:(nR+0.5)
        plot(ax(bi), [0.5 nR+0.5],[k k],'k-','LineWidth',0.3);
        plot(ax(bi), [k k],[0.5 nR+0.5],'k-','LineWidth',0.3);
    end
    hold(ax(bi), "off");
    split = 4;
    xline(ax(bi), split+0.5, 'k-', 'LineWidth', 2);
    yline(ax(bi), split+0.5, 'k-', 'LineWidth', 2);

end

% one shared colorbar for the whole figure
cb = colorbar(ax(1));
cb.Layout.Tile = "east";

cb.FontSize = fs_cb;
cb.LineWidth = 1.2;

ylabel(cb, "\beta(\Delta EEG FC)", ...
    "Interpreter","tex", ...
    "FontWeight","bold", ...
    "FontSize", fs_cb);


title(tl, ...
    sprintf("Spatially stratified dynamic EEG–fMRI coupling (masked p < %.2f)", pThresh), ...
    "FontWeight","bold", ...
    "FontSize", fs_super);


outPngAll = fullfile(figDir, "CouplingHeatmap_ALLBANDS_BETAdev_p" + strrep(num2str(pThresh),".","p") + ".png");
exportgraphics(fAll, outPngAll, "Resolution", 300);
close(fAll);

fprintf("[OK] Saved ONE combined figure here:\n%s\n", outPngAll);

end

function cmap = bluewhitered(m)
%BLUEWHITERED  Blue-white-red diverging colormap
if nargin < 1
    m = size(get(gcf,'colormap'),1);
end

bottom = [0 0 0.5];
middle = [1 1 1];
top    = [0.5 0 0];

% interpolate
cmap = zeros(m,3);
for i = 1:3
    cmap(:,i) = interp1([1 m/2 m],[bottom(i) middle(i) top(i)],1:m);
end
end

