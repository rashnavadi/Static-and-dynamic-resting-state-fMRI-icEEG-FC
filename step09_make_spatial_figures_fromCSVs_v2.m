
%RUN EXAMPLE
% outDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_spatial";
% figDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_spatial/figures";
function step09_make_spatial_figures_fromCSVs_v2(outDir, figDir)

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

    if isempty(regions_ref)
        regions_ref = regions;
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
% 2) Shared color limits
% -----------------------------
allVals = [];
for bi = 1:numel(betaMasked_all)
    M = betaMasked_all{bi};
    if isempty(M); continue; end
    allVals = [allVals; M(~isnan(M))]; %#ok<AGROW>
end
if isempty(allVals)
    clim = [-1 1];
else
    clim = [min(allVals) max(allVals)];
    if clim(1) == clim(2)
        clim = clim + [-1 1]*1e-6; % avoid flat caxis error
    end
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
        title(ax(bi), bandLabel, "Interpreter","tex", "FontWeight","bold");

        continue;
    end

    imagesc(ax(bi), M);
    axis(ax(bi), "image");
    colormap(ax(bi), parula);
    caxis(ax(bi), clim);

    set(ax(bi), "XTick", 1:numel(regions_ref), "XTickLabel", regions_ref, ...
               "YTick", 1:numel(regions_ref), "YTickLabel", regions_ref, ...
               "TickDir","out", "FontSize", 10);
    xtickangle(ax(bi), 45);

    bandLabel = bandLabelMap(band);
    title(ax(bi), bandLabel, "Interpreter","tex", "FontWeight","bold");

    % grid
    hold(ax(bi), "on");
    nR = numel(regions_ref);
    for k = 0.5:1:(nR+0.5)
        plot(ax(bi), [0.5 nR+0.5],[k k],'k-','LineWidth',0.3);
        plot(ax(bi), [k k],[0.5 nR+0.5],'k-','LineWidth',0.3);
    end
    hold(ax(bi), "off");
end

% one shared colorbar for the whole figure
cb = colorbar(ax(1));
cb.Layout.Tile = "east";
ylabel(cb, "\beta(\Delta EEG FC)", "FontWeight","bold");

title(tl, sprintf("Spatially stratified dynamic EEG–fMRI coupling (masked p<%.2f)", pThresh), ...
    "FontWeight","bold", "FontSize", 14);

outPngAll = fullfile(figDir, "CouplingHeatmap_ALLBANDS_BETAdev_p" + strrep(num2str(pThresh),".","p") + ".png");
exportgraphics(fAll, outPngAll, "Resolution", 300);
close(fAll);

fprintf("[OK] Saved ONE combined figure here:\n%s\n", outPngAll);

end
