function Fig_window_level_DeltaEEG_randomSlope_dynamicFC(outDir, outPathBase, corrMode)
% Fig_DeltaEEG_randomSlope_DEFAULT
% Single-panel figure for the DEFAULT window-level model (subject-specific random slope for ΔEEG)
% across EEG frequency bands.
%
% Uses (default-model output):
%   LMM_window_delta_randomSlopeSubject_DEFAULT_fdrAcrossBands.csv
%
% Figure:
%   x-axis: bands (δ, θ, α, β, low γ, high γ)
%   y-axis: β for ΔEEG term
%   Filled markers + (vertical) 95% CI
%   Asterisks: significance after BH-FDR across bands (default) or uncorrected p<0.05
%
% How to run:
% outDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_window";
% Fig_window_level_DeltaEEG_randomSlope_dynamicFC(outDir);

if nargin < 1 || strlength(string(outDir))==0
    error("outDir is required.");
end
if nargin < 2; outPathBase = ""; end
if nargin < 3 || strlength(string(corrMode))==0; corrMode = "fdr"; end
corrMode = lower(string(corrMode));

% --- input file (DEFAULT model) ---
inDefault = fullfile(outDir, "LMM_window_delta_randomSlopeSubject_DEFAULT_fdrAcrossBands.csv");
assert(exist(inDefault,"file")==2, "Missing file: %s", inDefault);

T = readtable(inDefault);

% --- band ordering ---
bandOrder  = ["delta","theta","alpha","beta","low_gamma","high_gamma"];
bandLabels = ["\delta","\theta","\alpha","\beta","low \gamma","high \gamma"]; % TeX-safe

T.band = categorical(string(T.band), cellstr(bandOrder), "Ordinal", true);
T = sortrows(T, "band");

% --- extract ΔEEG term ---
yb   = T.beta_dev;
yl   = T.ci_low_dev;
yh   = T.ci_high_dev;

punc = T.p_dev;
pfdr = T.p_dev_fdr;

% --- significance mask for stars ---
switch corrMode
    case "fdr"
        sig = pfdr < 0.05;
        sigLabel = "BH-FDR across bands (p_{FDR}<0.05)";
    otherwise
        sig = punc < 0.05;
        sigLabel = "Uncorrected p<0.05";
end

% --- y-limits from CIs ---
ylims = make_ylims([yl; yh; 0]);

% --- plot ---
fig = figure('Color','w','Position',[120 120 900 520], 'Renderer','opengl');
ax = axes(fig); hold(ax,'on');

title(ax, "\DeltaEEG \rightarrow windowed fMRI FC", ...
    'FontSize',17,'FontWeight','bold');

x = (1:numel(bandOrder))';

% Light alternating vertical bands (optional aesthetic)
for i = 1:numel(x)
    if mod(i,2)==0
        rectangle(ax,'Position',[x(i)-0.5, ylims(1), 1.0, diff(ylims)], ...
            'FaceColor',[0.97 0.97 0.97], 'EdgeColor','none');
    end
end

% Zero line
yline(ax,0,'--','LineWidth',1.2,'Color',[0.35 0.35 0.35]);

% Error bars (manual)
for i = 1:numel(x)
    line(ax,[x(i) x(i)],[yl(i) yh(i)],'LineWidth',2.6,'Color',[0.25 0.25 0.25]);
end

% Markers (filled)
scatter(ax, x, yb, 95, 'filled', ...
    'MarkerFaceColor',[0.85 0.10 0.10], 'MarkerEdgeColor','k','LineWidth',0.9);

% Stars above upper CI (if significant)
yrng = diff(ylims);
for i = 1:numel(x)
    if sig(i)
        text(ax, x(i), yh(i)+0.03*yrng, "*", ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',16,'FontWeight','bold','Color',[0.15 0.15 0.15]);
    end
end

% Axes formatting
ax.XLim = [0.5 numel(x)+0.5];
ax.YLim = ylims;
ax.XTick = x;
ax.XTickLabel = bandLabels;
ax.TickLabelInterpreter = 'tex';
ax.FontName = 'Arial';
ax.FontSize = 16;
ax.LineWidth = 1.1;
ax.TickDir = 'out';
ax.Box = 'off';
ylabel(ax, "LMM coefficient (\beta) for \DeltaEEG", 'FontSize',16);

grid(ax,'on');
ax.YGrid = 'on';
ax.XGrid = 'off';
ax.GridAlpha = 0.10;

% Significance note
text(ax, 0.02, 0.02, "Significance: " + sigLabel, ...
    'Units','normalized','FontSize',12,'Color',[0.35 0.35 0.35]);

% Export
if strlength(outPathBase)==0
    outPathBase = fullfile(outDir, "Fig_DeltaEEG_randomSlope_DEFAULT");
end

exportgraphics(fig, outPathBase + ".png", 'Resolution', 400);
try
    exportgraphics(fig, outPathBase + ".pdf", 'ContentType','vector');
catch
    exportgraphics(fig, outPathBase + ".pdf", 'ContentType','image', 'Resolution', 400);
end

end

% ---- helper: y-lims ----
function ylims = make_ylims(vals)
vals = vals(~isnan(vals));
yMin = min(vals); yMax = max(vals);
pad  = 0.12*(yMax - yMin + eps);
ylims = [yMin-pad, yMax+pad];
end
