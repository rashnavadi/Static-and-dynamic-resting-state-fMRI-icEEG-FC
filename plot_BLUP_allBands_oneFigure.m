% Run e.g.:
% plot_BLUP_subjects("/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_window/BLUPs/BLUP_deltaEEG_alpha.csv")
% plot_BLUP_subjects("/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_window/BLUPs/BLUP_deltaEEG_beta.csv")
% plot_BLUP_subjects("/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_window/BLUPs/BLUP_deltaEEG_delta.csv")
% plot_BLUP_subjects("/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_window/BLUPs/BLUP_deltaEEG_theta.csv")
% plot_BLUP_subjects("/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_window/BLUPs/BLUP_deltaEEG_low_gamma.csv")
% plot_BLUP_subjects("/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_window/BLUPs/BLUP_deltaEEG_high_gamma.csv")


function plot_BLUP_allBands_oneFigure(blupDirOrFiles, outPng)
% plot_BLUP_allBands_oneFigure
% One combined plot of subject-specific dynamic coupling (BLUP slopes) across bands.
%
% INPUTS
%   blupDirOrFiles: either a directory containing BLUP_deltaEEG_*.csv
%                   OR a string array/cell array of full CSV paths
%   outPng        : optional output PNG path (e.g., ".../BLUP_allBands.png")
%
% EXAMPLE (directory):
%   plot_BLUP_allBands_oneFigure("/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_window/BLUPs", "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_window/BLUPs/BLUP_allBands.png")

%
% EXAMPLE (explicit files):
%   files = [
%     "/.../BLUP_deltaEEG_delta.csv"
%     "/.../BLUP_deltaEEG_theta.csv"
%     "/.../BLUP_deltaEEG_alpha.csv"
%     "/.../BLUP_deltaEEG_beta.csv"
%     "/.../BLUP_deltaEEG_low_gamma.csv"
%     "/.../BLUP_deltaEEG_high_gamma.csv"
%   ];
%   plot_BLUP_allBands_oneFigure(files, "/.../BLUP_allBands.png")

if nargin < 2, outPng = ""; end

% --------- Resolve input files ----------
if isstring(blupDirOrFiles) || ischar(blupDirOrFiles)
    p = string(blupDirOrFiles);
    if isfolder(p)
        D = dir(fullfile(p, "BLUP_deltaEEG_*.csv"));
        files = string(fullfile({D.folder}, {D.name}));
    else
        files = string(blupDirOrFiles);
    end
else
    files = string(blupDirOrFiles);
end
if isempty(files), error("No BLUP_deltaEEG_*.csv files found."); end

% --------- Desired band order ----------
bandOrder = ["delta","theta","alpha","beta","low_gamma","high_gamma"];
bandLabel = containers.Map( ...
    ["delta","theta","alpha","beta","low_gamma","high_gamma"], ...
    ["\delta","\theta","\alpha","\beta","low \gamma","high \gamma"] );

% --------- Load + stack ----------
allBand = strings(0,1);
allSlope = [];
allSubj = strings(0,1);

for i = 1:numel(files)
    T = readtable(files(i));

    if ~any(strcmp(T.Properties.VariableNames, "term"))
        warning("No 'term' column in %s. Skipping.", files(i));
        continue;
    end

    % Keep only ΔEEG main effect
    T = T(strcmp(string(T.term), "delta_eeg_fc"), :);
    if isempty(T), continue; end

    if ~any(strcmp(T.Properties.VariableNames, "deltaEEG_subjectSlope"))
        error("Missing column 'deltaEEG_subjectSlope' in %s", files(i));
    end
    if ~any(strcmp(T.Properties.VariableNames, "band"))
        error("Missing column 'band' in %s", files(i));
    end
    if any(strcmp(T.Properties.VariableNames, "subject_id"))
        subj = string(T.subject_id);
    else
        subj = repmat("subj", height(T), 1);
    end

    allBand  = [allBand;  string(T.band)];
    allSlope = [allSlope; T.deltaEEG_subjectSlope];
    allSubj  = [allSubj;  subj];
end

% Keep only known bands
tf = ismember(allBand, bandOrder);
allBand  = allBand(tf);
allSlope = allSlope(tf);
allSubj  = allSubj(tf);

if isempty(allSlope)
    error("After filtering, no ΔEEG rows remain. Check term labels in your CSVs.");
end

% --------- Plot ----------
f = figure('Color','w','Position',[100 100 1100 520]);
ax = axes(f); hold(ax,'on');

% Band colors (muted, paper-friendly)
bandColors = containers.Map( ...
    ["delta","theta","alpha","beta","low_gamma","high_gamma"], ...
    { ...
      [0.40 0.60 0.80], ...
      [0.50 0.70 0.65], ...
      [0.70 0.60 0.80], ...
      [0.85 0.65 0.40], ...
      [0.85 0.45 0.45], ...
      [0.55 0.55 0.55] ...
    });

xPos = 1:numel(bandOrder);

% Symmetric y-lims around 0 (cleaner)
lim = max(abs(allSlope)) * 1.08;
if lim == 0, lim = 0.1; end

rng(0); % set once (don’t reset inside loop)

for b = 1:numel(bandOrder)

    band = bandOrder(b);
    idx = allBand == band;
    v = allSlope(idx);
    if isempty(v); continue; end

    baseColor = bandColors(band);
    fillColor = baseColor + (1 - baseColor)*0.55;  % lighter density color

    % ---- Density (half violin) ----
    [fden, yy] = ksdensity(v, 'NumPoints', 200);
    fden = fden / max(fden);
    width = 0.25;

    yy = yy(:);
    xx = (xPos(b) - width * fden(:));       % left half

    x1 = xPos(b) * ones(numel(yy),1);       % center line
    patch(ax, [x1; flipud(xx)], [yy; flipud(yy)], ...
        fillColor, 'EdgeColor','none', 'FaceAlpha',0.65);

    % ---- Jittered dots ----
    jitter = (rand(size(v)) - 0.5) * 0.18;
    scatter(ax, xPos(b) + jitter, v, 42, ...
        'filled', ...
        'MarkerFaceColor', baseColor, ...
        'MarkerEdgeColor', [0.15 0.15 0.15], ...
        'MarkerFaceAlpha', 0.92);

    % ---- Median + IQR ----
    med = median(v);
    q1 = prctile(v,25);
    q3 = prctile(v,75);

    plot(ax, [xPos(b)-0.18 xPos(b)+0.18], [med med], '-', ...
        'Color',[0.2 0.2 0.2], 'LineWidth', 2.2);
    plot(ax, [xPos(b) xPos(b)], [q1 q3], '-', ...
        'Color',[0.2 0.2 0.2], 'LineWidth', 2.2);
end

yline(ax, 0, '--', 'LineWidth', 1.6, 'Color', [0.45 0.45 0.45]);

% ---- Axis styling ----
ax.XLim = [0.5 numel(bandOrder)+0.5];
ax.YLim = [-lim lim];

ax.FontSize = 14;
ax.LineWidth = 1.2;
ax.TickDir = 'out';
box(ax,'off');
grid(ax,'on');
ax.GridAlpha = 0.15;

% X ticks (Greek)
xt = strings(size(bandOrder));
for b = 1:numel(bandOrder)
    xt(b) = bandLabel(bandOrder(b));
end
set(ax, 'XTick', xPos, 'XTickLabel', xt, 'TickLabelInterpreter','tex');

ylabel(ax, "\DeltaEEG \rightarrow fMRI coupling (subject-specific slope)", 'Interpreter','tex');
title(ax, "Subject-specific dynamic coupling across frequency bands", 'FontWeight','bold');

% Export
if strlength(outPng) > 0
    exportgraphics(f, outPng, 'Resolution', 300);
    fprintf("[OK] Saved: %s\n", outPng);
    close(f);
end
end


