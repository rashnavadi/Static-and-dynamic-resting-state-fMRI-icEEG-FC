function plot_dynamic_FC(subjectList, rootDir, runID, pairID, outPng)
% plot_dynamic_FC
% Recreates OLD from Kristina "Figure 5" style: dynamic fluctuations of windowed fMRI FC (dashed)
% vs band-limited icEEG FC (solid) for ONE given connection across time windows.
%
% INPUTS
%   subjectList : string array/cellstr, e.g. ["ICE013","ICE017","ICE045","ICE047"]
%   rootDir     : directory that contains subject folders with master_dynamicFC.csv
%                e.g. "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables"
%   runID       : run string (exact match), e.g. "Run1a"
%                (use "" to include all runs, but usually pick one)
%   pairID      : exact pair_id string, e.g. "DLA1-DLA2--DLH1-DLH2"
%   outPng      : optional output PNG path ("" to not save)
%
% EXAMPLE:
% subjects = ["ICE013","ICE017","ICE045","ICE047"];
% subjects = ["ICE060"];
% subjects = ["ICE017"];
% subjects = ["ICE044"];
% subjects = ["ICE024"];
% rootDir  = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables";
% runID    = "Run1";
% runID    = "Run1a";
% pairID   = "DLA1-DLA2--DLA3-DLA4";
% pairID   = "DRAIN5-DRAIN6--DRLOF3-DRLOF4";
% pairID   = "DLOF5-DLOF6--DROF3-DROF4";
% pairID   = "DLA2-DLA3--DRA7-DRA8";
% pairID   = "DLAC1-DLAC2--DLPIN5-DLPIN6";
% outPng   = fullfile(rootDir,"ICE024_dynamic_example.png");
% plot_dynamic_FC(subjects, rootDir, runID, pairID, outPng);

if nargin < 5, outPng = ""; end

% --- band order (matches your convention) ---
bandOrder = ["delta","theta","alpha","beta","low_gamma","high_gamma"];
bandLabelMap = containers.Map( ...
    ["delta","theta","alpha","beta","low_gamma","high_gamma"], ...
    ["\delta","\theta","\alpha","\beta","low \gamma","high \gamma"] );

nSubj = numel(subjectList);
nBand = numel(bandOrder);

% --- figure settings ---
fsTitle = 14;
fsTick  = 10;
fsLab   = 12;

f = figure("Color","w","Position",[50 50 1400 250*nSubj]);
tl = tiledlayout(f, nSubj, nBand, "Padding","compact", "TileSpacing","compact");
title(tl, "Dynamic fluctuations of windowed fMRI FC (dashed) and icEEG FC (solid)", ...
    "FontWeight","bold","FontSize",16);

% --- collect global y-limits for clean comparability ---
allY = [];

% pass 1: read + store per-subject tables (filtered) for speed and global y-lims
S = struct();
for si = 1:nSubj
    subj = string(subjectList(si));
    csvPath = fullfile(rootDir, subj, "master_dynamicFC.csv");
    if ~exist(csvPath,"file")
        warning("[SKIP] Missing: %s", csvPath);
        S(si).ok = false; %#ok<*AGROW>
        continue;
    end

    T = readtable(csvPath);

    % enforce string columns (safer)
    T.subject_id = string(T.subject_id);
    T.run_id     = string(T.run_id);
    T.pair_id    = string(T.pair_id);
    T.band       = string(T.band);

    % filter: subject, pair, (optional) run
    T = T(T.subject_id == subj & T.pair_id == string(pairID), :);
    if strlength(string(runID)) > 0
        T = T(T.run_id == string(runID), :);
    end

    if isempty(T)
        warning("[SKIP] No rows for %s (%s, %s)", subj, runID, pairID);
        S(si).ok = false;
        continue;
    end

    % sort by band then window
    T = sortrows(T, {'band','win_idx'});   % <-- contains chars (' ')

    S(si).ok = true;
    S(si).subj = subj;
    S(si).T = T;

    allY = [allY; T.fmri_fc_windowed; T.eeg_fc_windowed];
end

if isempty(allY)
    error("No data found after filtering. Check runID/pairID and file paths.");
end

% robust global y-limits (avoid a single extreme window dominating)
lo = prctile(allY, 2);
hi = prctile(allY, 98);
pad = 0.08 * (hi - lo);
if pad == 0, pad = 0.05; end
yLim = [lo-pad, hi+pad];

% pass 2: plot
for si = 1:nSubj
    for bi = 1:nBand
        ax = nexttile(tl);

        if ~S(si).ok
            axis(ax,"off");
            continue;
        end

        subj = S(si).subj;
        T = S(si).T;

        band = bandOrder(bi);
        Tb = T(T.band == band, :);

        if isempty(Tb)
            axis(ax,"off");
            continue;
        end

        x = Tb.win_idx;
        y_fmri = Tb.fmri_fc_windowed;
        y_eeg  = Tb.eeg_fc_windowed;

        % plot (match your description)
        plot(ax, x, y_eeg, "-",  "LineWidth", 1.8); hold(ax,"on");
        plot(ax, x, y_fmri, "--", "LineWidth", 1.8);

        % correlation (optional but very informative)
        r = corr(y_eeg, y_fmri, "Rows","complete");
        txt = sprintf("r = %.2f", r);

        % style
        ax.YLim = yLim;
        ax.FontSize = fsTick;
        ax.TickDir = "out";
        grid(ax,"on");
        ax.GridAlpha = 0.15;
        box(ax,"off");

        % titles: top row only
        if si == 1
            title(ax, bandLabelMap(band), "Interpreter","tex", ...
                "FontWeight","bold","FontSize",fsTitle);
        end

        % subject labels: first column only
        if bi == 1
            ylabel(ax, subj, "FontWeight","bold","FontSize",fsLab);
        end

        % x label: bottom row only
        if si == nSubj
            xlabel(ax, "window", "FontSize",fsLab);
        else
            ax.XTickLabel = [];
        end

        % show r in each panel (top-left corner)
        xText = min(x) + 0.05*(max(x)-min(x));
        yText = yLim(2) - 0.12*(yLim(2)-yLim(1));
        text(ax, xText, yText, txt, "FontSize",10, "FontWeight","bold");

        % legend only once
        if si == 1 && bi == nBand
            legend(ax, {"icEEG FC","fMRI FC"}, "Location","southoutside", "Box","off");
        end
    end
end

% save
if strlength(string(outPng)) > 0
    exportgraphics(f, outPng, "Resolution", 300);
    fprintf("[OK] Saved: %s\n", outPng);
end

end
