
% HOW TO RUN:
% outDir = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_pairlevel";
% Fig_pairLevel_trueStatic_vs_windowAvg_staticFC(outDir, fullfile(outDir,"Fig3A_trueStatic_vs_windowAvg"), "fdr");


function fig = Fig_pairLevel_trueStatic_vs_windowAvg_staticFC(outDir, outPathBase, corrMode)
% Fig3A final
% Panel A: trueStatic EEG → trueStatic fMRI
% Panel B: windowAvg EEG → windowAvg fMRI
%
% Reads:
%   Fig3A_LMM_trueStatic.csv
%   Fig3A_LMM_windowAvg.csv

if nargin < 2; outPathBase = ""; end
if nargin < 3 || strlength(string(corrMode))==0; corrMode = "fdr"; end
corrMode = lower(string(corrMode));

Ttrue = readtable(fullfile(outDir,"Fig3A_LMM_trueStatic.csv"));
Twin  = readtable(fullfile(outDir,"Fig3A_LMM_windowAvg.csv"));

% shared y-lims (CI + zero)
yMin = min([Ttrue.ci_low;  Twin.ci_low;  0]);
yMax = max([Ttrue.ci_high; Twin.ci_high; 0]);
pad  = 0.12 * (yMax - yMin + eps);
ylims = [yMin-pad, yMax+pad];

fig = figure('Color','w','Position',[80 80 980 420],'Renderer','opengl');
tlo = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','loose');
title(tlo,"Mean EEG FC predicting mean fMRI FC",'FontSize',17,'FontWeight','bold');

nexttile(tlo,1);
plot_coef_panel(Ttrue, "(A) trueStatic EEG \rightarrow trueStatic fMRI", ylims, corrMode);

nexttile(tlo,2);
plot_coef_panel(Twin, "(B) windowAvg EEG \rightarrow windowAvg fMRI", ylims, corrMode);

if strlength(string(outPathBase)) > 0
    exportgraphics(fig, outPathBase + ".png",'Resolution',400);
    try
        exportgraphics(fig, outPathBase + ".pdf",'ContentType','vector');
    catch
        exportgraphics(fig, outPathBase + ".pdf",'ContentType','image','Resolution',400);
    end
end
end

% ---------- helper ----------
function plot_coef_panel(T, panelTitle, ylims, corrMode)

bandOrder  = ["delta","theta","alpha","beta","low_gamma","high_gamma"];
bandLabels = ["\delta","\theta","\alpha","\beta","low \gamma","high \gamma"];

b = string(T.band);
b(b=="gammalow")  = "low_gamma";
b(b=="gammahigh") = "high_gamma";
T.band = categorical(b, cellstr(bandOrder), "Ordinal", true);
T = sortrows(T,"band");

x  = (1:height(T))';
yb = T.beta; yl = T.ci_low; yh = T.ci_high;

% significance
sig = false(height(T),1);
sigNote = "";
switch corrMode
    case "fdr"
        if ismember("p_fdr", T.Properties.VariableNames)
            sig = T.p_fdr < 0.05;
            sigNote = "Significance: BH-FDR across bands (p_{FDR}<0.05)";
        else
            sig = T.p_unc < 0.05;
            sigNote = "Significance: uncorrected p<0.05 (no p_fdr column)";
        end
    case "unc"
        sig = T.p_unc < 0.05;
        sigNote = "Significance: uncorrected p<0.05";
    otherwise
        error("corrMode must be 'fdr' or 'unc' for this figure.");
end

ax = gca; hold(ax,'on');

% highlight significant bands
for i=1:numel(x)
    if sig(i)
        rectangle(ax,'Position',[x(i)-0.45, ylims(1), 0.9, diff(ylims)], ...
            'FaceColor',[0.96 0.96 0.96],'EdgeColor','none');
    end
end

yline(ax,0,'--','LineWidth',1.2,'Color',[0.35 0.35 0.35]);

for i=1:numel(x)
    line(ax,[x(i) x(i)],[yl(i) yh(i)],'LineWidth',2.6,'Color',[0.35 0.35 0.35]);
end

scatter(ax,x(~sig),yb(~sig),80,'filled','MarkerFaceColor',[0.20 0.45 0.85], ...
    'MarkerEdgeColor','k','LineWidth',0.8);
scatter(ax,x(sig), yb(sig), 88,'filled','MarkerFaceColor',[0.90 0.20 0.20], ...
    'MarkerEdgeColor','k','LineWidth',0.8);

yrng = diff(ylims);
for i=1:numel(x)
    if sig(i)
        text(ax,x(i),yh(i)+0.03*yrng,"*", ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',16,'FontWeight','bold','Color',[0.15 0.15 0.15]);
    end
end

ax.XLim=[0.5 numel(x)+0.5];
ax.YLim=ylims;
ax.XTick=x;
ax.XTickLabel=bandLabels;
ax.TickLabelInterpreter='tex';
ax.FontName='Arial'; ax.FontSize=16; ax.LineWidth=1.1;
ax.TickDir='out'; ax.Box='off';

title(ax,panelTitle,'FontSize',14,'FontWeight','bold');
ylabel(ax,'LMM coefficient (\beta)','FontSize',14);

grid(ax,'on'); ax.YGrid='on'; ax.XGrid='off'; ax.GridAlpha=0.10;

text(ax,0.02,0.02,sigNote,'Units','normalized','FontSize',12,'Color',[0.35 0.35 0.35]);
end

