%% Supplementary figure: Motion robustness of EEG–fMRI coupling
% One-panel plot: beta(ΔEEG FC) across bands
% Low-motion vs moderate-motion with 95% CI

% -------- Paths --------
inCSV = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_motion/MotionRobustness_LMM_byGroup_thr0p5.csv";
outFig = "/Volumes/MIND/ICE/Tara/Kristina_paper/master_tables/GROUP_MASTER/group_stats_outputs_motion/Fig_Supp_MotionRobustness.png";

% -------- Load data --------
T = readtable(inCSV);

% Ensure correct ordering of bands
bandOrder = ["delta","theta","alpha","beta","low_gamma","high_gamma"];
bandLabels = ["\delta","\theta","\alpha","\beta","low \gamma","high \gamma"];

% Keep only dynamic term (already isolated in table)
T.band = categorical(T.band, bandOrder, "Ordinal", true);

% Split by motion group
Tlow  = T(T.motion_group=="low", :);
Tmod  = T(T.motion_group=="moderate", :);

% -------- Prepare figure --------
figure("Color","w","Position",[100 100 850 420]);
hold on;

x = 1:numel(bandOrder);
offset = 0.12;

% ---- LOW MOTION ----
y_low  = Tlow.beta_dev;
errL_low = y_low - Tlow.ci_low_dev;
errU_low = Tlow.ci_high_dev - y_low;

errorbar(x-offset, y_low, errL_low, errU_low, ...
    "o-", "LineWidth",1.6, "MarkerSize",7, ...
    "CapSize",10, "DisplayName","Low motion");

% ---- MODERATE MOTION ----
y_mod  = Tmod.beta_dev;
errL_mod = y_mod - Tmod.ci_low_dev;
errU_mod = Tmod.ci_high_dev - y_mod;

errorbar(x+offset, y_mod, errL_mod, errU_mod, ...
    "s--", "LineWidth",1.6, "MarkerSize",7, ...
    "CapSize",10, "DisplayName","Moderate motion");

% -------- Axes & labels --------
set(gca,"XTick",x,"XTickLabel",bandLabels,"FontSize",11);
ylabel("\beta (\Delta EEG FC)","FontWeight","bold");
xlabel("EEG frequency band","FontWeight","bold");

yline(0,"k:","LineWidth",1);

legend(["Low motion","Moderate motion"], ...
       "Location","northwest","Box","off");


title("Motion robustness of dynamic EEG–fMRI coupling","FontWeight","bold");

xlim([0.5 numel(x)+0.5]);
box off;
grid on;

% -------- Save --------
exportgraphics(gcf, outFig, "Resolution",300);
close(gcf);

fprintf("[OK] Saved motion robustness figure:\n%s\n", outFig);
