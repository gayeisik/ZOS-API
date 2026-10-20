clc;
clear;
close all;

%% Select folder (the Zemax folder that was used in EFL_Batch_Analysis.m)
folderPath = uigetdir(pwd, 'Select Zemax folder');

if isequal(folderPath, 0)
    return;
end

% FOV results are in the "FOV ANALİZİ" subfolder (char(304) = 'İ')
fovFolder = fullfile(folderPath, ['FOV ANAL' char(304) 'Z' char(304)]);
fovFile = fullfile(fovFolder, 'FOV_Results.xlsx');

if ~isfile(fovFile)
    error('FOV_Results.xlsx not found in:\n%s\nRun FOV_Analysis.m first.', fovFolder);
end

%% Read EFL and FOV results
try
    Data = readtable(fovFile, 'Sheet', 'FOV by File');
catch
    Data = readtable(fovFile);
end

fileNames = string(Data.FileName);
EFL_mm = Data.EFL_mm;
HFOV_deg = Data.HFOV_deg;
VFOV_deg = Data.VFOV_deg;

if isempty(EFL_mm)
    error('No data found in FOV_Results.xlsx.');
end

%% Ask user for target values and tolerances
prompt = {'Target EFL (mm):', ...
          'Target horizontal FOV (deg):', ...
          'Target vertical FOV (deg):', ...
          'EFL tolerance (+/- %):', ...
          'Horizontal FOV tolerance (+/- %):', ...
          'Vertical FOV tolerance (+/- %):'};
answer = inputdlg(prompt, 'Target Values', [1 50], {'', '', '', '2', '2', '2'});

if isempty(answer)
    return;
end

inputs = str2double(answer);

if any(~isfinite(inputs))
    error('All inputs must be numbers.');
end

targetEFL  = inputs(1);
targetHFOV = inputs(2);
targetVFOV = inputs(3);
tolEFL     = inputs(4);
tolHFOV    = inputs(5);
tolVFOV    = inputs(6);

if targetEFL == 0
    error('Target EFL cannot be zero.');
end
if targetHFOV <= 0 || targetVFOV <= 0
    error('Target FOV values must be positive.');
end
if any(inputs(4:6) < 0)
    error('Tolerances must be non-negative.');
end

%% Percentage error for each file
% error (%) = (actual - target) / target * 100   (signed: + larger, - smaller)
pctError = @(x, t) (x - t) ./ abs(t) * 100;

EFL_err_pct  = pctError(EFL_mm, targetEFL);
HFOV_err_pct = pctError(HFOV_deg, targetHFOV);
VFOV_err_pct = pctError(VFOV_deg, targetVFOV);

% Absolute FOV error in degrees (for reference)
HFOV_err_deg = HFOV_deg - targetHFOV;
VFOV_err_deg = VFOV_deg - targetVFOV;

% Tolerance check
EFL_pass  = abs(EFL_err_pct)  <= tolEFL;
HFOV_pass = abs(HFOV_err_pct) <= tolHFOV;
VFOV_pass = abs(VFOV_err_pct) <= tolVFOV;
All_pass  = EFL_pass & HFOV_pass & VFOV_pass;

passText = @(p) replace(string(p), ["true", "false"], ["PASS", "FAIL"]);

%% Statistics of the errors
statNames = ["Count"; "Mean error (%)"; "Median error (%)"; "Std error (%)"; ...
             "Min error (%)"; "Max error (%)"; "Max |error| (%)"; ...
             "95th percentile |error| (%)"; "Tolerance (+/- %)"; "Yield (%)"];

calcStats = @(e, tol, pass) [numel(e); mean(e); median(e); std(e); ...
                             min(e); max(e); max(abs(e)); ...
                             percentile95(abs(e)); tol; 100 * mean(pass)];

eflStats  = calcStats(EFL_err_pct,  tolEFL,  EFL_pass);
hfovStats = calcStats(HFOV_err_pct, tolHFOV, HFOV_pass);
vfovStats = calcStats(VFOV_err_pct, tolVFOV, VFOV_pass);
allYield  = 100 * mean(All_pass);

fprintf('Number of files: %d\n\n', numel(EFL_mm));
fprintf('%-28s %12s %12s %12s\n', '', 'EFL', 'HFOV', 'VFOV');
for i = 2:numel(statNames)
    fprintf('%-28s %12.4f %12.4f %12.4f\n', statNames(i), ...
        eflStats(i), hfovStats(i), vfovStats(i));
end
fprintf('\nOverall yield (EFL, HFOV and VFOV all in tolerance): %.2f %%\n', allYield);

%% Output folder: "HATA ANALİZİ"
outFolder = fullfile(folderPath, ['HATA ANAL' char(304) 'Z' char(304)]);

if ~isfolder(outFolder)
    mkdir(outFolder);
end

%% Histograms: percentage error distributions with tolerance limits
fig = figure('Name', 'Error Distribution', 'Color', 'w', 'Position', [100 100 1500 450]);

errData = {EFL_err_pct, HFOV_err_pct, VFOV_err_pct};
tolData = [tolEFL, tolHFOV, tolVFOV];
labels  = ["EFL", "Horizontal FOV", "Vertical FOV"];
yields  = [eflStats(end), hfovStats(end), vfovStats(end)];

for i = 1:3
    subplot(1, 3, i);
    histogram(errData{i}, 'BinMethod', 'auto');
    hold on;
    xline(-tolData(i), 'r--', 'LineWidth', 1.5);
    xline( tolData(i), 'r--', 'LineWidth', 1.5);
    xline(0, 'k-');
    hold off;
    xlabel(sprintf('%s error (%%)', labels(i)));
    ylabel('Number of optical systems');
    title(sprintf('%s - yield %.1f %%', labels(i), yields(i)));
    grid on;
end

sgtitle(sprintf('Percentage error distribution - %d files (red: tolerance limits)', numel(EFL_mm)));

imageFile = fullfile(outFolder, 'Error_Histogram.png');
try
    exportgraphics(fig, imageFile, 'Resolution', 300);   % R2020a+
catch
    saveas(fig, imageFile);
end

%% Export Excel
outputFile = fullfile(outFolder, 'Error_Results.xlsx');

% Remove old file so stale rows/sheets are not kept
if isfile(outputFile)
    delete(outputFile);
end

Results = table(fileNames, ...
    EFL_mm, EFL_err_pct, passText(EFL_pass), ...
    HFOV_deg, HFOV_err_pct, HFOV_err_deg, passText(HFOV_pass), ...
    VFOV_deg, VFOV_err_pct, VFOV_err_deg, passText(VFOV_pass), ...
    passText(All_pass), ...
    'VariableNames', {'FileName', ...
    'EFL_mm', 'EFL_err_pct', 'EFL_result', ...
    'HFOV_deg', 'HFOV_err_pct', 'HFOV_err_deg', 'HFOV_result', ...
    'VFOV_deg', 'VFOV_err_pct', 'VFOV_err_deg', 'VFOV_result', ...
    'Overall_result'});
writetable(Results, outputFile, 'Sheet', 'Errors by File');

Stats = table(statNames, eflStats, hfovStats, vfovStats, ...
    'VariableNames', {'Statistic', 'EFL', 'HFOV', 'VFOV'});
writetable(Stats, outputFile, 'Sheet', 'Statistics');

Inputs = table( ...
    ["Target EFL (mm)"; "Target HFOV (deg)"; "Target VFOV (deg)"; ...
     "EFL tolerance (+/- %)"; "HFOV tolerance (+/- %)"; "VFOV tolerance (+/- %)"; ...
     "Overall yield (%)"], ...
    [inputs; allYield], ...
    'VariableNames', {'Parameter', 'Value'});
writetable(Inputs, outputFile, 'Sheet', 'Inputs');

fprintf('\nHistogram saved to:\n%s\n', imageFile);
fprintf('Results saved to:\n%s\n', outputFile);


function p = percentile95(x)
% 95th percentile (nearest-rank), no Statistics Toolbox needed
x = sort(x(:));
p = x(max(1, ceil(0.95 * numel(x))));
end
