clc;
clear;
close all;

%% Select folder (the Zemax folder that was used in EFL_Batch_Analysis.m)
folderPath = uigetdir(pwd, 'Select Zemax folder');

if isequal(folderPath, 0)
    return;
end

% EFL results are in the "EFL ANALİZİ" subfolder (char(304) = 'İ')
eflFolder = fullfile(folderPath, ['EFL ANAL' char(304) 'Z' char(304)]);
eflFile = fullfile(eflFolder, 'EFL_Results.xlsx');

if ~isfile(eflFile)
    error('EFL_Results.xlsx not found in:\n%s\nRun EFL_Batch_Analysis.m first.', eflFolder);
end

%% Read EFL results
try
    EFLTable = readtable(eflFile, 'Sheet', 'EFL by File');
catch
    EFLTable = readtable(eflFile);
end

fileNames = string(EFLTable.FileName);
EFL_mm = EFLTable.EFL_mm;
status = string(EFLTable.Status);

valid = isfinite(EFL_mm) & status == "OK" & EFL_mm ~= 0;
fileNames = fileNames(valid);
EFL_mm = EFL_mm(valid);

if isempty(EFL_mm)
    error('No valid EFL values found in EFL_Results.xlsx.');
end

%% Ask user for detector information
prompt = {'Detector format - horizontal pixel count (e.g. 640):', ...
          'Detector format - vertical pixel count (e.g. 512):', ...
          'Pixel pitch (um):'};
answer = inputdlg(prompt, 'Detector Information', [1 50], {'640', '512', '12'});

if isempty(answer)
    return;
end

nPixH = str2double(answer{1});
nPixV = str2double(answer{2});
pitch_um = str2double(answer{3});

if ~(isfinite(nPixH) && nPixH > 0 && nPixH == round(nPixH))
    error('Horizontal pixel count must be a positive integer.');
end
if ~(isfinite(nPixV) && nPixV > 0 && nPixV == round(nPixV))
    error('Vertical pixel count must be a positive integer.');
end
if ~(isfinite(pitch_um) && pitch_um > 0)
    error('Pixel pitch must be a positive number.');
end

%% Detector size
detW_mm = nPixH * pitch_um / 1000;
detH_mm = nPixV * pitch_um / 1000;
detD_mm = hypot(detW_mm, detH_mm);

fprintf('Detector: %d x %d px, pitch %.3f um\n', nPixH, nPixV, pitch_um);
fprintf('Detector size: %.4f mm x %.4f mm (diagonal %.4f mm)\n', ...
    detW_mm, detH_mm, detD_mm);

%% Field of view for each EFL
% FOV = 2 * atan(detector size / (2 * EFL))
HFOV_deg = 2 * atand(detW_mm ./ (2 * abs(EFL_mm)));
VFOV_deg = 2 * atand(detH_mm ./ (2 * abs(EFL_mm)));

%% Output folder: "FOV ANALİZİ" (char(304) = 'İ', avoids .m file encoding issues)
outFolder = fullfile(folderPath, ['FOV ANAL' char(304) 'Z' char(304)]);

if ~isfolder(outFolder)
    mkdir(outFolder);
end

%% Statistical summary
statNames = ["Count"; "Mean"; "Median"; "Std"; "Min"; "Max"];
calcStats = @(x) [numel(x); mean(x); median(x); std(x); min(x); max(x)];
hStats = calcStats(HFOV_deg);
vStats = calcStats(VFOV_deg);

fprintf('\nNumber of valid files: %d\n', numel(EFL_mm));
fprintf('\n%-8s %14s %14s\n', '', 'HFOV (deg)', 'VFOV (deg)');
for i = 2:numel(statNames)
    fprintf('%-8s %14.4f %14.4f\n', statNames(i), hStats(i), vStats(i));
end

%% Histograms: HFOV and VFOV distributions
fig = figure('Name', 'FOV Distribution', 'Color', 'w', 'Position', [100 100 1100 450]);

subplot(1, 2, 1);
histogram(HFOV_deg, 'BinMethod', 'auto');
xlabel('Horizontal FOV (deg)');
ylabel('Number of optical systems');
title(sprintf('Horizontal FOV - %d files', numel(HFOV_deg)));
grid on;

subplot(1, 2, 2);
histogram(VFOV_deg, 'BinMethod', 'auto');
xlabel('Vertical FOV (deg)');
ylabel('Number of optical systems');
title(sprintf('Vertical FOV - %d files', numel(VFOV_deg)));
grid on;

sgtitle(sprintf('Detector %d x %d px, %.2f um pitch', nPixH, nPixV, pitch_um));

imageFile = fullfile(outFolder, 'FOV_Histogram.png');
try
    exportgraphics(fig, imageFile, 'Resolution', 300);   % R2020a+
catch
    saveas(fig, imageFile);
end

%% Export Excel
outputFile = fullfile(outFolder, 'FOV_Results.xlsx');

% Remove old file so stale rows/sheets are not kept
if isfile(outputFile)
    delete(outputFile);
end

Results = table(fileNames, EFL_mm, HFOV_deg, VFOV_deg, ...
    'VariableNames', {'FileName', 'EFL_mm', 'HFOV_deg', 'VFOV_deg'});
writetable(Results, outputFile, 'Sheet', 'FOV by File');

Stats = table(statNames, hStats, vStats, ...
    'VariableNames', {'Statistic', 'HFOV_deg', 'VFOV_deg'});
writetable(Stats, outputFile, 'Sheet', 'Statistics');

Detector = table( ...
    ["Horizontal pixels"; "Vertical pixels"; "Pixel pitch (um)"; ...
     "Detector width (mm)"; "Detector height (mm)"; "Detector diagonal (mm)"], ...
    [nPixH; nPixV; pitch_um; detW_mm; detH_mm; detD_mm], ...
    'VariableNames', {'Parameter', 'Value'});
writetable(Detector, outputFile, 'Sheet', 'Detector');

fprintf('\nHistogram saved to:\n%s\n', imageFile);
fprintf('Results saved to:\n%s\n', outputFile);
