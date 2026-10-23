clc;
clear;
close all;

%% User settings
filePrefix  = 'MC_T';             % only files starting with this are used
maxFreq     = 41.7;               % maximum spatial frequency (cyc/mm), e.g. Nyquist
nFreqPoints = 100;                % number of frequency steps from 0 to maxFreq
yieldLevels = [50 80 90 95];      % yield levels (%)
% FFT MTF sample size is set in the "Compute MTF" section (S_64x64 by default)

%% Select folder
folderPath = uigetdir(pwd, 'Select Zemax folder');

if isequal(folderPath, 0)
    return;
end

% Find ZMX files whose names start with the prefix
files = dir(fullfile(folderPath, [filePrefix '*.zmx']));
filesUpper = dir(fullfile(folderPath, [filePrefix '*.ZMX']));
files = [files; filesUpper];

% Windows dir() is case-insensitive; keep only names starting exactly with the prefix
files = files(startsWith(string({files.name}), filePrefix));

if isempty(files)
    error('No ZMX files starting with "%s" found.', filePrefix);
end

% Remove duplicate filenames, if any
[~, uniqueIdx] = unique(lower(string({files.name})), 'stable');
files = files(uniqueIdx);
nFiles = numel(files);

%% Initialize ZOS-API connection
% InitConnection and CleanupConnection are from the
% official Zemax MATLAB standalone sample files.

TheApplication = InitConnection();

if isempty(TheApplication)
    error('Could not connect to OpticStudio.');
end

%% Compute MTF (polychromatic FFT MTF, all fields, tangential and sagittal)
freqGrid = linspace(0, maxFreq, nFreqPoints + 1)';
nFreq = numel(freqGrid);

fileNames = string({files.name})';
status = strings(nFiles, 1);
MTF_T = [];          % nFiles x nFreq x nFields (tangential)
MTF_S = [];          % nFiles x nFreq x nFields (sagittal)
fieldLabels = strings(0, 1);
nFields = 0;
onAxisField = 1;
diffT = [];          % diffraction limit (tangential)
diffS = [];          % diffraction limit (sagittal)
diffSource = "";

fprintf('Computing MTF for %d files...\n', nFiles);
wb = waitbar(0, 'Computing MTF...', 'Name', 'MTF Analysis');

try
    TheSystem = TheApplication.PrimarySystem;

    for k = 1:nFiles

        filePath = fullfile(folderPath, files(k).name);
        mtf = [];

        try
            % Load optical system
            loaded = TheSystem.LoadFile(filePath, false);

            if ~loaded
                error('OpticStudio could not load this file.');
            end

            % FFT MTF analysis
            mtf = TheSystem.Analyses.New_FftMtf();
            settings = mtf.GetSettings();
            settings.MaximumFrequency = maxFreq;
            settings.SampleSize = ZOSAPI.Analysis.SampleSizes.S_64x64;
            settings.Wavelength.UseAllWavelengths();
            settings.Field.UseAllFields();
            try
                settings.ShowDiffractionLimit = true;
            catch
            end

            mtf.ApplyAndWaitForCompletion();
            results = mtf.GetResults();
            nSeries = results.NumberOfDataSeries;
            dataSeries = results.DataSeries;

            % Separate field series from the diffraction limit series
            isDiff = false(nSeries, 1);
            for s = 1:nSeries
                isDiff(s) = contains(lower(string(char(dataSeries(s).Description))), "diff");
            end
            fieldIdx = find(~isDiff);

            % First successful file defines the number of fields
            if nFields == 0
                nFields = numel(fieldIdx);
                MTF_T = nan(nFiles, nFreq, nFields);
                MTF_S = nan(nFiles, nFreq, nFields);
                fieldLabels = strings(nFields, 1);
                onAxisField = findOnAxisField(TheSystem);
            elseif numel(fieldIdx) ~= nFields
                error('Number of fields (%d) differs from first file (%d).', numel(fieldIdx), nFields);
            end

            for f = 1:nFields
                data = dataSeries(fieldIdx(f));
                [MTF_T(k, :, f), MTF_S(k, :, f)] = readSeries(data, freqGrid);

                if fieldLabels(f) == ""
                    fieldLabels(f) = strtrim(string(char(data.Description)));
                end
            end

            % Diffraction limit (taken once, from the first successful file)
            if isempty(diffT)
                diffIdx = find(isDiff, 1);
                if ~isempty(diffIdx)
                    [diffT, diffS] = readSeries(dataSeries(diffIdx), freqGrid);
                    diffSource = "OpticStudio FFT MTF";
                else
                    diffT = analyticDiffLimit(TheSystem, freqGrid);
                    diffS = diffT;
                    diffSource = "Analytical (circular aperture, working F/#)";
                end
            end

            if any(isnan(MTF_T(k, :, :)), 'all') || any(isnan(MTF_S(k, :, :)), 'all')
                error('MTF data contains invalid values.');
            end

            status(k) = "OK";

        catch ME
            status(k) = string(ME.message);
            fprintf('ERROR: %s - %s\n', files(k).name, ME.message);
        end

        if ~isempty(mtf)
            mtf.Close();
        end

        if isvalid(wb)
            waitbar(k / nFiles, wb, sprintf('Computing MTF... %d / %d', k, nFiles));
        end
    end

catch ME
    if isvalid(wb), close(wb); end
    CleanupConnection(TheApplication);
    rethrow(ME);
end

if isvalid(wb), close(wb); end
CleanupConnection(TheApplication);

%% Keep valid results only
valid = status == "OK";
nValid = sum(valid);

if nValid == 0
    error('No valid MTF results were obtained.');
end

MTF_T_valid = MTF_T(valid, :, :);
MTF_S_valid = MTF_S(valid, :, :);

% Clip small interpolation overshoots
MTF_T_valid = min(max(MTF_T_valid, 0), 1);
MTF_S_valid = min(max(MTF_S_valid, 0), 1);
diffT = min(max(diffT, 0), 1);
diffS = min(max(diffS, 0), 1);

for f = 1:nFields
    if fieldLabels(f) == ""
        fieldLabels(f) = sprintf('Field %d', f);
    end
end

fprintf('Valid files: %d / %d\n', nValid, nFiles);

%% Mean MTF and yield MTF curves
% Yield Y% MTF = the MTF value that at least Y% of the systems meet or exceed
nYield = numel(yieldLevels);
meanT = zeros(nFreq, nFields);
meanS = zeros(nFreq, nFields);
yieldT = zeros(nFreq, nYield, nFields);
yieldS = zeros(nFreq, nYield, nFields);

for f = 1:nFields
    meanT(:, f) = mean(MTF_T_valid(:, :, f), 1)';
    meanS(:, f) = mean(MTF_S_valid(:, :, f), 1)';
    for y = 1:nYield
        yieldT(:, y, f) = yieldValue(MTF_T_valid(:, :, f), yieldLevels(y))';
        yieldS(:, y, f) = yieldValue(MTF_S_valid(:, :, f), yieldLevels(y))';
    end
end

%% Output folder: "MTF ANALİZİ" (char(304) = 'İ', avoids .m file encoding issues)
outFolder = fullfile(folderPath, ['MTF ANAL' char(304) 'Z' char(304)]);

if ~isfolder(outFolder)
    mkdir(outFolder);
end

%% Plots: mean MTF + one figure per yield level (left: all fields, right: on-axis field)
% Solid: tangential, dashed: sagittal, black: diffraction limit
fieldColors = lines(nFields);
showDiffS = max(abs(diffT - diffS)) > 1e-6;

% Level 1: mean MTF, levels 2..end: yield MTF
plotT = cat(2, reshape(meanT, nFreq, 1, nFields), yieldT);
plotS = cat(2, reshape(meanS, nFreq, 1, nFields), yieldS);
levelNames = ["Mean MTF", compose("%d%% yield", yieldLevels)];
levelFiles = ["MTF_Mean", compose("MTF_Yield%d", yieldLevels)];
nLevels = numel(levelNames);
levelFigs = gobjects(nLevels, 1);

for y = 1:nLevels
    fig = figure('Name', sprintf('MTF - %s', levelNames(y)), 'Color', 'w', ...
        'Position', [100 100 1400 550]);

    panels = {1:nFields, onAxisField};
    panelTitles = ["All fields", "On-axis field: " + fieldLabels(onAxisField)];

    for p = 1:2
        subplot(1, 2, p);
        hold on;
        plot(freqGrid, diffT, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Diff. limit');
        if showDiffS
            plot(freqGrid, diffS, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Diff. limit S');
        end
        for f = panels{p}
            plot(freqGrid, plotT(:, y, f), '-', 'Color', fieldColors(f, :), ...
                'LineWidth', 1.2, 'DisplayName', sprintf('F%d T', f));
            plot(freqGrid, plotS(:, y, f), '--', 'Color', fieldColors(f, :), ...
                'LineWidth', 1.2, 'DisplayName', sprintf('F%d S', f));
        end
        hold off;
        xlim([0 maxFreq]);
        ylim([0 1]);
        xlabel('Spatial frequency (cyc/mm)');
        ylabel('Modulus of the OTF');
        title(panelTitles(p), 'Interpreter', 'none');
        legend('Location', 'northeast', 'NumColumns', 2);
        grid on;
        box on;
    end

    sgtitle(sprintf('Polychromatic FFT MTF - %s (%d files)', levelNames(y), nValid));
    levelFigs(y) = fig;

    imageFile = fullfile(outFolder, levelFiles(y) + ".png");
    try
        exportgraphics(fig, imageFile, 'Resolution', 300);   % R2020a+
    catch
        saveas(fig, imageFile);
    end
end

%% Summary at maximum frequency
yieldNames = compose("Yield%d", yieldLevels);
Summary = table(fieldLabels, meanT(end, :)', meanS(end, :)', ...
    'VariableNames', {'Field', 'Mean_T', 'Mean_S'});
for y = 1:nYield
    Summary.(char(yieldNames(y) + "_T")) = squeeze(yieldT(end, y, :));
    Summary.(char(yieldNames(y) + "_S")) = squeeze(yieldS(end, y, :));
end

%% Export Excel
outputFile = fullfile(outFolder, 'MTF_Results.xlsx');

% Remove old file so stale rows/sheets are not kept
if isfile(outputFile)
    delete(outputFile);
end

writetable(Summary, outputFile, 'Sheet', 'Summary at max freq');

for f = 1:nFields
    FieldTable = table(freqGrid, meanT(:, f), meanS(:, f), ...
        'VariableNames', {'Frequency_cyc_mm', 'Mean_T', 'Mean_S'});
    for y = 1:nYield
        FieldTable.(char(yieldNames(y) + "_T")) = yieldT(:, y, f);
        FieldTable.(char(yieldNames(y) + "_S")) = yieldS(:, y, f);
    end
    writetable(FieldTable, outputFile, 'Sheet', sprintf('Field %d', f));
end

FieldInfo = table((1:nFields)', fieldLabels, 'VariableNames', {'Field', 'Description'});
writetable(FieldInfo, outputFile, 'Sheet', 'Field Info');

DiffTable = table(freqGrid, diffT, diffS, ...
    'VariableNames', {'Frequency_cyc_mm', 'DiffLimit_T', 'DiffLimit_S'});
writetable(DiffTable, outputFile, 'Sheet', 'Diffraction Limit');

Files = table(fileNames, status, 'VariableNames', {'FileName', 'Status'});
writetable(Files, outputFile, 'Sheet', 'Files');

% Raw per-file MTF curves for further analysis
save(fullfile(outFolder, 'MTF_RawData.mat'), 'freqGrid', 'MTF_T', 'MTF_S', ...
    'diffT', 'diffS', 'fileNames', 'status', 'fieldLabels', 'yieldLevels');

%% PDF report: summary page + mean MTF page + one page per yield level
pdfFile = fullfile(outFolder, 'MTF_Report.pdf');

if isfile(pdfFile)
    delete(pdfFile);
end

% Summary page text
colNames = Summary.Properties.VariableNames(2:end);
header = sprintf('%-30s', 'Field');
header = [header, sprintf('%12s', colNames{:})];

report = {
    'MTF YIELD ANALYSIS REPORT'
    ''
    sprintf('Date:              %s', char(datetime('now', 'Format', 'dd.MM.yyyy HH:mm')))
    sprintf('Folder:            %s', folderPath)
    sprintf('File prefix:       %s', filePrefix)
    sprintf('Valid files:       %d / %d', nValid, nFiles)
    sprintf('Max frequency:     %.2f cyc/mm', maxFreq)
    sprintf('Yield levels:      %s %%', strjoin(string(yieldLevels), ', '))
    'MTF type:          Polychromatic FFT MTF (all wavelengths, all fields)'
    sprintf('On-axis field:     Field %d (%s)', onAxisField, char(fieldLabels(onAxisField)))
    sprintf('Diffraction limit: %s', char(diffSource))
    ''
    'Fields (plots: solid = tangential, dashed = sagittal)'
    ''
    sprintf('MTF at %.2f cyc/mm (T: tangential, S: sagittal)', maxFreq)
    sprintf('Yield X%% = MTF value that at least X%% of the systems meet or exceed')
    ''
    header
    repmat('-', 1, numel(header))
    };

fieldLines = compose("  F%d: %s", (1:nFields)', fieldLabels);
report = [report(1:14); cellstr(fieldLines); report(15:end)];

for f = 1:nFields
    label = sprintf('F%d: %s', f, char(fieldLabels(f)));
    label = label(1:min(end, 28));
    values = table2array(Summary(f, 2:end));
    report{end + 1, 1} = [sprintf('%-30s', label), sprintf('%12.4f', values)]; %#ok<SAGROW>
end
report{end + 1, 1} = sprintf('%-30s%12.4f%12.4f', 'Diffraction limit', diffT(end), diffS(end));

failed = fileNames(~valid);
if ~isempty(failed)
    report = [report; {''; sprintf('Failed files (%d):', numel(failed))}];
    shown = failed(1:min(end, 15));
    report = [report; cellstr("  " + shown)];
    if numel(failed) > numel(shown)
        report{end + 1, 1} = sprintf('  ... and %d more (see Files sheet in Excel)', numel(failed) - numel(shown));
    end
end

reportFig = figure('Name', 'MTF Report', 'Color', 'w', 'Position', [100 100 1200 850]);
ax = axes(reportFig, 'Position', [0 0 1 1], 'Visible', 'off');
text(ax, 0.04, 0.96, report, 'VerticalAlignment', 'top', ...
    'FontName', 'Courier New', 'FontSize', 9, 'Interpreter', 'none');

allFigs = [reportFig; levelFigs];

try
    % Multi-page PDF (exportgraphics with 'Append' needs R2021b+)
    exportgraphics(allFigs(1), pdfFile, 'ContentType', 'vector');
    for p = 2:numel(allFigs)
        exportgraphics(allFigs(p), pdfFile, 'ContentType', 'vector', 'Append', true);
    end
    fprintf('\nPDF report saved to:\n%s\n', pdfFile);
catch
    % Older MATLAB: one PDF file per page
    if isfile(pdfFile)
        delete(pdfFile);
    end
    for p = 1:numel(allFigs)
        set(allFigs(p), 'PaperOrientation', 'landscape', 'PaperPositionMode', 'auto');
        print(allFigs(p), fullfile(outFolder, sprintf('MTF_Report_page%d.pdf', p)), '-dpdf', '-bestfit');
    end
    fprintf('\nMulti-page PDF not supported in this MATLAB version;\n');
    fprintf('report saved as separate pages (MTF_Report_page*.pdf).\n');
end

fprintf('\nResults saved to:\n%s\n', outFolder);


function v = yieldValue(data, yieldPct)
% data: nSystems x nFreq. Returns, for each frequency, the MTF value that
% at least yieldPct % of the systems meet or exceed (no toolbox needed).
sorted = sort(data, 1, 'descend');
idx = max(1, ceil(yieldPct / 100 * size(data, 1)));
v = sorted(idx, :);
end

function [t, s] = readSeries(data, freqGrid)
% Reads one FFT MTF data series (column 1: tangential, column 2: sagittal)
% and interpolates it to the common frequency grid.
x = data.XData.Data.double;
y = data.YData.Data.double;

x = x(:);
if size(y, 1) ~= numel(x)
    y = y.';
end

t = interp1(x, y(:, 1), freqGrid, 'linear', 'extrap');
if size(y, 2) > 1
    s = interp1(x, y(:, 2), freqGrid, 'linear', 'extrap');
else
    s = t;
end
end

function idx = findOnAxisField(TheSystem)
% Returns the first field at (0, 0); field 1 if none is found
idx = 1;
try
    sysFields = TheSystem.SystemData.Fields;
    for i = 1:sysFields.NumberOfFields
        fld = sysFields.GetField(i);
        if fld.X == 0 && fld.Y == 0
            idx = i;
            return;
        end
    end
catch
end
end

function mtf = analyticDiffLimit(TheSystem, freqGrid)
% Polychromatic diffraction-limited MTF of a circular aperture (on axis),
% used only if OpticStudio does not return a diffraction limit series.
fno = TheSystem.MFE.GetOperandValue( ...
    ZOSAPI.Editors.MFE.MeritOperandType.WFNO, 0, 0, 0, 0, 0, 0, 0, 0);
wls = TheSystem.SystemData.Wavelengths;

mtf = zeros(size(freqGrid));
wSum = 0;
for i = 1:wls.NumberOfWavelengths
    wl = wls.GetWavelength(i);
    fc = 1 / (wl.Wavelength * 1e-3 * abs(fno));   % cutoff frequency (cyc/mm)
    nu = min(freqGrid / fc, 1);
    mtf = mtf + wl.Weight * (2 / pi) * (acos(nu) - nu .* sqrt(1 - nu.^2));
    wSum = wSum + wl.Weight;
end
mtf = mtf / wSum;
end

function app = InitConnection()

import System.Reflection.*;

% Find the installed version of OpticStudio.
zemaxData = winqueryreg('HKEY_CURRENT_USER', 'Software\Zemax', 'ZemaxRoot');
NetHelper = strcat(zemaxData, '\ZOS-API\Libraries\ZOSAPI_NetHelper.dll');
% Note -- uncomment the following line to use a custom NetHelper path
% NetHelper = 'C:\Users\gayeisik\Documents\Zemax\ZOS-API\Libraries\ZOSAPI_NetHelper.dll';

% This is the path to OpticStudio
NET.addAssembly(NetHelper);

success = ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize();
% Note -- uncomment the following line to use a custom initialization path
% success = ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize('C:\Program Files\Ansys Zemax OpticStudio 2025 R2.04');
if success == 1
    LogMessage(strcat('Found OpticStudio at: ', char(ZOSAPI_NetHelper.ZOSAPI_Initializer.GetZemaxDirectory())));
else
    LogMessage('NOT-Found OpticStudio');
    app = [];
    return;
end

% Now load the ZOS-API assemblies
NET.addAssembly(AssemblyName('ZOSAPI_Interfaces'));
NET.addAssembly(AssemblyName('ZOSAPI'));

% Create the initial connection class
TheConnection = ZOSAPI.ZOSAPI_Connection();

% NOTE - if this fails with a message like 'Unable to load one or more of
% the requested types', it is usually caused by try to connect to a 32-bit
% version of OpticStudio from a 64-bit version of MATLAB (or vice-versa).
TheConnection.ConnectionTimeoutSeconds = 20;
app = TheConnection.CreateNewApplication();

if isempty(app)
    HandleError('An unknown connection error occurred!');
end
if ~app.IsValidLicenseForAPI
    HandleError('License check failed!');
end

end

function LogMessage(msg)
disp(msg);
end

function HandleError(error)
ME = MException('zosapi:HandleError', error);
throw(ME);
end

function CleanupConnection(TheApplication)
% Note - this will close down the connection.
TheApplication.CloseApplication();
end
