clc;
clear;
close all;

%% User settings
filePrefix    = 'MC_T';             % only files starting with this are used
pixelPitch_um = 15;                 % detector pixel pitch (um)
nFreqPoints   = 100;                % number of frequency steps from 0 to maxFreq
yieldLevels   = [50 80 90 95];      % yield levels (%)
reportYield   = 90;                 % yield level shown in the per-error MTF plots
compErrors_um = [10 20 30 40 50];   % compensator setting errors (+/- um)
compSurface   = [];                 % compensator surface; [] = read from COMP operand of nominal file
nFilesMax     = Inf;                % use only the first N MC files (e.g. 50 for a quick test)
% FFT MTF sample size is set in the "Compute MTF" section (S_64x64 by default)
% Nominal file: the .zmx file in the folder whose name does NOT start with "MC_"
%
% For every MC file and every error d, the compensator thickness is set to
% (compensated value + d) and (compensated value - d); the worse of the two
% (lower mean MTF over all fields, T and S, all frequencies) is kept.

% Nyquist frequency of the detector, used as the maximum MTF frequency
maxFreq = 1000 / (2 * pixelPitch_um);   % cyc/mm

if ~ismember(reportYield, yieldLevels)
    error('reportYield (%g) must be one of yieldLevels.', reportYield);
end

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
files = files(1:min(end, nFilesMax));
nFiles = numel(files);

% Nominal file: .zmx file whose name does not start with "MC_"
allZmx = [dir(fullfile(folderPath, '*.zmx')); dir(fullfile(folderPath, '*.ZMX'))];
[~, uniqueIdx] = unique(lower(string({allZmx.name})), 'stable');
allZmx = allZmx(uniqueIdx);
nominalFiles = allZmx(~startsWith(string({allZmx.name}), 'MC_', 'IgnoreCase', true));

if isempty(nominalFiles)
    error('No nominal .zmx file (name not starting with "MC_") found in the folder.');
elseif isscalar(nominalFiles)
    nominalName = nominalFiles(1).name;
else
    [sel, ok] = listdlg('ListString', {nominalFiles.name}, 'SelectionMode', 'single', ...
        'PromptString', 'Select the nominal file:', 'ListSize', [400 200]);
    if ~ok
        return;
    end
    nominalName = nominalFiles(sel).name;
end

%% Initialize ZOS-API connection
% InitConnection and CleanupConnection are from the
% official Zemax MATLAB standalone sample files.

TheApplication = InitConnection();

if isempty(TheApplication)
    error('Could not connect to OpticStudio.');
end

%% Compute MTF for each MC file and each compensator error
freqGrid = linspace(0, maxFreq, nFreqPoints + 1)';
nFreq = numel(freqGrid);

nErr = numel(compErrors_um);
levelErrors_um = [0, compErrors_um];   % level 1: no error (compensated file as is)
nLevels = nErr + 1;

fileNames = string({files.name})';
status = strings(nFiles, 1);
MTF_T = [];          % nFiles x nFreq x nFields x nLevels (tangential)
MTF_S = [];          % nFiles x nFreq x nFields x nLevels (sagittal)
worstSign = nan(nFiles, nErr);   % +1: +d was worse, -1: -d was worse
compThickness = nan(nFiles, 1);  % compensated thickness of each MC file
fieldLabels = strings(0, 1);
nFields = 0;
diffT = [];
diffS = [];
diffSource = "";
paraxialEFL = NaN;
freqScale = 1;       % Zemax frequency -> cyc/mm (1000/|EFL| for afocal files)

fprintf('Computing MTF for %d files x %d compensator errors (+/-)...\n', nFiles, nErr);
wb = waitbar(0, 'Computing MTF...', 'Name', 'Compensator Error Analysis');

try
    TheSystem = TheApplication.PrimarySystem;

    %% Nominal file: afocal check, EFL, compensator, on-axis field
    if ~TheSystem.LoadFile(fullfile(folderPath, nominalName), false)
        error('OpticStudio could not load the nominal file: %s', nominalName);
    end

    nominalAfocal = logical(TheSystem.SystemData.Aperture.AFocalImageSpace);

    if nominalAfocal
        % EFL is read in the file's own (afocal) mode
        nominalEFL = TheSystem.MFE.GetOperandValue( ...
            ZOSAPI.Editors.MFE.MeritOperandType.EFFL, 0, 0, 0, 0, 0, 0, 0, 0);

        if ~isfinite(nominalEFL) || nominalEFL == 0
            error('Could not compute the EFL of the nominal file: %s', nominalName);
        end

        paraxialEFL = abs(nominalEFL);
        freqScale = 1000 / paraxialEFL;   % afocal MTF is in cyc/mrad
        fprintf('Nominal file %s is in afocal mode, |EFL| = %.4f mm\n', nominalName, paraxialEFL);
    else
        fprintf('Nominal file %s is not in afocal mode.\n', nominalName);
    end

    if isempty(compSurface)
        compSurface = findCompensator(TheSystem);
    end
    compSurfaceComment = string(char(TheSystem.LDE.GetSurfaceAt(compSurface).Comment));

    lensUnitsPerMm = unitsPerMm(TheSystem);
    onAxisField = findOnAxisField(TheSystem);

    fprintf('Compensator: thickness of surface %d %s\n', compSurface, compSurfaceComment);
    fprintf('Nyquist frequency: %.4f cyc/mm (pixel pitch %.2f um)\n', maxFreq, pixelPitch_um);

    tStart = tic;

    for k = 1:nFiles

        filePath = fullfile(folderPath, files(k).name);
        mtf = [];
        surf = [];

        try
            % Load optical system
            if ~TheSystem.LoadFile(filePath, false)
                error('OpticStudio could not load this file.');
            end

            surf = TheSystem.LDE.GetSurfaceAt(compSurface);
            t0 = surf.Thickness;
            compThickness(k) = t0;

            % FFT MTF analysis (re-applied after each compensator change)
            mtf = TheSystem.Analyses.New_FftMtf();
            settings = mtf.GetSettings();
            settings.MaximumFrequency = maxFreq / freqScale;   % cyc/mrad if afocal
            settings.SampleSize = ZOSAPI.Analysis.SampleSizes.S_64x64;
            settings.Wavelength.UseAllWavelengths();
            settings.Field.UseAllFields();
            try
                settings.ShowDiffractionLimit = isempty(diffT);
            catch
            end

            % Level 1: compensated file, no error
            [T0, S0, labels, dT, dS] = computeMtf(mtf, freqGrid, freqScale);

            % First successful file defines the number of fields
            if nFields == 0
                nFields = size(T0, 2);
                MTF_T = nan(nFiles, nFreq, nFields, nLevels);
                MTF_S = nan(nFiles, nFreq, nFields, nLevels);
                fieldLabels = labels;
            elseif size(T0, 2) ~= nFields
                error('Number of fields (%d) differs from first file (%d).', size(T0, 2), nFields);
            end

            if isempty(diffT)
                if ~isempty(dT)
                    diffT = dT;
                    diffS = dS;
                    diffSource = "OpticStudio FFT MTF";
                else
                    diffT = analyticDiffLimit(TheSystem, freqGrid);
                    diffS = diffT;
                    diffSource = "Analytical (circular aperture, working F/#)";
                end
                try
                    settings.ShowDiffractionLimit = false;
                catch
                end
            end

            MTF_T(k, :, :, 1) = T0;
            MTF_S(k, :, :, 1) = S0;

            % Levels 2..end: +/- compensator error, keep the worse one
            for e = 1:nErr
                d = compErrors_um(e) / 1000 * lensUnitsPerMm;   % um -> lens units

                surf.Thickness = t0 + d;
                [Tp, Sp] = computeMtf(mtf, freqGrid, freqScale);

                surf.Thickness = t0 - d;
                [Tm, Sm] = computeMtf(mtf, freqGrid, freqScale);

                if mean([Tp(:); Sp(:)]) <= mean([Tm(:); Sm(:)])
                    MTF_T(k, :, :, e + 1) = Tp;
                    MTF_S(k, :, :, e + 1) = Sp;
                    worstSign(k, e) = 1;
                else
                    MTF_T(k, :, :, e + 1) = Tm;
                    MTF_S(k, :, :, e + 1) = Sm;
                    worstSign(k, e) = -1;
                end
            end

            surf.Thickness = t0;

            if any(isnan(MTF_T(k, :, :, :)), 'all') || any(isnan(MTF_S(k, :, :, :)), 'all')
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

MTF_T_valid = min(max(MTF_T(valid, :, :, :), 0), 1);
MTF_S_valid = min(max(MTF_S(valid, :, :, :), 0), 1);
diffT = min(max(diffT, 0), 1);
diffS = min(max(diffS, 0), 1);

for f = 1:nFields
    if fieldLabels(f) == ""
        fieldLabels(f) = sprintf('Field %d', f);
    end
end

fprintf('Valid files: %d / %d (computation time: %.1f s)\n', nValid, nFiles, toc(tStart));

%% Mean MTF and yield MTF curves for each compensator error level
% Yield Y% MTF = the MTF value that at least Y% of the systems meet or exceed
nYield = numel(yieldLevels);
meanT = zeros(nFreq, nFields, nLevels);
meanS = zeros(nFreq, nFields, nLevels);
yieldT = zeros(nFreq, nYield, nFields, nLevels);
yieldS = zeros(nFreq, nYield, nFields, nLevels);

for lev = 1:nLevels
    for f = 1:nFields
        dataT = MTF_T_valid(:, :, f, lev);
        dataS = MTF_S_valid(:, :, f, lev);
        meanT(:, f, lev) = mean(dataT, 1)';
        meanS(:, f, lev) = mean(dataS, 1)';
        for y = 1:nYield
            yieldT(:, y, f, lev) = yieldValue(dataT, yieldLevels(y))';
            yieldS(:, y, f, lev) = yieldValue(dataS, yieldLevels(y))';
        end
    end
end

%% Sensitivity: yield MTF vs compensator error at Nyquist/2 and Nyquist
% Worst field and direction (min over all fields, T and S)
keyFreqs = [maxFreq / 2, maxFreq];
keyNames = ["Nyquist/2", "Nyquist"];
[~, keyIdx] = min(abs(freqGrid - keyFreqs), [], 1);

sensYield = zeros(nLevels, nYield, 2);    % worst field
sensMean = zeros(nLevels, 2);             % worst field, mean MTF
sensOnAxis = zeros(nLevels, 2);           % on-axis field, report yield
ry = find(yieldLevels == reportYield, 1);

for lev = 1:nLevels
    for q = 1:2
        i = keyIdx(q);
        for y = 1:nYield
            sensYield(lev, y, q) = min([squeeze(yieldT(i, y, :, lev)); squeeze(yieldS(i, y, :, lev))]);
        end
        sensMean(lev, q) = min([meanT(i, :, lev), meanS(i, :, lev)]);
        sensOnAxis(lev, q) = min(yieldT(i, ry, onAxisField, lev), yieldS(i, ry, onAxisField, lev));
    end
end

%% Output folder: "KOMPANZATÖR ANALİZİ" (char(214) = 'Ö', char(304) = 'İ')
outFolder = fullfile(folderPath, ['KOMPANZAT' char(214) 'R ANAL' char(304) 'Z' char(304)]);

if ~isfolder(outFolder)
    mkdir(outFolder);
end

%% Plot 1: sensitivity (MTF vs compensator error)
axPos = {[0.05 0.32 0.42 0.56], [0.54 0.32 0.42 0.56]};   % fixed axes boxes
yieldColors = lines(nYield);

sensFig = figure('Name', 'Compensator Sensitivity', 'Color', 'w', 'Position', [50 50 1500 750]);

for q = 1:2
    ax = axes(sensFig, 'Position', axPos{q});
    hold on;
    for y = 1:nYield
        plot(levelErrors_um, sensYield(:, y, q), '-o', 'Color', yieldColors(y, :), ...
            'LineWidth', 1.5, 'MarkerFaceColor', yieldColors(y, :), ...
            'DisplayName', sprintf('%d%% yield', yieldLevels(y)));
    end
    plot(levelErrors_um, sensMean(:, q), 'k--s', 'LineWidth', 1.5, 'DisplayName', 'Mean');
    hold off;
    xlim([0 max(levelErrors_um)]);
    ylim([0 1]);
    xticks(levelErrors_um);
    xlabel('Compensator error (\pm\mum, worse side)');
    ylabel('MTF (worst field, T/S)');
    title(sprintf('%s (%.2f cyc/mm)', keyNames(q), freqGrid(keyIdx(q))));
    grid on;
    box on;
    pbaspect(ax, [1.6 1 1]);

    lgd = legend(ax, 'NumColumns', 5, 'Location', 'none');
    lgd.Units = 'normalized';
    lp = lgd.Position;
    lgd.Position = [axPos{q}(1) + (axPos{q}(3) - lp(3)) / 2, 0.24 - lp(4), lp(3), lp(4)];
end

sgtitle(sprintf('MTF vs compensator error - %d files', nValid));
exportPng(sensFig, fullfile(outFolder, 'Compensator_Sensitivity.png'));

%% Plots 2..: report-yield MTF for each error level (left: all fields, right: on-axis)
% Solid: tangential, dashed: sagittal, black: diffraction limit
fieldColors = lines(nFields);
showDiffS = max(abs(diffT - diffS)) > 1e-6;
levelFigs = gobjects(nLevels, 1);

for lev = 1:nLevels
    if lev == 1
        levelName = 'No compensator error';
    else
        levelName = sprintf('Compensator error \\pm%g \\mum (worse side)', levelErrors_um(lev));
    end

    fig = figure('Name', sprintf('MTF - error %g um', levelErrors_um(lev)), 'Color', 'w', ...
        'Position', [50 50 1500 750]);

    panels = {1:nFields, onAxisField};
    panelTitles = ["All fields", "On-axis field: " + fieldLabels(onAxisField)];

    for p = 1:2
        ax = axes(fig, 'Position', axPos{p});
        hold on;
        plot(freqGrid, diffT, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Diff. limit');
        if showDiffS
            plot(freqGrid, diffS, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Diff. limit S');
        end
        for f = panels{p}
            plot(freqGrid, yieldT(:, ry, f, lev), '-', 'Color', fieldColors(f, :), ...
                'LineWidth', 1.2, 'DisplayName', sprintf('F%d T', f));
            plot(freqGrid, yieldS(:, ry, f, lev), '--', 'Color', fieldColors(f, :), ...
                'LineWidth', 1.2, 'DisplayName', sprintf('F%d S', f));
        end
        hold off;
        xlim([0 maxFreq]);
        ylim([0 1]);
        xlabel('Spatial frequency (cyc/mm)');
        ylabel('Modulus of the OTF');
        title(panelTitles(p), 'Interpreter', 'none');
        grid on;
        box on;
        pbaspect(ax, [1.6 1 1]);

        lgd = legend(ax, 'NumColumns', 6, 'Location', 'none');
        lgd.Units = 'normalized';
        lp = lgd.Position;
        lgd.Position = [axPos{p}(1) + (axPos{p}(3) - lp(3)) / 2, 0.24 - lp(4), lp(3), lp(4)];
    end

    sgtitle(sprintf('Polychromatic FFT MTF - %d%% yield - %s (%d files)', ...
        reportYield, levelName, nValid));
    levelFigs(lev) = fig;

    exportPng(fig, fullfile(outFolder, sprintf('MTF_Yield%d_CompError_%gum.png', ...
        reportYield, levelErrors_um(lev))));
end

%% Export Excel
outputFile = fullfile(outFolder, 'Compensator_Results.xlsx');

% Remove old file so stale rows/sheets are not kept
if isfile(outputFile)
    delete(outputFile);
end

% Sensitivity table
Sens = table(levelErrors_um', 'VariableNames', {'CompError_um'});
for q = 1:2
    tag = ["NyqHalf", "Nyq"];
    for y = 1:nYield
        Sens.(sprintf('Yield%d_%s', yieldLevels(y), tag(q))) = sensYield(:, y, q);
    end
    Sens.(sprintf('Mean_%s', tag(q))) = sensMean(:, q);
end
writetable(Sens, outputFile, 'Sheet', 'Sensitivity');

% Report-yield MTF curves for each error level
for lev = 1:nLevels
    LevelTable = table(freqGrid, 'VariableNames', {'Frequency_cyc_mm'});
    for f = 1:nFields
        LevelTable.(sprintf('F%d_T', f)) = yieldT(:, ry, f, lev);
        LevelTable.(sprintf('F%d_S', f)) = yieldS(:, ry, f, lev);
    end
    writetable(LevelTable, outputFile, 'Sheet', sprintf('Yield%d Err %gum', reportYield, levelErrors_um(lev)));
end

% Worse side (+1 / -1) for each file and error
WorstSide = table(fileNames, compThickness, 'VariableNames', {'FileName', 'CompThickness'});
for e = 1:nErr
    WorstSide.(sprintf('Err%gum', compErrors_um(e))) = worstSign(:, e);
end
WorstSide.Status = status;
writetable(WorstSide, outputFile, 'Sheet', 'Worse Side');

FieldInfo = table((1:nFields)', fieldLabels, 'VariableNames', {'Field', 'Description'});
writetable(FieldInfo, outputFile, 'Sheet', 'Field Info');

DiffTable = table(freqGrid, diffT, diffS, ...
    'VariableNames', {'Frequency_cyc_mm', 'DiffLimit_T', 'DiffLimit_S'});
writetable(DiffTable, outputFile, 'Sheet', 'Diffraction Limit');

% Raw per-file MTF curves for further analysis
save(fullfile(outFolder, 'Compensator_RawData.mat'), 'freqGrid', 'MTF_T', 'MTF_S', ...
    'levelErrors_um', 'worstSign', 'compThickness', 'compSurface', 'diffT', 'diffS', ...
    'fileNames', 'status', 'fieldLabels', 'yieldLevels', 'pixelPitch_um', 'maxFreq', ...
    'nominalName', 'nominalAfocal', 'paraxialEFL');

%% PDF report: summary page + sensitivity page + one page per error level
pdfFile = fullfile(outFolder, 'Compensator_Report.pdf');

% Remove old report files (single and per-page)
oldPages = dir(fullfile(outFolder, 'Compensator_Report_page*.pdf'));
for i = 1:numel(oldPages)
    delete(fullfile(outFolder, oldPages(i).name));
end
if isfile(pdfFile)
    delete(pdfFile);
end

% If the old report could not be deleted (e.g. open in a PDF viewer),
% save the new one under a time-stamped name
if isfile(pdfFile)
    pdfFile = fullfile(outFolder, ['Compensator_Report_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.pdf']);
    fprintf('\nCompensator_Report.pdf is open in another program; new report: %s\n', pdfFile);
end

% Summary page text
if nominalAfocal
    afocalText = sprintf(['Yes - MTF computed in cyc/mrad, frequency scaled ' ...
        'to cyc/mm by 1000/|EFL| (|EFL| = %.4f mm)'], paraxialEFL);
else
    afocalText = 'No';
end

header = sprintf('%-14s%-14s%14s%14s%14s%14s', 'Error (um)', 'Worse +/-', ...
    'Field Nyq/2', 'Field Nyq', 'Axis Nyq/2', 'Axis Nyq');

report = {
    'COMPENSATOR ERROR - MTF ANALYSIS REPORT'
    ''
    sprintf('Date:              %s', char(datetime('now', 'Format', 'dd.MM.yyyy HH:mm')))
    sprintf('Folder:            %s', folderPath)
    sprintf('File prefix:       %s', filePrefix)
    sprintf('Valid files:       %d / %d', nValid, nFiles)
    sprintf('Nominal file:      %s', nominalName)
    sprintf('Afocal mode:       %s', afocalText)
    sprintf('Compensator:       thickness of surface %d %s', compSurface, char(compSurfaceComment))
    sprintf('Errors:            +/- %s um (worse side kept for each file)', strjoin(string(compErrors_um), ', '))
    sprintf('Pixel pitch:       %.2f um', pixelPitch_um)
    sprintf('Max frequency:     %.2f cyc/mm (Nyquist)', maxFreq)
    sprintf('Yield levels:      %s %%', strjoin(string(yieldLevels), ', '))
    'MTF type:          Polychromatic FFT MTF (all wavelengths, all fields)'
    sprintf('Diffraction limit: %s', char(diffSource))
    sprintf('Yield X%% = MTF value that at least X%% of the systems meet or exceed')
    ''
    sprintf('%d%% yield MTF (Field: worst field and T/S, Axis: on-axis field, worse of T/S)', reportYield)
    ''
    header
    repmat('-', 1, numel(header))
    };

for lev = 1:nLevels
    if lev == 1
        sideText = '-';
    else
        e = lev - 1;
        sideText = sprintf('%d / %d', sum(worstSign(valid, e) == 1), sum(worstSign(valid, e) == -1));
    end
    report{end + 1, 1} = sprintf('%-14g%-14s%14.4f%14.4f%14.4f%14.4f', levelErrors_um(lev), sideText, ...
        sensYield(lev, ry, 1), sensYield(lev, ry, 2), sensOnAxis(lev, 1), sensOnAxis(lev, 2)); %#ok<SAGROW>
end

report = [report; {''; 'Fields (plots: solid = tangential, dashed = sagittal)'}];
report = [report; cellstr(compose("  F%d: %s", (1:nFields)', fieldLabels))];

failed = fileNames(~valid);
if ~isempty(failed)
    report = [report; {''; sprintf('Failed files (%d):', numel(failed))}];
    shown = failed(1:min(end, 10));
    report = [report; cellstr("  " + shown)];
    if numel(failed) > numel(shown)
        report{end + 1, 1} = sprintf('  ... and %d more (see Worse Side sheet in Excel)', numel(failed) - numel(shown));
    end
end

reportFig = figure('Name', 'Compensator Report', 'Color', 'w', 'Position', [50 50 1200 850]);
ax = axes(reportFig, 'Position', [0 0 1 1], 'Visible', 'off');
text(ax, 0.04, 0.97, report, 'VerticalAlignment', 'top', ...
    'FontName', 'Courier New', 'FontSize', 9, 'Interpreter', 'none');

allFigs = [reportFig; sensFig; levelFigs];

try
    % Multi-page PDF (exportgraphics with 'Append' needs R2021b+)
    exportgraphics(allFigs(1), pdfFile, 'ContentType', 'vector');
    for p = 2:numel(allFigs)
        exportgraphics(allFigs(p), pdfFile, 'ContentType', 'vector', 'Append', true);
    end
    fprintf('\nPDF report saved to:\n%s\n', pdfFile);
catch pdfErr
    % Older MATLAB: one PDF file per page
    fprintf('\nMulti-page PDF export failed: %s\n', pdfErr.message);
    if isfile(pdfFile)
        delete(pdfFile);
    end
    for p = 1:numel(allFigs)
        set(allFigs(p), 'PaperOrientation', 'landscape', 'PaperPositionMode', 'auto');
        print(allFigs(p), fullfile(outFolder, sprintf('Compensator_Report_page%d.pdf', p)), '-dpdf', '-bestfit');
    end
    fprintf('Report saved as separate pages (Compensator_Report_page*.pdf).\n');
end

fprintf('\nResults saved to:\n%s\n', outFolder);


function [T, S, labels, dT, dS] = computeMtf(mtf, freqGrid, freqScale)
% Runs the FFT MTF analysis on the current system and returns tangential
% and sagittal MTF (nFreq x nFields) on the common frequency grid, plus the
% diffraction limit if OpticStudio returned it.
mtf.ApplyAndWaitForCompletion();
results = mtf.GetResults();
nSeries = results.NumberOfDataSeries;
dataSeries = results.DataSeries;

isDiff = false(nSeries, 1);
for s = 1:nSeries
    isDiff(s) = contains(lower(string(char(dataSeries(s).Description))), "diff");
end
fieldIdx = find(~isDiff);

T = nan(numel(freqGrid), numel(fieldIdx));
S = T;
labels = strings(numel(fieldIdx), 1);
for f = 1:numel(fieldIdx)
    data = dataSeries(fieldIdx(f));
    [T(:, f), S(:, f)] = readSeries(data, freqGrid, freqScale);
    labels(f) = strtrim(string(char(data.Description)));
end

dT = [];
dS = [];
diffIdx = find(isDiff, 1);
if ~isempty(diffIdx)
    [dT, dS] = readSeries(dataSeries(diffIdx), freqGrid, freqScale);
end
end

function [t, s] = readSeries(data, freqGrid, freqScale)
% Reads one FFT MTF data series (column 1: tangential, column 2: sagittal),
% converts its frequencies to cyc/mm (freqScale) and interpolates it to the
% common frequency grid.
x = data.XData.Data.double;
y = data.YData.Data.double;

x = x(:) * freqScale;
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

function surfNo = findCompensator(TheSystem)
% Reads the compensator surface from the COMP operand in the Tolerance
% Data Editor. Only a single thickness compensator (code 0) is supported.
tde = TheSystem.TDE;
found = [];
for i = 1:tde.NumberOfOperands
    op = tde.GetOperandAt(i);
    if strcmpi(string(op.Type.ToString()), "COMP")
        if op.Param2 ~= 0
            error('Compensator in TDE row %d is not a thickness (code %d). Set compSurface manually.', i, op.Param2);
        end
        found(end + 1) = op.Param1; %#ok<AGROW>
    end
end

if isempty(found)
    error('No COMP operand found in the nominal file. Set compSurface in the user settings.');
elseif numel(found) > 1
    error('More than one COMP operand found (surfaces %s). Set compSurface in the user settings.', ...
        strjoin(string(found), ', '));
end
surfNo = found;
end

function f = unitsPerMm(TheSystem)
% Lens units per millimeter (compensator errors are given in um)
u = lower(string(TheSystem.SystemData.Units.LensUnits.ToString()));
if contains(u, "milli")
    f = 1;
elseif contains(u, "centi")
    f = 0.1;
elseif contains(u, "inch")
    f = 1 / 25.4;
elseif contains(u, "meter")
    f = 1e-3;
else
    f = 1;
end
end

function v = yieldValue(data, yieldPct)
% data: nSystems x nFreq. Returns, for each frequency, the MTF value that
% at least yieldPct % of the systems meet or exceed (no toolbox needed).
sorted = sort(data, 1, 'descend');
idx = max(1, ceil(yieldPct / 100 * size(data, 1)));
v = sorted(idx, :);
end

function exportPng(fig, imageFile)
try
    exportgraphics(fig, imageFile, 'Resolution', 300);   % R2020a+
catch
    saveas(fig, imageFile);
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
