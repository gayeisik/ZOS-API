clc;
clear;
close all;

%% User settings
filePrefix    = 'MC_T';           % only files starting with this are used
pixelPitch_um = 15;               % detector pixel pitch (um)
nFreqPoints   = 100;              % number of frequency steps from 0 to maxFreq
yieldLevels   = [50 80 90 95];    % yield levels (%)
fieldLimitPct = 75;               % max relative field used (% of detector half-diagonal)
                                  % (MTF can be measured up to ~70 % field)
% FFT MTF sample size is set in the "Compute MTF" section (S_64x64 by default)
% Nominal file: the .zmx file in the folder whose name does NOT start with "MC_"

% Nyquist frequency of the detector, used as the maximum MTF frequency
maxFreq = 1000 / (2 * pixelPitch_um);   % cyc/mm

%% Ask user for detector format
answer = inputdlg({'Detector format - horizontal pixel count (e.g. 640):', ...
                   'Detector format - vertical pixel count (e.g. 512):'}, ...
                  'Detector Format', [1 50], {'640', '512'});

if isempty(answer)
    return;
end

nPixH = str2double(answer{1});
nPixV = str2double(answer{2});

if ~(isfinite(nPixH) && nPixH > 0 && nPixH == round(nPixH)) || ...
   ~(isfinite(nPixV) && nPixV > 0 && nPixV == round(nPixV))
    error('Pixel counts must be positive integers.');
end

% Detector size; full field (100 %) = detector half-diagonal
detW_mm = nPixH * pixelPitch_um / 1000;
detH_mm = nPixV * pixelPitch_um / 1000;
halfDiag_mm = hypot(detW_mm, detH_mm) / 2;

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
paraxialEFL = NaN;   % focal length of the paraxial lens added to afocal files (mm)
nConverted = 0;      % number of afocal MC files converted with a paraxial lens

fprintf('Computing MTF for %d files...\n', nFiles);
wb = waitbar(0, 'Computing MTF...', 'Name', 'MTF Analysis');

try
    TheSystem = TheApplication.PrimarySystem;

    %% Nominal file: afocal check and EFL
    if ~TheSystem.LoadFile(fullfile(folderPath, nominalName), false)
        error('OpticStudio could not load the nominal file: %s', nominalName);
    end

    nominalAfocal = logical(TheSystem.SystemData.Aperture.AFocalImageSpace);

    % EFL is read in the file's own mode (also used for angle fields)
    nominalEFL = TheSystem.MFE.GetOperandValue( ...
        ZOSAPI.Editors.MFE.MeritOperandType.EFFL, 0, 0, 0, 0, 0, 0, 0, 0);

    if nominalAfocal
        if ~isfinite(nominalEFL) || nominalEFL == 0
            error('Could not compute the EFL of the nominal file: %s', nominalName);
        end

        paraxialEFL = abs(nominalEFL);
        fprintf('Nominal file %s is in afocal mode.\n', nominalName);
        fprintf('Paraxial lens focal length (|EFL|): %.4f mm\n', paraxialEFL);
    else
        fprintf('Nominal file %s is not in afocal mode.\n', nominalName);
    end

    fprintf('Nyquist frequency: %.4f cyc/mm (pixel pitch %.2f um)\n', maxFreq, pixelPitch_um);

    %% Field selection (done once, on the nominal file)
    [selFields, allRelPct, fieldTypeName] = selectFields(TheSystem, halfDiag_mm, nominalEFL, fieldLimitPct);
    nominalNumFields = TheSystem.SystemData.Fields.NumberOfFields;
    nominalOnAxis = findOnAxisField(TheSystem);
    fieldRelPct = allRelPct(selFields);

    fprintf('\nDetector: %d x %d px, %.4f x %.4f mm (half-diagonal %.4f mm)\n', ...
        nPixH, nPixV, detW_mm, detH_mm, halfDiag_mm);
    fprintf('Field type: %s, field limit: %g %% of full field\n', fieldTypeName, fieldLimitPct);
    for i = 1:nominalNumFields
        if ismember(i, selFields), mark = 'used'; else, mark = 'excluded'; end
        fprintf('  Field %2d: %6.1f %% -> %s\n', i, allRelPct(i), mark);
    end

    if isempty(selFields)
        error('No field of the nominal file is within %g %% of the full field.', fieldLimitPct);
    end

    for k = 1:nFiles

        filePath = fullfile(folderPath, files(k).name);
        mtf = [];

        try
            % Load optical system
            loaded = TheSystem.LoadFile(filePath, false);

            if ~loaded
                error('OpticStudio could not load this file.');
            end

            % Nominal is afocal (so MC files are too): add paraxial lens
            % before image, switch afocal mode off
            if nominalAfocal
                addParaxialFocus(TheSystem, paraxialEFL);
                nConverted = nConverted + 1;
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

            if numel(fieldIdx) ~= nominalNumFields
                error('Number of fields (%d) differs from nominal file (%d).', numel(fieldIdx), nominalNumFields);
            end

            % Keep only the fields selected on the nominal file
            fieldIdx = fieldIdx(selFields);

            if nFields == 0
                nFields = numel(selFields);
                MTF_T = nan(nFiles, nFreq, nFields);
                MTF_S = nan(nFiles, nFreq, nFields);
                fieldLabels = strings(nFields, 1);
                onAxisField = find(selFields == nominalOnAxis, 1);
                if isempty(onAxisField)
                    onAxisField = 1;
                end
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

%% Output folder: "MTF ANALİZİ - %75 FIELD" (char(304) = 'İ', avoids .m file encoding issues)
outFolder = fullfile(folderPath, ['MTF ANAL' char(304) 'Z' char(304) sprintf(' - %%%g FIELD', fieldLimitPct)]);

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
        'Position', [50 50 1500 750]);

    panels = {1:nFields, onAxisField};
    axPos = {[0.05 0.32 0.42 0.56], [0.54 0.32 0.42 0.56]};   % fixed axes boxes
    panelTitles = ["All fields", "On-axis field: " + fieldLabels(onAxisField)];

    for p = 1:2
        ax = axes(fig, 'Position', axPos{p});
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
        grid on;
        box on;

        % Same plot box size and aspect ratio (x:y = 1.6:1) for every plot
        pbaspect(ax, [1.6 1 1]);

        % Legend below the plot; placed manually so it does not resize the axes
        lgd = legend(ax, 'NumColumns', 6, 'Location', 'none');
        lgd.Units = 'normalized';
        lp = lgd.Position;
        lgd.Position = [axPos{p}(1) + (axPos{p}(3) - lp(3)) / 2, 0.24 - lp(4), lp(3), lp(4)];
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

%% Export Excel
outputFile = fullfile(outFolder, 'MTF_Results.xlsx');

% Remove old file so stale rows/sheets are not kept
if isfile(outputFile)
    delete(outputFile);
end

yieldNames = compose("Yield%d", yieldLevels);

for f = 1:nFields
    FieldTable = table(freqGrid, meanT(:, f), meanS(:, f), ...
        'VariableNames', {'Frequency_cyc_mm', 'Mean_T', 'Mean_S'});
    for y = 1:nYield
        FieldTable.(char(yieldNames(y) + "_T")) = yieldT(:, y, f);
        FieldTable.(char(yieldNames(y) + "_S")) = yieldS(:, y, f);
    end
    writetable(FieldTable, outputFile, 'Sheet', sprintf('Field %d', f));
end

FieldInfo = table((1:nFields)', selFields(:), fieldLabels, fieldRelPct(:), ...
    'VariableNames', {'Field', 'NominalFieldNo', 'Description', 'RelativeField_pct'});
writetable(FieldInfo, outputFile, 'Sheet', 'Field Info');

DiffTable = table(freqGrid, diffT, diffS, ...
    'VariableNames', {'Frequency_cyc_mm', 'DiffLimit_T', 'DiffLimit_S'});
writetable(DiffTable, outputFile, 'Sheet', 'Diffraction Limit');

Files = table(fileNames, status, 'VariableNames', {'FileName', 'Status'});
writetable(Files, outputFile, 'Sheet', 'Files');

% Raw per-file MTF curves for further analysis
save(fullfile(outFolder, 'MTF_RawData.mat'), 'freqGrid', 'MTF_T', 'MTF_S', ...
    'diffT', 'diffS', 'fileNames', 'status', 'fieldLabels', 'yieldLevels', ...
    'pixelPitch_um', 'maxFreq', 'nominalName', 'nominalAfocal', 'paraxialEFL', ...
    'nPixH', 'nPixV', 'fieldLimitPct', 'selFields', 'fieldRelPct');

%% PDF report: summary page + mean MTF page + one page per yield level
pdfFile = fullfile(outFolder, 'MTF_Report.pdf');

% Remove old report files (single and per-page)
oldPages = dir(fullfile(outFolder, 'MTF_Report_page*.pdf'));
for i = 1:numel(oldPages)
    delete(fullfile(outFolder, oldPages(i).name));
end
if isfile(pdfFile)
    delete(pdfFile);
end

% If the old report could not be deleted (e.g. open in a PDF viewer),
% save the new one under a time-stamped name
if isfile(pdfFile)
    pdfFile = fullfile(outFolder, ['MTF_Report_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.pdf']);
    fprintf('\nMTF_Report.pdf is open in another program; new report: %s\n', pdfFile);
end

% Summary page text
if nominalAfocal
    afocalText = sprintf(['Yes - paraxial lens (f = %.4f mm) added before image ' ...
        'in %d MC file(s), afocal mode off'], paraxialEFL, nConverted);
else
    afocalText = 'No';
end

report = {
    sprintf('MTF YIELD ANALYSIS REPORT - MAX FIELD %g %%', fieldLimitPct)
    ''
    sprintf('Date:              %s', char(datetime('now', 'Format', 'dd.MM.yyyy HH:mm')))
    sprintf('Folder:            %s', folderPath)
    sprintf('File prefix:       %s', filePrefix)
    sprintf('Valid files:       %d / %d', nValid, nFiles)
    sprintf('Detector format:   %d x %d px (%.4f x %.4f mm)', nPixH, nPixV, detW_mm, detH_mm)
    sprintf('Pixel pitch:       %.2f um', pixelPitch_um)
    sprintf('Field limit:       %g %% of full field (full field = detector half-diagonal %.4f mm)', fieldLimitPct, halfDiag_mm)
    sprintf('Fields used:       %d of %d nominal fields (field type: %s)', nFields, nominalNumFields, fieldTypeName)
    sprintf('Max frequency:     %.2f cyc/mm (Nyquist)', maxFreq)
    sprintf('Nominal file:      %s', nominalName)
    sprintf('Afocal mode:       %s', afocalText)
    sprintf('Yield levels:      %s %%', strjoin(string(yieldLevels), ', '))
    'MTF type:          Polychromatic FFT MTF (all wavelengths, selected fields)'
    sprintf('On-axis field:     Field %d (%s)', onAxisField, char(fieldLabels(onAxisField)))
    sprintf('Diffraction limit: %s', char(diffSource))
    ''
    'Fields (plots: solid = tangential, dashed = sagittal)'
    ''
    sprintf('Yield X%% = MTF value that at least X%% of the systems meet or exceed')
    };

fieldLines = compose("  F%d: %s  (nominal field %d, %.1f %% of full field)", ...
    (1:nFields)', fieldLabels, selFields(:), fieldRelPct(:));
idx = find(strcmp(report, 'Fields (plots: solid = tangential, dashed = sagittal)'), 1);
report = [report(1:idx); cellstr(fieldLines); report(idx + 1:end)];

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
catch pdfErr
    % Older MATLAB: one PDF file per page
    fprintf('\nMulti-page PDF export failed: %s\n', pdfErr.message);
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

function addParaxialFocus(TheSystem, focalLength)
% Inserts a paraxial lens before the image surface (focal length and
% thickness = focalLength) and switches afocal image space off, so the
% MTF is reported in cyc/mm. The file on disk is not changed.
lde = TheSystem.LDE;
imageIdx = lde.NumberOfSurfaces - 1;

lde.InsertNewSurfaceAt(imageIdx);
surf = lde.GetSurfaceAt(imageIdx);
surf.ChangeType(surf.GetSurfaceTypeSettings(ZOSAPI.Editors.LDE.SurfaceType.Paraxial));
surf.GetSurfaceCell(ZOSAPI.Editors.LDE.SurfaceColumn.Par1).DoubleValue = focalLength;
surf.Thickness = focalLength;
surf.Comment = 'Paraxial focus (MTF script)';

TheSystem.SystemData.Aperture.AFocalImageSpace = false;
end

function [sel, relPct, typeName] = selectFields(TheSystem, halfDiag_mm, efl, limitPct)
% Relative field (%) of each field = radial field / detector half-diagonal.
% Returns the field numbers within limitPct. Height fields (object or image
% height, X/Y in mm) are used directly as detector coordinates; angle
% fields are converted to image height with the EFL.
sysFields = TheSystem.SystemData.Fields;
typeName = string(char(sysFields.GetFieldType().ToString()));
n = sysFields.NumberOfFields;
r_mm = zeros(n, 1);

isAngle = contains(lower(typeName), "angle");
isHeight = contains(lower(typeName), "height");

if isAngle && (~isfinite(efl) || efl == 0)
    error('Angle fields need a valid EFL to convert to image height.');
end
if ~isAngle && ~isHeight
    error('Field type "%s" is not supported (use height or angle fields).', typeName);
end

for i = 1:n
    fld = sysFields.GetField(i);
    if isAngle
        r_mm(i) = abs(efl) * hypot(tand(fld.X), tand(fld.Y));
    else
        r_mm(i) = hypot(fld.X, fld.Y);
    end
end

relPct = 100 * r_mm / halfDiag_mm;
sel = find(relPct <= limitPct + 1e-9);
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
