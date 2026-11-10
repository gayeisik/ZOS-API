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
% (lower on-axis MTF: mean of T and S over all frequencies) is kept and the
% chosen signed error (e.g. +10 or -10 um) is reported.

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
worstSign = nan(nFiles, nErr);   % +1: +d was worse, -1: -d was worse (on-axis MTF)
compThickness = nan(nFiles, 1);  % compensated thickness of each MC file
fieldLabels = strings(0, 1);
nFields = 0;
diffT = [];
diffS = [];
diffSource = "";
paraxialEFL = NaN;
freqScale = 1;       % Zemax frequency -> cyc/mm (1000/|EFL| for afocal files)
checkPrinted = false;   % check output is printed for the first valid file only
maxChange = 0;          % largest MTF change caused by the errors in that file

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

    % Is the thickness set in the Lens Data Editor or by a multi-configuration
    % THIC operand? (An LDE value controlled by the MCE cannot be changed directly.)
    [compMceRow, compInfo] = findCompensatorControl(TheSystem, compSurface);

    lensUnitsPerMm = unitsPerMm(TheSystem);
    onAxisField = findOnAxisField(TheSystem);

    fprintf('Compensator: thickness of surface %d %s (%s)\n', compSurface, compSurfaceComment, compInfo);
    fprintf('Nyquist frequency: %.4f cyc/mm (pixel pitch %.2f um)\n', maxFreq, pixelPitch_um);

    tStart = tic;

    for k = 1:nFiles

        filePath = fullfile(folderPath, files(k).name);

        try
            % Load optical system
            if ~TheSystem.LoadFile(filePath, false)
                error('OpticStudio could not load this file.');
            end

            t0 = getCompensator(TheSystem, compSurface, compMceRow);
            compThickness(k) = t0;

            % Level 1: compensated file, no error (no optimization is run anywhere)
            [T0, S0, labels, dT, dS] = computeMtf(TheSystem, maxFreq, freqGrid, freqScale, isempty(diffT));

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
            end

            MTF_T(k, :, :, 1) = T0;
            MTF_S(k, :, :, 1) = S0;

            % Levels 2..end: +/- compensator error, keep the worse one
            for e = 1:nErr
                d = compErrors_um(e) / 1000 * lensUnitsPerMm;   % um -> lens units

                % Set the error, check that it was really applied, then compute MTF
                setCompensator(TheSystem, compSurface, compMceRow, t0 + d);
                [Tp, Sp] = computeMtf(TheSystem, maxFreq, freqGrid, freqScale, false);

                setCompensator(TheSystem, compSurface, compMceRow, t0 - d);
                [Tm, Sm] = computeMtf(TheSystem, maxFreq, freqGrid, freqScale, false);

                % Worse side decided by the on-axis field MTF
                onAxisPlus = mean([Tp(:, onAxisField); Sp(:, onAxisField)]);
                onAxisMinus = mean([Tm(:, onAxisField); Sm(:, onAxisField)]);

                if onAxisPlus <= onAxisMinus
                    MTF_T(k, :, :, e + 1) = Tp;
                    MTF_S(k, :, :, e + 1) = Sp;
                    worstSign(k, e) = 1;
                else
                    MTF_T(k, :, :, e + 1) = Tm;
                    MTF_S(k, :, :, e + 1) = Sm;
                    worstSign(k, e) = -1;
                end

                % Check output for the first file: does the MTF react to the error?
                if ~checkPrinted
                    if e == 1
                        fprintf('\nCheck on %s (compensator value %.6f):\n', files(k).name, t0);
                        fprintf('  on-axis mean MTF without error: %.4f\n', ...
                            mean([T0(:, onAxisField); S0(:, onAxisField)]));
                    end
                    fprintf('  +%g um: %.4f   -%g um: %.4f   -> applied %+g um\n', ...
                        compErrors_um(e), onAxisPlus, compErrors_um(e), onAxisMinus, ...
                        worstSign(k, e) * compErrors_um(e));
                    maxChange = max(maxChange, max(abs([Tp(:) - T0(:); Tm(:) - T0(:)])));
                end
            end

            if ~checkPrinted
                checkPrinted = true;
                if maxChange < 1e-6
                    warning(['The MTF did not change when the compensator thickness was changed. ' ...
                        'Surface %d is probably not the real compensator, or its thickness has no ' ...
                        'effect (e.g. it is in collimated space). Check compSurface.'], compSurface);
                end
                fprintf('\n');
            end

            setCompensator(TheSystem, compSurface, compMceRow, t0);

            if any(isnan(MTF_T(k, :, :, :)), 'all') || any(isnan(MTF_S(k, :, :, :)), 'all')
                error('MTF data contains invalid values.');
            end

            status(k) = "OK";

        catch ME
            status(k) = string(ME.message);
            fprintf('ERROR: %s - %s\n', files(k).name, ME.message);
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

%% Plot pages (A4 portrait, two plots per page, every plot the same size)
plotPos = {[0.12 0.63 0.80 0.27], [0.12 0.20 0.80 0.27]};   % top / bottom plot
fieldColors = lines(nFields);
showDiffS = max(abs(diffT - diffS)) > 1e-6;

% One page per error level: all fields (top), on-axis field (bottom)
levelFigs = gobjects(nLevels, 1);

for lev = 1:nLevels
    if lev == 1
        levelName = 'Kompanzatör hatası yok (0 µm)';
    else
        e = lev - 1;
        levelName = sprintf('Kompanzatör hatası ±%g µm (uygulanan: +%g → %d dosya, −%g → %d dosya)', ...
            compErrors_um(e), compErrors_um(e), sum(worstSign(valid, e) == 1), ...
            compErrors_um(e), sum(worstSign(valid, e) == -1));
    end

    fig = newA4Page();
    addPageTitle(fig, sprintf('%d%% yield MTF – %s', reportYield, levelName));

    panels = {1:nFields, onAxisField};
    panelTitles = ["All fields", "On-axis field: " + fieldLabels(onAxisField)];

    for p = 1:2
        ax = axes(fig, 'Position', plotPos{p}, 'FontSize', 9);
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
        placeLegend(ax, plotPos{p}, 6);
    end

    levelFigs(lev) = fig;
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

% Applied (worse-side) signed compensator error in um for each file, e.g. +10 or -10
WorstSide = table(fileNames, compThickness, 'VariableNames', {'FileName', 'CompThickness'});
for e = 1:nErr
    WorstSide.(sprintf('Err%gum_applied_um', compErrors_um(e))) = worstSign(:, e) * compErrors_um(e);
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

%% PDF report (A4): summary, explanation pages, one MTF page per error level
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

errListText = strjoin(string(compErrors_um), ', ');

if nominalAfocal
    afocalText = sprintf('Evet – MTF cyc/mrad olarak hesaplandı, 1000/|EFL| ile cyc/mm''ye çevrildi (|EFL| = %.4f mm)', paraxialEFL);
    afocalExplain = sprintf(['Bu analizde nominal dosya afokal moddadır: |EFL| = %.4f mm, ' ...
        'ölçek katsayısı 1000/|EFL| = %.4f.'], paraxialEFL, freqScale);
else
    afocalText = 'Hayır';
    afocalExplain = 'Bu analizde nominal dosya afokal modda değildir; MTF doğrudan cyc/mm olarak alınır.';
end

% --- Page: summary -------------------------------------------------------
summaryLines = {
    sprintf('Tarih:               %s', char(datetime('now', 'Format', 'dd.MM.yyyy HH:mm')))
    sprintf('Klasör:              %s', folderPath)
    sprintf('Nominal dosya:       %s', nominalName)
    sprintf('MC dosyaları:        %s*  (geçerli: %d / %d)', filePrefix, nValid, nFiles)
    sprintf('Afokal mod:          %s', afocalText)
    sprintf('Kompanzatör:         yüzey %d kalınlığı %s', compSurface, char(compSurfaceComment))
    sprintf('Kompanzatör hatası:  ±%s µm', errListText)
    sprintf('Piksel aralığı:      %.2f µm', pixelPitch_um)
    sprintf('Nyquist frekansı:    %.2f cyc/mm  (Nyquist/2 = %.2f cyc/mm)', maxFreq, maxFreq / 2)
    sprintf('Yield seviyeleri:    %s %%', strjoin(string(yieldLevels), ', '))
    'MTF tipi:            Polikromatik FFT MTF (tüm dalga boyları, tüm fieldler)'
    sprintf('Diffraction limit:   %s', char(diffSource))
    ''
    sprintf('%d%% YIELD MTF ÖZETİ', reportYield)
    'Field: en kötü field ve en kötü yön (T/S)    Axis: on-axis field, T/S''nin kötüsü'
    ''
    };

header = sprintf('%-10s%-26s%11s%11s%11s%11s', 'Hata(µm)', 'Uygulanan (dosya sayısı)', ...
    'Field N/2', 'Field N', 'Axis N/2', 'Axis N');
summaryLines = [summaryLines; {header; repmat('-', 1, numel(header))}];

for lev = 1:nLevels
    if lev == 1
        sideText = '-';
    else
        e = lev - 1;
        sideText = sprintf('+%g: %d   -%g: %d', compErrors_um(e), sum(worstSign(valid, e) == 1), ...
            compErrors_um(e), sum(worstSign(valid, e) == -1));
    end
    summaryLines{end + 1, 1} = sprintf('%-10g%-26s%11.4f%11.4f%11.4f%11.4f', levelErrors_um(lev), sideText, ...
        sensYield(lev, ry, 1), sensYield(lev, ry, 2), sensOnAxis(lev, 1), sensOnAxis(lev, 2)); %#ok<SAGROW>
end

summaryLines = [summaryLines; {''; 'N/2 = Nyquist/2, N = Nyquist'; ''; 'FIELDLER (grafiklerde düz çizgi = T, kesikli çizgi = S)'}];
summaryLines = [summaryLines; cellstr(compose("  F%d: %s", (1:nFields)', fieldLabels))];

failed = fileNames(~valid);
if ~isempty(failed)
    summaryLines = [summaryLines; {''; sprintf('Hesaplanamayan dosyalar (%d):', numel(failed))}];
    shown = failed(1:min(end, 8));
    summaryLines = [summaryLines; cellstr("  " + shown)];
    if numel(failed) > numel(shown)
        summaryLines{end + 1, 1} = sprintf('  ... ve %d dosya daha (Excel, Worse Side sayfası)', numel(failed) - numel(shown));
    end
end

summaryFig = textPage('KOMPANZATÖR HATASI – MTF ANALİZ RAPORU', summaryLines, 'Courier New', 8.5);

% --- Explanation pages ---------------------------------------------------
explain1 = {
    'AMAÇ'
    ['Monte Carlo (MC) tolerans analizinde her MC_T dosyası, toleranslar uygulandıktan sonra ' ...
     'kompanzatör en iyi konuma optimize edilerek kaydedilir. Gerçek sistemde ise kompanzatör bu ' ...
     'ideal konuma tam olarak getirilemez: motor adımı, mekanik boşluk (backlash), odak ayarındaki ' ...
     'ölçüm belirsizliği ve ısıl kaymalar küçük bir konumlandırma hatası bırakır. Bu analiz, ' ...
     'kompanzatör ±δ kadar hatalı konumlandığında MTF''in nasıl değiştiğini gösterir ve kompanzatör ' ...
     'için izin verilebilecek hata payının belirlenmesini sağlar.']
    ''
    'ADIM 1 – NOMİNAL DOSYA'
    sprintf(['Klasörde adı MC_ ile başlamayan dosya nominal tasarım kabul edilir (bu analizde: %s). ' ...
     'Bu dosyadan şunlar okunur:'], nominalName)
    ['•  Afokal mod: Sistem "Afocal Image Space" modundaysa Zemax MTF''i açısal frekansta (cyc/mrad) ' ...
     'verir. Bu durumda |EFL| okunur ve frekanslar  ν[cyc/mm] = ν[cyc/mrad] × 1000 / |EFL|  ile ' ...
     'cyc/mm''ye çevrilir. MTF değerleri değişmez, sadece frekans ekseni ölçeklenir. ' afocalExplain]
    sprintf(['•  Kompanzatör: Tolerans editöründeki COMP satırından okunan yüzeyin kalınlığıdır ' ...
     '(bu analizde: yüzey %d).'], compSurface)
    '•  On-axis field: Koordinatı (0, 0) olan fielddir.'
    ''
    'ADIM 2 – REFERANS MTF (0 µm)'
    sprintf(['Her MC_T dosyası kaydedildiği (kompanze edilmiş) haliyle açılır ve polikromatik FFT MTF ' ...
     'hesaplanır: tüm dalga boyları, tüm fieldler, tanjansiyel (T) ve sajital (S) ayrı ayrı, 0''dan ' ...
     'Nyquist frekansına (%.2f cyc/mm) kadar. Bu sonuç "0 µm" referansıdır ve kompanzatör hatası ' ...
     'olmayan durumu gösterir.'], maxFreq)
    ''
    'ADIM 3 – KOMPANZATÖR HATASI (±δ)'
    sprintf(['Her hata payı δ (%s µm) için kompanzatör kalınlığı önce +δ, sonra −δ kaydırılır ve her iki ' ...
     'durumda MTF yeniden hesaplanır. İki sonuçtan on-axis fieldin ortalama MTF''i (T ve S, tüm ' ...
     'frekanslar) daha düşük olan taraf "kötü taraf" kabul edilir ve o dosya için bu sonuç ' ...
     'kullanılır. Böylece her dosyada en kötü durum (worst case) alınmış olur.'], errListText)
    ['Her dosyada hangi tarafın (+δ veya −δ) uygulandığı Excel''deki "Worse Side" sayfasında, toplam ' ...
     'dosya sayıları ise özet tablosunda ve her MTF sayfasının başlığında verilir.']
    ['Hata eklendikten sonra hiçbir optimizasyon yapılmaz; MTF doğrudan hatalı kompanzatör konumunda ' ...
     'hesaplanır. Kod her değişiklikten sonra değerin sisteme gerçekten uygulandığını kontrol eder.']
    ['Not: Kötü taraf on-axis fielde göre seçildiği için bazı kenar fieldlerde diğer taraf daha kötü ' ...
     'olabilir. Hesaptan sonra kompanzatör dosyanın kendi değerine geri alınır; diskteki dosyalar ' ...
     'değiştirilmez.']
    ''
    'ADIM 4 – İSTATİSTİK: ORTALAMA VE YIELD MTF'
    'Her hata payı için tüm geçerli dosyaların MTF''leri her frekansta ayrı ayrı istatistiğe sokulur:'
    '•  Ortalama MTF: Dosyaların MTF değerlerinin ortalaması.'
    ['•  Yield %X MTF: Dosyaların en az %X''inin sağladığı (bu değere eşit veya üstünde kalan) MTF ' ...
     'değeri. Örneğin 200 dosyada %90 yield MTF = 0.25 ise sistemlerin en az 180''inin MTF''i 0.25 ' ...
     'veya üstündedir. Üretim açısından "sistemlerin %90''ı en az bu performansı verir" anlamına gelir.']
    ['Önemli: Yield eğrisi her frekansta ayrı hesaplanır. Bu yüzden %90 yield eğrisi tek bir gerçek ' ...
     'sistemin MTF eğrisi değildir; her frekanstaki alt sınırların birleşimidir.']
    };

explain2 = {
    'NYQUIST FREKANSI NEDİR?'
    sprintf(['Nyquist frekansı, dedektörün örnekleyebileceği en yüksek uzaysal frekanstır:  ' ...
     'ν_N = 1 / (2 × piksel aralığı). Bir çizgi çifti (bir açık + bir koyu çizgi) iki piksele denk ' ...
     'gelir. Bu analizde piksel aralığı %.2f µm olduğu için Nyquist = %.2f cyc/mm''dir. Nyquist''in ' ...
     'üzerindeki detaylar dedektör tarafından doğru örneklenemez (aliasing); bu yüzden MTF yalnızca ' ...
     'Nyquist''e kadar hesaplanır.'], pixelPitch_um, maxFreq)
    ['•  Nyquist''teki MTF: Sistemin, dedektörün çözebileceği en ince detayı ne kadar kontrastla ' ...
     'ilettiğini gösterir. Küçük hedeflerin ve ince çizgilerin seçilebilirliği buna bağlıdır.']
    ['•  Nyquist/2''deki MTF: Orta frekanslardaki kontrastı, yani görüntünün genel keskinliğini ve ' ...
     'algılanan netliğini temsil eder.']
    ['Kompanzatör hatası esas olarak odak kayması (defocus) üretir. Defocus yüksek frekansları daha çok ' ...
     'etkilediği için Nyquist''teki MTF genellikle Nyquist/2''dekinden daha hızlı düşer.']
    ''
    'KOMPANZATÖR HATA PAYI NASIL BELİRLENİR?'
    '1.  MTF şartını belirle, örneğin: "%90 yield''de Nyquist''te MTF ≥ 0.15".'
    ['2.  Özet sayfasındaki tabloda ve MTF sayfalarında, hata payı büyüdükçe MTF''in bu şartın altına ' ...
     'indiği ilk hata payını bul. Şartı hâlâ sağlayan en büyük hata payı, izin verilebilecek en büyük ' ...
     'kompanzatör hatasıdır. İki hata payı arasında doğrusal değişim kabul edilebilir.']
    ['3.  Güvenlik payı bırak: bulunan değerin biraz altını spesifikasyon olarak seç (diğer hata ' ...
     'kaynakları ve modelleme belirsizlikleri için).']
    ['4.  Bu değeri mekanizmanın gerçek kabiliyetiyle (motor adımı, boşluk, odak ayarı hassasiyeti, ' ...
     'sıcaklık) karşılaştır. Mekanizma bu hassasiyeti sağlayamıyorsa mekanizma ya da optik tasarım ' ...
     '(odak derinliği) iyileştirilmelidir.']
    ''
    'MTF EĞRİSİ SAYFALARI NASIL OKUNUR?'
    sprintf(['Her hata payı için bir sayfa vardır. Üstteki grafik tüm fieldlerin, alttaki grafik ' ...
     'on-axis fieldin %d%% yield MTF eğrisini gösterir. Düz çizgi tanjansiyel (T), kesikli çizgi ' ...
     'sajital (S); siyah çizgi diffraction limit, yani aberasyonsuz ideal sistemin teorik üst ' ...
     'sınırıdır. Eğriler diffraction limitine ne kadar yakınsa performans o kadar iyidir.'], reportYield)
    ['Hata payı arttıkça eğrilerin aşağı inmesi beklenir. Aynı fieldde T ve S eğrilerinin birbirinden ' ...
     'ayrılması astigmatizma / field curvature etkisini gösterir. Sayfa başlığında, o hata payında ' ...
     'kaç dosyaya +δ ve kaç dosyaya −δ uygulandığı yazar.']
    };

explainFig1 = textPage('1. ANALİZ YÖNTEMİ', explain1, 'Arial', 10);
explainFig2 = textPage('2. NYQUIST FREKANSI, HATA PAYI VE MTF SAYFALARI', explain2, 'Arial', 10);

% --- Assemble pages: frame, header and page numbers ----------------------
allFigs = [summaryFig; explainFig1; explainFig2; levelFigs];
nPages = numel(allFigs);
headerText = sprintf('Kompanzatör Hatası – MTF Analizi  |  %s', nominalName);

for p = 1:nPages
    addPageFrame(allFigs(p), headerText, p, nPages);
end

% PNG of the plot pages
for lev = 1:nLevels
    exportPng(levelFigs(lev), fullfile(outFolder, sprintf('MTF_Yield%d_CompError_%gum.png', ...
        reportYield, levelErrors_um(lev))));
end

try
    % Multi-page PDF (exportgraphics with 'Append' needs R2021b+)
    exportgraphics(allFigs(1), pdfFile, 'ContentType', 'vector');
    for p = 2:nPages
        exportgraphics(allFigs(p), pdfFile, 'ContentType', 'vector', 'Append', true);
    end
    fprintf('\nPDF report saved to:\n%s\n', pdfFile);
catch pdfErr
    % Older MATLAB: one A4 PDF file per page
    fprintf('\nMulti-page PDF export failed: %s\n', pdfErr.message);
    if isfile(pdfFile)
        delete(pdfFile);
    end
    for p = 1:nPages
        print(allFigs(p), fullfile(outFolder, sprintf('Compensator_Report_page%d.pdf', p)), '-dpdf');
    end
    fprintf('Report saved as separate pages (Compensator_Report_page*.pdf).\n');
end

close(allFigs);

fprintf('\nResults saved to:\n%s\n', outFolder);


function [T, S, labels, dT, dS] = computeMtf(TheSystem, maxFreq, freqGrid, freqScale, showDiffLimit)
% Runs a new FFT MTF analysis on the current system and returns tangential
% and sagittal MTF (nFreq x nFields) on the common frequency grid, plus the
% diffraction limit if requested and returned by OpticStudio. A new analysis
% is opened for every call so that no earlier result can be reused.
mtf = TheSystem.Analyses.New_FftMtf();
closeMtf = onCleanup(@() mtf.Close());

settings = mtf.GetSettings();
settings.MaximumFrequency = maxFreq / freqScale;   % cyc/mrad if afocal
settings.SampleSize = ZOSAPI.Analysis.SampleSizes.S_64x64;
settings.Wavelength.UseAllWavelengths();
settings.Field.UseAllFields();
try
    settings.ShowDiffractionLimit = showDiffLimit;
catch
end

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

function [mceRow, info] = findCompensatorControl(TheSystem, surfNo)
% Returns the Multi-Configuration Editor row of a THIC operand that controls
% the thickness of surfNo (0 if the thickness is set in the Lens Data Editor)
% and a short description for the log.
mceRow = 0;
info = "set in Lens Data Editor";

try
    mce = TheSystem.MCE;
    for i = 1:mce.NumberOfOperands
        op = mce.GetOperandAt(i);
        if strcmpi(string(op.Type.ToString()), "THIC") && op.Param1 == surfNo
            mceRow = i;
            info = sprintf("set by multi-configuration operand THIC (row %d), configuration %d of %d", ...
                i, mce.CurrentConfiguration, mce.NumberOfConfigurations);
            return;
        end
    end
catch
end

try
    solveType = string(TheSystem.LDE.GetSurfaceAt(surfNo).ThicknessCell.GetSolveData().Type.ToString());
    info = info + ", solve: " + solveType;
catch
end
end

function v = getCompensator(TheSystem, surfNo, mceRow)
if mceRow > 0
    mce = TheSystem.MCE;
    v = mce.GetOperandAt(mceRow).GetOperandCell(mce.CurrentConfiguration).DoubleValue;
else
    v = TheSystem.LDE.GetSurfaceAt(surfNo).Thickness;
end
end

function setCompensator(TheSystem, surfNo, mceRow, value)
% Sets the compensator thickness and checks that the system really took it
% (a pickup or solve on the thickness would silently restore the old value).
% The .NET object is stored in a variable first: MATLAB does not apply a
% property assignment written on a method call chain (obj.Method(x).Prop = v).
if mceRow > 0
    mce = TheSystem.MCE;
    compCell = mce.GetOperandAt(mceRow).GetOperandCell(mce.CurrentConfiguration);
    compCell.DoubleValue = value;
else
    surf = TheSystem.LDE.GetSurfaceAt(surfNo);
    surf.Thickness = value;
end

actual = TheSystem.LDE.GetSurfaceAt(surfNo).Thickness;
if abs(actual - value) > 1e-9
    error(['Compensator thickness of surface %d could not be changed (wanted %.6f, system has %.6f). ' ...
        'The thickness is probably controlled by a solve, pickup or multi-configuration operand.'], ...
        surfNo, value, actual);
end
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

function fig = newA4Page()
% Blank A4 portrait page (21 x 29.7 cm). Invisible, so its size is not
% limited by the screen.
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
    'Position', [1 1 21 29.7], 'PaperUnits', 'centimeters', ...
    'PaperSize', [21 29.7], 'PaperPosition', [0 0 21 29.7]);
end

function addPageTitle(fig, titleText)
annotation(fig, 'textbox', [0.06 0.915 0.88 0.04], 'String', titleText, ...
    'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'Interpreter', 'none');
end

function fig = textPage(titleText, lines, fontName, fontSize)
% A4 page with a title and wrapped text (one cell per paragraph/line)
fig = newA4Page();
annotation(fig, 'textbox', [0.08 0.89 0.84 0.05], 'String', titleText, ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'Interpreter', 'none', 'VerticalAlignment', 'top');
annotation(fig, 'textbox', [0.08 0.06 0.84 0.83], 'String', lines, ...
    'FontName', fontName, 'FontSize', fontSize, 'EdgeColor', 'none', ...
    'Interpreter', 'none', 'VerticalAlignment', 'top', 'FitBoxToText', 'off');
end

function placeLegend(ax, axPos, nCols)
% Legend centred below the axes; placed manually so it does not resize the axes
lgd = legend(ax, 'NumColumns', nCols, 'Location', 'none', 'FontSize', 8);
lgd.Units = 'normalized';
lp = lgd.Position;
lgd.Position = [axPos(1) + (axPos(3) - lp(3)) / 2, axPos(2) - 0.045 - lp(4), lp(3), lp(4)];
end

function addPageFrame(fig, headerText, pageNo, nPages)
% Header, footer with page number and an almost-white page border. The
% border keeps the whole A4 page when exportgraphics crops to the content.
ax = axes(fig, 'Position', [0 0 1 1], 'Visible', 'off', 'XLim', [0 1], 'YLim', [0 1]);
hold(ax, 'on');
plot(ax, [0.001 0.999 0.999 0.001 0.001], [0.001 0.001 0.999 0.999 0.001], ...
    '-', 'Color', [0.995 0.995 0.995], 'LineWidth', 0.5);
plot(ax, [0.06 0.94], [0.965 0.965], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
plot(ax, [0.06 0.94], [0.035 0.035], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
text(ax, 0.06, 0.970, headerText, 'FontSize', 8, 'Color', [0.4 0.4 0.4], ...
    'VerticalAlignment', 'bottom', 'Interpreter', 'none');
text(ax, 0.94, 0.030, sprintf('Sayfa %d / %d', pageNo, nPages), 'FontSize', 8, ...
    'Color', [0.4 0.4 0.4], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'Interpreter', 'none');
hold(ax, 'off');
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
