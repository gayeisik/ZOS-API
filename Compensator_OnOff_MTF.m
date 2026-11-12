%% KOMPANZATÖRÜN KATKISI - KOMPANZATÖRLÜ / KOMPANZATÖRSÜZ MTF (ZOS-API, MATLAB)
% NE YAPAR
%   Aynı Monte Carlo dosyalarında MTF'i iki durumda hesaplar ve karşılaştırır:
%     - KOMPANZATÖRLÜ  : dosyanın kendi hali (Zemax kompanzatörü optimize etmiş)
%     - KOMPANZATÖRSÜZ : aynı dosya, kompanzatör nominal konumuna geri alınmış
%   Böylece toleransların toplam zararı, kompanzatörün geri kazandırdığı kısım
%   ve kompanzasyondan sonra kalan kayıp ayrı ayrı görülür. Örneğin MC yalnızca
%   afokal grubun toleranslarıyla koşulduysa, afokal gruptan gelen hatanın ne
%   kadarının kompanzatörle düzeltildiği gösterilir.
%
% NASIL ÇALIŞIR
%   1. Klasör seçilir. Adı "MC_" ile başlamayan .zmx dosyası nominal tasarımdır.
%      Ondan afokal mod, |EFL|, on-axis field (0,0), nominal MTF ve
%      kompanzatörün NOMİNAL kalınlığı okunur (kompanzatör = tek variable
%      kalınlık; bulunamazsa tolerans editöründeki COMP satırı).
%   2. Her MC_T dosyası açılır. Kompanzatör o dosyanın içinde ayrı bulunur:
%      kalınlığı variable olan tek yüzey. (Toleranslama MC dosyalarına coordinate
%      break eklediği için yüzey numaraları nominal dosyadan farklıdır.) Variable
%      kalınlık sayısı 1 değilse dosya hata verir ve analize katılmaz.
%   3. Dosyanın kendi haliyle MTF hesaplanır (kompanzatörlü).
%   4. Kompanzatör kalınlığı nominal değere getirilir, değerin uygulandığı
%      kontrol edilir ve MTF yeniden hesaplanır (kompanzatörsüz). Optimizasyon
%      yapılmaz. Sonra kompanzatör dosyanın kendi değerine geri alınır; diskteki
%      dosyalar değişmez.
%   İki durum aynı bozulmuş sistemden çıktığı için dosyalar birebir eşleşir;
%   ikinci bir Monte Carlo koşusuna ve sabit seed'e gerek yoktur.
%
% HESAPLAR
%   MTF: polikromatik FFT MTF, tüm dalga boyları, tüm fieldler, T ve S ayrı,
%        0'dan Nyquist'e kadar. Nyquist = 1000 / (2 × piksel aralığı[µm]) [cyc/mm]
%   Afokal dosyalarda frekans ekseni ν[cyc/mm] = ν[cyc/mrad] × 1000 / |EFL|.
%   Kompanzatör hareketi = MC dosyasındaki kalınlık − nominal kalınlık.
%   Dosya başına en kötü MTF = tüm fieldlerin T ve S MTF'inin en düşüğü.
%   Yield %X MTF = dosyaların en az %X'inin sağladığı MTF (her frekansta ayrı).
%   Geri kazanım (%) = (yield_kompanzatörlü − yield_kompanzatörsüz) /
%                      (nominal − yield_kompanzatörsüz) × 100
%     %100: kompanzatör toleransların zararını tamamen geri kazandırıyor.
%     %0  : kompanzatör hiçbir şey kazandırmıyor.
%
% AYARLAR (User settings)
%   filePrefix, pixelPitch_um, nFreqPoints, yieldLevels, reportYield,
%   keyFreqFrac, nFilesMax
%
% ÇIKTILAR (seçilen klasörde "KOMPANZATÖR KATKI ANALİZİ" klasörüne)
%   Compensator_OnOff_Report.pdf : A4; özet ve sonuçlar, yöntem, yorum, grafikler
%   *.png                        : grafik sayfaları
%   Compensator_OnOff_Results.xlsx : dosya bazında sonuçlar, yield eğrileri
%   Compensator_OnOff_RawData.mat  : her dosyanın iki durumdaki ham MTF eğrileri

clc;
clear;
close all;

%% User settings
filePrefix    = 'MC_T';             % only files starting with this are used
pixelPitch_um = 15;                 % detector pixel pitch (um)
nFreqPoints   = 100;                % number of frequency steps from 0 to maxFreq
yieldLevels   = [50 80 90 95];      % yield levels (%)
reportYield   = 90;                 % yield level used in the plots and the summary
keyFreqFrac   = [0.5 1.0];          % summary frequencies as fractions of Nyquist (first one is used in the plots)
nFilesMax     = Inf;                % use only the first N MC files (e.g. 20 for a quick test)
% FFT MTF sample size is set in computeMtf (S_64x64 by default)
% Nominal file: the .zmx file in the folder whose name does NOT start with "MC_"

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

%% Compute MTF of every MC file with and without compensation
freqGrid = linspace(0, maxFreq, nFreqPoints + 1)';
nFreq = numel(freqGrid);

fileNames = string({files.name})';
status = strings(nFiles, 1);
MTF_T_on = [];   MTF_S_on = [];    % nFiles x nFreq x nFields, compensated (file as saved)
MTF_T_off = [];  MTF_S_off = [];   % nFiles x nFreq x nFields, compensator at nominal
compSurfUsed = nan(nFiles, 1);     % compensator surface found in each MC file
compThickness = nan(nFiles, 1);    % compensated thickness of each MC file
fieldLabels = strings(0, 1);
nFields = 0;
diffT = [];
diffS = [];
diffSource = "";
paraxialEFL = NaN;
freqScale = 1;                     % Zemax frequency -> cyc/mm (1000/|EFL| for afocal files)
checkPrinted = false;

fprintf('Computing MTF for %d files (with and without compensation)...\n', nFiles);
wb = waitbar(0, 'Computing MTF...', 'Name', 'Compensator On/Off Analysis');

try
    TheSystem = TheApplication.PrimarySystem;

    %% Nominal file: afocal check, EFL, on-axis field, nominal MTF and compensator
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

    onAxisField = findOnAxisField(TheSystem);

    [nomSurf, nomSource] = findNominalCompensator(TheSystem);
    tNominal = getCompensator(TheSystem, nomSurf);
    nomComment = string(char(TheSystem.LDE.GetSurfaceAt(nomSurf).Comment));
    fprintf('Nominal compensator: thickness of surface %d %s = %.6f (%s)\n', ...
        nomSurf, nomComment, tNominal, nomSource);

    [nomT, nomS, ~, dT, dS] = computeMtf(TheSystem, maxFreq, freqGrid, freqScale, true, 0);
    if ~isempty(dT)
        diffT = dT;
        diffS = dS;
        diffSource = "OpticStudio FFT MTF";
    else
        diffT = analyticDiffLimit(TheSystem, freqGrid);
        diffS = diffT;
        diffSource = "Analytical (circular aperture, working F/#)";
    end

    fprintf('Nyquist frequency: %.4f cyc/mm (pixel pitch %.2f um)\n', maxFreq, pixelPitch_um);

    tStart = tic;

    for k = 1:nFiles

        filePath = fullfile(folderPath, files(k).name);

        try
            % Load optical system
            if ~TheSystem.LoadFile(filePath, false)
                error('OpticStudio could not load this file.');
            end

            % Compensator of THIS file (numbers differ from the nominal file)
            surfNo = findVariableThickness(TheSystem);
            compSurfUsed(k) = surfNo;
            t0 = getCompensator(TheSystem, surfNo);
            compThickness(k) = t0;

            % 1) Compensated: the file as saved by the Monte Carlo analysis
            [Ton, Son, labels] = computeMtf(TheSystem, maxFreq, freqGrid, freqScale, false, 0);

            % 2) Not compensated: same file, compensator back at its nominal value
            setCompensator(TheSystem, surfNo, tNominal);
            [Toff, Soff] = computeMtf(TheSystem, maxFreq, freqGrid, freqScale, false, 0);
            setCompensator(TheSystem, surfNo, t0);

            % First successful file defines the number of fields
            if nFields == 0
                nFields = size(Ton, 2);
                if size(nomT, 2) ~= nFields
                    error('Number of fields (%d) differs from nominal file (%d).', nFields, size(nomT, 2));
                end
                MTF_T_on = nan(nFiles, nFreq, nFields);   MTF_S_on = MTF_T_on;
                MTF_T_off = MTF_T_on;                     MTF_S_off = MTF_T_on;
                fieldLabels = labels;
            elseif size(Ton, 2) ~= nFields || size(Toff, 2) ~= nFields
                error('Number of fields differs from first file (%d).', nFields);
            end

            MTF_T_on(k, :, :) = Ton;    MTF_S_on(k, :, :) = Son;
            MTF_T_off(k, :, :) = Toff;  MTF_S_off(k, :, :) = Soff;

            if any(isnan([Ton(:); Son(:); Toff(:); Soff(:)]))
                error('MTF data contains invalid values.');
            end

            % Check output for the first file
            if ~checkPrinted
                comment = string(char(TheSystem.LDE.GetSurfaceAt(surfNo).Comment));
                fprintf('\nCheck on %s: compensator = thickness of surface %d %s\n', ...
                    files(k).name, surfNo, comment);
                fprintf('  compensated value %.6f, nominal value %.6f, travel %+.4f\n', ...
                    t0, tNominal, t0 - tNominal);
                fprintf('  on-axis mean MTF: with compensation %.4f, without %.4f\n\n', ...
                    mean([Ton(:, onAxisField); Son(:, onAxisField)]), ...
                    mean([Toff(:, onAxisField); Soff(:, onAxisField)]));
                checkPrinted = true;
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

clip = @(x) min(max(x, 0), 1);   % clip small interpolation overshoots
Ton = clip(MTF_T_on(valid, :, :));    Son = clip(MTF_S_on(valid, :, :));
Toff = clip(MTF_T_off(valid, :, :));  Soff = clip(MTF_S_off(valid, :, :));
nomT = clip(nomT);  nomS = clip(nomS);
diffT = clip(diffT);  diffS = clip(diffS);
travel = compThickness(valid) - tNominal;

for f = 1:nFields
    if fieldLabels(f) == ""
        fieldLabels(f) = sprintf('Field %d', f);
    end
end

fprintf('Valid files: %d / %d (computation time: %.1f s)\n', nValid, nFiles, toc(tStart));

%% Statistics
% Per file and frequency: worst MTF over all fields and both directions, and on-axis
worstOn  = min(min(Ton, [], 3), min(Son, [], 3));     % nValid x nFreq
worstOff = min(min(Toff, [], 3), min(Soff, [], 3));
axisOn   = min(Ton(:, :, onAxisField), Son(:, :, onAxisField));
axisOff  = min(Toff(:, :, onAxisField), Soff(:, :, onAxisField));
nomWorst = min(min(nomT, [], 2), min(nomS, [], 2))';  % 1 x nFreq
nomAxis  = min(nomT(:, onAxisField), nomS(:, onAxisField))';

% Yield curves of every field (for the per-field plots)
nYield = numel(yieldLevels);
ry = find(yieldLevels == reportYield, 1);
yT_on = zeros(nFreq, nFields);  yS_on = yT_on;  yT_off = yT_on;  yS_off = yT_on;
for f = 1:nFields
    yT_on(:, f)  = yieldValue(Ton(:, :, f),  reportYield)';
    yS_on(:, f)  = yieldValue(Son(:, :, f),  reportYield)';
    yT_off(:, f) = yieldValue(Toff(:, :, f), reportYield)';
    yS_off(:, f) = yieldValue(Soff(:, :, f), reportYield)';
end

% Yield and mean curves of the worst field and the on-axis field
curve = @(d) struct('mean', mean(d, 1), 'yield', cell2mat(arrayfun(@(y) yieldValue(d, y), yieldLevels(:), 'UniformOutput', false)));
cWorstOn = curve(worstOn);   cWorstOff = curve(worstOff);
cAxisOn  = curve(axisOn);    cAxisOff  = curve(axisOff);

% Summary at the key frequencies
keyFreqs = keyFreqFrac(:)' * maxFreq;
keyIdx = zeros(size(keyFreqs));
for q = 1:numel(keyFreqs)
    [~, keyIdx(q)] = min(abs(freqGrid - keyFreqs(q)));
end
keyNames = compose("%g x Nyquist", keyFreqFrac(:)');
recovery = @(yOn, yOff, nom) 100 * (yOn - yOff) ./ max(nom - yOff, eps);

nKey = numel(keyIdx);
sumWorst = zeros(nKey, 6);   % nominal, mean on, mean off, yield on, yield off, recovery
sumAxis  = zeros(nKey, 6);
for q = 1:nKey
    i = keyIdx(q);
    sumWorst(q, :) = [nomWorst(i), cWorstOn.mean(i), cWorstOff.mean(i), ...
        cWorstOn.yield(ry, i), cWorstOff.yield(ry, i), ...
        recovery(cWorstOn.yield(ry, i), cWorstOff.yield(ry, i), nomWorst(i))];
    sumAxis(q, :) = [nomAxis(i), cAxisOn.mean(i), cAxisOff.mean(i), ...
        cAxisOn.yield(ry, i), cAxisOff.yield(ry, i), ...
        recovery(cAxisOn.yield(ry, i), cAxisOff.yield(ry, i), nomAxis(i))];
end

fprintf('\n%d%% yield MTF, worst field (T/S):\n', reportYield);
fprintf('%-16s %9s %9s %9s %11s\n', 'Frequency', 'Nominal', 'With', 'Without', 'Recovery');
for q = 1:nKey
    fprintf('%-16s %9.4f %9.4f %9.4f %10.1f%%\n', keyNames(q), sumWorst(q, 1), sumWorst(q, 4), sumWorst(q, 5), sumWorst(q, 6));
end
fprintf('Compensator travel: mean %+.4f, std %.4f, max |travel| %.4f (lens units)\n', ...
    mean(travel), std(travel), max(abs(travel)));

%% Output folder: "KOMPANZATÖR KATKI ANALİZİ" (char(214) = 'Ö', char(304) = 'İ')
outFolder = fullfile(folderPath, ['KOMPANZAT' char(214) 'R KATKI ANAL' char(304) 'Z' char(304)]);

if ~isfolder(outFolder)
    mkdir(outFolder);
end

%% Plot pages (A4 portrait, two plots per page, every plot the same size)
plotPos = {[0.12 0.63 0.80 0.27], [0.12 0.20 0.80 0.27]};   % top / bottom plot
fieldColors = lines(nFields);
showDiffS = max(abs(diffT - diffS)) > 1e-6;
colOn  = [0 0.45 0.74];
colOff = [0.85 0.33 0.10];
k1 = keyIdx(1);                     % frequency used in the per-file plots
k1Name = keyNames(1);

% Page: comparison of the yield curves (worst field top, on-axis bottom)
cmpFig = newA4Page();
addPageTitle(cmpFig, sprintf('Kompanzatörlü ve kompanzatörsüz MTF (%d dosya)', nValid));
pairs = {cWorstOn, cWorstOff, nomWorst, 'Worst field (T/S)'; cAxisOn, cAxisOff, nomAxis, 'On-axis field (worse of T/S)'};
for p = 1:2
    ax = axes(cmpFig, 'Position', plotPos{p}, 'FontSize', 9);
    hold on;
    plot(freqGrid, diffT, 'k:', 'LineWidth', 1.2, 'DisplayName', 'Diff. limit');
    plot(freqGrid, pairs{p, 3}, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Nominal');
    plot(freqGrid, pairs{p, 1}.yield(ry, :), '-', 'Color', colOn, 'LineWidth', 2, ...
        'DisplayName', sprintf('With comp. %d%% yield', reportYield));
    plot(freqGrid, pairs{p, 2}.yield(ry, :), '-', 'Color', colOff, 'LineWidth', 2, ...
        'DisplayName', sprintf('Without comp. %d%% yield', reportYield));
    plot(freqGrid, pairs{p, 1}.mean, '--', 'Color', colOn, 'LineWidth', 1.2, 'DisplayName', 'With comp. mean');
    plot(freqGrid, pairs{p, 2}.mean, '--', 'Color', colOff, 'LineWidth', 1.2, 'DisplayName', 'Without comp. mean');
    hold off;
    xlim([0 maxFreq]);
    ylim([0 1]);
    xlabel('Spatial frequency (cyc/mm)');
    ylabel('Modulus of the OTF');
    title(pairs{p, 4});
    grid on;
    box on;
    pbaspect(ax, [1.6 1 1]);
    placeLegend(ax, plotPos{p}, 3);
end

% Page: per-file comparison at the first key frequency
fileFig = newA4Page();
addPageTitle(fileFig, sprintf('Dosya bazında karşılaştırma (%s)', k1Name));

ax = axes(fileFig, 'Position', plotPos{1}, 'FontSize', 9);
hold on;
lim = [0, max([worstOn(:, k1); worstOff(:, k1); nomWorst(k1)]) * 1.05];
plot(lim, lim, 'k-', 'LineWidth', 1, 'DisplayName', 'No change (y = x)');
plot(worstOff(:, k1), worstOn(:, k1), '.', 'Color', colOn, 'MarkerSize', 12, 'DisplayName', 'MC file');
plot(nomWorst(k1), nomWorst(k1), 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y', 'DisplayName', 'Nominal');
hold off;
xlim(lim);
ylim(lim);
xlabel('MTF without compensation (worst field)');
ylabel('MTF with compensation (worst field)');
title(sprintf('Each point is one MC file (%s)', k1Name));
grid on;
box on;
pbaspect(ax, [1.6 1 1]);
placeLegend(ax, plotPos{1}, 3);

ax = axes(fileFig, 'Position', plotPos{2}, 'FontSize', 9);
plot(travel, worstOn(:, k1) - worstOff(:, k1), '.', 'Color', colOn, 'MarkerSize', 12, 'DisplayName', 'MC file');
xlabel('Compensator travel (lens units)');
ylabel('MTF gained by compensation');
title(sprintf('MTF gain vs compensator travel (%s, worst field)', k1Name));
grid on;
box on;
pbaspect(ax, [1.6 1 1]);
placeLegend(ax, plotPos{2}, 1);

% Pages: report-yield MTF of all fields, with and without compensation
stateFigs = gobjects(2, 1);
stateData = {yT_on, yS_on, 'kompanzatörlü'; yT_off, yS_off, 'kompanzatörsüz'};
for st = 1:2
    fig = newA4Page();
    addPageTitle(fig, sprintf('%d%% yield MTF – %s (%d dosya)', reportYield, stateData{st, 3}, nValid));
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
            plot(freqGrid, stateData{st, 1}(:, f), '-', 'Color', fieldColors(f, :), ...
                'LineWidth', 1.2, 'DisplayName', sprintf('F%d T', f));
            plot(freqGrid, stateData{st, 2}(:, f), '--', 'Color', fieldColors(f, :), ...
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
    stateFigs(st) = fig;
end

%% Export Excel
outputFile = fullfile(outFolder, 'Compensator_OnOff_Results.xlsx');

% Remove old file so stale rows/sheets are not kept
if isfile(outputFile)
    delete(outputFile);
end

ByFile = table(fileNames, compSurfUsed, compThickness, compThickness - tNominal, status, ...
    'VariableNames', {'FileName', 'CompSurface', 'CompThickness', 'CompTravel', 'Status'});
for q = 1:nKey
    tag = sprintf('N%03d', round(100 * keyFreqFrac(q)));
    colOnVals = nan(nFiles, 1);   colOnVals(valid) = worstOn(:, keyIdx(q));
    colOffVals = nan(nFiles, 1);  colOffVals(valid) = worstOff(:, keyIdx(q));
    ByFile.(['WorstMTF_with_' tag]) = colOnVals;
    ByFile.(['WorstMTF_without_' tag]) = colOffVals;
    ByFile.(['Gain_' tag]) = colOnVals - colOffVals;
end
writetable(ByFile, outputFile, 'Sheet', 'By File');

Summary = table([keyNames'; keyNames'], [repmat("Worst field", nKey, 1); repmat("On-axis", nKey, 1)], ...
    [sumWorst(:, 1); sumAxis(:, 1)], [sumWorst(:, 2); sumAxis(:, 2)], [sumWorst(:, 3); sumAxis(:, 3)], ...
    [sumWorst(:, 4); sumAxis(:, 4)], [sumWorst(:, 5); sumAxis(:, 5)], [sumWorst(:, 6); sumAxis(:, 6)], ...
    'VariableNames', {'Frequency', 'Field', 'Nominal', 'Mean_with', 'Mean_without', ...
    sprintf('Yield%d_with', reportYield), sprintf('Yield%d_without', reportYield), 'Recovery_pct'});
writetable(Summary, outputFile, 'Sheet', 'Summary');

Curves = table(freqGrid, nomWorst', nomAxis', diffT, 'VariableNames', {'Frequency_cyc_mm', 'Nominal_worst', 'Nominal_axis', 'DiffLimit_T'});
for y = 1:nYield
    Curves.(sprintf('Worst_Yield%d_with', yieldLevels(y)))    = cWorstOn.yield(y, :)';
    Curves.(sprintf('Worst_Yield%d_without', yieldLevels(y))) = cWorstOff.yield(y, :)';
    Curves.(sprintf('Axis_Yield%d_with', yieldLevels(y)))     = cAxisOn.yield(y, :)';
    Curves.(sprintf('Axis_Yield%d_without', yieldLevels(y)))  = cAxisOff.yield(y, :)';
end
Curves.Worst_Mean_with = cWorstOn.mean';   Curves.Worst_Mean_without = cWorstOff.mean';
Curves.Axis_Mean_with = cAxisOn.mean';     Curves.Axis_Mean_without = cAxisOff.mean';
writetable(Curves, outputFile, 'Sheet', 'Yield Curves');

FieldInfo = table((1:nFields)', fieldLabels, 'VariableNames', {'Field', 'Description'});
writetable(FieldInfo, outputFile, 'Sheet', 'Field Info');

% Raw per-file MTF curves for further analysis
save(fullfile(outFolder, 'Compensator_OnOff_RawData.mat'), 'freqGrid', 'MTF_T_on', 'MTF_S_on', ...
    'MTF_T_off', 'MTF_S_off', 'nomT', 'nomS', 'diffT', 'diffS', 'fileNames', 'status', ...
    'compSurfUsed', 'compThickness', 'tNominal', 'fieldLabels', 'pixelPitch_um', 'maxFreq', ...
    'nominalName', 'nominalAfocal', 'paraxialEFL');

%% PDF report (A4): summary, explanation pages, plot pages
pdfFile = fullfile(outFolder, 'Compensator_OnOff_Report.pdf');

% Remove old report files (single and per-page)
oldPages = dir(fullfile(outFolder, 'Compensator_OnOff_Report_page*.pdf'));
for i = 1:numel(oldPages)
    delete(fullfile(outFolder, oldPages(i).name));
end
if isfile(pdfFile)
    delete(pdfFile);
end

% If the old report could not be deleted (e.g. open in a PDF viewer),
% save the new one under a time-stamped name
if isfile(pdfFile)
    pdfFile = fullfile(outFolder, ['Compensator_OnOff_Report_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.pdf']);
    fprintf('\nCompensator_OnOff_Report.pdf is open in another program; new report: %s\n', pdfFile);
end

if nominalAfocal
    afocalText = sprintf('Evet – MTF cyc/mrad olarak hesaplandı, 1000/|EFL| ile cyc/mm''ye çevrildi (|EFL| = %.4f mm)', paraxialEFL);
else
    afocalText = 'Hayır';
end
compSurfText = char(strjoin(string(unique(compSurfUsed(valid))'), ', '));

% --- Page: summary -------------------------------------------------------
summaryLines = {
    sprintf('Tarih:                %s', char(datetime('now', 'Format', 'dd.MM.yyyy HH:mm')))
    sprintf('Klasör:               %s', folderPath)
    sprintf('Nominal dosya:        %s', nominalName)
    sprintf('MC dosyaları:         %s*  (geçerli: %d / %d)', filePrefix, nValid, nFiles)
    sprintf('Afokal mod:           %s', afocalText)
    sprintf('Kompanzatör:          nominal dosyada yüzey %d, nominal kalınlık %.6f', nomSurf, tNominal)
    sprintf('                      MC dosyalarında yüzey %s (tek variable kalınlık)', compSurfText)
    sprintf('Piksel aralığı:       %.2f µm   Nyquist: %.2f cyc/mm', pixelPitch_um, maxFreq)
    sprintf('Diffraction limit:    %s', char(diffSource))
    ''
    sprintf('%d%% YIELD MTF – EN KÖTÜ FIELD (T/S)', reportYield)
    };

header = sprintf('%-15s%10s%12s%12s%12s', 'Frekans', 'Nominal', 'Komp. var', 'Komp. yok', 'Geri kaz.');
summaryLines = [summaryLines; {header; repmat('-', 1, numel(header))}];
for q = 1:nKey
    summaryLines{end + 1, 1} = sprintf('%-15s%10.4f%12.4f%12.4f%11.1f%%', keyNames(q), ...
        sumWorst(q, 1), sumWorst(q, 4), sumWorst(q, 5), sumWorst(q, 6)); %#ok<SAGROW>
end
summaryLines = [summaryLines; {''; sprintf('%d%% YIELD MTF – ON-AXIS FIELD (T/S''nin kötüsü)', reportYield); header; repmat('-', 1, numel(header))}];
for q = 1:nKey
    summaryLines{end + 1, 1} = sprintf('%-15s%10.4f%12.4f%12.4f%11.1f%%', keyNames(q), ...
        sumAxis(q, 1), sumAxis(q, 4), sumAxis(q, 5), sumAxis(q, 6)); %#ok<SAGROW>
end
summaryLines = [summaryLines; {
    ''
    'Geri kaz. = (komp. var - komp. yok) / (nominal - komp. yok): toleransların zararının'
    '            kompanzatörle geri kazanılan yüzdesi.'
    ''
    'KOMPANZATÖR HAREKETİ (MC değeri - nominal değer, lens birimi)'
    sprintf('  ortalama %+.4f   std %.4f   en büyük |hareket| %.4f', mean(travel), std(travel), max(abs(travel)))
    ''
    'FIELDLER'
    }];
summaryLines = [summaryLines; cellstr(compose("  F%d: %s", (1:nFields)', fieldLabels))];

failed = fileNames(~valid);
if ~isempty(failed)
    summaryLines = [summaryLines; {''; sprintf('Hesaplanamayan dosyalar (%d):', numel(failed))}];
    shown = failed(1:min(end, 8));
    summaryLines = [summaryLines; cellstr("  " + shown)];
    if numel(failed) > numel(shown)
        summaryLines{end + 1, 1} = sprintf('  ... ve %d dosya daha (Excel, By File sayfası)', numel(failed) - numel(shown));
    end
end

summaryFig = textPage('KOMPANZATÖRÜN KATKISI – MTF ANALİZ RAPORU', summaryLines, 'Courier New', 8.5);

% --- Explanation pages ---------------------------------------------------
explain1 = {
    'AMAÇ'
    ['Bu analiz, kompanzatörün toleranslardan gelen hatanın ne kadarını düzelttiğini gösterir. Aynı Monte ' ...
     'Carlo dosyalarında MTF iki durumda hesaplanır: kompanzatör ayarlıyken ve kompanzatör hiç ' ...
     'oynatılmamış (nominal konumunda) iken. Monte Carlo yalnızca afokal grubun toleranslarıyla ' ...
     'koşulduysa sonuç, afokal gruptan gelen hatanın objektifteki kompanzatörle ne kadar düzeltilebildiğini ' ...
     'doğrudan gösterir. Bu karşılaştırma odak hatasıyla sınırlı değildir; toleransların ürettiği tüm ' ...
     'aberasyonlar (odak, koma, astigmatizma, küresel) içindedir.']
    ''
    'ADIM 1 – NOMİNAL DOSYA'
    sprintf(['Klasörde adı MC_ ile başlamayan dosya nominal tasarımdır (bu analizde: %s). Bu dosyadan ' ...
     'afokal mod, on-axis field, nominal MTF ve kompanzatörün nominal kalınlığı (%.6f) okunur.'], nominalName, tNominal)
    ''
    'ADIM 2 – İKİ DURUMDA MTF'
    ['Her MC dosyası, Zemax Monte Carlo analizinin rastgele bozduğu ve kompanzatörünü optimize ettiği ' ...
     'bir sistemdir. Her dosyada:']
    '•  Kompanzatörlü: dosyanın kendi haliyle MTF hesaplanır.'
    ['•  Kompanzatörsüz: kompanzatör kalınlığı nominal değere geri alınır ve MTF yeniden hesaplanır. ' ...
     'Bu, aynı bozulmuş sistemin kompanzatör hiç ayarlanmamış halidir. Optimizasyon yapılmaz.']
    ['Kompanzatör her MC dosyasında kalınlığı variable olan tek yüzey olarak bulunur; toleranslama MC ' ...
     'dosyalarına coordinate break eklediği için yüzey numaraları nominal dosyadan farklıdır. İki durum aynı ' ...
     'dosyadan çıktığı için birebir eşleşir; ikinci bir Monte Carlo koşusuna gerek yoktur. Diskteki ' ...
     'dosyalar değiştirilmez.']
    sprintf(['MTF: polikromatik FFT MTF, tüm dalga boyları, tüm fieldler, T ve S ayrı, 0''dan Nyquist''e ' ...
     '(%.2f cyc/mm) kadar.'], maxFreq)
    ''
    'ADIM 3 – İSTATİSTİK'
    ['•  En kötü field: her dosyada, her frekansta tüm fieldlerin T ve S MTF''inin en düşüğü. Sistemin en ' ...
     'zayıf noktasını temsil eder.']
    '•  On-axis: on-axis fieldin T ve S MTF''inin düşüğü.'
    sprintf(['•  Yield %%%d MTF: dosyaların en az %%%d''inin sağladığı MTF. Her frekansta ayrı hesaplanır; ' ...
     'tek bir sistemin eğrisi değildir.'], reportYield, reportYield)
    ['•  Geri kazanım = (kompanzatörlü − kompanzatörsüz) / (nominal − kompanzatörsüz) × 100. Toleransların ' ...
     'yaptığı zararın yüzde kaçının kompanzatörle geri kazanıldığını gösterir.']
    };

explain2 = {
    'SONUÇLAR NASIL YORUMLANIR?'
    ['•  Kompanzatörsüz yield nominalden çok düşük, kompanzatörlü yield nominale yakın (geri kazanım ' ...
     'yüksek): toleransların hatası ağırlıklı olarak odak türündendir ve kompanzatör bunu düzeltiyor.']
    ['•  Her iki yield de düşük (geri kazanım düşük): hata ağırlıklı olarak kompanzatörün düzeltemediği ' ...
     'aberasyonlardan (koma, astigmatizma, küresel) geliyor. Toleransları sıkılaştırmak ya da başka bir ' ...
     'ayar (yanal ayar, clocking) düşünmek gerekir.']
    ['•  Her iki yield de nominale yakın: toleranslar zaten zararsız; kompanzatöre bu toleranslar için ' ...
     'ihtiyaç yok.']
    ['•  Kompanzatörlü yield ile nominal arasındaki fark, kompanzasyondan sonra kalan kayıptır ve ' ...
     'kompanzatörün düzeltemediği kısmı gösterir.']
    ['Geri kazanım, kompanzatörsüz yield nominale çok yakınsa (zarar çok küçükse) anlamsızlaşır; o durumda ' ...
     'yield değerlerinin kendisine bakılmalıdır.']
    ''
    'GRAFİKLER NASIL OKUNUR?'
    sprintf(['•  Kompanzatörlü ve kompanzatörsüz MTF: Üstte en kötü field, altta on-axis. Siyah düz çizgi ' ...
     'nominal, siyah noktalı çizgi diffraction limit. Mavi kompanzatörlü, turuncu kompanzatörsüz; kalın düz ' ...
     'çizgi %%%d yield, kesikli çizgi ortalama. Mavi ile turuncu arasındaki fark kompanzatörün katkısıdır; ' ...
     'mavi ile siyah arasındaki fark kalan kayıptır.'], reportYield)
    ['•  Dosya bazında karşılaştırma (üst): Her nokta bir MC dosyası. Yatay eksen kompanzatörsüz, dikey eksen ' ...
     'kompanzatörlü MTF. Siyah çizginin (y = x) üstündeki noktalar kompanzatörün iyileştirdiği dosyalardır; ' ...
     'çizgiden ne kadar uzaksa kazanç o kadar büyüktür. Yıldız nominal sistemdir.']
    ['•  Kazanç – kompanzatör hareketi (alt): Kompanzatör ne kadar çok hareket etmek zorunda kaldıysa kazanç ' ...
     'o kadar büyükse, kompanzatörün düzelttiği şey odak türü hatadır (afokallik bozulması). Noktalar ' ...
     'hareketten bağımsız dağılmışsa kazancın kaynağı başka bir etkidir.']
    sprintf(['•  Son iki sayfa: %%%d yield MTF eğrileri, kompanzatörlü ve kompanzatörsüz. Üstte tüm fieldler, ' ...
     'altta on-axis. Düz çizgi T, kesikli çizgi S, siyah diffraction limit.'], reportYield)
    ''
    'BU ANALİZİN GÖSTERMEDİKLERİ'
    ['•  Sonuç, Monte Carlo''ya konan toleranslara bağlıdır. Yalnızca afokal grubun toleranslarıyla koşulduysa ' ...
     'sadece afokal grup hakkında konuşur; objektifin toleransları ve gruplar arası hizalama ayrıca ele alınmalıdır.']
    ['•  Kompanzatör hareketinin mekanik stroğa sığıp sığmadığı burada değerlendirilmez; hareket değerleri ' ...
     'özet sayfasında ve Excel''de verilir.']
    '•  Sıcaklık değişimi ve kompanzatörün ayar hatası bu analizde yoktur.'
    };

explainFig1 = textPage('1. ANALİZ YÖNTEMİ', explain1, 'Arial', 10);
explainFig2 = textPage('2. SONUÇLARIN YORUMU', explain2, 'Arial', 10);

% --- Assemble pages: frame, header and page numbers ----------------------
plotFigs = [cmpFig; fileFig; stateFigs];
plotNames = ["Compensator_OnOff_Yield", "Compensator_OnOff_PerFile", ...
    sprintf("MTF_Yield%d_with_comp", reportYield), sprintf("MTF_Yield%d_without_comp", reportYield)];
allFigs = [summaryFig; explainFig1; explainFig2; plotFigs];
nPages = numel(allFigs);
headerText = sprintf('Kompanzatörün Katkısı – MTF Analizi  |  %s', nominalName);

for p = 1:nPages
    addPageFrame(allFigs(p), headerText, p, nPages);
end

% PNG of the plot pages
for p = 1:numel(plotFigs)
    exportPng(plotFigs(p), fullfile(outFolder, plotNames(p) + ".png"));
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
        print(allFigs(p), fullfile(outFolder, sprintf('Compensator_OnOff_Report_page%d.pdf', p)), '-dpdf');
    end
    fprintf('Report saved as separate pages (Compensator_OnOff_Report_page*.pdf).\n');
end

close(allFigs);

fprintf('\nResults saved to:\n%s\n', outFolder);


function [surfNo, source] = findNominalCompensator(TheSystem)
% Compensator of the NOMINAL file: the only variable thickness; if there is
% not exactly one, the COMP operand of the Tolerance Data Editor.
lde = TheSystem.LDE;
found = [];
for i = 0:lde.NumberOfSurfaces - 1
    solveType = string(lde.GetSurfaceAt(i).ThicknessCell.GetSolveData().Type.ToString());
    if solveType == "Variable"
        found(end + 1) = i; %#ok<AGROW>
    end
end
if isscalar(found)
    surfNo = found;
    source = "only variable thickness";
    return;
end

comp = [];
try
    tde = TheSystem.TDE;
    for i = 1:tde.NumberOfOperands
        op = tde.GetOperandAt(i);
        if strcmpi(string(op.Type.ToString()), "COMP") && double(op.Param2) == 0
            comp(end + 1) = double(op.Param1); %#ok<AGROW>
        end
    end
catch
end
if isscalar(comp)
    surfNo = comp;
    source = "COMP operand";
    return;
end

error(['Compensator of the nominal file could not be identified: %d variable thicknesses, ' ...
    '%d thickness COMP operands. Make the compensator the only variable thickness.'], numel(found), numel(comp));
end

function [T, S, labels, dT, dS] = computeMtf(TheSystem, maxFreq, freqGrid, freqScale, showDiffLimit, fieldNo)
% Runs a new FFT MTF analysis on the current system and returns tangential
% and sagittal MTF (nFreq x nFields) on the common frequency grid, plus the
% diffraction limit if requested and returned by OpticStudio. A new analysis
% is opened for every call so that no earlier result can be reused.
% fieldNo = 0: all fields; fieldNo > 0: only that field (one column).
mtf = TheSystem.Analyses.New_FftMtf();
closeMtf = onCleanup(@() mtf.Close());

settings = mtf.GetSettings();
settings.MaximumFrequency = maxFreq / freqScale;   % cyc/mrad if afocal
settings.SampleSize = ZOSAPI.Analysis.SampleSizes.S_64x64;
settings.Wavelength.UseAllWavelengths();
if fieldNo > 0
    settings.Field.SetFieldNumber(fieldNo);
else
    settings.Field.UseAllFields();
end
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

if fieldNo > 0 && numel(fieldIdx) ~= 1
    error('Single-field MTF (field %d) returned %d data series instead of 1.', fieldNo, numel(fieldIdx));
end

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

function surfNo = findVariableThickness(TheSystem)
% Returns the surface whose thickness is a variable: the compensator of a
% Monte Carlo file. Exactly one variable thickness is expected.
lde = TheSystem.LDE;
found = [];
for i = 0:lde.NumberOfSurfaces - 1
    solveType = string(lde.GetSurfaceAt(i).ThicknessCell.GetSolveData().Type.ToString());
    if solveType == "Variable"
        found(end + 1) = i; %#ok<AGROW>
    end
end

if isempty(found)
    error('No variable thickness found: the compensator cannot be identified.');
elseif numel(found) > 1
    error('%d variable thicknesses found (surfaces %s): exactly one (the compensator) is expected.', ...
        numel(found), strjoin(string(found), ', '));
end
surfNo = found;
end

function v = getCompensator(TheSystem, surfNo)
v = TheSystem.LDE.GetSurfaceAt(surfNo).Thickness;
end

function setCompensator(TheSystem, surfNo, value)
% Sets the compensator thickness and checks that the system really took it.
% The .NET object is stored in a variable first: MATLAB does not apply a
% property assignment written on a method call chain (obj.Method(x).Prop = v).
surf = TheSystem.LDE.GetSurfaceAt(surfNo);
surf.Thickness = value;

actual = TheSystem.LDE.GetSurfaceAt(surfNo).Thickness;
if abs(actual - value) > 1e-9
    error('Compensator thickness of surface %d could not be changed (wanted %.6f, system has %.6f).', ...
        surfNo, value, actual);
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
% Arial is set explicitly: the default font has no Turkish glyphs in the PDF
annotation(fig, 'textbox', [0.06 0.915 0.88 0.04], 'String', titleText, ...
    'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'Interpreter', 'none');
end

function fig = textPage(titleText, lines, fontName, fontSize)
% A4 page with a title and wrapped text (one cell per paragraph/line)
fig = newA4Page();
annotation(fig, 'textbox', [0.08 0.89 0.84 0.05], 'String', titleText, ...
    'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
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
text(ax, 0.06, 0.970, headerText, 'FontName', 'Arial', 'FontSize', 8, 'Color', [0.4 0.4 0.4], ...
    'VerticalAlignment', 'bottom', 'Interpreter', 'none');
text(ax, 0.94, 0.030, sprintf('Sayfa %d / %d', pageNo, nPages), 'FontName', 'Arial', 'FontSize', 8, ...
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
