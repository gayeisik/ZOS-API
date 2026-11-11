%% MTF YIELD ANALİZİ (ZOS-API, MATLAB)
% NE YAPAR
%   Monte Carlo dosyalarının (adı "MC_T" ile başlayan .zmx) her biri için MTF
%   hesaplar ve tüm dosyalar üzerinden ortalama MTF ile yield MTF eğrilerini
%   çıkarır. Zemax'in Monte Carlo yield mantığıyla aynıdır.
%
% NASIL ÇALIŞIR
%   1. Klasör seçilir. Adı "MC_" ile başlamayan .zmx dosyası nominal tasarım
%      kabul edilir (birden fazlaysa listeden seçilir).
%   2. Nominal dosyada "Afocal Image Space" modu kontrol edilir. Afokalse EFFL
%      okunur ve |EFL| kaydedilir.
%   3. Her MC_T dosyası açılır ve polikromatik FFT MTF hesaplanır: tüm dalga
%      boyları, tüm fieldler, tanjansiyel (T) ve sajital (S) ayrı ayrı, 0'dan
%      Nyquist frekansına kadar. Dosyalar değiştirilmez.
%   4. Diffraction limit eğrisi ilk dosyadan alınır.
%
% HESAPLAR
%   Nyquist frekansı = 1000 / (2 × piksel aralığı[µm])  [cyc/mm]  -> maksimum frekans
%   Afokal dosyalarda Zemax MTF'i cyc/mrad verir. Frekans ekseni
%       ν[cyc/mm] = ν[cyc/mrad] × 1000 / |EFL|
%   ile cyc/mm'ye çevrilir; MTF değerleri değişmez.
%   Ortalama MTF  : her frekansta dosyaların MTF ortalaması.
%   Yield %X MTF  : her frekansta, dosyaların en az %X'inin sağladığı (eşit veya
%                   üstünde kaldığı) MTF değeri. Eğri her frekansta ayrı
%                   hesaplanır; tek bir sistemin eğrisi değildir.
%
% AYARLAR (User settings)
%   filePrefix, pixelPitch_um, nFreqPoints, yieldLevels
%
% ÇIKTILAR (seçilen klasörde "MTF ANALİZİ" klasörüne)
%   MTF_Report.pdf   : A4; özet, yöntem ve yorum açıklamaları, ortalama MTF ve her
%                      yield seviyesi için bir sayfa (üstte tüm fieldler, altta
%                      on-axis, siyah: diff. limit)
%   MTF_Mean.png, MTF_Yield*.png : aynı grafikler
%   MTF_Results.xlsx : her field için frekansa göre ortalama ve yield MTF,
%                      field bilgisi, diffraction limit, dosya durumları
%   MTF_RawData.mat  : her dosyanın ham MTF eğrileri

clc;
clear;
close all;

%% User settings
filePrefix    = 'MC_T';           % only files starting with this are used
pixelPitch_um = 15;               % detector pixel pitch (um)
nFreqPoints   = 100;              % number of frequency steps from 0 to maxFreq
yieldLevels   = [50 80 90 95];    % yield levels (%)
% FFT MTF sample size is set in the "Compute MTF" section (S_64x64 by default)
% Nominal file: the .zmx file in the folder whose name does NOT start with "MC_"

% Nyquist frequency of the detector, used as the maximum MTF frequency
maxFreq = 1000 / (2 * pixelPitch_um);   % cyc/mm

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
paraxialEFL = NaN;   % |EFL| of the nominal file, used for afocal frequency scaling (mm)
freqScale = 1;       % Zemax frequency -> cyc/mm (1000/|EFL| for afocal files)

fprintf('Computing MTF for %d files...\n', nFiles);
wb = waitbar(0, 'Computing MTF...', 'Name', 'MTF Analysis');

try
    TheSystem = TheApplication.PrimarySystem;

    %% Nominal file: afocal check and EFL
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
        fprintf('Nominal file %s is in afocal mode.\n', nominalName);
        fprintf('|EFL|: %.4f mm, MTF frequency cyc/mrad -> cyc/mm (x 1000/|EFL|)\n', paraxialEFL);

        % Afocal MTF is in cyc/mrad: nu_mm = nu_mrad * 1000 / |EFL|
        freqScale = 1000 / paraxialEFL;
    else
        fprintf('Nominal file %s is not in afocal mode.\n', nominalName);
    end

    fprintf('Nyquist frequency: %.4f cyc/mm (pixel pitch %.2f um)\n', maxFreq, pixelPitch_um);

    tStart = tic;

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
            settings.MaximumFrequency = maxFreq / freqScale;   % cyc/mrad if afocal
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
                [MTF_T(k, :, f), MTF_S(k, :, f)] = readSeries(data, freqGrid, freqScale);

                if fieldLabels(f) == ""
                    fieldLabels(f) = strtrim(string(char(data.Description)));
                end
            end

            % Diffraction limit (taken once, from the first successful file)
            if isempty(diffT)
                diffIdx = find(isDiff, 1);
                if ~isempty(diffIdx)
                    [diffT, diffS] = readSeries(dataSeries(diffIdx), freqGrid, freqScale);
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

fprintf('Valid files: %d / %d (MTF computation time: %.1f s)\n', nValid, nFiles, toc(tStart));

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

%% Plot pages (A4 portrait, two plots per page, every plot the same size)
% Top: all fields, bottom: on-axis field
% Solid: tangential, dashed: sagittal, black: diffraction limit
plotPos = {[0.12 0.63 0.80 0.27], [0.12 0.20 0.80 0.27]};   % top / bottom plot
fieldColors = lines(nFields);
showDiffS = max(abs(diffT - diffS)) > 1e-6;

% Level 1: mean MTF, levels 2..end: yield MTF
plotT = cat(2, reshape(meanT, nFreq, 1, nFields), yieldT);
plotS = cat(2, reshape(meanS, nFreq, 1, nFields), yieldS);
levelNames = ["Ortalama MTF", compose("%d%% yield MTF", yieldLevels)];
levelFiles = ["MTF_Mean", compose("MTF_Yield%d", yieldLevels)];
nLevels = numel(levelNames);
levelFigs = gobjects(nLevels, 1);

for y = 1:nLevels
    fig = newA4Page();
    addPageTitle(fig, sprintf('%s (%d dosya)', levelNames(y), nValid));

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
        pbaspect(ax, [1.6 1 1]);
        placeLegend(ax, plotPos{p}, 6);
    end

    levelFigs(y) = fig;
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

FieldInfo = table((1:nFields)', fieldLabels, 'VariableNames', {'Field', 'Description'});
writetable(FieldInfo, outputFile, 'Sheet', 'Field Info');

DiffTable = table(freqGrid, diffT, diffS, ...
    'VariableNames', {'Frequency_cyc_mm', 'DiffLimit_T', 'DiffLimit_S'});
writetable(DiffTable, outputFile, 'Sheet', 'Diffraction Limit');

Files = table(fileNames, status, 'VariableNames', {'FileName', 'Status'});
writetable(Files, outputFile, 'Sheet', 'Files');

% Raw per-file MTF curves for further analysis
save(fullfile(outFolder, 'MTF_RawData.mat'), 'freqGrid', 'MTF_T', 'MTF_S', ...
    'diffT', 'diffS', 'fileNames', 'status', 'fieldLabels', 'yieldLevels', ...
    'pixelPitch_um', 'maxFreq', 'nominalName', 'nominalAfocal', 'paraxialEFL');

%% PDF report (A4): summary, explanation pages, mean MTF page, one page per yield level
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

yieldListText = strjoin(string(yieldLevels), ', ');

if nominalAfocal
    afocalText = sprintf('Evet – MTF cyc/mrad olarak hesaplandı, 1000/|EFL| ile cyc/mm''ye çevrildi (|EFL| = %.4f mm)', paraxialEFL);
    afocalPara = sprintf(['Nominal dosya afokal moddadır. Bu modda Zemax MTF''i açısal frekansta (cyc/mrad) ' ...
        'verir. Nominal dosyadan |EFL| = %.4f mm okunmuş ve frekanslar  ν[cyc/mm] = ν[cyc/mrad] × 1000 / |EFL|  ' ...
        'ile cyc/mm''ye çevrilmiştir (katsayı %.4f). MTF değerleri değişmez, sadece frekans ekseni ölçeklenir. ' ...
        'Bu, sistemin sonuna odak uzaklığı |EFL| olan ideal bir lens koymakla eşdeğerdir.'], paraxialEFL, freqScale);
else
    afocalText = 'Hayır';
    afocalPara = 'Nominal dosya afokal modda değildir; MTF doğrudan cyc/mm olarak hesaplanır.';
end

% --- Page: summary -------------------------------------------------------
summaryLines = {
    sprintf('Tarih:               %s', char(datetime('now', 'Format', 'dd.MM.yyyy HH:mm')))
    sprintf('Klasör:              %s', folderPath)
    sprintf('Nominal dosya:       %s', nominalName)
    sprintf('MC dosyaları:        %s*  (geçerli: %d / %d)', filePrefix, nValid, nFiles)
    sprintf('Afokal mod:          %s', afocalText)
    sprintf('Piksel aralığı:      %.2f µm', pixelPitch_um)
    sprintf('Nyquist frekansı:    %.2f cyc/mm  (maksimum MTF frekansı)', maxFreq)
    sprintf('Yield seviyeleri:    %s %%', yieldListText)
    'MTF tipi:            Polikromatik FFT MTF (tüm dalga boyları, tüm fieldler)'
    sprintf('On-axis field:       F%d (%s)', onAxisField, char(fieldLabels(onAxisField)))
    sprintf('Diffraction limit:   %s', char(diffSource))
    ''
    'RAPORUN İÇERİĞİ'
    '  Sayfa 2-3 : Analiz yöntemi ve grafiklerin yorumu'
    sprintf('  Sayfa 4   : Ortalama MTF')
    sprintf('  Sayfa 5-%d : Yield MTF (%s %%)', 4 + numel(yieldLevels), yieldListText)
    ''
    'FIELDLER (grafiklerde düz çizgi = T, kesikli çizgi = S)'
    };

fieldLines = compose("  F%d: %s", (1:nFields)', fieldLabels);
summaryLines = [summaryLines; cellstr(fieldLines)];

failed = fileNames(~valid);
if ~isempty(failed)
    summaryLines = [summaryLines; {''; sprintf('Hesaplanamayan dosyalar (%d):', numel(failed))}];
    shown = failed(1:min(end, 10));
    summaryLines = [summaryLines; cellstr("  " + shown)];
    if numel(failed) > numel(shown)
        summaryLines{end + 1, 1} = sprintf('  ... ve %d dosya daha (Excel, Files sayfası)', numel(failed) - numel(shown));
    end
end

summaryFig = textPage('MTF YIELD ANALİZ RAPORU', summaryLines, 'Courier New', 8.5);

% --- Explanation pages ---------------------------------------------------
explain1 = {
    'AMAÇ'
    ['Bu analiz, Monte Carlo (MC) tolerans analizinin kaydettiği MC_T dosyalarının MTF performansını ' ...
     'istatistiksel olarak değerlendirir. Her MC dosyası, toleranslar uygulanmış ve kompanzatörü ' ...
     'optimize edilmiş olası bir üretim sistemini temsil eder. Amaç, üretilecek sistemlerin belirli ' ...
     'bir oranının (yield) en az hangi MTF değerini sağlayacağını görmektir.']
    ''
    'ADIM 1 – NOMİNAL DOSYA'
    sprintf(['Klasörde adı MC_ ile başlamayan dosya nominal tasarım kabul edilir (bu analizde: %s). ' ...
     'Bu dosyada sistemin "Afocal Image Space" modunda olup olmadığına bakılır.'], nominalName)
    afocalPara
    ''
    'ADIM 2 – MTF HESABI'
    sprintf(['Her MC_T dosyası açılır ve polikromatik FFT MTF hesaplanır: tüm dalga boyları, ' ...
     'tüm fieldler, tanjansiyel (T) ve sajital (S) ayrı ayrı, 0''dan Nyquist frekansına ' ...
     '(%.2f cyc/mm) kadar. Diskteki dosyalar değiştirilmez. Açılamayan veya MTF''i hesaplanamayan ' ...
     'dosyalar özet sayfasında listelenir ve istatistiğe katılmaz.'], maxFreq)
    sprintf(['Diffraction limit eğrisi (aberasyonsuz ideal sistemin teorik üst sınırı) ilk dosyadan ' ...
     'bir kez alınır; MC dosyaları arasında pratikte değişmez. Kaynak: %s.'], char(diffSource))
    ''
    'ADIM 3 – İSTATİSTİK: ORTALAMA VE YIELD MTF'
    'Tüm geçerli dosyaların MTF''leri her frekansta ayrı ayrı istatistiğe sokulur:'
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
     'üzerindeki detaylar dedektör tarafından doğru örneklenemez (aliasing); bu yüzden grafikler ' ...
     'Nyquist''e kadar çizilir.'], pixelPitch_um, maxFreq)
    ['•  Nyquist''teki MTF (grafiğin en sağı): Sistemin, dedektörün çözebileceği en ince detayı ne kadar ' ...
     'kontrastla ilettiğini gösterir. Küçük hedeflerin ve ince çizgilerin seçilebilirliği buna bağlıdır.']
    ['•  Nyquist/2''deki MTF (grafiğin ortası): Orta frekanslardaki kontrastı, yani görüntünün genel ' ...
     'keskinliğini ve algılanan netliğini temsil eder.']
    ''
    'GRAFİK SAYFALARI NASIL OKUNUR?'
    ['İlk grafik sayfası ortalama MTF''i, sonraki sayfalar her yield seviyesinin MTF eğrilerini gösterir. ' ...
     'Her sayfada üstteki grafik tüm fieldleri, alttaki grafik yalnızca on-axis fieldi içerir.']
    '•  Yatay eksen: uzaysal frekans (cyc/mm), 0''dan Nyquist''e kadar. Dikey eksen: MTF (0-1).'
    ['•  Her renk bir fielddir (F1, F2, ...; liste özet sayfasındadır). Düz çizgi tanjansiyel (T), ' ...
     'kesikli çizgi sajital (S) yöndür.']
    ['•  Siyah çizgi diffraction limittir. Hiçbir eğri bunun üzerine çıkamaz; eğriler bu çizgiye ne kadar ' ...
     'yakınsa performans o kadar iyidir.']
    ''
    'SONUÇLAR NASIL YORUMLANIR?'
    ['•  %50 yield medyan sistemi temsil eder: sistemlerin yarısı bundan iyi, yarısı kötüdür. Yüksek ' ...
     'yield (%90, %95) eğrileri üretimde neredeyse tüm sistemlerin sağlayacağı alt sınırı gösterir ve ' ...
     'her zaman daha alttadır. Spesifikasyon karşılaştırması genellikle yüksek yield eğrisiyle yapılır.']
    ['•  Ortalama MTF ile yüksek yield eğrisi arasındaki fark, toleransların performansı ne kadar ' ...
     'dağıttığını gösterir. Fark büyükse tasarım toleranslara hassastır; toleransları sıkılaştırmak veya ' ...
     'kompanzasyonu iyileştirmek gerekebilir.']
    ['•  Fieldler arasındaki fark: Kenar fieldlerin eğrileri on-axis''e göre belirgin düşükse performansı ' ...
     'kenar fieldler sınırlıyor demektir. Aynı fieldde T ve S eğrilerinin ayrılması astigmatizma / ' ...
     'field curvature etkisini gösterir.']
    ['•  Bir MTF şartını kontrol etmek için (örneğin "%90 yield''de Nyquist''te MTF ≥ 0.15") ilgili yield ' ...
     'sayfasında o frekanstaki en düşük eğriye bakılır. Sayısal değerler MTF_Results.xlsx dosyasında ' ...
     'her field için frekansa göre verilmiştir.']
    };

explainFig1 = textPage('1. ANALİZ YÖNTEMİ', explain1, 'Arial', 10);
explainFig2 = textPage('2. NYQUIST FREKANSI VE GRAFİKLERİN YORUMU', explain2, 'Arial', 10);

% --- Assemble pages: frame, header and page numbers ----------------------
allFigs = [summaryFig; explainFig1; explainFig2; levelFigs];
nPages = numel(allFigs);
headerText = sprintf('MTF Yield Analizi  |  %s', nominalName);

for p = 1:nPages
    addPageFrame(allFigs(p), headerText, p, nPages);
end

% PNG of the plot pages
for y = 1:nLevels
    exportPng(levelFigs(y), fullfile(outFolder, levelFiles(y) + ".png"));
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
        print(allFigs(p), fullfile(outFolder, sprintf('MTF_Report_page%d.pdf', p)), '-dpdf');
    end
    fprintf('Report saved as separate pages (MTF_Report_page*.pdf).\n');
end

close(allFigs);

fprintf('\nResults saved to:\n%s\n', outFolder);


function v = yieldValue(data, yieldPct)
% data: nSystems x nFreq. Returns, for each frequency, the MTF value that
% at least yieldPct % of the systems meet or exceed (no toolbox needed).
sorted = sort(data, 1, 'descend');
idx = max(1, ceil(yieldPct / 100 * size(data, 1)));
v = sorted(idx, :);
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
