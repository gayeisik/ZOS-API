%% YIELD'E GÖRE EFL ANALİZİ (ZOS-API, MATLAB)
% NE YAPAR
%   EFL dağılımını tüm Monte Carlo dosyalarından değil, MTF'e göre belirlenen
%   yield kümesinden çıkarır. Örneğin %90 yield için: MTF'i en iyi olan %90
%   dosya seçilir (MTF'i en kötü %10 dışarıda bırakılır) ve EFL dağılımı ile
%   istatistikleri yalnızca bu dosyalardan hesaplanır. Karşılaştırma için tüm
%   dosyaların dağılımı da verilir.
%
% NASIL ÇALIŞIR
%   1. Klasör seçilir; adı filePrefix ("MC_T") ile başlayan .zmx dosyaları bulunur.
%   2. Her dosya için iki değer gerekir: EFL ve bir MTF ölçütü.
%      useSavedResults = true ise ve klasörde önceki analizlerin sonuçları
%      varsa onlar kullanılır (OpticStudio açılmaz, saniyeler sürer):
%        - MTF : "MTF ANALİZİ/MTF_RawData.mat"   (MTF_Yield_Analysis.m çıktısı)
%        - EFL : "EFL ANALİZİ/EFL_Results.xlsx"  (EFL_Batch_Analysis.m çıktısı)
%      Eksik olan değer OpticStudio ile hesaplanır: EFL için EFFL operandı, MTF
%      için polikromatik FFT MTF (tüm dalga boyları, tüm fieldler, T ve S).
%      Nominal dosya (adı "MC_" ile başlamayan .zmx) afokal moddaysa MTF
%      frekansı  ν[cyc/mm] = ν[cyc/mrad] × 1000 / |EFL|  ile çevrilir.
%   3. Dosyalar MTF ölçütüne göre sıralanır ve en iyi %yieldPct seçilir.
%
% HESAPLAR
%   Nyquist frekansı = 1000 / (2 × piksel aralığı[µm])  [cyc/mm]
%   MTF ölçütü (dosya başına tek sayı) = mtfFreqFrac × Nyquist frekansında,
%     seçilen fieldlerin (mtfFields; boş = tüm fieldler) tanjansiyel ve sajital
%     MTF değerlerinin EN DÜŞÜĞÜ. Yani dosyanın en kötü field ve yöndeki MTF'i.
%   Yield eşiği = dosyaların en az %yieldPct'inin sağladığı (eşit veya üstünde
%     kaldığı) MTF ölçütü değeri. Ölçütü bu eşiğe eşit veya üstünde olan
%     dosyalar yield kümesidir.
%   EFL istatistikleri (sayı, ortalama, medyan, standart sapma, min, max) hem
%     tüm geçerli dosyalar hem de yield kümesi için hesaplanır.
%
% AYARLAR (User settings)
%   filePrefix, pixelPitch_um, yieldPct, mtfFreqFrac, mtfFields, useSavedResults
%
% ÇIKTILAR (seçilen klasörde "EFL ANALİZİ - %90 YIELD" klasörüne)
%   EFL_Yield_Histogram.png : tüm dosyaların ve yield kümesinin EFL histogramı
%   EFL_Yield_Results.xlsx  : "EFL by File" (dosya adı, EFL, MTF ölçütü, sıra,
%                             yield kümesinde mi), "Statistics", "Settings"
%   Command Window          : eşik değeri ve iki kümenin EFL istatistikleri

clc;
clear;
close all;

%% User settings
filePrefix      = 'MC_T';   % only files starting with this are used
pixelPitch_um   = 15;       % detector pixel pitch (um)
yieldPct        = 90;       % yield level (%): best yieldPct % of the files by MTF are kept
mtfFreqFrac     = 0.5;      % MTF criterion frequency as a fraction of Nyquist (0.5 = Nyquist/2)
mtfFields       = [];       % fields used for the criterion; [] = all fields, e.g. 1 = field 1 only
useSavedResults = true;     % true: use MTF_RawData.mat / EFL_Results.xlsx if they exist in the folder
% FFT MTF sample size is set in computeFileMtf (S_64x64 by default)

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
fileNames = string({files.name})';

EFL_mm = nan(nFiles, 1);
mtfMetric = nan(nFiles, 1);      % worst-field MTF at the criterion frequency
mtfFreq = NaN;                   % criterion frequency (cyc/mm)
mtfSource = "";
eflSource = "";

%% Saved results of earlier analyses (no OpticStudio needed)
analysisSuffix = ['ANAL' char(304) 'Z' char(304)];   % "ANALİZİ" (char(304) = 'İ')
mtfMatFile = fullfile(folderPath, ['MTF ' analysisSuffix], 'MTF_RawData.mat');
eflXlsFile = fullfile(folderPath, ['EFL ' analysisSuffix], 'EFL_Results.xlsx');

if useSavedResults && isfile(mtfMatFile)
    saved = load(mtfMatFile);
    if isfield(saved, 'maxFreq')
        savedNyquist = saved.maxFreq;
    else
        savedNyquist = max(saved.freqGrid);
    end
    mtfFreq = mtfFreqFrac * savedNyquist;

    % MTF at the criterion frequency: 1 x nSavedFiles x nFields -> nSavedFiles x nFields
    savedT = permute(interp1(saved.freqGrid, permute(saved.MTF_T, [2 1 3]), mtfFreq), [2 3 1]);
    savedS = permute(interp1(saved.freqGrid, permute(saved.MTF_S, [2 1 3]), mtfFreq), [2 3 1]);
    savedMetric = fileMetric(savedT, savedS, mtfFields);
    savedMetric(string(saved.status) ~= "OK") = NaN;

    [found, loc] = ismember(lower(fileNames), lower(string(saved.fileNames)));
    mtfMetric(found) = savedMetric(loc(found));
    mtfSource = "MTF_RawData.mat";
    fprintf('MTF: %d / %d files read from %s\n', sum(isfinite(mtfMetric)), nFiles, mtfMatFile);

    if isfield(saved, 'pixelPitch_um') && saved.pixelPitch_um ~= pixelPitch_um
        warning(['Saved MTF data was computed for a pixel pitch of %g um (Nyquist %.2f cyc/mm), ' ...
            'not %g um. The saved Nyquist frequency is used.'], saved.pixelPitch_um, savedNyquist, pixelPitch_um);
    end
end

if useSavedResults && isfile(eflXlsFile)
    try
        EFLTable = readtable(eflXlsFile, 'Sheet', 'EFL by File');
    catch
        EFLTable = readtable(eflXlsFile);
    end
    okRows = string(EFLTable.Status) == "OK" & isfinite(EFLTable.EFL_mm);
    [found, loc] = ismember(lower(fileNames), lower(string(EFLTable.FileName)));
    found(found) = okRows(loc(found));
    EFL_mm(found) = EFLTable.EFL_mm(loc(found));
    eflSource = "EFL_Results.xlsx";
    fprintf('EFL: %d / %d files read from %s\n', sum(isfinite(EFL_mm)), nFiles, eflXlsFile);
end

needMtf = any(~isfinite(mtfMetric));
needEfl = any(~isfinite(EFL_mm));

%% Compute what is missing with OpticStudio
if needMtf || needEfl
    % If the MTF has to be computed, it is computed for ALL files so that every
    % file uses the same frequency and settings.
    if needMtf
        mtfMetric(:) = NaN;
        mtfFreq = mtfFreqFrac * 1000 / (2 * pixelPitch_um);   % fraction of Nyquist (cyc/mm)
        mtfSource = "OpticStudio FFT MTF";
    end
    if needEfl
        eflSource = "OpticStudio EFFL";
    end

    % Nominal file (name does not start with "MC_"): afocal check for the MTF frequency
    freqScale = 1;   % Zemax frequency -> cyc/mm (1000/|EFL| for afocal files)
    nominalName = '';
    if needMtf
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
    end

    TheApplication = InitConnection();

    if isempty(TheApplication)
        error('Could not connect to OpticStudio.');
    end

    fprintf('Computing %s for %d files...\n', strjoin([repmat("EFL", 1, needEfl), repmat("MTF", 1, needMtf)], ' and '), nFiles);
    wb = waitbar(0, 'Computing...', 'Name', 'EFL Yield Analysis');

    try
        TheSystem = TheApplication.PrimarySystem;

        if needMtf
            if ~TheSystem.LoadFile(fullfile(folderPath, nominalName), false)
                error('OpticStudio could not load the nominal file: %s', nominalName);
            end

            if logical(TheSystem.SystemData.Aperture.AFocalImageSpace)
                % EFL is read in the file's own (afocal) mode
                nominalEFL = TheSystem.MFE.GetOperandValue( ...
                    ZOSAPI.Editors.MFE.MeritOperandType.EFFL, 0, 0, 0, 0, 0, 0, 0, 0);

                if ~isfinite(nominalEFL) || nominalEFL == 0
                    error('Could not compute the EFL of the nominal file: %s', nominalName);
                end

                freqScale = 1000 / abs(nominalEFL);   % afocal MTF is in cyc/mrad
                fprintf('Nominal file %s is in afocal mode, |EFL| = %.4f mm\n', nominalName, abs(nominalEFL));
            else
                fprintf('Nominal file %s is not in afocal mode.\n', nominalName);
            end
        end

        tStart = tic;

        for k = 1:nFiles
            try
                if ~TheSystem.LoadFile(fullfile(folderPath, files(k).name), false)
                    error('OpticStudio could not load this file.');
                end

                if needEfl && ~isfinite(EFL_mm(k))
                    EFL_mm(k) = TheSystem.MFE.GetOperandValue( ...
                        ZOSAPI.Editors.MFE.MeritOperandType.EFFL, 0, 0, 0, 0, 0, 0, 0, 0);
                end

                if needMtf
                    [T, S] = computeFileMtf(TheSystem, mtfFreq, freqScale);
                    mtfMetric(k) = fileMetric(T, S, mtfFields);
                end
            catch ME
                fprintf('ERROR: %s - %s\n', files(k).name, ME.message);
            end

            if isvalid(wb)
                waitbar(k / nFiles, wb, sprintf('Computing... %d / %d', k, nFiles));
            end
        end

        fprintf('Computation time: %.1f s\n', toc(tStart));

    catch ME
        if isvalid(wb), close(wb); end
        CleanupConnection(TheApplication);
        rethrow(ME);
    end

    if isvalid(wb), close(wb); end
    CleanupConnection(TheApplication);
end

%% Valid files and yield set
valid = isfinite(EFL_mm) & isfinite(mtfMetric);
nValid = sum(valid);

if nValid == 0
    error('No file has both a valid EFL and a valid MTF.');
end

% Yield threshold: the MTF value that at least yieldPct % of the valid files meet or exceed
sortedMetric = sort(mtfMetric(valid), 'descend');
threshold = sortedMetric(max(1, ceil(yieldPct / 100 * nValid)));

inYield = valid & mtfMetric >= threshold;
nYield = sum(inYield);

% Rank by MTF (1 = best) among the valid files
mtfRank = nan(nFiles, 1);
[~, order] = sort(mtfMetric(valid), 'descend');
validIdx = find(valid);
mtfRank(validIdx(order)) = (1:nValid)';

%% Statistical summary
statNames = ["Count"; "Mean"; "Median"; "Std"; "Min"; "Max"];
calcStats = @(x) [numel(x); mean(x); median(x); std(x); min(x); max(x)];
allStats = calcStats(EFL_mm(valid));
yieldStats = calcStats(EFL_mm(inYield));

if isempty(mtfFields)
    fieldText = 'all fields';
else
    fieldText = ['field(s) ' char(strjoin(string(mtfFields), ', '))];
end

fprintf('\nMTF criterion: worst MTF (T/S) of %s at %.2f cyc/mm (%g x Nyquist)\n', ...
    fieldText, mtfFreq, mtfFreqFrac);
fprintf('Yield %g %%: MTF criterion >= %.4f  ->  %d of %d valid files\n', ...
    yieldPct, threshold, nYield, nValid);

fprintf('\n%-8s %16s %16s\n', 'EFL', 'All files', sprintf('%g%% yield', yieldPct));
fprintf('%-8s %16d %16d\n', statNames(1), allStats(1), yieldStats(1));
for i = 2:numel(statNames)
    fprintf('%-8s %13.4f mm %13.4f mm\n', statNames(i), allStats(i), yieldStats(i));
end

%% Output folder: "EFL ANALİZİ - %90 YIELD"
outFolder = fullfile(folderPath, ['EFL ' analysisSuffix sprintf(' - %%%g YIELD', yieldPct)]);

if ~isfolder(outFolder)
    mkdir(outFolder);
end

%% Histogram: EFL distribution of all files and of the yield set (same bins)
fig = figure('Name', 'EFL Distribution - Yield', 'Color', 'w', 'Position', [100 100 900 520]);

[~, binEdges] = histcounts(EFL_mm(valid));
histogram(EFL_mm(valid), binEdges, 'FaceColor', [0.75 0.75 0.75], ...
    'DisplayName', sprintf('All files (%d)', nValid));
hold on;
histogram(EFL_mm(inYield), binEdges, 'FaceColor', [0 0.45 0.74], 'FaceAlpha', 0.85, ...
    'DisplayName', sprintf('%g%% yield (%d)', yieldPct, nYield));
hold off;

xlabel('Effective Focal Length (mm)');
ylabel('Number of optical systems');
title(sprintf('EFL distribution - %g%% yield by MTF at %g x Nyquist (%s)', ...
    yieldPct, mtfFreqFrac, fieldText));
legend('Location', 'best');
grid on;

imageFile = fullfile(outFolder, 'EFL_Yield_Histogram.png');
try
    exportgraphics(fig, imageFile, 'Resolution', 300);   % R2020a+
catch
    saveas(fig, imageFile);
end

%% Export Excel
outputFile = fullfile(outFolder, 'EFL_Yield_Results.xlsx');

% Remove old file so stale rows/sheets are not kept
if isfile(outputFile)
    delete(outputFile);
end

Results = table(fileNames, EFL_mm, mtfMetric, mtfRank, inYield, ...
    'VariableNames', {'FileName', 'EFL_mm', 'MTF_criterion', 'MTF_rank', 'InYieldSet'});
writetable(Results, outputFile, 'Sheet', 'EFL by File');

Stats = table(statNames, allStats, yieldStats, ...
    'VariableNames', {'Statistic', 'EFL_mm_AllFiles', 'EFL_mm_YieldSet'});
writetable(Stats, outputFile, 'Sheet', 'Statistics');

Settings = table( ...
    ["Yield level (%)"; "MTF criterion frequency (cyc/mm)"; "Fraction of Nyquist"; ...
     "MTF criterion threshold"; "Valid files"; "Files in yield set"], ...
    [yieldPct; mtfFreq; mtfFreqFrac; threshold; nValid; nYield], ...
    'VariableNames', {'Parameter', 'Value'});
writetable(Settings, outputFile, 'Sheet', 'Settings');

Sources = table(["MTF criterion fields"; "MTF source"; "EFL source"], ...
    [string(fieldText); mtfSource; eflSource], 'VariableNames', {'Parameter', 'Value'});
writetable(Sources, outputFile, 'Sheet', 'Sources');

fprintf('\nHistogram saved to:\n%s\n', imageFile);
fprintf('Results saved to:\n%s\n', outputFile);


function m = fileMetric(T, S, fieldList)
% Worst MTF of each file: minimum over the chosen fields of tangential and
% sagittal MTF. T, S: nFiles x nFields. fieldList = [] uses all fields.
if ~isempty(fieldList)
    if any(fieldList < 1 | fieldList > size(T, 2))
        error('mtfFields contains a field number outside 1..%d.', size(T, 2));
    end
    T = T(:, fieldList);
    S = S(:, fieldList);
end
m = min([T, S], [], 2);
end

function [T, S] = computeFileMtf(TheSystem, mtfFreq, freqScale)
% Polychromatic FFT MTF of the current system at one frequency (cyc/mm) for
% all fields. T, S: 1 x nFields (tangential, sagittal).
mtf = TheSystem.Analyses.New_FftMtf();
closeMtf = onCleanup(@() mtf.Close());

settings = mtf.GetSettings();
settings.MaximumFrequency = 1.5 * mtfFreq / freqScale;   % cyc/mrad if afocal
settings.SampleSize = ZOSAPI.Analysis.SampleSizes.S_64x64;
settings.Wavelength.UseAllWavelengths();
settings.Field.UseAllFields();
try
    settings.ShowDiffractionLimit = false;
catch
end

mtf.ApplyAndWaitForCompletion();
results = mtf.GetResults();
nSeries = results.NumberOfDataSeries;
dataSeries = results.DataSeries;

% Skip a diffraction limit series, if OpticStudio returns one
isDiff = false(nSeries, 1);
for i = 1:nSeries
    isDiff(i) = contains(lower(string(char(dataSeries(i).Description))), "diff");
end
fieldIdx = find(~isDiff);

T = nan(1, numel(fieldIdx));
S = T;
for f = 1:numel(fieldIdx)
    data = dataSeries(fieldIdx(f));
    x = data.XData.Data.double;
    y = data.YData.Data.double;

    x = x(:) * freqScale;   % cyc/mm
    if size(y, 1) ~= numel(x)
        y = y.';
    end

    T(f) = interp1(x, y(:, 1), mtfFreq);
    if size(y, 2) > 1
        S(f) = interp1(x, y(:, 2), mtfFreq);
    else
        S(f) = T(f);
    end
end
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
