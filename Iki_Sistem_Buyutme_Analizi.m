%% İKİ SİSTEM BÜYÜTME TOLERANS ANALİZİ (ZOS-API, MATLAB)
% NE YAPAR
%   Birlikte çalışan iki sistemin (1: MWIR objektif + dedektör, 2: OLED +
%   gündüz kanalına aktaran optik) Monte Carlo tolerans dosyalarından toplam
%   büyütmenin istatistiğini çıkarır. Her sistemde MTF'e göre %yieldPct yield
%   kümesine giren dosyalar seçilir; 1. sistemin yield kümesindeki HER EFL,
%   2. sistemin yield kümesindeki HER EFL ile eşleştirilir (tam çaprazlama,
%   n1 x n2 çift). Her çiftin büyütmesi hesaplanır ve büyütme isterine
%   (magTarget ± magTol) göre istatistik verilir.
%   Başka bir scriptin çıktısı kullanılmaz; EFL ve MTF bu scriptte
%   OpticStudio ile hesaplanır.
%
% NASIL ÇALIŞIR
%   1. Kullanıcı önce 1. sistemin (MWIR), sonra 2. sistemin (OLED) Monte Carlo
%      klasörünü seçer. Her klasörde adı filePrefix ("MC_T") ile başlayan .zmx
%      dosyaları MC dosyalarıdır; adı "MC_" ile başlamayan .zmx dosyası nominal
%      tasarımdır (birden fazlaysa listeden seçilir).
%   2. Her sistem için OpticStudio ile:
%        - Nominal dosyanın |EFL|'si okunur ve afokal mod kontrol edilir.
%        - Her MC dosyasının EFL'si (EFFL operandı) okunur. EFL HER ZAMAN
%          mutlak değer olarak kullanılır.
%        - Her MC dosyasının polikromatik FFT MTF'i (tüm dalga boyları, tüm
%          fieldler, T ve S) o sistemin Nyquist'inin mtfFreqFrac katında
%          hesaplanır. Afokal dosyalarda Zemax cyc/mrad verir; frekans
%          ν[cyc/mm] = ν[cyc/mrad] × 1000 / |EFL_nominal|  ile çevrilir.
%   3. Her sistemde dosyalar MTF ölçütüne göre sıralanır, en iyi %yieldPct
%      seçilir (yield kümesi).
%   4. İki yield kümesi çaprazlanır ve büyütme istatistiği çıkarılır.
%
% HESAPLAR
%   Nyquist frekansı = 1000 / (2 × piksel aralığı[µm])  [cyc/mm], her sistem
%     kendi piksel aralığıyla (MWIR 15 µm, OLED 5.65 µm).
%   MTF ölçütü (dosya başına tek sayı) = mtfFreqFrac × Nyquist frekansında,
%     seçilen fieldlerin tanjansiyel ve sajital MTF'lerinin EN DÜŞÜĞÜ.
%   Yield eşiği = dosyaların en az %yieldPct'inin sağladığı MTF ölçütü değeri.
%   Büyütme:
%       M = (|EFL_1| / |EFL_2|) × (k × p_2 / p_1)
%     p_1 = dedektör piksel aralığı, p_2 = OLED piksel aralığı,
%     k   = bir dedektör pikselinin OLED'de kapladığı piksel sayısı
%           (512 piksellik görüntü 1024 piksele basılıyor -> k = 2).
%     Nominal: M = 1 için |EFL_2| = 0.7533 × |EFL_1|.
%   Geçme koşulu: |M - magTarget| <= magTol.
%   İstatistikler:
%     - Geçen çift oranı (rastgele eşleşen iki sistemin geçme olasılığı)
%     - M'nin ortalama, std, medyan, min, max ve yüzdelikleri
%     - Alt / üst sınır dışında kalan oranlar
%     - Cp, Cpk ve normal dağılım varsayımıyla tahmini geçme oranı
%     - Ortalama M hedefe kaydırılsaydı geçme oranı (nominal EFL düzeltme
%       potansiyeli)
%     - Büyütme hatasına katkı: var(ln M) = var(ln EFL_1) + var(ln EFL_2);
%       her sistemin payı yüzde olarak
%     - Dosya bazında geçme oranı: bir sistemin her dosyası diğer sistemin
%       dosyalarının yüzde kaçıyla isteri sağlıyor
%     - Nominal EFL'ler için diğer sistemin EFL'sine izin verilen aralık
%
% AYARLAR (User settings)
%   sysName, pixelPitch_um, nPixels, pixelScale, filePrefix, yieldPct,
%   mtfFreqFrac, mtfFields, magTarget, magTol, writePairsSheet, nHistBins
%
% ÇIKTILAR (1. sistemin klasöründe "BÜYÜTME ANALİZİ - %90 YIELD" klasörüne)
%   Magnification_Histogram.png : tüm çiftlerin büyütme dağılımı ve tolerans bandı
%   EFL_Pair_Map.png            : EFL_1 - EFL_2 düzleminde geçen / kalan çiftler
%   PassRate_by_File.png        : her dosyanın diğer sistemle geçme oranı
%   EFL_Distributions.png       : her sistemin EFL dağılımı (tümü ve yield kümesi)
%   Magnification_Results.xlsx  : "Summary", "Percentiles", sistem başına dosya
%                                 tabloları, "Pairs" (tüm çiftler), "Settings"
%   Magnification_Data.mat      : tüm ham veriler

clc;
clear;
close all;

%% User settings
sysName         = ["MWIR", "OLED"];  % system 1, system 2
pixelPitch_um   = [15, 5.65];        % pixel pitch of detector / OLED (um)
nPixels         = [512, 1024];       % pixel count of detector / OLED
pixelScale      = nPixels(2) / nPixels(1);   % k: OLED pixels per detector pixel
filePrefix      = 'MC_T';            % only files starting with this are used
yieldPct        = 90;                % yield level (%): best yieldPct % of each system's files by MTF
mtfFreqFrac     = 0.5;               % MTF criterion frequency as a fraction of each system's Nyquist
mtfFields       = {[], []};          % fields used for the MTF criterion per system; [] = all fields
magTarget       = 1.00;              % required magnification
magTol          = 0.03;              % allowed deviation: |M - magTarget| <= magTol
writePairsSheet = true;              % write every pair to the "Pairs" sheet (Excel limit 1048575 rows)
nHistBins       = 60;                % histogram bins of the magnification
% FFT MTF sample size is set in computeFileMtf (S_64x64 by default)

magConst = pixelScale * pixelPitch_um(2) / pixelPitch_um(1);   % M = magConst * |EFL_1| / |EFL_2|
magLow   = magTarget - magTol;
magHigh  = magTarget + magTol;

%% Select the two Monte Carlo folders
sys = struct([]);

for s = 1:2
    folderPath = uigetdir(pwd, sprintf('Select the %s Monte Carlo folder (system %d)', sysName(s), s));
    if isequal(folderPath, 0)
        return;
    end

    if s == 2 && strcmpi(folderPath, sys(1).folder)
        error('The same folder was selected for both systems.');
    end

    [mcFiles, nominalName] = findFiles(folderPath, filePrefix, sysName(s));
    if isempty(nominalName)
        return;   % nominal file selection cancelled
    end

    sys(s).name        = sysName(s);
    sys(s).folder      = folderPath;
    sys(s).pitch       = pixelPitch_um(s);
    sys(s).fileNames   = mcFiles;
    sys(s).nominalName = nominalName;

    fprintf('%s: %d MC files, nominal file %s\n', sysName(s), numel(mcFiles), nominalName);
end

%% Compute EFL and MTF of both systems with OpticStudio
TheApplication = InitConnection();

if isempty(TheApplication)
    error('Could not connect to OpticStudio.');
end

tStart = tic;

try
    TheSystem = TheApplication.PrimarySystem;

    for s = 1:2
        [sys(s).EFL, sys(s).mtfMetric, sys(s).nominalEFL, sys(s).isAfocal, sys(s).mtfFreq] = ...
            computeSystem(TheSystem, sys(s), mtfFreqFrac, mtfFields{s});
    end
catch ME
    CleanupConnection(TheApplication);
    rethrow(ME);
end

CleanupConnection(TheApplication);
fprintf('Computation time: %.1f s\n', toc(tStart));

%% Yield set of each system
for s = 1:2
    [sys(s).valid, sys(s).inYield, sys(s).threshold, sys(s).mtfRank] = ...
        selectYield(sys(s).EFL, sys(s).mtfMetric, yieldPct);

    if ~any(sys(s).inYield)
        error('%s: no file has both a valid EFL and a valid MTF.', sys(s).name);
    end

    fprintf('%s yield %g %%: MTF >= %.4f at %.2f cyc/mm -> %d of %d valid files\n', ...
        sys(s).name, yieldPct, sys(s).threshold, sys(s).mtfFreq, sum(sys(s).inYield), sum(sys(s).valid));
end

%% Cross-match: every yield EFL of system 1 with every yield EFL of system 2
E1 = sys(1).EFL(sys(1).inYield);           % n1 x 1 (absolute values)
E2 = sys(2).EFL(sys(2).inYield).';         % 1 x n2
names1 = sys(1).fileNames(sys(1).inYield);
names2 = sys(2).fileNames(sys(2).inYield);
n1 = numel(E1);
n2 = numel(E2);

M = magConst * E1 ./ E2;                    % n1 x n2 magnification matrix
pass = abs(M - magTarget) <= magTol + 1e-12;

mAll   = M(:);
nPairs = numel(mAll);
nPass  = sum(pass(:));
passPct = 100 * nPass / nPairs;

%% Statistics of the magnification
mMean = mean(mAll);
mStd  = std(mAll);
mMed  = median(mAll);
mMin  = min(mAll);
mMax  = max(mAll);
belowPct = 100 * mean(mAll < magLow - 1e-12);
abovePct = 100 * mean(mAll > magHigh + 1e-12);

pctLevels = [0.135 1 2.5 5 25 50 75 95 97.5 99 99.865]';
pctValues = quantileSorted(sort(mAll), pctLevels / 100);

Cp  = (magHigh - magLow) / (6 * mStd);
Cpk = min(magHigh - mMean, mMean - magLow) / (3 * mStd);
normPhi = @(x) 0.5 * erfc(-x / sqrt(2));
normPassPct = 100 * (normPhi((magHigh - mMean) / mStd) - normPhi((magLow - mMean) / mStd));

% If the mean magnification were moved to the target (nominal EFL correction)
mCentered = mAll * (magTarget / mMean);
centeredPassPct = 100 * mean(abs(mCentered - magTarget) <= magTol + 1e-12);

% Nominal magnification and EFL windows
mNominal = magConst * sys(1).nominalEFL / sys(2).nominalEFL;
efl2Window = magConst * sys(1).nominalEFL ./ [magHigh, magLow];   % allowed |EFL_2| for nominal |EFL_1|
efl1Window = sys(2).nominalEFL / magConst * [magLow, magHigh];    % allowed |EFL_1| for nominal |EFL_2|

% Contribution of each system to the magnification error: var(ln M) = var(ln E1) + var(ln E2)
varLn = [var(log(E1)), var(log(E2))];
contribPct = 100 * varLn / sum(varLn);

% Pass rate of each file against all files of the other system
passRate1 = 100 * mean(pass, 2);           % n1 x 1
passRate2 = 100 * mean(pass, 1).';         % n2 x 1

sys(1).passRate = nan(numel(sys(1).fileNames), 1);
sys(1).passRate(sys(1).inYield) = passRate1;
sys(2).passRate = nan(numel(sys(2).fileNames), 1);
sys(2).passRate(sys(2).inYield) = passRate2;

anyMatch = [100 * mean(passRate1 > 0), 100 * mean(passRate2 > 0)];
allMatch = [100 * mean(passRate1 == 100), 100 * mean(passRate2 == 100)];

%% Command window summary
fprintf('\n===== MAGNIFICATION ANALYSIS =====\n');
fprintf('M = (|EFL_%s| / |EFL_%s|) x (%g x %g um / %g um) = %.6f x |EFL_1| / |EFL_2|\n', ...
    sys(1).name, sys(2).name, pixelScale, pixelPitch_um(2), pixelPitch_um(1), magConst);
fprintf('Requirement: M = %.4f +/- %.4f  [%.4f, %.4f]\n', magTarget, magTol, magLow, magHigh);
fprintf('Nominal: |EFL_1| = %.4f mm, |EFL_2| = %.4f mm, M = %.5f\n', ...
    sys(1).nominalEFL, sys(2).nominalEFL, mNominal);
fprintf('Pairs: %d x %d = %d\n', n1, n2, nPairs);
fprintf('\nPASS: %d of %d pairs = %.2f %%\n', nPass, nPairs, passPct);
fprintf('Below %.4f: %.2f %%   Above %.4f: %.2f %%\n', magLow, belowPct, magHigh, abovePct);
fprintf('\nM mean %.5f  std %.5f  median %.5f  min %.5f  max %.5f\n', mMean, mStd, mMed, mMin, mMax);
fprintf('Cp %.3f  Cpk %.3f  normal-fit pass %.2f %%\n', Cp, Cpk, normPassPct);
fprintf('Pass rate if mean M were centered at %.4f: %.2f %%\n', magTarget, centeredPassPct);
fprintf('\nContribution to magnification variance: %s %.1f %%, %s %.1f %%\n', ...
    sys(1).name, contribPct(1), sys(2).name, contribPct(2));
fprintf('Allowed |EFL_%s| for nominal |EFL_%s|: %.4f - %.4f mm\n', ...
    sys(2).name, sys(1).name, efl2Window(1), efl2Window(2));
fprintf('Allowed |EFL_%s| for nominal |EFL_%s|: %.4f - %.4f mm\n', ...
    sys(1).name, sys(2).name, efl1Window(1), efl1Window(2));
for s = 1:2
    fprintf('%s files matching at least one / all files of the other system: %.1f %% / %.1f %%\n', ...
        sys(s).name, anyMatch(s), allMatch(s));
end

%% Output folder: "BÜYÜTME ANALİZİ - %90 YIELD" in the system 1 folder
outName = ['B' char(220) 'Y' char(220) 'TME ANAL' char(304) 'Z' char(304) sprintf(' - %%%g YIELD', yieldPct)];
outFolder = fullfile(sys(1).folder, outName);

if ~isfolder(outFolder)
    mkdir(outFolder);
end

colPass = [0.47 0.67 0.19];
colFail = [0.85 0.33 0.10];
colBlue = [0 0.45 0.74];
eflLabel = @(s) sprintf('|EFL| %s (mm)', sys(s).name);

%% Figure 1: magnification histogram with tolerance band
fig1 = figure('Name', 'Magnification Distribution', 'Color', 'w', 'Position', [100 100 950 540]);
h = histogram(mAll, nHistBins, 'FaceColor', colBlue, 'DisplayName', sprintf('All pairs (%d)', nPairs));
hold on;
yl = ylim;
band = patch([magLow magHigh magHigh magLow], [0 0 yl(2) yl(2)], colPass, ...
    'FaceAlpha', 0.18, 'EdgeColor', 'none', 'DisplayName', sprintf('Requirement %.2f \\pm %.2f', magTarget, magTol));
uistack(band, 'bottom');
xline(magTarget, '-', 'Color', colPass, 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(magLow, '--', 'Color', colPass, 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline(magHigh, '--', 'Color', colPass, 'LineWidth', 1.2, 'HandleVisibility', 'off');
hNom = xline(mNominal, ':', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', sprintf('Nominal M = %.4f', mNominal));
hMean = xline(mMean, '-.', 'Color', colFail, 'LineWidth', 1.5, 'DisplayName', sprintf('Mean M = %.4f', mMean));
hold off;
xlim([min(mMin, magLow) max(mMax, magHigh)] + [-1 1] * 0.05 * (max(mMax, magHigh) - min(mMin, magLow)));
xlabel('Magnification M');
ylabel('Number of pairs');
title(sprintf('Magnification of %s x %s pairs (%g%% MTF yield sets): PASS %.2f%%', ...
    sys(1).name, sys(2).name, yieldPct, passPct));
legend([h band hNom hMean], 'Location', 'best');
grid on;
saveFigure(fig1, fullfile(outFolder, 'Magnification_Histogram.png'));

%% Figure 2: EFL pair map (pass / fail) with the requirement band
fig2 = figure('Name', 'EFL Pair Map', 'Color', 'w', 'Position', [120 120 900 700]);
[G1, G2] = ndgrid(E1, E2);
plot(G1(~pass), G2(~pass), '.', 'Color', colFail, 'MarkerSize', 5, 'DisplayName', 'Fail');
hold on;
plot(G1(pass), G2(pass), '.', 'Color', colPass, 'MarkerSize', 5, 'DisplayName', 'Pass');
e1Line = linspace(min(E1), max(E1), 50);
plot(e1Line, magConst * e1Line / magTarget, '-k', 'LineWidth', 1.2, 'DisplayName', sprintf('M = %.2f', magTarget));
plot(e1Line, magConst * e1Line / magLow, '--k', 'LineWidth', 1, 'DisplayName', sprintf('M = %.2f', magLow));
plot(e1Line, magConst * e1Line / magHigh, '--k', 'LineWidth', 1, 'DisplayName', sprintf('M = %.2f', magHigh));
plot(sys(1).nominalEFL, sys(2).nominalEFL, 'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'y', 'DisplayName', 'Nominal');
hold off;
xlabel(eflLabel(1));
ylabel(eflLabel(2));
title(sprintf('EFL pairs: %d of %d pass (%.2f%%)', nPass, nPairs, passPct));
legend('Location', 'best');
grid on;
saveFigure(fig2, fullfile(outFolder, 'EFL_Pair_Map.png'));

%% Figure 3: pass rate of each file against the other system
fig3 = figure('Name', 'Pass Rate by File', 'Color', 'w', 'Position', [140 140 950 700]);
passRates = {passRate1, passRate2};
eflSets = {E1, E2(:)};
for s = 1:2
    subplot(2, 1, s);
    plot(eflSets{s}, passRates{s}, 'o', 'Color', colBlue, 'MarkerFaceColor', colBlue, 'MarkerSize', 4);
    hold on;
    xline(sys(s).nominalEFL, ':k', 'Nominal', 'LineWidth', 1.2);
    hold off;
    ylim([-5 105]);
    xlabel(eflLabel(s));
    ylabel(sprintf('Pass rate with %s files (%%)', sys(3 - s).name));
    title(sprintf('%s files: %.1f%% match at least one, %.1f%% match all %s files', ...
        sys(s).name, anyMatch(s), allMatch(s), sys(3 - s).name));
    grid on;
end
saveFigure(fig3, fullfile(outFolder, 'PassRate_by_File.png'));

%% Figure 4: EFL distribution of each system (all valid files and yield set)
fig4 = figure('Name', 'EFL Distributions', 'Color', 'w', 'Position', [160 160 950 700]);
for s = 1:2
    subplot(2, 1, s);
    eflValid = sys(s).EFL(sys(s).valid);
    [~, binEdges] = histcounts(eflValid);
    histogram(eflValid, binEdges, 'FaceColor', [0.75 0.75 0.75], ...
        'DisplayName', sprintf('All valid files (%d)', numel(eflValid)));
    hold on;
    histogram(sys(s).EFL(sys(s).inYield), binEdges, 'FaceColor', colBlue, 'FaceAlpha', 0.85, ...
        'DisplayName', sprintf('%g%% yield (%d)', yieldPct, sum(sys(s).inYield)));
    xline(sys(s).nominalEFL, ':k', 'LineWidth', 1.5, 'DisplayName', sprintf('Nominal %.4f mm', sys(s).nominalEFL));
    hold off;
    xlabel(eflLabel(s));
    ylabel('Number of systems');
    title(sprintf('%s: |EFL| distribution, %g%% yield by MTF at %.2f cyc/mm', ...
        sys(s).name, yieldPct, sys(s).mtfFreq));
    legend('Location', 'best');
    grid on;
end
saveFigure(fig4, fullfile(outFolder, 'EFL_Distributions.png'));

%% Export Excel
outputFile = fullfile(outFolder, 'Magnification_Results.xlsx');

if isfile(outputFile)
    delete(outputFile);
end

summaryNames = [ ...
    "Pairs (n1 x n2)"; "Passing pairs"; "Pass rate"; "Below lower limit"; "Above upper limit"; ...
    "Mean M"; "Std M"; "Median M"; "Min M"; "Max M"; "Nominal M"; ...
    "Cp"; "Cpk"; "Normal-fit pass rate"; "Pass rate if mean M centered at target"; ...
    "Variance contribution " + sys(1).name; "Variance contribution " + sys(2).name; ...
    sys(1).name + " files matching >= 1 " + sys(2).name + " file"; ...
    sys(1).name + " files matching all " + sys(2).name + " files"; ...
    sys(2).name + " files matching >= 1 " + sys(1).name + " file"; ...
    sys(2).name + " files matching all " + sys(1).name + " files"; ...
    "Nominal |EFL| " + sys(1).name; "Nominal |EFL| " + sys(2).name; ...
    "Allowed |EFL| " + sys(2).name + " min (nominal " + sys(1).name + ")"; ...
    "Allowed |EFL| " + sys(2).name + " max (nominal " + sys(1).name + ")"; ...
    "Allowed |EFL| " + sys(1).name + " min (nominal " + sys(2).name + ")"; ...
    "Allowed |EFL| " + sys(1).name + " max (nominal " + sys(2).name + ")"];
summaryValues = [nPairs; nPass; passPct; belowPct; abovePct; ...
    mMean; mStd; mMed; mMin; mMax; mNominal; ...
    Cp; Cpk; normPassPct; centeredPassPct; contribPct(1); contribPct(2); ...
    anyMatch(1); allMatch(1); anyMatch(2); allMatch(2); ...
    sys(1).nominalEFL; sys(2).nominalEFL; efl2Window(1); efl2Window(2); efl1Window(1); efl1Window(2)];
summaryUnits = ["count"; "count"; "%"; "%"; "%"; ""; ""; ""; ""; ""; ""; ""; ""; "%"; "%"; "%"; "%"; ...
    "%"; "%"; "%"; "%"; "mm"; "mm"; "mm"; "mm"; "mm"; "mm"];
writetable(table(summaryNames, summaryValues, summaryUnits, ...
    'VariableNames', {'Parameter', 'Value', 'Unit'}), outputFile, 'Sheet', 'Summary');

writetable(table(pctLevels, pctValues, 'VariableNames', {'Percentile', 'Magnification'}), ...
    outputFile, 'Sheet', 'Percentiles');

for s = 1:2
    fileTable = table(sys(s).fileNames, sys(s).EFL, sys(s).mtfMetric, sys(s).mtfRank, ...
        sys(s).inYield, sys(s).passRate, ...
        'VariableNames', {'FileName', 'EFL_abs_mm', 'MTF_criterion', 'MTF_rank', 'InYieldSet', ...
        char("PassRate_with_" + sys(3 - s).name + "_pct")});
    writetable(fileTable, outputFile, 'Sheet', char(sys(s).name + " Files"));
end

if writePairsSheet
    if nPairs <= 1048575
        [I, J] = ndgrid(1:n1, 1:n2);
        pairTable = table(names1(I(:)), names2(J(:)), E1(I(:)), E2(J(:)).', mAll, pass(:), ...
            'VariableNames', {char("File_" + sys(1).name), char("File_" + sys(2).name), ...
            char("EFL_" + sys(1).name + "_mm"), char("EFL_" + sys(2).name + "_mm"), ...
            'Magnification', 'Pass'});
        writetable(pairTable, outputFile, 'Sheet', 'Pairs');
    else
        warning('%d pairs exceed the Excel row limit; the "Pairs" sheet is not written (see the .mat file).', nPairs);
    end
end

settingNames = [ ...
    "System 1"; "System 1 folder"; "System 1 nominal file"; "System 1 afocal"; ...
    "System 1 pixel pitch (um)"; "System 1 pixels"; "System 1 MTF frequency (cyc/mm)"; ...
    "System 1 MTF threshold"; "System 1 MTF fields"; "System 1 valid / yield files"; ...
    "System 2"; "System 2 folder"; "System 2 nominal file"; "System 2 afocal"; ...
    "System 2 pixel pitch (um)"; "System 2 pixels"; "System 2 MTF frequency (cyc/mm)"; ...
    "System 2 MTF threshold"; "System 2 MTF fields"; "System 2 valid / yield files"; ...
    "File prefix"; "Yield level (%)"; "MTF frequency (fraction of Nyquist)"; ...
    "Pixel scale k (OLED px per detector px)"; "Magnification constant k*p2/p1"; ...
    "Magnification formula"; "Target magnification"; "Tolerance (+/-)"];
settingValues = strings(numel(settingNames), 1);
for s = 1:2
    o = (s - 1) * 10;
    settingValues(o + (1:10)) = [sys(s).name; string(sys(s).folder); string(sys(s).nominalName); ...
        string(sys(s).isAfocal); string(sys(s).pitch); string(nPixels(s)); string(sys(s).mtfFreq); ...
        string(sys(s).threshold); fieldText(mtfFields{s}); ...
        sprintf('%d / %d', sum(sys(s).valid), sum(sys(s).inYield))];
end
settingValues(21:end) = [string(filePrefix); string(yieldPct); string(mtfFreqFrac); ...
    string(pixelScale); string(magConst); "M = k*p2/p1 * |EFL_1| / |EFL_2|"; ...
    string(magTarget); string(magTol)];
writetable(table(settingNames, settingValues, 'VariableNames', {'Parameter', 'Value'}), ...
    outputFile, 'Sheet', 'Settings');

%% Save raw data
save(fullfile(outFolder, 'Magnification_Data.mat'), 'sys', 'E1', 'E2', 'names1', 'names2', 'M', 'pass', ...
    'magConst', 'magTarget', 'magTol', 'yieldPct', 'mtfFreqFrac', 'pixelPitch_um', 'nPixels', 'pixelScale');

fprintf('\nResults saved to:\n%s\n', outFolder);


%% Local functions
function [mcNames, nominalName] = findFiles(folderPath, filePrefix, sysLabel)
% MC files (name starts with filePrefix) and the nominal file (name does not
% start with "MC_") of one folder. nominalName = '' if the selection is cancelled.
allZmx = [dir(fullfile(folderPath, '*.zmx')); dir(fullfile(folderPath, '*.ZMX'))];
if isempty(allZmx)
    error('%s: no .zmx file found in %s', sysLabel, folderPath);
end

% Windows dir() is case-insensitive; remove duplicate names
[~, uniqueIdx] = unique(lower(string({allZmx.name})), 'stable');
allZmx = allZmx(uniqueIdx);
allNames = string({allZmx.name})';

mcNames = allNames(startsWith(allNames, filePrefix));
if isempty(mcNames)
    error('%s: no .zmx file starting with "%s" found in %s', sysLabel, filePrefix, folderPath);
end

nominalFiles = allNames(~startsWith(allNames, 'MC_', 'IgnoreCase', true));
if isempty(nominalFiles)
    error('%s: no nominal .zmx file (name not starting with "MC_") found in %s', sysLabel, folderPath);
elseif isscalar(nominalFiles)
    nominalName = char(nominalFiles);
else
    [sel, ok] = listdlg('ListString', cellstr(nominalFiles), 'SelectionMode', 'single', ...
        'PromptString', sprintf('Select the %s nominal file:', sysLabel), 'ListSize', [400 200]);
    if ~ok
        nominalName = '';
        return;
    end
    nominalName = char(nominalFiles(sel));
end
end

function [EFL, metric, nominalEFL, isAfocal, mtfFreq] = computeSystem(TheSystem, S, mtfFreqFrac, mtfFields)
% |EFL| and MTF criterion of every MC file of one system.
nFiles = numel(S.fileNames);
EFL = nan(nFiles, 1);
metric = nan(nFiles, 1);
mtfFreq = mtfFreqFrac * 1000 / (2 * S.pitch);   % cyc/mm

% Nominal file: |EFL| and afocal mode
if ~TheSystem.LoadFile(char(fullfile(S.folder, S.nominalName)), false)
    error('%s: OpticStudio could not load the nominal file %s', S.name, S.nominalName);
end
isAfocal = logical(TheSystem.SystemData.Aperture.AFocalImageSpace);
nominalEFL = abs(readEfl(TheSystem));
if ~isfinite(nominalEFL) || nominalEFL == 0
    error('%s: could not compute the EFL of the nominal file %s', S.name, S.nominalName);
end

freqScale = 1;   % Zemax frequency -> cyc/mm
if isAfocal
    freqScale = 1000 / nominalEFL;   % afocal MTF is in cyc/mrad
end
fprintf('%s nominal %s: |EFL| = %.4f mm, afocal = %d, MTF criterion at %.2f cyc/mm\n', ...
    S.name, S.nominalName, nominalEFL, isAfocal, mtfFreq);

wb = waitbar(0, sprintf('%s: computing...', S.name), 'Name', 'Magnification Analysis');
closeWb = onCleanup(@() delete(wb(isvalid(wb))));

for k = 1:nFiles
    try
        if ~TheSystem.LoadFile(char(fullfile(S.folder, S.fileNames(k))), false)
            error('OpticStudio could not load this file.');
        end
        EFL(k) = abs(readEfl(TheSystem));
        [T, Sg] = computeFileMtf(TheSystem, mtfFreq, freqScale);
        metric(k) = fileMetric(T, Sg, mtfFields);
    catch ME
        fprintf('ERROR: %s %s - %s\n', S.name, S.fileNames(k), ME.message);
    end

    if isvalid(wb)
        waitbar(k / nFiles, wb, sprintf('%s: computing... %d / %d', S.name, k, nFiles));
    end
end
end

function efl = readEfl(TheSystem)
efl = TheSystem.MFE.GetOperandValue(ZOSAPI.Editors.MFE.MeritOperandType.EFFL, 0, 0, 0, 0, 0, 0, 0, 0);
end

function [valid, inYield, threshold, mtfRank] = selectYield(EFL, metric, yieldPct)
% Yield set: files whose MTF criterion is at or above the value that at least
% yieldPct % of the valid files meet. mtfRank: 1 = best MTF.
valid = isfinite(EFL) & EFL > 0 & isfinite(metric);
nValid = sum(valid);
mtfRank = nan(numel(EFL), 1);

if nValid == 0
    inYield = false(size(valid));
    threshold = NaN;
    return;
end

sortedMetric = sort(metric(valid), 'descend');
threshold = sortedMetric(max(1, ceil(yieldPct / 100 * nValid)));
inYield = valid & metric >= threshold;

[~, order] = sort(metric(valid), 'descend');
validIdx = find(valid);
mtfRank(validIdx(order)) = (1:nValid)';
end

function q = quantileSorted(x, p)
% Linear-interpolation quantiles of sorted data x (no toolbox needed).
n = numel(x);
if n == 1
    q = repmat(x, size(p));
    return;
end
pos = 1 + p * (n - 1);
q = interp1((1:n)', x(:), pos);
end

function t = fieldText(fieldList)
if isempty(fieldList)
    t = "all fields";
else
    t = "field(s) " + strjoin(string(fieldList), ", ");
end
end

function saveFigure(fig, fileName)
try
    exportgraphics(fig, fileName, 'Resolution', 300);   % R2020a+
catch
    saveas(fig, fileName);
end
end

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
