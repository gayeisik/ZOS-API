%% İKİ SİSTEM BÜYÜTME TOLERANS ANALİZİ (ZOS-API, MATLAB)
% NE YAPAR
%   Birlikte çalışan iki sistemin (1: MWIR objektif + dedektör, 2: OLED +
%   gündüz kanalına aktaran optik) Monte Carlo tolerans dosyalarından toplam
%   büyütmenin istatistiğini çıkarır. Her sistemde MTF'e göre %yieldPct yield
%   kümesine giren dosyalar seçilir; 1. sistemin yield kümesindeki HER EFL,
%   2. sistemin yield kümesindeki HER EFL ile eşleştirilir (tam çaprazlama,
%   n1 x n2 çift). Her çiftin büyütmesi hesaplanır ve büyütme isterine
%   (magTarget ± magTol) göre istatistik verilir. magTol birden fazla değer
%   olabilir (örn. [0.01 0.03]); her tolerans için sonuçlar ayrı verilir.
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
%   Geçme koşulu: |M - magTarget| <= magTol, magTol'daki her değer için ayrı.
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
%   Magnification_Results.xlsx  : "Summary", "By Tolerance" (her tolerans bir
%                                 sütun), "Percentiles", sistem başına dosya
%                                 tabloları, "Pairs" (tüm çiftler), "Settings"
%   Magnification_Data.mat      : tüm ham veriler
%   *.fig                       : grafiklerin düzenlenebilir MATLAB dosyaları
%   Analysis_Summary.txt        : özet ve sonuç / değerlendirme metni
%   Analysis_Log.txt            : Command Window çıktısının tamamı
%   Magnification_Report.pdf    : A4 rapor; özet, otomatik sonuç / değerlendirme,
%                                 yöntem, istatistiklerin anlamı ve her grafik için
%                                 "nasıl okunur / nasıl yorumlanır" açıklaması

clc;
clear;
close all;

% Command Window log; copied to the output folder at the end
diary off;
logFile = fullfile(tempdir, 'Magnification_Analysis_Log.txt');
if isfile(logFile)
    delete(logFile);
end
diary(logFile);

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
magTol          = [0.01 0.03];       % allowed deviation(s): |M - magTarget| <= magTol; one or more values
writePairsSheet = true;              % write every pair to the "Pairs" sheet (Excel limit 1048575 rows)
nHistBins       = 60;                % histogram bins of the magnification
% FFT MTF sample size is set in computeFileMtf (S_64x64 by default)

magConst = pixelScale * pixelPitch_um(2) / pixelPitch_um(1);   % M = magConst * |EFL_1| / |EFL_2|
magTol   = unique(abs(magTol(:)).');      % tightest tolerance first
nTol     = numel(magTol);
if nTol == 0 || any(magTol == 0)
    error('magTol must contain one or more positive values.');
end

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
mAll   = M(:);
nPairs = numel(mAll);

%% Statistics of the magnification (independent of the tolerance)
mMean = mean(mAll);
mStd  = std(mAll);
mMed  = median(mAll);
mMin  = min(mAll);
mMax  = max(mAll);

pctLevels = [0.135 1 2.5 5 25 50 75 95 97.5 99 99.865]';
pctValues = quantileSorted(sort(mAll), pctLevels / 100);

mNominal = magConst * sys(1).nominalEFL / sys(2).nominalEFL;

% Contribution of each system to the magnification error: var(ln M) = var(ln E1) + var(ln E2)
varLn = [var(log(E1)), var(log(E2))];
contribPct = 100 * varLn / sum(varLn);

%% Statistics for each tolerance
for t = nTol:-1:1
    tol(t) = tolStats(M, magTarget, magTol(t), mMean, mStd, magConst, ...
        sys(1).nominalEFL, sys(2).nominalEFL); %#ok<SAGROW>
end

% Pass rate of each file against all files of the other system (one column per tolerance)
for s = 1:2
    sys(s).passRate = nan(numel(sys(s).fileNames), nTol);
end
for t = 1:nTol
    sys(1).passRate(sys(1).inYield, t) = tol(t).passRate1;
    sys(2).passRate(sys(2).inYield, t) = tol(t).passRate2;
end

%% Command window summary
fprintf('\n===== MAGNIFICATION ANALYSIS =====\n');
fprintf('M = (|EFL_%s| / |EFL_%s|) x (%g x %g um / %g um) = %.6f x |EFL_1| / |EFL_2|\n', ...
    sys(1).name, sys(2).name, pixelScale, pixelPitch_um(2), pixelPitch_um(1), magConst);
fprintf('Nominal: |EFL_1| = %.4f mm, |EFL_2| = %.4f mm, M = %.5f\n', ...
    sys(1).nominalEFL, sys(2).nominalEFL, mNominal);
fprintf('Pairs: %d x %d = %d\n', n1, n2, nPairs);
fprintf('M mean %.5f  std %.5f  median %.5f  min %.5f  max %.5f\n', mMean, mStd, mMed, mMin, mMax);
fprintf('Contribution to magnification variance: %s %.1f %%, %s %.1f %%\n', ...
    sys(1).name, contribPct(1), sys(2).name, contribPct(2));

for t = 1:nTol
    T = tol(t);
    fprintf('\n--- Requirement M = %.4f +/- %.4f  [%.4f, %.4f] ---\n', magTarget, T.tol, T.low, T.high);
    fprintf('PASS: %d of %d pairs = %.2f %%   (below %.2f %%, above %.2f %%)\n', ...
        T.nPass, nPairs, T.passPct, T.belowPct, T.abovePct);
    fprintf('Cp %.3f  Cpk %.3f  normal-fit pass %.2f %%  centered-mean pass %.2f %%\n', ...
        T.Cp, T.Cpk, T.normPassPct, T.centeredPassPct);
    fprintf('Allowed |EFL_%s| for nominal |EFL_%s|: %.4f - %.4f mm\n', ...
        sys(2).name, sys(1).name, T.efl2Window(1), T.efl2Window(2));
    fprintf('Allowed |EFL_%s| for nominal |EFL_%s|: %.4f - %.4f mm\n', ...
        sys(1).name, sys(2).name, T.efl1Window(1), T.efl1Window(2));
    for s = 1:2
        fprintf('%s files matching at least one / all files of the other system: %.1f %% / %.1f %%\n', ...
            sys(s).name, T.anyMatch(s), T.allMatch(s));
    end
end

%% Output folder: "BÜYÜTME ANALİZİ - %90 YIELD" in the system 1 folder
outName = ['B' char(220) 'Y' char(220) 'TME ANAL' char(304) 'Z' char(304) sprintf(' - %%%g YIELD', yieldPct)];
outFolder = fullfile(sys(1).folder, outName);

if ~isfolder(outFolder)
    mkdir(outFolder);
end

colFail = [0.85 0.33 0.10];
colBlue = [0 0.45 0.74];
% One color per tolerance, tightest first
tolPalette = [0.10 0.50 0.10; 0.55 0.75 0.20; 0.93 0.69 0.13; 0.49 0.18 0.56; 0.30 0.75 0.93];
tolColors = tolPalette(mod(0:nTol - 1, size(tolPalette, 1)) + 1, :);
tolMarkers = {'o', 's', '^', 'd', 'v'};
eflLabel = @(s) sprintf('|EFL| %s (mm)', sys(s).name);
passTitle = strjoin(arrayfun(@(T) sprintf('\\pm%g: %.2f%%', T.tol, T.passPct), tol, 'UniformOutput', false), ',  ');

%% Figure 1: magnification histogram with the tolerance bands
fig1 = figure('Name', 'Magnification Distribution', 'Color', 'w', 'Position', [100 100 950 540]);
h = histogram(mAll, nHistBins, 'FaceColor', colBlue, 'DisplayName', sprintf('All pairs (%d)', nPairs));
hold on;
yl = ylim;
bands = gobjects(nTol, 1);
for t = nTol:-1:1   % widest band first, so the tighter bands are drawn on top of it
    bands(t) = patch([tol(t).low tol(t).high tol(t).high tol(t).low], [0 0 yl(2) yl(2)], tolColors(t, :), ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
        'DisplayName', sprintf('%.2f \\pm %g  (PASS %.2f%%)', magTarget, tol(t).tol, tol(t).passPct));
    xline(tol(t).low, '--', 'Color', tolColors(t, :), 'LineWidth', 1.2, 'HandleVisibility', 'off');
    xline(tol(t).high, '--', 'Color', tolColors(t, :), 'LineWidth', 1.2, 'HandleVisibility', 'off');
end
uistack(bands, 'bottom');
xline(magTarget, '-', 'Color', tolColors(1, :), 'LineWidth', 1.5, 'HandleVisibility', 'off');
hNom = xline(mNominal, ':', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', sprintf('Nominal M = %.4f', mNominal));
hMean = xline(mMean, '-.', 'Color', colFail, 'LineWidth', 1.5, 'DisplayName', sprintf('Mean M = %.4f', mMean));
hold off;
xMin = min(mMin, tol(end).low);
xMax = max(mMax, tol(end).high);
xlim([xMin xMax] + [-1 1] * 0.05 * (xMax - xMin));
xlabel('Magnification M');
ylabel('Number of pairs');
title(sprintf('Magnification of %s x %s pairs (%g%% MTF yield sets)  |  PASS %s', ...
    sys(1).name, sys(2).name, yieldPct, passTitle));
legend([h; bands; hNom; hMean], 'Location', 'best');
grid on;
saveFigure(fig1, fullfile(outFolder, 'Magnification_Histogram.png'));

%% Figure 2: EFL pair map, colored by the tightest tolerance each pair meets
fig2 = figure('Name', 'EFL Pair Map', 'Color', 'w', 'Position', [120 120 900 700]);
[G1, G2] = ndgrid(E1, E2);
passLevel = zeros(size(M));                 % 0 = fails every tolerance
for t = nTol:-1:1
    passLevel(tol(t).pass) = t;             % tightest passed tolerance wins
end
plot(G1(passLevel == 0), G2(passLevel == 0), '.', 'Color', colFail, 'MarkerSize', 6, ...
    'DisplayName', sprintf('Fail (outside \\pm%g)', tol(end).tol));
hold on;
for t = 1:nTol
    if t == 1
        label = sprintf('Pass \\pm%g', tol(t).tol);
    else
        label = sprintf('Pass \\pm%g only', tol(t).tol);
    end
    plot(G1(passLevel == t), G2(passLevel == t), '.', 'Color', tolColors(t, :), 'MarkerSize', 6, ...
        'DisplayName', label);
end
e1Line = linspace(min(E1), max(E1), 50);
plot(e1Line, magConst * e1Line / magTarget, '-k', 'LineWidth', 1.2, 'DisplayName', sprintf('M = %.2f', magTarget));
for t = 1:nTol
    plot(e1Line, magConst * e1Line / tol(t).low, '--', 'Color', tolColors(t, :), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('M = %.2f \\pm %g', magTarget, tol(t).tol));
    plot(e1Line, magConst * e1Line / tol(t).high, '--', 'Color', tolColors(t, :), 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
end
plot(sys(1).nominalEFL, sys(2).nominalEFL, 'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'y', 'DisplayName', 'Nominal');
hold off;
xlabel(eflLabel(1));
ylabel(eflLabel(2));
title(sprintf('EFL pairs (%d)  |  PASS %s', nPairs, passTitle));
legend('Location', 'best');
grid on;
saveFigure(fig2, fullfile(outFolder, 'EFL_Pair_Map.png'));

%% Figure 3: pass rate of each file against the other system, one series per tolerance
fig3 = figure('Name', 'Pass Rate by File', 'Color', 'w', 'Position', [140 140 950 700]);
eflSets = {E1, E2(:)};
for s = 1:2
    subplot(2, 1, s);
    hold on;
    for t = 1:nTol
        if s == 1
            rates = tol(t).passRate1;
        else
            rates = tol(t).passRate2;
        end
        plot(eflSets{s}, rates, tolMarkers{mod(t - 1, numel(tolMarkers)) + 1}, 'Color', tolColors(t, :), ...
            'MarkerFaceColor', tolColors(t, :), 'MarkerSize', 5, ...
            'DisplayName', sprintf('\\pm%g: %.0f%% match \\geq1, %.0f%% match all', ...
            tol(t).tol, tol(t).anyMatch(s), tol(t).allMatch(s)));
    end
    xline(sys(s).nominalEFL, ':k', 'Nominal', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    hold off;
    ylim([-5 105]);
    xlabel(eflLabel(s));
    ylabel(sprintf('Pass rate with %s files (%%)', sys(3 - s).name));
    title(sprintf('%s files: pass rate against the %s yield set', sys(s).name, sys(3 - s).name));
    legend('Location', 'best');
    grid on;
    box on;
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

tolTag = @(x) strrep(sprintf('%g', x), '.', 'p');   % 0.01 -> "0p01" for Excel column names

% Summary: values that do not depend on the tolerance
summaryNames = [ ...
    "Pairs (n1 x n2)"; "Mean M"; "Std M"; "Median M"; "Min M"; "Max M"; "Nominal M"; ...
    "Variance contribution " + sys(1).name; "Variance contribution " + sys(2).name; ...
    "Nominal |EFL| " + sys(1).name; "Nominal |EFL| " + sys(2).name];
summaryValues = [nPairs; mMean; mStd; mMed; mMin; mMax; mNominal; contribPct(1); contribPct(2); ...
    sys(1).nominalEFL; sys(2).nominalEFL];
summaryUnits = ["count"; ""; ""; ""; ""; ""; ""; "%"; "%"; "mm"; "mm"];
writetable(table(summaryNames, summaryValues, summaryUnits, ...
    'VariableNames', {'Parameter', 'Value', 'Unit'}), outputFile, 'Sheet', 'Summary');

% By Tolerance: one column per tolerance
byTolNames = [ ...
    "Tolerance (+/-)"; "Lower limit"; "Upper limit"; "Passing pairs"; "Pass rate"; ...
    "Below lower limit"; "Above upper limit"; "Cp"; "Cpk"; "Normal-fit pass rate"; ...
    "Pass rate if mean M centered at target"; ...
    sys(1).name + " files matching >= 1 " + sys(2).name + " file"; ...
    sys(1).name + " files matching all " + sys(2).name + " files"; ...
    sys(2).name + " files matching >= 1 " + sys(1).name + " file"; ...
    sys(2).name + " files matching all " + sys(1).name + " files"; ...
    "Allowed |EFL| " + sys(2).name + " min (nominal " + sys(1).name + ")"; ...
    "Allowed |EFL| " + sys(2).name + " max (nominal " + sys(1).name + ")"; ...
    "Allowed |EFL| " + sys(1).name + " min (nominal " + sys(2).name + ")"; ...
    "Allowed |EFL| " + sys(1).name + " max (nominal " + sys(2).name + ")"];
byTolUnits = [""; ""; ""; "count"; "%"; "%"; "%"; ""; ""; "%"; "%"; "%"; "%"; "%"; "%"; "mm"; "mm"; "mm"; "mm"];
byTol = table(byTolNames, byTolUnits, 'VariableNames', {'Parameter', 'Unit'});
for t = 1:nTol
    T = tol(t);
    byTol.(['Tol_' tolTag(T.tol)]) = [T.tol; T.low; T.high; T.nPass; T.passPct; T.belowPct; T.abovePct; ...
        T.Cp; T.Cpk; T.normPassPct; T.centeredPassPct; T.anyMatch(1); T.allMatch(1); T.anyMatch(2); ...
        T.allMatch(2); T.efl2Window(1); T.efl2Window(2); T.efl1Window(1); T.efl1Window(2)];
end
writetable(byTol, outputFile, 'Sheet', 'By Tolerance');

writetable(table(pctLevels, pctValues, 'VariableNames', {'Percentile', 'Magnification'}), ...
    outputFile, 'Sheet', 'Percentiles');

for s = 1:2
    fileTable = table(sys(s).fileNames, sys(s).EFL, sys(s).mtfMetric, sys(s).mtfRank, sys(s).inYield, ...
        'VariableNames', {'FileName', 'EFL_abs_mm', 'MTF_criterion', 'MTF_rank', 'InYieldSet'});
    for t = 1:nTol
        fileTable.(sprintf('PassRate_with_%s_tol_%s_pct', sys(3 - s).name, tolTag(magTol(t)))) = sys(s).passRate(:, t);
    end
    writetable(fileTable, outputFile, 'Sheet', char(sys(s).name + " Files"));
end

if writePairsSheet
    if nPairs <= 1048575
        [I, J] = ndgrid(1:n1, 1:n2);
        pairTable = table(names1(I(:)), names2(J(:)), E1(I(:)), E2(J(:)).', mAll, ...
            'VariableNames', {char("File_" + sys(1).name), char("File_" + sys(2).name), ...
            char("EFL_" + sys(1).name + "_mm"), char("EFL_" + sys(2).name + "_mm"), 'Magnification'});
        for t = 1:nTol
            pairTable.(['Pass_tol_' tolTag(magTol(t))]) = tol(t).pass(:);
        end
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
    "Magnification formula"; "Target magnification"; "Tolerances (+/-)"];
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
    string(magTarget); strjoin(string(magTol), ", ")];
writetable(table(settingNames, settingValues, 'VariableNames', {'Parameter', 'Value'}), ...
    outputFile, 'Sheet', 'Settings');

%% Save raw data
save(fullfile(outFolder, 'Magnification_Data.mat'), 'sys', 'E1', 'E2', 'names1', 'names2', 'M', 'tol', ...
    'magConst', 'magTarget', 'magTol', 'yieldPct', 'mtfFreqFrac', 'pixelPitch_um', 'nPixels', 'pixelScale');

%% PDF report (A4): summary, assessment, method, interpretation, one page per figure
pdfFile = fullfile(outFolder, 'Magnification_Report.pdf');

% Remove old report files (single and per-page)
oldPages = dir(fullfile(outFolder, 'Magnification_Report_page*.pdf'));
for i = 1:numel(oldPages)
    delete(fullfile(outFolder, oldPages(i).name));
end
if isfile(pdfFile)
    delete(pdfFile);
end

% If the old report could not be deleted (e.g. open in a PDF viewer),
% save the new one under a time-stamped name
if isfile(pdfFile)
    pdfFile = fullfile(outFolder, ['Magnification_Report_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.pdf']);
    fprintf('\nMagnification_Report.pdf is open in another program; new report: %s\n', pdfFile);
end

n1Name = char(sys(1).name);
n2Name = char(sys(2).name);
reqTextOf = @(T) sprintf('%.4f ± %g  (%.4f – %.4f)', magTarget, T.tol, T.low, T.high);
tolListText = strjoin(arrayfun(@(x) sprintf('±%g', x), magTol, 'UniformOutput', false), ' ve ');
reqText = sprintf('%.4f, tolerans %s', magTarget, tolListText);
afocalWord = ["Hayır", "Evet"];

% --- Page: summary -------------------------------------------------------
summaryLines = {
    sprintf('Tarih:                 %s', char(datetime('now', 'Format', 'dd.MM.yyyy HH:mm')))
    sprintf('Büyütme isteri:        M = %s', reqText)
    sprintf('Formül:                M = (%g × %.2f µm / %.2f µm) × |EFL_%s| / |EFL_%s|', ...
        pixelScale, pixelPitch_um(2), pixelPitch_um(1), n1Name, n2Name)
    sprintf('                         = %.6f × |EFL_%s| / |EFL_%s|', magConst, n1Name, n2Name)
    sprintf('Yield seviyesi:        %%%g  (MTF, Nyquist × %g)', yieldPct, mtfFreqFrac)
    ''
    };
for s = 1:2
    summaryLines = [summaryLines; {
        sprintf('SİSTEM %d – %s', s, sys(s).name)
        sprintf('  Klasör:              %s', sys(s).folder)
        sprintf('  Nominal dosya:       %s   |EFL| = %.4f mm', sys(s).nominalName, sys(s).nominalEFL)
        sprintf('  Afokal mod:          %s', afocalWord(sys(s).isAfocal + 1))
        sprintf('  Piksel:              %d × %.2f µm', nPixels(s), sys(s).pitch)
        sprintf('  MTF ölçütü:          %s, %.2f cyc/mm, eşik %.4f', fieldText(mtfFields{s}), sys(s).mtfFreq, sys(s).threshold)
        sprintf('  Dosyalar:            %d MC, %d geçerli, %d yield kümesinde', ...
            numel(sys(s).fileNames), sum(sys(s).valid), sum(sys(s).inYield))
        sprintf('  |EFL| yield kümesi:  ort. %.4f mm, std %.4f mm, %.4f – %.4f mm', ...
            mean(sys(s).EFL(sys(s).inYield)), std(sys(s).EFL(sys(s).inYield)), ...
            min(sys(s).EFL(sys(s).inYield)), max(sys(s).EFL(sys(s).inYield)))
        ''}]; %#ok<AGROW>
end
summaryLines = [summaryLines; {
    'SONUÇLAR (toleranstan bağımsız)'
    sprintf('  Çift sayısı:         %d × %d = %d', n1, n2, nPairs)
    sprintf('  Nominal M:           %.5f', mNominal)
    sprintf('  Ortalama M:          %.5f      Std M: %.5f', mMean, mStd)
    sprintf('  Medyan M:            %.5f      Min / Max: %.5f / %.5f', mMed, mMin, mMax)
    sprintf('  %%1 – %%99 aralığı:    %.5f – %.5f', pctValues(pctLevels == 1), pctValues(pctLevels == 99))
    sprintf('  Varyans katkısı:     %s %%%.1f, %s %%%.1f', n1Name, contribPct(1), n2Name, contribPct(2))
    ''
    }];
for t = 1:nTol
    T = tol(t);
    summaryLines = [summaryLines; {
        sprintf('TOLERANS ±%g  (M = %.4f – %.4f)', T.tol, T.low, T.high)
        sprintf('  Geçen çift:          %d   -> GEÇME ORANI %%%.2f', T.nPass, T.passPct)
        sprintf('  Alt / üst dışında:   %%%.2f / %%%.2f', T.belowPct, T.abovePct)
        sprintf('  Cp / Cpk:            %.3f / %.3f', T.Cp, T.Cpk)
        sprintf('  Normal dağ. tahmini: %%%.2f    Ortalama merkezde: %%%.2f', T.normPassPct, T.centeredPassPct)
        sprintf('  İzinli |EFL_%s|:    %.4f – %.4f mm  (nominal %s ile)', n2Name, T.efl2Window(1), T.efl2Window(2), n1Name)
        sprintf('  İzinli |EFL_%s|:    %.4f – %.4f mm  (nominal %s ile)', n1Name, T.efl1Window(1), T.efl1Window(2), n2Name)
        sprintf('  %s dosyaları:      en az bir %s ile %%%.1f, tümüyle %%%.1f', n1Name, n2Name, T.anyMatch(1), T.allMatch(1))
        sprintf('  %s dosyaları:      en az bir %s ile %%%.1f, tümüyle %%%.1f', n2Name, n1Name, T.anyMatch(2), T.allMatch(2))
        ''}]; %#ok<AGROW>
end

failedText = {};
for s = 1:2
    failed = sys(s).fileNames(~sys(s).valid);
    if ~isempty(failed)
        shown = failed(1:min(end, 5));
        failedText = [failedText; {sprintf('%s: hesaplanamayan dosyalar (%d): %s', ...
            sys(s).name, numel(failed), strjoin(shown, ', '))}]; %#ok<AGROW>
        if numel(failed) > numel(shown)
            failedText{end} = [failedText{end} sprintf(' ... (+%d)', numel(failed) - numel(shown))];
        end
    end
end
if ~isempty(failedText)
    summaryLines = [summaryLines; failedText];
end

summaryFig = textPage('BÜYÜTME TOLERANS ANALİZ RAPORU', summaryLines, 'Courier New', 8 - 0.5 * (nTol > 2));

% --- Page: assessment (generated from the results) -----------------------
assess = {'GENEL SONUÇ'};
for t = 1:nTol
    T = tol(t);
    assess{end + 1, 1} = sprintf(['•  ±%g (%s): %d çiftin %d''i isteri sağlıyor, geçme oranı %%%.2f. ' ...
        'Yield kümesinden rastgele seçilen bir %s ile bir %s birlikte kullanıldığında isteri sağlama ' ...
        'olasılığı yaklaşık %%%.1f''dir.'], T.tol, reqTextOf(T), nPairs, T.nPass, T.passPct, ...
        n1Name, n2Name, T.passPct); %#ok<SAGROW>
end
if nTol > 1
    assess{end + 1, 1} = sprintf(['Tolerans ±%g''den ±%g''ye daraltıldığında geçme oranı %%%.2f''den ' ...
        '%%%.2f''ye düşüyor. Dağılım aynı kaldığı için aradaki fark, M''si iki bant arasında kalan çiftlerdir.'], ...
        tol(end).tol, tol(1).tol, tol(end).passPct, tol(1).passPct);
end

assess = [assess; {''; 'NOMİNAL TASARIM'}];
nomDev = mNominal - magTarget;
nomText = sprintf('Nominal büyütme M = %.5f, hedeften sapma %+.5f.', mNominal, nomDev);
for t = 1:nTol
    if abs(nomDev) > tol(t).tol
        nomText = [nomText sprintf(' ±%g için nominal bile istere uymuyor; önce nominal EFL oranı düzeltilmeli.', ...
            tol(t).tol)]; %#ok<AGROW>
    else
        nomText = [nomText sprintf(' ±%g bandının %%%.0f''i nominalde kullanılmış.', ...
            tol(t).tol, 100 * abs(nomDev) / tol(t).tol)]; %#ok<AGROW>
    end
end
assess{end + 1, 1} = nomText;

assess = [assess; {''; 'ORTALAMA KAYMASI MI, SAÇILIM MI?'}];
shift = mMean - magTarget;
assess{end + 1, 1} = sprintf('Ortalama M = %.5f, hedeften sapma %+.5f. Std M = %.5f.', mMean, shift, mStd);
for t = 1:nTol
    T = tol(t);
    if abs(shift) > 0.2 * T.tol
        assess{end + 1, 1} = sprintf(['•  ±%g: kayma bandın %%%.0f''i kadar, belirgin. Ortalama hedefe ' ...
            'getirilirse (nominal EFL oranı düzeltilerek) geçme oranı %%%.2f''den %%%.2f''ye çıkar.'], ...
            T.tol, 100 * abs(shift) / T.tol, T.passPct, T.centeredPassPct); %#ok<SAGROW>
    else
        assess{end + 1, 1} = sprintf(['•  ±%g: kayma küçük (bandın %%%.0f''i). Ortalamayı merkeze almak geçme ' ...
            'oranını yalnızca %%%.2f''den %%%.2f''ye çıkarır; kayıp esas olarak saçılımdan geliyor.'], ...
            T.tol, 100 * abs(shift) / T.tol, T.passPct, T.centeredPassPct); %#ok<SAGROW>
    end
end

assess = [assess; {''; 'SÜREÇ YETERLİLİĞİ (Cpk)'}];
for t = 1:nTol
    T = tol(t);
    if T.Cpk >= 1.33
        cpkText = 'yeterli (≥ 1.33): dağılım banda rahatça sığıyor';
    elseif T.Cpk >= 1
        cpkText = 'sınırda (1 – 1.33): dağılım banda sığıyor ama pay az';
    else
        cpkText = 'yetersiz (< 1): dağılımın kuyrukları bandın dışına taşıyor';
    end
    assess{end + 1, 1} = sprintf(['•  ±%g: Cp = %.3f, Cpk = %.3f -> %s. Normal dağılım tahmini %%%.2f ' ...
        '(ölçülen %%%.2f).'], T.tol, T.Cp, T.Cpk, cpkText, T.normPassPct, T.passPct); %#ok<SAGROW>
end
assess{end + 1, 1} = 'Normal tahmin ile ölçülen oran belirgin farklıysa dağılım normal değildir; ölçülen oran esas alınmalı.';

assess = [assess; {''; 'HANGİ SİSTEM BASKIN?'}];
[~, dom] = max(contribPct);
assess{end + 1, 1} = sprintf(['Büyütme saçılımının %%%.1f''i %s, %%%.1f''i %s sisteminin EFL saçılımından geliyor ' ...
    '(toleranstan bağımsız). Geçme oranını artırmak için EFL toleranslarını sıkılaştırmada öncelik %s sistemindedir.'], ...
    contribPct(1), n1Name, contribPct(2), n2Name, sys(dom).name);

assess = [assess; {''; 'İZİN VERİLEN EFL ARALIKLARI'}];
for t = 1:nTol
    T = tol(t);
    assess{end + 1, 1} = sprintf(['•  ±%g: %s nominal (%.4f mm) ise %s |EFL| %.4f – %.4f mm; %s nominal ' ...
        '(%.4f mm) ise %s |EFL| %.4f – %.4f mm olmalı.'], T.tol, n1Name, sys(1).nominalEFL, n2Name, ...
        T.efl2Window(1), T.efl2Window(2), n2Name, sys(2).nominalEFL, n1Name, T.efl1Window(1), T.efl1Window(2)); %#ok<SAGROW>
end

assess = [assess; {''; 'EŞLEŞTİRME (SEÇİCİ MONTAJ) POTANSİYELİ'}];
for t = 1:nTol
    T = tol(t);
    assess{end + 1, 1} = sprintf(['•  ±%g: %s dosyalarının %%%.1f''i en az bir %s ile, %%%.1f''i tümüyle; ' ...
        '%s dosyalarının %%%.1f''i en az bir %s ile, %%%.1f''i tümüyle isteri sağlıyor.'], T.tol, ...
        n1Name, T.anyMatch(1), n2Name, T.allMatch(1), n2Name, T.anyMatch(2), n1Name, T.allMatch(2)); %#ok<SAGROW>
end
assess{end + 1, 1} = sprintf(['"Tümüyle" oranı, diğer sistemin saçılımı banda sığmıyorsa düşük çıkar (saçılımı büyük ' ...
    'olan sistemin karşısındaki dosyalar "tümü" testini geçemez). "En az bir" oranı "tümü" oranından belirgin ' ...
    'büyükse, EFL ölçülüp sistemler eşleştirilerek (seçici montaj) geçme oranı artırılabilir.']);

if min(n1, n2) < 50
    assess = [assess; {''; 'UYARI – ÖRNEKLEM BÜYÜKLÜĞÜ'}];
    assess{end + 1, 1} = sprintf(['Yield kümelerinde %d ve %d dosya var. Çift sayısı (%d) büyük görünse de ' ...
        'çiftler birbirinden bağımsız değildir (her dosya %d veya %d çiftte tekrar kullanılır); sonucun ' ...
        'güvenilirliği çift sayısıyla değil dosya sayısıyla sınırlıdır. Güvenilir bir oran için sistem başına ' ...
        'en az 50 – 100 MC dosyası önerilir.'], n1, n2, nPairs, n2, n1);
end

assessFig = textPage('1. SONUÇ VE DEĞERLENDİRME', assess, 'Arial', 10 - 0.5 * (nTol > 1));

% --- Page: method --------------------------------------------------------
method = {
    'AMAÇ'
    sprintf(['%s ve %s birlikte çalışır: %s görüntüsü %s üzerinde gösterilir ve gündüz kanalına aktarılır. ' ...
     'Toplam büyütmenin %s olması istenir. Her iki sistemin üretim toleransları EFL''lerini değiştirir; bu ' ...
     'analiz, toleranslı iki sistem birleştirildiğinde büyütmenin isteri hangi olasılıkla sağladığını bulur.'], ...
     n1Name, n2Name, n1Name, n2Name, reqText)
    ''
    'ADIM 1 – DOSYALAR'
    sprintf(['Her sistem için ayrı bir Monte Carlo klasörü seçilir. Adı %s ile başlayan .zmx dosyaları MC ' ...
     'dosyalarıdır (her biri olası bir üretilmiş sistem); adı MC_ ile başlamayan dosya nominal tasarımdır.'], filePrefix)
    ''
    'ADIM 2 – EFL VE MTF'
    ['Her dosya OpticStudio''da açılır. EFL, EFFL operandıyla okunur ve her zaman mutlak değeri (|EFL|) ' ...
     'kullanılır. MTF polikromatik FFT MTF''dir (tüm dalga boyları, tüm fieldler, T ve S).']
    sprintf(['MTF ölçütü her sistemin kendi Nyquist frekansının %g katında alınır: %s için %.2f cyc/mm ' ...
     '(%.2f µm piksel), %s için %.2f cyc/mm (%.2f µm piksel). Dosya başına tek sayı, seçilen fieldlerde T ' ...
     've S MTF''lerinin en düşüğüdür. Afokal dosyalarda Zemax MTF''i cyc/mrad verir; frekans ' ...
     'ν[cyc/mm] = ν[cyc/mrad] × 1000 / |EFL_nominal| ile çevrilir.'], mtfFreqFrac, n1Name, sys(1).mtfFreq, ...
     sys(1).pitch, n2Name, sys(2).mtfFreq, sys(2).pitch)
    ''
    'ADIM 3 – YIELD KÜMESİ'
    sprintf(['Her sistemde dosyalar MTF ölçütüne göre sıralanır ve en iyi %%%g seçilir. Yalnızca bu dosyalar ' ...
     'çaprazlamaya girer: MTF şartını sağlamayan sistemler zaten üretimde elenecektir.'], yieldPct)
    ''
    'ADIM 4 – ÇAPRAZLAMA VE BÜYÜTME'
    sprintf(['%s yield kümesindeki her dosya, %s yield kümesindeki her dosyayla eşleştirilir: %d × %d = %d çift. ' ...
     'Bu, olası tüm (%s, %s) ikililerinin listesidir; bir %s her %s ile denenir.'], ...
     n1Name, n2Name, n1, n2, nPairs, n1Name, n2Name, n1Name, n2Name)
    sprintf(['Bir hedef, %s''da |EFL_1| × tanθ büyüklüğünde görüntü oluşturur; bu %.2f µm''lik dedektör ' ...
     'pikselleriyle örneklenir ve her dedektör pikseli OLED''de %g piksele (%.2f µm) basılır. OLED görüntüsü ' ...
     '|EFL_2| ile açıya çevrilir. Çıkış açısının giriş açısına oranı büyütmedir:'], ...
     n1Name, pixelPitch_um(1), pixelScale, pixelPitch_um(2))
    sprintf('        M = (%g × %.2f / %.2f) × |EFL_1| / |EFL_2| = %.6f × |EFL_1| / |EFL_2|', ...
        pixelScale, pixelPitch_um(2), pixelPitch_um(1), magConst)
    sprintf(['Bir çift, her tolerans için ayrı değerlendirilir: |M − %.2f| ≤ tolerans ise o toleransta geçer ' ...
     '(bu analizde %s). Büyütme dağılımı toleranstan bağımsızdır; yalnızca kabul bandı değişir.'], magTarget, tolListText)
    ''
    'NEDEN n1 × n2 ÇİFT (n! DEĞİL)?'
    ['n1 × n2, birbirinden farklı tüm ikililerin sayısıdır. n! ise n sistemin hepsini aynı anda birebir ' ...
     'dağıtmanın yol sayısıdır; her dağıtım bu ikililerden n tanesini seçer ve yeni bir ikili üretmez. Her ' ...
     'ikili tüm dağıtımlarda eşit sayıda yer aldığından, geçme oranı iki yolla da aynı çıkar.']
    };
methodFig = textPage('2. ANALİZ YÖNTEMİ', method, 'Arial', 10);

% --- Page: meaning of the statistics -------------------------------------
meaning = {
    'GEÇME ORANI'
    ['Kriteri sağlayan çiftlerin tüm çiftlere oranı (her tolerans için ayrı). Yield kümelerinden rastgele ' ...
     'seçilen iki sistemin birlikte isteri sağlama olasılığının tahminidir. Alt / üst sınır dışında kalan ' ...
     'oranlar, büyütmenin hangi yönde sorun çıkardığını gösterir.']
    ''
    'ORTALAMA, STD, YÜZDELİKLER'
    ['Ortalama M''nin hedeften farkı sistematik kaymadır; std saçılımdır. Bunlar toleranstan bağımsızdır. ' ...
     '%1 – %99 aralığı, sistemlerin %98''inin içinde kaldığı büyütme aralığıdır; bu aralık bir tolerans ' ...
     'bandının içindeyse tasarım o tolerans için güvenlidir. Tüm yüzdelikler Excel''de "Percentiles" sayfasındadır.']
    ''
    'Cp VE Cpk'
    ['Cp = bant genişliği / (6σ): dağılımın genişliği banda sığıyor mu? Cpk = ortalamadan en yakın sınıra ' ...
     'uzaklık / (3σ): kayma da hesaba katılır. Cpk ≥ 1.33 yeterli, 1 – 1.33 sınırda, < 1 yetersiz kabul ' ...
     'edilir. Cp büyük ama Cpk küçükse sorun kaymadır, ikisi de küçükse sorun saçılımdır. Bant daraldıkça ' ...
     'Cp ve Cpk aynı oranda küçülür. Normal dağılım varsayar; az dosyayla dikkatli yorumlanmalıdır.']
    ''
    'ORTALAMA MERKEZDE GEÇME ORANI'
    ['Tüm büyütmeler, ortalama tam hedefe gelecek şekilde ölçeklenirse elde edilecek geçme oranı. Bu, ' ...
     'nominal EFL oranını düzelterek kazanılabilecek en fazla iyileşmeyi gösterir. Mevcut oranla farkı ' ...
     'küçükse nominal düzeltme işe yaramaz, toleranslar sıkılaşmalıdır. Dar toleransta aynı kayma daha ' ...
     'büyük kayıp yaratır.']
    ''
    'VARYANS KATKISI'
    ['M = sabit × EFL_1 / EFL_2 olduğundan ln M = ln EFL_1 − ln EFL_2 ve var(ln M) = var(ln EFL_1) + ' ...
     'var(ln EFL_2). Her terimin toplama oranı, o sistemin büyütme saçılımına katkısıdır (yaklaşık olarak ' ...
     'EFL''nin göreli std''sinin karesi). Katkısı büyük sistemin EFL toleransı önceliklidir.']
    ''
    'DOSYA BAZINDA GEÇME ORANI'
    ['Bir dosyanın diğer sistemin dosyalarının yüzde kaçıyla isteri sağladığı. %100: hangi eşle birleşirse ' ...
     'birleşsin çalışır. %0: hiçbiriyle çalışmaz. Oranların EFL''ye göre nasıl düştüğü, o sistem için ' ...
     'kabul edilebilir EFL sınırını verir. Excel''de sistem başına dosya tablolarında her tolerans için verilmiştir.']
    ''
    'İZİN VERİLEN EFL ARALIĞI'
    ['Bir sistem tam nominal üretilirse, isteri sağlamak için diğer sistemin |EFL|''sinin kalması gereken ' ...
     'aralık. Diğer sistem de toleranslı olduğundan gerçek pay bu aralıktan daha dardır.']
    };
meaningFig = textPage('3. İSTATİSTİKLERİN ANLAMI', meaning, 'Arial', 10);

% --- Figure pages --------------------------------------------------------
plotFigs = gobjects(4, 1);

plotFigs(1) = figurePage(fig1, '4. BÜYÜTME DAĞILIMI', {
    'NASIL OKUNUR'
    sprintf(['•  Her çubuk, büyütmesi o aralıkta olan %s – %s çifti sayısıdır. Renkli bantlar tolerans ' ...
     'isterleridir (%s); iç içe bantlarda en koyu ve dar olan en sıkı toleranstır. Bir bandın içindeki ' ...
     'çubuklar o toleransta geçen çiftlerdir; her toleransın geçme oranı lejantta verilmiştir.'], n1Name, n2Name, tolListText)
    ['•  Siyah noktalı çizgi nominal M, turuncu kesikli çizgi tüm çiftlerin ortalama M''sidir. İkisi ' ...
     'farklıysa toleranslar büyütmeyi sistematik olarak bir yöne itiyor demektir.']
    'NASIL YORUMLANIR'
    ['•  Dağılım bir bandın ortasında ve o banttan dar: tasarım o tolerans için güvenli. Geniş bantta ' ...
     'güvenli olup dar bantta taşıyorsa, dar isterin karşılanması için saçılım azaltılmalıdır.']
    ['•  Dağılım bir kenara yığılmış veya tek taraftan taşıyor: sorun kaymadır; nominal EFL oranı ' ...
     'düzeltilerek çözülür ("ortalama merkezde geçme oranı" kazancı gösterir).']
    ['•  Dağılım bandın iki tarafından da taşıyor: sorun saçılımdır; EFL toleransları sıkılaştırılmalı ' ...
     '(öncelik varyans katkısı büyük olan sistemde) veya seçici montaj yapılmalıdır.']
    }, false);

plotFigs(2) = figurePage(fig2, '5. EFL ÇİFT HARİTASI', {
    'NASIL OKUNUR'
    sprintf(['•  Her nokta bir çifttir: yatay eksen %s |EFL|, dikey eksen %s |EFL|. Nokta rengi, çiftin ' ...
     'sağladığı en sıkı toleranstır (koyu yeşil = en dar bantta geçti, açık renk = yalnızca daha geniş ' ...
     'bantta geçti); kırmızı = hiçbir toleransı sağlamadı. Sarı yıldız nominal tasarımdır.'], n1Name, n2Name)
    sprintf(['•  Düz siyah çizgi M = %.2f (ideal oran); her toleransın sınırları kendi renginde kesikli ' ...
     'çizgilerdir. Sol üstte M küçük, sağ altta M büyüktür.'], magTarget)
    'NASIL YORUMLANIR'
    sprintf(['•  Nokta bulutu yatayda genişse saçılım %s''dan, dikeyde genişse %s''den gelir. Bulut kabul ' ...
     'bölgesinin bir tarafına kaymışsa nominal oran düzeltilmelidir.'], n1Name, n2Name)
    ['•  Kabul bölgesi çizgilerine paralel bir şerit içinde kalan EFL çiftleri birlikte çalışır. Bulut ' ...
     'bölgeden taşıyorsa, EFL ölçülüp sistemler bu harita üzerinden eşleştirilebilir (seçici montaj).']
    }, true);

plotFigs(3) = figurePage(fig3, '6. DOSYA BAZINDA GEÇME ORANI', {
    'NASIL OKUNUR'
    ['•  Her nokta bir dosyadır; her tolerans ayrı renk ve işaretle gösterilir. Yatay eksen dosyanın ' ...
     '|EFL|''si, dikey eksen diğer sistemin yield kümesindeki dosyaların yüzde kaçıyla isteri sağladığıdır. ' ...
     'Noktalı çizgi nominal |EFL|''dir. Lejantta her tolerans için "en az birine" ve "hepsine" uyan ' ...
     'dosyaların oranı verilmiştir.']
    'NASIL YORUMLANIR'
    ['•  %100 olan dosyalar diğer sistemin her dosyasıyla çalışır; %0 olanlar hiçbiriyle çalışmaz. ' ...
     'Genellikle nominal EFL çevresinde yüksek, uçlarda düşük bir tepe şekli görülür; dar toleransta tepe ' ...
     'daralır ve alçalır. Tepe nominalden kaymışsa nominal oran düzeltilmelidir.']
    ['•  Oranın hızla düştüğü EFL değerleri, o sistem için kabul edilebilir EFL sınırlarıdır; EFL toleransı ' ...
     'veya EFL''ye göre sınıflama (seçici montaj) için kullanılabilir.']
    ['•  "Hepsine uyan" oranı, diğer sistemin saçılımı banda sığmıyorsa sıfıra iner; saçılımı büyük olan ' ...
     'sistemin karşısındaki dosyalar bu testi geçemez.']
    }, false);

plotFigs(4) = figurePage(fig4, '7. EFL DAĞILIMLARI', {
    'NASIL OKUNUR'
    sprintf(['•  Gri: tüm geçerli MC dosyalarının |EFL| dağılımı. Mavi: MTF''e göre %%%g yield kümesi ' ...
     '(çaprazlamada kullanılan dosyalar). Noktalı çizgi nominal |EFL|''dir. Bu grafik toleranstan bağımsızdır.'], yieldPct)
    'NASIL YORUMLANIR'
    ['•  Mavi histogram griden belirgin farklıysa MTF ile EFL ilişkilidir: MTF''i kötü olup elenen ' ...
     'dosyalar EFL''de de uçtadır.']
    ['•  Dağılımın nominale göre kayması, toleransların EFL''yi sistematik olarak değiştirdiğini gösterir; ' ...
     'bu kayma doğrudan büyütmenin ortalamasına yansır.']
    ['•  İki sistemin göreli genişliği (std / ortalama EFL) büyütme saçılımına katkılarını belirler; ' ...
     'sayısal pay Sonuç sayfasında verilmiştir.']
    }, false);

% --- Assemble pages: frame, header and page numbers ----------------------
allFigs = [summaryFig; assessFig; methodFig; meaningFig; plotFigs];
nPages = numel(allFigs);
headerText = sprintf('Büyütme Tolerans Analizi  |  %s: %s  |  %s: %s', ...
    n1Name, sys(1).nominalName, n2Name, sys(2).nominalName);

for p = 1:nPages
    addPageFrame(allFigs(p), headerText, p, nPages);
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
        print(allFigs(p), fullfile(outFolder, sprintf('Magnification_Report_page%d.pdf', p)), '-dpdf');
    end
    fprintf('Report saved as separate pages (Magnification_Report_page*.pdf).\n');
end

close(allFigs);

%% Text summary of the results (same text as the summary and assessment pages)
summaryFile = fullfile(outFolder, 'Analysis_Summary.txt');
fid = fopen(summaryFile, 'w', 'n', 'UTF-8');
if fid > 0
    fprintf(fid, '%s\r\n', 'BÜYÜTME TOLERANS ANALİZ ÖZETİ', '');
    fprintf(fid, '%s\r\n', summaryLines{:});
    fprintf(fid, '\r\n%s\r\n\r\n', 'SONUÇ VE DEĞERLENDİRME');
    fprintf(fid, '%s\r\n', assess{:});
    fclose(fid);
end

%% List of saved outputs; the Command Window log is copied last
fprintf('\nResults saved to:\n%s\n', outFolder);
saved = dir(outFolder);
saved = saved(~[saved.isdir]);
fprintf('  %s\n', saved.name);
fprintf('  Analysis_Log.txt\n');

diary off;
try
    copyfile(logFile, fullfile(outFolder, 'Analysis_Log.txt'));
catch logErr
    fprintf('Command Window log could not be copied: %s\n', logErr.message);
end


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

function T = tolStats(M, magTarget, tolValue, mMean, mStd, magConst, nominalEFL1, nominalEFL2)
% Pass / fail statistics of the magnification matrix M (n1 x n2) for one
% tolerance: |M - magTarget| <= tolValue.
T.tol  = tolValue;
T.low  = magTarget - tolValue;
T.high = magTarget + tolValue;

mAll = M(:);
T.pass     = abs(M - magTarget) <= tolValue + 1e-12;
T.nPass    = sum(T.pass(:));
T.passPct  = 100 * T.nPass / numel(mAll);
T.belowPct = 100 * mean(mAll < T.low - 1e-12);
T.abovePct = 100 * mean(mAll > T.high + 1e-12);

T.Cp  = (T.high - T.low) / (6 * mStd);
T.Cpk = min(T.high - mMean, mMean - T.low) / (3 * mStd);
normPhi = @(x) 0.5 * erfc(-x / sqrt(2));
T.normPassPct = 100 * (normPhi((T.high - mMean) / mStd) - normPhi((T.low - mMean) / mStd));

% If the mean magnification were moved to the target (nominal EFL correction)
T.centeredPassPct = 100 * mean(abs(mAll * (magTarget / mMean) - magTarget) <= tolValue + 1e-12);

% Allowed |EFL| of one system when the other one is nominal
T.efl2Window = magConst * nominalEFL1 ./ [T.high, T.low];
T.efl1Window = nominalEFL2 / magConst * [T.low, T.high];

% Pass rate of each file against all files of the other system
T.passRate1 = 100 * mean(T.pass, 2);       % n1 x 1
T.passRate2 = 100 * mean(T.pass, 1).';     % n2 x 1
T.anyMatch = [100 * mean(T.passRate1 > 0), 100 * mean(T.passRate2 > 0)];
T.allMatch = [100 * mean(T.passRate1 == 100), 100 * mean(T.passRate2 == 100)];
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

function fig = newA4Page()
% Blank A4 portrait page (21 x 29.7 cm). Invisible, so its size is not
% limited by the screen.
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
    'Position', [1 1 21 29.7], 'PaperUnits', 'centimeters', ...
    'PaperSize', [21 29.7], 'PaperPosition', [0 0 21 29.7]);
end

function fig = textPage(titleText, lines, fontName, fontSize)
% A4 page with a title and wrapped text (one cell per paragraph/line).
% Arial / Courier New are set explicitly: the default font has no Turkish glyphs in the PDF
fig = newA4Page();
annotation(fig, 'textbox', [0.08 0.89 0.84 0.05], 'String', titleText, ...
    'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'Interpreter', 'none', 'VerticalAlignment', 'top');
annotation(fig, 'textbox', [0.08 0.06 0.84 0.83], 'String', lines, ...
    'FontName', fontName, 'FontSize', fontSize, 'EdgeColor', 'none', ...
    'Interpreter', 'none', 'VerticalAlignment', 'top', 'FitBoxToText', 'off');
end

function fig = figurePage(srcFig, titleText, captionLines, squareAxes)
% A4 page: title, the axes of srcFig (copied with their legends, stacked
% top to bottom) and an explanation text below them.
fig = newA4Page();
annotation(fig, 'textbox', [0.08 0.89 0.84 0.05], 'String', titleText, ...
    'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'Interpreter', 'none', 'VerticalAlignment', 'top');

srcAx = findobj(srcFig, 'Type', 'axes');
yPos = arrayfun(@(a) a.Position(2), srcAx);
[~, order] = sort(yPos, 'descend');
srcAx = srcAx(order);

% Plot area: y = 0.33 ... 0.87 of the page, the text below it
nAx = numel(srcAx);
areaBottom = 0.34;
areaTop = 0.86;
slot = (areaTop - areaBottom) / nAx;
for i = 1:nAx
    a = srcAx(i);
    if ~isempty(a.Legend)
        copied = copyobj([a.Legend, a], fig);
    else
        copied = copyobj(a, fig);
    end
    newAx = findobj(copied, 'Type', 'axes');
    newAx.Units = 'normalized';
    newAx.Position = [0.13, areaTop - i * slot + 0.05, 0.80, slot - 0.08];
    newAx.FontSize = 8;
    newAx.Title.FontSize = 9;
    if squareAxes
        pbaspect(newAx, [1 1 1]);
    end
end

annotation(fig, 'textbox', [0.08 0.05 0.84 0.25], 'String', captionLines, ...
    'FontName', 'Arial', 'FontSize', 9.5, 'EdgeColor', 'none', ...
    'Interpreter', 'none', 'VerticalAlignment', 'top', 'FitBoxToText', 'off');
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

function saveFigure(fig, fileName)
% PNG (300 dpi) and an editable MATLAB .fig with the same name
try
    exportgraphics(fig, fileName, 'Resolution', 300);   % R2020a+
catch
    saveas(fig, fileName);
end
[folder, name] = fileparts(fileName);
savefig(fig, fullfile(folder, [name '.fig']));
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
