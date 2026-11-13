%% EŞLEŞTİRMELİ MONTAJ ANALİZİ (ZOS-API, MATLAB)
% NE YAPAR
%   İki sistemin (1: MWIR, 2: kolimatör + OLED) MTF'e göre %yieldPct yield
%   kümelerindeki sistemler, büyütme isterini (magTarget ± magTol) sağlayacak
%   şekilde birbirleriyle eşleştirilir. Her sistem en fazla bir kez kullanılır.
%   Soru: yield kümelerindeki N_MWIR ve N_OLED sistemden en fazla kaç ürün
%   çıkar? Örneğin 400 MC dosyasının %90 yield kümesi 360 sistemse, 360 MWIR
%   ve 360 OLED'den en fazla kaç ürün çıktığı bulunur.
%   Başka bir scriptin çıktısı kullanılmaz; EFL ve MTF bu scriptte
%   OpticStudio ile hesaplanır.
%
% NASIL ÇALIŞIR
%   1. Kullanıcı önce 1. sistemin (MWIR), sonra 2. sistemin (OLED) Monte Carlo
%      klasörünü seçer. Adı filePrefix ("MC_T") ile başlayan .zmx dosyaları MC
%      dosyalarıdır; adı "MC_" ile başlamayan .zmx dosyası nominal tasarımdır
%      (birden fazlaysa listeden seçilir).
%   2. Her sistem için OpticStudio ile her MC dosyasının |EFL|'si (EFFL
%      operandı, HER ZAMAN mutlak değer) ve polikromatik FFT MTF'i (tüm dalga
%      boyları, tüm fieldler, T ve S; sistemin Nyquist'inin mtfFreqFrac
%      katında) hesaplanır. Afokal dosyalarda frekans
%      ν[cyc/mm] = ν[cyc/mrad] × 1000 / |EFL_nominal|  ile çevrilir.
%   3. Her sistemde MTF'e göre en iyi %yieldPct dosya seçilir (yield kümesi).
%   4. Her tolerans için en iyi eşleştirme bulunur.
%
% HESAPLAR
%   Büyütme:  M = (k × p_2 / p_1) × |EFL_1| / |EFL_2|
%     p_1 = dedektör piksel aralığı, p_2 = OLED piksel aralığı,
%     k   = bir dedektör pikselinin OLED'de kapladığı piksel sayısı (1024/512 = 2).
%   Bir çift, |M - magTarget| <= magTol ise uyumludur.
%   En iyi eşleştirme (maksimum eşleştirme):
%     Her MWIR için uyumlu OLED EFL'leri sürekli bir aralıktır ve bu aralık
%     MWIR EFL'si büyüdükçe kayar. Bu yüzden şu kural KESİN en büyük ürün
%     sayısını verir: MWIR'lar EFL'ye göre küçükten büyüğe alınır, her birine
%     uyumlu ve henüz kullanılmamış en küçük EFL'li OLED verilir.
%     Aynı ürün sayısını veren eşleştirmeler arasından toplam |M - hedef| en
%     küçük olan seçilir (matchpairs, R2019a+; yoksa yukarıdaki kural).
%   Eşleştirme verimi    = ürün sayısı / min(N_MWIR, N_OLED)
%   Toplam verim         = ürün sayısı / MC dosya sayısı (MTF yield + eşleştirme)
%   Rastgele montaj      = eşleştirme yapılmadan birebir rastgele dağıtımda
%                          beklenen ürün sayısı = min(N_MWIR, N_OLED) × uyumlu
%                          çiftlerin oranı (karşılaştırma için)
%   Ortalama merkezde    = ortalama büyütme hedefe kaydırılsaydı (nominal EFL
%                          oranı düzeltilerek) eşleştirmeyle çıkacak ürün sayısı
%
% AYARLAR (User settings)
%   sysName, pixelPitch_um, nPixels, pixelScale, filePrefix, yieldPct,
%   mtfFreqFrac, mtfFields, magTarget, magTol, nHistBins
%
% ÇIKTILAR (1. sistemin klasöründe "EŞLEŞTİRME ANALİZİ - %90 YIELD" klasörüne)
%   Products.png               : rastgele / eşleştirmeli / merkezlenmiş ürün sayısı
%   Matched_Pairs_Map.png      : eşleştirilen çiftler, EFL_1 - EFL_2 düzleminde
%   Unmatched_EFL.png          : açıkta kalan sistemlerin EFL'leri
%   Matched_Magnification.png  : eşleştirilen çiftlerin büyütme dağılımı
%   *.fig                      : grafiklerin düzenlenebilir MATLAB dosyaları
%   Matching_Results.xlsx      : "Summary" (her tolerans bir sütun),
%                                "Matched tol ..." (ürün listesi), "Unmatched
%                                tol ..." (açıkta kalanlar), sistem başına dosya
%                                tabloları (her dosyanın eşi), "Settings"
%   Matching_Report.pdf        : A4 rapor; özet, sonuç, yöntem, grafikler ve
%                                nasıl yorumlanacakları
%   Matching_Data.mat          : tüm ham veriler
%   Analysis_Summary.txt       : özet ve sonuç metni
%   Analysis_Log.txt           : Command Window çıktısının tamamı

clc;
clear;
close all;

% Command Window log; copied to the output folder at the end
diary off;
logFile = fullfile(tempdir, 'Matching_Analysis_Log.txt');
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
nHistBins       = 40;                % histogram bins
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


%% Yield sets
E1 = sys(1).EFL(sys(1).inYield);           % N_MWIR x 1 (absolute values)
E2 = sys(2).EFL(sys(2).inYield).';         % 1 x N_OLED
names1 = sys(1).fileNames(sys(1).inYield);
names2 = sys(2).fileNames(sys(2).inYield);
idx1 = find(sys(1).inYield);
idx2 = find(sys(2).inYield);
n1 = numel(E1);
n2 = numel(E2);
nMax = min(n1, n2);                                              % largest possible number of products
nMcRef = min(numel(sys(1).fileNames), numel(sys(2).fileNames));  % MC files per system

M = magConst * E1 ./ E2;                    % n1 x n2 magnification of every possible pair
mMean = mean(M(:));
mNominal = magConst * sys(1).nominalEFL / sys(2).nominalEFL;
Mcentered = M * (magTarget / mMean);       % mean moved to the target (nominal EFL ratio corrected)

%% Best matching for each tolerance
for t = nTol:-1:1
    R = struct();
    R.tol  = magTol(t);
    R.low  = magTarget - magTol(t);
    R.high = magTarget + magTol(t);

    pass = abs(M - magTarget) <= magTol(t) + 1e-12;
    R.compatPct = 100 * mean(pass(:));                  % compatible share of all possible pairs
    R.randomProducts = R.compatPct / 100 * nMax;        % expected with random one-to-one assembly

    [R.pairs, R.method] = bestMatching(M, pass, E1, E2, magTarget);
    R.nMatched = size(R.pairs, 1);
    R.matchPct = 100 * R.nMatched / nMax;
    R.overallPct = 100 * R.nMatched / nMcRef;
    R.matchM = M(sub2ind(size(M), R.pairs(:, 1), R.pairs(:, 2)));
    R.unmatched1 = setdiff((1:n1).', R.pairs(:, 1));
    R.unmatched2 = setdiff((1:n2).', R.pairs(:, 2));

    % Unmatched systems on the low / high EFL side of their yield set
    R.unmatchedLow  = [sum(E1(R.unmatched1) < median(E1)), sum(E2(R.unmatched2) < median(E2))];
    R.unmatchedHigh = [numel(R.unmatched1), numel(R.unmatched2)] - R.unmatchedLow;

    % Products if the mean magnification were at the target
    passC = abs(Mcentered - magTarget) <= magTol(t) + 1e-12;
    R.nMatchedCentered = size(bestMatching(Mcentered, passC, E1, E2, magTarget), 1);

    res(t) = R; %#ok<SAGROW>
end

% Partner of every file for each tolerance ("" = unmatched or not in the yield set)
for s = 1:2
    sys(s).partner = strings(numel(sys(s).fileNames), nTol);
end
for t = 1:nTol
    sys(1).partner(idx1(res(t).pairs(:, 1)), t) = names2(res(t).pairs(:, 2));
    sys(2).partner(idx2(res(t).pairs(:, 2)), t) = names1(res(t).pairs(:, 1));
end

%% Command window summary
fprintf('\n===== MATCHED ASSEMBLY =====\n');
fprintf('M = %.6f x |EFL_%s| / |EFL_%s|   nominal M = %.5f, mean M of all pairs = %.5f\n', ...
    magConst, sys(1).name, sys(2).name, mNominal, mMean);
fprintf('Yield sets: %d %s, %d %s  ->  at most %d products (MC files: %d)\n', ...
    n1, sys(1).name, n2, sys(2).name, nMax, nMcRef);
for t = 1:nTol
    R = res(t);
    fprintf('\n--- M = %.4f +/- %g  [%.4f, %.4f] ---\n', magTarget, R.tol, R.low, R.high);
    fprintf('Matched assembly : %d products = %.1f %% of %d   (overall %.1f %% of %d MC)\n', ...
        R.nMatched, R.matchPct, nMax, R.overallPct, nMcRef);
    fprintf('Random assembly  : %.1f products expected (%.2f %% of all pairs compatible)\n', ...
        R.randomProducts, R.compatPct);
    fprintf('Mean centered    : %d products\n', R.nMatchedCentered);
    fprintf('Unmatched        : %d %s (low/high EFL %d/%d), %d %s (low/high EFL %d/%d)\n', ...
        numel(R.unmatched1), sys(1).name, R.unmatchedLow(1), R.unmatchedHigh(1), ...
        numel(R.unmatched2), sys(2).name, R.unmatchedLow(2), R.unmatchedHigh(2));
    fprintf('Method           : %s\n', R.method);
end

%% Output folder: "EŞLEŞTİRME ANALİZİ - %90 YIELD" in the system 1 folder
outName = ['E' char(350) 'LE' char(350) 'T' char(304) 'RME ANAL' char(304) 'Z' char(304) ...
    sprintf(' - %%%g YIELD', yieldPct)];
outFolder = fullfile(sys(1).folder, outName);

if ~isfolder(outFolder)
    mkdir(outFolder);
end

colFail = [0.85 0.33 0.10];
colBlue = [0 0.45 0.74];
colGray = [0.75 0.75 0.75];
tolPalette = [0.10 0.50 0.10; 0.55 0.75 0.20; 0.93 0.69 0.13; 0.49 0.18 0.56; 0.30 0.75 0.93];
tolColors = tolPalette(mod(0:nTol - 1, size(tolPalette, 1)) + 1, :);
tolMarkers = {'o', 's', '^', 'd', 'v'};
eflLabel = @(s) sprintf('|EFL| %s (mm)', sys(s).name);
tolLabels = arrayfun(@(x) sprintf('\\pm%g', x), magTol, 'UniformOutput', false);

%% Figure 1: number of products
fig1 = figure('Name', 'Products', 'Color', 'w', 'Position', [100 100 950 560]);
barData = [[res.randomProducts].', [res.nMatched].', [res.nMatchedCentered].'];
if nTol == 1
    hb = bar([1 2], [barData; nan nan nan]);   % a single group needs a dummy row
    xlim([0.5 1.5]);
else
    hb = bar(1:nTol, barData);
end
hb(1).FaceColor = colGray;
hb(2).FaceColor = colBlue;
hb(3).FaceColor = [0.47 0.67 0.19];
hb(1).DisplayName = 'Random assembly (expected)';
hb(2).DisplayName = 'Matched assembly (best)';
hb(3).DisplayName = 'Matched, mean M centered';
hold on;
yline(nMax, '--k', sprintf('Max. %d products', nMax), 'LineWidth', 1.2, ...
    'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
try   % value labels on the bars (R2019b+)
    for b = 1:3
        for t = 1:nTol
            text(hb(b).XEndPoints(t), hb(b).YEndPoints(t), sprintf('%.0f', barData(t, b)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
        end
    end
catch
end
hold off;
set(gca, 'XTick', 1:nTol, 'XTickLabel', tolLabels);
ylim([0 nMax * 1.18]);
xlabel('Magnification tolerance');
ylabel('Number of products');
title(sprintf('Products from %d %s + %d %s (%g%% MTF yield sets)', n1, sys(1).name, n2, sys(2).name, yieldPct));
legend(hb, 'Location', 'southeast');
grid on;
saveFigure(fig1, fullfile(outFolder, 'Products.png'));

%% Figure 2: matched pairs in the EFL plane, one panel per tolerance
fig2 = figure('Name', 'Matched Pairs', 'Color', 'w', 'Position', [120 60 900 300 + 300 * nTol]);
e1Line = linspace(min(E1), max(E1), 50);
for t = 1:nTol
    R = res(t);
    subplot(nTol, 1, t);
    hold on;
    plot(e1Line, magConst * e1Line / magTarget, '-k', 'LineWidth', 1, 'DisplayName', sprintf('M = %.2f', magTarget));
    plot(e1Line, magConst * e1Line / R.low, '--', 'Color', tolColors(t, :), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('M = %.2f \\pm %g', magTarget, R.tol));
    plot(e1Line, magConst * e1Line / R.high, '--', 'Color', tolColors(t, :), 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
    if R.nMatched > 0
        plot(E1(R.pairs(:, 1)), E2(R.pairs(:, 2)), tolMarkers{mod(t - 1, numel(tolMarkers)) + 1}, ...
            'Color', tolColors(t, :), 'MarkerFaceColor', tolColors(t, :), 'MarkerSize', 4, ...
            'DisplayName', sprintf('Matched pairs (%d)', R.nMatched));
    end
    plot(sys(1).nominalEFL, sys(2).nominalEFL, 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y', 'DisplayName', 'Nominal');
    hold off;
    xlabel(eflLabel(1));
    ylabel(eflLabel(2));
    title(sprintf('\\pm%g: %d of %d possible products', R.tol, R.nMatched, nMax));
    legend('Location', 'best');
    grid on;
    box on;
end
saveFigure(fig2, fullfile(outFolder, 'Matched_Pairs_Map.png'));

%% Figure 3: EFL of the unmatched systems
fig3 = figure('Name', 'Unmatched Systems', 'Color', 'w', 'Position', [140 140 950 700]);
eflSets = {E1, E2(:)};
for s = 1:2
    subplot(2, 1, s);
    [~, binEdges] = histcounts(eflSets{s}, nHistBins);
    histogram(eflSets{s}, binEdges, 'FaceColor', colGray, ...
        'DisplayName', sprintf('%g%% yield set (%d)', yieldPct, numel(eflSets{s})));
    hold on;
    for t = nTol:-1:1
        if s == 1
            un = E1(res(t).unmatched1);
        else
            un = E2(res(t).unmatched2);
        end
        histogram(un, binEdges, 'FaceColor', tolColors(t, :), 'FaceAlpha', 0.85, ...
            'DisplayName', sprintf('Unmatched at \\pm%g (%d)', res(t).tol, numel(un)));
    end
    xline(sys(s).nominalEFL, ':k', 'LineWidth', 1.5, 'DisplayName', sprintf('Nominal %.4f mm', sys(s).nominalEFL));
    xline(median(eflSets{s}), '--k', 'LineWidth', 1, 'DisplayName', 'Median of the yield set');
    hold off;
    xlabel(eflLabel(s));
    ylabel('Number of systems');
    title(sprintf('%s: systems left without a partner', sys(s).name));
    legend('Location', 'best');
    grid on;
end
saveFigure(fig3, fullfile(outFolder, 'Unmatched_EFL.png'));

%% Figure 4: magnification of the matched pairs vs. all possible pairs
fig4 = figure('Name', 'Matched Magnification', 'Color', 'w', 'Position', [160 160 950 540]);
allM = M(:);
mEdges = linspace(min([allM; magTarget - magTol(end)]), max([allM; magTarget + magTol(end)]), 2 * nHistBins + 1);
histogram(allM, mEdges, 'Normalization', 'probability', 'FaceColor', colGray, 'EdgeColor', 'none', ...
    'DisplayName', sprintf('All possible pairs (random assembly, %d)', numel(allM)));
hold on;
for t = 1:nTol
    if res(t).nMatched > 0
        histogram(res(t).matchM, mEdges, 'Normalization', 'probability', 'DisplayStyle', 'stairs', ...
            'EdgeColor', tolColors(t, :), 'LineWidth', 1.8, ...
            'DisplayName', sprintf('Matched pairs \\pm%g (%d)', res(t).tol, res(t).nMatched));
    end
    xline(res(t).low, '--', 'Color', tolColors(t, :), 'LineWidth', 1.1, 'HandleVisibility', 'off');
    xline(res(t).high, '--', 'Color', tolColors(t, :), 'LineWidth', 1.1, 'HandleVisibility', 'off');
end
xline(magTarget, '-k', 'HandleVisibility', 'off');
xline(mMean, '-.', 'Color', colFail, 'LineWidth', 1.3, 'DisplayName', sprintf('Mean M of all pairs = %.4f', mMean));
hold off;
xlabel('Magnification M');
ylabel('Fraction of pairs');
title('Magnification: matched pairs vs. random assembly');
legend('Location', 'best');
grid on;
saveFigure(fig4, fullfile(outFolder, 'Matched_Magnification.png'));

%% Export Excel
outputFile = fullfile(outFolder, 'Matching_Results.xlsx');

if isfile(outputFile)
    delete(outputFile);
end

tolTag = @(x) strrep(sprintf('%g', x), '.', 'p');   % 0.01 -> "0p01" for sheet / column names

summaryNames = [ ...
    "Tolerance (+/-)"; "Lower limit M"; "Upper limit M"; ...
    sys(1).name + " systems in yield set"; sys(2).name + " systems in yield set"; ...
    "Maximum possible products"; "MC files per system"; ...
    "Matched assembly: products"; "Matched assembly: products / maximum"; ...
    "Overall: products / MC files"; "Random assembly: expected products"; ...
    "Compatible share of all possible pairs"; "Matched, mean M centered: products"; ...
    "Unmatched " + sys(1).name; "Unmatched " + sys(1).name + " (EFL below median)"; ...
    "Unmatched " + sys(1).name + " (EFL above median)"; ...
    "Unmatched " + sys(2).name; "Unmatched " + sys(2).name + " (EFL below median)"; ...
    "Unmatched " + sys(2).name + " (EFL above median)"; ...
    "Matched pairs: min M"; "Matched pairs: max M"; "Matched pairs: mean |M - target|"; ...
    "Nominal M"; "Mean M of all pairs"];
summaryUnits = [""; ""; ""; "count"; "count"; "count"; "count"; "count"; "%"; "%"; "count"; "%"; "count"; ...
    "count"; "count"; "count"; "count"; "count"; "count"; ""; ""; ""; ""; ""];
summaryTable = table(summaryNames, summaryUnits, 'VariableNames', {'Parameter', 'Unit'});
for t = 1:nTol
    R = res(t);
    summaryTable.(['Tol_' tolTag(R.tol)]) = [R.tol; R.low; R.high; n1; n2; nMax; nMcRef; ...
        R.nMatched; R.matchPct; R.overallPct; R.randomProducts; R.compatPct; R.nMatchedCentered; ...
        numel(R.unmatched1); R.unmatchedLow(1); R.unmatchedHigh(1); ...
        numel(R.unmatched2); R.unmatchedLow(2); R.unmatchedHigh(2); ...
        statOrNaN(@min, R.matchM); statOrNaN(@max, R.matchM); statOrNaN(@mean, abs(R.matchM - magTarget)); ...
        mNominal; mMean];
end
writetable(summaryTable, outputFile, 'Sheet', 'Summary');

for t = 1:nTol
    R = res(t);
    p1 = R.pairs(:, 1);
    p2 = R.pairs(:, 2);
    matchTable = table((1:R.nMatched).', names1(p1), names2(p2), E1(p1), E2(p2).', R.matchM, R.matchM - magTarget, ...
        'VariableNames', {'Product', char("File_" + sys(1).name), char("File_" + sys(2).name), ...
        char("EFL_" + sys(1).name + "_mm"), char("EFL_" + sys(2).name + "_mm"), 'Magnification', 'Deviation'});
    writetable(matchTable, outputFile, 'Sheet', ['Matched tol ' tolTag(R.tol)]);

    unSystem = [repmat(sys(1).name, numel(R.unmatched1), 1); repmat(sys(2).name, numel(R.unmatched2), 1)];
    unFile = [names1(R.unmatched1); names2(R.unmatched2)];
    unEfl = [E1(R.unmatched1); E2(R.unmatched2).'];
    writetable(table(unSystem, unFile, unEfl, 'VariableNames', {'System', 'FileName', 'EFL_abs_mm'}), ...
        outputFile, 'Sheet', ['Unmatched tol ' tolTag(R.tol)]);
end

for s = 1:2
    fileTable = table(sys(s).fileNames, sys(s).EFL, sys(s).mtfMetric, sys(s).mtfRank, sys(s).inYield, ...
        'VariableNames', {'FileName', 'EFL_abs_mm', 'MTF_criterion', 'MTF_rank', 'InYieldSet'});
    for t = 1:nTol
        fileTable.(sprintf('Partner_%s_tol_%s', sys(3 - s).name, tolTag(magTol(t)))) = sys(s).partner(:, t);
    end
    writetable(fileTable, outputFile, 'Sheet', char(sys(s).name + " Files"));
end

settingNames = [ ...
    "System 1"; "System 1 folder"; "System 1 nominal file"; "System 1 nominal |EFL| (mm)"; "System 1 afocal"; ...
    "System 1 pixel pitch (um)"; "System 1 MTF frequency (cyc/mm)"; "System 1 MTF threshold"; "System 1 MTF fields"; ...
    "System 2"; "System 2 folder"; "System 2 nominal file"; "System 2 nominal |EFL| (mm)"; "System 2 afocal"; ...
    "System 2 pixel pitch (um)"; "System 2 MTF frequency (cyc/mm)"; "System 2 MTF threshold"; "System 2 MTF fields"; ...
    "File prefix"; "Yield level (%)"; "MTF frequency (fraction of Nyquist)"; ...
    "Pixel scale k (OLED px per detector px)"; "Magnification constant k*p2/p1"; ...
    "Magnification formula"; "Target magnification"; "Tolerances (+/-)"; "Matching method"];
settingValues = strings(numel(settingNames), 1);
for s = 1:2
    o = (s - 1) * 9;
    settingValues(o + (1:9)) = [sys(s).name; string(sys(s).folder); string(sys(s).nominalName); ...
        string(sys(s).nominalEFL); string(sys(s).isAfocal); string(sys(s).pitch); string(sys(s).mtfFreq); ...
        string(sys(s).threshold); fieldText(mtfFields{s})];
end
settingValues(19:end) = [string(filePrefix); string(yieldPct); string(mtfFreqFrac); ...
    string(pixelScale); string(magConst); "M = k*p2/p1 * |EFL_1| / |EFL_2|"; ...
    string(magTarget); strjoin(string(magTol), ", "); res(1).method];
writetable(table(settingNames, settingValues, 'VariableNames', {'Parameter', 'Value'}), ...
    outputFile, 'Sheet', 'Settings');

%% Save raw data
save(fullfile(outFolder, 'Matching_Data.mat'), 'sys', 'E1', 'E2', 'names1', 'names2', 'M', 'res', ...
    'magConst', 'magTarget', 'magTol', 'yieldPct', 'mtfFreqFrac', 'pixelPitch_um', 'nPixels', 'pixelScale');

%% PDF report (A4)
pdfFile = fullfile(outFolder, 'Matching_Report.pdf');

oldPages = dir(fullfile(outFolder, 'Matching_Report_page*.pdf'));
for i = 1:numel(oldPages)
    delete(fullfile(outFolder, oldPages(i).name));
end
if isfile(pdfFile)
    delete(pdfFile);
end
if isfile(pdfFile)   % still open in a PDF viewer
    pdfFile = fullfile(outFolder, ['Matching_Report_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.pdf']);
    fprintf('\nMatching_Report.pdf is open in another program; new report: %s\n', pdfFile);
end

n1Name = char(sys(1).name);
n2Name = char(sys(2).name);
tolListText = strjoin(arrayfun(@(x) sprintf('±%g', x), magTol, 'UniformOutput', false), ' ve ');
afocalWord = ["Hayır", "Evet"];

% --- Page: summary -------------------------------------------------------
summaryLines = {
    sprintf('Tarih:                 %s', char(datetime('now', 'Format', 'dd.MM.yyyy HH:mm')))
    sprintf('Büyütme isteri:        M = %.4f, tolerans %s', magTarget, tolListText)
    sprintf('Formül:                M = (%g × %.2f µm / %.2f µm) × |EFL_%s| / |EFL_%s| = %.6f × oran', ...
        pixelScale, pixelPitch_um(2), pixelPitch_um(1), n1Name, n2Name, magConst)
    sprintf('Yield seviyesi:        %%%g  (MTF, Nyquist × %g)', yieldPct, mtfFreqFrac)
    ''
    };
for s = 1:2
    nYield = sum(sys(s).inYield);
    summaryLines = [summaryLines; {
        sprintf('SİSTEM %d – %s', s, sys(s).name)
        sprintf('  Klasör:              %s', sys(s).folder)
        sprintf('  Nominal dosya:       %s   |EFL| = %.4f mm', sys(s).nominalName, sys(s).nominalEFL)
        sprintf('  Afokal mod:          %s', afocalWord(sys(s).isAfocal + 1))
        sprintf('  MTF ölçütü:          %s, %.2f cyc/mm, eşik %.4f', fieldText(mtfFields{s}), sys(s).mtfFreq, sys(s).threshold)
        sprintf('  Dosyalar:            %d MC, %d geçerli, %d yield kümesinde', ...
            numel(sys(s).fileNames), sum(sys(s).valid), nYield)
        sprintf('  |EFL| yield kümesi:  ort. %.4f mm, std %.4f mm, %.4f – %.4f mm', ...
            mean(sys(s).EFL(sys(s).inYield)), std(sys(s).EFL(sys(s).inYield)), ...
            min(sys(s).EFL(sys(s).inYield)), max(sys(s).EFL(sys(s).inYield)))
        ''}]; %#ok<AGROW>
end
summaryLines = [summaryLines; {
    'GENEL'
    sprintf('  Yield kümeleri:      %d %s, %d %s  ->  en fazla %d ürün', n1, n1Name, n2, n2Name, nMax)
    sprintf('  Nominal M:           %.5f', mNominal)
    sprintf('  Tüm çiftlerin ort. M: %.5f', mMean)
    ''
    }];
for t = 1:nTol
    R = res(t);
    summaryLines = [summaryLines; {
        sprintf('TOLERANS ±%g  (M = %.4f – %.4f)', R.tol, R.low, R.high)
        sprintf('  EŞLEŞTİRMELİ MONTAJ: %d / %d ürün  ->  %%%.1f', R.nMatched, nMax, R.matchPct)
        sprintf('  Toplam verim:        %d / %d MC dosyası  ->  %%%.1f', R.nMatched, nMcRef, R.overallPct)
        sprintf('  Rastgele montaj:     %.1f ürün (tüm çiftlerin %%%.2f''si uyumlu)', R.randomProducts, R.compatPct)
        sprintf('  Ortalama merkezde:   %d ürün', R.nMatchedCentered)
        sprintf('  Açıkta kalan %s:   %d  (EFL küçük / büyük: %d / %d)', n1Name, numel(R.unmatched1), ...
            R.unmatchedLow(1), R.unmatchedHigh(1))
        sprintf('  Açıkta kalan %s:   %d  (EFL küçük / büyük: %d / %d)', n2Name, numel(R.unmatched2), ...
            R.unmatchedLow(2), R.unmatchedHigh(2))
        sprintf('  Eşleşen çiftlerin M: %.4f – %.4f', statOrNaN(@min, R.matchM), statOrNaN(@max, R.matchM))
        ''}]; %#ok<AGROW>
end
summaryFig = textPage('EŞLEŞTİRMELİ MONTAJ ANALİZ RAPORU', summaryLines, 'Courier New', 8);

% --- Page: results and assessment ----------------------------------------
assess = {
    'SORU'
    sprintf(['MTF''e göre %%%g yield kümelerinde %d %s ve %d %s var. EFL''ler ölçülüp her %s''a uygun bir %s ' ...
     'seçilerek montaj yapılırsa (her sistem en fazla bir kez), büyütme isterini sağlayan en fazla kaç ürün çıkar?'], ...
     yieldPct, n1, n1Name, n2, n2Name, n1Name, n2Name)
    ''
    'SONUÇ'
    };
for t = 1:nTol
    R = res(t);
    assess{end + 1, 1} = sprintf(['•  ±%g: en iyi eşleştirmeyle %d ürün çıkar (en fazla %d ürünün %%%.1f''i). ' ...
        'MC dosya sayısına göre toplam verim %d / %d = %%%.1f. Eşleştirme yapılmadan rastgele montajda ' ...
        'beklenen ürün sayısı %.1f; eşleştirmenin kazancı %.1f ürün.'], R.tol, R.nMatched, nMax, R.matchPct, ...
        R.nMatched, nMcRef, R.overallPct, R.randomProducts, R.nMatched - R.randomProducts); %#ok<SAGROW>
end

assess = [assess; {''; 'AÇIKTA KALAN SİSTEMLER'}];
for t = 1:nTol
    R = res(t);
    assess{end + 1, 1} = sprintf(['•  ±%g: %d %s ve %d %s eşsiz kalıyor. %s''lardan %d''i EFL''si küçük, %d''i ' ...
        'büyük tarafta; %s''lerden %d''i küçük, %d''i büyük tarafta (yield kümesinin medyanına göre).'], ...
        R.tol, numel(R.unmatched1), n1Name, numel(R.unmatched2), n2Name, n1Name, R.unmatchedLow(1), ...
        R.unmatchedHigh(1), n2Name, R.unmatchedLow(2), R.unmatchedHigh(2)); %#ok<SAGROW>
end
assess{end + 1, 1} = ['Açıkta kalanlar bir tarafta toplanıyorsa sebep ortalama büyütmenin hedeften kaymasıdır; ' ...
    'iki tarafa dağılıyorsa sebep EFL saçılımıdır (dağılımın uçlarındaki sistemlere uygun eş yoktur). ' ...
    'Dosya listesi Excel''de "Unmatched" sayfalarındadır.'];

assess = [assess; {''; 'NOMİNAL ORAN DÜZELTİLİRSE'}];
assess{end + 1, 1} = sprintf(['Nominal büyütme %.5f, tüm olası çiftlerin ortalama büyütmesi %.5f (hedeften sapma %+.5f). ' ...
    'Ortalama hedefe getirilseydi (bir sistemin nominal EFL''si düzeltilerek):'], mNominal, mMean, mMean - magTarget);
for t = 1:nTol
    R = res(t);
    assess{end + 1, 1} = sprintf('•  ±%g: %d yerine %d ürün (%+d).', R.tol, R.nMatched, R.nMatchedCentered, ...
        R.nMatchedCentered - R.nMatched); %#ok<SAGROW>
end
assess{end + 1, 1} = ['Kazanç büyükse eşleştirmeden önce nominal EFL oranını düzeltmek en etkili adımdır; ' ...
    'küçükse ürün sayısını EFL saçılımı sınırlıyordur.'];

assess = [assess; {''; 'VARSAYIMLAR'}];
assess{end + 1, 1} = ['•  EFL''ler hatasız ölçülüyor ve yield kümesinin tamamı aynı anda eşleştiriliyor; sonuç bu ' ...
    'koşullarda ulaşılabilecek en fazla ürün sayısıdır. Ölçüm belirsizliği veya küçük partiler (örneğin 10''ar ' ...
    'sistemlik) ürün sayısını düşürür.'];
assess{end + 1, 1} = ['•  Sonuç tek bir Monte Carlo örneğine aittir; MC dosya sayısı azsa ürün sayısı birkaç adet ' ...
    'oynayabilir.'];
assessFig = textPage('1. SONUÇ VE DEĞERLENDİRME', assess, 'Arial', 10);

% --- Page: method --------------------------------------------------------
method = {
    'ADIM 1 – DOSYALAR, EFL VE MTF'
    sprintf(['Her sistem için ayrı bir Monte Carlo klasörü seçilir. Adı %s ile başlayan .zmx dosyaları MC dosyalarıdır; ' ...
     'adı MC_ ile başlamayan dosya nominal tasarımdır. Her dosyanın EFL''si EFFL operandıyla okunur ve her zaman ' ...
     'mutlak değeri kullanılır. MTF polikromatik FFT MTF''dir; ölçüt her sistemin kendi Nyquist frekansının %g ' ...
     'katında (%s: %.2f cyc/mm, %s: %.2f cyc/mm), seçilen fieldlerde T ve S MTF''lerinin en düşüğüdür.'], ...
     filePrefix, mtfFreqFrac, n1Name, sys(1).mtfFreq, n2Name, sys(2).mtfFreq)
    ''
    'ADIM 2 – YIELD KÜMELERİ'
    sprintf(['Her sistemde dosyalar MTF ölçütüne göre sıralanır ve en iyi %%%g seçilir: %d %s ve %d %s. ' ...
     'MTF şartını sağlamayan sistemler üretimde elenecekleri için eşleştirmeye girmez. Bu yüzden en fazla ' ...
     'min(%d, %d) = %d ürün çıkabilir.'], yieldPct, n1, n1Name, n2, n2Name, n1, n2, nMax)
    ''
    'ADIM 3 – UYUMLU ÇİFTLER'
    sprintf(['Her olası (%s, %s) çifti için büyütme M = %.6f × |EFL_1| / |EFL_2| hesaplanır. |M − %.2f| ≤ tolerans ' ...
     'olan çiftler uyumludur.'], n1Name, n2Name, magConst, magTarget)
    ''
    'ADIM 4 – EN İYİ EŞLEŞTİRME'
    ['Amaç, her sistemi en fazla bir kez kullanarak en çok sayıda uyumlu çift kurmaktır (maksimum eşleştirme). ' ...
     'Bir ' n1Name '''a uyan ' n2Name ' EFL''leri sürekli bir aralıktır ve bu aralık ' n1Name ' EFL''si büyüdükçe ' ...
     'yukarı kayar. Bu yapıda şu kural kesin en büyük ürün sayısını verir:']
    ['     ' n1Name '''lar EFL''ye göre küçükten büyüğe alınır; her birine uyumlu ve henüz kullanılmamış ' ...
     'en küçük EFL''li ' n2Name ' verilir.']
    ['Neden en iyisi: küçük EFL''li ' n1Name '''ın uyabileceği ' n2Name '''lar, daha büyük EFL''li ' n1Name ...
     '''lara da uyuyorsa, ona en küçüğünü vermek sonrakilerin seçeneklerini hiç azaltmaz. Aynı ürün sayısını ' ...
     'veren eşleştirmeler arasından toplam |M − hedef| en küçük olan seçilir (MATLAB matchpairs).']
    ''
    'ADIM 5 – KARŞILAŞTIRMALAR'
    sprintf(['•  Rastgele montaj: eşleştirme yapılmadan birebir rastgele dağıtımda beklenen ürün sayısı = %d × ' ...
     '(uyumlu çiftlerin oranı).'], nMax)
    ['•  Ortalama merkezde: tüm büyütmeler ortalama tam hedefe gelecek şekilde ölçeklenip (nominal EFL oranı ' ...
     'düzeltilmiş gibi) eşleştirme tekrarlanır.']
    '•  Toplam verim: ürün sayısı / MC dosya sayısı; MTF elemesi ve eşleştirme kaybı birlikte.'
    };
methodFig = textPage('2. YÖNTEM', method, 'Arial', 10);

% --- Figure pages --------------------------------------------------------
plotFigs = gobjects(4, 1);
plotFigs(1) = figurePage(fig1, '3. ÜRÜN SAYISI', {
    'NASIL OKUNUR'
    sprintf(['•  Her tolerans için üç çubuk: gri = eşleştirme yapılmadan rastgele montajda beklenen ürün, mavi = ' ...
     'en iyi eşleştirmeyle çıkan ürün, yeşil = ortalama büyütme hedefe getirilmiş olsaydı eşleştirmeyle çıkacak ' ...
     'ürün. Kesikli çizgi olası en fazla ürün sayısıdır (%d).'], nMax)
    'NASIL YORUMLANIR'
    '•  Mavi çubuk kesikli çizgiye ne kadar yakınsa eşleştirme o kadar başarılıdır.'
    '•  Gri ile mavi arasındaki fark, EFL ölçüp eşleştirmenin kazancıdır.'
    ['•  Mavi ile yeşil arasındaki fark, nominal EFL oranını düzeltmenin kazancıdır. Bu fark büyükse önce ' ...
     'nominal oran düzeltilmelidir.']
    }, false);
plotFigs(2) = figurePage(fig2, '4. EŞLEŞTİRİLEN ÇİFTLER', {
    'NASIL OKUNUR'
    sprintf(['•  Her panel bir toleranstır. Her nokta bir üründür: yatay eksen %s |EFL|, dikey eksen ona ' ...
     'eşleştirilen %s |EFL|. Düz çizgi M = %.2f, kesikli çizgiler tolerans sınırlarıdır; tüm noktalar bu ' ...
     'şeridin içindedir. Sarı yıldız nominal tasarımdır.'], n1Name, n2Name, magTarget)
    'NASIL YORUMLANIR'
    ['•  Noktalar şeridin bir kenarına yığılıyorsa eşleştirme zorlanıyor demektir; büyük ihtimalle ortalama ' ...
     'büyütme hedeften kaymıştır.']
    ['•  Şeridin uçlarında nokta yoksa, EFL''si uçta olan sistemlere uygun eş bulunamamıştır (açıkta kalanlar ' ...
     'sonraki sayfada).']
    }, false);
plotFigs(3) = figurePage(fig3, '5. AÇIKTA KALAN SİSTEMLER', {
    'NASIL OKUNUR'
    ['•  Gri histogram yield kümesindeki tüm sistemlerin |EFL''si, renkli histogramlar her toleransta eşsiz ' ...
     'kalan sistemlerdir. Kesikli çizgi yield kümesinin medyanı, noktalı çizgi nominal |EFL''dir.']
    'NASIL YORUMLANIR'
    ['•  Açıkta kalanlar tek bir uçta toplanıyorsa sebep ortalama kaymadır; nominal oran düzeltilerek o ' ...
     'sistemler de eşleştirilebilir.']
    ['•  İki uçta da açıkta kalan varsa sebep saçılımdır; o sistemin EFL toleransı sıkılaştırılmalı veya ' ...
     'diğer sistemden daha fazla stok tutulmalıdır.']
    }, false);
plotFigs(4) = figurePage(fig4, '6. EŞLEŞTİRİLEN ÇİFTLERİN BÜYÜTMESİ', {
    'NASIL OKUNUR'
    ['•  Gri histogram tüm olası çiftlerin büyütme dağılımıdır, yani eşleştirme yapılmadan rastgele montajda ' ...
     'görülecek dağılım. Renkli basamak çizgileri eşleştirilen çiftlerin büyütmesidir. İkisi de oran olarak ' ...
     'çizilmiştir (toplam = 1). Kesikli çizgiler tolerans sınırlarıdır.']
    'NASIL YORUMLANIR'
    ['•  Eşleştirme, büyütmeyi tolerans bandının içine toplar. Eşleştirilen çiftler hedef çizgisi çevresinde ' ...
     'yoğunsa ürünlerin büyütmesi hedefe yakındır; bir sınıra yakın toplanıyorsa pay azdır.']
    }, false);

% --- Assemble pages ------------------------------------------------------
allFigs = [summaryFig; assessFig; methodFig; plotFigs];
nPages = numel(allFigs);
headerText = sprintf('Eşleştirmeli Montaj Analizi  |  %s: %s  |  %s: %s', ...
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
    fprintf('\nMulti-page PDF export failed: %s\n', pdfErr.message);
    if isfile(pdfFile)
        delete(pdfFile);
    end
    for p = 1:nPages
        print(allFigs(p), fullfile(outFolder, sprintf('Matching_Report_page%d.pdf', p)), '-dpdf');
    end
    fprintf('Report saved as separate pages (Matching_Report_page*.pdf).\n');
end
close(allFigs);

%% Text summary
summaryFile = fullfile(outFolder, 'Analysis_Summary.txt');
fid = fopen(summaryFile, 'w', 'n', 'UTF-8');
if fid > 0
    fprintf(fid, '%s\r\n', 'EŞLEŞTİRMELİ MONTAJ ANALİZ ÖZETİ', '');
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
function [pairs, method] = bestMatching(M, pass, E1, E2, magTarget)
% Largest set of disjoint (system 1, system 2) pairs that meet the requirement.
% pass(i, j): pair i-j meets |M - target| <= tol. For each system 1 file the
% acceptable system 2 EFLs form one interval that moves up with EFL_1, so the
% greedy rule below is exact: take system 1 files in increasing EFL and give
% each the unused acceptable system 2 file with the smallest EFL.
[~, rowOrder] = sort(E1(:));
[~, colOrder] = sort(E2(:));
used = false(numel(E2), 1);
pairs = zeros(0, 2);
for i = rowOrder.'
    ok = pass(i, colOrder).' & ~used(colOrder);
    cand = colOrder(ok);
    if ~isempty(cand)
        used(cand(1)) = true;
        pairs(end + 1, :) = [i, cand(1)]; %#ok<AGROW>
    end
end
method = "greedy (maximum count)";

% Among the maximum matchings, take the one with the smallest total |M - target|
% (matchpairs, R2019a+). An unmatched file costs 1 >> any |M - target|, so the
% number of pairs is maximised first; incompatible pairs get a prohibitive cost.
try
    cost = abs(M - magTarget);
    cost(~pass) = 1e6;
    mp = matchpairs(cost, 1);
    if size(mp, 1) == size(pairs, 1) && all(pass(sub2ind(size(M), mp(:, 1), mp(:, 2))))
        pairs = mp;
        method = "maximum matching, minimum total |M - target| (matchpairs)";
    end
catch
end

% Sort the products by EFL of system 1
[~, order] = sort(E1(pairs(:, 1)));
pairs = pairs(order, :);
end

function v = statOrNaN(fun, x)
% fun(x), or NaN for an empty x
if isempty(x)
    v = NaN;
else
    v = fun(x);
end
end

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
