%% DİSTORSİYONLU BÜYÜTME ANALİZİ (ZOS-API, MATLAB)
% NE YAPAR
%   Birlikte çalışan iki sistemin (1: MWIR objektif + dedektör, 2: OLED +
%   kolimatör) toplam büyütmesini, distorsiyonu da dahil ederek, görüş
%   alanının birkaç noktasında (normFields, örn. 0, 0.5, 0.7, 1 normalize alan)
%   hesaplar. Her sistemde MTF'e göre %yieldPct yield kümesine giren dosyalar
%   seçilir; 1. sistemin yield kümesindeki her dosya 2. sistemin yield
%   kümesindeki her dosyayla eşleştirilir (n1 x n2 çift). Her çiftin on-axis
%   (0) ve off-axis alanlardaki büyütmesi hesaplanır ve alan bazında büyütme
%   değerleri / dağılımları verilir. Geçti / kaldı değerlendirmesi yapılmaz.
%   On-axis (0) paraksiyel EFL ile hesaplanır; distorsiyonsuz önceki
%   analizle (Iki_Sistem_Buyutme_Analizi.m) aynı büyütmeyi verir.
%   Başka bir scriptin çıktısı kullanılmaz; her şey bu scriptte OpticStudio
%   ile hesaplanır.
%
% NASIL ÇALIŞIR
%   1. Kullanıcı önce 1. sistemin (MWIR), sonra 2. sistemin (OLED) Monte Carlo
%      klasörünü seçer. Adı filePrefix ("MC_T") ile başlayan .zmx dosyaları MC
%      dosyalarıdır; adı "MC_" ile başlamayan .zmx dosyası nominal tasarımdır.
%   2. Her dosya için OpticStudio ile:
%        - Paraksiyel |EFL| (EFFL operandı, HER ZAMAN mutlak değer)
%        - Her normalize alanda (Hy) gerçek ana ışının (Px = Py = 0, birincil
%          dalga boyu) yüksekliği ve açısı, buradan etkin EFL:
%              f_eff(h) = |y(h)| / |tan θ(h)|      (h = 0 için paraksiyel EFL)
%          Dosya afokal modda değilse (nesne sonsuzda, görüntü düzlemde):
%              y = görüntü yüzeyinde REAY, tan θ = nesne uzayında REAB / REAC
%              (yüzey 0'daki doğrultu kosinüsleri)
%          Dosya afokal moddaysa (nesne = OLED, çıkış açısal):
%              y = nesne yüzeyinde REAY, tan θ = görüntü uzayında REAB / REAC
%        - Distorsiyon: D(h) = f_eff(h) / EFL_paraksiyel − 1  (F-Tan(θ)
%          distorsiyonu; Zemax'in Field Curvature / Distortion grafiğiyle
%          karşılaştırılabilir)
%        - Polikromatik FFT MTF (yield seçimi için, büyütme scriptiyle aynı)
%   3. Her sistemde MTF'e göre en iyi %yieldPct dosya seçilir.
%   4. Her çift ve her alan için büyütme hesaplanır.
%
% HESAPLAR
%   c = k × p_2 / p_1  (k = 1024 / 512 = 2 OLED pikseli / dedektör pikseli)
%   MWIR alanı h: dedektörde yükseklik y_1(h), giriş açısı tan θ_1(h)
%   OLED'deki yükseklik:            y_OLED = c × y_1(h)
%   Kolimatörün o yükseklikteki etkin EFL'si f_2(y_OLED): 2. sistemin alan
%     noktalarından  f_2(y) = EFL_2 + a·y² + b·y⁴  (dosya başına en küçük
%     kareler uyumu, y = 0'da tam olarak paraksiyel EFL) ile bulunur.
%   Çıkış açısı:                    tan θ_2 = y_OLED / f_2(y_OLED)
%   Büyütme:                        M(h) = tan θ_2 / tan θ_1(h)
%                                        = c × f_1,eff(h) / f_2(y_OLED)
%   Yaklaşık olarak M(h) ≈ M(0) × (1 + D_1(h) − D_2), yani kenar alanda
%   büyütme iki sistemin distorsiyon farkı kadar değişir.
%   Off-axis büyütme değişimi:      ΔM(h) = M(h) / M(0) − 1  (çift başına, %)
%   Her alan için: nominal çiftin M'si; tüm çiftlerin M ortalaması, std,
%   min, max ve yüzdelikleri; ΔM'nin aynı istatistikleri.
%
% AYARLAR (User settings)
%   sysName, pixelPitch_um, nPixels, pixelScale, filePrefix, yieldPct,
%   mtfFreqFrac, mtfFields, normFields, magTarget (yalnız referans çizgisi),
%   writePairsSheet, nHistBins
%
% ÇIKTILAR (1. sistemin klasöründe "BÜYÜTME-DİSTORSİYON ANALİZİ - %90 YIELD")
%   Distortion_by_Field.png      : iki sistemin alana göre distorsiyonu
%                                  (nominal, yield kümesi ortalaması ve aralığı)
%   Magnification_by_Field.png   : alana göre büyütme (ortalama, %1-%99, min-max)
%   Magnification_Change.png     : on-axis'e göre büyütme değişimi ΔM(h)
%   Magnification_Histograms.png : her alanın büyütme dağılımı
%   *.fig                        : grafiklerin düzenlenebilir MATLAB dosyaları
%   Distortion_Magnification_Results.xlsx : "Summary", "Field Stats",
%                                  "Percentiles", "Distortion", sistem başına
%                                  dosya tabloları, "Pairs" (her çiftin her
%                                  alandaki M'si), "Settings"
%   Distortion_Magnification_Report.pdf : A4 rapor; özet, değerlendirme, yöntem,
%                                  grafikler ve nasıl yorumlanacakları
%   Distortion_Magnification_Data.mat : tüm ham veriler
%   Analysis_Summary.txt, Analysis_Log.txt

clc;
clear;
close all;

% Command Window log; copied to the output folder at the end
diary off;
logFile = fullfile(tempdir, 'Distortion_Magnification_Log.txt');
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
normFields      = [0 0.5 0.7 1];     % normalized fields (Hy) of system 1 where M is evaluated
magTarget       = 1.00;              % target magnification (drawn as a reference line only)
writePairsSheet = true;              % write every pair to the "Pairs" sheet (Excel limit 1048575 rows)
nHistBins       = 60;                % histogram bins
% FFT MTF sample size is set in computeFileMtf (S_64x64 by default)

magConst = pixelScale * pixelPitch_um(2) / pixelPitch_um(1);   % c: OLED height = c * detector height
normFields = unique([0, abs(normFields(:)).']);   % the center (0) is always included
if any(normFields > 1)
    error('normFields must be between 0 and 1.');
end
nF = numel(normFields);
fieldNames = compose("%.2f", normFields);

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

%% Compute EFL, field data and MTF of both systems with OpticStudio
TheApplication = InitConnection();

if isempty(TheApplication)
    error('Could not connect to OpticStudio.');
end

tStart = tic;

try
    TheSystem = TheApplication.PrimarySystem;

    for s = 1:2
        % System 1 is evaluated at normFields; system 2 at the same normalized
        % fields (used for its distortion curve and the f_2(y) fit)
        [sys(s).EFL, sys(s).fEff, sys(s).yH, sys(s).tH, sys(s).mtfMetric, sys(s).nominal, ...
            sys(s).isAfocal, sys(s).mtfFreq] = computeSystemField(TheSystem, sys(s), normFields, ...
            mtfFreqFrac, mtfFields{s});
    end
catch ME
    CleanupConnection(TheApplication);
    rethrow(ME);
end

CleanupConnection(TheApplication);
fprintf('Computation time: %.1f s\n', toc(tStart));

%% Yield set of each system (files with a valid EFL, valid field data and valid MTF)
for s = 1:2
    fieldOk = all(isfinite(sys(s).fEff), 2) & all(sys(s).fEff > 0, 2);
    eflForYield = sys(s).EFL;
    eflForYield(~fieldOk) = NaN;
    [sys(s).valid, sys(s).inYield, sys(s).threshold, sys(s).mtfRank] = ...
        selectYield(eflForYield, sys(s).mtfMetric, yieldPct);

    if ~any(sys(s).inYield)
        error('%s: no file has a valid EFL, valid field data and a valid MTF.', sys(s).name);
    end

    sys(s).D = sys(s).fEff ./ sys(s).EFL - 1;                  % distortion per file and field
    sys(s).nominal.D = sys(s).nominal.fEff / sys(s).nominal.EFL - 1;

    fprintf('%s yield %g %%: MTF >= %.4f at %.2f cyc/mm -> %d of %d valid files\n', ...
        sys(s).name, yieldPct, sys(s).threshold, sys(s).mtfFreq, sum(sys(s).inYield), sum(sys(s).valid));
    fprintf('%s nominal distortion: %s\n', sys(s).name, ...
        strjoin(compose("h=%.2f: %+.3f %%", normFields(:), 100 * sys(s).nominal.D(:)).', ', '));
end


%% Field-dependent magnification of every pair
in1 = sys(1).inYield;
in2 = sys(2).inYield;
names1 = sys(1).fileNames(in1);
names2 = sys(2).fileNames(in2);
n1 = sum(in1);
n2 = sum(in2);
nPairs = n1 * n2;

E1 = sys(1).EFL(in1);                 % n1 x 1 paraxial |EFL|
E2 = sys(2).EFL(in2);                 % n2 x 1
F1 = sys(1).fEff(in1, :);             % n1 x nF effective EFL of system 1
Y1 = sys(1).yH(in1, :);               % n1 x nF detector height of system 1

% f_2(y) = EFL_2 + a*y^2 + b*y^4 for every system 2 file
coef2 = zeros(n2, 2);
fitErr2 = zeros(n2, 1);
Y2 = sys(2).yH(in2, :);
F2pts = sys(2).fEff(in2, :);
for j = 1:n2
    [coef2(j, :), fitErr2(j)] = fitFieldEfl(Y2(j, :), F2pts(j, :), E2(j));
end
[nomCoef2, nomFitErr2] = fitFieldEfl(sys(2).nominal.yH, sys(2).nominal.fEff, sys(2).nominal.EFL);

% OLED height of every system 1 field point, and f_2 there: n1 x n2 x nF
Yo = magConst * reshape(Y1, n1, 1, nF);
F2 = E2.' + coef2(:, 1).' .* Yo.^2 + coef2(:, 2).' .* Yo.^4;
M = magConst * reshape(F1, n1, 1, nF) ./ F2;
M(:, :, normFields == 0) = magConst * E1 ./ E2.';     % center: paraxial (same as the analysis without distortion)

% Nominal pair
yoNom = magConst * sys(1).nominal.yH;
f2Nom = sys(2).nominal.EFL + nomCoef2(1) * yoNom.^2 + nomCoef2(2) * yoNom.^4;
mNomField = magConst * sys(1).nominal.fEff ./ f2Nom;
mNomField(normFields == 0) = magConst * sys(1).nominal.EFL / sys(2).nominal.EFL;
dDistNom = mNomField / mNomField(1) - 1;              % magnification change caused by the nominal distortion

% Extrapolation check: OLED heights beyond the largest sampled height of system 2
extrapPct = 100 * (max(Yo(:)) / min(sys(2).yH(in2, end)) - 1);
if extrapPct > 2
    warning(['The largest OLED height reached by system 1 is %.1f %% above the largest sampled height of ' ...
        'system 2; f_2(y) is extrapolated there.'], extrapPct);
end
if max(fitErr2) > 1e-3 * median(E2)
    warning('f_2(y) fit error up to %.4f mm (relative %.2e).', max(fitErr2), max(fitErr2) / median(E2));
end

%% Statistics of the magnification at every field
dM = 100 * (M ./ M(:, :, 1) - 1);                    % change from on-axis (%), n1 x n2 x nF
pctLevels = [1 5 50 95 99];
nP = numel(pctLevels);
fieldStats = zeros(nF, 4);                            % M: mean std min max
fieldPct = zeros(nF, nP);                             % M percentiles
changeStats = zeros(nF, 4);                           % dM (%): mean std min max
changePct = zeros(nF, 2);                             % dM (%): P1 P99
for f = 1:nF
    mf = M(:, :, f);
    mf = mf(:);
    df = dM(:, :, f);
    df = df(:);
    fieldStats(f, :) = [mean(mf), std(mf), min(mf), max(mf)];
    fieldPct(f, :) = quantileSorted(sort(mf), pctLevels / 100);
    changeStats(f, :) = [mean(df), std(df), min(df), max(df)];
    changePct(f, :) = quantileSorted(sort(df), [1 99] / 100);
end
meanShiftEdge = fieldStats(end, 1) - fieldStats(1, 1);
stdIncreaseEdge = fieldStats(end, 2) - fieldStats(1, 2);

%% Command window summary
fprintf('\n===== MAGNIFICATION BY FIELD (%g %% MTF yield sets) =====\n', yieldPct);
fprintf('c = k x p2 / p1 = %g x %g / %g = %.6f;  pairs: %d x %d = %d\n', ...
    pixelScale, pixelPitch_um(2), pixelPitch_um(1), magConst, n1, n2, nPairs);
fprintf('%-6s %10s %10s %9s %10s %10s %10s %10s %10s %10s\n', 'Field', 'Nominal M', 'Mean M', 'Std M', ...
    'Min M', 'Max M', 'dM mean %', 'dM min %', 'dM max %', 'Nom. dM %');
for f = 1:nF
    fprintf('%-6s %10.5f %10.5f %9.5f %10.5f %10.5f %+10.3f %+10.3f %+10.3f %+10.3f\n', fieldNames(f), ...
        mNomField(f), fieldStats(f, 1), fieldStats(f, 2), fieldStats(f, 3), fieldStats(f, 4), ...
        changeStats(f, 1), changeStats(f, 3), changeStats(f, 4), 100 * dDistNom(f));
end

%% Output folder: "BÜYÜTME-DİSTORSİYON ANALİZİ - %90 YIELD" in the system 1 folder
outName = ['B' char(220) 'Y' char(220) 'TME-D' char(304) 'STORS' char(304) 'YON ANAL' char(304) 'Z' char(304) ...
    sprintf(' - %%%g YIELD', yieldPct)];
outFolder = fullfile(sys(1).folder, outName);
if ~isfolder(outFolder)
    mkdir(outFolder);
end

colBlue = [0 0.45 0.74];
fieldColors = lines(nF);
fieldLabel = sprintf('Normalized field of %s (0 = on-axis)', sys(1).name);

%% Figure 1: distortion by field
fig1 = figure('Name', 'Distortion by Field', 'Color', 'w', 'Position', [100 100 950 700]);
for s = 1:2
    subplot(2, 1, s);
    Dy = 100 * sys(s).D(sys(s).inYield, :);
    hold on;
    fill([normFields, fliplr(normFields)], [min(Dy, [], 1), fliplr(max(Dy, [], 1))], colBlue, ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', 'Yield set min – max');
    plot(normFields, mean(Dy, 1), '-o', 'Color', colBlue, 'LineWidth', 1.5, 'MarkerFaceColor', colBlue, ...
        'DisplayName', 'Yield set mean');
    plot(normFields, 100 * sys(s).nominal.D, '--sk', 'LineWidth', 1.3, 'DisplayName', 'Nominal');
    yline(0, ':k', 'HandleVisibility', 'off');
    hold off;
    xlabel(sprintf('Normalized field of %s', sys(s).name));
    ylabel('Distortion F-Tan(\theta) (%)');
    title(sprintf('%s distortion: nominal %+.3f %% at full field', sys(s).name, 100 * sys(s).nominal.D(end)));
    legend('Location', 'best');
    grid on;
    box on;
end
saveFigure(fig1, fullfile(outFolder, 'Distortion_by_Field.png'));

%% Figure 2: magnification by field
fig2 = figure('Name', 'Magnification by Field', 'Color', 'w', 'Position', [120 120 950 560]);
hold on;
fill([normFields, fliplr(normFields)], [fieldStats(:, 3).', fliplr(fieldStats(:, 4).')], colBlue, ...
    'FaceAlpha', 0.12, 'EdgeColor', 'none', 'DisplayName', 'All pairs min – max');
fill([normFields, fliplr(normFields)], [fieldPct(:, 1).', fliplr(fieldPct(:, end).')], colBlue, ...
    'FaceAlpha', 0.30, 'EdgeColor', 'none', 'DisplayName', 'All pairs 1 % – 99 %');
plot(normFields, fieldStats(:, 1), '-o', 'Color', colBlue, 'LineWidth', 1.6, 'MarkerFaceColor', colBlue, ...
    'DisplayName', 'Mean of all pairs');
plot(normFields, mNomField, '--sk', 'LineWidth', 1.3, 'DisplayName', 'Nominal pair');
yline(magTarget, ':k', sprintf('M = %.2f', magTarget), 'LineWidth', 1, 'HandleVisibility', 'off');
for f = 1:nF
    text(normFields(f), fieldStats(f, 4), sprintf('%.4f', fieldStats(f, 1)), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', colBlue);
end
hold off;
xlim([-0.03 1.03]);
xlabel(fieldLabel);
ylabel('Magnification M');
title(sprintf('Magnification on-axis and off-axis (%d pairs, %g%% MTF yield sets)', nPairs, yieldPct));
legend('Location', 'best');
grid on;
box on;
saveFigure(fig2, fullfile(outFolder, 'Magnification_by_Field.png'));

%% Figure 3: change of magnification from on-axis
fig3 = figure('Name', 'Magnification Change', 'Color', 'w', 'Position', [140 140 950 560]);
hold on;
fill([normFields, fliplr(normFields)], [changeStats(:, 3).', fliplr(changeStats(:, 4).')], colBlue, ...
    'FaceAlpha', 0.12, 'EdgeColor', 'none', 'DisplayName', 'All pairs min – max');
fill([normFields, fliplr(normFields)], [changePct(:, 1).', fliplr(changePct(:, 2).')], colBlue, ...
    'FaceAlpha', 0.30, 'EdgeColor', 'none', 'DisplayName', 'All pairs 1 % – 99 %');
plot(normFields, changeStats(:, 1), '-o', 'Color', colBlue, 'LineWidth', 1.6, 'MarkerFaceColor', colBlue, ...
    'DisplayName', 'Mean of all pairs');
plot(normFields, 100 * dDistNom, '--sk', 'LineWidth', 1.3, 'DisplayName', 'Nominal pair');
yline(0, ':k', 'HandleVisibility', 'off');
hold off;
xlim([-0.03 1.03]);
xlabel(fieldLabel);
ylabel('\DeltaM = M(h) / M(0) - 1  (%)');
title('Change of magnification from on-axis (effect of distortion)');
legend('Location', 'best');
grid on;
box on;
saveFigure(fig3, fullfile(outFolder, 'Magnification_Change.png'));

%% Figure 4: magnification histogram of every field
fig4 = figure('Name', 'Magnification Histograms', 'Color', 'w', 'Position', [160 160 950 540]);
mEdges = linspace(min(M(:)), max(M(:)), nHistBins + 1);
hold on;
for f = 1:nF
    mf = M(:, :, f);
    histogram(mf(:), mEdges, 'DisplayStyle', 'stairs', 'EdgeColor', fieldColors(f, :), 'LineWidth', 1.6, ...
        'DisplayName', sprintf('Field %s: mean %.4f, std %.4f', fieldNames(f), fieldStats(f, 1), fieldStats(f, 2)));
end
if magTarget >= mEdges(1) && magTarget <= mEdges(end)
    xline(magTarget, ':k', sprintf('M = %.2f', magTarget), 'HandleVisibility', 'off');
end
hold off;
xlabel('Magnification M');
ylabel('Number of pairs');
title('Magnification distribution on-axis and off-axis');
legend('Location', 'best');
grid on;
box on;
saveFigure(fig4, fullfile(outFolder, 'Magnification_Histograms.png'));

%% Export Excel
outputFile = fullfile(outFolder, 'Distortion_Magnification_Results.xlsx');
if isfile(outputFile)
    delete(outputFile);
end
fieldTag = @(h) strrep(sprintf('%.2f', h), '.', 'p');  % 0.70 -> "0p70" for column names

summaryNames = ["Pairs (n1 x n2)"; sys(1).name + " files in yield set"; sys(2).name + " files in yield set"; ...
    "c = k*p2/p1"; "Nominal M on-axis"; "Nominal M at full field"; ...
    "Nominal magnification change on-axis -> full field"; ...
    "Mean M on-axis"; "Mean M at full field"; "Mean M shift on-axis -> full field"; ...
    "Std M on-axis"; "Std M at full field"; ...
    "Nominal distortion " + sys(1).name + " at full field"; "Nominal distortion " + sys(2).name + " at full field"; ...
    "Largest OLED height reached / largest sampled height - 1"; "Max f_2(y) fit error"];
summaryValues = [nPairs; n1; n2; magConst; mNomField(1); mNomField(end); 100 * dDistNom(end); ...
    fieldStats(1, 1); fieldStats(end, 1); meanShiftEdge; fieldStats(1, 2); fieldStats(end, 2); ...
    100 * sys(1).nominal.D(end); 100 * sys(2).nominal.D(end); extrapPct; max(fitErr2)];
summaryUnits = ["count"; "count"; "count"; ""; ""; ""; "%"; ""; ""; ""; ""; ""; "%"; "%"; "%"; "mm"];
writetable(table(summaryNames, summaryValues, summaryUnits, 'VariableNames', {'Parameter', 'Value', 'Unit'}), ...
    outputFile, 'Sheet', 'Summary');

writetable(table(normFields(:), sys(1).nominal.yH(:), atand(sys(1).nominal.tH(:)), yoNom(:), ...
    mNomField(:), fieldStats(:, 1), fieldStats(:, 2), fieldStats(:, 3), fieldStats(:, 4), ...
    100 * dDistNom(:), changeStats(:, 1), changeStats(:, 2), changeStats(:, 3), changeStats(:, 4), ...
    'VariableNames', {'NormField', 'Nominal_detector_height_mm', 'Nominal_field_angle_deg', ...
    'Nominal_OLED_height_mm', 'Nominal_M', 'Mean_M', 'Std_M', 'Min_M', 'Max_M', ...
    'Nominal_dM_pct', 'Mean_dM_pct', 'Std_dM_pct', 'Min_dM_pct', 'Max_dM_pct'}), ...
    outputFile, 'Sheet', 'Field Stats');

pctTable = table(pctLevels(:), 'VariableNames', {'Percentile'});
for f = 1:nF
    pctTable.(['M_field_' fieldTag(normFields(f))]) = fieldPct(f, :).';
end
writetable(pctTable, outputFile, 'Sheet', 'Percentiles');

distTable = table(normFields(:), 'VariableNames', {'NormField'});
for s = 1:2
    Dy = 100 * sys(s).D(sys(s).inYield, :);
    distTable.(char(sys(s).name + "_nominal_pct")) = 100 * sys(s).nominal.D(:);
    distTable.(char(sys(s).name + "_mean_pct")) = mean(Dy, 1).';
    distTable.(char(sys(s).name + "_std_pct")) = std(Dy, 0, 1).';
    distTable.(char(sys(s).name + "_min_pct")) = min(Dy, [], 1).';
    distTable.(char(sys(s).name + "_max_pct")) = max(Dy, [], 1).';
    distTable.(char(sys(s).name + "_nominal_height_mm")) = sys(s).nominal.yH(:);
end
writetable(distTable, outputFile, 'Sheet', 'Distortion');

for s = 1:2
    fileTable = table(sys(s).fileNames, sys(s).EFL, sys(s).mtfMetric, sys(s).mtfRank, sys(s).inYield, ...
        'VariableNames', {'FileName', 'EFL_paraxial_abs_mm', 'MTF_criterion', 'MTF_rank', 'InYieldSet'});
    for f = 1:nF
        fileTable.(['fEff_' fieldTag(normFields(f)) '_mm']) = sys(s).fEff(:, f);
    end
    for f = 1:nF
        fileTable.(['Distortion_' fieldTag(normFields(f)) '_pct']) = 100 * sys(s).D(:, f);
    end
    writetable(fileTable, outputFile, 'Sheet', char(sys(s).name + " Files"));
end

if writePairsSheet
    if nPairs <= 1048575
        [I, J] = ndgrid(1:n1, 1:n2);
        pairTable = table(names1(I(:)), names2(J(:)), E1(I(:)), E2(J(:)), ...
            'VariableNames', {char("File_" + sys(1).name), char("File_" + sys(2).name), ...
            char("EFL_" + sys(1).name + "_mm"), char("EFL_" + sys(2).name + "_mm")});
        for f = 1:nF
            mf = M(:, :, f);
            pairTable.(['M_field_' fieldTag(normFields(f))]) = mf(:);
        end
        for f = 2:nF
            df = dM(:, :, f);
            pairTable.(['dM_field_' fieldTag(normFields(f)) '_pct']) = df(:);
        end
        writetable(pairTable, outputFile, 'Sheet', 'Pairs');
    else
        warning('%d pairs exceed the Excel row limit; the "Pairs" sheet is not written (see the .mat file).', nPairs);
    end
end

settingNames = ["System 1"; "System 1 folder"; "System 1 nominal file"; "System 1 afocal"; ...
    "System 1 MTF frequency (cyc/mm)"; "System 1 MTF threshold"; "System 1 MTF fields"; ...
    "System 2"; "System 2 folder"; "System 2 nominal file"; "System 2 afocal"; ...
    "System 2 MTF frequency (cyc/mm)"; "System 2 MTF threshold"; "System 2 MTF fields"; ...
    "File prefix"; "Yield level (%)"; "MTF frequency (fraction of Nyquist)"; "Pixel pitch (um)"; ...
    "Pixels"; "Pixel scale k"; "c = k*p2/p1"; "Normalized fields"; "Magnification formula"];
settingValues = strings(numel(settingNames), 1);
for s = 1:2
    o = (s - 1) * 7;
    settingValues(o + (1:7)) = [sys(s).name; string(sys(s).folder); string(sys(s).nominalName); ...
        string(sys(s).isAfocal); string(sys(s).mtfFreq); string(sys(s).threshold); fieldText(mtfFields{s})];
end
settingValues(15:end) = [string(filePrefix); string(yieldPct); string(mtfFreqFrac); ...
    strjoin(string(pixelPitch_um), " / "); strjoin(string(nPixels), " / "); string(pixelScale); ...
    string(magConst); strjoin(fieldNames, ", "); ...
    "M(h) = c * f1_eff(h) / f2(c * y1(h)),  f2(y) = EFL2 + a*y^2 + b*y^4"];
writetable(table(settingNames, settingValues, 'VariableNames', {'Parameter', 'Value'}), ...
    outputFile, 'Sheet', 'Settings');

%% Save raw data
save(fullfile(outFolder, 'Distortion_Magnification_Data.mat'), 'sys', 'names1', 'names2', 'M', 'dM', ...
    'normFields', 'fieldStats', 'fieldPct', 'changeStats', 'mNomField', 'coef2', 'magConst', 'yieldPct', ...
    'mtfFreqFrac', 'pixelPitch_um', 'nPixels', 'pixelScale');

%% PDF report (A4)
pdfFile = fullfile(outFolder, 'Distortion_Magnification_Report.pdf');
oldPages = dir(fullfile(outFolder, 'Distortion_Magnification_Report_page*.pdf'));
for i = 1:numel(oldPages)
    delete(fullfile(outFolder, oldPages(i).name));
end
if isfile(pdfFile)
    delete(pdfFile);
end
if isfile(pdfFile)   % still open in a PDF viewer
    pdfFile = fullfile(outFolder, ['Distortion_Magnification_Report_' ...
        char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.pdf']);
    fprintf('\nThe old report is open in another program; new report: %s\n', pdfFile);
end

n1Name = char(sys(1).name);
n2Name = char(sys(2).name);
modeText = ["nesne sonsuzda: y = görüntüde REAY, açı = nesne uzayında REAB/REAC", ...
            "afokal: y = nesnede (OLED) REAY, açı = görüntü uzayında REAB/REAC"];

% --- Page: summary -------------------------------------------------------
summaryLines = {
    sprintf('Tarih:                 %s', char(datetime('now', 'Format', 'dd.MM.yyyy HH:mm')))
    sprintf('Alanlar (%s):        %s  (0 = on-axis)', n1Name, strjoin(fieldNames, ', '))
    sprintf('c = k × p2 / p1:       %g × %.2f / %.2f = %.6f', pixelScale, pixelPitch_um(2), pixelPitch_um(1), magConst)
    sprintf('Yield seviyesi:        %%%g  (MTF, Nyquist × %g)', yieldPct, mtfFreqFrac)
    sprintf('Çift sayısı:           %d × %d = %d', n1, n2, nPairs)
    ''
    };
for s = 1:2
    summaryLines = [summaryLines; {
        sprintf('SİSTEM %d – %s', s, sys(s).name)
        sprintf('  Klasör:              %s', sys(s).folder)
        sprintf('  Nominal dosya:       %s   |EFL| = %.4f mm', sys(s).nominalName, sys(s).nominal.EFL)
        sprintf('  Okuma:               %s', modeText(sys(s).isAfocal + 1))
        sprintf('  Dosyalar:            %d MC, %d geçerli, %d yield kümesinde', ...
            numel(sys(s).fileNames), sum(sys(s).valid), sum(sys(s).inYield))
        sprintf('  Nominal distorsiyon: %s', strjoin(compose("%s: %+.3f%%", fieldNames(:), 100 * sys(s).nominal.D(:)).', ', '))
        ''}]; %#ok<AGROW>
end
summaryLines = [summaryLines; {
    'BÜYÜTME M (tüm çiftler)'
    '  Alan   Nominal M   Ort. M     Std M     Min M     Max M     %1 M      %99 M'
    }];
for f = 1:nF
    summaryLines{end + 1, 1} = sprintf('  %-5s  %9.5f  %9.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f', fieldNames(f), ...
        mNomField(f), fieldStats(f, 1), fieldStats(f, 2), fieldStats(f, 3), fieldStats(f, 4), ...
        fieldPct(f, 1), fieldPct(f, end)); %#ok<SAGROW>
end
summaryLines = [summaryLines; {
    ''
    'ON-AXIS''E GÖRE DEĞİŞİM ΔM = M(h) / M(0) − 1  (%)'
    '  Alan   Nominal    Ortalama   Std       Min        Max'
    }];
for f = 2:nF
    summaryLines{end + 1, 1} = sprintf('  %-5s  %+8.3f  %+9.3f  %8.3f  %+9.3f  %+9.3f', fieldNames(f), ...
        100 * dDistNom(f), changeStats(f, 1), changeStats(f, 2), changeStats(f, 3), changeStats(f, 4)); %#ok<SAGROW>
end
summaryFig = textPage('DİSTORSİYONLU BÜYÜTME ANALİZ RAPORU', summaryLines, 'Courier New', 8);

% --- Page: assessment ----------------------------------------------------
assess = {
    'ON-AXIS VE OFF-AXIS BÜYÜTME'
    sprintf(['%%%g MTF yield kümelerindeki %d çiftin on-axis büyütmesi ortalama %.5f (std %.5f, %.5f – %.5f). ' ...
     'Tam alanda (%s) ortalama %.5f (std %.5f, %.5f – %.5f).'], yieldPct, nPairs, fieldStats(1, 1), ...
     fieldStats(1, 2), fieldStats(1, 3), fieldStats(1, 4), fieldNames(end), fieldStats(end, 1), ...
     fieldStats(end, 2), fieldStats(end, 3), fieldStats(end, 4))
    ''
    'NOMİNAL DİSTORSİYONUN ETKİSİ'
    sprintf(['Tam alanda nominal distorsiyon: %s %+.3f %%, %s %+.3f %% (her biri kendi tam alanında). Nominal ' ...
     'çiftin büyütmesi on-axis %.5f, tam alanda %.5f; distorsiyon tam alanda büyütmeyi %+.3f %% değiştiriyor.'], ...
     n1Name, 100 * sys(1).nominal.D(end), n2Name, 100 * sys(2).nominal.D(end), mNomField(1), mNomField(end), ...
     100 * dDistNom(end))
    };
if sign(sys(1).nominal.D(end)) == sign(sys(2).nominal.D(end)) && sys(1).nominal.D(end) ~= 0
    assess{end + 1, 1} = ['İki sistemin distorsiyonu aynı işaretli; büyütmede birbirini kısmen götürüyor. ' ...
        'Off-axis büyütmeyi belirleyen iki distorsiyonun farkıdır.'];
else
    assess{end + 1, 1} = ['İki sistemin distorsiyonu zıt işaretli; büyütmedeki etkileri toplanıyor. ' ...
        'Distorsiyonların eşleştirilmesi (aynı işaret ve büyüklük) off-axis büyütme farkını azaltır.'];
end
assess = [assess; {
    ''
    'TOLERANSLARIN ETKİSİ'
    sprintf(['Çiftlerin ortalama büyütmesi on-axis''ten tam alana %+.5f değişiyor (ortalama ΔM %+.3f %%). Bu kısım ' ...
     'nominal distorsiyon farkından gelir ve tüm ürünlerde aynıdır; tasarımda distorsiyonlar eşleştirilerek ' ...
     'veya görüntü işlemede distorsiyon düzeltmesiyle giderilebilir.'], meanShiftEdge, changeStats(end, 1))
    sprintf(['ΔM''nin çiftler arasındaki std''si tam alanda %.3f %% (%.3f %% – %+.3f %%). Bu, toleransların ' ...
     'distorsiyonu dosyadan dosyaya değiştirmesinden gelir; ürüne özeldir ve sabit bir düzeltmeyle giderilemez. ' ...
     'Büyütmenin std''si on-axis %.5f, tam alanda %.5f (%+.5f).'], changeStats(end, 2), changeStats(end, 3), ...
     changeStats(end, 4), fieldStats(1, 2), fieldStats(end, 2), stdIncreaseEdge)
    ''
    'KONTROL'
    ['Okunan nominal distorsiyonu, nominal dosyalarda OpticStudio''nun Field Curvature / Distortion analizindeki ' ...
     'F-Tan(Theta) distorsiyonuyla karşılaştırın (özet sayfasında her alan için verilmiştir). Değerler ' ...
     'uyuşmuyorsa alan tanımı veya ışın okuma yönü (afokal / sonlu) kontrol edilmelidir.']
    }];
if extrapPct > 2
    assess{end + 1, 1} = sprintf(['UYARI: %s''ın tam alanı OLED''de, %s''in örneklenen en büyük yüksekliğinin %%%.1f ' ...
        'üzerine çıkıyor; orada kolimatörün etkin EFL''si uyum eğrisinden dışarı doğru tahmin edildi.'], ...
        n1Name, n2Name, extrapPct);
end
if min(n1, n2) < 50
    assess = [assess; {''; 'UYARI – ÖRNEKLEM BÜYÜKLÜĞÜ'}];
    assess{end + 1, 1} = sprintf(['Yield kümelerinde %d ve %d dosya var; istatistiklerin güvenilirliği çift ' ...
        'sayısıyla değil dosya sayısıyla sınırlıdır. Sistem başına en az 50 – 100 MC dosyası önerilir.'], n1, n2);
end
assessFig = textPage('1. DEĞERLENDİRME', assess, 'Arial', 10);

% --- Page: method --------------------------------------------------------
method = {
    'NEDEN DİSTORSİYON?'
    ['Paraksiyel EFL yalnızca on-axis büyütmeyi verir. Distorsiyonlu bir sistemde gerçek görüntü yüksekliği ' ...
     'y = f · tan θ · (1 + D(θ)) olur; bu yüzden büyütme alan boyunca değişir. İki sistem arka arkaya çalıştığı ' ...
     'için off-axis büyütme yaklaşık iki sistemin distorsiyon farkı kadar değişir: M(θ) ≈ M(0) × (1 + D_1 − D_2).']
    ''
    'ADIM 1 – ETKİN EFL'
    sprintf(['Her MC dosyasında normalize alanlarda (%s) birincil dalga boyunda gerçek ana ışın izlenir ve ' ...
     'f_eff(h) = |y| / |tan θ| hesaplanır. Nesnesi sonsuzda olan dosyada (afokal değil) y görüntü yüzeyinde ' ...
     'REAY, tan θ nesne uzayında REAB / REAC''dır. Afokal dosyada (nesne = OLED) y nesne yüzeyinde REAY, tan θ ' ...
     'görüntü uzayında REAB / REAC''dır. On-axis''te (h = 0) paraksiyel EFL (EFFL) kullanılır. Distorsiyon ' ...
     'D(h) = f_eff(h) / EFL − 1''dir.'], strjoin(fieldNames, ', '))
    ''
    'ADIM 2 – YIELD KÜMELERİ'
    sprintf(['Her sistemde MTF''e göre en iyi %%%g dosya seçilir (büyütme scriptiyle aynı ölçüt). Alan ışınları ' ...
     'izlenemeyen dosyalar geçersiz sayılır.'], yieldPct)
    ''
    'ADIM 3 – ALANA GÖRE BÜYÜTME'
    sprintf(['%s''ın h alanında dedektördeki yükseklik y_1(h), giriş açısı θ_1(h)''dir. Görüntü OLED''e c = %.6f ' ...
     'ölçeğiyle basılır: y_OLED = c × y_1(h). Kolimatörün etkin EFL''si yüksekliğe bağlıdır; %s''in alan ' ...
     'noktalarından her dosya için f_2(y) = EFL_2 + a·y² + b·y⁴ eğrisi uydurulur (y = 0''da tam olarak ' ...
     'paraksiyel EFL). Çıkış açısı tan θ_2 = y_OLED / f_2(y_OLED), büyütme M(h) = tan θ_2 / tan θ_1(h) = ' ...
     'c × f_1,eff(h) / f_2(y_OLED)''dir. İki sistemin alanları farklı olsa da her nokta OLED üzerindeki gerçek ' ...
     'yükseklikte eşleştirilir.'], n1Name, magConst, n2Name)
    ''
    'ADIM 4 – ÇAPRAZLAMA'
    sprintf(['%s yield kümesindeki her dosya %s yield kümesindeki her dosyayla eşleştirilir (%d × %d = %d çift) ' ...
     've her çiftin her alandaki büyütmesi hesaplanır. Her alan için büyütmenin istatistikleri ile on-axis''e ' ...
     'göre değişimi ΔM(h) = M(h) / M(0) − 1 verilir.'], n1Name, n2Name, n1, n2, nPairs)
    ''
    'VARSAYIMLAR'
    ['•  Sistemler dönel simetrik kabul edilir ve alanlar y ekseninde (Hy) örneklenir. Ana ışın birincil ' ...
     'dalga boyunda izlenir; renk farkı (lateral renk) dahil değildir.']
    '•  Görüntü işlemede distorsiyon düzeltmesi yapılmadığı varsayılır.'
    };
methodFig = textPage('2. YÖNTEM', method, 'Arial', 10);

% --- Figure pages --------------------------------------------------------
plotFigs = gobjects(4, 1);
plotFigs(1) = figurePage(fig2, '3. ON-AXIS VE OFF-AXIS BÜYÜTME', {
    'NASIL OKUNUR'
    sprintf(['•  Yatay eksen %s''ın normalize alanıdır (0 = on-axis). Dikey eksen büyütmedir. Mavi çizgi tüm ' ...
     'çiftlerin ortalaması (değeri noktaların üstünde yazılı), koyu mavi alan çiftlerin %%1 – %%99''u, açık mavi ' ...
     'alan en küçük ve en büyük değerdir. Kesikli siyah çizgi nominal çifttir, noktalı çizgi M = %.2f.'], n1Name, magTarget)
    'NASIL YORUMLANIR'
    ['•  Ortalama çizginin on-axis''ten kenara eğimi distorsiyon farkının etkisidir; tüm ürünlerde aynıdır ve ' ...
     'tasarımla veya görüntü işlemeyle düzeltilebilir.']
    '•  Bantların kenara doğru genişlemesi toleransların distorsiyonu değiştirmesinden gelir; ürüne özeldir.'
    }, false);
plotFigs(2) = figurePage(fig3, '4. ON-AXIS''E GÖRE BÜYÜTME DEĞİŞİMİ', {
    'NASIL OKUNUR'
    ['•  Her çiftin off-axis büyütmesinin kendi on-axis büyütmesine göre yüzde değişimidir: ' ...
     'ΔM(h) = M(h) / M(0) − 1. On-axis''te tanım gereği sıfırdır. Bantlar ve çizgiler önceki grafikle aynı anlamdadır.']
    'NASIL YORUMLANIR'
    ['•  Bu grafik EFL saçılımını ayıklar ve yalnız distorsiyonun etkisini gösterir. Ortalama çizgi nominal ' ...
     'distorsiyon farkını, bant genişliği toleransların distorsiyona etkisini verir.']
    ['•  Örneğin tam alanda ortalama −%0.8 ise, on-axis büyütmesi 1.000 olan bir ürünün kenarda büyütmesi ' ...
     'yaklaşık 0.992''dir.']
    }, false);
plotFigs(3) = figurePage(fig4, '5. ALANLARIN BÜYÜTME DAĞILIMI', {
    'NASIL OKUNUR'
    ['•  Her renkli basamak çizgisi bir alanın büyütme dağılımıdır (tüm çiftler). Lejantta her alanın ' ...
     'ortalaması ve std''si yazılıdır.']
    'NASIL YORUMLANIR'
    ['•  Off-axis dağılımları on-axis''e göre kaymışsa sebep nominal distorsiyon farkıdır; genişlemişse sebep ' ...
     'toleransların distorsiyona etkisidir.']
    }, false);
plotFigs(4) = figurePage(fig1, '6. ALANA GÖRE DİSTORSİYON', {
    'NASIL OKUNUR'
    ['•  Her panel bir sistemdir. Yatay eksen o sistemin normalize alanı, dikey eksen F-Tan(θ) ' ...
     'distorsiyonudur (%). Kesikli siyah çizgi nominal tasarım, mavi çizgi yield kümesinin ortalaması, ' ...
     'açık mavi alan yield kümesindeki en küçük ve en büyük değerlerdir.']
    'NASIL YORUMLANIR'
    ['•  İki sistemin tam alandaki nominal distorsiyonları aynı işaretli ve yakın büyüklükteyse büyütmeye ' ...
     'etkileri birbirini götürür; zıt işaretliyse toplanır.']
    ['•  Açık mavi alanın genişliği, toleransların distorsiyonu ne kadar değiştirdiğidir. Geniş ise off-axis ' ...
     'büyütme üründen ürüne değişir.']
    }, false);

% --- Assemble pages ------------------------------------------------------
allFigs = [summaryFig; assessFig; methodFig; plotFigs];
nPages = numel(allFigs);
headerText = sprintf('Distorsiyonlu Büyütme Analizi  |  %s: %s  |  %s: %s', ...
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
        print(allFigs(p), fullfile(outFolder, sprintf('Distortion_Magnification_Report_page%d.pdf', p)), '-dpdf');
    end
    fprintf('Report saved as separate pages.\n');
end
close(allFigs);

%% Text summary
summaryFile = fullfile(outFolder, 'Analysis_Summary.txt');
fid = fopen(summaryFile, 'w', 'n', 'UTF-8');
if fid > 0
    fprintf(fid, '%s\r\n', 'DİSTORSİYONLU BÜYÜTME ANALİZ ÖZETİ', '');
    fprintf(fid, '%s\r\n', summaryLines{:});
    fprintf(fid, '\r\n%s\r\n\r\n', 'DEĞERLENDİRME');
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
function [EFL, fEff, yH, tH, metric, nominal, isAfocal, mtfFreq] = computeSystemField(TheSystem, S, ...
    normFields, mtfFreqFrac, mtfFields)
% Paraxial |EFL|, field data (effective EFL, height, tan of angle at each
% normalized field) and MTF criterion of every MC file of one system, and the
% same field data of the nominal file.
nFiles = numel(S.fileNames);
nF = numel(normFields);
EFL = nan(nFiles, 1);
fEff = nan(nFiles, nF);
yH = nan(nFiles, nF);
tH = nan(nFiles, nF);
metric = nan(nFiles, 1);
mtfFreq = mtfFreqFrac * 1000 / (2 * S.pitch);   % cyc/mm

% Nominal file
if ~TheSystem.LoadFile(char(fullfile(S.folder, S.nominalName)), false)
    error('%s: OpticStudio could not load the nominal file %s', S.name, S.nominalName);
end
isAfocal = logical(TheSystem.SystemData.Aperture.AFocalImageSpace);
nominal.EFL = abs(readEfl(TheSystem));
if ~isfinite(nominal.EFL) || nominal.EFL == 0
    error('%s: could not compute the EFL of the nominal file %s', S.name, S.nominalName);
end
checkObjectDistance(TheSystem, isAfocal, S.name);
[nominal.fEff, nominal.yH, nominal.tH] = fieldEfl(TheSystem, normFields, isAfocal, nominal.EFL);
if any(~isfinite(nominal.fEff))
    error('%s: the chief ray could not be traced at every field of the nominal file.', S.name);
end

freqScale = 1;   % Zemax frequency -> cyc/mm
if isAfocal
    freqScale = 1000 / nominal.EFL;   % afocal MTF is in cyc/mrad
end
fprintf('%s nominal %s: |EFL| = %.4f mm, afocal = %d, MTF criterion at %.2f cyc/mm\n', ...
    S.name, S.nominalName, nominal.EFL, isAfocal, mtfFreq);

wb = waitbar(0, sprintf('%s: computing...', S.name), 'Name', 'Distortion Magnification Analysis');
closeWb = onCleanup(@() delete(wb(isvalid(wb))));

for k = 1:nFiles
    try
        if ~TheSystem.LoadFile(char(fullfile(S.folder, S.fileNames(k))), false)
            error('OpticStudio could not load this file.');
        end
        EFL(k) = abs(readEfl(TheSystem));
        [fEff(k, :), yH(k, :), tH(k, :)] = fieldEfl(TheSystem, normFields, isAfocal, EFL(k));
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

function [fEff, yH, tH] = fieldEfl(TheSystem, normFields, isAfocal, efl0)
% Real chief ray (Hx = 0, Hy = h, Px = Py = 0, primary wavelength) at each
% normalized field h. yH = |height|, tH = |tan(angle)|, fEff = yH / tH.
%   Not afocal (object at infinity): height on the image surface, angle in
%     object space (direction cosines at surface 0).
%   Afocal (object = OLED): height on the object surface, angle in image space.
% At h = 0 the paraxial EFL efl0 is used.
imgSurf = TheSystem.LDE.NumberOfSurfaces - 1;
wave = primaryWave(TheSystem);
nF = numel(normFields);
fEff = nan(1, nF);
yH = zeros(1, nF);
tH = zeros(1, nF);

if isAfocal
    hSurf = 0;
    aSurf = imgSurf;
else
    hSurf = imgSurf;
    aSurf = 0;
end

for f = 1:nF
    h = normFields(f);
    if h == 0
        fEff(f) = efl0;
        continue;
    end
    y = TheSystem.MFE.GetOperandValue(ZOSAPI.Editors.MFE.MeritOperandType.REAY, hSurf, wave, 0, h, 0, 0, 0, 0);
    m = TheSystem.MFE.GetOperandValue(ZOSAPI.Editors.MFE.MeritOperandType.REAB, aSurf, wave, 0, h, 0, 0, 0, 0);
    n = TheSystem.MFE.GetOperandValue(ZOSAPI.Editors.MFE.MeritOperandType.REAC, aSurf, wave, 0, h, 0, 0, 0, 0);
    yH(f) = abs(y);
    tH(f) = abs(m / n);
    if yH(f) > 0 && tH(f) > 0 && isfinite(yH(f)) && isfinite(tH(f))
        fEff(f) = yH(f) / tH(f);
    end
end
end

function w = primaryWave(TheSystem)
% Number of the primary wavelength (1 if it cannot be read)
w = 1;
try
    wls = TheSystem.SystemData.Wavelengths;
    for i = 1:wls.NumberOfWavelengths
        if wls.GetWavelength(i).IsPrimary
            w = i;
            return;
        end
    end
catch
end
end

function checkObjectDistance(TheSystem, isAfocal, sysLabel)
% Warn if the object distance does not match the way the field data is read
try
    objThickness = TheSystem.LDE.GetSurfaceAt(0).Thickness;
    atInfinity = ~isfinite(objThickness) || abs(objThickness) > 1e9;
    if isAfocal && atInfinity
        warning(['%s is in afocal mode with the object at infinity (afocal telescope). The field data is read ' ...
            'as object height -> image angle, which is not meaningful here.'], sysLabel);
    elseif ~isAfocal && ~atInfinity
        warning(['%s is not in afocal mode and its object is at a finite distance. The field angle is read from ' ...
            'the chief ray direction in object space; check that this is the intended model.'], sysLabel);
    end
catch
end
end

function [coef, maxErr] = fitFieldEfl(y, fEff, efl0)
% Effective EFL as a function of height: fEff(y) = efl0 + a*y^2 + b*y^4 (least
% squares over the field points with y > 0; exact paraxial EFL at y = 0).
use = y > 0 & isfinite(fEff);
coef = [0 0];
maxErr = 0;
if ~any(use)
    return;
end
y = y(use);
dfe = fEff(use) - efl0;
if nnz(use) >= 2
    A = [y(:).^2, y(:).^4];
else
    A = y(:).^2;
end
c = A \ dfe(:);
coef(1:numel(c)) = c;
maxErr = max(abs(A * c - dfe(:)));
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
