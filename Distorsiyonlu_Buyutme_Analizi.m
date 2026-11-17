%% DİSTORSİYONLU BÜYÜTME ANALİZİ (ZOS-API, MATLAB)
% NE YAPAR
%   Birlikte çalışan iki sistemin (1: MWIR objektif + dedektör, 2: OLED +
%   kolimatör) toplam büyütmesini, distorsiyonu da dahil ederek, görüş
%   alanının birkaç noktasında (normFields, örn. 0, 0.5, 0.7, 1 normalize alan)
%   hesaplar. Her sistemde MTF'e göre %yieldPct yield kümesine giren dosyalar
%   seçilir; 1. sistemin yield kümesindeki her dosya 2. sistemin yield
%   kümesindeki her dosyayla eşleştirilir (n1 x n2 çift). Her çiftin on-axis
%   (0) ve off-axis alanlardaki büyütmesi, y ve x eksenlerinde ayrı ayrı
%   hesaplanır ve alan / eksen bazında büyütme değerleri / dağılımları ile
%   iki eksen arasındaki fark (tilt / decenter kaynaklı asimetri) verilir.
%   Geçti / kaldı değerlendirmesi yapılmaz.
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
%        - Her normalize alanda ve her eksende (fieldAxes) gerçek ana ışının
%          (Px = Py = 0, birincil dalga boyu) yüksekliği ve açısı:
%              y ekseni: (Hx, Hy) = (0, h), yükseklik REAY, açı REAB / REAC
%              x ekseni: (Hx, Hy) = (h, 0), yükseklik REAX, açı REAA / REAC
%          Dosya afokal modda değilse (nesne sonsuzda, görüntü düzlemde):
%              yükseklik görüntü yüzeyinde, açı nesne uzayında (yüzey 0)
%          Dosya afokal moddaysa (nesne = OLED, çıkış açısal):
%              yükseklik nesne yüzeyinde, açı görüntü uzayında
%          Tilt / decenter toleransları on-axis ana ışını da kaydırır
%          (boresight). Bu kayma büyütmeyi değiştirmediği için etkin EFL
%          on-axis ana ışına göre farkla hesaplanır:
%              f_eff(h) = |y(h) − y(0)| / |tan θ(h) − tan θ(0)|
%          (h = 0 için paraksiyel EFL)
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
%   Her alan ve eksen için: nominal çiftin M'si; tüm çiftlerin M ortalaması,
%   std, min, max ve yüzdelikleri; ΔM'nin aynı istatistikleri.
%   Eksen asimetrisi (aynı çift, aynı alan):  M_x / M_y − 1  (%)
%
% AYARLAR (User settings)
%   sysName, pixelPitch_um, nPixels, pixelScale, filePrefix, yieldPct,
%   mtfFreqFrac, mtfFields, normFields, fieldAxes, magTarget (yalnız referans çizgisi),
%   writePairsSheet, nHistBins
%
% ÇIKTILAR (1. sistemin klasöründe "BÜYÜTME-DİSTORSİYON ANALİZİ - %90 YIELD")
%   Distortion_by_Field.png      : iki sistemin alana göre distorsiyonu
%                                  (nominal, yield kümesi ortalaması ve aralığı)
%   Magnification_by_Field.png   : alana göre büyütme (ortalama, %1-%99, min-max)
%   Magnification_Change.png     : on-axis'e göre büyütme değişimi ΔM(h)
%   Magnification_Histograms.png : her alanın büyütme dağılımı (eksen başına)
%   Axis_Asymmetry.png           : x – y ekseni büyütme farkı
%   *.fig                        : grafiklerin düzenlenebilir MATLAB dosyaları
%   Distortion_Magnification_Results.xlsx : "Summary", "Field Stats",
%                                  "Percentiles", "Distortion", "Axis
%                                  Asymmetry", sistem başına
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
normFields      = [0 0.5 0.7 1];     % normalized fields of system 1 where M is evaluated
fieldAxes       = ["y", "x"];        % field axes: "y" = (Hx, Hy) = (0, h), "x" = (h, 0)
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
fieldAxes = lower(string(fieldAxes));
nAx = numel(fieldAxes);
if nAx == 0 || any(~ismember(fieldAxes, ["x", "y"]))
    error('fieldAxes must contain "y" and / or "x".');
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

%% Compute EFL, field data (both axes) and MTF of both systems with OpticStudio
TheApplication = InitConnection();

if isempty(TheApplication)
    error('Could not connect to OpticStudio.');
end

tStart = tic;

try
    TheSystem = TheApplication.PrimarySystem;

    for s = 1:2
        [sys(s).EFL, sys(s).fEff, sys(s).yH, sys(s).tH, sys(s).mtfMetric, sys(s).nominal, ...
            sys(s).isAfocal, sys(s).mtfFreq] = computeSystemField(TheSystem, sys(s), normFields, ...
            fieldAxes, mtfFreqFrac, mtfFields{s});
    end
catch ME
    CleanupConnection(TheApplication);
    rethrow(ME);
end

CleanupConnection(TheApplication);
fprintf('Computation time: %.1f s\n', toc(tStart));

%% Yield set of each system (files with a valid EFL, valid field data and valid MTF)
for s = 1:2
    nFiles = numel(sys(s).fileNames);
    fe = reshape(sys(s).fEff, nFiles, []);
    fieldOk = all(isfinite(fe), 2) & all(fe > 0, 2);
    eflForYield = sys(s).EFL;
    eflForYield(~fieldOk) = NaN;
    [sys(s).valid, sys(s).inYield, sys(s).threshold, sys(s).mtfRank] = ...
        selectYield(eflForYield, sys(s).mtfMetric, yieldPct);

    if ~any(sys(s).inYield)
        error('%s: no file has a valid EFL, valid field data and a valid MTF.', sys(s).name);
    end

    sys(s).D = sys(s).fEff ./ sys(s).EFL - 1;                  % distortion: file x field x axis
    sys(s).nominal.D = sys(s).nominal.fEff / sys(s).nominal.EFL - 1;

    fprintf('%s yield %g %%: MTF >= %.4f at %.2f cyc/mm -> %d of %d valid files\n', ...
        sys(s).name, yieldPct, sys(s).threshold, sys(s).mtfFreq, sum(sys(s).inYield), sum(sys(s).valid));
    for a = 1:nAx
        fprintf('%s nominal distortion (%s): %s\n', sys(s).name, fieldAxes(a), ...
            strjoin(compose("h=%.2f: %+.3f %%", normFields(:), 100 * sys(s).nominal.D(1, :, a).').', ', '));
    end
end

%% Field-dependent magnification of every pair, for every axis
in1 = sys(1).inYield;
in2 = sys(2).inYield;
names1 = sys(1).fileNames(in1);
names2 = sys(2).fileNames(in2);
n1 = sum(in1);
n2 = sum(in2);
nPairs = n1 * n2;

E1 = sys(1).EFL(in1);                 % n1 x 1 paraxial |EFL|
E2 = sys(2).EFL(in2);                 % n2 x 1

M = zeros(n1, n2, nF, nAx);           % magnification: pair x field x axis
coef2 = zeros(n2, 2, nAx);            % f_2(y) = EFL_2 + a*y^2 + b*y^4 of every system 2 file and axis
fitErr2 = zeros(n2, nAx);
mNomField = zeros(nF, nAx);
yoNom = zeros(nF, nAx);
extrapPct = zeros(1, nAx);
for a = 1:nAx
    F1 = sys(1).fEff(in1, :, a);      % n1 x nF effective EFL of system 1
    Y1 = sys(1).yH(in1, :, a);        % n1 x nF detector height (from the on-axis chief ray)
    Y2 = sys(2).yH(in2, :, a);
    F2pts = sys(2).fEff(in2, :, a);
    for j = 1:n2
        [coef2(j, :, a), fitErr2(j, a)] = fitFieldEfl(Y2(j, :), F2pts(j, :), E2(j));
    end
    nomCoef2 = fitFieldEfl(sys(2).nominal.yH(1, :, a), sys(2).nominal.fEff(1, :, a), sys(2).nominal.EFL);

    % OLED height of every system 1 field point, and f_2 there: n1 x n2 x nF
    Yo = magConst * reshape(Y1, n1, 1, nF);
    F2 = E2.' + coef2(:, 1, a).' .* Yo.^2 + coef2(:, 2, a).' .* Yo.^4;
    Ma = magConst * reshape(F1, n1, 1, nF) ./ F2;
    Ma(:, :, 1) = magConst * E1 ./ E2.';                % on-axis: paraxial (same as the analysis without distortion)
    M(:, :, :, a) = Ma;

    % Nominal pair
    yoNom(:, a) = magConst * sys(1).nominal.yH(1, :, a).';
    f2Nom = sys(2).nominal.EFL + nomCoef2(1) * yoNom(:, a).^2 + nomCoef2(2) * yoNom(:, a).^4;
    mNomField(:, a) = magConst * sys(1).nominal.fEff(1, :, a).' ./ f2Nom;
    mNomField(1, a) = magConst * sys(1).nominal.EFL / sys(2).nominal.EFL;

    % Extrapolation check: OLED heights beyond the largest sampled height of system 2
    extrapPct(a) = 100 * (max(Yo(:)) / min(Y2(:, end)) - 1);
    if extrapPct(a) > 2
        warning(['%s axis: the largest OLED height reached by system 1 is %.1f %% above the largest sampled ' ...
            'height of system 2; f_2(y) is extrapolated there.'], fieldAxes(a), extrapPct(a));
    end
end
if max(fitErr2(:)) > 1e-3 * median(E2)
    warning('f_2(y) fit error up to %.4f mm (relative %.2e).', max(fitErr2(:)), max(fitErr2(:)) / median(E2));
end
dDistNom = mNomField ./ mNomField(1, :) - 1;            % nF x nAx, change caused by the nominal distortion

%% Statistics of the magnification at every field and axis
dM = 100 * (M ./ M(:, :, 1, :) - 1);                    % change from on-axis (%)
pctLevels = [1 5 50 95 99];
nP = numel(pctLevels);
fieldStats = zeros(nF, 4, nAx);                          % M: mean std min max
fieldPct = zeros(nF, nP, nAx);                           % M percentiles
changeStats = zeros(nF, 4, nAx);                         % dM (%): mean std min max
changePct = zeros(nF, 2, nAx);                           % dM (%): P1 P99
for a = 1:nAx
    for f = 1:nF
        mf = M(:, :, f, a);
        mf = mf(:);
        df = dM(:, :, f, a);
        df = df(:);
        fieldStats(f, :, a) = [mean(mf), std(mf), min(mf), max(mf)];
        fieldPct(f, :, a) = quantileSorted(sort(mf), pctLevels / 100);
        changeStats(f, :, a) = [mean(df), std(df), min(df), max(df)];
        changePct(f, :, a) = quantileSorted(sort(df), [1 99] / 100);
    end
end

% x - y asymmetry: M of axis 2 relative to axis 1 for the same pair and field (%)
hasAniso = nAx >= 2;
if hasAniso
    aniso = 100 * (M(:, :, :, 2) ./ M(:, :, :, 1) - 1);
    anisoStats = zeros(nF, 6);                           % mean std min max P1 P99
    for f = 1:nF
        af = aniso(:, :, f);
        af = af(:);
        q = quantileSorted(sort(af), [1 99] / 100);
        anisoStats(f, :) = [mean(af), std(af), min(af), max(af), q(1), q(2)];
    end
    anisoNom = 100 * (mNomField(:, 2) ./ mNomField(:, 1) - 1);
end

%% Command window summary
fprintf('\n===== MAGNIFICATION BY FIELD (%g %% MTF yield sets) =====\n', yieldPct);
fprintf('c = k x p2 / p1 = %g x %g / %g = %.6f;  pairs: %d x %d = %d\n', ...
    pixelScale, pixelPitch_um(2), pixelPitch_um(1), magConst, n1, n2, nPairs);
for a = 1:nAx
    fprintf('\nAxis %s:\n', fieldAxes(a));
    fprintf('%-6s %10s %10s %9s %10s %10s %10s %10s %10s %10s\n', 'Field', 'Nominal M', 'Mean M', 'Std M', ...
        'Min M', 'Max M', 'dM mean %', 'dM min %', 'dM max %', 'Nom. dM %');
    for f = 1:nF
        fprintf('%-6s %10.5f %10.5f %9.5f %10.5f %10.5f %+10.3f %+10.3f %+10.3f %+10.3f\n', fieldNames(f), ...
            mNomField(f, a), fieldStats(f, 1, a), fieldStats(f, 2, a), fieldStats(f, 3, a), fieldStats(f, 4, a), ...
            changeStats(f, 1, a), changeStats(f, 3, a), changeStats(f, 4, a), 100 * dDistNom(f, a));
    end
end
if hasAniso
    fprintf('\nAsymmetry M_%s / M_%s - 1 (%%):\n', fieldAxes(2), fieldAxes(1));
    fprintf('%-6s %10s %10s %9s %10s %10s\n', 'Field', 'Nominal', 'Mean', 'Std', 'Min', 'Max');
    for f = 1:nF
        fprintf('%-6s %+10.4f %+10.4f %9.4f %+10.4f %+10.4f\n', fieldNames(f), anisoNom(f), anisoStats(f, 1), ...
            anisoStats(f, 2), anisoStats(f, 3), anisoStats(f, 4));
    end
end

%% Output folder: "BÜYÜTME-DİSTORSİYON ANALİZİ - %90 YIELD" in the system 1 folder
outName = ['B' char(220) 'Y' char(220) 'TME-D' char(304) 'STORS' char(304) 'YON ANAL' char(304) 'Z' char(304) ...
    sprintf(' - %%%g YIELD', yieldPct)];
outFolder = fullfile(sys(1).folder, outName);
if ~isfolder(outFolder)
    mkdir(outFolder);
end

axisPalette = [0 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];
axisColors = axisPalette(mod(0:nAx - 1, size(axisPalette, 1)) + 1, :);
axisMarkers = {'o', 's', '^'};
fieldColors = lines(nF);
fieldLabel = sprintf('Normalized field of %s (0 = on-axis)', sys(1).name);
axLabel = @(a) sprintf('%s axis', fieldAxes(a));

%% Figure 1: distortion by field
fig1 = figure('Name', 'Distortion by Field', 'Color', 'w', 'Position', [100 100 950 700]);
for s = 1:2
    subplot(2, 1, s);
    hold on;
    for a = 1:nAx
        Dy = 100 * sys(s).D(sys(s).inYield, :, a);
        fill([normFields, fliplr(normFields)], [min(Dy, [], 1), fliplr(max(Dy, [], 1))], axisColors(a, :), ...
            'FaceAlpha', 0.12, 'EdgeColor', 'none', 'DisplayName', [axLabel(a) ': yield set min – max']);
        plot(normFields, mean(Dy, 1), ['-' axisMarkers{a}], 'Color', axisColors(a, :), 'LineWidth', 1.5, ...
            'MarkerFaceColor', axisColors(a, :), 'DisplayName', [axLabel(a) ': yield set mean']);
        plot(normFields, 100 * sys(s).nominal.D(1, :, a), ['--' axisMarkers{a}], 'Color', 'k', 'LineWidth', 1.1, ...
            'DisplayName', [axLabel(a) ': nominal']);
    end
    yline(0, ':k', 'HandleVisibility', 'off');
    hold off;
    xlabel(sprintf('Normalized field of %s', sys(s).name));
    ylabel('Distortion F-Tan(\theta) (%)');
    title(sprintf('%s distortion: nominal %+.3f %% at full field (%s)', sys(s).name, ...
        100 * sys(s).nominal.D(1, end, 1), fieldAxes(1)));
    legend('Location', 'best');
    grid on;
    box on;
end
saveFigure(fig1, fullfile(outFolder, 'Distortion_by_Field.png'));

%% Figure 2: magnification by field
fig2 = figure('Name', 'Magnification by Field', 'Color', 'w', 'Position', [120 120 950 560]);
hold on;
for a = 1:nAx
    fill([normFields, fliplr(normFields)], [fieldPct(:, 1, a).', fliplr(fieldPct(:, end, a).')], axisColors(a, :), ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none', 'DisplayName', [axLabel(a) ': 1 % – 99 %']);
    plot(normFields, fieldStats(:, 3, a), ':', 'Color', axisColors(a, :), 'HandleVisibility', 'off');
    plot(normFields, fieldStats(:, 4, a), ':', 'Color', axisColors(a, :), 'DisplayName', [axLabel(a) ': min / max']);
    plot(normFields, fieldStats(:, 1, a), ['-' axisMarkers{a}], 'Color', axisColors(a, :), 'LineWidth', 1.6, ...
        'MarkerFaceColor', axisColors(a, :), 'DisplayName', [axLabel(a) ': mean']);
    plot(normFields, mNomField(:, a), ['--' axisMarkers{a}], 'Color', 'k', 'LineWidth', 1.1, ...
        'DisplayName', [axLabel(a) ': nominal pair']);
end
yline(magTarget, ':k', sprintf('M = %.2f', magTarget), 'LineWidth', 1, 'HandleVisibility', 'off');
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
for a = 1:nAx
    fill([normFields, fliplr(normFields)], [changePct(:, 1, a).', fliplr(changePct(:, 2, a).')], axisColors(a, :), ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none', 'DisplayName', [axLabel(a) ': 1 % – 99 %']);
    plot(normFields, changeStats(:, 3, a), ':', 'Color', axisColors(a, :), 'HandleVisibility', 'off');
    plot(normFields, changeStats(:, 4, a), ':', 'Color', axisColors(a, :), 'DisplayName', [axLabel(a) ': min / max']);
    plot(normFields, changeStats(:, 1, a), ['-' axisMarkers{a}], 'Color', axisColors(a, :), 'LineWidth', 1.6, ...
        'MarkerFaceColor', axisColors(a, :), 'DisplayName', [axLabel(a) ': mean']);
    plot(normFields, 100 * dDistNom(:, a), ['--' axisMarkers{a}], 'Color', 'k', 'LineWidth', 1.1, ...
        'DisplayName', [axLabel(a) ': nominal pair']);
end
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

%% Figure 4: magnification histogram of every field, one panel per axis
fig4 = figure('Name', 'Magnification Histograms', 'Color', 'w', 'Position', [160 60 950 300 + 300 * nAx]);
mEdges = linspace(min(M(:)), max(M(:)), nHistBins + 1);
for a = 1:nAx
    subplot(nAx, 1, a);
    hold on;
    for f = 1:nF
        mf = M(:, :, f, a);
        histogram(mf(:), mEdges, 'DisplayStyle', 'stairs', 'EdgeColor', fieldColors(f, :), 'LineWidth', 1.6, ...
            'DisplayName', sprintf('Field %s: mean %.4f, std %.4f', fieldNames(f), fieldStats(f, 1, a), ...
            fieldStats(f, 2, a)));
    end
    if magTarget >= mEdges(1) && magTarget <= mEdges(end)
        xline(magTarget, ':k', sprintf('M = %.2f', magTarget), 'HandleVisibility', 'off');
    end
    hold off;
    xlabel('Magnification M');
    ylabel('Number of pairs');
    title(sprintf('Magnification distribution, %s', axLabel(a)));
    legend('Location', 'best');
    grid on;
    box on;
end
saveFigure(fig4, fullfile(outFolder, 'Magnification_Histograms.png'));

%% Figure 5: x - y asymmetry
if hasAniso
    fig5 = figure('Name', 'Axis Asymmetry', 'Color', 'w', 'Position', [180 180 950 560]);
    hold on;
    fill([normFields, fliplr(normFields)], [anisoStats(:, 3).', fliplr(anisoStats(:, 4).')], axisColors(2, :), ...
        'FaceAlpha', 0.12, 'EdgeColor', 'none', 'DisplayName', 'All pairs min – max');
    fill([normFields, fliplr(normFields)], [anisoStats(:, 5).', fliplr(anisoStats(:, 6).')], axisColors(2, :), ...
        'FaceAlpha', 0.30, 'EdgeColor', 'none', 'DisplayName', 'All pairs 1 % – 99 %');
    plot(normFields, anisoStats(:, 1), '-o', 'Color', axisColors(2, :), 'LineWidth', 1.6, ...
        'MarkerFaceColor', axisColors(2, :), 'DisplayName', 'Mean of all pairs');
    plot(normFields, anisoNom, '--sk', 'LineWidth', 1.1, 'DisplayName', 'Nominal pair');
    yline(0, ':k', 'HandleVisibility', 'off');
    hold off;
    xlim([-0.03 1.03]);
    xlabel(fieldLabel);
    ylabel(sprintf('M_%s / M_%s - 1  (%%)', fieldAxes(2), fieldAxes(1)));
    title(sprintf('Asymmetry between the %s and %s axes (tilt / decenter effect)', fieldAxes(2), fieldAxes(1)));
    legend('Location', 'best');
    grid on;
    box on;
    saveFigure(fig5, fullfile(outFolder, 'Axis_Asymmetry.png'));
end

%% Export Excel
outputFile = fullfile(outFolder, 'Distortion_Magnification_Results.xlsx');
if isfile(outputFile)
    delete(outputFile);
end
fieldTag = @(h) strrep(sprintf('%.2f', h), '.', 'p');  % 0.70 -> "0p70" for column names

summaryNames = ["Pairs (n1 x n2)"; sys(1).name + " files in yield set"; sys(2).name + " files in yield set"; ...
    "c = k*p2/p1"; "Nominal M on-axis"; "Mean M on-axis"; "Std M on-axis"];
summaryValues = [nPairs; n1; n2; magConst; mNomField(1, 1); fieldStats(1, 1, 1); fieldStats(1, 2, 1)];
summaryUnits = ["count"; "count"; "count"; ""; ""; ""; ""];
for a = 1:nAx
    tag = " (" + fieldAxes(a) + " axis)";
    summaryNames = [summaryNames; "Nominal M at full field" + tag; ...
        "Nominal magnification change on-axis -> full field" + tag; "Mean M at full field" + tag; ...
        "Std M at full field" + tag; "Mean dM at full field" + tag; ...
        "Nominal distortion " + sys(1).name + " at full field" + tag; ...
        "Nominal distortion " + sys(2).name + " at full field" + tag; ...
        "Largest OLED height reached / largest sampled height - 1" + tag]; %#ok<AGROW>
    summaryValues = [summaryValues; mNomField(end, a); 100 * dDistNom(end, a); fieldStats(end, 1, a); ...
        fieldStats(end, 2, a); changeStats(end, 1, a); 100 * sys(1).nominal.D(1, end, a); ...
        100 * sys(2).nominal.D(1, end, a); extrapPct(a)]; %#ok<AGROW>
    summaryUnits = [summaryUnits; ""; "%"; ""; ""; "%"; "%"; "%"; "%"]; %#ok<AGROW>
end
if hasAniso
    summaryNames = [summaryNames; "Asymmetry M_x/M_y - 1 at full field: nominal"; ...
        "Asymmetry at full field: mean"; "Asymmetry at full field: std"; ...
        "Asymmetry at full field: min"; "Asymmetry at full field: max"];
    summaryValues = [summaryValues; anisoNom(end); anisoStats(end, 1); anisoStats(end, 2); ...
        anisoStats(end, 3); anisoStats(end, 4)];
    summaryUnits = [summaryUnits; "%"; "%"; "%"; "%"; "%"];
end
summaryNames = [summaryNames; "Max f_2(y) fit error"];
summaryValues = [summaryValues; max(fitErr2(:))];
summaryUnits = [summaryUnits; "mm"];
writetable(table(summaryNames, summaryValues, summaryUnits, 'VariableNames', {'Parameter', 'Value', 'Unit'}), ...
    outputFile, 'Sheet', 'Summary');

[FF, AA] = ndgrid(1:nF, 1:nAx);
fsT = table(fieldAxes(AA(:)).', normFields(FF(:)).', 'VariableNames', {'Axis', 'NormField'});
fsT.Nominal_detector_height_mm = reshape(sys(1).nominal.yH(1, :, :), [], 1);
fsT.Nominal_OLED_height_mm = yoNom(:);
fsT.Nominal_M = mNomField(:);
mS = reshape(permute(fieldStats, [1 3 2]), [], 4);
cS = reshape(permute(changeStats, [1 3 2]), [], 4);
fsT.Mean_M = mS(:, 1); fsT.Std_M = mS(:, 2); fsT.Min_M = mS(:, 3); fsT.Max_M = mS(:, 4);
fsT.Nominal_dM_pct = 100 * dDistNom(:);
fsT.Mean_dM_pct = cS(:, 1); fsT.Std_dM_pct = cS(:, 2); fsT.Min_dM_pct = cS(:, 3); fsT.Max_dM_pct = cS(:, 4);
writetable(fsT, outputFile, 'Sheet', 'Field Stats');

pctTable = table(pctLevels(:), 'VariableNames', {'Percentile'});
for a = 1:nAx
    for f = 1:nF
        pctTable.(['M_' char(fieldAxes(a)) '_field_' fieldTag(normFields(f))]) = fieldPct(f, :, a).';
    end
end
writetable(pctTable, outputFile, 'Sheet', 'Percentiles');

distTable = table(fieldAxes(AA(:)).', normFields(FF(:)).', 'VariableNames', {'Axis', 'NormField'});
for s = 1:2
    Dy = 100 * sys(s).D(sys(s).inYield, :, :);           % nYield x nF x nAx
    Dy = reshape(permute(Dy, [2 3 1]), nF * nAx, []);   % (field, axis) x file
    distTable.(char(sys(s).name + "_nominal_pct")) = reshape(100 * sys(s).nominal.D(1, :, :), [], 1);
    distTable.(char(sys(s).name + "_mean_pct")) = mean(Dy, 2);
    distTable.(char(sys(s).name + "_std_pct")) = std(Dy, 0, 2);
    distTable.(char(sys(s).name + "_min_pct")) = min(Dy, [], 2);
    distTable.(char(sys(s).name + "_max_pct")) = max(Dy, [], 2);
    distTable.(char(sys(s).name + "_nominal_height_mm")) = reshape(sys(s).nominal.yH(1, :, :), [], 1);
end
writetable(distTable, outputFile, 'Sheet', 'Distortion');

if hasAniso
    writetable(table(normFields(:), anisoNom, anisoStats(:, 1), anisoStats(:, 2), anisoStats(:, 3), ...
        anisoStats(:, 4), anisoStats(:, 5), anisoStats(:, 6), 'VariableNames', {'NormField', 'Nominal_pct', ...
        'Mean_pct', 'Std_pct', 'Min_pct', 'Max_pct', 'P1_pct', 'P99_pct'}), outputFile, 'Sheet', 'Axis Asymmetry');
end

for s = 1:2
    fileTable = table(sys(s).fileNames, sys(s).EFL, sys(s).mtfMetric, sys(s).mtfRank, sys(s).inYield, ...
        'VariableNames', {'FileName', 'EFL_paraxial_abs_mm', 'MTF_criterion', 'MTF_rank', 'InYieldSet'});
    for a = 1:nAx
        for f = 1:nF
            fileTable.(['fEff_' char(fieldAxes(a)) '_' fieldTag(normFields(f)) '_mm']) = sys(s).fEff(:, f, a);
        end
    end
    for a = 1:nAx
        for f = 1:nF
            fileTable.(['Distortion_' char(fieldAxes(a)) '_' fieldTag(normFields(f)) '_pct']) = 100 * sys(s).D(:, f, a);
        end
    end
    writetable(fileTable, outputFile, 'Sheet', char(sys(s).name + " Files"));
end

if writePairsSheet
    if nPairs <= 1048575
        fprintf('Writing %d pairs to the "Pairs" sheet (this can take a few minutes)...\n', nPairs);
        [I, J] = ndgrid(1:n1, 1:n2);
        pairTable = table(names1(I(:)), names2(J(:)), E1(I(:)), E2(J(:)), ...
            'VariableNames', {char("File_" + sys(1).name), char("File_" + sys(2).name), ...
            char("EFL_" + sys(1).name + "_mm"), char("EFL_" + sys(2).name + "_mm")});
        for a = 1:nAx
            for f = 1:nF
                mf = M(:, :, f, a);
                pairTable.(['M_' char(fieldAxes(a)) '_' fieldTag(normFields(f))]) = mf(:);
            end
        end
        for a = 1:nAx
            for f = 2:nF
                df = dM(:, :, f, a);
                pairTable.(['dM_' char(fieldAxes(a)) '_' fieldTag(normFields(f)) '_pct']) = df(:);
            end
        end
        if hasAniso
            for f = 2:nF
                af = aniso(:, :, f);
                pairTable.(['Asym_' fieldTag(normFields(f)) '_pct']) = af(:);
            end
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
    "Pixels"; "Pixel scale k"; "c = k*p2/p1"; "Normalized fields"; "Field axes"; "Magnification formula"];
settingValues = strings(numel(settingNames), 1);
for s = 1:2
    o = (s - 1) * 7;
    settingValues(o + (1:7)) = [sys(s).name; string(sys(s).folder); string(sys(s).nominalName); ...
        string(sys(s).isAfocal); string(sys(s).mtfFreq); string(sys(s).threshold); fieldText(mtfFields{s})];
end
settingValues(15:end) = [string(filePrefix); string(yieldPct); string(mtfFreqFrac); ...
    strjoin(string(pixelPitch_um), " / "); strjoin(string(nPixels), " / "); string(pixelScale); ...
    string(magConst); strjoin(fieldNames, ", "); strjoin(fieldAxes, ", "); ...
    "M(h) = c * f1_eff(h) / f2(c * y1(h)),  f2(y) = EFL2 + a*y^2 + b*y^4,  f_eff from the on-axis chief ray"];
writetable(table(settingNames, settingValues, 'VariableNames', {'Parameter', 'Value'}), ...
    outputFile, 'Sheet', 'Settings');

%% Save raw data
save(fullfile(outFolder, 'Distortion_Magnification_Data.mat'), 'sys', 'names1', 'names2', 'M', 'dM', ...
    'normFields', 'fieldAxes', 'fieldStats', 'fieldPct', 'changeStats', 'mNomField', 'coef2', 'magConst', ...
    'yieldPct', 'mtfFreqFrac', 'pixelPitch_um', 'nPixels', 'pixelScale', '-v7.3');

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
axisList = strjoin(fieldAxes, ' ve ');
modeText = ["nesne sonsuzda: y = görüntüde REAY/REAX, açı = nesne uzayında REAB/REAC, REAA/REAC", ...
            "afokal: y = nesnede (OLED) REAY/REAX, açı = görüntü uzayında REAB/REAC, REAA/REAC"];

% --- Page: summary -------------------------------------------------------
summaryLines = {
    sprintf('Tarih:                 %s', char(datetime('now', 'Format', 'dd.MM.yyyy HH:mm')))
    sprintf('Alanlar (%s):        %s  (0 = on-axis), eksenler: %s', n1Name, strjoin(fieldNames, ', '), axisList)
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
        }]; %#ok<AGROW>
    for a = 1:nAx
        summaryLines{end + 1, 1} = sprintf('  Nominal dist. (%s):  %s', fieldAxes(a), ...
            strjoin(compose("%s: %+.3f%%", fieldNames(:), 100 * sys(s).nominal.D(1, :, a).').', ', ')); %#ok<SAGROW>
    end
    summaryLines{end + 1, 1} = ''; %#ok<SAGROW>
end
for a = 1:nAx
    summaryLines = [summaryLines; {
        sprintf('BÜYÜTME M – %s EKSENİ (tüm çiftler)', upper(fieldAxes(a)))
        '  Alan   Nominal M   Ort. M     Std M     Min M     Max M     ΔM ort. %'
        }]; %#ok<AGROW>
    for f = 1:nF
        summaryLines{end + 1, 1} = sprintf('  %-5s  %9.5f  %9.5f  %8.5f  %8.5f  %8.5f  %+8.3f', fieldNames(f), ...
            mNomField(f, a), fieldStats(f, 1, a), fieldStats(f, 2, a), fieldStats(f, 3, a), fieldStats(f, 4, a), ...
            changeStats(f, 1, a)); %#ok<SAGROW>
    end
    summaryLines{end + 1, 1} = ''; %#ok<SAGROW>
end
if hasAniso
    summaryLines = [summaryLines; {
        sprintf('ASİMETRİ M_%s / M_%s − 1 (%%)', fieldAxes(2), fieldAxes(1))
        '  Alan   Nominal    Ortalama   Std       Min        Max'
        }];
    for f = 2:nF
        summaryLines{end + 1, 1} = sprintf('  %-5s  %+8.4f  %+9.4f  %8.4f  %+9.4f  %+9.4f', fieldNames(f), ...
            anisoNom(f), anisoStats(f, 1), anisoStats(f, 2), anisoStats(f, 3), anisoStats(f, 4)); %#ok<SAGROW>
    end
end
summaryFig = textPage('DİSTORSİYONLU BÜYÜTME ANALİZ RAPORU', summaryLines, 'Courier New', 7.5);

% --- Page: assessment ----------------------------------------------------
assess = {'ON-AXIS VE OFF-AXIS BÜYÜTME'};
assess{end + 1, 1} = sprintf(['%%%g MTF yield kümelerindeki %d çiftin on-axis büyütmesi ortalama %.5f (std %.5f, ' ...
    '%.5f – %.5f).'], yieldPct, nPairs, fieldStats(1, 1, 1), fieldStats(1, 2, 1), fieldStats(1, 3, 1), fieldStats(1, 4, 1));
for a = 1:nAx
    assess{end + 1, 1} = sprintf(['•  %s ekseni, tam alan (%s): ortalama %.5f (std %.5f, %.5f – %.5f); on-axis''e göre ' ...
        'ortalama değişim %+.3f %%, çiftler arasında %+.3f %% ile %+.3f %% arası.'], fieldAxes(a), fieldNames(end), ...
        fieldStats(end, 1, a), fieldStats(end, 2, a), fieldStats(end, 3, a), fieldStats(end, 4, a), ...
        changeStats(end, 1, a), changeStats(end, 3, a), changeStats(end, 4, a)); %#ok<SAGROW>
end
assess = [assess; {''; 'NOMİNAL DİSTORSİYONUN ETKİSİ'}];
assess{end + 1, 1} = sprintf(['Tam alanda nominal distorsiyon (%s ekseni): %s %+.3f %%, %s %+.3f %%. Nominal çiftin ' ...
    'büyütmesi on-axis %.5f, tam alanda %.5f; distorsiyon tam alanda büyütmeyi %+.3f %% değiştiriyor.'], ...
    fieldAxes(1), n1Name, 100 * sys(1).nominal.D(1, end, 1), n2Name, 100 * sys(2).nominal.D(1, end, 1), ...
    mNomField(1, 1), mNomField(end, 1), 100 * dDistNom(end, 1));
if sign(sys(1).nominal.D(1, end, 1)) == sign(sys(2).nominal.D(1, end, 1)) && sys(1).nominal.D(1, end, 1) ~= 0
    assess{end + 1, 1} = ['İki sistemin distorsiyonu aynı işaretli; büyütmede birbirini kısmen götürüyor. ' ...
        'Off-axis büyütmeyi belirleyen iki distorsiyonun farkıdır.'];
else
    assess{end + 1, 1} = ['İki sistemin distorsiyonu zıt işaretli; büyütmedeki etkileri toplanıyor. ' ...
        'Distorsiyonların eşleştirilmesi (aynı işaret ve büyüklük) off-axis büyütme farkını azaltır.'];
end
assess = [assess; {''; 'TOLERANSLARIN ETKİSİ'}];
assess{end + 1, 1} = sprintf(['Ortalama ΔM (%s, tam alan) %+.3f %%: bu kısım nominal distorsiyon farkından gelir ve tüm ' ...
    'ürünlerde aynıdır; tasarımda veya görüntü işlemede düzeltilebilir. ΔM''nin çiftler arasındaki std''si %.3f %%: ' ...
    'toleransların distorsiyonu dosyadan dosyaya değiştirmesinden gelir ve ürüne özeldir.'], fieldAxes(1), ...
    changeStats(end, 1, 1), changeStats(end, 2, 1));
if hasAniso
    assess = [assess; {''; sprintf('%s – %s ASİMETRİSİ (TILT / DECENTER)', upper(fieldAxes(2)), upper(fieldAxes(1)))}];
    assess{end + 1, 1} = sprintf(['Nominal tasarım simetrik olduğu için nominal çiftte asimetri %+.4f %%. Toleranslı ' ...
        'çiftlerde tam alanda M_%s / M_%s − 1 ortalama %+.4f %%, std %.4f %%, %+.4f %% ile %+.4f %% arası (%%1 – %%99: ' ...
        '%+.4f %% … %+.4f %%). Bu fark tilt ve decenter toleranslarının görüntüyü bir yönde diğerinden farklı ' ...
        'büyütmesinden gelir.'], anisoNom(end), fieldAxes(2), fieldAxes(1), anisoStats(end, 1), anisoStats(end, 2), ...
        anisoStats(end, 3), anisoStats(end, 4), anisoStats(end, 5), anisoStats(end, 6));
end
assess = [assess; {
    ''
    'KONTROL'
    ['Okunan nominal distorsiyonu, nominal dosyalarda OpticStudio''nun Field Curvature / Distortion analizindeki ' ...
     'F-Tan(Theta) distorsiyonuyla karşılaştırın (özet sayfasında her alan ve eksen için verilmiştir).']
    }];
if any(extrapPct > 2)
    assess{end + 1, 1} = sprintf(['UYARI: %s''ın tam alanı OLED''de, %s''in örneklenen en büyük yüksekliğinin %%%.1f ' ...
        'üzerine çıkıyor; orada kolimatörün etkin EFL''si uyum eğrisinden dışarı doğru tahmin edildi.'], ...
        n1Name, n2Name, max(extrapPct));
end
if min(n1, n2) < 50
    assess = [assess; {''; 'UYARI – ÖRNEKLEM BÜYÜKLÜĞÜ'}];
    assess{end + 1, 1} = sprintf(['Yield kümelerinde %d ve %d dosya var; istatistiklerin güvenilirliği çift ' ...
        'sayısıyla değil dosya sayısıyla sınırlıdır. Sistem başına en az 50 – 100 MC dosyası önerilir.'], n1, n2);
end
assessFig = textPage('1. DEĞERLENDİRME', assess, 'Arial', 9.5);

% --- Page: method --------------------------------------------------------
method = {
    'NEDEN DİSTORSİYON?'
    ['Paraksiyel EFL yalnızca on-axis büyütmeyi verir. Distorsiyonlu bir sistemde gerçek görüntü yüksekliği ' ...
     'y = f · tan θ · (1 + D(θ)) olur; bu yüzden büyütme alan boyunca değişir: M(θ) ≈ M(0) × (1 + D_1 − D_2).']
    ''
    'ADIM 1 – ETKİN EFL (İKİ EKSEN)'
    sprintf(['Her MC dosyasında normalize alanlarda (%s) ve %s eksenlerinde birincil dalga boyunda gerçek ana ışın ' ...
     'izlenir (y ekseni: Hx = 0, Hy = h, REAY ve REAB / REAC; x ekseni: Hx = h, Hy = 0, REAX ve REAA / REAC). ' ...
     'Toleranslı dosyalarda tilt / decenter on-axis ana ışını da eksenden kaydırır (boresight); bu kayma büyütmeyi ' ...
     'değiştirmediği için etkin EFL on-axis ana ışına göre farkla hesaplanır: ' ...
     'f_eff(h) = |y(h) − y(0)| / |tan θ(h) − tan θ(0)|. Nesnesi sonsuzda olan dosyada y görüntü yüzeyinde, açı ' ...
     'nesne uzayında; afokal dosyada (nesne = OLED) y nesne yüzeyinde, açı görüntü uzayında okunur. On-axis''te ' ...
     'paraksiyel EFL (EFFL) kullanılır. Distorsiyon D(h) = f_eff(h) / EFL − 1''dir.'], strjoin(fieldNames, ', '), axisList)
    ''
    'ADIM 2 – YIELD KÜMELERİ'
    sprintf(['Her sistemde MTF''e göre en iyi %%%g dosya seçilir. Alan ışınları izlenemeyen dosyalar geçersiz ' ...
     'sayılır.'], yieldPct)
    ''
    'ADIM 3 – ALANA GÖRE BÜYÜTME'
    sprintf(['Her eksen için ayrı: %s''ın h alanında dedektördeki yükseklik y_1(h) OLED''e c = %.6f ölçeğiyle basılır: ' ...
     'y_OLED = c × y_1(h). %s''in aynı eksendeki alan noktalarından her dosya için f_2(y) = EFL_2 + a·y² + b·y⁴ ' ...
     'eğrisi uydurulur. Büyütme M(h) = c × f_1,eff(h) / f_2(y_OLED)''dir.'], n1Name, magConst, n2Name)
    ''
    'ADIM 4 – ÇAPRAZLAMA VE ASİMETRİ'
    sprintf(['%s yield kümesindeki her dosya %s yield kümesindeki her dosyayla eşleştirilir (%d × %d = %d çift). Her ' ...
     'çift, alan ve eksen için büyütme; on-axis''e göre değişim ΔM(h) = M(h) / M(0) − 1 ve aynı çiftin iki eksen ' ...
     'arasındaki farkı (asimetri) M_x / M_y − 1 hesaplanır.'], n1Name, n2Name, n1, n2, nPairs)
    ''
    'VARSAYIMLAR'
    ['•  Her eksende alanın yalnız pozitif tarafı örneklenir. Ana ışın birincil dalga boyunda izlenir; lateral ' ...
     'renk dahil değildir. Görüntü işlemede distorsiyon düzeltmesi yapılmadığı ve MWIR''ın x / y eksenlerinin ' ...
     'OLED''in x / y eksenlerine karşılık geldiği varsayılır.']
    };
methodFig = textPage('2. YÖNTEM', method, 'Arial', 9.5);

% --- Figure pages --------------------------------------------------------
plotFigs = figurePage(fig2, '3. ON-AXIS VE OFF-AXIS BÜYÜTME', {
    'NASIL OKUNUR'
    sprintf(['•  Yatay eksen %s''ın normalize alanıdır (0 = on-axis). Her renk bir eksendir (%s). Kalın çizgi tüm ' ...
     'çiftlerin ortalaması, renkli bant çiftlerin %%1 – %%99''u, noktalı çizgiler en küçük ve en büyük değer, ' ...
     'siyah kesikli çizgi nominal çifttir.'], n1Name, axisList)
    'NASIL YORUMLANIR'
    ['•  Ortalama çizginin on-axis''ten kenara eğimi distorsiyon farkının etkisidir. İki eksenin ortalamaları ' ...
     'ayrışıyorsa toleranslar ortalamada bir yönü diğerinden farklı etkiliyordur.']
    }, false);
plotFigs(end + 1, 1) = figurePage(fig3, '4. ON-AXIS''E GÖRE BÜYÜTME DEĞİŞİMİ', {
    'NASIL OKUNUR'
    ['•  Her çiftin off-axis büyütmesinin kendi on-axis büyütmesine göre yüzde değişimidir: ' ...
     'ΔM(h) = M(h) / M(0) − 1. EFL saçılımını ayıklar ve yalnız distorsiyonun etkisini gösterir.']
    'NASIL YORUMLANIR'
    ['•  Ortalama çizgi nominal distorsiyon farkını, bant genişliği toleransların distorsiyona etkisini verir. ' ...
     'Örneğin tam alanda ortalama −%0.8 ise on-axis büyütmesi 1.000 olan ürünün kenarda büyütmesi ≈ 0.992''dir.']
    }, false);
if hasAniso
    plotFigs(end + 1, 1) = figurePage(fig5, sprintf('5. %s – %s ASİMETRİSİ', upper(fieldAxes(2)), upper(fieldAxes(1))), {
        'NASIL OKUNUR'
        sprintf(['•  Aynı çiftin aynı alanda %s eksenindeki büyütmesinin %s eksenindekine oranı − 1 (%%). Nominal ' ...
         'simetrik tasarımda sıfırdır; on-axis''te tanım gereği sıfırdır.'], fieldAxes(2), fieldAxes(1))
        'NASIL YORUMLANIR'
        ['•  Bant genişliği tilt / decenter toleranslarının yarattığı asimetridir: görüntü bir yönde diğerinden ' ...
         'farklı büyütülür (örneğin bir kare hafif dikdörtgen görünür). Bu ürüne özeldir.']
        ['•  Ortalama sıfırdan belirgin uzaksa toleranslar sistematik olarak bir yönü tercih ediyordur ' ...
         '(örneğin yalnız bir eksende tanımlı tilt toleransları).']
        }, false);
end
plotFigs(end + 1, 1) = figurePage(fig4, sprintf('%d. ALANLARIN BÜYÜTME DAĞILIMI', 4 + hasAniso + 1), {
    'NASIL OKUNUR'
    '•  Her panel bir eksen, her renkli basamak çizgisi bir alanın büyütme dağılımıdır (tüm çiftler).'
    'NASIL YORUMLANIR'
    ['•  Off-axis dağılımları on-axis''e göre kaymışsa sebep nominal distorsiyon farkıdır; genişlemişse sebep ' ...
     'toleransların distorsiyona etkisidir.']
    }, false);
plotFigs(end + 1, 1) = figurePage(fig1, sprintf('%d. ALANA GÖRE DİSTORSİYON', 4 + hasAniso + 2), {
    'NASIL OKUNUR'
    ['•  Her panel bir sistemdir; her renk bir eksen. Kalın çizgi yield kümesinin ortalaması, renkli bant en küçük ' ...
     've en büyük değer, siyah kesikli çizgi nominal tasarımdır (her eksen kendi işaretiyle).']
    'NASIL YORUMLANIR'
    ['•  İki sistemin tam alandaki nominal distorsiyonları aynı işaretli ve yakın büyüklükteyse büyütmeye ' ...
     'etkileri birbirini götürür; zıt işaretliyse toplanır. Bantların genişliği toleransların etkisidir.']
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
    normFields, fieldAxes, mtfFreqFrac, mtfFields)
% Paraxial |EFL|, field data (effective EFL, height and tan of the angle
% measured from the on-axis chief ray, at each normalized field and axis) and
% MTF criterion of every MC file of one system, and the field data of the
% nominal file. fEff, yH, tH: nFiles x nF x nAxes.
nFiles = numel(S.fileNames);
nF = numel(normFields);
nAx = numel(fieldAxes);
EFL = nan(nFiles, 1);
fEff = nan(nFiles, nF, nAx);
yH = nan(nFiles, nF, nAx);
tH = nan(nFiles, nF, nAx);
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
nominal.fEff = nan(1, nF, nAx);
nominal.yH = nan(1, nF, nAx);
nominal.tH = nan(1, nF, nAx);
for a = 1:nAx
    [nominal.fEff(1, :, a), nominal.yH(1, :, a), nominal.tH(1, :, a)] = ...
        fieldEfl(TheSystem, normFields, isAfocal, nominal.EFL, fieldAxes(a));
end
if any(~isfinite(nominal.fEff(:)))
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
        for a = 1:nAx
            [fEff(k, :, a), yH(k, :, a), tH(k, :, a)] = fieldEfl(TheSystem, normFields, isAfocal, EFL(k), fieldAxes(a));
        end
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

function [fEff, dH, dT] = fieldEfl(TheSystem, normFields, isAfocal, efl0, axisName)
% Real chief ray (Px = Py = 0, primary wavelength) at each normalized field h
% along one axis: "y" -> (Hx, Hy) = (0, h), height REAY, direction REAB / REAC;
% "x" -> (Hx, Hy) = (h, 0), height REAX, direction REAA / REAC.
%   Not afocal (object at infinity): height on the image surface, angle in
%     object space (direction cosines at surface 0).
%   Afocal (object = OLED): height on the object surface, angle in image space.
% Height and tan(angle) are measured from the on-axis chief ray (h = 0), which
% removes the boresight shift of tilted / decentered files:
%   dH = |y(h) - y(0)|, dT = |tan(h) - tan(0)|, fEff = dH / dT.
% normFields(1) must be 0; there the paraxial EFL efl0 is used.
imgSurf = TheSystem.LDE.NumberOfSurfaces - 1;
wave = primaryWave(TheSystem);
nF = numel(normFields);

if isAfocal
    hSurf = 0;
    aSurf = imgSurf;
else
    hSurf = imgSurf;
    aSurf = 0;
end
isX = strcmpi(axisName, "x");

hv = nan(1, nF);
tv = nan(1, nF);
for f = 1:nF
    if isX
        hx = normFields(f);
        hy = 0;
        p = TheSystem.MFE.GetOperandValue(ZOSAPI.Editors.MFE.MeritOperandType.REAX, hSurf, wave, hx, hy, 0, 0, 0, 0);
        d = TheSystem.MFE.GetOperandValue(ZOSAPI.Editors.MFE.MeritOperandType.REAA, aSurf, wave, hx, hy, 0, 0, 0, 0);
    else
        hx = 0;
        hy = normFields(f);
        p = TheSystem.MFE.GetOperandValue(ZOSAPI.Editors.MFE.MeritOperandType.REAY, hSurf, wave, hx, hy, 0, 0, 0, 0);
        d = TheSystem.MFE.GetOperandValue(ZOSAPI.Editors.MFE.MeritOperandType.REAB, aSurf, wave, hx, hy, 0, 0, 0, 0);
    end
    n = TheSystem.MFE.GetOperandValue(ZOSAPI.Editors.MFE.MeritOperandType.REAC, aSurf, wave, hx, hy, 0, 0, 0, 0);
    hv(f) = p;
    tv(f) = d / n;
end

dH = abs(hv - hv(1));
dT = abs(tv - tv(1));
fEff = nan(1, nF);
fEff(1) = efl0;
ok = dH > 0 & dT > 0 & isfinite(dH) & isfinite(dT);
ok(1) = false;
fEff(ok) = dH(ok) ./ dT(ok);
dH(1) = 0;
dT(1) = 0;
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
