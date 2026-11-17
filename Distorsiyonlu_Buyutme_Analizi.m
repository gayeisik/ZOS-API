%% DİSTORSİYONLU BÜYÜTME ANALİZİ (ZOS-API, MATLAB)
% NE YAPAR
%   Birlikte çalışan iki sistemin (1: MWIR objektif + dedektör, 2: OLED +
%   kolimatör) toplam büyütmesini, distorsiyonu da dahil ederek, görüş
%   alanının birkaç noktasında (normFields, örn. 0, 0.5, 0.7, 1 normalize alan)
%   hesaplar. Her sistemde MTF'e göre %yieldPct yield kümesine giren dosyalar
%   seçilir; 1. sistemin yield kümesindeki her dosya 2. sistemin yield
%   kümesindeki her dosyayla eşleştirilir (n1 x n2 çift). Her çiftin on-axis
%   (0) ve off-axis alanlardaki büyütmesi, y ve x eksenlerinde ayrı ayrı
%   hesaplanır ve alan / eksen bazında büyütme değerleri / dağılımları verilir.
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
%   Her alan ve eksen için: nominal çiftin M'si; tüm çiftlerin M ortalaması,
%   std, min, max ve yüzdelikleri.
%
% AYARLAR (User settings)
%   sysName, pixelPitch_um, nPixels, pixelScale, filePrefix, yieldPct,
%   mtfFreqFrac, mtfFields, normFields, fieldAxes, magTarget (yalnız referans çizgisi),
%   writePairsSheet, nHistBins
%
% ÇIKTILAR (1. sistemin klasöründe "BÜYÜTME-DİSTORSİYON ANALİZİ - %90 YIELD")
%   Yalnızca büyütme M ve dağılımı verilir (distorsiyon arka planda kullanılır).
%   Magnification_Histograms_y.png / _x.png : her alanın büyütme histogramı
%                                  (ortalama, std, min, max, %1 – %99, nominal)
%   Magnification_by_Field.png   : özet; alana göre ortalama, %1 – %99, min – max
%   *.fig                        : grafiklerin düzenlenebilir MATLAB dosyaları
%   Magnification_by_Field_Results.xlsx : "Magnification" (alan / eksen başına
%                                  istatistik), sistem başına dosya tabloları,
%                                  "Pairs" (her çiftin her alandaki M'si), "Settings"
%   Magnification_by_Field_Report.pdf : sonuç tablosu, nasıl okunur, grafikler
%   Magnification_by_Field_Data.mat : ham veriler
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

%% Statistics of the magnification at every field and axis
pctLevels = [1 5 50 95 99];
nP = numel(pctLevels);
fieldStats = zeros(nF, 4, nAx);                          % M: mean std min max
fieldPct = zeros(nF, nP, nAx);                           % M percentiles
for a = 1:nAx
    for f = 1:nF
        mf = M(:, :, f, a);
        mf = mf(:);
        fieldStats(f, :, a) = [mean(mf), std(mf), min(mf), max(mf)];
        fieldPct(f, :, a) = quantileSorted(sort(mf), pctLevels / 100);
    end
end
iP1 = find(pctLevels == 1);
iP99 = find(pctLevels == 99);
fieldKind = repmat("off-axis", 1, nF);
fieldKind(normFields == 0) = "on-axis";

%% Command window summary
fprintf('\n===== MAGNIFICATION (%g %% MTF yield sets, %d x %d = %d pairs) =====\n', yieldPct, n1, n2, nPairs);
for a = 1:nAx
    fprintf('\n%s axis\n', fieldAxes(a));
    fprintf('%-6s %-9s %10s %10s %9s %10s %10s %10s %10s\n', 'Field', '', 'Nominal', 'Mean', 'Std', ...
        'Min', 'Max', '1 %', '99 %');
    for f = 1:nF
        fprintf('%-6s %-9s %10.5f %10.5f %9.5f %10.5f %10.5f %10.5f %10.5f\n', fieldNames(f), fieldKind(f), ...
            mNomField(f, a), fieldStats(f, 1, a), fieldStats(f, 2, a), fieldStats(f, 3, a), fieldStats(f, 4, a), ...
            fieldPct(f, iP1, a), fieldPct(f, iP99, a));
    end
end

%% Output folder: "BÜYÜTME-DİSTORSİYON ANALİZİ - %90 YIELD" in the system 1 folder
outName = ['B' char(220) 'Y' char(220) 'TME-D' char(304) 'STORS' char(304) 'YON ANAL' char(304) 'Z' char(304) ...
    sprintf(' - %%%g YIELD', yieldPct)];
outFolder = fullfile(sys(1).folder, outName);
if ~isfolder(outFolder)
    mkdir(outFolder);
end

colBlue = [0 0.45 0.74];
axisPalette = [0 0.45 0.74; 0.85 0.33 0.10];
axisColors = axisPalette(mod(0:nAx - 1, 2) + 1, :);
mEdges = linspace(min(M(:)), max(M(:)), nHistBins + 1);    % same bins in every panel
nCols = min(nF, 2);
nRows = ceil(nF / nCols);

%% Figures 1..nAx: magnification histogram of every field (one figure per axis)
histFigs = gobjects(nAx, 1);
for a = 1:nAx
    histFigs(a) = figure('Name', sprintf('Magnification %s axis', fieldAxes(a)), 'Color', 'w', ...
        'Position', [100 + 20 * a, 60, 1000, 380 * nRows]);
    for f = 1:nF
        subplot(nRows, nCols, f);
        mf = M(:, :, f, a);
        histogram(mf(:), mEdges, 'FaceColor', axisColors(a, :), 'EdgeColor', 'none');
        hold on;
        yl = ylim;
        patch([fieldPct(f, iP1, a) fieldPct(f, iP99, a) fieldPct(f, iP99, a) fieldPct(f, iP1, a)], ...
            [0 0 yl(2) yl(2)], [0.5 0.5 0.5], 'FaceAlpha', 0.10, 'EdgeColor', 'none');
        xline(fieldStats(f, 1, a), '-k', 'LineWidth', 1.8);
        xline(mNomField(f, a), '--k', 'LineWidth', 1.1);
        if magTarget >= mEdges(1) && magTarget <= mEdges(end)
            xline(magTarget, ':', 'Color', [0.8 0 0], 'LineWidth', 1.2);
        end
        hold off;
        xlim([mEdges(1) mEdges(end)]);
        ylim(yl);
        text(0.98, 0.96, sprintf(['Mean     %.5f\nStd      %.5f\nMin      %.5f\nMax      %.5f\n' ...
            '1 %% – 99 %%\n%.5f – %.5f\nNominal  %.5f'], fieldStats(f, 1, a), fieldStats(f, 2, a), ...
            fieldStats(f, 3, a), fieldStats(f, 4, a), fieldPct(f, iP1, a), fieldPct(f, iP99, a), mNomField(f, a)), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'FontName', 'Courier New', 'FontSize', 8, 'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7]);
        xlabel('Magnification M');
        ylabel('Number of pairs');
        title(sprintf('Field %s (%s), %s axis', fieldNames(f), fieldKind(f), fieldAxes(a)));
        grid on;
        box on;
    end
    saveFigure(histFigs(a), fullfile(outFolder, sprintf('Magnification_Histograms_%s.png', fieldAxes(a))));
end

%% Figure: magnification by field (overview)
figOverview = figure('Name', 'Magnification by Field', 'Color', 'w', 'Position', [140 140 1000 560]);
hold on;
dx = 0.012 * ((1:nAx) - (nAx + 1) / 2);              % small offset so the axes do not overlap
hLeg = gobjects(0);
for a = 1:nAx
    xa = normFields + dx(a);
    for f = 1:nF
        plot([xa(f) xa(f)], [fieldStats(f, 3, a) fieldStats(f, 4, a)], '-', 'Color', axisColors(a, :), 'LineWidth', 0.8);
        plot([xa(f) xa(f)], [fieldPct(f, iP1, a) fieldPct(f, iP99, a)], '-', 'Color', axisColors(a, :), 'LineWidth', 5);
    end
    hLeg(end + 1) = plot(xa, fieldStats(:, 1, a), 'o-', 'Color', axisColors(a, :), 'LineWidth', 1.5, ...
        'MarkerFaceColor', 'w', 'MarkerSize', 7, 'DisplayName', sprintf('%s axis: mean', fieldAxes(a))); %#ok<SAGROW>
    for f = 1:nF
        text(xa(f) + 0.012, fieldStats(f, 1, a), sprintf('%.4f', fieldStats(f, 1, a)), 'FontSize', 8, ...
            'Color', axisColors(a, :), 'VerticalAlignment', 'middle');
    end
end
hLeg(end + 1) = plot(normFields, mNomField(:, 1), 'sk--', 'LineWidth', 1, 'DisplayName', 'Nominal pair');
hLeg(end + 1) = plot(nan, nan, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 5, 'DisplayName', 'Thick bar: 1 % – 99 % of pairs');
hLeg(end + 1) = plot(nan, nan, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8, 'DisplayName', 'Thin bar: min – max');
yline(magTarget, ':', 'Color', [0.8 0 0], 'LineWidth', 1.2, 'HandleVisibility', 'off');
hold off;
xlim([-0.06 1.1]);
set(gca, 'XTick', normFields);
xlabel(sprintf('Normalized field of %s (0 = on-axis)', sys(1).name));
ylabel('Magnification M');
title(sprintf('Magnification on-axis and off-axis (%d pairs, %g%% MTF yield sets)', nPairs, yieldPct));
legend(hLeg, 'Location', 'best');
grid on;
box on;
saveFigure(figOverview, fullfile(outFolder, 'Magnification_by_Field.png'));

%% Export Excel
outputFile = fullfile(outFolder, 'Magnification_by_Field_Results.xlsx');
if isfile(outputFile)
    delete(outputFile);
end
fieldTag = @(h) strrep(sprintf('%.2f', h), '.', 'p');  % 0.70 -> "0p70" for column names

[FF, AA] = ndgrid(1:nF, 1:nAx);
statT = table(fieldAxes(AA(:)).', normFields(FF(:)).', fieldKind(FF(:)).', 'VariableNames', {'Axis', 'NormField', 'Type'});
mS = reshape(permute(fieldStats, [1 3 2]), [], 4);
pS = reshape(permute(fieldPct, [1 3 2]), [], nP);
statT.Nominal_M = mNomField(:);
statT.Mean_M = mS(:, 1);
statT.Std_M = mS(:, 2);
statT.Min_M = mS(:, 3);
statT.Max_M = mS(:, 4);
for k = 1:nP
    statT.(sprintf('P%g_M', pctLevels(k))) = pS(:, k);
end
writetable(statT, outputFile, 'Sheet', 'Magnification');

for s = 1:2
    fileTable = table(sys(s).fileNames, sys(s).EFL, sys(s).mtfMetric, sys(s).mtfRank, sys(s).inYield, ...
        'VariableNames', {'FileName', 'EFL_paraxial_abs_mm', 'MTF_criterion', 'MTF_rank', 'InYieldSet'});
    writetable(fileTable, outputFile, 'Sheet', char(sys(s).name + " Files"));
end

if writePairsSheet
    if nPairs <= 1048575
        fprintf('Writing %d pairs to the "Pairs" sheet (this can take a few minutes)...\n', nPairs);
        [I, J] = ndgrid(1:n1, 1:n2);
        pairTable = table(names1(I(:)), names2(J(:)), ...
            'VariableNames', {char("File_" + sys(1).name), char("File_" + sys(2).name)});
        for a = 1:nAx
            for f = 1:nF
                mf = M(:, :, f, a);
                pairTable.(['M_' char(fieldAxes(a)) '_' fieldTag(normFields(f))]) = mf(:);
            end
        end
        writetable(pairTable, outputFile, 'Sheet', 'Pairs');
    else
        warning('%d pairs exceed the Excel row limit; the "Pairs" sheet is not written (see the .mat file).', nPairs);
    end
end

settingNames = ["System 1 folder"; "System 1 nominal file"; "System 1 valid / yield files"; ...
    "System 2 folder"; "System 2 nominal file"; "System 2 valid / yield files"; ...
    "Yield level (%)"; "c = k*p2/p1"; "Normalized fields"; "Field axes"];
settingValues = [string(sys(1).folder); string(sys(1).nominalName); ...
    sprintf('%d / %d', sum(sys(1).valid), n1); string(sys(2).folder); string(sys(2).nominalName); ...
    sprintf('%d / %d', sum(sys(2).valid), n2); string(yieldPct); string(magConst); ...
    strjoin(fieldNames, ", "); strjoin(fieldAxes, ", ")];
writetable(table(settingNames, settingValues, 'VariableNames', {'Parameter', 'Value'}), ...
    outputFile, 'Sheet', 'Settings');

%% Save raw data
save(fullfile(outFolder, 'Magnification_by_Field_Data.mat'), 'names1', 'names2', 'M', 'normFields', ...
    'fieldAxes', 'fieldStats', 'fieldPct', 'pctLevels', 'mNomField', 'magConst', 'yieldPct', '-v7.3');

%% PDF report (A4)
pdfFile = fullfile(outFolder, 'Magnification_by_Field_Report.pdf');
oldPages = dir(fullfile(outFolder, 'Magnification_by_Field_Report_page*.pdf'));
for i = 1:numel(oldPages)
    delete(fullfile(outFolder, oldPages(i).name));
end
if isfile(pdfFile)
    delete(pdfFile);
end
if isfile(pdfFile)   % still open in a PDF viewer
    pdfFile = fullfile(outFolder, ['Magnification_by_Field_Report_' ...
        char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.pdf']);
    fprintf('\nThe old report is open in another program; new report: %s\n', pdfFile);
end

n1Name = char(sys(1).name);
n2Name = char(sys(2).name);

% --- Page: results -------------------------------------------------------
summaryLines = {
    sprintf('Tarih:        %s', char(datetime('now', 'Format', 'dd.MM.yyyy HH:mm')))
    sprintf('%s:         %s  (%d / %d dosya %%%g MTF yield kümesinde)', n1Name, sys(1).nominalName, n1, ...
        numel(sys(1).fileNames), yieldPct)
    sprintf('%s:         %s  (%d / %d dosya %%%g MTF yield kümesinde)', n2Name, sys(2).nominalName, n2, ...
        numel(sys(2).fileNames), yieldPct)
    sprintf('Çift sayısı:  %d × %d = %d  (her %s her %s ile)', n1, n2, nPairs, n1Name, n2Name)
    ''
    };
for a = 1:nAx
    summaryLines = [summaryLines; {
        sprintf('BÜYÜTME M – %s EKSENİ', upper(fieldAxes(a)))
        '  Alan              Nominal   Ortalama   Std       Min       Max       %1        %99'
        }]; %#ok<AGROW>
    for f = 1:nF
        summaryLines{end + 1, 1} = sprintf('  %-5s %-10s  %8.5f  %8.5f  %7.5f  %8.5f  %8.5f  %8.5f  %8.5f', ...
            fieldNames(f), char(fieldKind(f)), mNomField(f, a), fieldStats(f, 1, a), fieldStats(f, 2, a), ...
            fieldStats(f, 3, a), fieldStats(f, 4, a), fieldPct(f, iP1, a), fieldPct(f, iP99, a)); %#ok<SAGROW>
    end
    summaryLines{end + 1, 1} = ''; %#ok<SAGROW>
end
summaryLines = [summaryLines; {
    'TABLO NASIL OKUNUR'
    '  Alan      : 0 = görüntünün ortası (on-axis), 1 = görüntünün kenarı.'
    '  Nominal   : nominal MWIR + nominal kolimatör çiftinin büyütmesi.'
    '  Ortalama  : tüm çiftlerin büyütmesinin ortalaması.'
    '  Std       : çiftler arasındaki saçılım (standart sapma).'
    '  Min / Max : en küçük ve en büyük büyütme.'
    '  %1 / %99  : çiftlerin %98''i bu iki değerin arasındadır.'
    }];
summaryFig = textPage('ALANA GÖRE BÜYÜTME – SONUÇLAR', summaryLines, 'Courier New', 8);

% --- Page: how to read ---------------------------------------------------
readText = {
    'NE HESAPLANDI?'
    sprintf(['%%%g MTF yield kümelerindeki her %s, her %s ile eşleştirildi (%d çift). Her çiftin büyütmesi ' ...
     'görüntünün ortasında (alan 0, on-axis) ve kenara doğru (alan %s, off-axis) ayrı ayrı hesaplandı. ' ...
     'Kenardaki büyütme, lenslerin distorsiyonunun etkisini de içerir; ortadaki büyütme EFL oranından gelir.'], ...
     yieldPct, n1Name, n2Name, nPairs, strjoin(fieldNames(2:end), ', '))
    sprintf('Eksenler: %s (görüntünün dikey ve yatay yönü).', strjoin(fieldAxes, ', '))
    ''
    'HİSTOGRAMLAR NASIL OKUNUR?'
    ['•  Her panel bir alandır. Çubuklar, büyütmesi o aralıkta olan çift sayısını gösterir. Bütün panellerde yatay ' ...
     'eksen aynıdır; böylece paneller doğrudan karşılaştırılabilir.']
    '•  Siyah düz çizgi: ortalama. Siyah kesikli çizgi: nominal çift. Kırmızı noktalı çizgi: M = 1.'
    '•  Gri gölge: çiftlerin %98''inin bulunduğu aralık (%1 – %99).'
    '•  Sağ üstteki kutuda sayısal değerler vardır.'
    ''
    'ÖZET GRAFİK NASIL OKUNUR?'
    ['•  Yatay eksen alan (0 = orta, 1 = kenar). Her alan için daire ortalamayı, kalın çubuk çiftlerin %98''ini, ' ...
     'ince çubuk en küçük – en büyük değeri gösterir. Her renk bir eksendir.']
    ''
    'NELERE BAKMALI?'
    ['•  Ortalamanın alanla değişimi: ortadan kenara ortalama büyütme ne kadar değişiyor? Bu değişim bütün ' ...
     'ürünlerde aynıdır (tasarımdan gelir).']
    ['•  Dağılımın genişliği (std, %1 – %99): üründen ürüne ne kadar fark var? Kenarda genişliyorsa toleranslar ' ...
     'kenardaki büyütmeyi ortadakinden daha çok oynatıyor demektir.']
    ['•  İki eksenin karşılaştırması: aynı alanda y ve x büyütmeleri farklıysa görüntü bir yönde diğerinden ' ...
     'farklı büyütülüyor demektir (tilt / decenter toleranslarından).']
    };
if min(n1, n2) < 50
    readText = [readText; {''; sprintf(['UYARI: Yield kümelerinde %d ve %d dosya var; istatistikler dosya sayısıyla ' ...
        'sınırlıdır. Sistem başına en az 50 – 100 MC dosyası önerilir.'], n1, n2)}];
end
readFig = textPage('SONUÇLAR NASIL OKUNUR?', readText, 'Arial', 10);

% --- Figure pages --------------------------------------------------------
plotFigs = figurePage(figOverview, 'ÖZET: ALANA GÖRE BÜYÜTME', {
    '•  Daire: ortalama (değeri yanında yazılı). Kalın çubuk: çiftlerin %98''i. İnce çubuk: en küçük – en büyük.'
    '•  Siyah kesikli çizgi: nominal çift. Kırmızı noktalı çizgi: M = 1.'
    }, false);
for a = 1:nAx
    plotFigs(end + 1, 1) = figurePage(histFigs(a), sprintf('BÜYÜTME DAĞILIMI – %s EKSENİ', upper(fieldAxes(a))), {
        '•  Her panel bir alan. Çubuklar: o büyütmedeki çift sayısı. Siyah düz çizgi: ortalama.'
        '•  Siyah kesikli: nominal çift. Kırmızı noktalı: M = 1. Gri gölge: çiftlerin %98''i (%1 – %99).'
        }, false); %#ok<AGROW>
end

% --- Assemble pages ------------------------------------------------------
allFigs = [summaryFig; readFig; plotFigs];
nPages = numel(allFigs);
headerText = sprintf('Alana Göre Büyütme  |  %s: %s  |  %s: %s', n1Name, sys(1).nominalName, n2Name, sys(2).nominalName);
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
        print(allFigs(p), fullfile(outFolder, sprintf('Magnification_by_Field_Report_page%d.pdf', p)), '-dpdf');
    end
    fprintf('Report saved as separate pages.\n');
end
close(allFigs);

%% Text summary
summaryFile = fullfile(outFolder, 'Analysis_Summary.txt');
fid = fopen(summaryFile, 'w', 'n', 'UTF-8');
if fid > 0
    fprintf(fid, '%s\r\n', 'ALANA GÖRE BÜYÜTME – SONUÇLAR', '');
    fprintf(fid, '%s\r\n', summaryLines{:});
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
% A4 page: title, the axes of srcFig (copied with their legends, in the same
% row / column layout as in srcFig) and an explanation text below them.
fig = newA4Page();
annotation(fig, 'textbox', [0.08 0.89 0.84 0.05], 'String', titleText, ...
    'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'Interpreter', 'none', 'VerticalAlignment', 'top');

srcAx = findobj(srcFig, 'Type', 'axes');
xPos = round(arrayfun(@(a) a.Position(1), srcAx), 2);
yPos = round(arrayfun(@(a) a.Position(2), srcAx), 2);
rowVals = unique(yPos, 'sorted');
rowVals = flipud(rowVals(:));                 % top row first
colVals = unique(xPos, 'sorted');
nRows = numel(rowVals);
nCols = numel(colVals);

% Plot area: up to y = 0.86 of the page; the text below it (less room for a short text)
if numel(captionLines) <= 3
    areaBottom = 0.17;
    captionBox = [0.08 0.04 0.84 0.11];
else
    areaBottom = 0.34;
    captionBox = [0.08 0.05 0.84 0.25];
end
areaTop = 0.86;
rowH = (areaTop - areaBottom) / nRows;
colW = 0.84 / nCols;
for i = 1:numel(srcAx)
    a = srcAx(i);
    r = find(rowVals == yPos(i), 1);
    c = find(colVals == xPos(i), 1);
    if ~isempty(a.Legend)
        copied = copyobj([a.Legend, a], fig);
    else
        copied = copyobj(a, fig);
    end
    newAx = findobj(copied, 'Type', 'axes');
    newAx.Units = 'normalized';
    newAx.Position = [0.08 + (c - 1) * colW + 0.06, areaTop - r * rowH + 0.05, colW - 0.08, rowH - 0.08];
    newAx.FontSize = 8;
    newAx.Title.FontSize = 9;
    if squareAxes
        pbaspect(newAx, [1 1 1]);
    end
end

annotation(fig, 'textbox', captionBox, 'String', captionLines, ...
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
