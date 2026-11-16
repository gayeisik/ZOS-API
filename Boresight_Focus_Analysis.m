%% BORESIGHT ANALİZİ – SONSUZ ODAK → YAKIN ODAK (ZOS-API, MATLAB)
% NE YAPAR
%   Sistem sonsuz odaktan yakın odağa geçtiğinde (yalnızca kompanzatör lens
%   hareket eder) on-axis boresight'ın dedektörde kaç PİKSEL kaydığını her
%   Monte Carlo (MC) dosyası için hesaplar. Kompanzatörün yakın odaktaki
%   konumlandırma (encoder) hatası ±δ da eklenebilir.
%
% GİRDİLER (tek klasör)
%   - Nominal yakın odak dosyası: adı "MC_" ile başlamayan .zmx. En az iki
%     konfigürasyon (biri sonsuz odak, biri yakın odak). Yakın odak bir
%     paraxial yüzeyle simüle edilir; konfigürasyonlar arasında yalnızca
%     kompanzatör kalınlığı değişir.
%   - MC dosyaları: adı filePrefix ile başlayan .zmx. Afokal mod, sonsuz odak,
%     tek konfigürasyon (tolerans analizinin kaydettiği haliyle).
%
% OTOMATİK TESPİT (yorum/comment alanı kullanılmaz)
%   Nominal yakın odak dosyasında, MCE'den:
%     - Yakın odak paraxial yüzeyi: IGNR satırının yüzeyi. Ignore olduğu (1)
%       konfigürasyon sonsuz odak, ignore olmadığı (0) konfigürasyon yakın odaktır.
%     - Kompanzatör: değeri sonsuz ve yakın odak arasında değişen THIC satırının
%       yüzeyi (sonsuz odakta 0, yakın odakta 0 olmayan).
%   Her MC dosyasında (yüzey numaraları nominalden farklı olabilir, çünkü
%   toleranslama coordinate break ekler; bu yüzden her dosyada ayrı bulunur):
%     - Kompanzatör: kalınlığı variable olan tek yüzey.
%     - Yakın odak paraxial yüzeyi görüntü yüzeyinden hemen önce eklenir.
%     - |EFL|: dosyanın kendi EFL'i (afokal modda EFFL).
%   Tespit edilemezse ya da birden fazla aday varsa kullanıcıya liste ile
%   sorulur. Ayarlarda değer verilirse tespit yapılmaz, o değer kullanılır.
%
% NASIL ÇALIŞIR
%   1. Nominal dosyadan:
%        Δkomp  = T_komp(yakın config) − T_komp(sonsuz config)
%        f_yakın = yakın odak paraxial yüzeyinin focal length'i (yakın config)
%   2. Her MC dosyasında:
%      a. Sonsuz odak: dosya olduğu gibi. On-axis ana ışının (Hx=Hy=0,
%         Px=Py=0, ana dalga boyu) görüntü uzayındaki yönü okunur:
%         tan θx = L/N, tan θy = M/N (REAA, REAB, REAC).
%      b. Yakın odak: görüntü yüzeyinden hemen önce paraxial yüzey eklenir
%         (thickness 0, focal length f_yakın) ve kompanzatör kalınlığı
%         T0 + Δkomp (± δ) yapılır. Yön tekrar okunur.
%      c. Boresight hatası [piksel] = (tan θ_yakın − tan θ_sonsuz) × |EFL| / piksel
%         Bu, dedektör tarafında f = |EFL| paraxial lens ile yapılan hesabın
%         aynısıdır (y = f·tan θ; ışının yüksekliğinden bağımsız).
%      d. Eklenen yüzey silinir, kompanzatör geri alınır. Diske kaydedilmez.
%   3. Encoder hatası δ için +δ ve −δ ikisi de hesaplanır; boresight kaymasını
%      büyüten taraf (worst case) raporlanır. δ = 0 her zaman hesaplanır.
%   Hiçbir adımda optimizasyon yapılmaz.
%
% KONTROLLER (ilk geçerli MC dosyasında, Command Window)
%   - Kolimasyon: yakın odak durumunda çıkış kolime olmalıdır. Paraxial
%     focal length'in işareti ters çevrildiğinde çıkış belirgin şekilde daha
%     iyi kolime oluyorsa uyarı verilir (işaret ya da Δkomp yönü şüpheli).
%   - Kompanzatör değişikliğinin sisteme gerçekten uygulandığı kontrol edilir.
%
% ÇIKTILAR (seçilen klasörde "BORESIGHT ANALİZİ" klasörüne)
%   Boresight_Results.xlsx : Ozet (seviye başına ortalama, std, yield, max),
%                            Dosyalar (dosya başına dx, dy, büyüklük, uygulanan
%                            işaret), Ayarlar (okunan/seçilen tüm değerler)
%   Boresight_CDF.png      : boresight hatasının kümülatif dağılımı
%   Boresight_Scatter.png  : dx–dy dağılımı (her seviye için)
%   Boresight_RawData.mat  : ham sonuçlar

clc;
clear;
close all;

%% Kullanıcı ayarları
filePrefix    = 'MC_T';          % yalnızca bu önekle başlayan dosyalar kullanılır
pixelPitch_um = 15;              % dedektör piksel boyutu (um)
compErrors_um = [8 16];          % yakın odakta kompanzatör (encoder) hatası, +/- um
                                 % (0 um = hatasız durum her zaman hesaplanır)
yieldLevels   = [50 80 90 95];   % yield seviyeleri (%): dosyaların %Y'si bu değerin altında
nFilesMax     = Inf;             % hızlı deneme için yalnızca ilk N MC dosyası

% Konfigürasyon ve yüzeyler: [] = otomatik (açıklama yukarıda), gerekirse sorulur
infConfig      = [];   % nominal dosyada sonsuz odak konfigürasyonu
nearConfig     = [];   % nominal dosyada yakın odak konfigürasyonu
nomCompSurface = [];   % nominal dosyada kompanzatör yüzeyi
nomNearSurface = [];   % nominal dosyada yakın odak paraxial yüzeyi
mcCompSurface  = [];   % MC dosyalarında kompanzatör yüzeyi ([] = her dosyada variable kalınlık)
mcInsertAfter  = [];   % MC'de paraxial BU yüzeyden sonra eklenir ([] = görüntüden hemen önce)
efl_mm         = [];   % |EFL| (mm) ([] = her MC dosyasının kendi EFL'i)

%% Klasör seçimi
folderPath = uigetdir(pwd, 'Klasörü seçin (nominal yakın odak dosyası + MC dosyaları)');
if isequal(folderPath, 0)
    return;
end

files = [dir(fullfile(folderPath, [filePrefix '*.zmx'])); dir(fullfile(folderPath, [filePrefix '*.ZMX']))];
files = files(startsWith(string({files.name}), filePrefix));
[~, uniqueIdx] = unique(lower(string({files.name})), 'stable');
files = files(uniqueIdx);
files = files(1:min(end, nFilesMax));
nFiles = numel(files);
if nFiles == 0
    error('"%s" ile başlayan .zmx dosyası bulunamadı.', filePrefix);
end

allZmx = [dir(fullfile(folderPath, '*.zmx')); dir(fullfile(folderPath, '*.ZMX'))];
[~, uniqueIdx] = unique(lower(string({allZmx.name})), 'stable');
allZmx = allZmx(uniqueIdx);
nominalFiles = allZmx(~startsWith(string({allZmx.name}), 'MC_', 'IgnoreCase', true));

if isempty(nominalFiles)
    error('Klasörde nominal yakın odak dosyası (adı "MC_" ile başlamayan .zmx) yok.');
elseif isscalar(nominalFiles)
    nominalName = nominalFiles(1).name;
else
    [sel, ok] = listdlg('ListString', {nominalFiles.name}, 'SelectionMode', 'single', ...
        'PromptString', 'Nominal YAKIN ODAK dosyasını seçin:', 'ListSize', [400 200]);
    if ~ok
        return;
    end
    nominalName = nominalFiles(sel).name;
end

%% ZOS-API bağlantısı
% InitConnection ve CleanupConnection resmi Zemax MATLAB örneklerinden alınmıştır.
TheApplication = InitConnection();
if isempty(TheApplication)
    error('OpticStudio''ya bağlanılamadı.');
end

levelErrors_um = [0, compErrors_um];
nLevels = numel(levelErrors_um);
fileNames = string({files.name})';
status = strings(nFiles, 1);
dxPix = nan(nFiles, nLevels);        % boresight kayması x (piksel)
dyPix = nan(nFiles, nLevels);        % boresight kayması y (piksel)
appliedSign = zeros(nFiles, nLevels);% +1 / -1: uygulanan encoder hatası işareti (0: hatasız)
compT0 = nan(nFiles, 1);             % MC dosyasındaki kompanzatör kalınlığı
tanInf = nan(nFiles, 2);             % sonsuz odakta ana ışın yönü (tan θx, tan θy)
compSurfUsed = nan(nFiles, 1);       % MC dosyasında kullanılan kompanzatör yüzeyi
eflUsed = nan(nFiles, 1);            % MC dosyasında kullanılan |EFL| (mm)

try
    TheSystem = TheApplication.PrimarySystem;

    %% 1. Nominal yakın odak dosyası: Δkomp ve yakın odak focal length
    loadFile(TheSystem, fullfile(folderPath, nominalName));
    nomUnits = unitsPerMm(TheSystem);
    mce = TheSystem.MCE;
    nConf = mce.NumberOfConfigurations;
    if nConf < 2
        error('Nominal dosyada %d konfigürasyon var; en az 2 gerekli (sonsuz ve yakın odak).', nConf);
    end

    confItems = cellstr(compose('Config %d', (1:nConf)'));

    % MCE'den otomatik tespit: IGNR -> yakın odak paraxial yüzeyi ve konfigürasyonlar
    [thicOps, ignrOps] = readMceOperands(TheSystem);
    if isempty(nomNearSurface)
        isSwitch = false(numel(ignrOps), 1);
        for i = 1:numel(ignrOps)
            v = ignrOps(i).val;
            isSwitch(i) = any(v ~= 0) && any(v == 0);   % bazı config'lerde ignore, bazılarında değil
        end
        cand = ignrOps(isSwitch);
        if numel(cand) > 1
            names = arrayfun(@(r) sprintf('Yüzey %d: IGNR %s', r.surf, mat2str(r.val)), cand, 'UniformOutput', false);
            cand = cand(askList(names, 'Birden fazla IGNR satırı var. Yakın odak paraxial yüzeyi hangisi?', 1, 'Nominal dosya'));
        end
        if isscalar(cand)
            nomNearSurface = cand.surf;
            if isempty(infConfig)
                infConfig = find(cand.val ~= 0, 1);   % ignore -> sonsuz odak
            end
            if isempty(nearConfig)
                nearList = find(cand.val == 0);       % ignore değil -> yakın odak
                if isscalar(nearList)
                    nearConfig = nearList;
                else
                    nearConfig = nearList(askList(confItems(nearList), ...
                        'Birden fazla yakın odak konfigürasyonu var. Hangisi?', 1, 'Nominal dosya'));
                end
            end
        end
    end

    if isempty(infConfig)
        infConfig = askList(confItems, 'SONSUZ odak konfigürasyonu:', 1, 'Nominal dosya');
    end
    if isempty(nearConfig)
        nearConfig = askList(confItems, 'YAKIN odak konfigürasyonu:', min(2, nConf), 'Nominal dosya');
    end
    if infConfig == nearConfig
        error('Sonsuz ve yakın odak konfigürasyonları aynı olamaz.');
    end

    mce.SetCurrentConfiguration(infConfig);
    [tInfAll, fInfAll] = readThicknessFocal(TheSystem);
    mce.SetCurrentConfiguration(nearConfig);
    [tNearAll, fNearAll] = readThicknessFocal(TheSystem);

    nNomSurf = numel(tInfAll);
    changed = abs(tNearAll - tInfAll) > 1e-12;
    isParax = ~isnan(fNearAll) | ~isnan(fInfAll);

    % MCE'den otomatik tespit: değeri değişen THIC satırı -> kompanzatör
    if isempty(nomCompSurface)
        isMove = false(numel(thicOps), 1);
        for i = 1:numel(thicOps)
            isMove(i) = abs(thicOps(i).val(nearConfig) - thicOps(i).val(infConfig)) > 1e-12;
        end
        cand = thicOps(isMove);
        if numel(cand) > 1
            names = arrayfun(@(r) sprintf('Yüzey %d: THIC %.6g -> %.6g', r.surf, r.val(infConfig), ...
                r.val(nearConfig)), cand, 'UniformOutput', false);
            cand = cand(askList(names, 'Birden fazla THIC satırı değişiyor. Kompanzatör hangisi?', 1, 'Nominal dosya'));
        end
        if isscalar(cand)
            nomCompSurface = cand.surf;
        end
    end

    % Tespit edilemediyse sor: kalınlığı değişen yüzeyler "*" ile işaretlenir
    if isempty(nomCompSurface)
        marks = repmat("", nNomSurf, 1);
        marks(changed) = "  * kalınlık değişiyor";
        items = surfaceItems(TheSystem, marks);
        init = find(changed, 1);
        if isempty(init), init = 1; end
        nomCompSurface = askList(items, ...
            {'Nominal dosyada KOMPANZATÖR yüzeyi (kalınlığı hareket eden yüzey).', ...
             '* = kalınlığı sonsuz ve yakın odak arasında değişen yüzeyler.', ...
             'Pickup ile bağlı yüzeyi değil, asıl hareket eden yüzeyi seçin.'}, init, 'Nominal dosya') - 1;
    end
    checkSurface(nomCompSurface, nNomSurf, 'nomCompSurface');

    % Yakın odak paraxial yüzeyi tespit edilemediyse sor
    if isempty(nomNearSurface)
        marks = repmat("", nNomSurf, 1);
        marks(isParax) = "  * paraxial";
        items = surfaceItems(TheSystem, marks);
        cand = find(isParax);
        if isempty(cand)
            init = 1;
        else
            [~, k] = max(abs(fNearAll(cand)));   % öneri: yakın config'te en büyük |f|
            init = cand(k);
        end
        nomNearSurface = askList(items, ...
            {'Nominal dosyada YAKIN ODAK paraxial yüzeyi.', ...
             '* = paraxial yüzeyler (değerler yakın odak konfigürasyonundan).'}, init, 'Nominal dosya') - 1;
    end
    checkSurface(nomNearSurface, nNomSurf, 'nomNearSurface');

    dComp_mm = (tNearAll(nomCompSurface + 1) - tInfAll(nomCompSurface + 1)) / nomUnits;
    fNear_mm = fNearAll(nomNearSurface + 1) / nomUnits;

    if abs(dComp_mm) < 1e-9
        error('Yüzey %d kalınlığı iki konfigürasyonda aynı: Δkomp = 0. Kompanzatör yüzeyini kontrol edin.', nomCompSurface);
    end
    if isnan(fNear_mm)
        error('Yüzey %d paraxial yüzey değil (yakın odak konfigürasyonunda).', nomNearSurface);
    end
    if fNear_mm == 0
        error('Yüzey %d focal length = 0 (yakın odak konfigürasyonunda).', nomNearSurface);
    end

    fprintf('Nominal dosya     : %s\n', nominalName);
    fprintf('Konfigürasyonlar  : sonsuz = %d, yakın = %d\n', infConfig, nearConfig);
    fprintf('Kompanzatör       : yüzey %d, Δkomp = %+.6f mm (yakın − sonsuz)\n', nomCompSurface, dComp_mm);
    fprintf('Yakın odak lensi  : yüzey %d, focal length = %.6g mm\n', nomNearSurface, fNear_mm);

    %% 2. MC dosyaları: kompanzatör ve |EFL| her dosyada ayrı bulunur
    fprintf('MC dosyaları      : %d dosya\n', nFiles);
    if isempty(mcCompSurface)
        fprintf('MC kompanzatör    : her dosyada kalınlığı variable olan yüzey\n');
    else
        checkSurface(mcCompSurface, Inf, 'mcCompSurface');
        fprintf('MC kompanzatör    : yüzey %d (ayardan)\n', mcCompSurface);
    end
    if isempty(mcInsertAfter)
        fprintf('Paraxial ekleme   : görüntü yüzeyinden hemen önce\n');
    else
        checkSurface(mcInsertAfter, Inf, 'mcInsertAfter');
        fprintf('Paraxial ekleme   : yüzey %d''den sonra (ayardan)\n', mcInsertAfter);
    end
    if isempty(efl_mm)
        fprintf('|EFL|             : her MC dosyasının kendi EFL''i\n');
    elseif ~isfinite(efl_mm) || efl_mm <= 0
        error('Geçersiz efl_mm: %g', efl_mm);
    else
        fprintf('|EFL|             : %.4f mm (ayardan)\n', efl_mm);
    end
    fprintf('Piksel boyutu     : %.2f um\n', pixelPitch_um);
    fprintf('Encoder hatası    : %s um\n\n', mat2str(compErrors_um));

    %% 3. Her MC dosyası için boresight hesabı
    checkPrinted = false;
    tStart = tic;
    wb = waitbar(0, 'Boresight hesaplanıyor...', 'Name', 'Boresight Analizi');

    for k = 1:nFiles
        try
            loadFile(TheSystem, fullfile(folderPath, files(k).name));
            lde = TheSystem.LDE;
            if ~logical(TheSystem.SystemData.Aperture.AFocalImageSpace)
                error('Dosya afokal modda değil.');
            end

            u = unitsPerMm(TheSystem);
            wave = primaryWave(TheSystem);

            % Bu dosyanın kompanzatörü ve |EFL|'i
            if isempty(mcCompSurface)
                compSurf = findVariableThickness(TheSystem);
            else
                compSurf = mcCompSurface;
            end
            compSurfUsed(k) = compSurf;
            if isempty(efl_mm)
                eflK = abs(TheSystem.MFE.GetOperandValue( ...
                    ZOSAPI.Editors.MFE.MeritOperandType.EFFL, 0, 0, 0, 0, 0, 0, 0, 0)) / u;
                if ~isfinite(eflK) || eflK == 0
                    error('EFL hesaplanamadı.');
                end
            else
                eflK = efl_mm;
            end
            eflUsed(k) = eflK;
            pixPerTan = eflK * 1000 / pixelPitch_um;   % tan θ -> piksel

            % a. Sonsuz odak (dosya olduğu gibi)
            t0 = lde.GetSurfaceAt(compSurf).Thickness;
            compT0(k) = t0 / u;
            b0 = rayTan(TheSystem, lde.NumberOfSurfaces - 1, wave, 0);
            tanInf(k, :) = b0;
            if ~checkPrinted
                cInf = collimation(TheSystem, lde.NumberOfSurfaces - 1, wave);
            end

            % b. Yakın odak: paraxial ekle, kompanzatörü kaydır
            if isempty(mcInsertAfter)
                insIdx = lde.NumberOfSurfaces - 1;   % görüntüden hemen önce
            else
                insIdx = mcInsertAfter + 1;
            end
            insertParaxial(TheSystem, insIdx, fNear_mm * u);
            compIdx = compSurf + (insIdx <= compSurf);   % ekleme önceyse numara kayar
            img = TheSystem.LDE.NumberOfSurfaces - 1;

            for lev = 1:nLevels
                d = levelErrors_um(lev) / 1000;   % mm
                if d == 0
                    signs = 0;
                else
                    signs = [1 -1];
                end
                best = -Inf;
                for s = signs
                    setThickness(TheSystem, compIdx, t0 + (dComp_mm + s * d) * u);
                    b = rayTan(TheSystem, img, wave, 0);
                    dxy = (b - b0) * pixPerTan;
                    if hypot(dxy(1), dxy(2)) > best
                        best = hypot(dxy(1), dxy(2));
                        dxPix(k, lev) = dxy(1);
                        dyPix(k, lev) = dxy(2);
                        appliedSign(k, lev) = s;
                    end
                end

                % Kontrol: yakın odakta çıkış kolime mi, focal length işareti doğru mu?
                if ~checkPrinted && lev == 1
                    setThickness(TheSystem, compIdx, t0 + dComp_mm * u);
                    cNear = collimation(TheSystem, img, wave);
                    setFocal(TheSystem, insIdx, -fNear_mm * u);
                    cFlip = collimation(TheSystem, img, wave);
                    setFocal(TheSystem, insIdx, fNear_mm * u);
                    fprintf('Kontrol (%s):\n', files(k).name);
                    fprintf('  Kompanzatör yüzey %d: %.6f -> %.6f mm (yakın odak)\n', ...
                        compIdx, t0 / u, t0 / u + dComp_mm);
                    fprintf('  Kolimasyon hatası (Py = +/-0.7 ışın açı farkı): sonsuz %.3g, yakın %.3g, ', cInf, cNear);
                    fprintf('yakın (f işareti ters) %.3g rad\n', cFlip);
                    fprintf('  Boresight kayması (0 um): dx = %+.4f, dy = %+.4f piksel\n\n', dxPix(k, 1), dyPix(k, 1));
                    if isfinite(cNear) && isfinite(cFlip) && cFlip < 0.5 * cNear
                        warning(['Focal length işareti ters çevrildiğinde çıkış daha iyi kolime oluyor ' ...
                            '(%.3g < %.3g rad). Yakın odak paraxial focal length işaretini veya Δkomp ' ...
                            'yönünü kontrol edin.'], cFlip, cNear);
                    end
                    checkPrinted = true;
                end
            end

            % c. Geri al (diske kaydedilmez)
            setThickness(TheSystem, compIdx, t0);
            TheSystem.LDE.RemoveSurfaceAt(insIdx);

            if any(isnan(dxPix(k, :))) || any(isnan(dyPix(k, :)))
                error('Geçersiz boresight değeri.');
            end
            status(k) = "OK";
        catch ME
            status(k) = string(ME.message);
            fprintf('HATA: %s - %s\n', files(k).name, ME.message);
        end

        if isvalid(wb)
            waitbar(k / nFiles, wb, sprintf('Boresight hesaplanıyor... %d / %d', k, nFiles));
        end
    end
    if isvalid(wb), close(wb); end

catch ME
    if exist('wb', 'var') && isvalid(wb), close(wb); end
    CleanupConnection(TheApplication);
    rethrow(ME);
end

CleanupConnection(TheApplication);

%% Geçerli sonuçlar ve istatistik
valid = status == "OK";
nValid = sum(valid);
if nValid == 0
    error('Geçerli sonuç alınamadı.');
end
fprintf('Geçerli dosya: %d / %d (süre: %.1f s)\n\n', nValid, nFiles, toc(tStart));

magPix = hypot(dxPix, dyPix);
magV = magPix(valid, :);
nYield = numel(yieldLevels);
yieldPix = zeros(nLevels, nYield);
for lev = 1:nLevels
    sortedMag = sort(magV(:, lev), 'ascend');
    for y = 1:nYield
        yieldPix(lev, y) = sortedMag(max(1, ceil(yieldLevels(y) / 100 * nValid)));
    end
end

levelNames = strings(nLevels, 1);
for lev = 1:nLevels
    if levelErrors_um(lev) == 0
        levelNames(lev) = "Encoder hatası yok";
    else
        levelNames(lev) = sprintf('Encoder ±%g µm', levelErrors_um(lev));
    end
end

Ozet = table(levelErrors_um', mean(magV, 1)', std(magV, 0, 1)', max(magV, [], 1)', ...
    mean(abs(dxPix(valid, :)), 1)', mean(abs(dyPix(valid, :)), 1)', ...
    'VariableNames', {'EncoderHatasi_um', 'Ortalama_piksel', 'Std_piksel', 'Max_piksel', ...
    'Ortalama_abs_dx_piksel', 'Ortalama_abs_dy_piksel'});
for y = 1:nYield
    Ozet.(sprintf('Yield%d_piksel', yieldLevels(y))) = yieldPix(:, y);
end

fprintf('Boresight hatası (sonsuz -> yakın odak), piksel = %.2f um, %d dosya:\n', pixelPitch_um, nValid);
disp(Ozet);
fprintf('Yield %%Y: dosyaların %%Y''sinde boresight hatası bu değerin altında.\n');

%% Çıktı klasörü: "BORESIGHT ANALİZİ" (char(304) = 'İ')
outFolder = fullfile(folderPath, ['BORESIGHT ANAL' char(304) 'Z' char(304)]);
if ~isfolder(outFolder)
    mkdir(outFolder);
end

%% Grafikler
colors = lines(nLevels);

% Kümülatif dağılım
figCdf = figure('Color', 'w', 'Position', [100 100 800 500]);
hold on;
for lev = 1:nLevels
    sortedMag = sort(magV(:, lev), 'ascend');
    stairs([0; sortedMag], [0; (1:nValid)' / nValid * 100], 'LineWidth', 1.5, ...
        'Color', colors(lev, :), 'DisplayName', levelNames(lev));
end
hold off;
grid on; box on;
xlabel('Boresight hatası (piksel)');
ylabel('Dosyaların yüzdesi (%) – hata bu değerin altında');
title(sprintf('Sonsuz → yakın odak boresight hatası (%d MC dosyası, piksel %.0f µm)', nValid, pixelPitch_um));
legend('Location', 'southeast');
ylim([0 100]);
set(findall(figCdf, '-property', 'FontName'), 'FontName', 'Arial');
exportPng(figCdf, fullfile(outFolder, 'Boresight_CDF.png'));

% dx - dy dağılımı
figSc = figure('Color', 'w', 'Position', [150 150 400 * nLevels 430]);
lim = max(magV(:)) * 1.1;
if lim == 0, lim = 1; end
for lev = 1:nLevels
    subplot(1, nLevels, lev);
    plot(dxPix(valid, lev), dyPix(valid, lev), '.', 'Color', colors(lev, :), 'MarkerSize', 10);
    hold on;
    th = linspace(0, 2 * pi, 200);
    r = yieldPix(lev, end);
    plot(r * cos(th), r * sin(th), 'k--', 'LineWidth', 1);
    hold off;
    axis equal; grid on; box on;
    xlim([-lim lim]); ylim([-lim lim]);
    xlabel('dx (piksel)'); ylabel('dy (piksel)');
    title({levelNames(lev), sprintf('kesikli: yield %d%% = %.3f piksel', yieldLevels(end), r)});
end
set(findall(figSc, '-property', 'FontName'), 'FontName', 'Arial');
exportPng(figSc, fullfile(outFolder, 'Boresight_Scatter.png'));

%% Excel
outputFile = fullfile(outFolder, 'Boresight_Results.xlsx');
if isfile(outputFile)
    delete(outputFile);
end
writetable(Ozet, outputFile, 'Sheet', 'Ozet');

Dosyalar = table(fileNames, status, compSurfUsed, compT0, eflUsed, tanInf(:, 1), tanInf(:, 2), ...
    'VariableNames', {'Dosya', 'Durum', 'Kompanzator_Yuzey', 'Kompanzator_T0_mm', 'EFL_mm', ...
    'Sonsuz_tanThetaX', 'Sonsuz_tanThetaY'});
for lev = 1:nLevels
    tag = sprintf('%gum', levelErrors_um(lev));
    Dosyalar.(['dx_piksel_' tag]) = dxPix(:, lev);
    Dosyalar.(['dy_piksel_' tag]) = dyPix(:, lev);
    Dosyalar.(['Hata_piksel_' tag]) = magPix(:, lev);
    if levelErrors_um(lev) ~= 0
        Dosyalar.(['UygulananIsaret_' tag]) = appliedSign(:, lev);
    end
end
writetable(Dosyalar, outputFile, 'Sheet', 'Dosyalar');

Ayarlar = table( ...
    ["Klasor"; "Nominal dosya"; "Sonsuz odak config"; "Yakin odak config"; ...
     "Nominal kompanzator yuzeyi"; "Delta komp (mm)"; "Nominal yakin odak paraxial yuzeyi"; ...
     "Yakin odak focal length (mm)"; "MC kompanzator yuzeyi"; "Paraxial eklendigi yer"; ...
     "|EFL| (mm)"; "Piksel (um)"; "Encoder hatalari (um)"; "MC dosya sayisi"; "Gecerli dosya"], ...
    [string(folderPath); string(nominalName); string(infConfig); string(nearConfig); ...
     string(nomCompSurface); sprintf('%.6f', dComp_mm); string(nomNearSurface); ...
     sprintf('%.6g', fNear_mm); settingText(mcCompSurface, "her dosyada variable kalinlik"); ...
     settingText(mcInsertAfter, "goruntuden hemen once"); ...
     sprintf('%.4f - %.4f', min(eflUsed(valid)), max(eflUsed(valid))); string(pixelPitch_um); ...
     string(mat2str(compErrors_um)); string(nFiles); string(nValid)], ...
    'VariableNames', {'Ayar', 'Deger'});
writetable(Ayarlar, outputFile, 'Sheet', 'Ayarlar');

save(fullfile(outFolder, 'Boresight_RawData.mat'), 'fileNames', 'status', 'levelErrors_um', ...
    'dxPix', 'dyPix', 'magPix', 'appliedSign', 'compT0', 'compSurfUsed', 'eflUsed', 'tanInf', 'dComp_mm', 'fNear_mm', ...
    'efl_mm', 'pixelPitch_um', 'infConfig', 'nearConfig', 'nomCompSurface', 'nomNearSurface', ...
    'mcCompSurface', 'mcInsertAfter', 'nominalName');

fprintf('\nÇıktılar: %s\n', outFolder);

%% ===================== YARDIMCI FONKSİYONLAR =====================

function loadFile(TheSystem, filePath)
if ~TheSystem.LoadFile(filePath, false)
    error('OpticStudio dosyayı açamadı: %s', filePath);
end
end

function [t, f] = readThicknessFocal(TheSystem)
% Aktif konfigürasyondaki her yüzeyin kalınlığı ve (paraxial ise) focal length'i.
% Dizi indeksi = yüzey numarası + 1.
lde = TheSystem.LDE;
n = lde.NumberOfSurfaces;
t = nan(n, 1);
f = nan(n, 1);
for i = 0:n - 1
    s = lde.GetSurfaceAt(i);
    t(i + 1) = s.Thickness;
    if isParaxialType(s)
        f(i + 1) = s.GetSurfaceCell(ZOSAPI.Editors.LDE.SurfaceColumn.Par1).DoubleValue;   % Focal Length
    end
end
end

function [thic, ignr] = readMceOperands(TheSystem)
% MCE'deki THIC ve IGNR satırları: yüzey numarası (Param1) ve her
% konfigürasyondaki değer.
mce = TheSystem.MCE;
nConf = mce.NumberOfConfigurations;
thic = struct('surf', {}, 'val', {});
ignr = struct('surf', {}, 'val', {});
for i = 1:mce.NumberOfOperands
    op = mce.GetOperandAt(i);
    typ = upper(strtrim(string(char(op.Type.ToString()))));
    if typ ~= "THIC" && typ ~= "IGNR"
        continue;
    end
    v = nan(1, nConf);
    for c = 1:nConf
        v(c) = cellValue(op.GetOperandCell(c));
    end
    r = struct('surf', double(op.Param1), 'val', v);
    if typ == "THIC"
        thic(end + 1) = r; %#ok<AGROW>
    else
        ignr(end + 1) = r; %#ok<AGROW>
    end
end
end

function v = cellValue(h)
% Hücre tipi Double, Integer ya da Boolean olabilir: önce metin değerinden oku
v = NaN;
try, v = str2double(strrep(char(h.Value), ',', '.')); catch, end
if isnan(v)
    try, v = double(h.IntegerValue); catch
        try, v = h.DoubleValue; catch, end
    end
end
if isnan(v)
    try
        txt = lower(strtrim(char(h.Value)));
        if any(strcmp(txt, {'true', 'on', 'yes'})), v = 1; end
        if any(strcmp(txt, {'false', 'off', 'no'})), v = 0; end
    catch
    end
end
end

function surfNo = findVariableThickness(TheSystem)
% Kalınlığı variable olan tek yüzey: MC dosyasının kompanzatörü
lde = TheSystem.LDE;
found = [];
for i = 0:lde.NumberOfSurfaces - 1
    solveType = string(lde.GetSurfaceAt(i).ThicknessCell.GetSolveData().Type.ToString());
    if solveType == "Variable"
        found(end + 1) = i; %#ok<AGROW>
    end
end
if isempty(found)
    error('Variable kalınlık yok: kompanzatör bulunamadı (mcCompSurface ayarı ile verilebilir).');
elseif numel(found) > 1
    error('%d variable kalınlık var (yüzey %s): tek kompanzatör bekleniyor (mcCompSurface ayarı ile verilebilir).', ...
        numel(found), mat2str(found));
end
surfNo = found;
end

function s = settingText(value, autoText)
if isempty(value)
    s = string(autoText);
else
    s = string(value);
end
end

function tf = isParaxialType(s)
name = lower(string(char(s.TypeName)));
tf = contains(name, "paraxial") && ~contains(name, "xy");
end

function items = surfaceItems(TheSystem, marks)
% Liste kutusunda gösterilecek satırlar (aktif konfigürasyonun değerleriyle)
lde = TheSystem.LDE;
n = lde.NumberOfSurfaces;
items = cell(n, 1);
for i = 0:n - 1
    s = lde.GetSurfaceAt(i);
    solve = string(s.ThicknessCell.GetSolveData().Type.ToString());
    if solve == "Fixed"
        solve = "";
    else
        solve = " (" + solve + ")";
    end
    extra = "";
    if isParaxialType(s)
        extra = sprintf('  f = %.6g', s.GetSurfaceCell(ZOSAPI.Editors.LDE.SurfaceColumn.Par1).DoubleValue);
    end
    items{i + 1} = char(sprintf("%3d | %s | %s | R = %.6g | T = %.6f%s | %s%s%s", i, ...
        string(char(s.TypeName)), strtrim(string(char(s.Comment))), s.Radius, s.Thickness, solve, ...
        strtrim(string(char(s.Material))), extra, marks(i + 1)));
end
end

function idx = askList(items, prompt, initial, titleText)
[idx, ok] = listdlg('ListString', items, 'SelectionMode', 'single', 'PromptString', prompt, ...
    'InitialValue', initial, 'ListSize', [750 420], 'Name', titleText);
if ~ok
    error('Seçim iptal edildi.');
end
end

function checkSurface(surfNo, nSurf, name)
if ~isscalar(surfNo) || surfNo < 0 || surfNo > nSurf - 1 || surfNo ~= round(surfNo)
    error('%s = %g geçersiz (0..%d olmalı).', name, surfNo, nSurf - 1);
end
end

function insertParaxial(TheSystem, idx, focalLength)
% idx konumuna paraxial yüzey ekler (thickness 0). Dosya diske kaydedilmez.
lde = TheSystem.LDE;
lde.InsertNewSurfaceAt(idx);
surf = lde.GetSurfaceAt(idx);
surf.ChangeType(surf.GetSurfaceTypeSettings(ZOSAPI.Editors.LDE.SurfaceType.Paraxial));
surf.Thickness = 0;
surf.Comment = 'Yakin odak (boresight script)';
setFocal(TheSystem, idx, focalLength);
end

function setFocal(TheSystem, idx, focalLength)
surf = TheSystem.LDE.GetSurfaceAt(idx);
focalCell = surf.GetSurfaceCell(ZOSAPI.Editors.LDE.SurfaceColumn.Par1);
focalCell.DoubleValue = focalLength;   % değişken üzerinden: zincirli atama uygulanmaz
actual = TheSystem.LDE.GetSurfaceAt(idx).GetSurfaceCell(ZOSAPI.Editors.LDE.SurfaceColumn.Par1).DoubleValue;
if abs(actual - focalLength) > 1e-9 * max(1, abs(focalLength))
    error('Paraxial focal length ayarlanamadı (istenen %.6g, sistemdeki %.6g).', focalLength, actual);
end
end

function setThickness(TheSystem, surfNo, value)
% Kalınlığı ayarlar ve sistemin gerçekten aldığını kontrol eder.
surf = TheSystem.LDE.GetSurfaceAt(surfNo);
surf.Thickness = value;
actual = TheSystem.LDE.GetSurfaceAt(surfNo).Thickness;
if abs(actual - value) > 1e-9
    error('Yüzey %d kalınlığı değiştirilemedi (istenen %.6f, sistemdeki %.6f). Kalınlık solve/pickup ile bağlı olabilir.', ...
        surfNo, value, actual);
end
end

function t = rayTan(TheSystem, surf, wave, py)
% On-axis (Hx = Hy = 0) ışının surf yüzeyinden sonraki yönü: [L/N, M/N] = [tan θx, tan θy]
mfe = TheSystem.MFE;
L = mfe.GetOperandValue(ZOSAPI.Editors.MFE.MeritOperandType.REAA, surf, wave, 0, 0, 0, py, 0, 0);
M = mfe.GetOperandValue(ZOSAPI.Editors.MFE.MeritOperandType.REAB, surf, wave, 0, 0, 0, py, 0, 0);
N = mfe.GetOperandValue(ZOSAPI.Editors.MFE.MeritOperandType.REAC, surf, wave, 0, 0, 0, py, 0, 0);
if ~all(isfinite([L M N])) || abs(N) < 1e-12
    error('Işın izlenemedi (yüzey %d, Py = %g).', surf, py);
end
t = [L / N, M / N];
end

function c = collimation(TheSystem, surf, wave)
% Kolimasyon göstergesi: Py = +0.7 ve -0.7 ışınlarının açı farkı (rad).
% Kolime çıkışta ~0; defokusta büyür.
try
    a = rayTan(TheSystem, surf, wave, 0.7);
    b = rayTan(TheSystem, surf, wave, -0.7);
    c = abs(atan(a(2)) - atan(b(2)));
catch
    c = NaN;
end
end

function w = primaryWave(TheSystem)
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

function f = unitsPerMm(TheSystem)
% 1 mm kaç lens birimi
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

function exportPng(fig, imageFile)
try
    exportgraphics(fig, imageFile, 'Resolution', 300);   % R2020a+
catch
    saveas(fig, imageFile);
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
