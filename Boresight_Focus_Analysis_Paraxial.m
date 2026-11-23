%% BORESIGHT ANALİZİ (PARAXIAL / AFOKAL OLMAYAN) – SONSUZ ODAK → YAKIN ODAK (ZOS-API, MATLAB)
% NE YAPAR
%   Boresight_Focus_Analysis.m'nin afokal OLMAYAN MC dosyaları için olan
%   karşılığı. Sistem sonsuz odaktan yakın odağa geçtiğinde (yalnızca
%   kompanzatör lens hareket eder) on-axis boresight'ın dedektörde kaç PİKSEL
%   kaydığını her Monte Carlo (MC) dosyası için hesaplar. Kompanzatörün yakın
%   odaktaki konumlandırma (encoder) hatası ±δ da eklenebilir.
%
% GİRDİLER (tek klasör)
%   - Nominal yakın odak dosyası: adı "MC_" ile başlamayan .zmx. En az iki
%     konfigürasyon (biri sonsuz odak, biri yakın odak). Yakın odak bir
%     paraxial yüzeyle simüle edilir; konfigürasyonlar arasında yalnızca
%     kompanzatör kalınlığı değişir.
%   - MC dosyaları: adı filePrefix ile başlayan .zmx. Afokal mod KAPALI,
%     sonsuz odak, tek konfigürasyon. Görüntü yüzeyinden hemen önce odaklayıcı
%     bir paraxial yüzey vardır: focal length = thickness = |EFL|, yani görüntü
%     yüzeyi bu lensin odak düzlemindedir (eşdeğer dedektör).
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
%     - Odaklayıcı paraxial yüzey: görüntü yüzeyinden önceki son yüzey.
%     - Yakın odak paraxial yüzeyi, odaklayıcı paraxial yüzeyin hemen önüne
%       eklenir (nominal yakın odak dosyasındaki düzen).
%   Tespit edilemezse ya da birden fazla aday varsa kullanıcıya liste ile
%   sorulur. Ayarlarda değer verilirse tespit yapılmaz, o değer kullanılır.
%
% NASIL ÇALIŞIR
%   1. Nominal dosyadan:
%        Δkomp  = T_komp(yakın config) − T_komp(sonsuz config)
%        f_yakın = yakın odak paraxial yüzeyinin focal length'i (yakın config)
%   2. Her MC dosyasında:
%      a. Sonsuz odak: dosya olduğu gibi. On-axis ana ışının (Hx=Hy=0,
%         Px=Py=0, ana dalga boyu) görüntü yüzeyindeki konumu okunur:
%         X, Y (REAX, REAY).
%      b. Yakın odak: odaklayıcı paraxial yüzeyin hemen önüne paraxial yüzey
%         eklenir (thickness 0, focal length f_yakın) ve kompanzatör kalınlığı
%         T0 + Δkomp (± δ) yapılır. Konum tekrar okunur.
%      c. Boresight hatası [piksel] = (konum_yakın − konum_sonsuz) / piksel
%         Görüntü yüzeyi odaklayıcı lensin odak düzleminde olduğu için konum
%         y = f·tan θ'dır; ışının lense çarptığı yükseklikten bağımsızdır.
%         Bu yüzden sonuç afokal koddaki açı × |EFL| hesabıyla aynıdır.
%      d. Eklenen yüzey silinir, kompanzatör geri alınır. Diske kaydedilmez.
%   3. Encoder hatası δ için +δ ve −δ ikisi de hesaplanır; boresight kaymasını
%      büyüten taraf (worst case) raporlanır. δ = 0 her zaman hesaplanır.
%   Hiçbir adımda optimizasyon yapılmaz.
%
% KONTROLLER (Command Window ve PDF)
%   - Odak düzlemi: odaklayıcı paraxial yüzeyin thickness'i focal length'ine
%     eşit olmalıdır. Eşit değilse ışının yüksekliği sonuca karışır; dosya
%     başına kaydedilir ve uyarı verilir.
%   - Odak: ilk geçerli MC dosyasında Py = +0.7 ve −0.7 ışınlarının görüntü
%     yüzeyindeki konum farkına bakılır (odakta ~0). Yakın odak lensinin
%     focal length işareti ters çevrildiğinde belirgin şekilde daha iyi odak
%     oluyorsa uyarı verilir (işaret ya da Δkomp yönü şüpheli).
%   - Kompanzatör değişikliğinin sisteme gerçekten uygulandığı kontrol edilir.
%
% ÇIKTILAR (seçilen klasörde "BORESIGHT ANALİZİ - PARAXIAL" klasörüne; afokal kodun
%           çıktılarının üzerine yazılmaz)
%   Boresight_Results.xlsx : Ozet (seviye başına ortalama, std, yield, max),
%                            Dosyalar (dosya başına dx, dy, büyüklük, uygulanan
%                            işaret), Ayarlar (okunan/seçilen tüm değerler)
%   Boresight_Report.pdf   : A4; özet tablosu ve kontrol sonucu, yöntem, kapsam ve
%                            yorum açıklamaları, grafik sayfaları (kümülatif dağılım,
%                            encoder hatasına göre yield, dx–dy dağılımı)
%   Boresight_Grafik_*.png : grafik sayfaları
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
mcFocusSurface = [];   % MC'de odaklayıcı paraxial (f = |EFL|) yüzeyi ([] = görüntüden önceki son yüzey)
                       % Yakın odak paraxial yüzeyi bu yüzeyin hemen önüne eklenir.

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
posInf = nan(nFiles, 2);             % sonsuz odakta ana ışının görüntüdeki konumu (X, Y) mm
compSurfUsed = nan(nFiles, 1);       % MC dosyasında kullanılan kompanzatör yüzeyi
focusSurfUsed = nan(nFiles, 1);      % MC dosyasındaki odaklayıcı paraxial yüzey
focusF = nan(nFiles, 1);             % odaklayıcı paraxial focal length (mm)
focusT = nan(nFiles, 1);             % odaklayıcı paraxial thickness (mm)

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
    % (kaynak metinleri PDF raporu için)
    confSource = "ayar / kullanıcı seçimi";
    nearSource = "ayar / kullanıcı seçimi";
    compSource = "ayar / kullanıcı seçimi";
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
            nearSource = "MCE, IGNR satırı";
            if isempty(infConfig)
                confSource = "MCE, IGNR satırı (ignore = sonsuz odak)";
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
            compSource = "MCE, THIC satırı";
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

    %% 2. MC dosyaları: kompanzatör ve odaklayıcı paraxial yüzey her dosyada ayrı bulunur
    fprintf('MC dosyaları      : %d dosya\n', nFiles);
    if isempty(mcCompSurface)
        fprintf('MC kompanzatör    : her dosyada kalınlığı variable olan yüzey\n');
    else
        checkSurface(mcCompSurface, Inf, 'mcCompSurface');
        fprintf('MC kompanzatör    : yüzey %d (ayardan)\n', mcCompSurface);
    end
    if isempty(mcFocusSurface)
        fprintf('Odaklayıcı lens   : görüntüden önceki son yüzey (paraxial olmalı)\n');
    else
        checkSurface(mcFocusSurface, Inf, 'mcFocusSurface');
        fprintf('Odaklayıcı lens   : yüzey %d (ayardan)\n', mcFocusSurface);
    end
    fprintf('Yakın odak lensi  : odaklayıcı paraxial yüzeyin hemen önüne eklenir\n');
    fprintf('Piksel boyutu     : %.2f um\n', pixelPitch_um);
    fprintf('Encoder hatası    : %s um\n\n', mat2str(compErrors_um));

    %% 3. Her MC dosyası için boresight hesabı
    checkPrinted = false;
    checkFile = "";          % kontrolün yapıldığı dosya (PDF raporu için)
    cInf = NaN; cNear = NaN; cFlip = NaN;
    signWarning = false;
    checkShift = [NaN NaN];
    focusWarned = false;     % odak düzlemi uyarısı yalnızca bir kez
    tStart = tic;
    wb = waitbar(0, 'Boresight hesaplanıyor...', 'Name', 'Boresight Analizi');

    for k = 1:nFiles
        try
            loadFile(TheSystem, fullfile(folderPath, files(k).name));
            lde = TheSystem.LDE;
            if logical(TheSystem.SystemData.Aperture.AFocalImageSpace)
                error('Dosya afokal modda. Afokal MC dosyaları için Boresight_Focus_Analysis.m kullanın.');
            end

            u = unitsPerMm(TheSystem);
            wave = primaryWave(TheSystem);
            pixPerUnit = 1000 / (u * pixelPitch_um);   % lens birimi -> piksel

            % Bu dosyanın kompanzatörü
            if isempty(mcCompSurface)
                compSurf = findVariableThickness(TheSystem);
            else
                compSurf = mcCompSurface;
            end
            compSurfUsed(k) = compSurf;

            % Odaklayıcı paraxial yüzey (f = |EFL|): görüntü yüzeyi onun odak düzleminde olmalı
            if isempty(mcFocusSurface)
                focSurf = lde.NumberOfSurfaces - 2;   % görüntüden önceki son yüzey
            else
                focSurf = mcFocusSurface;
            end
            focSurfObj = lde.GetSurfaceAt(focSurf);
            if ~isParaxialType(focSurfObj)
                error(['Yüzey %d paraxial değil (%s). Görüntüden önceki yüzey odaklayıcı paraxial ' ...
                    'olmalı ya da mcFocusSurface ayarı verilmeli.'], focSurf, char(focSurfObj.TypeName));
            end
            focusSurfUsed(k) = focSurf;
            focusF(k) = focSurfObj.GetSurfaceCell(ZOSAPI.Editors.LDE.SurfaceColumn.Par1).DoubleValue / u;
            focusT(k) = focSurfObj.Thickness / u;
            if focSurf ~= lde.NumberOfSurfaces - 2 && ~focusWarned
                warning(['Odaklayıcı paraxial yüzey (%d) ile görüntü yüzeyi arasında başka yüzeyler var; ' ...
                    'görüntü yüzeyinin odak düzleminde olduğunu kontrol edin.'], focSurf);
                focusWarned = true;
            end
            if abs(focusT(k) - focusF(k)) > 1e-6 * max(1, abs(focusF(k))) && ~focusWarned
                warning(['%s: odaklayıcı paraxial yüzeyin thickness''i (%.6f mm) focal length''ine ' ...
                    '(%.6f mm) eşit değil. Görüntü yüzeyi odak düzleminde değil; ışının yüksekliği ' ...
                    'sonuca karışır.'], files(k).name, focusT(k), focusF(k));
                focusWarned = true;
            end

            % a. Sonsuz odak (dosya olduğu gibi)
            t0 = lde.GetSurfaceAt(compSurf).Thickness;
            compT0(k) = t0 / u;
            img0 = lde.NumberOfSurfaces - 1;
            b0 = rayPos(TheSystem, img0, wave, 0);
            posInf(k, :) = b0 / u;
            if ~checkPrinted
                cInf = focusSpread(TheSystem, img0, wave) / u;
            end

            % b. Yakın odak: odaklayıcı paraxial yüzeyin önüne paraxial ekle, kompanzatörü kaydır
            insIdx = focSurf;                    % yeni yüzey odaklayıcı lensin hemen önüne gelir
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
                    b = rayPos(TheSystem, img, wave, 0);
                    dxy = (b - b0) * pixPerUnit;
                    if hypot(dxy(1), dxy(2)) > best
                        best = hypot(dxy(1), dxy(2));
                        dxPix(k, lev) = dxy(1);
                        dyPix(k, lev) = dxy(2);
                        appliedSign(k, lev) = s;
                    end
                end

                % Kontrol: yakın odakta görüntü odakta mı, focal length işareti doğru mu?
                if ~checkPrinted && lev == 1
                    setThickness(TheSystem, compIdx, t0 + dComp_mm * u);
                    cNear = focusSpread(TheSystem, img, wave) / u;
                    setFocal(TheSystem, insIdx, -fNear_mm * u);
                    cFlip = focusSpread(TheSystem, img, wave) / u;
                    setFocal(TheSystem, insIdx, fNear_mm * u);
                    checkFile = string(files(k).name);
                    checkShift = [dxPix(k, 1), dyPix(k, 1)];
                    signWarning = isfinite(cNear) && isfinite(cFlip) && cFlip < 0.5 * cNear;
                    fprintf('Kontrol (%s):\n', files(k).name);
                    fprintf('  Odaklayıcı paraxial: yüzey %d, f = %.6f mm, thickness = %.6f mm\n', ...
                        focSurf, focusF(k), focusT(k));
                    fprintf('  Kompanzatör yüzey %d: %.6f -> %.6f mm (yakın odak)\n', ...
                        compIdx, t0 / u, t0 / u + dComp_mm);
                    fprintf('  Odak hatası (Py = +/-0.7 ışın konum farkı): sonsuz %.3g, yakın %.3g, ', cInf, cNear);
                    fprintf('yakın (f işareti ters) %.3g mm\n', cFlip);
                    fprintf('  Boresight kayması (0 um): dx = %+.4g, dy = %+.4g piksel\n\n', dxPix(k, 1), dyPix(k, 1));
                    if signWarning
                        warning(['Focal length işareti ters çevrildiğinde görüntü daha iyi odaklanıyor ' ...
                            '(%.3g < %.3g mm). Yakın odak paraxial focal length işaretini veya Δkomp ' ...
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

%% Çıktı klasörü: "BORESIGHT ANALİZİ - PARAXIAL" (char(304) = 'İ')
outFolder = fullfile(folderPath, ['BORESIGHT ANAL' char(304) 'Z' char(304) ' - PARAXIAL']);
if ~isfolder(outFolder)
    mkdir(outFolder);
end

%% Grafik sayfaları (A4 dikey, sayfa başına iki grafik)
plotPos = {[0.14 0.58 0.74 0.30], [0.14 0.15 0.74 0.30]};   % üst / alt grafik
colors = lines(nLevels);
plotFigs = gobjects(0, 1);
maxMag = max(magV(:));
if maxMag == 0, maxMag = 1; end

% Sayfa: kümülatif dağılım (üst) ve encoder hatasına göre yield (alt)
fig = newA4Page();
addPageTitle(fig, sprintf('Boresight hatası dağılımı – sonsuz → yakın odak (%d MC dosyası)', nValid));

ax = axes(fig, 'Position', plotPos{1}, 'FontSize', 9, 'FontName', 'Arial');
hold(ax, 'on');
for lev = 1:nLevels
    sortedMag = sort(magV(:, lev), 'ascend');
    stairs(ax, [0; sortedMag], [0; (1:nValid)' / nValid * 100], 'LineWidth', 1.5, ...
        'Color', colors(lev, :), 'DisplayName', levelNames(lev));
end
hold(ax, 'off');
grid(ax, 'on'); box(ax, 'on');
xlim(ax, [0 maxMag * 1.05]);
ylim(ax, [0 100]);
xlabel(ax, 'Boresight hatası (piksel)');
ylabel(ax, 'Dosyaların yüzdesi (%)');
title(ax, 'Kümülatif dağılım: dosyaların ne kadarında hata bu değerin altında');
legend(ax, 'Location', 'southeast', 'FontSize', 8);

ax = axes(fig, 'Position', plotPos{2}, 'FontSize', 9, 'FontName', 'Arial');
hold(ax, 'on');
yColors = lines(nYield);
for y = 1:nYield
    plot(ax, levelErrors_um, yieldPix(:, y), '-o', 'Color', yColors(y, :), 'LineWidth', 1.5, ...
        'MarkerFaceColor', yColors(y, :), 'DisplayName', sprintf('Yield %d%%', yieldLevels(y)));
end
plot(ax, levelErrors_um, max(magV, [], 1), 'k--s', 'LineWidth', 1, 'DisplayName', 'Maksimum');
hold(ax, 'off');
grid(ax, 'on'); box(ax, 'on');
xPad = max(1, 0.05 * (max(levelErrors_um) - min(levelErrors_um)));
xlim(ax, [min(levelErrors_um) - xPad, max(levelErrors_um) + xPad]);
xticks(ax, levelErrors_um);
ylim(ax, [0 maxMag * 1.1]);
xlabel(ax, 'Encoder hatası ±δ (µm)');
ylabel(ax, 'Boresight hatası (piksel)');
title(ax, 'Encoder hatasına göre boresight hatası');
legend(ax, 'Location', 'northwest', 'FontSize', 8);
plotFigs(end + 1, 1) = fig;

% Sayfalar: dx - dy dağılımı (sayfa başına iki seviye)
scPos = {[0.22 0.55 0.56 0.33], [0.22 0.12 0.56 0.33]};
lim = maxMag * 1.15;
th = linspace(0, 2 * pi, 200);
for lev = 1:nLevels
    slot = mod(lev - 1, 2) + 1;
    if slot == 1
        fig = newA4Page();
        addPageTitle(fig, 'Boresight kayma yönü (dx – dy): her nokta bir MC dosyası');
        plotFigs(end + 1, 1) = fig; %#ok<SAGROW>
    end
    ax = axes(fig, 'Position', scPos{slot}, 'FontSize', 9, 'FontName', 'Arial');
    plot(ax, dxPix(valid, lev), dyPix(valid, lev), '.', 'Color', colors(lev, :), 'MarkerSize', 10);
    hold(ax, 'on');
    r = yieldPix(lev, end);
    plot(ax, r * cos(th), r * sin(th), 'k--', 'LineWidth', 1);
    plot(ax, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1);
    hold(ax, 'off');
    grid(ax, 'on'); box(ax, 'on');
    xlim(ax, [-lim lim]); ylim(ax, [-lim lim]);
    pbaspect(ax, [1 1 1]);
    xlabel(ax, 'dx (piksel)'); ylabel(ax, 'dy (piksel)');
    title(ax, {char(levelNames(lev)), sprintf('kesikli daire: yield %d%% = %.3f piksel', yieldLevels(end), r)});
end

%% Odaklayıcı paraxial yüzey özeti (Excel ve PDF için)
fMin = min(focusF(valid));
fMax = max(focusF(valid));
if abs(fMax - fMin) < 5e-5
    fText = string(sprintf('%.4f', fMin));
else
    fText = string(sprintf('%.4f - %.4f', fMin, fMax));
end
focusMismatch = valid & abs(focusT - focusF) > 1e-6 * max(1, abs(focusF));
if any(focusMismatch)
    focusOkText = string(sprintf('HAYIR: %d dosyada thickness ~= focal length', sum(focusMismatch)));
else
    focusOkText = "evet (tum gecerli dosyalarda)";
end

%% Excel
outputFile = fullfile(outFolder, 'Boresight_Results.xlsx');
if isfile(outputFile)
    delete(outputFile);
end
writetable(Ozet, outputFile, 'Sheet', 'Ozet');

Dosyalar = table(fileNames, status, compSurfUsed, compT0, focusSurfUsed, focusF, focusT, ...
    posInf(:, 1), posInf(:, 2), ...
    'VariableNames', {'Dosya', 'Durum', 'Kompanzator_Yuzey', 'Kompanzator_T0_mm', ...
    'Odaklayici_Yuzey', 'Odaklayici_f_mm', 'Odaklayici_Thickness_mm', 'Sonsuz_X_mm', 'Sonsuz_Y_mm'});
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
     "Yakin odak focal length (mm)"; "MC kompanzator yuzeyi"; "MC odaklayici paraxial yuzey"; ...
     "Odaklayici focal length (mm)"; "Odaklayici thickness = focal length?"; "Piksel (um)"; ...
     "Encoder hatalari (um)"; "MC dosya sayisi"; "Gecerli dosya"], ...
    [string(folderPath); string(nominalName); string(infConfig); string(nearConfig); ...
     string(nomCompSurface); sprintf('%.6f', dComp_mm); string(nomNearSurface); ...
     sprintf('%.6g', fNear_mm); settingText(mcCompSurface, "her dosyada variable kalinlik"); ...
     settingText(mcFocusSurface, "goruntuden onceki son yuzey"); fText; focusOkText; ...
     string(pixelPitch_um); string(mat2str(compErrors_um)); string(nFiles); string(nValid)], ...
    'VariableNames', {'Ayar', 'Deger'});
writetable(Ayarlar, outputFile, 'Sheet', 'Ayarlar');

save(fullfile(outFolder, 'Boresight_RawData.mat'), 'fileNames', 'status', 'levelErrors_um', ...
    'dxPix', 'dyPix', 'magPix', 'appliedSign', 'compT0', 'compSurfUsed', 'focusSurfUsed', 'focusF', ...
    'focusT', 'posInf', 'dComp_mm', 'fNear_mm', 'pixelPitch_um', 'infConfig', 'nearConfig', ...
    'nomCompSurface', 'nomNearSurface', 'mcCompSurface', 'mcFocusSurface', 'nominalName');

%% PDF raporu (A4): özet, yöntem, yorum, grafik sayfaları
pdfFile = fullfile(outFolder, 'Boresight_Report.pdf');

oldPages = dir(fullfile(outFolder, 'Boresight_Report_page*.pdf'));
for i = 1:numel(oldPages)
    delete(fullfile(outFolder, oldPages(i).name));
end
if isfile(pdfFile)
    delete(pdfFile);
end
% Eski rapor silinemediyse (ör. PDF okuyucuda açık) yeni rapor zaman damgalı adla kaydedilir
if isfile(pdfFile)
    pdfFile = fullfile(outFolder, ['Boresight_Report_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.pdf']);
    fprintf('\nBoresight_Report.pdf başka bir programda açık; yeni rapor: %s\n', pdfFile);
end

errListText = char(strjoin(string(compErrors_um), ', '));
mcCompText = char(strjoin(string(unique(compSurfUsed(valid))'), ', '));
focSurfText = char(strjoin(string(unique(focusSurfUsed(valid))'), ', '));
insertText = 'odaklayıcı paraxial yüzeyin hemen önüne';
yieldTop = yieldLevels(end);

% --- Sayfa: özet ---------------------------------------------------------
summaryLines = {
    sprintf('Tarih:                  %s', char(datetime('now', 'Format', 'dd.MM.yyyy HH:mm')))
    sprintf('Klasör:                 %s', folderPath)
    sprintf('Nominal yakın odak:     %s', nominalName)
    sprintf('MC dosyaları:           %s*  (geçerli: %d / %d)', filePrefix, nValid, nFiles)
    sprintf('Konfigürasyonlar:       sonsuz = %d, yakın = %d  (%s)', infConfig, nearConfig, char(confSource))
    sprintf('Kompanzatör (nominal):  yüzey %d, Δkomp = %+.4f mm  (%s)', nomCompSurface, dComp_mm, char(compSource))
    sprintf('Yakın odak lensi:       yüzey %d, f = %.6g mm  (%s)', nomNearSurface, fNear_mm, char(nearSource))
    sprintf('MC kompanzatör:         yüzey %s', mcCompText)
    sprintf('MC odaklayıcı lens:     yüzey %s, f = %s mm, thickness = f: %s', focSurfText, char(fText), char(focusOkText))
    sprintf('MC yakın odak lensi:    %s, thickness 0', insertText)
    'Mod:                    afokal kapalı, konum (REAX, REAY) odak düzleminde okunur'
    sprintf('Piksel boyutu:          %.2f µm', pixelPitch_um)
    sprintf('Encoder hatası:         ±%s µm (kötü taraf)', errListText)
    'Boresight:              on-axis ana ışın (Hx = Hy = 0, Px = Py = 0), ana dalga boyu'
    ''
    'BORESIGHT HATASI ÖZETİ (piksel) – sonsuz odak → yakın odak'
    ''
    };

header = sprintf('%-11s%10s%10s', 'Hata (µm)', 'Ortalama', 'Std');
for y = 1:nYield
    header = [header sprintf('%10s', sprintf('Y%d%%', yieldLevels(y)))]; %#ok<AGROW>
end
header = [header sprintf('%10s', 'Max')];
summaryLines = [summaryLines; {header; repmat('-', 1, numel(header))}];
for lev = 1:nLevels
    row = sprintf('%-11s%10.4f%10.4f', sprintf('±%g', levelErrors_um(lev)), Ozet.Ortalama_piksel(lev), Ozet.Std_piksel(lev));
    for y = 1:nYield
        row = [row sprintf('%10.4f', yieldPix(lev, y))]; %#ok<AGROW>
    end
    summaryLines{end + 1, 1} = [row sprintf('%10.4f', Ozet.Max_piksel(lev))]; %#ok<SAGROW>
end
summaryLines = [summaryLines; {
    ''
    'Y%: dosyaların %Y''sinde boresight hatası bu değerin altındadır.'
    sprintf('Örnek: ±0 µm satırında Y%d%% = %.4f piksel -> dosyaların %%%d''ında hata %.4f pikselden küçük.', ...
        yieldTop, yieldPix(1, end), yieldTop, yieldPix(1, end))
    ''
    sprintf('KONTROL (%s)', char(checkFile))
    sprintf('  Odak hatası: sonsuz %.3g mm, yakın %.3g mm, yakın (f işareti ters) %.3g mm', cInf, cNear, cFlip)
    }];
if signWarning
    summaryLines{end + 1, 1} = ['  UYARI: Focal length işareti ters çevrildiğinde görüntü daha iyi odaklanıyor. ' ...
        'Yakın odak lensinin işaretini veya Δkomp yönünü kontrol edin.'];
else
    summaryLines{end + 1, 1} = '  Yakın odak lensinin işareti ve Δkomp yönü tutarlı.';
end
summaryLines{end + 1, 1} = sprintf('  Boresight kayması (0 µm): dx = %+.4f, dy = %+.4f piksel', checkShift(1), checkShift(2));

failed = fileNames(~valid);
if ~isempty(failed)
    summaryLines = [summaryLines; {''; sprintf('Hesaplanamayan dosyalar (%d):', numel(failed))}];
    shown = find(~valid);
    shown = shown(1:min(end, 8));
    for i = shown'
        summaryLines{end + 1, 1} = sprintf('  %s: %s', fileNames(i), status(i)); %#ok<SAGROW>
    end
    if sum(~valid) > numel(shown)
        summaryLines{end + 1, 1} = sprintf('  ... ve %d dosya daha (Excel, Dosyalar sayfası)', sum(~valid) - numel(shown));
    end
end

summaryFig = textPage('BORESIGHT HATASI (PARAXIAL) – SONSUZ → YAKIN ODAK ANALİZ RAPORU', summaryLines, 'Courier New', 8.5);

% --- Sayfa: yöntem -------------------------------------------------------
explain1 = {
    'AMAÇ'
    ['Sistem sonsuz odaktan yakın odağa geçtiğinde yalnızca kompanzatör lens eksen boyunca hareket eder. ' ...
     'Nominal (hatasız) sistemde bu hareket görüş eksenini (boresight) değiştirmez. Üretilmiş bir ' ...
     'sistemde ise lensler decenter/tilt hatalı olduğu için kompanzatörün hareketi ve yakın hedefin ' ...
     'paralaksı görüş eksenini kaydırır: sonsuz odakta merkez pikselin baktığı hedef, yakın odakta ' ...
     'başka bir piksele düşer. Bu analiz, Monte Carlo (MC) dosyalarının her biri için bu kaymayı ' ...
     'piksel cinsinden hesaplar.']
    ''
    'ADIM 1 – NOMİNAL YAKIN ODAK DOSYASI'
    sprintf(['Nominal dosyada (%s) iki konfigürasyon vardır. Yakın odak, görüntü tarafına eklenen bir ' ...
     'paraxial yüzeyle simüle edilir ve konfigürasyonlar arasında yalnızca kompanzatör kalınlığı ' ...
     'değişir. Bu dosyadan iki değer okunur:'], nominalName)
    sprintf('•  Δkomp = kompanzatör kalınlığı (yakın) − (sonsuz) = %+.4f mm (yüzey %d).', dComp_mm, nomCompSurface)
    sprintf('•  Yakın odak lensinin focal length''i f = %.6g mm (yüzey %d, yakın odak konfigürasyonunda).', fNear_mm, nomNearSurface)
    ['Yakın odak paraxial yüzeyi MCE''deki IGNR satırından (ignore = sonsuz odak), kompanzatör ise iki ' ...
     'konfigürasyonda değeri değişen THIC satırından bulunur. Yorum (comment) alanı kullanılmaz.']
    ''
    'ADIM 2 – SONSUZ ODAK (HER MC DOSYASI)'
    ['MC dosyası afokal mod kapalı, kaydedildiği haliyle açılır (toleranslar uygulanmış, kompanzatör ' ...
     'sonsuz odağa optimize edilmiş). Görüntü yüzeyinden hemen önce odaklayıcı bir paraxial yüzey ' ...
     'vardır (focal length = thickness = |EFL|); görüntü yüzeyi bu lensin odak düzlemindedir. ' ...
     'Kompanzatör, dosyada kalınlığı variable olan yüzeydir (toleranslama coordinate break eklediği ' ...
     'için yüzey numarası her dosyada ayrıca bulunur). On-axis ana ışının (Hx = Hy = 0, Px = Py = 0, ' ...
     'ana dalga boyu) görüntü yüzeyindeki konumu okunur: X = REAX, Y = REAY.']
    ''
    'ADIM 3 – YAKIN ODAK (AYNI MC DOSYASI)'
    sprintf(['Aynı dosyaya, nominal dosyadaki gibi bir yakın odak paraxial yüzeyi eklenir (%s, ' ...
     'thickness 0, f = %.6g mm) ve kompanzatör kalınlığı Δkomp kadar değiştirilir. Böylece hatalar ' ...
     'aynen korunur, yalnızca odak değişir. Ana ışının konumu tekrar okunur. Hiçbir adımda optimizasyon ' ...
     'yapılmaz; hesaptan sonra eklenen yüzey silinir, kompanzatör geri alınır ve diskteki dosyalar ' ...
     'değiştirilmez.'], insertText, fNear_mm)
    ''
    'ADIM 4 – PİKSELE ÇEVİRME'
    '   Boresight hatası [piksel] = (REAY_yakın − REAY_sonsuz) / piksel boyutu   (x için REAX)'
    ['Görüntü yüzeyi odaklayıcı lensin odak düzleminde olduğu için konum y = f · tan θ''dır ve ışının ' ...
     'lense çarptığı yükseklikten bağımsızdır; paralel kayma elenir. Sonuç, afokal dosyada açıyı ' ...
     'okuyup |EFL| ile çarpmakla aynıdır. Bunun için odaklayıcı lensin thickness''i focal length''ine ' ...
     'eşit olmalıdır; kod bunu her dosyada kontrol eder. dx, dy ve büyüklük √(dx² + dy²) ayrı ayrı ' ...
     'kaydedilir.']
    ''
    'ADIM 5 – ENCODER (KOMPANZATÖR KONUM) HATASI'
    sprintf(['Yakın odakta kompanzatör, mekanizmanın konumlandırma hatası kadar yanlış konumlanabilir. ' ...
     'Her hata payı δ (%s µm) için kompanzatör Δkomp + δ ve Δkomp − δ konumlarına getirilir; ' ...
     'boresight hatasını büyüten taraf (worst case) o dosya için kullanılır. Uygulanan işaret ' ...
     'Excel''deki Dosyalar sayfasındadır.'], errListText)
    ''
    'ADIM 6 – İSTATİSTİK'
    ['Her hata payı için tüm geçerli dosyaların boresight hatası istatistiğe sokulur: ortalama, ' ...
     'standart sapma, maksimum ve yield değerleri. Yield %Y değeri, dosyaların %Y''sinin altında ' ...
     'kaldığı boresight hatasıdır; üretim açısından "sistemlerin %Y''sinde boresight kayması bu ' ...
     'değerden küçüktür" anlamına gelir.']
    };

% --- Sayfa: yorum --------------------------------------------------------
explain2 = {
    'KONTROL: YAKIN ODAK DOĞRU SİMÜLE EDİLDİ Mİ?'
    ['İlk geçerli MC dosyasında görüntünün odakta olup olmadığına bakılır: Py = +0.7 ve Py = −0.7 ' ...
     'ışınlarının görüntü yüzeyindeki konum farkı (odak hatası). Yakın odakta bu değer sonsuz ' ...
     'odaktakine yakın olmalıdır. Ayrıca yakın odak lensinin focal length işareti ters çevrilerek ' ...
     'aynı değer hesaplanır; ters işaret belirgin şekilde daha iyi odak veriyorsa (yarısından küçük) ' ...
     'uyarı verilir. Sistem dedektörden sahneye doğru tasarlandığında, sonlu uzaklıktaki bir sahne ' ...
     'noktası çıkışta yakınsayan bir ışın demektir; bunu kolime eden lensin işareti, nominal ' ...
     'dosyadaki Δkomp yönüyle tutarlı olmalıdır.']
    ''
    'KONTROL: GÖRÜNTÜ YÜZEYİ ODAK DÜZLEMİNDE Mİ?'
    ['Odaklayıcı paraxial yüzeyin thickness''i (t) focal length''ine (f) eşit değilse dedektördeki ' ...
     'konum y_det = t · tan θ + y · (1 − t/f) olur ve ışının lense çarptığı yükseklik (y) sonuca ' ...
     'karışır; paralel kayma görüntü kayması gibi görünür. Özet sayfasındaki "thickness = f" satırı ' ...
     'bunu gösterir. Tolerans analizinde bu thickness kompanzatör olarak kullanılmamalıdır.']
    ''
    'BU ANALİZ NEYİ KAPSAR?'
    ['Tolerans analizinin MC dosyalarına uyguladığı tüm hatalar (lens decenter/tilt, yüzey hataları, ' ...
     'boşluklar vb.) sonsuz ve yakın odakta AYNI kalır. Hesaplanan kayma, bu hatalı sistemde ' ...
     'kompanzatörün eksenel hareketinden ve yakın hedefin paralaksından kaynaklanan boresight ' ...
     'değişimidir. Bu kayma her seferinde aynı olduğu için (sistematik) odak konumuna göre bir ' ...
     'kalibrasyon tablosuyla yazılımda düzeltilebilir.']
    ''
    'BU ANALİZ NEYİ KAPSAMAZ?'
    ['•  Kızak hatası: Lensin kızak boyunca hareket ederken farklı decenter/tilt alması (kızağın ' ...
     'optik eksene göre açısı ve düzlüğü). MC dosyalarında kompanzatörün hatası her iki odakta ' ...
     'aynıdır. Bu etki için MCE''de konfigürasyon başına decenter/tilt (TMCO) toleransı gerekir.']
    ['•  Gidiş-dönüş (tekrarlanabilirlik) hatası: Lens yakın odaktan sonsuza döndüğünde ilk konumuna ' ...
     'tam oturmaması (boşluk, histerezis). Bu hata rastgeledir ve kalibre edilemez.']
    ['•  Yalnızca on-axis ana ışın kullanılır; diğer fieldler ve görüntü centroid''i (koma etkisi) ' ...
     'bu sürümde yoktur.']
    ''
    'GRAFİKLER NASIL OKUNUR?'
    ['•  Kümülatif dağılım: Yatay eksende boresight hatası (piksel), dikey eksende bu değerin altında ' ...
     'kalan dosyaların yüzdesi. Eğri sola ne kadar yakınsa sistem o kadar iyidir. Eğrinin %90 ' ...
     'çizgisini kestiği yer, yield %90 değeridir.']
    ['•  Encoder hatasına göre boresight hatası: Her yield seviyesinin ve maksimum değerin hata ' ...
     'payıyla değişimi. Eğri yataysa encoder hatası boresight''ı pratikte etkilemiyordur.']
    sprintf(['•  dx – dy dağılımı: Her nokta bir MC dosyasının kayma vektörüdür. Noktaların merkez ' ...
     'etrafında rastgele dağılması beklenir; tek bir yönde toplanması sistematik bir hatayı ' ...
     'gösterir. Kesikli daire yield %%%d değeridir.'], yieldTop)
    ''
    'SPESİFİKASYON NASIL KONUR?'
    '1.  Kabul edilebilir boresight hatasını piksel cinsinden belirle (ör. 0.5 piksel).'
    ['2.  Özet tablosunda ilgili yield seviyesindeki (ör. Y90%) değeri bu sınırla karşılaştır. ' ...
     'Encoder hatası satırları, mekanizmanın konumlandırma hassasiyetinin sonuca etkisini gösterir.']
    ['3.  Sınır aşılıyorsa: kalibrasyonla düzeltme (sistematik kısım), hareketli lensin ve ' ...
     'arkasındaki lenslerin decenter/tilt toleranslarının sıkılaştırılması ya da mekanizmanın ' ...
     'iyileştirilmesi değerlendirilir. Hangi toleransın en çok katkı verdiğini görmek için Zemax''te ' ...
     'boresight farkını kriter alan bir sensitivity analizi yapılabilir.']
    };

explainFig1 = textPage('1. ANALİZ YÖNTEMİ', explain1, 'Arial', 9.5);
explainFig2 = textPage('2. KONTROL, KAPSAM VE YORUM', explain2, 'Arial', 9.5);

% --- Sayfaları birleştir: çerçeve, başlık ve sayfa numaraları -------------
allFigs = [summaryFig; explainFig1; explainFig2; plotFigs];
nPages = numel(allFigs);
headerText = sprintf('Boresight Hatası (paraxial, REAX/REAY) – Sonsuz → Yakın Odak  |  %s', nominalName);
for p = 1:nPages
    addPageFrame(allFigs(p), headerText, p, nPages);
end

% Grafik sayfalarının PNG'leri
for p = 1:numel(plotFigs)
    exportPng(plotFigs(p), fullfile(outFolder, sprintf('Boresight_Grafik_%d.png', p)));
end

try
    % Çok sayfalı PDF (exportgraphics 'Append' R2021b+ gerektirir)
    exportgraphics(allFigs(1), pdfFile, 'ContentType', 'vector');
    for p = 2:nPages
        exportgraphics(allFigs(p), pdfFile, 'ContentType', 'vector', 'Append', true);
    end
    fprintf('\nPDF raporu: %s\n', pdfFile);
catch pdfErr
    % Eski MATLAB: her sayfa ayrı A4 PDF
    fprintf('\nÇok sayfalı PDF oluşturulamadı: %s\n', pdfErr.message);
    if isfile(pdfFile)
        delete(pdfFile);
    end
    for p = 1:nPages
        print(allFigs(p), fullfile(outFolder, sprintf('Boresight_Report_page%d.pdf', p)), '-dpdf');
    end
    fprintf('Rapor ayrı sayfalar olarak kaydedildi (Boresight_Report_page*.pdf).\n');
end

close(allFigs);

fprintf('\nÇıktılar: %s\n', outFolder);

%% ===================== YARDIMCI FONKSİYONLAR =====================

function fig = newA4Page()
% Boş A4 dikey sayfa (21 x 29.7 cm). Görünmez; boyutu ekranla sınırlanmaz.
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
    'Position', [1 1 21 29.7], 'PaperUnits', 'centimeters', ...
    'PaperSize', [21 29.7], 'PaperPosition', [0 0 21 29.7]);
end

function addPageTitle(fig, titleText)
% Arial açıkça verilir: varsayılan fontta PDF'te Türkçe karakterler yok
annotation(fig, 'textbox', [0.06 0.915 0.88 0.04], 'String', titleText, ...
    'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'Interpreter', 'none');
end

function fig = textPage(titleText, lines, fontName, fontSize)
% Başlık ve kaydırılmış metin içeren A4 sayfa (her hücre bir paragraf/satır)
fig = newA4Page();
annotation(fig, 'textbox', [0.08 0.89 0.84 0.05], 'String', titleText, ...
    'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'Interpreter', 'none', 'VerticalAlignment', 'top');
annotation(fig, 'textbox', [0.08 0.06 0.84 0.83], 'String', lines, ...
    'FontName', fontName, 'FontSize', fontSize, 'EdgeColor', 'none', ...
    'Interpreter', 'none', 'VerticalAlignment', 'top', 'FitBoxToText', 'off');
end

function addPageFrame(fig, headerText, pageNo, nPages)
% Üst bilgi, sayfa numaralı alt bilgi ve neredeyse beyaz sayfa çerçevesi. Çerçeve,
% exportgraphics içeriğe göre kırptığında tüm A4 sayfanın korunmasını sağlar.
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

function p = rayPos(TheSystem, surf, wave, py)
% On-axis (Hx = Hy = 0) ışının surf yüzeyindeki konumu: [X, Y] (lens birimi), REAX ve REAY.
% REAC yalnızca ışının izlenebildiğini kontrol etmek için okunur (izlenemezse 0 döner).
mfe = TheSystem.MFE;
X = mfe.GetOperandValue(ZOSAPI.Editors.MFE.MeritOperandType.REAX, surf, wave, 0, 0, 0, py, 0, 0);
Y = mfe.GetOperandValue(ZOSAPI.Editors.MFE.MeritOperandType.REAY, surf, wave, 0, 0, 0, py, 0, 0);
N = mfe.GetOperandValue(ZOSAPI.Editors.MFE.MeritOperandType.REAC, surf, wave, 0, 0, 0, py, 0, 0);
if ~all(isfinite([X Y N])) || abs(N) < 1e-12
    error('Işın izlenemedi (yüzey %d, Py = %g).', surf, py);
end
p = [X, Y];
end

function c = focusSpread(TheSystem, surf, wave)
% Odak göstergesi: Py = +0.7 ve -0.7 ışınlarının görüntü yüzeyindeki y konum
% farkı (lens birimi). Odakta ~0; defokusta büyür.
try
    a = rayPos(TheSystem, surf, wave, 0.7);
    b = rayPos(TheSystem, surf, wave, -0.7);
    c = abs(a(2) - b(2));
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
