%% BORESIGHT ANALİZİ – SONSUZ ODAK → YAKIN ODAK (ZOS-API, MATLAB)
% NE YAPAR
%   Sistem sonsuz odaktan yakın odağa geçtiğinde (yalnızca kompanzatör lens
%   hareket eder) görüntünün dedektörde kaç PİKSEL kaydığını her Monte Carlo
%   (MC) dosyası için hesaplar:
%     - On-axis (merkez piksel) = boresight, ayrıca dosyadaki tüm fieldler
%     - İki yöntemle: ana ışın (chief ray) ve centroid (tüm dalga boyları)
%     - Kompanzatörün yakın odaktaki konumlandırma (encoder) hatası ±δ ile
%   Ayrıca nominal tasarımın kendi kayması (odak nefesi) referans olarak verilir.
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
%      Nominal tasarımın kendi kayması: iki konfigürasyonda görüntü
%      yüzeyindeki ışın konumları (X, Y) farkı / piksel boyutu.
%   2. Her MC dosyasında (batch ray trace, gerçek ışınlar):
%      a. Sonsuz odak: dosya olduğu gibi. Her field için ana ışın (Px=Py=0,
%         ana dalga boyu) ve centroid ışınları (pupil ızgarası, tüm dalga
%         boyları) izlenir; görüntü uzayındaki yön tan θx = L/N, tan θy = M/N.
%         Centroid = ışınların tan θ değerlerinin dalga boyu ağırlıklı ortalaması.
%      b. Yakın odak: görüntü yüzeyinden hemen önce paraxial yüzey eklenir
%         (thickness 0, focal length f_yakın) ve kompanzatör kalınlığı
%         T0 + Δkomp (± δ) yapılır. Işınlar tekrar izlenir.
%      c. Kayma [piksel] = (tan θ_yakın − tan θ_sonsuz) × |EFL| / piksel
%         Bu, dedektör tarafında f = |EFL| paraxial lens ile yapılan hesabın
%         aynısıdır (y = f·tan θ; ışının yüksekliğinden bağımsız).
%      d. Eklenen yüzey silinir, kompanzatör geri alınır. Diske kaydedilmez.
%   3. Encoder hatası δ için +δ ve −δ ikisi de hesaplanır; her field ve her
%      yöntem için kaymayı büyüten taraf (worst case) raporlanır. δ = 0 her
%      zaman hesaplanır. Hiçbir adımda optimizasyon yapılmaz.
%
% KONTROLLER (ilk geçerli MC dosyasında, Command Window ve PDF)
%   - Kolimasyon: yakın odak durumunda çıkış kolime olmalıdır. Paraxial
%     focal length'in işareti ters çevrildiğinde çıkış belirgin şekilde daha
%     iyi kolime oluyorsa uyarı verilir (işaret ya da Δkomp yönü şüpheli).
%   - Kompanzatör değişikliğinin sisteme gerçekten uygulandığı kontrol edilir.
%
% ÇIKTILAR (seçilen klasörde "BORESIGHT ANALİZİ" klasörüne)
%   Boresight_Report.pdf   : A4; özet tabloları ve kontrol sonucu, yöntem,
%                            hesaplama detayları (operandlar, formüller),
%                            kapsam ve yorum, grafik sayfaları
%   Boresight_Results.xlsx : Ozet, Nominal, Dosyalar, Ayarlar sayfaları
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
fieldsToUse   = [];              % [] = dosyadaki tüm fieldler; ör. [1 3 5] = yalnızca bunlar
                                 % (on-axis (0,0) her zaman 1. sırada hesaplanır)
pupilGrid     = 12;              % centroid için pupil ızgarası (çap boyunca ışın sayısı)

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
methodNames = ["Centroid", "Ana ışın"];   % m = 1: centroid, m = 2: ana ışın
methodCols = {3:4, 1:2};                 % traceFields çıktısındaki sütunlar
fileNames = string({files.name})';
status = strings(nFiles, 1);
DX = [];                             % kayma x (piksel): dosya x seviye x field x yöntem
DY = [];                             % kayma y (piksel)
SGN = [];                            % uygulanan encoder hatası işareti (+1 / -1, 0: hatasız)
compT0 = nan(nFiles, 1);             % MC dosyasındaki kompanzatör kalınlığı (mm)
compSurfUsed = nan(nFiles, 1);       % MC dosyasında kullanılan kompanzatör yüzeyi
eflUsed = nan(nFiles, 1);            % MC dosyasında kullanılan |EFL| (mm)
nF = 0;                              % field sayısı (ilk geçerli MC dosyasından)
fieldLabels = strings(0, 1);
fieldTypeName = "";

% Centroid için pupil ızgarası: hücre merkezleri, birim çemberin içi
g = ((1:pupilGrid) - 0.5) / pupilGrid * 2 - 1;
[PX, PY] = meshgrid(g, g);
inside = PX .^ 2 + PY .^ 2 <= 1;
pupilPx = PX(inside);
pupilPy = PY(inside);

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
            [~, i] = max(abs(fNearAll(cand)));   % öneri: yakın config'te en büyük |f|
            init = cand(i);
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

    % Nominal tasarımın kendi kayması (odak nefesi): iki konfigürasyonda görüntü
    % yüzeyindeki ışın konumları. Afokal değilse konum (X, Y), afokalse yön (tan θ).
    nomAfocal = logical(TheSystem.SystemData.Aperture.AFocalImageSpace);
    [nomHx, nomHy, nomFieldLabels] = fieldList(TheSystem, fieldsToUse);
    [nomWIdx, nomWW, nomWave] = waveList(TheSystem);
    nomImg = TheSystem.LDE.NumberOfSurfaces - 1;
    mce.SetCurrentConfiguration(infConfig);
    Rinf = traceFields(TheSystem, nomImg, nomHx, nomHy, nomWave, nomWIdx, nomWW, pupilPx, pupilPy, nomAfocal);
    mce.SetCurrentConfiguration(nearConfig);
    Rnear = traceFields(TheSystem, nomImg, nomHx, nomHy, nomWave, nomWIdx, nomWW, pupilPx, pupilPy, nomAfocal);
    if nomAfocal
        nomFactor = abs(TheSystem.MFE.GetOperandValue( ...
            ZOSAPI.Editors.MFE.MeritOperandType.EFFL, 0, 0, 0, 0, 0, 0, 0, 0)) / nomUnits * 1000 / pixelPitch_um;
    else
        nomFactor = 1000 / (nomUnits * pixelPitch_um);   % lens birimi -> piksel
    end
    Dnom = (Rnear - Rinf) * nomFactor;
    nomNF = numel(nomHx);
    nomDX = zeros(nomNF, 2);
    nomDY = zeros(nomNF, 2);
    for m = 1:2
        nomDX(:, m) = Dnom(:, methodCols{m}(1));
        nomDY(:, m) = Dnom(:, methodCols{m}(2));
    end
    nomMag = hypot(nomDX, nomDY);
    fprintf('Nominal tasarım   : on-axis kayma centroid %.3g, ana ışın %.3g piksel (≈ 0 beklenir)\n', ...
        nomMag(1, 1), nomMag(1, 2));

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
    fprintf('Encoder hatası    : %s um\n', mat2str(compErrors_um));
    fprintf('Centroid ışınları : %d pupil noktası x dalga boyu sayısı (her field)\n\n', numel(pupilPx));

    %% 3. Her MC dosyası için kayma hesabı
    checkPrinted = false;
    checkFile = "";          % kontrolün yapıldığı dosya (PDF raporu için)
    cInf = NaN; cNear = NaN; cFlip = NaN;
    signWarning = false;
    checkShift = nan(2, 2);  % on-axis, 0 um: [dx dy] x yöntem
    diagLines = {};          % ilk dosyanın ham değerleri (Command Window ve PDF)
    tStart = tic;
    wb = waitbar(0, 'Boresight hesaplanıyor...', 'Name', 'Boresight Analizi');

    for k = 1:nFiles
        insIdx = [];
        try
            loadFile(TheSystem, fullfile(folderPath, files(k).name));
            lde = TheSystem.LDE;
            if ~logical(TheSystem.SystemData.Aperture.AFocalImageSpace)
                error('Dosya afokal modda değil.');
            end

            u = unitsPerMm(TheSystem);
            [wIdx, wW, wave] = waveList(TheSystem);
            [hx, hy, labels] = fieldList(TheSystem, fieldsToUse);

            % İlk geçerli dosya field sayısını belirler
            if nF == 0
                nF = numel(hx);
                fieldLabels = labels;
                fieldTypeName = string(char(TheSystem.SystemData.Fields.GetFieldType().ToString()));
                DX = nan(nFiles, nLevels, nF, 2);
                DY = nan(nFiles, nLevels, nF, 2);
                SGN = zeros(nFiles, nLevels, nF, 2);
            elseif numel(hx) ~= nF
                error('Field sayısı (%d) ilk dosyadan (%d) farklı.', numel(hx), nF);
            end

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
            img0 = lde.NumberOfSurfaces - 1;
            R0 = traceFields(TheSystem, img0, hx, hy, wave, wIdx, wW, pupilPx, pupilPy, true);
            if ~checkPrinted
                cInf = collimation(TheSystem, img0, wave);
                opInf = rayTan(TheSystem, img0, wave, 0);   % batch'ten bağımsız kontrol (REAA/REAB/REAC)
                nSurfBefore = lde.NumberOfSurfaces;
            end

            % b. Yakın odak: paraxial ekle, kompanzatörü kaydır
            if isempty(mcInsertAfter)
                insIdx = img0;                       % görüntüden hemen önce
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
                for s = signs
                    setThickness(TheSystem, compIdx, t0 + (dComp_mm + s * d) * u);
                    R = traceFields(TheSystem, img, hx, hy, wave, wIdx, wW, pupilPx, pupilPy, true);
                    D = (R - R0) * pixPerTan;
                    % Her field ve yöntem için kaymayı büyüten taraf tutulur
                    for m = 1:2
                        c = methodCols{m};
                        for f = 1:nF
                            cur = hypot(DX(k, lev, f, m), DY(k, lev, f, m));
                            new = hypot(D(f, c(1)), D(f, c(2)));
                            if isnan(cur) || new > cur
                                DX(k, lev, f, m) = D(f, c(1));
                                DY(k, lev, f, m) = D(f, c(2));
                                SGN(k, lev, f, m) = s;
                            end
                        end
                    end
                end

                % Kontrol: yakın odakta çıkış kolime mi, focal length işareti doğru mu?
                if ~checkPrinted && lev == 1
                    setThickness(TheSystem, compIdx, t0 + dComp_mm * u);
                    opNear = rayTan(TheSystem, img, wave, 0);
                    insSurf = TheSystem.LDE.GetSurfaceAt(insIdx);
                    diagLines = {
                        sprintf('  |EFL| = %.6g mm, tan θ -> piksel çarpanı = %.6g', eflK, pixPerTan)
                        sprintf('  Yüzey sayısı: %d -> %d, eklenen yüzey %d: %s, f = %.6g', nSurfBefore, ...
                            TheSystem.LDE.NumberOfSurfaces, insIdx, char(insSurf.TypeName), ...
                            insSurf.GetSurfaceCell(ZOSAPI.Editors.LDE.SurfaceColumn.Par1).DoubleValue)
                        sprintf('  Kompanzatör yüzey %d kalınlığı: %.6f -> %.6f (okunan)', compIdx, t0, ...
                            TheSystem.LDE.GetSurfaceAt(compIdx).Thickness)
                        sprintf('  Ana ışın tan θ (batch)  sonsuz: (%.6e, %.6e)  yakın: (%.6e, %.6e)', ...
                            R0(1, 1), R0(1, 2), R(1, 1), R(1, 2))
                        sprintf('  Ana ışın tan θ (REAA)   sonsuz: (%.6e, %.6e)  yakın: (%.6e, %.6e)', ...
                            opInf(1), opInf(2), opNear(1), opNear(2))
                        sprintf('  Centroid tan θ (batch)  sonsuz: (%.6e, %.6e)  yakın: (%.6e, %.6e)', ...
                            R0(1, 3), R0(1, 4), R(1, 3), R(1, 4))
                        };
                    cNear = collimation(TheSystem, img, wave);
                    setFocal(TheSystem, insIdx, -fNear_mm * u);
                    cFlip = collimation(TheSystem, img, wave);
                    setFocal(TheSystem, insIdx, fNear_mm * u);
                    checkFile = string(files(k).name);
                    for m = 1:2
                        checkShift(:, m) = [DX(k, 1, 1, m); DY(k, 1, 1, m)];
                    end
                    signWarning = isfinite(cNear) && isfinite(cFlip) && cFlip < 0.5 * cNear;
                    fprintf('Kontrol (%s):\n', files(k).name);
                    fprintf('  Kompanzatör yüzey %d: %.6f -> %.6f mm (yakın odak)\n', ...
                        compIdx, t0 / u, t0 / u + dComp_mm);
                    fprintf('  Kolimasyon hatası (Py = +/-0.7 ışın açı farkı): sonsuz %.3g, yakın %.3g, ', cInf, cNear);
                    fprintf('yakın (f işareti ters) %.3g rad\n', cFlip);
                    fprintf('  On-axis kayma (0 um): centroid dx = %+.4g, dy = %+.4g | ana ışın dx = %+.4g, dy = %+.4g piksel\n', ...
                        checkShift(1, 1), checkShift(2, 1), checkShift(1, 2), checkShift(2, 2));
                    fprintf('  Tanı değerleri:\n');
                    fprintf('%s\n', diagLines{:});
                    fprintf('\n');
                    if signWarning
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

            if any(isnan(DX(k, :, :, :)), 'all') || any(isnan(DY(k, :, :, :)), 'all')
                error('Geçersiz kayma değeri (ışın izlenemedi ya da tüm centroid ışınları vignette oldu).');
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

MAG = hypot(DX, DY);
Mv = MAG(valid, :, :, :);
nYield = numel(yieldLevels);
statMean = reshape(mean(Mv, 1), nLevels, nF, 2);
statStd = reshape(std(Mv, 0, 1), nLevels, nF, 2);
statMax = reshape(max(Mv, [], 1), nLevels, nF, 2);
statYield = zeros(nLevels, nF, 2, nYield);
sortedMv = sort(Mv, 1, 'ascend');
for y = 1:nYield
    idx = max(1, ceil(yieldLevels(y) / 100 * nValid));
    statYield(:, :, :, y) = reshape(sortedMv(idx, :, :, :), nLevels, nF, 2);
end
yieldTop = yieldLevels(end);

% Nominal tasarım kayması MC fieldleriyle eşleşiyor mu (aynı field sayısı)?
nomMatch = nomNF == nF;
if ~nomMatch
    warning('Nominal dosyada %d, MC dosyalarında %d field var: nominal kayma field bazında eşleştirilmedi.', nomNF, nF);
end

levelNames = strings(nLevels, 1);
for lev = 1:nLevels
    if levelErrors_um(lev) == 0
        levelNames(lev) = "Encoder hatası yok";
    else
        levelNames(lev) = sprintf('Encoder ±%g µm', levelErrors_um(lev));
    end
end

% Özet tablosu (uzun biçim): seviye x field x yöntem
nRow = nLevels * nF * 2;
oEnc = zeros(nRow, 1); oField = strings(nRow, 1); oMethod = strings(nRow, 1);
oMean = zeros(nRow, 1); oStd = zeros(nRow, 1); oMax = zeros(nRow, 1);
oYield = zeros(nRow, nYield); oNom = nan(nRow, 1);
r = 0;
for m = 1:2
    for f = 1:nF
        for lev = 1:nLevels
            r = r + 1;
            oEnc(r) = levelErrors_um(lev);
            oField(r) = fieldLabels(f);
            oMethod(r) = methodNames(m);
            oMean(r) = statMean(lev, f, m);
            oStd(r) = statStd(lev, f, m);
            oMax(r) = statMax(lev, f, m);
            oYield(r, :) = squeeze(statYield(lev, f, m, :))';
            if nomMatch
                oNom(r) = nomMag(f, m);
            end
        end
    end
end
Ozet = table(oMethod, oField, oEnc, oMean, oStd, oMax, oNom, ...
    'VariableNames', {'Yontem', 'Field', 'EncoderHatasi_um', 'Ortalama_piksel', 'Std_piksel', ...
    'Max_piksel', 'NominalTasarim_piksel'});
for y = 1:nYield
    Ozet.(sprintf('Yield%d_piksel', yieldLevels(y))) = oYield(:, y);
end

fprintf('On-axis boresight hatası (sonsuz -> yakın odak), piksel = %.2f um, %d dosya:\n', pixelPitch_um, nValid);
disp(Ozet(oField == fieldLabels(1), :));
fprintf('Yield %%Y: dosyaların %%Y''sinde kayma bu değerin altında. Tüm fieldler Excel''de.\n');

%% Çıktı klasörü: "BORESIGHT ANALİZİ" (char(304) = 'İ')
outFolder = fullfile(folderPath, ['BORESIGHT ANAL' char(304) 'Z' char(304)]);
if ~isfolder(outFolder)
    mkdir(outFolder);
end

%% Grafik sayfaları (A4 dikey, sayfa başına iki grafik)
plotPos = {[0.14 0.58 0.74 0.30], [0.14 0.15 0.74 0.30]};   % üst / alt grafik
colors = lines(max(nLevels, 2));
plotFigs = gobjects(0, 1);
axMag = MAG(valid, :, 1, :);   % on-axis
maxMag = max(axMag(:));
if ~isfinite(maxMag) || maxMag == 0, maxMag = 1; end

% Sayfa: on-axis boresight – kümülatif dağılım (üst), encoder hatasına göre (alt)
fig = newA4Page();
addPageTitle(fig, sprintf('On-axis boresight hatası – sonsuz → yakın odak (%d MC dosyası)', nValid));

ax = axes(fig, 'Position', plotPos{1}, 'FontSize', 9, 'FontName', 'Arial');
hold(ax, 'on');
lineStyles = {'-', '--'};
for m = 1:2
    for lev = 1:nLevels
        sortedMag = sort(MAG(valid, lev, 1, m), 'ascend');
        stairs(ax, [0; sortedMag], [0; (1:nValid)' / nValid * 100], lineStyles{m}, 'LineWidth', 1.4, ...
            'Color', colors(lev, :), 'DisplayName', sprintf('%s – %s', methodNames(m), levelNames(lev)));
    end
end
hold(ax, 'off');
grid(ax, 'on'); box(ax, 'on');
xlim(ax, [0 maxMag * 1.05]);
ylim(ax, [0 100]);
xlabel(ax, 'Boresight hatası (piksel)');
ylabel(ax, 'Dosyaların yüzdesi (%)');
title(ax, 'Kümülatif dağılım (düz: centroid, kesikli: ana ışın)');
legend(ax, 'Location', 'southeast', 'FontSize', 7);

ax = axes(fig, 'Position', plotPos{2}, 'FontSize', 9, 'FontName', 'Arial');
hold(ax, 'on');
yColors = lines(nYield);
for y = 1:nYield
    plot(ax, levelErrors_um, statYield(:, 1, 1, y), '-o', 'Color', yColors(y, :), 'LineWidth', 1.5, ...
        'MarkerFaceColor', yColors(y, :), 'DisplayName', sprintf('Centroid Y%d%%', yieldLevels(y)));
end
plot(ax, levelErrors_um, statYield(:, 1, 2, end), 'k--o', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('Ana ışın Y%d%%', yieldTop));
plot(ax, levelErrors_um, statMax(:, 1, 1), 'k:s', 'LineWidth', 1, 'DisplayName', 'Centroid maksimum');
hold(ax, 'off');
grid(ax, 'on'); box(ax, 'on');
xPad = max(1, 0.05 * (max(levelErrors_um) - min(levelErrors_um)));
xlim(ax, [min(levelErrors_um) - xPad, max(levelErrors_um) + xPad]);
xticks(ax, levelErrors_um);
ylim(ax, [0 maxMag * 1.1]);
xlabel(ax, 'Encoder hatası ±δ (µm)');
ylabel(ax, 'Boresight hatası (piksel)');
title(ax, 'Encoder hatasına göre on-axis boresight hatası');
legend(ax, 'Location', 'northwest', 'FontSize', 8);
plotFigs(end + 1, 1) = fig;

% Sayfa: fieldlere göre yield (üst: centroid, alt: ana ışın)
if nF > 1
    fig = newA4Page();
    addPageTitle(fig, sprintf('Fieldlere göre kayma – yield %d%% (siyah: nominal tasarım)', yieldTop));
    for m = 1:2
        ax = axes(fig, 'Position', plotPos{m}, 'FontSize', 9, 'FontName', 'Arial');
        vals = reshape(statYield(:, :, m, end), nLevels, nF)';   % nF x nLevels
        hb = bar(ax, 1:nF, vals, 'grouped');
        for lev = 1:nLevels
            hb(lev).FaceColor = colors(lev, :);
            hb(lev).DisplayName = levelNames(lev);
        end
        if nomMatch
            hold(ax, 'on');
            plot(ax, 1:nF, nomMag(:, m), 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 6, ...
                'DisplayName', 'Nominal tasarım');
            hold(ax, 'off');
        end
        grid(ax, 'on'); box(ax, 'on');
        xticks(ax, 1:nF);
        xticklabels(ax, cellstr(fieldLabels));
        xtickangle(ax, 30);
        ylabel(ax, sprintf('Y%d%% kayma (piksel)', yieldTop));
        title(ax, char(methodNames(m)));
        legend(ax, 'Location', 'northwest', 'FontSize', 7);
    end
    plotFigs(end + 1, 1) = fig;
end

% Sayfalar: on-axis dx - dy dağılımı (sayfa başına iki seviye, centroid)
scPos = {[0.22 0.55 0.56 0.33], [0.22 0.12 0.56 0.33]};
lim = maxMag * 1.15;
th = linspace(0, 2 * pi, 200);
for lev = 1:nLevels
    slot = mod(lev - 1, 2) + 1;
    if slot == 1
        fig = newA4Page();
        addPageTitle(fig, 'On-axis kayma yönü (dx – dy, centroid): her nokta bir MC dosyası');
        plotFigs(end + 1, 1) = fig; %#ok<SAGROW>
    end
    ax = axes(fig, 'Position', scPos{slot}, 'FontSize', 9, 'FontName', 'Arial');
    plot(ax, DX(valid, lev, 1, 1), DY(valid, lev, 1, 1), '.', 'Color', colors(lev, :), 'MarkerSize', 10);
    hold(ax, 'on');
    rr = statYield(lev, 1, 1, end);
    plot(ax, rr * cos(th), rr * sin(th), 'k--', 'LineWidth', 1);
    plot(ax, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1);
    hold(ax, 'off');
    grid(ax, 'on'); box(ax, 'on');
    xlim(ax, [-lim lim]); ylim(ax, [-lim lim]);
    pbaspect(ax, [1 1 1]);
    xlabel(ax, 'dx (piksel)'); ylabel(ax, 'dy (piksel)');
    title(ax, {char(levelNames(lev)), sprintf('kesikli daire: yield %d%% = %.3g piksel', yieldTop, rr)});
end

%% Excel
outputFile = fullfile(outFolder, 'Boresight_Results.xlsx');
if isfile(outputFile)
    delete(outputFile);
end
writetable(Ozet, outputFile, 'Sheet', 'Ozet');

Nominal = table(nomFieldLabels, nomDX(:, 1), nomDY(:, 1), nomMag(:, 1), nomDX(:, 2), nomDY(:, 2), nomMag(:, 2), ...
    'VariableNames', {'Field', 'Centroid_dx_piksel', 'Centroid_dy_piksel', 'Centroid_piksel', ...
    'AnaIsin_dx_piksel', 'AnaIsin_dy_piksel', 'AnaIsin_piksel'});
writetable(Nominal, outputFile, 'Sheet', 'Nominal');

% Dosyalar: her satır bir dosya x field
nR = nFiles * nF;
Dosyalar = table(repelem(fileNames, nF), repelem(status, nF), repmat(fieldLabels, nFiles, 1), ...
    repelem(compSurfUsed, nF), repelem(compT0, nF), repelem(eflUsed, nF), ...
    'VariableNames', {'Dosya', 'Durum', 'Field', 'Kompanzator_Yuzey', 'Kompanzator_T0_mm', 'EFL_mm'});
methodTags = ["Centroid", "AnaIsin"];
for m = 1:2
    for lev = 1:nLevels
        tag = sprintf('%s_%sum', methodTags(m), strrep(sprintf('%g', levelErrors_um(lev)), '.', 'p'));
        dxv = reshape(permute(DX(:, lev, :, m), [3 1 2 4]), nR, 1);   % dosya içinde field sırası
        dyv = reshape(permute(DY(:, lev, :, m), [3 1 2 4]), nR, 1);
        Dosyalar.([tag '_dx_piksel']) = dxv;
        Dosyalar.([tag '_dy_piksel']) = dyv;
        Dosyalar.([tag '_piksel']) = hypot(dxv, dyv);
        if levelErrors_um(lev) ~= 0
            Dosyalar.([tag '_Isaret']) = reshape(permute(SGN(:, lev, :, m), [3 1 2 4]), nR, 1);
        end
    end
end
writetable(Dosyalar, outputFile, 'Sheet', 'Dosyalar');

Ayarlar = table( ...
    ["Klasor"; "Nominal dosya"; "Sonsuz odak config"; "Yakin odak config"; ...
     "Nominal kompanzator yuzeyi"; "Delta komp (mm)"; "Nominal yakin odak paraxial yuzeyi"; ...
     "Yakin odak focal length (mm)"; "MC kompanzator yuzeyi"; "Paraxial eklendigi yer"; ...
     "|EFL| (mm)"; "Piksel (um)"; "Encoder hatalari (um)"; "Field tipi"; "Field sayisi"; ...
     "Centroid pupil noktasi (her dalga boyu)"; "MC dosya sayisi"; "Gecerli dosya"], ...
    [string(folderPath); string(nominalName); string(infConfig); string(nearConfig); ...
     string(nomCompSurface); sprintf('%.6f', dComp_mm); string(nomNearSurface); ...
     sprintf('%.6g', fNear_mm); settingText(mcCompSurface, "her dosyada variable kalinlik"); ...
     settingText(mcInsertAfter, "goruntuden hemen once"); ...
     sprintf('%.4f - %.4f', min(eflUsed(valid)), max(eflUsed(valid))); string(pixelPitch_um); ...
     string(mat2str(compErrors_um)); fieldTypeName; string(nF); string(numel(pupilPx)); ...
     string(nFiles); string(nValid)], ...
    'VariableNames', {'Ayar', 'Deger'});
writetable(Ayarlar, outputFile, 'Sheet', 'Ayarlar');

save(fullfile(outFolder, 'Boresight_RawData.mat'), 'fileNames', 'status', 'levelErrors_um', ...
    'methodNames', 'fieldLabels', 'DX', 'DY', 'MAG', 'SGN', 'compT0', 'compSurfUsed', 'eflUsed', ...
    'statMean', 'statStd', 'statMax', 'statYield', 'yieldLevels', 'nomFieldLabels', 'nomDX', 'nomDY', ...
    'nomMag', 'dComp_mm', 'fNear_mm', 'efl_mm', 'pixelPitch_um', 'infConfig', 'nearConfig', ...
    'nomCompSurface', 'nomNearSurface', 'mcCompSurface', 'mcInsertAfter', 'nominalName');

%% PDF raporu (A4)
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
eflMin = min(eflUsed(valid));
eflMax = max(eflUsed(valid));
if abs(eflMax - eflMin) < 5e-5
    eflText = sprintf('%.4f mm', eflMin);
else
    eflText = sprintf('%.4f – %.4f mm (her dosyanın kendi EFL''i)', eflMin, eflMax);
end
if isempty(mcInsertAfter)
    insertText = 'görüntü yüzeyinden hemen önce';
else
    insertText = sprintf('yüzey %d''den sonra', mcInsertAfter);
end
nPupil = numel(pupilPx);

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
    sprintf('MC yakın odak lensi:    %s, thickness 0', insertText)
    sprintf('|EFL|:                  %s', eflText)
    sprintf('Piksel boyutu:          %.2f µm', pixelPitch_um)
    sprintf('Encoder hatası:         ±%s µm (kötü taraf)', errListText)
    sprintf('Fieldler:               %d (%s), on-axis 1. sırada', nF, char(fieldTypeName))
    sprintf('Centroid:               %d pupil noktası x tüm dalga boyları (ağırlıklı)', nPupil)
    ''
    };

for m = 1:2
    summaryLines = [summaryLines; {sprintf('ON-AXIS BORESIGHT HATASI – %s (piksel)', upper(char(methodNames(m))))}]; %#ok<AGROW>
    summaryLines = [summaryLines; statTableLines(levelErrors_um, statMean(:, 1, m), statStd(:, 1, m), ...
        reshape(statYield(:, 1, m, :), nLevels, nYield), statMax(:, 1, m), yieldLevels); {''}]; %#ok<AGROW>
end

if nF > 1
    summaryLines = [summaryLines; {sprintf('FIELDLER – YIELD %d%% KAYMA, CENTROID (piksel)', yieldTop)}];
    hdr = sprintf('%-24s', 'Field');
    for lev = 1:nLevels
        hdr = [hdr sprintf('%10s', sprintf('±%g µm', levelErrors_um(lev)))]; %#ok<AGROW>
    end
    hdr = [hdr sprintf('%10s', 'Nominal')];
    summaryLines = [summaryLines; {hdr; repmat('-', 1, numel(hdr))}];
    for f = 1:nF
        row = sprintf('%-24s', char(fieldLabels(f)));
        for lev = 1:nLevels
            row = [row sprintf('%10.3g', statYield(lev, f, 1, end))]; %#ok<AGROW>
        end
        if nomMatch
            row = [row sprintf('%10.3g', nomMag(f, 1))]; %#ok<AGROW>
        else
            row = [row sprintf('%10s', '-')]; %#ok<AGROW>
        end
        summaryLines{end + 1, 1} = row; %#ok<SAGROW>
    end
    summaryLines = [summaryLines; {'Nominal = hatasız tasarımın kendi kayması (odak nefesi). Ana ışın değerleri Excel''de.'; ''}];
end

summaryLines = [summaryLines; {
    sprintf('Y%%: dosyaların %%Y''sinde kayma bu değerin altındadır. Örnek: ±0 µm, centroid Y%d%% = %.3g piksel.', ...
        yieldTop, statYield(1, 1, 1, end))
    sprintf('Nominal tasarım on-axis kayması: centroid %.3g, ana ışın %.3g piksel (≈ 0 beklenir).', nomMag(1, 1), nomMag(1, 2))
    ''
    sprintf('KONTROL (%s)', char(checkFile))
    sprintf('  Kolimasyon hatası: sonsuz %.3g rad, yakın %.3g rad, yakın (f işareti ters) %.3g rad', cInf, cNear, cFlip)
    }];
if signWarning
    summaryLines{end + 1, 1} = ['  UYARI: Focal length işareti ters çevrildiğinde çıkış daha iyi kolime oluyor. ' ...
        'Yakın odak lensinin işaretini veya Δkomp yönünü kontrol edin.'];
else
    summaryLines{end + 1, 1} = '  Yakın odak lensinin işareti ve Δkomp yönü tutarlı.';
end
summaryLines{end + 1, 1} = sprintf('  On-axis kayma (0 µm): centroid (%+.3g, %+.3g), ana ışın (%+.3g, %+.3g) piksel', ...
    checkShift(1, 1), checkShift(2, 1), checkShift(1, 2), checkShift(2, 2));
summaryLines = [summaryLines; diagLines(:)];

if any(~valid)
    summaryLines = [summaryLines; {''; sprintf('Hesaplanamayan dosyalar (%d):', sum(~valid))}];
    shown = find(~valid);
    shown = shown(1:min(end, 6));
    for i = shown'
        summaryLines{end + 1, 1} = sprintf('  %s: %s', fileNames(i), status(i)); %#ok<SAGROW>
    end
    if sum(~valid) > numel(shown)
        summaryLines{end + 1, 1} = sprintf('  ... ve %d dosya daha (Excel, Dosyalar sayfası)', sum(~valid) - numel(shown));
    end
end

summaryFig = textPage('BORESIGHT HATASI – SONSUZ → YAKIN ODAK ANALİZ RAPORU', summaryLines, 'Courier New', 8);

% --- Sayfa: yöntem -------------------------------------------------------
explain1 = {
    'AMAÇ'
    ['Sistem sonsuz odaktan yakın odağa geçtiğinde yalnızca kompanzatör lens eksen boyunca hareket eder. ' ...
     'Hatasız sistemde bu hareket görüş eksenini (boresight) değiştirmez. Üretilmiş bir sistemde ise ' ...
     'lensler decenter/tilt hatalı olduğu için kompanzatörün hareketi ve yakın hedefin paralaksı ' ...
     'görüş eksenini kaydırır: sonsuz odakta merkez pikselin baktığı hedef, yakın odakta başka bir ' ...
     'piksele düşer. Bu analiz, Monte Carlo (MC) dosyalarının her biri için bu kaymayı piksel ' ...
     'cinsinden hesaplar: on-axis (boresight) ve diğer fieldler için, ana ışın ve centroid ile.']
    ''
    'ADIM 1 – NOMİNAL YAKIN ODAK DOSYASI'
    sprintf(['Nominal dosyada (%s) iki konfigürasyon vardır. Yakın odak, görüntü tarafına eklenen bir ' ...
     'paraxial yüzeyle simüle edilir ve konfigürasyonlar arasında yalnızca kompanzatör kalınlığı ' ...
     'değişir. Bu dosyadan şunlar okunur:'], nominalName)
    sprintf('•  Δkomp = kompanzatör kalınlığı (yakın) − (sonsuz) = %+.4f mm (yüzey %d).', dComp_mm, nomCompSurface)
    sprintf('•  Yakın odak lensinin focal length''i f = %.6g mm (yüzey %d, yakın odak konfigürasyonunda).', fNear_mm, nomNearSurface)
    ['•  Nominal tasarımın kendi kayması (odak nefesi): iki konfigürasyonda görüntü yüzeyindeki ışın ' ...
     'konumlarının farkı. On-axis''te ≈ 0 olmalıdır; kenar fieldlerde tasarımdan gelen bir değer olabilir.']
    ''
    'ADIM 2 – SONSUZ ODAK (HER MC DOSYASI)'
    ['MC dosyası afokal modda, kaydedildiği haliyle açılır (toleranslar uygulanmış, kompanzatör sonsuz ' ...
     'odağa optimize edilmiş). Her field için ana ışın ve centroid ışınları izlenir ve ışınların ' ...
     'görüntü uzayındaki yönü (tan θx, tan θy) kaydedilir.']
    ''
    'ADIM 3 – YAKIN ODAK (AYNI MC DOSYASI)'
    sprintf(['Aynı dosyaya, nominal dosyadaki gibi bir yakın odak paraxial yüzeyi eklenir (%s, ' ...
     'thickness 0, f = %.6g mm) ve kompanzatör kalınlığı Δkomp kadar değiştirilir. Böylece hatalar ' ...
     'aynen korunur, yalnızca odak değişir. Işınlar tekrar izlenir. Hiçbir adımda optimizasyon ' ...
     'yapılmaz; hesaptan sonra eklenen yüzey silinir, kompanzatör geri alınır ve diskteki dosyalar ' ...
     'değiştirilmez.'], insertText, fNear_mm)
    ''
    'ADIM 4 – PİKSELE ÇEVİRME'
    '   Kayma [piksel] = (tan θ_yakın − tan θ_sonsuz) × |EFL| / piksel boyutu'
    'x ve y ayrı ayrı hesaplanır; büyüklük √(dx² + dy²). Formülün çıkarılışı hesaplama detaylarında.'
    ''
    'ADIM 5 – ENCODER (KOMPANZATÖR KONUM) HATASI'
    sprintf(['Yakın odakta kompanzatör, mekanizmanın konumlandırma hatası kadar yanlış konumlanabilir. ' ...
     'Her hata payı δ (%s µm) için kompanzatör Δkomp + δ ve Δkomp − δ konumlarına getirilir; her ' ...
     'field ve her yöntem için kaymayı büyüten taraf (worst case) kullanılır.'], errListText)
    ''
    'ADIM 6 – İSTATİSTİK'
    ['Her hata payı, field ve yöntem için tüm geçerli dosyaların kayması istatistiğe sokulur: ' ...
     'ortalama, standart sapma, maksimum ve yield değerleri. Yield %Y, dosyaların %Y''sinin altında ' ...
     'kaldığı kaymadır; "sistemlerin %Y''sinde kayma bu değerden küçüktür" anlamına gelir.']
    };

% --- Sayfa: hesaplama detayları I ----------------------------------------
detail1 = {
    'MCE''DEN OKUNANLAR (NOMİNAL DOSYA)'
    ['MCE''deki tüm satırlar dolaşılır; THIC ve IGNR satırlarının Param1 değeri (yüzey numarası) ve her ' ...
     'konfigürasyondaki değeri (GetOperandCell) okunur.']
    '•  IGNR: bazı konfigürasyonlarda 1, bazılarında 0 olan satır → yakın odak paraxial yüzeyi.'
    '      IGNR = 1 (ignore) olan konfigürasyon → sonsuz odak; IGNR = 0 → yakın odak.'
    '•  THIC: iki konfigürasyonda değeri farklı olan satır → kompanzatör yüzeyi.'
    ['Sonra MCE.SetCurrentConfiguration ile her iki konfigürasyona geçilir ve LDE''den kalınlıklar ' ...
     've paraxial yüzeylerin focal length''i (Par1 sütunu) okunur:']
    '      Δkomp = T_komp(yakın) − T_komp(sonsuz)        f_yakın = Par1(yakın odak paraxial yüzeyi)'
    ''
    'MC DOSYASINDA OKUNANLAR'
    ['•  Kompanzatör: kalınlık hücresinin solve tipi Variable olan tek yüzey. Toleranslama coordinate ' ...
     'break eklediği için numara her dosyada ayrıca bulunur.']
    '•  |EFL|: EFFL operandı (afokal modda sistemin paraxial odak uzaklığı).'
    '•  Dalga boyları ve ağırlıkları (SystemData.Wavelengths); ana dalga boyu IsPrimary.'
    ['•  Fieldler (SystemData.Fields) normalize edilir: radyal normalizasyonda Hx = X / r_max, ' ...
     'Hy = Y / r_max (r_max = en büyük √(X² + Y²)); dikdörtgen normalizasyonda Hx = X / |X|max, ' ...
     'Hy = Y / |Y|max. On-axis (Hx = Hy = 0) her zaman 1. field olarak eklenir.']
    ''
    'IŞIN İZLEME (BATCH RAY TRACE, GERÇEK IŞINLAR)'
    ['Her durum (sonsuz odak, yakın odak ± δ) için tek bir batch ray trace (CreateNormUnpol, ' ...
     'RaysType.Real) görüntü yüzeyine kadar çalıştırılır. Her field için:']
    '•  Ana ışın: Px = Py = 0, ana dalga boyu. Stop''un merkezinden geçen ışındır.'
    sprintf(['•  Centroid ışınları: pupil, çap boyunca %d hücreye bölünür; birim çemberin içinde kalan ' ...
     '%d hücre merkezi (Px, Py) her dalga boyu için izlenir. Hücre merkezleri kullanıldığı için pupil ' ...
     'kenarındaki ışınlar (vignetting) dışarıda kalır.'], pupilGrid, nPupil)
    ['Her ışın için görüntü yüzeyinden sonraki yön kosinüsleri (L, M, N) okunur. Hata kodu ya da ' ...
     'vignette kodu sıfır olmayan ışınlar centroid''e katılmaz.']
    '      tan θx = L / N          tan θy = M / N'
    ['Dedektörden sahneye trace edilen sistemde on-axis nokta merkez pikseldir; onun ana ışınının ' ...
     'çıkış yönü, merkez pikselin sahnede baktığı yöndür (boresight).']
    ''
    'CENTROID'
    '      tan θc = Σ w(λ) · tan θ(ışın) / Σ w(λ)       (tüm geçerli centroid ışınları, tüm dalga boyları)'
    ['w(λ) dalga boyu ağırlığıdır; pupil üniform aydınlatılmış kabul edilir (her ışın eşit alanı ' ...
     'temsil eder). Dedektörde y = f · tan θ olduğundan (sonraki sayfa), açı uzayındaki bu ortalama ' ...
     'doğrudan dedektördeki görüntü ağırlık merkezine karşılık gelir. Centroid, koma ve renk ' ...
     'etkilerini de içerdiği için bir izleyicinin (tracker) ölçeceği konuma ana ışından daha yakındır.']
    };

% --- Sayfa: hesaplama detayları II ---------------------------------------
detail2 = {
    'PİKSELE ÇEVİRME – NEDEN |EFL| İLE ÇARPILIR?'
    ['Nominal dosyadaki gibi, dedektör tarafında odak uzaklığı f olan ideal bir lens ve onun f kadar ' ...
     'gerisinde dedektör düşünülür. Lense y yüksekliğinde, u açısıyla (u = tan θ) gelen bir ışın için:']
    '      Lens sonrası açı:      u'' = u − y / f'
    '      Dedektördeki konum:    y_det = y + f · u'' = y + f · u − y = f · u'
    ['Dedektördeki konum yalnızca açıya bağlıdır, ışının yüksekliğine bağlı değildir. Bu yüzden afokal ' ...
     'moddaki açıyı |EFL| ile çarpmak, f = |EFL| paraxial lensle yapılan hesabın aynısıdır:']
    '      dx = (tan θx_yakın − tan θx_sonsuz) × |EFL| [mm] × 1000 / piksel [µm]'
    '      dy = (tan θy_yakın − tan θy_sonsuz) × |EFL| [mm] × 1000 / piksel [µm]'
    '      Kayma = √(dx² + dy²)'
    sprintf(['Örnek: Δtan θ = 10 µrad, |EFL| = %.1f mm, piksel = %.0f µm → 10e-6 × %.1f × 1000 / %.0f = ' ...
     '%.3f piksel.'], eflMin, pixelPitch_um, eflMin, pixelPitch_um, 10e-6 * eflMin * 1000 / pixelPitch_um)
    ''
    'YAKIN ODAK LENSİNİN ETKİSİ (PARALAKS)'
    ['Eklenen paraxial yüzey ışının açısını u'' = u − y / f_yakın kadar değiştirir. Hatasız sistemde ' ...
     'on-axis ana ışın eksendedir (y = 0) ve etki yoktur. MC dosyasında ana ışın eksenden biraz ' ...
     'kaçık çıktığı için yakın hedefte −y / f_yakın kadar bir paralaks açısı oluşur. Ayrıca ' ...
     'kompanzatörün eksenel hareketi, decenter/tilt hatalı sistemde ışının çıkış yönünü değiştirir. ' ...
     'İki etki de doğrudan ışın izlemeden gelir; formülle yaklaşım yapılmaz.']
    ''
    'NOMİNAL TASARIM KAYMASI (ODAK NEFESİ)'
    ['Nominal yakın odak dosyasında aynı ışınlar her iki konfigürasyonda görüntü yüzeyine kadar ' ...
     'izlenir ve konum farkı (X, Y) / piksel boyutu alınır (dosya afokalse yön farkı × |EFL|). ' ...
     'On-axis''te simetri nedeniyle ≈ 0 çıkmalıdır; kenar fieldlerdeki değer tasarımın kendisinden ' ...
     'gelir (odaklamayla büyütme/distorsiyon değişimi) ve MC sonuçlarına da dahildir.']
    ''
    'ENCODER HATASI'
    ['Her δ için iki durum hesaplanır: T_komp = T0 + Δkomp + δ ve T0 + Δkomp − δ. Her field ve yöntem ' ...
     'için büyüklüğü daha büyük olan taraf alınır; uygulanan işaret Excel''de (_Isaret sütunları). ' ...
     'Sonsuz odakta kompanzatör MC dosyasındaki optimize değerinde (T0) kalır.']
    ''
    'KOLİMASYON KONTROLÜ'
    '      Kolimasyon hatası = | atan(tan θy(Py = +0.7)) − atan(tan θy(Py = −0.7)) |   [rad], on-axis'
    ['Kolime çıkışta ≈ 0, defokusta büyür. Sonsuz odak, yakın odak (f) ve yakın odak (−f) için ' ...
     'hesaplanır; −f, kolimasyon hatasını yarıdan fazla azaltıyorsa (cFlip < 0.5 · cNear) uyarı verilir.']
    ''
    'İSTATİSTİK'
    '      Ortalama = Σ kayma / N        Std = √( Σ (kayma − ortalama)² / (N − 1) )'
    ['      Yield %Y = değerler küçükten büyüğe sıralanır, ceil(Y/100 × N)''inci değer alınır ' ...
     '(ör. N = 200, Y90 → 180. değer). Küçük değer iyidir.']
    ''
    'BİRİMLER'
    ['Lens birimi mm değilse (cm, inç, m) kalınlık, focal length ve EFL mm''ye çevrilir; MC dosyasına ' ...
     'uygulanırken o dosyanın birimine geri çevrilir. Encoder hatası µm olarak girilir.']
    };

% --- Sayfa: kapsam ve yorum ----------------------------------------------
explain2 = {
    'KONTROL: YAKIN ODAK DOĞRU SİMÜLE EDİLDİ Mİ?'
    ['Yakın odakta çıkışın kolimasyon hatası sonsuz odaktakine yakın olmalıdır. Sistem dedektörden ' ...
     'sahneye doğru tasarlandığında, sonlu uzaklıktaki bir sahne noktası çıkışta yakınsayan bir ışın ' ...
     'demektir; bunu kolime eden lensin işareti nominal dosyadaki Δkomp yönüyle tutarlı olmalıdır. ' ...
     'Özet sayfasındaki kontrol satırı bunu gösterir.']
    ''
    'BU ANALİZ NEYİ KAPSAR?'
    ['Tolerans analizinin MC dosyalarına uyguladığı tüm hatalar (lens decenter/tilt, yüzey hataları, ' ...
     'boşluklar vb.) sonsuz ve yakın odakta AYNI kalır. Hesaplanan kayma, bu hatalı sistemde ' ...
     'kompanzatörün eksenel hareketinden ve yakın hedefin paralaksından kaynaklanan görüntü kaymasıdır; ' ...
     'kenar fieldlerde tasarımın kendi odak nefesi de dahildir. Bu kayma her seferinde aynı olduğu ' ...
     'için (sistematik) odak konumuna göre bir kalibrasyon tablosuyla yazılımda düzeltilebilir.']
    ''
    'BU ANALİZ NEYİ KAPSAMAZ?'
    ['•  Kızak hatası: Lensin kızak boyunca hareket ederken farklı decenter/tilt alması (kızağın optik ' ...
     'eksene göre açısı ve düzlüğü). MC dosyalarında kompanzatörün hatası her iki odakta aynıdır.']
    ['•  Gidiş-dönüş (tekrarlanabilirlik) hatası: Lens yakın odaktan sonsuza döndüğünde ilk konumuna ' ...
     'tam oturmaması (boşluk, histerezis).']
    ['•  Piksele çevirmede ideal (distorsiyonsuz) dedektör eşlemesi y = |EFL| · tan θ kullanılır. ' ...
     'Kenar fieldlerde distorsiyon varsa yerel ölçek biraz farklıdır; küçük kaymalar için etkisi ' ...
     'ihmal edilebilir.']
    ''
    'ANA IŞIN MI, CENTROID Mİ?'
    ['Ana ışın tek bir ışındır ve geometrik görüş eksenini verir. Centroid ise pupil boyunca ve tüm ' ...
     'dalga boylarında görüntünün ağırlık merkezidir; koma ve lateral renk sapması varsa ana ışından ' ...
     'farklıdır. Hedef takibi ve boresight ölçümü görüntünün ağırlık merkezine göre yapıldığından ' ...
     'spesifikasyon için centroid önerilir. İkisi arasındaki büyük fark, hatalı sistemde koma veya ' ...
     'renk sapmasının belirgin olduğunu gösterir.']
    ''
    'GRAFİKLER NASIL OKUNUR?'
    ['•  Kümülatif dağılım: yatay eksende on-axis kayma, dikey eksende bu değerin altında kalan ' ...
     'dosyaların yüzdesi. Eğri sola ne kadar yakınsa sistem o kadar iyidir.']
    ['•  Encoder hatasına göre: yield seviyelerinin hata payıyla değişimi. Eğri yataysa encoder ' ...
     'hatası boresight''ı pratikte etkilemiyordur.']
    sprintf(['•  Fieldlere göre: her field için yield %%%d kayma; siyah işaret nominal tasarımın ' ...
     'kendi kayması. Çubuk ile işaret arasındaki fark toleranslardan gelir.'], yieldTop)
    sprintf(['•  dx – dy dağılımı: her nokta bir MC dosyasının on-axis kayma vektörü. Tek yönde ' ...
     'toplanma sistematik bir hatayı gösterir. Kesikli daire yield %%%d değeridir.'], yieldTop)
    ''
    'SPESİFİKASYON NASIL KONUR?'
    '1.  Kabul edilebilir boresight kaymasını piksel cinsinden belirle (ör. 0.5 piksel).'
    ['2.  Özet tablosunda ilgili yield seviyesindeki (ör. Y90%) centroid değerini bu sınırla ' ...
     'karşılaştır. Encoder satırları mekanizma hassasiyetinin etkisini gösterir.']
    ['3.  Sınır aşılıyorsa: kalibrasyonla düzeltme (sistematik kısım), hareketli lensin ve arkasındaki ' ...
     'lenslerin decenter/tilt toleranslarının sıkılaştırılması ya da mekanizmanın iyileştirilmesi ' ...
     'değerlendirilir.']
    };

explainFig1 = textPage('1. ANALİZ YÖNTEMİ', explain1, 'Arial', 9.5);
detailFig1 = textPage('2. HESAPLAMA DETAYLARI – OKUNAN DEĞERLER VE IŞIN İZLEME', detail1, 'Arial', 9);
detailFig2 = textPage('3. HESAPLAMA DETAYLARI – FORMÜLLER VE KONTROLLER', detail2, 'Arial', 9);
explainFig2 = textPage('4. KAPSAM VE YORUM', explain2, 'Arial', 9.5);

% --- Sayfaları birleştir: çerçeve, başlık ve sayfa numaraları -------------
allFigs = [summaryFig; explainFig1; detailFig1; detailFig2; explainFig2; plotFigs];
nPages = numel(allFigs);
headerText = sprintf('Boresight Hatası – Sonsuz → Yakın Odak  |  %s', nominalName);
for p = 1:nPages
    addPageFrame(allFigs(p), headerText, p, nPages);
end

% Grafik sayfalarının PNG'leri
oldPng = dir(fullfile(outFolder, 'Boresight_Grafik_*.png'));
for i = 1:numel(oldPng)
    delete(fullfile(outFolder, oldPng(i).name));
end
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

function lines = statTableLines(levels, mn, sd, yl, mx, yieldLevels)
% Özet sayfası için sabit genişlikli tablo satırları (seviye başına bir satır)
header = sprintf('%-11s%10s%10s', 'Hata (µm)', 'Ortalama', 'Std');
for y = 1:numel(yieldLevels)
    header = [header sprintf('%10s', sprintf('Y%d%%', yieldLevels(y)))]; %#ok<AGROW>
end
header = [header sprintf('%10s', 'Max')];
lines = {header; repmat('-', 1, numel(header))};
for lev = 1:numel(levels)
    row = sprintf('%-11s%10.3g%10.3g', sprintf('±%g', levels(lev)), mn(lev), sd(lev));
    for y = 1:numel(yieldLevels)
        row = [row sprintf('%10.3g', yl(lev, y))]; %#ok<AGROW>
    end
    lines{end + 1, 1} = [row sprintf('%10.3g', mx(lev))]; %#ok<AGROW>
end
end

function R = traceFields(TheSystem, surf, hx, hy, chiefWave, waves, weights, px, py, useAngle)
% Her field için ana ışın ve centroid (batch ray trace, gerçek ışınlar).
% R(f, :) = [ana ışın x, ana ışın y, centroid x, centroid y]
% useAngle = true : x, y = tan θx = L/N, tan θy = M/N (afokal)
% useAngle = false: x, y = görüntü yüzeyindeki konum X, Y (lens birimi)
nF = numel(hx);
nP = numel(px);
nW = numel(waves);
nRays = nF * (1 + nP * nW);
opd = ZOSAPI.Tools.RayTrace.OPDMode.None;

rt = TheSystem.Tools.OpenBatchRayTrace();
if isempty(rt)
    error('Batch ray trace açılamadı (başka bir araç açık olabilir).');
end
try
    normUnpol = rt.CreateNormUnpol(nRays, ZOSAPI.Tools.RayTrace.RaysType.Real, surf);
    for f = 1:nF
        normUnpol.AddRay(chiefWave, hx(f), hy(f), 0, 0, opd);
        for w = 1:nW
            for p = 1:nP
                normUnpol.AddRay(waves(w), hx(f), hy(f), px(p), py(p), opd);
            end
        end
    end
    rt.RunAndWaitForCompletion();
    normUnpol.StartReadingResults();

    R = nan(nF, 4);
    for f = 1:nF
        [ok, v] = readRay(normUnpol, useAngle, false);
        if ok
            R(f, 1:2) = v;
        end
        sx = 0; sy = 0; sw = 0;
        for w = 1:nW
            for p = 1:nP
                [ok, v] = readRay(normUnpol, useAngle, true);
                if ok
                    sx = sx + weights(w) * v(1);
                    sy = sy + weights(w) * v(2);
                    sw = sw + weights(w);
                end
            end
        end
        if sw > 0
            R(f, 3:4) = [sx sy] / sw;
        end
    end
catch ME
    rt.Close();
    rethrow(ME);
end
rt.Close();
end

function [ok, v] = readRay(normUnpol, useAngle, rejectVignetted)
[success, ~, errCode, vigCode, X, Y, ~, L, M, N] = normUnpol.ReadNextResult();
ok = logical(success) && errCode == 0 && (~rejectVignetted || vigCode == 0);
v = [NaN NaN];
if ~ok
    return;
end
if useAngle
    if abs(N) < 1e-12
        ok = false;
        return;
    end
    v = [L / N, M / N];
else
    v = [X, Y];
end
end

function [hx, hy, labels] = fieldList(TheSystem, fieldsToUse)
% Normalize field koordinatları. On-axis (0, 0) her zaman 1. sıradadır;
% dosyadaki (0, 0) fieldi tekrar eklenmez.
F = TheSystem.SystemData.Fields;
n = F.NumberOfFields;
X = zeros(n, 1);
Y = zeros(n, 1);
for i = 1:n
    fld = F.GetField(i);
    X(i) = fld.X;
    Y(i) = fld.Y;
end
rect = false;
try
    rect = contains(lower(string(F.Normalization.ToString())), "rect");
catch
end
if rect
    mx = max(abs(X)); my = max(abs(Y));
    HX = X / max(mx, eps) .* (mx > 0);
    HY = Y / max(my, eps) .* (my > 0);
else
    rMax = max(hypot(X, Y));
    HX = X / max(rMax, eps) .* (rMax > 0);
    HY = Y / max(rMax, eps) .* (rMax > 0);
end

if isempty(fieldsToUse)
    sel = 1:n;
else
    sel = fieldsToUse(:)';
    if any(sel < 1 | sel > n | sel ~= round(sel))
        error('fieldsToUse geçersiz: dosyada %d field var.', n);
    end
end

hx = 0; hy = 0; labels = "On-axis (0, 0)";
for i = sel
    if X(i) == 0 && Y(i) == 0
        continue;
    end
    hx(end + 1, 1) = HX(i); %#ok<AGROW>
    hy(end + 1, 1) = HY(i); %#ok<AGROW>
    labels(end + 1, 1) = sprintf('F%d (%g, %g)', i, X(i), Y(i)); %#ok<AGROW>
end
end

function [idx, w, primary] = waveList(TheSystem)
% Ağırlığı sıfır olmayan dalga boyları, ağırlıkları ve ana dalga boyu
wls = TheSystem.SystemData.Wavelengths;
n = wls.NumberOfWavelengths;
idx = (1:n)';
w = zeros(n, 1);
primary = 1;
for i = 1:n
    wl = wls.GetWavelength(i);
    w(i) = wl.Weight;
    try
        if wl.IsPrimary
            primary = i;
        end
    catch
    end
end
keep = w > 0;
idx = idx(keep);
w = w(keep);
if isempty(idx)
    idx = primary;
    w = 1;
end
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
