%% HIZALAMA TOLERANS ANALIZI (ZOS-API, MATLAB)
% Amac: hangi lensin hizalamasi (decenter/tilt) MTF yield'ini ne kadar
% etkiliyor, toleransi sikilastirinca yield ne kadar artiyor?
%
% Kod adim adim yaziliyor. Su an sadece ADIM 1 var.
%
%  ADIM 1 - OKU VE GOSTER (dosyada hicbir sey degistirilmez)
%    .zmx dosyasi acilir ve su bilgiler cekilir:
%      - Sistem: lens birimi, afokal mod, dalga boylari, alanlar
%      - LDE   : her yuzey (tip, yorum, R, kalinlik, malzeme, solve tipi,
%                coordinate break parametreleri, global z konumu)
%      - TDE   : her satir (operand, parametreler, min/max, hangi STAT'in
%                altinda oldugu), satirin sinifi (Uretim / Hizalama /
%                Bosluk / Kompanzator / Kontrol) ve okunur birimi
%                (decenter um, tilt arcmin)
%      - Hizalama hiyerarsisi: TEDX/TEDY/TETX/TETY satirlarinin yuzey
%                araliklari karsilastirilir. Baska araligi icine alan aralik
%                "Grup" (hucre/govde), almayan "Eleman" sayilir.
%    Supheli girisler icin uyari verilir (ornegin 0.5 dereceden buyuk tilt).
%    Ciktilar: komut penceresi + <cikti>_adim1.xlsx + <cikti>_adim1.mat
%
% Notlar
%  - Zemax'te TETX/TETY/TSTX/TSTY DERECE, TEDX/TEDY/TSDX/TSDY/TIRX/TIRY
%    LENS BIRIMI (mm) cinsindendir.
%  - TETX/TETY elemani, araligin ILK yuzeyinin tepe noktasi etrafinda dondurur
%    (Zemax varsayilani). Pivot sutunu bu yuzeyi gosterir.
%  - Turkce karakter kullanilmadi (eski MATLAB surumlerinde kodlama sorunu olmasin diye).

clc; clear; close all;

%% ===================== KULLANICI AYARLARI =====================
ayar.dosya = 'C:\...\tasarim.zmx';   % TDE'si dolu nominal tasarim
ayar.cikti = 'hizalama';             % cikti dosyalarinin on adi

% Uyari esikleri (sadece kontrol icin, analiz degil)
ayar.uyari.tiltDeg   = 0.5;   % [derece] bundan buyuk tilt -> "arcmin yerine derece mi girildi?"
ayar.uyari.decMm     = 0.5;   % [mm] bundan buyuk decenter -> supheli

% "Yol olarak kopyala" ile gelen cift tirnaklari temizle
ayar.dosya = strtrim(strrep(ayar.dosya, '"', ''));
if ~isfile(ayar.dosya), error('Dosya bulunamadi: %s', ayar.dosya); end

%% ===================== BAGLANTI =====================
TheApplication = baglan();
TheSystem = TheApplication.PrimarySystem;
if ~TheSystem.LoadFile(ayar.dosya, false)
    error('Dosya acilamadi: %s', ayar.dosya);
end

%% ===================== ADIM 1: OKU VE GOSTER =====================
fprintf('\n================ ADIM 1: DOSYA OKUMA ================\n');
fprintf('Dosya: %s\n', ayar.dosya);

% ---------- 1a. Sistem ----------
SD = TheSystem.SystemData;
birimAd = char(SD.Units.LensUnits.ToString());
umCarpan = birimdenUm(birimAd);          % 1 lens birimi kac um
afokal = false;
try, afokal = logical(SD.Aperture.AFocalImageSpace); catch, end

nDalga = SD.Wavelengths.NumberOfWavelengths;
dalga = nan(nDalga,1); dalgaW = nan(nDalga,1);
for i = 1:nDalga
    w = SD.Wavelengths.GetWavelength(i);
    dalga(i) = w.Wavelength; dalgaW(i) = w.Weight;
end

F = SD.Fields;
alanTip = char(F.GetFieldType().ToString());
nAlan = F.NumberOfFields;
alanX = nan(nAlan,1); alanY = nan(nAlan,1); alanW = nan(nAlan,1);
for i = 1:nAlan
    f = F.GetField(i);
    alanX(i) = f.X; alanY(i) = f.Y; alanW(i) = f.Weight;
end

fprintf('\n--- SISTEM ---\n');
fprintf('Lens birimi      : %s (1 birim = %g um)\n', birimAd, umCarpan);
fprintf('Afocal Image Space: %s\n', evetHayir(afokal));
fprintf('Dalga boylari    : %s um\n', strjoin(compose('%.3f', dalga'), ', '));
fprintf('Alan tipi        : %s | Alan sayisi: %d\n', alanTip, nAlan);
Talan = table((1:nAlan)', alanX, alanY, alanW, 'VariableNames', {'Alan', 'X', 'Y', 'Agirlik'});
disp(Talan);

uyarilar = {};
if all(alanX == 0) && all(alanY >= 0)
    uyarilar{end+1} = ['Alanlarin hepsi +y ekseninde. Decenter/tilt simetriyi bozar: ' ...
        'analizde -y ve +/-x alanlarini da eklemek gerekebilir.'];
end

% ---------- 1b. LDE ----------
LDE = TheSystem.LDE;
nYuzey = LDE.NumberOfSurfaces;
Y = struct('No', {}, 'Tip', {}, 'Yorum', {}, 'R', {}, 'Kalinlik', {}, 'KalinlikSolve', {}, ...
           'Malzeme', {}, 'YariCap', {}, 'CB', {}, 'GlobalZ', {});
for i = 0:nYuzey-1
    s = LDE.GetSurfaceAt(i);
    y.No       = i;
    y.Tip      = char(s.TypeName);
    y.Yorum    = strtrim(char(s.Comment));
    y.R        = s.Radius;
    y.Kalinlik = s.Thickness;
    y.KalinlikSolve = char(s.ThicknessCell.GetSolveData().Type.ToString());
    y.Malzeme  = strtrim(char(s.Material));
    y.YariCap  = s.SemiDiameter;
    y.CB       = '';
    if contains(lower(y.Tip), 'coordinate')
        p = nan(1,6);
        for k = 1:6
            p(k) = hucreSayi(s.GetSurfaceCell(parSutun(k)));   % Par6 (sira) tamsayi hucre
        end
        y.CB = sprintf('dX=%g dY=%g tX=%g tY=%g tZ=%g sira=%g', p);
    end
    y.GlobalZ = globalZ(LDE, i);
    Y(end+1) = orderfields(y, Y); %#ok<SAGROW>
end
Tlde = struct2table(Y, 'AsArray', true);
fprintf('\n--- LDE (%d yuzey) ---\n', nYuzey);
disp(Tlde(:, {'No', 'Tip', 'Yorum', 'R', 'Kalinlik', 'KalinlikSolve', 'Malzeme', 'GlobalZ', 'CB'}));

% ---------- 1c. TDE ----------
TDE = TheSystem.TDE;
nOp = TDE.NumberOfOperands;
if nOp == 0, error('TDE bos. Dosyayi TDE dolu haliyle .zmx olarak kaydet.'); end

statAd = 'Normal 2sigma (varsayilan)';
R = struct('Satir', {}, 'Operand', {}, 'P1', {}, 'P2', {}, 'P3', {}, 'P4', {}, ...
           'Min', {}, 'Max', {}, 'Birim', {}, 'Okunur', {}, 'Sinif', {}, ...
           'AitOldugu', {}, 'Dagilim', {}, 'Yorum', {});
okuYolu = ''; okuHata = ''; dokumYapildi = false;
for i = 1:nOp
    op = TDE.GetOperandAt(i);
    r.Satir   = i;
    r.Operand = upper(strtrim(char(op.Type.ToString())));
    prm = nan(1,4);
    for k = 1:4
        [prm(k), yol, hata] = parametreOku(op, k);
        if ~isempty(yol) && isempty(okuYolu), okuYolu = yol; end
        if ~isempty(hata), okuHata = hata; end
    end
    % Yuzey parametresi olmasi gereken bir satirda Param1 okunamadiysa teshis yaz (bir kez)
    if isnan(prm(1)) && ~dokumYapildi && ~any(strcmp(r.Operand, {'STAT','TWAV','SAVE','SEED','TOFF','CEXP'}))
        fprintf('\n! Satir %d (%s): parametreler okunamadi. Son hata: %s\n', i, r.Operand, okuHata);
        tdeSatirDokum(op);
        dokumYapildi = true;
    end
    r.P1 = prm(1); r.P2 = prm(2); r.P3 = prm(3); r.P4 = prm(4);
    [r.Min, r.Max] = minMaxOku(op);
    r.Yorum = '';
    try, r.Yorum = strtrim(char(op.Comment)); catch, end

    if strcmp(r.Operand, 'STAT')
        statAd = statMetin(prm(1), prm(2));
    end
    [r.Sinif, tur] = siniflandir(r.Operand, prm, Y);
    [r.Birim, r.Okunur] = okunurBirim(tur, r.Min, r.Max, birimAd, umCarpan);
    r.AitOldugu = aitOldugu(r.Operand, tur, prm, Y);
    if any(strcmp(r.Sinif, {'Kontrol', 'Kompanzator'}))
        r.Dagilim = '-';
    else
        r.Dagilim = statAd;
    end
    R(end+1) = orderfields(r, R); %#ok<SAGROW>

    % --- Birim kontrolleri ---
    buyuk = max(abs([r.Min r.Max]));
    if strcmp(tur, 'aci') && buyuk > ayar.uyari.tiltDeg
        uyarilar{end+1} = sprintf(['Satir %d %s: %.4g derece buyuk. Arcmin yerine derece mi girildi? ' ...
            '(1 arcmin = 0.0167 derece)'], i, r.Operand, buyuk); %#ok<SAGROW>
    end
    if any(strcmp(r.Operand, {'TEDX', 'TEDY', 'TSDX', 'TSDY'})) && buyuk > ayar.uyari.decMm
        uyarilar{end+1} = sprintf('Satir %d %s: %.4g %s buyuk. Kontrol et.', i, r.Operand, buyuk, birimAd); %#ok<SAGROW>
    end
    if ~any(strcmp(r.Sinif, {'Kontrol', 'Kompanzator'})) && ~isnan(r.Min) && ~isnan(r.Max) && ...
            abs(r.Min + r.Max) > 1e-12 * max(1, buyuk)
        uyarilar{end+1} = sprintf('Satir %d %s: min (%g) ve max (%g) simetrik degil.', i, r.Operand, r.Min, r.Max); %#ok<SAGROW>
    end
end
Ttde = struct2table(R, 'AsArray', true);
if isempty(okuYolu)
    uyarilar{end+1} = ['TDE parametreleri (yuzey numaralari) hicbir yoldan okunamadi. ' ...
        'Yukaridaki teshis ciktisini gonder.'];
else
    fprintf('TDE parametreleri su yolla okundu: %s\n', okuYolu);
end
fprintf('\n--- TDE (%d satir) ---\n', nOp);
disp(Ttde(:, {'Satir', 'Operand', 'P1', 'P2', 'P3', 'Min', 'Max', 'Okunur', 'Sinif', 'AitOldugu', 'Dagilim'}));

% Sinif ozeti
[sinifAd, ~, idx] = unique(Ttde.Sinif);
fprintf('Sinif ozeti: ');
for k = 1:numel(sinifAd)
    fprintf('%s = %d   ', sinifAd{k}, sum(idx == k));
end
fprintf('\n');

% ---------- 1d. Kompanzator ----------
fprintf('\n--- KOMPANZATOR ---\n');
kRows = find(strcmp(Ttde.Operand, 'COMP'));
if isempty(kRows)
    fprintf('TDE''de COMP satiri yok.\n');
    uyarilar{end+1} = 'TDE''de COMP satiri yok: analiz kompanzatorsuz kosulur.';
end
for k = kRows'
    sNo = Ttde.P1(k);
    yy = Y([Y.No] == sNo);
    if isempty(yy), fprintf('Satir %d: yuzey %g bulunamadi\n', k, sNo); continue; end
    fprintf('Satir %d: yuzey %d (%s) | kod %g | kalinlik %.4f %s | strok %+g ... %+g\n', ...
        k, sNo, bosIse(yy.Yorum, '-'), Ttde.P2(k), yy.Kalinlik, birimAd, Ttde.Min(k), Ttde.Max(k));
end

% ---------- 1e. Hizalama hiyerarsisi ----------
hz = find(strcmp(Ttde.Sinif, 'Hizalama'));
fprintf('\n--- HIZALAMA HIYERARSISI ---\n');
if isempty(hz)
    fprintf('Hizalama satiri (TEDX/TEDY/TETX/TETY) bulunamadi.\n');
    Thiz = table();
else
    nanSatir = hz(isnan(Ttde.P1(hz)) | isnan(Ttde.P2(hz)));
    if ~isempty(nanSatir)
        uyarilar{end+1} = sprintf(['%d hizalama satirinda yuzey araligi okunamadi (NaN), ' ...
            'hiyerarsiye katilmadi: satir %s'], numel(nanSatir), strjoin(compose('%d', nanSatir'), ','));
        hz = setdiff(hz, nanSatir);
    end
end
if isempty(hz)
    if exist('nanSatir', 'var') && ~isempty(nanSatir)
        fprintf('Hizalama satirlarinin yuzey araliklari okunamadi, hiyerarsi kurulamadi.\n');
    end
    Thiz = table();
else
    aralik = unique([Ttde.P1(hz) Ttde.P2(hz)], 'rows');
    nA = size(aralik,1);
    derinlik = zeros(nA,1); altSayi = zeros(nA,1);
    for a = 1:nA
        for b = 1:nA
            if a == b, continue; end
            icinde = aralik(b,1) >= aralik(a,1) && aralik(b,2) <= aralik(a,2);
            kesisir = aralik(b,1) <= aralik(a,2) && aralik(b,2) >= aralik(a,1);
            if icinde
                altSayi(a) = altSayi(a) + 1;
                derinlik(b) = derinlik(b) + 1;
            elseif kesisir && a < b && ~(aralik(a,1) >= aralik(b,1) && aralik(a,2) <= aralik(b,2))
                uyarilar{end+1} = sprintf(['Hizalama araliklari S%d-S%d ve S%d-S%d kismen kesisiyor ' ...
                    '(biri digerini tam icermiyor). Hiyerarsi belirsiz.'], aralik(a,:), aralik(b,:)); %#ok<SAGROW>
            end
        end
    end

    H = struct('Seviye', {}, 'Ad', {}, 'Yuzeyler', {}, 'Tip', {}, 'Decenter', {}, 'Tilt', {}, ...
               'Pivot', {}, 'PivotZ', {}, 'Satirlar', {});
    [~, sira] = sortrows([aralik(:,1) -aralik(:,2)]);
    for a = sira'
        h.Seviye = derinlik(a);
        h.Ad = aralikAdi(aralik(a,:), Y);
        h.Yuzeyler = sprintf('S%d-S%d', aralik(a,:));
        if altSayi(a) > 0, h.Tip = sprintf('Grup (%d alt aralik)', altSayi(a)); else, h.Tip = 'Eleman'; end
        bu = hz(Ttde.P1(hz) == aralik(a,1) & Ttde.P2(hz) == aralik(a,2));
        h.Decenter = ozetle(Ttde, bu, {'TEDX', 'TEDY'});
        h.Tilt     = ozetle(Ttde, bu, {'TETX', 'TETY', 'TETZ'});
        yy = Y([Y.No] == aralik(a,1));
        if isempty(yy)
            h.Pivot = sprintf('S%d (LDE''de yok)', aralik(a,1));
            h.PivotZ = NaN;
        else
            h.Pivot = sprintf('S%d %s', aralik(a,1), bosIse(yy.Yorum, ''));
            h.PivotZ = yy.GlobalZ;
        end
        h.Satirlar = strjoin(compose('%d', bu'), ',');
        H(end+1) = orderfields(h, H); %#ok<SAGROW>

        fprintf('%s%s [%s] %s | dec: %s | tilt: %s | pivot: %s\n', repmat('    ', 1, h.Seviye), ...
            h.Ad, h.Yuzeyler, h.Tip, h.Decenter, h.Tilt, h.Pivot);
    end
    Thiz = struct2table(H, 'AsArray', true);
end

% ---------- 1f. Uyarilar ----------
fprintf('\n--- UYARILAR ---\n');
if isempty(uyarilar)
    fprintf('Uyari yok.\n');
else
    uyarilar = unique(uyarilar, 'stable');
    for k = 1:numel(uyarilar), fprintf('! %s\n', uyarilar{k}); end
end

% ---------- 1g. Kaydet ----------
xlsx = [ayar.cikti '_adim1.xlsx'];
if isfile(xlsx), delete(xlsx); end
Tsis = table({'Dosya'; 'Lens birimi'; 'Afocal Image Space'; 'Dalga boylari [um]'; 'Alan tipi'}, ...
    {ayar.dosya; birimAd; evetHayir(afokal); strjoin(compose('%.4f', dalga'), ', '); alanTip}, ...
    'VariableNames', {'Bilgi', 'Deger'});
writetable(Tsis, xlsx, 'Sheet', 'Sistem');
writetable(Talan, xlsx, 'Sheet', 'Alanlar');
writetable(Tlde, xlsx, 'Sheet', 'LDE');
writetable(Ttde, xlsx, 'Sheet', 'TDE');
if ~isempty(Thiz), writetable(Thiz, xlsx, 'Sheet', 'Hizalama'); end
writetable(table(uyarilar(:), 'VariableNames', {'Uyari'}), xlsx, 'Sheet', 'Uyarilar');

yapi.dosya = ayar.dosya; yapi.birim = birimAd; yapi.umCarpan = umCarpan; yapi.afokal = afokal;
yapi.dalga = dalga; yapi.alanlar = Talan; yapi.lde = Tlde; yapi.tde = Ttde; yapi.hizalama = Thiz;
yapi.uyarilar = uyarilar;
save([ayar.cikti '_adim1.mat'], 'yapi');
fprintf('\nKaydedildi: %s ve %s_adim1.mat\n', xlsx, ayar.cikti);
fprintf('Dosyada hicbir degisiklik yapilmadi.\n');

%% ===================== YARDIMCI FONKSIYONLAR =====================
function app = baglan()
    import System.Reflection.*;
    zemaxData = winqueryreg('HKEY_CURRENT_USER', 'Software\Zemax', 'ZemaxRoot');
    NetHelper = strcat(zemaxData, '\ZOS-API\Libraries\ZOSAPI_NetHelper.dll');
    NET.addAssembly(NetHelper);
    if ~ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize()
    % Birden fazla surum kuruluysa ustteki satir yerine:
    % if ~ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize('C:\Program Files\Ansys Zemax OpticStudio 2025 R2.04')
        error('OpticStudio kurulumu bulunamadi.');
    end
    disp(['OpticStudio: ' char(ZOSAPI_NetHelper.ZOSAPI_Initializer.GetZemaxDirectory())]);
    NET.addAssembly(AssemblyName('ZOSAPI_Interfaces'));
    NET.addAssembly(AssemblyName('ZOSAPI'));
    conn = ZOSAPI.ZOSAPI_Connection();
    conn.ConnectionTimeoutSeconds = 20;
    app = conn.CreateNewApplication();
    if isempty(app), error('OpticStudio''ya baglanilamadi.'); end
    if ~app.IsValidLicenseForAPI, error('Lisans kontrolu basarisiz.'); end
end

function c = birimdenUm(ad)
    ad = lower(ad);
    if contains(ad, 'milli'),      c = 1e3;
    elseif contains(ad, 'centi'),  c = 1e4;
    elseif contains(ad, 'inch'),   c = 25400;
    elseif contains(ad, 'meter'),  c = 1e6;
    else,                          c = NaN;
    end
end

function t = evetHayir(b)
    if b, t = 'Evet'; else, t = 'Hayir'; end
end

function t = bosIse(s, yedek)
    if isempty(s), t = yedek; else, t = s; end
end

function c = parSutun(k)
    switch k
        case 1, c = ZOSAPI.Editors.LDE.SurfaceColumn.Par1;
        case 2, c = ZOSAPI.Editors.LDE.SurfaceColumn.Par2;
        case 3, c = ZOSAPI.Editors.LDE.SurfaceColumn.Par3;
        case 4, c = ZOSAPI.Editors.LDE.SurfaceColumn.Par4;
        case 5, c = ZOSAPI.Editors.LDE.SurfaceColumn.Par5;
        case 6, c = ZOSAPI.Editors.LDE.SurfaceColumn.Par6;
    end
end

function v = hucreSayi(h)
    % Hucre tipi Double ya da Integer olabilir: once metin degerinden oku
    v = NaN;
    try, v = str2double(char(h.Value)); catch, end
    if isnan(v)
        try, v = double(h.IntegerValue); catch
            try, v = h.DoubleValue; catch, end
        end
    end
end

function z = globalZ(LDE, i)
    % Yuzeyin global koordinattaki z konumu (pivot yerini gormek icin)
    z = NaN;
    try
        [ok, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, Z] = LDE.GetGlobalMatrix(i);
        if ok, z = Z; end
    catch
    end
end

function [v, yol, hata] = parametreOku(op, k)
    % TDE parametresini (Surf1, Surf2, Code...) sirayla uc yoldan okumayi dener:
    %   1) op.Param1..Param4 ozelligi
    %   2) op.GetOperandCell(TolerancingColumn.ParamN)
    %   3) op.GetCellAt(k)  (sutun 0 = Type, 1.. = parametreler)
    % yol : hangi yolun calistigi, hata: hicbiri calismazsa son hata mesaji
    v = NaN; yol = ''; hata = '';
    ad = sprintf('Param%d', k);
    try
        v = double(op.(ad));
        if ~isnan(v), yol = ['op.' ad]; return; end
    catch e, hata = e.message;
    end
    try
        c = ZOSAPI.Editors.TDE.TolerancingColumn.(ad);
        v = hucreSayi(op.GetOperandCell(c));
        if ~isnan(v), yol = ['GetOperandCell(' ad ')']; return; end
    catch e, hata = e.message;
    end
    try
        v = hucreSayi(op.GetCellAt(k));
        if ~isnan(v), yol = sprintf('GetCellAt(%d)', k); return; end
    catch e, hata = e.message;
    end
    v = NaN;
end

function tdeSatirDokum(op)
    % Parametre okunamazsa teshis icin satirin hucrelerini (baslik = deger) yazar
    fprintf('  Teshis - TDE satiri hucreleri:\n');
    for c = 0:9
        try
            h = op.GetCellAt(c);
            bas = ''; try, bas = char(h.Header); catch, end
            fprintf('    GetCellAt(%d): baslik="%s" deger="%s" tip=%s\n', c, bas, ...
                char(h.Value), char(h.DataType.ToString()));
        catch e
            fprintf('    GetCellAt(%d): %s\n', c, e.message);
            break;
        end
    end
end

function [mn, mx] = minMaxOku(op)
    mn = NaN; mx = NaN;
    try
        mn = op.Min; mx = op.Max;
    catch
        try
            mn = op.GetOperandCell(ZOSAPI.Editors.TDE.TolerancingColumn.Min).DoubleValue;
            mx = op.GetOperandCell(ZOSAPI.Editors.TDE.TolerancingColumn.Max).DoubleValue;
        catch
        end
    end
end

function t = statMetin(tip, sigma)
    switch tip
        case 0, t = sprintf('Normal %gsigma', sigma);
        case 1, t = 'Uniform';
        case 2, t = 'Parabolic';
        case 3, t = 'Kullanici tanimli';
        otherwise, t = sprintf('STAT %g', tip);
    end
end

function [sinif, tur] = siniflandir(op, prm, Y)
    % tur: okunur birimi belirler ('uzunluk', 'aci', 'fringe', 'diger')
    uzunluk = {'TEDX','TEDY','TSDX','TSDY','TIRX','TIRY','TTHI','TRAD','TUDX','TUDY'};
    aci     = {'TETX','TETY','TETZ','TSTX','TSTY','TUTX','TUTY','TUTZ'};
    fringe  = {'TFRN','TIRR','TEZI','TEXI'};
    if any(strcmp(op, uzunluk)), tur = 'uzunluk';
    elseif any(strcmp(op, aci)), tur = 'aci';
    elseif any(strcmp(op, fringe)), tur = 'fringe';
    else, tur = 'diger';
    end

    if any(strcmp(op, {'TEDX','TEDY','TETX','TETY','TETZ','TUDX','TUDY','TUTX','TUTY','TUTZ'}))
        sinif = 'Hizalama';
    elseif strcmp(op, 'TTHI')
        % Hava boslugu mekanik/hizalama, cam kalinligi uretim
        yy = Y([Y.No] == prm(1));
        if isempty(yy) || isempty(yy.Malzeme) || strcmpi(yy.Malzeme, 'AIR')
            sinif = 'Bosluk';
        else
            sinif = 'Uretim';
        end
    elseif any(strcmp(op, {'COMP','CPAR','CMCO'}))
        sinif = 'Kompanzator';
    elseif any(strcmp(op, {'STAT','TWAV','SAVE','SEED','TOFF','CEXP','TMCO','TMCC'}))
        sinif = 'Kontrol';
    elseif any(strcmp(op, {'TRAD','TCUR','TFRN','TCON','TIND','TABB','TSDX','TSDY','TSTX','TSTY', ...
                           'TIRX','TIRY','TIRR','TEZI','TEXI','TPAR','TEXA','TNPS','TCMU','TCIO','TCEO'}))
        sinif = 'Uretim';
    else
        sinif = 'Diger';
    end
end

function [birim, metin] = okunurBirim(tur, mn, mx, birimAd, umCarpan)
    buyuk = max(abs([mn mx]));
    simetrik = abs(mn + mx) <= 1e-12 * max(1, buyuk);
    switch tur
        case 'uzunluk'
            birim = birimAd;
            if simetrik, metin = sprintf('+/-%.4g um', buyuk * umCarpan);
            else, metin = sprintf('%.4g ... %.4g um', mn * umCarpan, mx * umCarpan); end
        case 'aci'
            birim = 'derece';
            if simetrik, metin = sprintf('+/-%.4g arcmin (%.4g arcsec)', buyuk * 60, buyuk * 3600);
            else, metin = sprintf('%.4g ... %.4g arcmin', mn * 60, mx * 60); end
        case 'fringe'
            birim = 'fringe';
            metin = sprintf('%.4g ... %.4g fringe', mn, mx);
        otherwise
            birim = '-';
            if isnan(buyuk), metin = '-'; else, metin = sprintf('%.4g ... %.4g', mn, mx); end
    end
end

function ad = aitOldugu(op, tur, prm, Y)
    if any(strcmp(op, {'STAT','TWAV','SAVE','SEED','TOFF','CEXP','TMCO','TMCC'}))
        ad = '-'; return;
    end
    aralikli = {'TEDX','TEDY','TETX','TETY','TETZ','TUDX','TUDY','TUTX','TUTY','TUTZ'};
    if any(strcmp(op, aralikli))
        ad = aralikAdi(prm(1:2), Y);
    else
        yy = Y([Y.No] == prm(1));
        if isempty(yy), ad = '-';
        else, ad = sprintf('S%d %s', prm(1), bosIse(yy.Yorum, ''));
        end
    end
    if strcmp(tur, 'diger') && isnan(prm(1)), ad = '-'; end
end

function ad = aralikAdi(ar, Y)
    % Ad: araligin ilk yuzeyinin yorumu; bossa araliktaki ilk dolu yorum; o da yoksa S#-S#
    ad = '';
    yy = Y([Y.No] == ar(1));
    if ~isempty(yy) && ~isempty(yy.Yorum), ad = yy.Yorum; end
    if isempty(ad)
        for k = ar(1):ar(2)
            yy = Y([Y.No] == k);
            if ~isempty(yy) && ~isempty(yy.Yorum), ad = yy.Yorum; break; end
        end
    end
    if isempty(ad), ad = sprintf('S%d-S%d', ar(1), ar(2)); end
end

function t = ozetle(T, satirlar, ops)
    p = {};
    for k = satirlar'
        if any(strcmp(T.Operand{k}, ops))
            p{end+1} = sprintf('%s %s', T.Operand{k}, T.Okunur{k}); %#ok<AGROW>
        end
    end
    if isempty(p), t = '-'; else, t = strjoin(p, '; '); end
end
