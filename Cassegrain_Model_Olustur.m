%% CASSEGRAIN MODEL OLUSTURUCU (ZOS-API, MATLAB)
% NE YAPAR
%   Kullanicidan isterleri (EFL, F#, obscuration, back focus, dalga boyu,
%   alan, boyut sinirlari) ve Cassegrain turunu/turlerini alir. Birinci
%   derece (paraksiyel) geometriyi ve aynalarin conic sabitlerini formullerle
%   hesaplar, isterleri kontrol eder, ardindan OpticStudio'da iki aynali
%   modeli kurar, istenirse optimize eder ve .zmx olarak kaydeder.
%   Birden fazla tur secilirse her biri ayri dosyaya kaydedilir ve sonuclar
%   bir tabloda karsilastirilir.
%
% NASIL CALISIR
%   1. Isterler bir pencereden girilir (varsayilanlar ortak gereksinim seti).
%   2. Tur(ler) secilir: Klasik Cassegrain, Ritchey-Chretien, Dall-Kirkham,
%      Pressman-Carmichael.
%   3. Kayit klasoru secilir.
%   4. Her tur icin: hesap -> kontrol -> LDE kurulumu -> (optimizasyon) ->
%      sonuc okuma -> kayit.
%   Ciktilar: <klasor>\Cassegrain_<tur>.zmx, Cassegrain_Sonuclar.xlsx,
%   Command Window ozeti.
%
% HESAPLAR (ince ayna, tum mesafeler pozitif buyukluk olarak)
%   D   = f / F#                     giris acikligi (M1 capi)
%   b   = k f                        M2 -> odak mesafesi
%   d   = b - e                      aynalar arasi mesafe
%   f1  = (k f - e) / (1 - k)        M1 odak uzakligi
%   m   = f / f1                     ikincil buyutme
%   f2  = (f1 - d) b / (b - (f1 - d))  M2 odak uzakligi (dis bukey)
%   R1  = -2 f1 ,  R2 = -2 f2        Zemax isaretiyle yaricaplar
%   DM2 = k D + 2 d tan(theta)       alan dahil M2 capi
%   beta = e / f1
%   Klasik : K1 = -1, K2 = -((m+1)/(m-1))^2
%   RC     : K1 = -1 - 2(1+beta) / (m^2 (m-beta))
%            K2 = -((m+1)/(m-1))^2 - 2m(m+1) / ((m-beta)(m-1)^3)
%   DK, PC : kuresel = 0 kosulu
%            K1 + 1 = k (m-1)^3 / m^3 * [K2 + ((m+1)/(m-1))^2]
%            DK: K2 = 0 -> K1 ;  PC: K1 = 0 -> K2
%   Formuller 3. derece (Seidel) sonuclaridir. Hizli M1'de yuksek dereceden
%   kuresel aberasyon kalir; bu yuzden conic'ler optimizasyonda serbest
%   birakilir (Klasik haric: parabol + hiperbol eksende tam kusursuzdur).
%
% ZEMAX MODELI (sequential)
%   0 OBJ   sonsuz
%   1 Dummy  t = d,  Circular Obscuration (M2 + baffle golgesi)
%   2 M1     STOP, R1, t = -d, MIRROR, K1, Circular Aperture (min = M1 deligi)
%   3 M2     R2, t = b, MIRROR, K2, Circular Aperture
%   4 IMA
%   Aperture: Entrance Pupil Diameter = D. Alanlar: 0, 0.7, 1.0 x yari alan.
%   Ray aiming: Real.
%   Optimizasyon degiskenleri: goruntu mesafesi + turun serbest conic'leri.
%   R1, R2, d sabit kalir, boylece EFL ve obscuration birinci derecede korunur.
%   Merit: RMS wavefront, centroid, Gaussian quadrature, obscuration orani ile.
%
% NOTLAR
%   - Raporlanan RMS WFE / spot degerleri RWCE / RSCE operandlarindandir
%     (Gaussian quadrature). Obscuration'i tam hesaba katmazlar; kesin degerler
%     icin Zemax'ta Wavefront Map ve Spot Diagram'a bakilmali.
%   - Model bir baslangic noktasidir; tolerans ve katadioptrik asama sonra.
%   - Turkce karakter kullanilmadi (eski MATLAB surumlerinde kodlama sorunu olmasin diye).

clc; clear; close all;

%% ===================== ISTERLER =====================
prompt = { ...
    'EFL (mm):', ...
    'F# :', ...
    'Maks. obscuration orani (cap, baffle dahil, 0-1):', ...
    'Tasarim eksenel obscuration k (M2 eksen isin capi / D, 0-1):', ...
    'Back focus e: odagin M1 tepesinin arkasindaki mesafesi (mm):', ...
    'Calisma dalga boyu (um):', ...
    'Yari alan acisi (derece):', ...
    'Maks. cap (mm):', ...
    'Maks. uzunluk (mm):', ...
    'M2 baffle / tutucu radyal payi (mm):', ...
    'Ek uzunluk payi: M2 arkasi + dedektor + M1 kalinligi (mm):', ...
    'Optimizasyon yapilsin mi? (1 = evet, 0 = hayir):', ...
    'Optimizasyon cekirdek sayisi:'};
varsayilan = {'400', '2.3', '0.30', '0.27', '20', '0.632', '0.5', ...
              '180', '160', '1', '15', '1', '8'};
cevap = inputdlg(prompt, 'Cassegrain isterleri', [1 70], varsayilan);
if isempty(cevap), return; end
v = str2double(cevap);
if any(~isfinite(v)), error('Tum alanlara sayi girilmeli.'); end

ist.f        = v(1);
ist.Fno      = v(2);
ist.kMax     = v(3);
ist.k        = v(4);
ist.e        = v(5);
ist.lambda   = v(6);
ist.theta    = v(7);
ist.Dmax     = v(8);
ist.Lmax     = v(9);
ist.bafflePay = v(10);
ist.ekUzunluk = v(11);
ist.optimize = v(12) ~= 0;
ist.cekirdek = max(1, round(v(13)));

if ist.f <= 0 || ist.Fno <= 0, error('EFL ve F# pozitif olmali.'); end
if ist.k <= 0 || ist.k >= 1,   error('k 0 ile 1 arasinda olmali.'); end
if ist.e < 0,                  error('Back focus e negatif olamaz.'); end
if ist.theta < 0,              error('Yari alan acisi negatif olamaz.'); end

turAdlari = {'Klasik Cassegrain', 'Ritchey-Chretien', 'Dall-Kirkham', 'Pressman-Carmichael'};
turKod    = {'Klasik', 'RC', 'DK', 'PC'};
[secim, ok] = listdlg('ListString', turAdlari, 'SelectionMode', 'multiple', ...
    'InitialValue', 2, 'Name', 'Cassegrain turu', ...
    'PromptString', 'Tur(ler) secin (Ctrl ile coklu):', 'ListSize', [260 90]);
if ~ok || isempty(secim), return; end

klasor = uigetdir(pwd, 'Zemax dosyalarinin kaydedilecegi klasoru secin');
if isequal(klasor, 0), return; end

%% ===================== BIRINCI DERECE HESAP =====================
g = ilkDerece(ist);
kontrolEt(ist, g);
ozetYaz(ist, g);

%% ===================== ZEMAX =====================
TheApplication = baglan();
TheSystem = TheApplication.PrimarySystem;

sonuc = struct([]);
for s = secim
    tur = turKod{s};
    [K1, K2] = conicHesap(tur, g);
    fprintf('\n================ %s ================\n', turAdlari{s});
    fprintf('Formul conic: K1 = %.4f, K2 = %.4f\n', K1, K2);

    modelKur(TheSystem, ist, g, K1, K2, tur, turAdlari{s});

    efl0 = operand(TheSystem, 'EFFL', 0, 0, 0, 0);
    fprintf('Paraksiyel EFL (Zemax): %.4f mm (hedef %.4f)\n', efl0, ist.f);
    if abs(abs(efl0) - ist.f) > 1e-3 * ist.f
        warning('Zemax EFL hedeften %%0.1''den fazla farkli. Model kontrol edilmeli.');
    end

    if ist.optimize
        meritKur(TheSystem, ist, g);
        mf0 = TheSystem.MFE.CalculateMeritFunction();
        optimizeEt(TheSystem, ist.cekirdek);
        mf1 = TheSystem.MFE.CalculateMeritFunction();
        fprintf('Merit: %.6f -> %.6f\n', mf0, mf1);
    end

    r = sonucOku(TheSystem, ist, g);
    r.Tur = string(turAdlari{s});
    r.K1_formul = K1;  r.K2_formul = K2;

    dosya = fullfile(klasor, sprintf('Cassegrain_%s.zmx', tur));
    TheSystem.SaveAs(dosya);
    r.Dosya = string(dosya);
    fprintf('Kaydedildi: %s\n', dosya);

    sonuc = [sonuc; r]; %#ok<AGROW>
end

TheApplication.CloseApplication();

%% ===================== KARSILASTIRMA =====================
T = struct2table(sonuc);
T = T(:, {'Tur', 'K1_formul', 'K2_formul', 'K1', 'K2', 'EFL_mm', 'WFNO', ...
          'GoruntuMesafesi_mm', 'e_mm', ...
          'WFE0_dalga', 'WFE07_dalga', 'WFE1_dalga', ...
          'Spot0_um', 'Spot07_um', 'Spot1_um', 'Dosya'});

fprintf('\n================ KARSILASTIRMA ================\n');
fprintf('Alanlar: 0, %.3f, %.3f derece. WFE: RMS dalga (centroid), Spot: RMS yaricap (um)\n', ...
    0.7 * ist.theta, ist.theta);
fprintf('%-20s %9s %9s %9s %9s %9s %9s %9s %9s\n', 'Tur', 'K1', 'K2', ...
    'WFE 0', 'WFE 0.7', 'WFE 1', 'Spot 0', 'Spot 0.7', 'Spot 1');
for i = 1:height(T)
    fprintf('%-20s %9.4f %9.4f %9.4f %9.4f %9.4f %9.2f %9.2f %9.2f\n', T.Tur(i), ...
        T.K1(i), T.K2(i), T.WFE0_dalga(i), T.WFE07_dalga(i), T.WFE1_dalga(i), ...
        T.Spot0_um(i), T.Spot07_um(i), T.Spot1_um(i));
end
fprintf('Kirinim siniri: RMS WFE <= 0.07 dalga (Marechal), Airy yaricapi %.2f um\n', g.airy_um);

xlsx = fullfile(klasor, 'Cassegrain_Sonuclar.xlsx');
writetable(T, xlsx, 'Sheet', 'Turler');
writetable(ilkDereceTablo(ist, g), xlsx, 'Sheet', 'Birinci derece');
fprintf('\nSonuclar: %s\n', xlsx);

%% ===================== YARDIMCI FONKSIYONLAR =====================
function g = ilkDerece(ist)
    f = ist.f;  k = ist.k;  e = ist.e;
    g.D    = f / ist.Fno;
    g.b    = k * f;
    g.d    = g.b - e;
    g.f1   = (k * f - e) / (1 - k);
    g.m    = f / g.f1;
    g.f2   = (g.f1 - g.d) * g.b / (g.b - (g.f1 - g.d));
    g.R1   = -2 * g.f1;
    g.R2   = -2 * g.f2;
    g.beta = e / g.f1;
    g.M1F  = g.f1 / g.D;
    th = tand(ist.theta);
    g.yImg    = f * th;                         % goruntu yuksekligi (yari)
    g.DM2     = k * g.D + 2 * g.d * th;         % alan dahil M2 capi
    g.rObs    = g.DM2 / 2 + ist.bafflePay;      % golge yaricapi (M2 + baffle)
    g.kAlan   = 2 * g.rObs / g.D;               % baffle dahil obscuration orani
    % M1 deligi: M2'den odaga giden isin demetinin M1 duzlemindeki yaricapi + 1 mm
    g.rDelik  = min(k * g.D / 2 * e / g.b + max(g.d, f) * th + 1, g.rObs);
    g.L       = g.d + e + ist.ekUzunluk;        % tahmini toplam uzunluk
    g.airy_um = 1.22 * ist.lambda * ist.Fno;
    g.odakKaymasi = g.m^2 + 1;                  % odak kaymasi / despace
    g.odakDerinligi_um = 2 * ist.lambda * ist.Fno^2;
end

function kontrolEt(ist, g)
    if g.d <= 0
        error('d = k f - e <= 0. Back focus e cok buyuk (e < k f = %.2f mm olmali).', ist.k * ist.f);
    end
    if g.f1 <= g.d
        error('f1 <= d: M2 odaktan sonra kaliyor (Gregory geometrisi). Girdiler kontrol edilmeli.');
    end
    if g.m <= 1
        error('Ikincil buyutme m <= 1. Cassegrain icin m > 1 olmali.');
    end
    if g.kAlan > ist.kMax
        warning('Baffle dahil obscuration %.3f > maks %.3f. k veya baffle payi kucultulmeli.', ...
            g.kAlan, ist.kMax);
    end
    if g.D > ist.Dmax
        warning('Giris acikligi D = %.2f mm > maks cap %.2f mm.', g.D, ist.Dmax);
    elseif ist.Dmax - g.D < 10
        fprintf('Not: capta kenar basina yalnizca %.1f mm mekanik pay kaliyor.\n', (ist.Dmax - g.D) / 2);
    end
    if g.L > ist.Lmax
        warning('Tahmini uzunluk %.1f mm > maks %.1f mm.', g.L, ist.Lmax);
    end
    if g.M1F < 0.7
        warning('M1 F/%.2f cok hizli: uretim ve hizalama hassasiyeti cok yuksek olur.', g.M1F);
    end
end

function ozetYaz(ist, g)
    fprintf('================ BIRINCI DERECE ================\n');
    fprintf('Giris acikligi D        : %.3f mm\n', g.D);
    fprintf('M1 odak uzakligi f1     : %.3f mm (M1 F/%.3f)\n', g.f1, g.M1F);
    fprintf('Ikincil buyutme m       : %.4f\n', g.m);
    fprintf('M2 odak uzakligi f2     : %.3f mm\n', g.f2);
    fprintf('R1 / R2                 : %.3f / %.3f mm\n', g.R1, g.R2);
    fprintf('Aynalar arasi d         : %.3f mm\n', g.d);
    fprintf('M2 -> odak b            : %.3f mm\n', g.b);
    fprintf('Back focus e            : %.3f mm\n', ist.e);
    fprintf('Goruntu yari yuksekligi : %.4f mm\n', g.yImg);
    fprintf('M2 capi (alan dahil)    : %.3f mm\n', g.DM2);
    fprintf('Obscuration (baffle ile): %.3f (maks %.3f)\n', g.kAlan, ist.kMax);
    fprintf('M1 delik yaricapi       : %.3f mm\n', g.rDelik);
    fprintf('Tahmini uzunluk         : %.1f mm (maks %.1f)\n', g.L, ist.Lmax);
    fprintf('Airy yaricapi           : %.3f um\n', g.airy_um);
    fprintf('Odak derinligi          : +/- %.2f um\n', g.odakDerinligi_um);
    fprintf('Odak kaymasi / despace  : %.2f um/um (m^2+1)\n', g.odakKaymasi);
end

function T = ilkDereceTablo(ist, g)
    ad = {'EFL_mm'; 'F#'; 'D_mm'; 'k'; 'e_mm'; 'd_mm'; 'b_mm'; 'f1_mm'; 'M1_F#'; ...
          'm'; 'f2_mm'; 'R1_mm'; 'R2_mm'; 'DM2_mm'; 'Obscuration_baffle_ile'; ...
          'M1_delik_yaricap_mm'; 'Uzunluk_tahmini_mm'; 'Airy_yaricap_um'; ...
          'Odak_derinligi_pm_um'; 'Odak_kaymasi_um_per_um'};
    deger = [ist.f; ist.Fno; g.D; ist.k; ist.e; g.d; g.b; g.f1; g.M1F; g.m; g.f2; ...
             g.R1; g.R2; g.DM2; g.kAlan; g.rDelik; g.L; g.airy_um; ...
             g.odakDerinligi_um; g.odakKaymasi];
    T = table(ad, deger, 'VariableNames', {'Buyukluk', 'Deger'});
end

function [K1, K2] = conicHesap(tur, g)
    m = g.m;  k = 1 - g.d / g.f1;  be = g.beta;
    A = ((m + 1) / (m - 1))^2;
    switch tur
        case 'Klasik'
            K1 = -1;
            K2 = -A;
        case 'RC'
            K1 = -1 - 2 * (1 + be) / (m^2 * (m - be));
            K2 = -A - 2 * m * (m + 1) / ((m - be) * (m - 1)^3);
        case 'DK'
            K2 = 0;
            K1 = k * (m - 1)^3 / m^3 * (K2 + A) - 1;
        case 'PC'
            K1 = 0;
            K2 = (K1 + 1) * m^3 / (k * (m - 1)^3) - A;
    end
end

function modelKur(TheSystem, ist, g, K1, K2, tur, turAdi)
    TheSystem.New(false);
    SD = TheSystem.SystemData;
    try, SD.TitleNotes.Title = sprintf('%s f=%g F/%g', turAdi, ist.f, ist.Fno); catch, end

    SD.Aperture.ApertureType = ZOSAPI.SystemData.ZemaxApertureType.EntrancePupilDiameter;
    SD.Aperture.ApertureValue = g.D;
    try, SD.RayAiming.RayAiming = ZOSAPI.SystemData.RayAimingMethod.Real; catch, end

    SD.Wavelengths.GetWavelength(1).Wavelength = ist.lambda;
    SD.Wavelengths.GetWavelength(1).Weight = 1;

    F = SD.Fields;
    F.GetField(1).Y = 0;
    if ist.theta > 0
        F.AddField(0, 0.7 * ist.theta, 1);
        F.AddField(0, ist.theta, 1);
    end

    LDE = TheSystem.LDE;               % yeni sistem: 0 OBJ, 1 STO, 2 IMA
    LDE.InsertNewSurfaceAt(2);
    LDE.InsertNewSurfaceAt(2);         % -> 0 OBJ, 1, 2, 3, 4 IMA

    s1 = LDE.GetSurfaceAt(1);
    s1.Comment = 'M2 golgesi (dummy)';
    s1.Thickness = g.d;
    aciklikAyarla(s1, 'obscuration', 0, g.rObs);

    s2 = LDE.GetSurfaceAt(2);
    s2.Comment = 'M1';
    s2.Radius = g.R1;
    s2.Thickness = -g.d;
    s2.Material = 'MIRROR';
    s2.Conic = K1;
    durdurmaYap(LDE, s2, 2);
    aciklikAyarla(s2, 'acik', g.rDelik, g.D / 2 + 0.5);

    s3 = LDE.GetSurfaceAt(3);
    s3.Comment = 'M2';
    s3.Radius = g.R2;
    s3.Thickness = g.b;
    s3.Material = 'MIRROR';
    s3.Conic = K2;
    aciklikAyarla(s3, 'acik', 0, g.DM2 / 2 + 0.5);

    % Optimizasyon degiskenleri
    s3.ThicknessCell.MakeSolveVariable();          % goruntu mesafesi (odak)
    switch tur
        case 'RC'
            s2.ConicCell.MakeSolveVariable();
            s3.ConicCell.MakeSolveVariable();
        case 'DK'
            s2.ConicCell.MakeSolveVariable();
        case 'PC'
            s3.ConicCell.MakeSolveVariable();
    end
end

function durdurmaYap(LDE, surf, idx)
    try
        surf.IsStop = true;
    catch
        try
            LDE.StopSurface = idx;
        catch
            warning('STOP yuzeyi M1''e tasinamadi. Zemax''ta elle yapilmali (yuzey %d).', idx);
        end
    end
end

function aciklikAyarla(surf, tip, rMin, rMax)
    AD = surf.ApertureData;
    if strcmp(tip, 'obscuration')
        tt = ZOSAPI.Editors.LDE.SurfaceApertureTypes.CircularObscuration;
        adlar = {'S_CircularObscuration_', 'x_S_CircularObscuration', 'S_CircularObscuration'};
    else
        tt = ZOSAPI.Editors.LDE.SurfaceApertureTypes.CircularAperture;
        adlar = {'S_CircularAperture_', 'x_S_CircularAperture', 'S_CircularAperture'};
    end
    ayar = AD.CreateApertureTypeSettings(tt);
    ic = [];
    for i = 1:numel(adlar)
        try, ic = ayar.(adlar{i}); break; catch, end
    end
    if isempty(ic), ic = ayar; end
    try
        ic.MinimumRadius = rMin;
        ic.MaximumRadius = rMax;
    catch
        warning('Yuzey "%s" icin aciklik yaricaplari yazilamadi; Zemax''ta elle girilmeli (min %.3f, max %.3f).', ...
            char(surf.Comment), rMin, rMax);
    end
    AD.ChangeApertureTypeSettings(ayar);
end

function meritKur(TheSystem, ist, g)
    MFE = TheSystem.MFE;
    kuruldu = false;
    try
        w = MFE.SEQOptimizationWizard;
        w.Type = 0;          % RMS
        w.Data = 0;          % Wavefront
        w.Reference = 0;     % Centroid
        w.Ring = 3;          % 4 halka
        w.Arm = 0;           % 6 kol
        try, w.Obscuration = min(g.kAlan, 0.95); catch, end
        try, w.IsGlassUsed = false; catch, end
        try, w.IsAirUsed = false; catch, end
        try
            w.Apply();
        catch
            w.CommonSettings.OK();
        end
        kuruldu = MFE.NumberOfOperands > 1;
    catch
    end
    if ~kuruldu
        % Yedek: her alan icin RMS wavefront (centroid) operandi
        warning('Optimization Wizard kullanilamadi, RWCE operandlari ile merit kuruluyor.');
        hy = [0 0.7 1];
        if ist.theta == 0, hy = 0; end
        for i = 1:numel(hy)
            op = MFE.AddOperand();
            op.ChangeType(ZOSAPI.Editors.MFE.MeritOperandType.RWCE);
            op.GetOperandCell(ZOSAPI.Editors.MFE.MeritColumn.Param1).IntegerValue = 6;
            op.GetOperandCell(ZOSAPI.Editors.MFE.MeritColumn.Param2).IntegerValue = 1;
            op.GetOperandCell(ZOSAPI.Editors.MFE.MeritColumn.Param4).DoubleValue = hy(i);
            op.Target = 0;
            op.Weight = 1;
        end
    end
    % EFL izleme satiri (agirlik 0; R1, R2, d sabit oldugundan EFL degismez)
    op = MFE.InsertNewOperandAt(1);
    op.ChangeType(ZOSAPI.Editors.MFE.MeritOperandType.EFFL);
    op.Target = ist.f;
    op.Weight = 0;
end

function optimizeEt(TheSystem, cekirdek)
    opt = TheSystem.Tools.OpenLocalOptimization();
    opt.Algorithm = ZOSAPI.Tools.Optimization.OptimizationAlgorithm.DampedLeastSquares;
    opt.Cycles = ZOSAPI.Tools.Optimization.OptimizationCycles.Automatic;
    opt.NumberOfCores = cekirdek;
    opt.RunAndWaitForCompletion();
    opt.Close();
end

function r = sonucOku(TheSystem, ist, g)
    LDE = TheSystem.LDE;
    r.K1 = LDE.GetSurfaceAt(2).Conic;
    r.K2 = LDE.GetSurfaceAt(3).Conic;
    r.GoruntuMesafesi_mm = LDE.GetSurfaceAt(3).Thickness;
    r.e_mm = r.GoruntuMesafesi_mm - g.d;
    r.EFL_mm = operand(TheSystem, 'EFFL', 0, 0, 0, 0);
    r.WFNO = operand(TheSystem, 'WFNO', 0, 0, 0, 0);
    hy = [0 0.7 1];
    wfe = nan(1, 3);  spot = nan(1, 3);
    for i = 1:3
        if ist.theta == 0 && hy(i) > 0, continue; end
        wfe(i)  = operand(TheSystem, 'RWCE', 6, 1, 0, hy(i));
        spot(i) = 1000 * operand(TheSystem, 'RSCE', 6, 1, 0, hy(i));
    end
    r.WFE0_dalga = wfe(1);  r.WFE07_dalga = wfe(2);  r.WFE1_dalga = wfe(3);
    r.Spot0_um = spot(1);   r.Spot07_um = spot(2);   r.Spot1_um = spot(3);
    fprintf('EFL %.4f mm, WFNO %.4f, K1 %.4f, K2 %.4f, e %.3f mm\n', ...
        r.EFL_mm, r.WFNO, r.K1, r.K2, r.e_mm);
    fprintf('RMS WFE (dalga) 0/0.7/1 alan: %.4f / %.4f / %.4f\n', wfe);
end

function val = operand(TheSystem, ad, i1, i2, hx, hy)
    tip = ZOSAPI.Editors.MFE.MeritOperandType.(ad);
    val = TheSystem.MFE.GetOperandValue(tip, i1, i2, hx, hy, 0, 0, 0, 0);
end

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
