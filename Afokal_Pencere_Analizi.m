%% AFOKAL PENCERE ANALIZI (ZOS-API, MATLAB)
% Afokal on grup ile objektif arasindaki kolime bosluga konmus bir Paraxial
% yuzeyle yapay "afokallik hatasi" (vergens, phi) uretir. Her phi icin:
%   1) Kompanzator bosluguna (T_on) nominalden baslayarak local optimizasyon yapar
%   2) Kompanzator hareketini (dz) ve Nyquist'teki FFT MTF'yi (T ve S, her alan) okur
% Sonunda kompanzatorun tolere edebildigi phi penceresini bulur, grafik ve CSV uretir.
%
% DOSYA HAZIRLIGI (ozet): Paraxial yuzey eklenmis (f = 1e10 mm), kolime
% boslukta (afokal grup ile kompanzator arasinda). Kompanzator dummy yontemi:
% lensten onceki dummy kalinligi degisken (nominal 0), lensten sonraki dummy
% pickup -1. Baska degisken yok, tek konfigurasyon, default merit function hazir.
% Birimler mm olmali.
%
% Not: Turkce karakter kullanilmadi (eski MATLAB surumlerinde kodlama sorunu olmasin diye).

clear; clc;

%% ===================== KULLANICI AYARLARI =====================
ayar.mod          = 'standalone'; % 'standalone': dosyayi arka planda kendisi acar (diger scriptlerinle ayni yontem)
                                  % 'extension' : acik OpticStudio'ya baglan (once Programming > Interactive Extension)
ayar.dosya        = 'C:\Tasarimlar\MWIR_B_paraxial.zmx';  % sadece 'standalone' icin

ayar.paraxSurf    = 12;   % Paraxial yuzeyin numarasi
ayar.kompOnSurf   = 12;   % Kompanzatoru hareket ettiren DEGISKEN kalinlik
                          % (dummy yontemi: lensten onceki dummy yuzey, nominal 0)
ayar.kompArkaSurf = 14;   % Pickup -1 olan kalinlik (dummy yontemi: lensten sonraki dummy yuzey)

ayar.ters         = true;  % true: tasarim TERS (dedektor = obje, cikis sonsuzda, Afocal Image Space)
                           %       MTF cycles/mrad olarak, acisal Nyquist'te okunur
                           % false: normal yon, MTF lp/mm olarak goruntu duzleminde
ayar.efl          = 100;   % [mm] objektifin EFL'si (sadece ters=true icin; acisal Nyquist hesabi)

ayar.phiD         = linspace(-0.05, 0.05, 11);  % taranacak vergens [diyoptri = 1/m]; 0 otomatik eklenir
                                                % 0.05 D = f 20 m, 0.01 D = f 100 m
ayar.pikselUm     = 15;    % dedektor piksel araligi [um] -> Nyquist = 1000/(2*piksel) lp/mm
ayar.strok        = 0.30;  % kompanzator icin kullanilabilir strok, +/- [mm] (termal pay dusulmus)
ayar.mtfDusus     = 0.05;  % kabul edilen goreli MTF dususu (0.05 = nominalin %5'i)
ayar.cekirdek     = 8;     % optimizasyonda kullanilacak cekirdek sayisi
ayar.cikti        = 'afokal_pencere_sonuc.csv';

ayar.dosya = strtrim(strrep(ayar.dosya, '"', ''));   % "Yol olarak kopyala" tirnaklarini temizle

%% ===================== BAGLANTI =====================
TheApplication = baglan(ayar);
TheSystem = TheApplication.PrimarySystem;
if strcmp(ayar.mod, 'standalone')
    TheSystem.LoadFile(ayar.dosya, false);
end
LDE = TheSystem.LDE;

%% ===================== KONTROLLER =====================
sP    = LDE.GetSurfaceAt(ayar.paraxSurf);
sOn   = LDE.GetSurfaceAt(ayar.kompOnSurf);
sArka = LDE.GetSurfaceAt(ayar.kompArkaSurf);

fprintf('Acilan dosya: %s\n', char(TheSystem.SystemFile));
if ~paraxialMi(sP)
    % Dosyadaki tum Paraxial yuzeyleri bul ve listele
    bulunan = [];
    for i = 1:LDE.NumberOfSurfaces - 1
        if paraxialMi(LDE.GetSurfaceAt(i)), bulunan(end+1) = i; end %#ok<AGROW>
    end
    if isempty(bulunan)
        error(['Acilan dosyada hic Paraxial yuzey yok. ayar.dosya dogru dosyayi mi gosteriyor, ' ...
               've paraxial ekledikten sonra dosyayi KAYDETTIN mi?']);
    else
        error('Yuzey %d Paraxial degil (tip: %s). Paraxial yuzey(ler): %s -> ayar.paraxSurf''i buna gore duzelt.', ...
              ayar.paraxSurf, char(sP.TypeName), mat2str(bulunan));
    end
end
if ~contains(lower(char(TheSystem.SystemData.Units.LensUnits.ToString())), 'millimeter')
    warning('Lens birimi mm degil. Kod mm varsayiyor.');
end
if ~contains(char(sArka.ThicknessCell.GetSolveData().Type.ToString()), 'Pickup')
    warning('Yuzey %d kalinliginda pickup solve yok. Kompanzator hareket ederken toplam uzunluk korunmayacak.', ...
            ayar.kompArkaSurf);
end
if ~ayar.ters
    nSurf = LDE.NumberOfSurfaces;               % object dahil
    sGoruntuOncesi = LDE.GetSurfaceAt(nSurf - 2);
    solveTipi = char(sGoruntuOncesi.ThicknessCell.GetSolveData().Type.ToString());
    if ~any(strcmp(solveTipi, {'Fixed', 'Variable'}))
        warning('Goruntu oncesi kalinlikta "%s" solve var. Odak otomatik duzeltilir, test anlamsizlasir. Kaldir!', solveTipi);
    end
else
    try
        if ~TheSystem.SystemData.Aperture.AFocalImageSpace
            warning('ters=true ama Afocal Image Space kapali. System Explorer > Aperture''dan ac.');
        end
    catch
        % ozellik adi surume gore farkli olabilir; elle kontrol et
    end
end
if TheSystem.MCE.NumberOfConfigurations > 1
    warning('Birden fazla konfigurasyon var. Merit function hepsini kapsiyorsa sonuclar karisir. Tek konfigurasyon onerilir.');
end
if TheSystem.MFE.NumberOfOperands < 1
    error('Merit function bos. Optimization Wizard ile default merit function olustur.');
end

%% ===================== HAZIRLIK =====================
parCell = sP.GetSurfaceCell(ZOSAPI.Editors.LDE.SurfaceColumn.Par1);  % Paraxial: Par1 = Focal Length
f_orj   = parCell.DoubleValue;
T_on0   = sOn.Thickness;
T_arka0 = sArka.Thickness;

if ayar.ters
    fNyq  = ayar.efl / (2 * ayar.pikselUm);         % cycles/mrad  (= EFL[mm] / (2 * p[um]))
    birim = 'cycles/mrad';
else
    fNyq  = 1000 / (2 * ayar.pikselUm);             % lp/mm
    birim = 'lp/mm';
end
nAlan = TheSystem.SystemData.Fields.NumberOfFields;
phiD  = unique([ayar.phiD(:); 0]);                  % sirali, 0 dahil
n     = numel(phiD);
i0    = find(phiD == 0);

fprintf('Nyquist: %.3f %s | Alan sayisi: %d | Nokta sayisi: %d\n', fNyq, birim, nAlan, n);
fprintf('Nominal T_on = %.4f mm, T_arka = %.4f mm\n\n', T_on0, T_arka0);

TheSystem.Tools.RemoveAllVariables();
sOn.ThicknessCell.MakeSolveVariable();

dz   = nan(n, 1);  mf = nan(n, 1);
mtfT = nan(n, nAlan);  mtfS = nan(n, nAlan);

%% ===================== TARAMA =====================
try
    for k = 1:n
        if phiD(k) == 0
            f = 1e10;                     % pratikte sifir guc
        else
            f = 1000 / phiD(k);           % D -> mm
        end
        parCell.DoubleValue = f;
        sOn.Thickness = T_on0;            % her adimda nominalden basla

        optimizeEt(TheSystem, ayar.cekirdek);

        dz(k) = sOn.Thickness - T_on0;
        mf(k) = TheSystem.MFE.CalculateMeritFunction();
        [mtfT(k,:), mtfS(k,:)] = mtfOku(TheSystem, fNyq, nAlan);

        fprintf('phi = %+7.4f D (f = %10.3g mm) | dz = %+8.4f mm | MF = %.5g | MTF_T(eksen) = %.3f\n', ...
                phiD(k), f, dz(k), mf(k), mtfT(k,1));
    end
catch hata
    geriYukle(sP, parCell, f_orj, sOn, T_on0);
    rethrow(hata);
end
geriYukle(sP, parCell, f_orj, sOn, T_on0);

%% ===================== DEGERLENDIRME =====================
dzRel = dz - dz(i0);                      % nominal tam optimum degilse bile dogru fark
nT = mtfT ./ mtfT(i0,:);                  % nominale gore normalize MTF
nS = mtfS ./ mtfS(i0,:);

okMTF   = all(nT >= (1 - ayar.mtfDusus), 2) & all(nS >= (1 - ayar.mtfDusus), 2);
okStrok = abs(dzRel) <= ayar.strok;
ok      = okMTF & okStrok;

lo = i0; while lo > 1 && ok(lo-1), lo = lo - 1; end
hi = i0; while hi < n && ok(hi+1), hi = hi + 1; end

egim = polyfit(phiD, dzRel, 1);           % mm / D

fprintf('\n===== SONUC =====\n');
fprintf('Kompanzator duyarliligi: %.4f mm/D (dogrusal fit)\n', egim(1));
fprintf('Kabul penceresi (taranan noktalarda): %+.4f D ... %+.4f D\n', phiD(lo), phiD(hi));
if phiD(lo) ~= 0, fprintf('  -> negatif tarafta |f| >= %.1f m\n', 1/abs(phiD(lo))); end
if phiD(hi) ~= 0, fprintf('  -> pozitif tarafta |f| >= %.1f m\n', 1/abs(phiD(hi))); end
if lo == 1 || hi == n
    fprintf('UYARI: Pencere tarama sinirina dayandi. ayar.phiD araligini genislet.\n');
end
sinirlayan = '';
if any(~okStrok), sinirlayan = [sinirlayan 'strok ']; end
if any(~okMTF),   sinirlayan = [sinirlayan 'MTF'];    end
if ~isempty(sinirlayan), fprintf('Pencereyi sinirlayan: %s\n', sinirlayan); end

%% ===================== CSV =====================
tablo = table(phiD, 1000./phiD, dzRel, mf, okStrok, okMTF, ...
    'VariableNames', {'phi_D', 'f_mm', 'dz_komp_mm', 'MeritFunction', 'Strok_OK', 'MTF_OK'});
for a = 1:nAlan
    tablo.(sprintf('MTF_T_alan%d', a)) = mtfT(:,a);
    tablo.(sprintf('MTF_S_alan%d', a)) = mtfS(:,a);
end
writetable(tablo, ayar.cikti);
fprintf('Sonuclar yazildi: %s\n', ayar.cikti);

%% ===================== GRAFIKLER =====================
figure('Name', 'Kompanzator hareketi');
plot(phiD, dzRel, 'o-', 'LineWidth', 1.5); hold on; grid on;
yline( ayar.strok, 'r--', 'strok siniri');
yline(-ayar.strok, 'r--');
xline(phiD(lo), 'g:'); xline(phiD(hi), 'g:');
xlabel('\phi [D]'); ylabel('\Deltaz_{komp} [mm]');
title(sprintf('Kompanzator hareketi (egim %.3f mm/D)', egim(1)));

figure('Name', 'Normalize MTF');
renk = lines(nAlan); hold on; grid on;
for a = 1:nAlan
    plot(phiD, nT(:,a), '-o', 'Color', renk(a,:), 'DisplayName', sprintf('Alan %d T', a));
    plot(phiD, nS(:,a), '--s', 'Color', renk(a,:), 'DisplayName', sprintf('Alan %d S', a));
end
yline(1 - ayar.mtfDusus, 'r--', 'esik', 'HandleVisibility', 'off');
xline(phiD(lo), 'g:', 'HandleVisibility', 'off'); xline(phiD(hi), 'g:', 'HandleVisibility', 'off');
xlabel('\phi [D]'); ylabel('MTF / MTF_{nominal}');
title(sprintf('Kompanzasyon sonrasi MTF @ %.3f %s', fNyq, birim));
legend('Location', 'best');

if strcmp(ayar.mod, 'standalone')
    TheApplication.CloseApplication();
end

%% ===================== YARDIMCI FONKSIYONLAR =====================
function app = baglan(ayar)
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
    if strcmp(ayar.mod, 'extension')
        app = conn.ConnectAsExtension(0);
    else
        app = conn.CreateNewApplication();
    end
    if isempty(app)
        error('OpticStudio''ya baglanilamadi.');
    end
    if ~app.IsValidLicenseForAPI
        if strcmp(ayar.mod, 'extension')
            error(['Extension baglantisi kurulamadi. OpticStudio''da Programming > Interactive Extension''a ' ...
                   'basip "Waiting for connection" yazisini gordukten sonra calistir; ya da ayar.mod = ''standalone'' kullan.']);
        else
            error('Lisans kontrolu basarisiz (standalone).');
        end
    end
end

function e = paraxialMi(surf)
    e = contains(lower(char(surf.TypeName)), 'paraxial') || ...
        contains(lower(char(surf.Type.ToString())), 'paraxial');
end

function optimizeEt(TheSystem, cekirdek)
    opt = TheSystem.Tools.OpenLocalOptimization();
    opt.Algorithm = ZOSAPI.Tools.Optimization.OptimizationAlgorithm.DampedLeastSquares;
    opt.Cycles = ZOSAPI.Tools.Optimization.OptimizationCycles.Automatic;
    opt.NumberOfCores = cekirdek;
    opt.RunAndWaitForCompletion();
    opt.Close();
end

function [T, S] = mtfOku(TheSystem, fNyq, nAlan)
    a = TheSystem.Analyses.New_FftMtf();
    s = a.GetSettings();
    s.MaximumFrequency = 1.5 * fNyq;
    s.SampleSize = ZOSAPI.Analysis.SampleSizes.S_64x64;
    try, s.ShowDiffractionLimit = false; catch, end
    a.ApplyAndWaitForCompletion();
    r = a.GetResults();
    if r.NumberOfDataSeries ~= nAlan
        warning('MTF veri serisi sayisi (%d) alan sayisina (%d) esit degil. Sira kontrol edilmeli.', ...
                r.NumberOfDataSeries, nAlan);
    end
    T = nan(1, nAlan); S = nan(1, nAlan);
    for i = 0:min(r.NumberOfDataSeries, nAlan) - 1
        ds = r.GetDataSeries(i);
        x = ds.XData.Data.double;
        y = ds.YData.Data.double;     % sutun 1: tanjansiyel, sutun 2: sagital
        T(i+1) = interp1(x, y(:,1), fNyq);
        S(i+1) = interp1(x, y(:,2), fNyq);
    end
    a.Close();
end

function geriYukle(~, parCell, f_orj, sOn, T_on0)
    parCell.DoubleValue = f_orj;
    sOn.Thickness = T_on0;
    sOn.ThicknessCell.MakeSolveFixed();
end
