%% AFOKAL TOLERANS ANALIZI (ZOS-API, MATLAB)
% "Afokal grubun bozulmasini kompanzatorle tolere edebilir miyim?" sorusunu
% iki parcada cevaplar:
%
%  PARCA 1 - KAPASITE (Dosya B, paraxial yuzeyli dosya)
%    Kolime bosluktaki Paraxial yuzeyle yapay vergens (phi) verilir, her adimda
%    kompanzator optimize edilir. Cikti: kompanzator duyarliligi [mm/D] ve
%    kabul penceresi [D] (strok + MTF kosulu).
%
%  PARCA 2 - IHTIYAC (Dosya A'nin Monte Carlo dosyalari)
%    Zemax'te sadece afokal grubun ic toleranslariyla kosulan Monte Carlo'nun
%    kaydettigi MC_T*.zmx dosyalari tek tek acilir, kompanzatorun ne kadar
%    hareket ettigi okunur ve Parca 1'deki egimle vergense cevrilir.
%
%  KARAR: Gercek vergens dagilimi kabul penceresinin icinde mi?
%
% Yuzey numaralari iki dosyada farkli olabilir: her dosyanin ayari ayri.
% 0 birakilan numaralar otomatik bulunur ve ekrana yazdirilir.
%
% Not: Turkce karakter kullanilmadi (eski MATLAB surumlerinde kodlama sorunu olmasin diye).

clc; clear; close all;

%% ===================== KULLANICI AYARLARI =====================
% ---------- Genel ----------
ayar.ters      = true;   % true: ters tasarim (Afocal Image Space, MTF cycles/mrad)
ayar.efl       = 239.5;  % [mm] tum sistemin EFL'si (sadece ters=true icin)
ayar.pikselUm  = 15;     % dedektor piksel araligi [um]
ayar.strok     = 0.50;   % [mm] kompanzatorun toplam mekanik stroku, +/-
ayar.termal    = 0.00;   % [mm] bu stroktan termal duzeltme icin ayrilan pay, +/-
ayar.mtfDusus  = 0.05;   % Parca 1'de kabul edilen goreli MTF dususu (0.05 = %5)
ayar.cekirdek  = 8;      % optimizasyon cekirdek sayisi
ayar.cikti     = 'afokal_tolerans_sonuc';   % cikti dosyalarinin on adi (.csv eklenir)

% ---------- PARCA 1: Dosya B (paraxial tarama) ----------
ayar.B.calistir   = true;   % false: tarama yapma, asagidaki egim/pencere degerlerini kullan
ayar.B.dosya      = 'C:\...\DOSYA_B_paraxial.zmx';
ayar.B.paraxSurf  = 0;      % 0 = otomatik (dosyadaki tek Paraxial yuzey)
ayar.B.kompOnSurf = 0;      % 0 = otomatik (dosyadaki tek DEGISKEN kalinlik = on dummy)
ayar.B.phiD       = linspace(-0.5, 0.5, 21);   % taranacak vergens [D]
% calistir = false ise kullanilacak degerler (onceki kosudan):
ayar.B.egim       = -1.311;          % [mm/D]
ayar.B.pencere    = [-0.05 0.05];    % [D]
ayar.B.pencereAlt = true;            % true: pencere tarama sinirina dayanmisti (gercek pencere daha genis)

% ---------- PARCA 2: Dosya A (Monte Carlo dosyalari) ----------
ayar.A.calistir   = true;
ayar.A.dosya      = 'C:\...\DOSYA_A_tolerans.zmx';  % MC'nin kosuldugu NOMINAL dosya
ayar.A.mcKlasor   = '';        % bos = Dosya A'nin klasoru
ayar.A.mcOnek     = 'MC_T';    % Zemax'in kaydettigi MC dosyalarinin on eki
ayar.A.kompOnSurf = 0;         % 0 = otomatik (TDE'deki COMP operandinin yuzeyi)

% "Yol olarak kopyala" ile gelen cift tirnaklari temizle
ayar.B.dosya = strtrim(strrep(ayar.B.dosya, '"', ''));
ayar.A.dosya = strtrim(strrep(ayar.A.dosya, '"', ''));
ayar.A.mcKlasor = strtrim(strrep(ayar.A.mcKlasor, '"', ''));
if ayar.B.calistir && ~isfile(ayar.B.dosya), error('Dosya B bulunamadi: %s', ayar.B.dosya); end
if ayar.A.calistir && ~isfile(ayar.A.dosya), error('Dosya A bulunamadi: %s', ayar.A.dosya); end

%% ===================== BAGLANTI =====================
TheApplication = baglan();
TheSystem = TheApplication.PrimarySystem;
izinStrok = ayar.strok - ayar.termal;   % toleranslara kalan strok
if izinStrok <= 0, error('ayar.termal, ayar.strok''tan kucuk olmali.'); end

%% ===================== PARCA 1: KAPASITE =====================
if ayar.B.calistir
    fprintf('\n========== PARCA 1: PARAXIAL TARAMA (Dosya B) ==========\n');
    if ~TheSystem.LoadFile(ayar.B.dosya, false)
        error('Dosya B acilamadi: %s', ayar.B.dosya);
    end
    LDE = TheSystem.LDE;

    % --- Yuzey numaralari ---
    pS = ayar.B.paraxSurf;
    if pS == 0
        bul = paraxialBul(LDE);
        if numel(bul) ~= 1
            error('Dosya B''de %d Paraxial yuzey bulundu (%s). ayar.B.paraxSurf''i elle gir.', numel(bul), mat2str(bul));
        end
        pS = bul;
    elseif ~paraxialMi(LDE.GetSurfaceAt(pS))
        error('Dosya B yuzey %d Paraxial degil. Paraxial yuzeyler: %s', pS, mat2str(paraxialBul(LDE)));
    end
    kS = ayar.B.kompOnSurf;
    if kS == 0
        bul = degiskenKalinlikBul(LDE);
        if numel(bul) ~= 1
            error(['Dosya B''de %d degisken kalinlik var (%s). Sadece kompanzatorun on dummy''si degisken olmali, ' ...
                   'ya da ayar.B.kompOnSurf''i elle gir.'], numel(bul), mat2str(bul));
        end
        kS = bul;
    end
    fprintf('Dosya B: Paraxial yuzey = %d | Kompanzator (on dummy) = %d\n', pS, kS);
    if isempty(pickupBul(LDE))
        warning('Dosya B''de pickup solve''li kalinlik bulunamadi. Kompanzatorun arka dummy''si pickup -1 olmali.');
    end

    % --- Hazirlik ---
    sP = LDE.GetSurfaceAt(pS);
    sK = LDE.GetSurfaceAt(kS);
    parCell = sP.GetSurfaceCell(ZOSAPI.Editors.LDE.SurfaceColumn.Par1);   % Focal Length
    T0 = sK.Thickness;
    TheSystem.Tools.RemoveAllVariables();
    sK.ThicknessCell.MakeSolveVariable();

    [fNyq, birim] = nyquist(ayar);
    nAlan = TheSystem.SystemData.Fields.NumberOfFields;
    phiD = unique([ayar.B.phiD(:); 0]);
    n = numel(phiD);  i0 = find(phiD == 0);
    dzB = nan(n,1);  mtfT = nan(n,nAlan);  mtfS = nan(n,nAlan);
    fprintf('Nyquist: %.3f %s | Alan: %d | Nokta: %d\n\n', fNyq, birim, nAlan, n);

    % --- Tarama ---
    for k = 1:n
        if phiD(k) == 0, f = 1e10; else, f = 1000 / phiD(k); end
        parCell.DoubleValue = f;
        sK.Thickness = T0;
        optimizeEt(TheSystem, ayar.cekirdek);
        dzB(k) = sK.Thickness - T0;
        [mtfT(k,:), mtfS(k,:)] = mtfOku(TheSystem, fNyq, nAlan);
        fprintf('phi = %+7.4f D | dz = %+8.4f mm | MTF_T(alan1) = %.3f\n', phiD(k), dzB(k), mtfT(k,1));
    end

    % --- Pencere ---
    dzB = dzB - dzB(i0);
    nT = mtfT ./ mtfT(i0,:);  nS = mtfS ./ mtfS(i0,:);
    okMTF   = all(nT >= 1 - ayar.mtfDusus, 2) & all(nS >= 1 - ayar.mtfDusus, 2);
    okStrok = abs(dzB) <= izinStrok;
    ok = okMTF & okStrok;
    lo = i0; while lo > 1 && ok(lo-1), lo = lo - 1; end
    hi = i0; while hi < n && ok(hi+1), hi = hi + 1; end

    p = polyfit(phiD, dzB, 1);
    egim = p(1);
    pencere = [phiD(lo) phiD(hi)];
    pencereAlt = (lo == 1 || hi == n);

    fprintf('\nKompanzator duyarliligi: %.4f mm/D\n', egim);
    fprintf('Kabul penceresi: %+.4f ... %+.4f D', pencere(1), pencere(2));
    if pencereAlt, fprintf('  (tarama sinirina dayandi: gercek pencere DAHA GENIS)'); end
    fprintf('\n');
    if lo > 1 || hi < n
        sebepler = {'strok', 'MTF'};
        fprintf('Pencereyi sinirlayan: %s\n', strjoin(sebepler([any(~okStrok) any(~okMTF)]), ' + '));
    end

    T1 = table(phiD, dzB, okStrok, okMTF, 'VariableNames', {'phi_D', 'dz_komp_mm', 'Strok_OK', 'MTF_OK'});
    for a = 1:nAlan
        T1.(sprintf('MTFn_T_alan%d', a)) = nT(:,a);
        T1.(sprintf('MTFn_S_alan%d', a)) = nS(:,a);
    end
    writetable(T1, [ayar.cikti '_parca1.csv']);

    figure('Name', 'Parca 1 - Kompanzator hareketi');
    plot(phiD, dzB, 'o-', 'LineWidth', 1.5); hold on; grid on;
    yline(izinStrok, 'r--', 'strok (termal dusulmus)'); yline(-izinStrok, 'r--');
    xline(pencere(1), 'g:'); xline(pencere(2), 'g:');
    xlabel('\phi [D]'); ylabel('\Deltaz_{komp} [mm]');
    title(sprintf('Kapasite: kompanzator hareketi (egim %.3f mm/D)', egim));

    figure('Name', 'Parca 1 - Normalize MTF');
    hold on; grid on; renk = lines(nAlan);
    for a = 1:nAlan
        plot(phiD, nT(:,a), '-o',  'Color', renk(a,:), 'DisplayName', sprintf('Alan %d T', a));
        plot(phiD, nS(:,a), '--s', 'Color', renk(a,:), 'DisplayName', sprintf('Alan %d S', a));
    end
    yline(1 - ayar.mtfDusus, 'r--', 'esik', 'HandleVisibility', 'off');
    xlabel('\phi [D]'); ylabel('MTF / MTF_{nominal}');
    title(sprintf('Kapasite: kompanzasyon sonrasi MTF @ %.3f %s', fNyq, birim));
    legend('Location', 'bestoutside');
else
    egim = ayar.B.egim;  pencere = ayar.B.pencere;  pencereAlt = ayar.B.pencereAlt;
    fprintf('\nParca 1 atlandi. Egim = %.4f mm/D, pencere = [%+.4f %+.4f] D\n', egim, pencere);
end

%% ===================== PARCA 2: IHTIYAC =====================
if ~ayar.A.calistir
    fprintf('\nParca 2 atlandi.\n');
    TheApplication.CloseApplication();
    return;
end
fprintf('\n========== PARCA 2: MONTE CARLO DOSYALARI (Dosya A) ==========\n');
if ~TheSystem.LoadFile(ayar.A.dosya, false)
    error('Dosya A acilamadi: %s', ayar.A.dosya);
end

kA = ayar.A.kompOnSurf;
if kA == 0
    bul = compBul(TheSystem);
    if numel(bul) ~= 1
        error(['Dosya A''nin TDE''sinde %d kalinlik kompanzatoru (COMP) bulundu (%s). ' ...
               'Tek COMP olmali, ya da ayar.A.kompOnSurf''i elle gir.'], numel(bul), mat2str(bul));
    end
    kA = bul;
end
Tnom = TheSystem.LDE.GetSurfaceAt(kA).Thickness;
fprintf('Dosya A: Kompanzator (on dummy) = %d | nominal kalinlik = %.4f mm\n', kA, Tnom);

klasor = ayar.A.mcKlasor;
if isempty(klasor), klasor = fileparts(ayar.A.dosya); end
mc = mcDosyalari(klasor, ayar.A.mcOnek);
nMC = numel(mc);
fprintf('%d MC dosyasi bulundu: %s\n\n', nMC, klasor);

dzA = nan(nMC,1);
for k = 1:nMC
    yol = fullfile(klasor, mc{k});
    if ~TheSystem.LoadFile(yol, false)
        warning('Acilamadi: %s', mc{k}); continue;
    end
    if kA > TheSystem.LDE.NumberOfSurfaces - 1
        warning('%s: yuzey sayisi farkli, atlandi.', mc{k}); continue;
    end
    dzA(k) = TheSystem.LDE.GetSurfaceAt(kA).Thickness - Tnom;
    if mod(k, 50) == 0 || k == nMC, fprintf('  %d / %d okundu\n', k, nMC); end
end
TheApplication.CloseApplication();

gecerli = ~isnan(dzA);
dzA  = dzA(gecerli);  mcAd = mc(gecerli);
phiA = dzA / egim;                 % Parca 1 ile ayni isaret kurali
absP = abs(phiA);

st.ort   = mean(phiA);
st.sigma = std(phiA);
st.p95   = yuzdelik(absP, 95);
st.max   = max(absP);
st.icinde = mean(phiA >= pencere(1) & phiA <= pencere(2)) * 100;
st.dzMax = max(abs(dzA));

fprintf('\n===== GERCEK VERGENS (%d MC kosusu) =====\n', numel(phiA));
fprintf('Ortalama       : %+.4f D\n', st.ort);
fprintf('Std sapma (1s) : %.4f D   (~%.0f m)\n', st.sigma, 1/max(st.sigma, eps));
fprintf('%%95 |phi|      : %.4f D   (~%.0f m)\n', st.p95, 1/max(st.p95, eps));
fprintf('En kotu |phi|  : %.4f D   (~%.0f m)\n', st.max, 1/max(st.max, eps));
fprintf('Kompanzator en buyuk hareket: %.4f mm (+ termal %.2f = %.4f / strok %.2f)\n', ...
        st.dzMax, ayar.termal, st.dzMax + ayar.termal, ayar.strok);
fprintf('Pencere icinde kalan kosu orani: %%%.1f\n', st.icinde);

%% ===================== KARAR =====================
fprintf('\n===== KARAR =====\n');
strokOK = st.dzMax + ayar.termal <= ayar.strok;
if st.icinde == 100 && strokOK
    fprintf('EVET: Tum MC kosularinda afokallik bozulmasi kompanzatorun kapasitesi icinde.\n');
elseif st.icinde >= 95 && strokOK
    fprintf('BUYUK OLCUDE EVET: Kosularin %%%.1f''i pencere icinde. Disaridakileri kontrol et.\n', st.icinde);
elseif pencereAlt
    fprintf(['BELIRSIZ: Bazi kosular bilinen pencerenin disinda ama pencere tarama sinirina dayanmisti.\n' ...
             'Parca 1''de ayar.B.phiD araligini genisletip tekrar calistir.\n']);
else
    fprintf(['HAYIR: Afokallik bozulmasi kompanzatorun kapasitesini asiyor (pencere icinde %%%.1f).\n' ...
             'Sensitivity raporunda en buyuk kompanzator degerini veren toleranslari sikilastir ya da stroku artir.\n'], st.icinde);
end
if ~strokOK
    fprintf('UYARI: Kompanzator hareketi + termal pay stroku asiyor.\n');
end

T2 = table(string(mcAd(:)), dzA, phiA, 'VariableNames', {'Dosya', 'dz_komp_mm', 'phi_esdeger_D'});
writetable(T2, [ayar.cikti '_parca2.csv']);
fprintf('\nCSV dosyalari yazildi: %s_parca*.csv\n', ayar.cikti);

%% ===================== GRAFIKLER (PARCA 2) =====================
figure('Name', 'Parca 2 - Gercek vergens dagilimi');
histogram(phiA, max(10, round(sqrt(numel(phiA))))); hold on; grid on;
xline(pencere(1), 'g-', 'LineWidth', 2); xline(pencere(2), 'g-', 'kabul penceresi', 'LineWidth', 2);
xline(-st.p95, 'm--'); xline(st.p95, 'm--', '%95');
xlabel('\phi_{esdeger} [D]'); ylabel('Kosu sayisi');
title(sprintf('Ihtiyac vs kapasite: %%%.1f pencere icinde', st.icinde));

figure('Name', 'Parca 2 - Kompanzator hareketi');
plot(dzA, '.', 'MarkerSize', 10); hold on; grid on;
yline(izinStrok, 'r--', 'strok (termal dusulmus)'); yline(-izinStrok, 'r--');
xlabel('MC kosusu'); ylabel('\Deltaz_{komp} [mm]');
title('Her MC kosusunda kompanzator hareketi');

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

function [fNyq, birim] = nyquist(ayar)
    if ayar.ters
        fNyq = ayar.efl / (2 * ayar.pikselUm);  birim = 'cycles/mrad';
    else
        fNyq = 1000 / (2 * ayar.pikselUm);      birim = 'lp/mm';
    end
end

function e = paraxialMi(surf)
    e = contains(lower(char(surf.TypeName)), 'paraxial') || ...
        contains(lower(char(surf.Type.ToString())), 'paraxial');
end

function n = paraxialBul(LDE)
    n = [];
    for i = 1:LDE.NumberOfSurfaces - 1
        if paraxialMi(LDE.GetSurfaceAt(i)), n(end+1) = i; end %#ok<AGROW>
    end
end

function n = degiskenKalinlikBul(LDE)
    n = [];
    for i = 0:LDE.NumberOfSurfaces - 1
        t = char(LDE.GetSurfaceAt(i).ThicknessCell.GetSolveData().Type.ToString());
        if strcmp(t, 'Variable'), n(end+1) = i; end %#ok<AGROW>
    end
end

function n = pickupBul(LDE)
    n = [];
    for i = 0:LDE.NumberOfSurfaces - 1
        t = char(LDE.GetSurfaceAt(i).ThicknessCell.GetSolveData().Type.ToString());
        if contains(t, 'Pickup'), n(end+1) = i; end %#ok<AGROW>
    end
end

function n = compBul(TheSystem)
    % TDE'deki COMP operandlarindan kalinlik kompanzatorlerini bulur (Int2/kod = 0)
    n = [];
    try
        TDE = TheSystem.TDE;
        for i = 1:TDE.NumberOfOperands
            op = TDE.GetOperandAt(i);
            if strcmpi(char(op.Type.ToString()), 'COMP')
                yuzey = op.GetOperandCell(ZOSAPI.Editors.TDE.TolerancingColumn.Param1).IntegerValue;
                kod   = op.GetOperandCell(ZOSAPI.Editors.TDE.TolerancingColumn.Param2).IntegerValue;
                if kod == 0, n(end+1) = yuzey; end %#ok<AGROW>
            end
        end
    catch hata
        warning('COMP otomatik okunamadi (%s). ayar.A.kompOnSurf''i elle gir.', hata.message);
    end
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

function adlar = mcDosyalari(klasor, onek)
    d = [dir(fullfile(klasor, [onek '*.zmx'])); dir(fullfile(klasor, [onek '*.ZMX']))];
    if isempty(d), error('"%s" ile baslayan MC dosyasi bulunamadi: %s', onek, klasor); end
    adlar = {d.name};
    adlar = adlar(startsWith(string(adlar), onek));          % Windows dir buyuk/kucuk harf duyarsiz
    [~, u] = unique(lower(string(adlar)), 'stable');
    adlar = sort(adlar(u));
end

function v = yuzdelik(x, p)
    x = sort(x(:));
    if isempty(x), v = NaN; return; end
    v = interp1(linspace(0, 100, numel(x)), x, p);
end
