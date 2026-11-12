%% AFOKAL TOLERANS ANALIZI (ZOS-API, MATLAB)
% "Afokal grubun bozulmasini kompanzatorle tolere edebilir miyim?" sorusunu
% iki parcada cevaplar. Tum sonuclar KOMPANZATOR HAREKETI [mm] cinsindendir.
%
%  PARCA 1 - KAPASITE (Dosya B, paraxial yuzeyli dosya)
%    Kolime bosluktaki Paraxial yuzeyin gucu adim adim degistirilerek yapay
%    afokallik hatasi verilir, her adimda kompanzator optimize edilir ve
%    kompanzatorun hareketi (dz) ile Nyquist/2'deki MTF okunur.
%    Cikti: kabul penceresi [mm] = kompanzatorun, hem strok hem MTF kosulunu
%    saglayarak duzeltebildigi hareket araligi. MTF kosulu: Nyquist/2'deki
%    (ayar.mtfFrek) MTF, hatasiz durumdaki degerin %5'inden (ayar.mtfDusus)
%    fazla dusmemeli. Nyquist'te MTF cok dusuk oldugu icin olcut olarak
%    kullanilmaz.
%
%  PARCA 2 - IHTIYAC (Dosya A'nin Monte Carlo dosyalari)
%    Zemax'te sadece afokal grubun ic toleranslariyla kosulan Monte Carlo'nun
%    kaydettigi MC_T*.zmx dosyalari tek tek acilir ve kompanzatorun ne kadar
%    hareket ettigi okunur (dz = MC degeri - nominal deger) [mm].
%    Kompanzator her MC dosyasinda ayri bulunur (tek DEGISKEN kalinlik):
%    toleranslama MC dosyalarina coordinate break ekledigi icin yuzey
%    numaralari nominal dosyadan farklidir.
%
%  KARAR: MC dosyalarindaki kompanzator hareketleri kabul penceresinin icinde mi?
%
% DOSYALAR: Kod calisinca bir klasor sectirir ve dosyalari oradan kendisi bulur:
%    - MC dosyalari : adi ayar.A.mcOnek ("MC_T") ile baslayan .zmx dosyalari
%    - Dosya B      : adi "MC_" ile baslamayan ve Paraxial yuzey iceren .zmx
%    - Dosya A      : adi "MC_" ile baslamayan ve Paraxial yuzey icermeyen .zmx
%                     (MC'nin kosuldugu nominal dosya)
%    Ayni turden birden fazla dosya varsa listeden sectirilir.
%    Ciktilar secilen klasorde "AFOKAL TOLERANS ANALIZI" klasorune yazilir.
%
% 0 birakilan yuzey numaralari otomatik bulunur ve ekrana yazdirilir.
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
ayar.mtfFrek   = 0.5;    % MTF dususune bakilan frekans, Nyquist'in kati (0.5 = Nyquist/2)
ayar.cekirdek  = 8;      % optimizasyon cekirdek sayisi
ayar.cikti     = 'afokal_tolerans_sonuc';   % cikti dosyalarinin on adi

% ---------- PARCA 1: Dosya B (paraxial tarama) ----------
ayar.B.calistir   = true;   % true : tarama yap, pencereyi BU tasarimdan hesapla (yeni tasarimda hep true)
                            % false: tarama yapma, onceki kosunun kaydettigi sonucu kullan
ayar.B.paraxSurf  = 0;      % 0 = otomatik (dosyadaki tek Paraxial yuzey)
ayar.B.kompOnSurf = 0;      % 0 = otomatik (dosyadaki tek DEGISKEN kalinlik = on dummy)
ayar.B.guc        = linspace(-0.5, 0.5, 21);   % Paraxial yuzeye verilecek guc adimlari [1/m]
                            % (sadece tarama girdisi: f_paraxial = 1000/guc mm. Sonuclar mm cinsindendir.
                            %  Pencere tarama sinirina dayanirsa araligi genislet.)
% calistir = false iken: Parca 1 her kosuda sonucunu cikti klasorundeki <cikti>_parca1.mat
% dosyasina kaydeder, kod onu okur. Elle vermek istersen (bos = .mat dosyasindan oku):
ayar.B.pencere    = [];     % [mm] kompanzator hareketi penceresi, ornek [-0.20 0.20]
ayar.B.pencereAlt = false;  % true: pencere tarama sinirina dayanmisti

% ---------- PARCA 2: Dosya A (Monte Carlo dosyalari) ----------
ayar.A.calistir   = true;
ayar.A.mcOnek     = 'MC_T';    % Zemax'in kaydettigi MC dosyalarinin on eki
ayar.A.kompOnSurf = 0;         % NOMINAL dosyadaki kompanzator yuzeyi; 0 = otomatik (TDE'deki COMP operandi)
% MC dosyalarinda kompanzator HER DOSYADA ayri bulunur (tek DEGISKEN kalinlik): toleranslama
% MC dosyalarina coordinate break ekledigi icin yuzey numaralari nominal dosyadan farklidir.

izinStrok = ayar.strok - ayar.termal;   % toleranslara kalan strok
if izinStrok <= 0, error('ayar.termal, ayar.strok''tan kucuk olmali.'); end

%% ===================== KLASOR SECIMI =====================
klasor = uigetdir(pwd, 'Zemax klasorunu sec (Dosya A, Dosya B ve MC dosyalari)');
if isequal(klasor, 0), return; end

ciktiKlasor = fullfile(klasor, 'AFOKAL TOLERANS ANALIZI');
if ~isfolder(ciktiKlasor), mkdir(ciktiKlasor); end
cikti = fullfile(ciktiKlasor, ayar.cikti);

% Adi "MC_" ile baslamayan .zmx dosyalari: Dosya A / Dosya B adaylari
adaylar = zmxDosyalari(klasor);
adaylar = adaylar(~startsWith(string(adaylar), 'MC_', 'IgnoreCase', true));
if isempty(adaylar)
    error('Klasorde adi "MC_" ile baslamayan .zmx dosyasi yok: %s', klasor);
end

if ayar.A.calistir
    mc = mcDosyalari(klasor, ayar.A.mcOnek);   % yoksa hata verir
end

%% ===================== BAGLANTI =====================
TheApplication = baglan();
TheSystem = TheApplication.PrimarySystem;

% Adaylari sinifla: Paraxial yuzey iceren = Dosya B, icermeyen = Dosya A
paraxVar = false(size(adaylar));
acildi = false(size(adaylar));
for i = 1:numel(adaylar)
    if ~TheSystem.LoadFile(fullfile(klasor, adaylar{i}), false)
        warning('Acilamadi, atlandi: %s', adaylar{i}); continue;
    end
    acildi(i) = true;
    paraxVar(i) = ~isempty(paraxialBul(TheSystem.LDE));
end
adaylar = adaylar(acildi);  paraxVar = paraxVar(acildi);

dosyaB = '';  dosyaA = '';
if ayar.B.calistir
    dosyaB = dosyaSec(adaylar(paraxVar), 'Dosya B (Paraxial yuzeyli dosya)', ...
        'Klasorde Paraxial yuzey iceren, adi "MC_" ile baslamayan .zmx dosyasi yok (Dosya B).');
    fprintf('Dosya B: %s\n', dosyaB);
end
if ayar.A.calistir
    dosyaA = dosyaSec(adaylar(~paraxVar), 'Dosya A (MC''nin kosuldugu nominal dosya)', ...
        'Klasorde Paraxial yuzey icermeyen, adi "MC_" ile baslamayan .zmx dosyasi yok (Dosya A).');
    fprintf('Dosya A: %s\n', dosyaA);
end

%% ===================== PARCA 1: KAPASITE =====================
if ayar.B.calistir
    fprintf('\n========== PARCA 1: PARAXIAL TARAMA (Dosya B) ==========\n');
    if ~TheSystem.LoadFile(fullfile(klasor, dosyaB), false)
        error('Dosya B acilamadi: %s', dosyaB);
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
    fOlc = ayar.mtfFrek * fNyq;   % MTF kosulunun bakildigi frekans
    nAlan = TheSystem.SystemData.Fields.NumberOfFields;
    guc = unique([ayar.B.guc(:); 0]);          % [1/m], 0 = hatasiz durum
    n = numel(guc);  i0 = find(guc == 0);
    fPar = 1000 ./ guc;  fPar(i0) = 1e10;      % Paraxial yuzeyin odak uzakligi [mm]
    dzB = nan(n,1);  mtfT = nan(n,nAlan);  mtfS = nan(n,nAlan);
    fprintf('Nyquist: %.3f %s | MTF kosulu: %.3f %s (%g x Nyquist) | Alan: %d | Nokta: %d\n\n', ...
        fNyq, birim, fOlc, birim, ayar.mtfFrek, nAlan, n);

    % --- Tarama ---
    for k = 1:n
        parCell.DoubleValue = fPar(k);
        sK.Thickness = T0;            % her adimda nominalden basla
        optimizeEt(TheSystem, ayar.cekirdek);
        dzB(k) = sK.Thickness - T0;
        [mtfT(k,:), mtfS(k,:)] = mtfOku(TheSystem, fOlc, nAlan);
        fprintf('adim %2d / %d | f_paraxial = %+12.1f mm | dz = %+8.4f mm | MTF_T(alan1) = %.3f\n', ...
            k, n, fPar(k), dzB(k), mtfT(k,1));
    end

    % --- Pencere (kompanzator hareketi cinsinden) ---
    dzB = dzB - dzB(i0);
    nT = mtfT ./ mtfT(i0,:);  nS = mtfS ./ mtfS(i0,:);
    okMTF   = all(nT >= 1 - ayar.mtfDusus, 2) & all(nS >= 1 - ayar.mtfDusus, 2);
    okStrok = abs(dzB) <= izinStrok;
    ok = okMTF & okStrok;
    lo = i0; while lo > 1 && ok(lo-1), lo = lo - 1; end
    hi = i0; while hi < n && ok(hi+1), hi = hi + 1; end

    pencere = [min(dzB(lo:hi)) max(dzB(lo:hi))];   % [mm]
    pencereAlt = (lo == 1 || hi == n);

    fprintf('\nKabul penceresi (kompanzator hareketi): %+.4f ... %+.4f mm', pencere(1), pencere(2));
    if pencereAlt, fprintf('  (tarama sinirina dayandi: gercek pencere DAHA GENIS)'); end
    fprintf('\n');
    if lo > 1 || hi < n
        sebepler = {'strok', 'MTF'};
        fprintf('Pencereyi sinirlayan: %s\n', strjoin(sebepler([any(~okStrok) any(~okMTF)]), ' + '));
    end

    T1 = table(dzB, fPar, okStrok, okMTF, 'VariableNames', {'dz_komp_mm', 'f_paraxial_mm', 'Strok_OK', 'MTF_OK'});
    for a = 1:nAlan
        T1.(sprintf('MTFn_T_alan%d', a)) = nT(:,a);
        T1.(sprintf('MTFn_S_alan%d', a)) = nS(:,a);
    end
    T1 = sortrows(T1, 'dz_komp_mm');
    writetable(T1, [cikti '_parca1.csv']);
    kaynakDosya = dosyaB;
    save([cikti '_parca1.mat'], 'pencere', 'pencereAlt', 'kaynakDosya');

    [dzSirali, sira] = sort(dzB);
    figure('Name', 'Parca 1 - Normalize MTF');
    hold on; grid on; renk = lines(nAlan);
    for a = 1:nAlan
        plot(dzSirali, nT(sira,a), '-o',  'Color', renk(a,:), 'DisplayName', sprintf('Alan %d T', a));
        plot(dzSirali, nS(sira,a), '--s', 'Color', renk(a,:), 'DisplayName', sprintf('Alan %d S', a));
    end
    yline(1 - ayar.mtfDusus, 'r--', 'MTF esigi', 'HandleVisibility', 'off');
    xline(izinStrok, 'r--', 'strok', 'HandleVisibility', 'off');
    xline(-izinStrok, 'r--', 'HandleVisibility', 'off');
    xline(pencere(1), 'g-', 'LineWidth', 2, 'HandleVisibility', 'off');
    xline(pencere(2), 'g-', 'kabul penceresi', 'LineWidth', 2, 'HandleVisibility', 'off');
    xlabel('Kompanzator hareketi \Deltaz_{komp} [mm]'); ylabel('MTF / MTF_{nominal}');
    title(sprintf('Kapasite: kompanzasyon sonrasi MTF @ %.3f %s (%g x Nyquist)', fOlc, birim, ayar.mtfFrek));
    legend('Location', 'bestoutside');
    resimKaydet(gcf, [cikti '_parca1_MTF.png']);
else
    if ~isempty(ayar.B.pencere)
        pencere = sort(ayar.B.pencere);  pencereAlt = ayar.B.pencereAlt;
    else
        matDosya = [cikti '_parca1.mat'];
        if ~isfile(matDosya)
            error(['Parca 1 atlandi ama %s yok. Bu tasarim icin once ayar.B.calistir = true ile calistir ' ...
                   '(ya da ayar.B.pencere''yi elle doldur).'], matDosya);
        end
        k1 = load(matDosya);
        if isfield(k1, 'egim')
            error(['%s eski surumle (diyoptri cinsinden) kaydedilmis. ' ...
                   'ayar.B.calistir = true ile Parca 1''i tekrar calistir.'], matDosya);
        end
        pencere = k1.pencere;  pencereAlt = k1.pencereAlt;
        fprintf('\nParca 1 sonucu okundu: %s (kaynak: %s)\n', matDosya, k1.kaynakDosya);
    end
    fprintf('\nParca 1 atlandi. Kabul penceresi = [%+.4f %+.4f] mm\n', pencere);
end

%% ===================== PARCA 2: IHTIYAC =====================
if ~ayar.A.calistir
    fprintf('\nParca 2 atlandi.\n');
    TheApplication.CloseApplication();
    return;
end
fprintf('\n========== PARCA 2: MONTE CARLO DOSYALARI (Dosya A) ==========\n');
if ~TheSystem.LoadFile(fullfile(klasor, dosyaA), false)
    error('Dosya A acilamadi: %s', dosyaA);
end

kA = ayar.A.kompOnSurf;
if kA == 0
    bul = compBul(TheSystem);
    if numel(bul) ~= 1
        error(['Dosya A''da (nominal) kompanzator otomatik bulunamadi: %d aday (%s). ' ...
               'ayar.A.kompOnSurf''a NOMINAL dosyadaki kompanzator yuzey numarasini elle gir.'], numel(bul), mat2str(bul));
    end
    kA = bul;
end
Tnom = TheSystem.LDE.GetSurfaceAt(kA).Thickness;
fprintf('Dosya A (nominal): Kompanzator (on dummy) = %d | nominal kalinlik = %.4f mm\n', kA, Tnom);

nMC = numel(mc);
fprintf('%d MC dosyasi bulundu: %s\n\n', nMC, klasor);

dzA = nan(nMC,1);
kMCilk = NaN;   % MC dosyalarinda bulunan kompanzator yuzeyi (ilk gecerli dosya)
for k = 1:nMC
    yol = fullfile(klasor, mc{k});
    if ~TheSystem.LoadFile(yol, false)
        warning('Acilamadi: %s', mc{k}); continue;
    end
    % MC dosyasinda kompanzator = tek degisken kalinlik (yuzey numarasi nominalden farkli)
    kMC = degiskenKalinlikBul(TheSystem.LDE);
    if numel(kMC) ~= 1
        warning('%s: %d degisken kalinlik bulundu (%s), tek kompanzator bekleniyordu; atlandi.', ...
            mc{k}, numel(kMC), mat2str(kMC));
        continue;
    end
    if isnan(kMCilk), kMCilk = kMC; end
    dzA(k) = TheSystem.LDE.GetSurfaceAt(kMC).Thickness - Tnom;
    if mod(k, 50) == 0 || k == nMC, fprintf('  %d / %d okundu\n', k, nMC); end
end
TheApplication.CloseApplication();
fprintf('MC dosyalarinda kompanzator yuzeyi (tek degisken kalinlik) = %d\n', kMCilk);

gecerli = ~isnan(dzA);
if ~any(gecerli)
    error('Hicbir MC dosyasinda tek degisken kalinlik (kompanzator) bulunamadi.');
end
dzA  = dzA(gecerli);  mcAd = mc(gecerli);
icinde = dzA >= pencere(1) & dzA <= pencere(2);   % Parca 1 ile ayni isaret kurali
absDz = abs(dzA);

st.ort    = mean(dzA);
st.sigma  = std(dzA);
st.p95    = yuzdelik(absDz, 95);
st.max    = max(absDz);
st.icinde = mean(icinde) * 100;

fprintf('\n===== KOMPANZATOR HAREKETI (%d MC kosusu) =====\n', numel(dzA));
fprintf('Kabul penceresi   : %+.4f ... %+.4f mm\n', pencere(1), pencere(2));
fprintf('Ortalama          : %+.4f mm\n', st.ort);
fprintf('Std sapma (1s)    : %.4f mm\n', st.sigma);
fprintf('%%95 |dz|          : %.4f mm\n', st.p95);
fprintf('En buyuk |dz|     : %.4f mm (+ termal %.2f = %.4f / strok %.2f)\n', ...
        st.max, ayar.termal, st.max + ayar.termal, ayar.strok);
fprintf('Pencere icinde kalan kosu orani: %%%.1f\n', st.icinde);

%% ===================== KARAR =====================
fprintf('\n===== KARAR =====\n');
strokOK = st.max + ayar.termal <= ayar.strok;
if st.icinde == 100 && strokOK
    fprintf('EVET: Tum MC kosularinda kompanzator hareketi kabul penceresinin icinde.\n');
elseif st.icinde >= 95 && strokOK
    fprintf('BUYUK OLCUDE EVET: Kosularin %%%.1f''i pencere icinde. Disaridakileri kontrol et.\n', st.icinde);
elseif pencereAlt
    fprintf(['BELIRSIZ: Bazi kosular bilinen pencerenin disinda ama pencere tarama sinirina dayanmisti.\n' ...
             'Parca 1''de ayar.B.guc araligini genisletip tekrar calistir.\n']);
else
    fprintf(['HAYIR: Gereken kompanzator hareketi kompanzatorun kapasitesini asiyor (pencere icinde %%%.1f).\n' ...
             'Sensitivity raporunda en buyuk kompanzator degerini veren toleranslari sikilastir ya da stroku artir.\n'], st.icinde);
end
if ~strokOK
    fprintf('UYARI: Kompanzator hareketi + termal pay stroku asiyor.\n');
end

T2 = table(string(mcAd(:)), dzA, icinde, 'VariableNames', {'Dosya', 'dz_komp_mm', 'Pencere_Icinde'});
writetable(T2, [cikti '_parca2.csv']);
fprintf('\nCiktilar yazildi: %s\n', ciktiKlasor);

%% ===================== GRAFIKLER (PARCA 2) =====================
figure('Name', 'Parca 2 - Kompanzator hareketi dagilimi');
histogram(dzA, max(10, round(sqrt(numel(dzA))))); hold on; grid on;
xline(pencere(1), 'g-', 'LineWidth', 2); xline(pencere(2), 'g-', 'kabul penceresi', 'LineWidth', 2);
xline(-st.p95, 'm--'); xline(st.p95, 'm--', '%95');
xlabel('Kompanzator hareketi \Deltaz_{komp} [mm]'); ylabel('Kosu sayisi');
title(sprintf('Ihtiyac vs kapasite: %%%.1f pencere icinde', st.icinde));
resimKaydet(gcf, [cikti '_parca2_dagilim.png']);

figure('Name', 'Parca 2 - Kompanzator hareketi');
plot(dzA, '.', 'MarkerSize', 10); hold on; grid on;
yline(izinStrok, 'r--', 'strok (termal dusulmus)'); yline(-izinStrok, 'r--');
yline(pencere(1), 'g-', 'LineWidth', 1.5); yline(pencere(2), 'g-', 'kabul penceresi', 'LineWidth', 1.5);
xlabel('MC kosusu'); ylabel('\Deltaz_{komp} [mm]');
title('Her MC kosusunda kompanzator hareketi');
resimKaydet(gcf, [cikti '_parca2_hareket.png']);

%% ===================== YARDIMCI FONKSIYONLAR =====================
function adlar = zmxDosyalari(klasor)
    % Klasordeki tum .zmx dosyalari (ayni dosya iki kez sayilmaz)
    d = [dir(fullfile(klasor, '*.zmx')); dir(fullfile(klasor, '*.ZMX'))];
    adlar = {d.name};
    [~, u] = unique(lower(string(adlar)), 'stable');
    adlar = sort(adlar(u));
end

function ad = dosyaSec(adaylar, baslik, yokMesaji)
    % Tek aday varsa onu, birden fazlaysa listeden secileni dondurur
    if isempty(adaylar)
        error('%s', yokMesaji);
    elseif isscalar(adaylar)
        ad = adaylar{1};
    else
        [sec, tamam] = listdlg('ListString', adaylar, 'SelectionMode', 'single', ...
            'PromptString', [baslik ' sec:'], 'ListSize', [450 200]);
        if ~tamam, error('Dosya secilmedi: %s', baslik); end
        ad = adaylar{sec};
    end
end

function resimKaydet(fig, dosya)
    try
        exportgraphics(fig, dosya, 'Resolution', 300);   % R2020a+
    catch
        saveas(fig, dosya);
    end
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
    % TDE'deki COMP operandlarindan kalinlik kompanzatorlerini bulur (kod = 0).
    % Bu numara sadece NOMINAL dosya icin gecerlidir; MC dosyalarinda kompanzator
    % her dosyada tek degisken kalinlik olarak ayri bulunur.
    n = [];
    try
        TDE = TheSystem.TDE;
        for i = 1:TDE.NumberOfOperands
            op = TDE.GetOperandAt(i);
            if strcmpi(char(op.Type.ToString()), 'COMP')
                if double(op.Param2) == 0, n(end+1) = double(op.Param1); end %#ok<AGROW>
            end
        end
    catch hata
        warning('COMP otomatik okunamadi (%s).', hata.message);
        n = [];
    end
    % Yedek: COMP okunamadiysa nominal dosyadaki tek degisken kalinlik
    if isempty(n)
        d = degiskenKalinlikBul(TheSystem.LDE);
        if isscalar(d)
            n = d;
            fprintf('COMP operandi okunamadi; nominal dosyadaki tek degisken kalinlik kullanildi (yuzey %d).\n', d);
        end
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
