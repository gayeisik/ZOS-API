%% AFOKAL TOLERANS ANALIZI - V2 (ZOS-API, MATLAB)
% "Afokal grubun (extender) bozulmasini objektifteki kompanzatorle tolere
% edebilir miyim?" sorusunu uc ayri kontrolle cevaplar. Kompanzatorle ilgili
% tum sonuclar KOMPANZATOR HAREKETI [mm] cinsindendir.
%
%  PARCA 1 - KAPASITE (Dosya B, paraxial yuzeyli dosya)
%    Kolime bosluktaki Paraxial yuzeyin gucu adim adim degistirilerek yapay
%    afokallik hatasi (saf odak hatasi) verilir. Her adimda kompanzator
%    optimize edilir; kompanzatorun hareketi (dz) ve MTF okunur.
%    Cikti: kabul penceresi [mm] = kompanzatorun, hem strok hem MTF kosulunu
%    saglayarak duzeltebildigi hareket araligi.
%
%  PARCA 2 - IHTIYAC VE SONUC (Dosya A'nin Monte Carlo dosyalari)
%    Her MC_T*.zmx dosyasi acilir ve iki sey okunur:
%      - kompanzatorun hareketi: dz = MC degeri - nominal deger [mm]
%      - dosyanin GERCEK MTF'i (kompanzasyon sonrasi kalan performans)
%    Kompanzator her MC dosyasinda ayri bulunur (tek DEGISKEN kalinlik):
%    toleranslama MC dosyalarina coordinate break ekledigi icin yuzey
%    numaralari nominal dosyadan farklidir.
%
%  KARAR - uc kontrol:
%    1) STROK    : gereken hareket + termal pay stroka sigiyor mu?
%    2) KAPASITE : gereken hareket Parca 1 penceresinin icinde mi? (kompanzator
%                  o kadar hareket ettiginde kendi basina MTF'i bozmuyor mu?)
%    3) MTF      : MC dosyalarinin kompanzasyon sonrasi MTF'i sarti sagliyor mu?
%    1 ve 2 gerekli kosullardir; asil sonucu 3 verir. Kompanzator sadece odak
%    hatasini duzeltir; toleranslarin urettigi koma/astigmatizma 3'te gorunur.
%
%  V1'DEN FARKLARI
%    - MTF'e tek frekans yerine birden fazla frekansta bakilir (varsayilan
%      Nyquist/2 ve Nyquist). Nyquist'te MTF dusuk oldugu icin tek basina
%      guvenilir bir olcut degildir.
%    - Parca 1 MTF kosuluna sadece Nyquist/2'de bakilir (ayar.B.kosulFrek);
%      Nyquist sadece raporlanir. Izin verilen dusus = max(goreli dusus x
%      nominal MTF, mutlak alt sinir). Mutlak alt sinir, dusuk MTF'te hesap
%      gurultusunun pencereyi belirlemesini onler.
%    - Parca 2'de MC dosyalarinin gercek MTF'i hesaplanir ve karara girer.
%    - Karar uc kontrolu ayri ayri raporlar.
%
% ONEMLI: Sonuc, Monte Carlo'ya hangi toleranslari koyduguna baglidir. Sadece
%    extender'in etkisini gormek icin MC sadece extender toleranslariyla
%    kosulmus olmalidir; tum sistem toleranslariyla kosulduysa kompanzator
%    hareketi ve MTF objektifin hatalarini da icerir.
%
% DOSYALAR: Kod calisinca bir klasor sectirir ve dosyalari oradan kendisi bulur:
%    - MC dosyalari : adi ayar.A.mcOnek ("MC_T") ile baslayan .zmx dosyalari
%    - Dosya B      : adi "MC_" ile baslamayan ve Paraxial yuzey iceren .zmx
%    - Dosya A      : adi "MC_" ile baslamayan ve Paraxial yuzey icermeyen .zmx
%                     (MC'nin kosuldugu nominal dosya)
%    Ayni turden birden fazla dosya varsa listeden sectirilir.
%    Ciktilar secilen klasorde "AFOKAL TOLERANS ANALIZI V2" klasorune yazilir.
%
% DOSYA B HAZIRLIGI: Kolime bosluga Paraxial yuzey eklenmis (f = 1e10 mm).
%    Kompanzatorun on kalinligi degisken, arka kalinligi pickup -1. Baska
%    degisken yok. Merit function hazir (optimizasyon onunla yapilir).
%
% Not: Turkce karakter kullanilmadi (eski MATLAB surumlerinde kodlama sorunu olmasin diye).

clc; clear; close all;

%% ===================== KULLANICI AYARLARI =====================
% ---------- Genel ----------
ayar.ters        = true;   % true: ters tasarim (Afocal Image Space, MTF cycles/mrad)
ayar.efl         = 239.5;  % [mm] tum sistemin EFL'si (sadece ters=true icin)
ayar.pikselUm    = 15;     % dedektor piksel araligi [um]
ayar.strok       = 0.50;   % [mm] kompanzatorun toplam mekanik stroku, +/-
ayar.termal      = 0.00;   % [mm] bu stroktan termal duzeltme icin ayrilan pay, +/-
ayar.mtfFrek     = [0.5 1.0];   % MTF'e bakilacak frekanslar (Nyquist'in kati)
ayar.cekirdek    = 8;      % optimizasyon cekirdek sayisi
ayar.cikti       = 'afokal_tolerans_v2';   % cikti dosyalarinin on adi

% ---------- PARCA 1: Dosya B (paraxial tarama) ----------
ayar.B.calistir    = true;   % true : tarama yap, pencereyi BU tasarimdan hesapla
                             % false: tarama yapma, onceki kosunun kaydettigi sonucu kullan
ayar.B.paraxSurf   = 0;      % 0 = otomatik (dosyadaki tek Paraxial yuzey)
ayar.B.kompOnSurf  = 0;      % 0 = otomatik (dosyadaki tek DEGISKEN kalinlik = on dummy)
ayar.B.guc         = linspace(-0.5, 0.5, 21);   % Paraxial yuzeye verilecek guc adimlari [1/m]
                             % (sadece tarama girdisi: f_paraxial = 1000/guc mm. Sonuclar mm cinsindendir.
                             %  Pencere tarama sinirina dayanirsa araligi genislet.)
ayar.B.kosulFrek   = 0.5;    % MTF dusus kosulunun bakildigi frekans(lar), Nyquist'in kati (0.5 = Nyquist/2).
                             % ayar.mtfFrek icinde olmali. Diger frekanslar sadece raporlanir.
ayar.B.mtfDusus    = 0.05;   % kabul edilen GORELI MTF dususu (0.05 = nominal MTF'in %5'i)
ayar.B.mtfDususMin = 0.01;   % kabul edilen en kucuk MUTLAK MTF dususu (gurultuye karsi alt sinir)
% calistir = false iken: Parca 1 her kosuda sonucunu cikti klasorundeki <cikti>_parca1.mat
% dosyasina kaydeder, kod onu okur. Elle vermek istersen (bos = .mat dosyasindan oku):
ayar.B.pencere     = [];     % [mm] kompanzator hareketi penceresi, ornek [-0.20 0.20]
ayar.B.pencereAlt  = false;  % true: pencere tarama sinirina dayanmisti

% ---------- PARCA 2: Dosya A (Monte Carlo dosyalari) ----------
ayar.A.calistir    = true;
ayar.A.mcOnek      = 'MC_T';    % Zemax'in kaydettigi MC dosyalarinin on eki
ayar.A.kompOnSurf  = 0;         % NOMINAL dosyadaki kompanzator yuzeyi; 0 = otomatik (TDE'deki COMP operandi)
ayar.A.mtfHesapla  = true;      % true: her MC dosyasinin MTF'ini hesapla (kontrol 3). false: sadece hareket (hizli)
ayar.A.mtfSart     = [];        % her frekans icin en dusuk kabul edilen MTF (en kotu alan, T/S),
                                % ayar.mtfFrek ile ayni sirada, ornek [0.30 0.10]. [] = sart yok, sadece raporla
ayar.A.verimHedef  = 90;        % [%] sartlari saglamasi istenen MC dosyasi orani

izinStrok = ayar.strok - ayar.termal;   % toleranslara kalan strok
if izinStrok <= 0, error('ayar.termal, ayar.strok''tan kucuk olmali.'); end
nFrek = numel(ayar.mtfFrek);
kosulVar = ismembertol(ayar.mtfFrek(:), ayar.B.kosulFrek(:), 1e-9);   % kosula giren frekanslar
if ~any(kosulVar) || nnz(kosulVar) ~= numel(unique(ayar.B.kosulFrek))
    error('ayar.B.kosulFrek icindeki her frekans ayar.mtfFrek icinde de olmali.');
end
if ~isempty(ayar.A.mtfSart) && numel(ayar.A.mtfSart) ~= nFrek
    error('ayar.A.mtfSart, ayar.mtfFrek ile ayni sayida deger icermeli (%d).', nFrek);
end

%% ===================== KLASOR SECIMI =====================
klasor = uigetdir(pwd, 'Zemax klasorunu sec (Dosya A, Dosya B ve MC dosyalari)');
if isequal(klasor, 0), return; end

ciktiKlasor = fullfile(klasor, 'AFOKAL TOLERANS ANALIZI V2');
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

[fNyq, birim] = nyquist(ayar);
frek = ayar.mtfFrek(:)' * fNyq;                       % MTF'e bakilan frekanslar
frekAd = compose("%.2g x Nyquist", ayar.mtfFrek(:)');  % grafik ve tablo etiketleri
fprintf('Nyquist: %.3f %s | MTF frekanslari: %s\n', fNyq, birim, strjoin(frekAd, ', '));

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

    nAlan = TheSystem.SystemData.Fields.NumberOfFields;
    guc = unique([ayar.B.guc(:); 0]);          % [1/m], 0 = hatasiz durum
    n = numel(guc);  i0 = find(guc == 0);
    fPar = 1000 ./ guc;  fPar(i0) = 1e10;      % Paraxial yuzeyin odak uzakligi [mm]
    dzB = nan(n,1);
    mtfT = nan(n, nFrek, nAlan);  mtfS = nan(n, nFrek, nAlan);
    fprintf('Alan: %d | Tarama noktasi: %d\n\n', nAlan, n);

    % --- Tarama ---
    for k = 1:n
        parCell.DoubleValue = fPar(k);
        sK.Thickness = T0;            % her adimda nominalden basla
        optimizeEt(TheSystem, ayar.cekirdek);
        dzB(k) = sK.Thickness - T0;
        [mT, mS] = mtfOkuCok(TheSystem, frek, nAlan);
        mtfT(k,:,:) = mT;  mtfS(k,:,:) = mS;
        fprintf('adim %2d / %d | f_paraxial = %+12.1f mm | dz = %+8.4f mm | en kotu MTF: %s\n', ...
            k, n, fPar(k), dzB(k), strjoin(compose("%.3f", min([mT mS], [], 2)'), ' / '));
    end

    % --- Pencere (kompanzator hareketi cinsinden) ---
    dzB = dzB - dzB(i0);
    nomT = reshape(mtfT(i0,:,:), nFrek, nAlan);
    nomS = reshape(mtfS(i0,:,:), nFrek, nAlan);
    % Izin verilen dusus: goreli dusus, ama en az mutlak alt sinir kadar
    izinT = max(ayar.B.mtfDusus * nomT, ayar.B.mtfDususMin);
    izinS = max(ayar.B.mtfDusus * nomS, ayar.B.mtfDususMin);

    okMTF = false(n,1);
    mtfMin = nan(n, nFrek);                    % en kotu alan ve yon (T/S)
    for k = 1:n
        vT = reshape(mtfT(k,:,:), nFrek, nAlan);
        vS = reshape(mtfS(k,:,:), nFrek, nAlan);
        % Kosul sadece ayar.B.kosulFrek frekanslarinda aranir (varsayilan Nyquist/2)
        dT = nomT(kosulVar,:) - vT(kosulVar,:);
        dS = nomS(kosulVar,:) - vS(kosulVar,:);
        okMTF(k) = all(dT <= izinT(kosulVar,:), 'all') && all(dS <= izinS(kosulVar,:), 'all');
        mtfMin(k,:) = min([vT vS], [], 2)';
    end
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
    for q = 1:nFrek
        T1.(sprintf('MTFmin_%s', frekEtiket(ayar.mtfFrek(q)))) = mtfMin(:,q);
    end
    for q = 1:nFrek
        for a = 1:nAlan
            T1.(sprintf('MTF_T_alan%d_%s', a, frekEtiket(ayar.mtfFrek(q)))) = mtfT(:,q,a);
            T1.(sprintf('MTF_S_alan%d_%s', a, frekEtiket(ayar.mtfFrek(q)))) = mtfS(:,q,a);
        end
    end
    T1 = sortrows(T1, 'dz_komp_mm');
    writetable(T1, [cikti '_parca1.csv']);
    kaynakDosya = dosyaB;  surum = 2;
    save([cikti '_parca1.mat'], 'pencere', 'pencereAlt', 'kaynakDosya', 'surum');

    [dzSirali, sira] = sort(dzB);
    figure('Name', 'Parca 1 - Kompanzasyon sonrasi MTF');
    hold on; grid on; renk = lines(nFrek);
    for q = 1:nFrek
        plot(dzSirali, mtfMin(sira,q), '-o', 'Color', renk(q,:), 'LineWidth', 1.5, ...
            'DisplayName', char(frekAd(q)));
    end
    xline(izinStrok, 'r--', 'strok', 'HandleVisibility', 'off');
    xline(-izinStrok, 'r--', 'HandleVisibility', 'off');
    xline(pencere(1), 'g-', 'LineWidth', 2, 'HandleVisibility', 'off');
    xline(pencere(2), 'g-', 'kabul penceresi', 'LineWidth', 2, 'HandleVisibility', 'off');
    xlabel('Kompanzator hareketi \Deltaz_{komp} [mm]'); ylabel('MTF (en kotu alan, T/S)');
    title('Kapasite: saf odak hatasi kompanze edildikten sonra kalan MTF');
    legend('Location', 'best');
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
        pencere = k1.pencere;  pencereAlt = k1.pencereAlt;
        fprintf('\nParca 1 sonucu okundu: %s (kaynak: %s)\n', matDosya, k1.kaynakDosya);
    end
    fprintf('\nParca 1 atlandi. Kabul penceresi = [%+.4f %+.4f] mm\n', pencere);
end

%% ===================== PARCA 2: IHTIYAC VE SONUC =====================
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
nAlanA = TheSystem.SystemData.Fields.NumberOfFields;
fprintf('Dosya A (nominal): Kompanzator (on dummy) = %d | nominal kalinlik = %.4f mm\n', kA, Tnom);

% Nominal dosyanin MTF'i (karsilastirma icin)
nomMin = nan(1, nFrek);
if ayar.A.mtfHesapla
    [mT, mS] = mtfOkuCok(TheSystem, frek, nAlanA);
    nomMin = min([mT mS], [], 2)';
    fprintf('Nominal MTF (en kotu alan, T/S): %s\n', ...
        strjoin(compose("%s = %.3f", frekAd(:), nomMin(:))', ' | '));
end

nMC = numel(mc);
fprintf('%d MC dosyasi bulundu: %s\n\n', nMC, klasor);

dzA = nan(nMC,1);
mcMin = nan(nMC, nFrek);   % her MC dosyasinin MTF'i (en kotu alan, T/S)
kMCilk = NaN;              % MC dosyalarinda bulunan kompanzator yuzeyi (ilk gecerli dosya)
tBas = tic;
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

    if ayar.A.mtfHesapla
        [mT, mS] = mtfOkuCok(TheSystem, frek, nAlanA);
        mcMin(k,:) = min([mT mS], [], 2)';
    end
    if mod(k, 25) == 0 || k == nMC, fprintf('  %d / %d tamamlandi (%.0f s)\n', k, nMC, toc(tBas)); end
end
TheApplication.CloseApplication();
fprintf('MC dosyalarinda kompanzator yuzeyi (tek degisken kalinlik) = %d\n', kMCilk);

gecerli = ~isnan(dzA);
if ayar.A.mtfHesapla, gecerli = gecerli & all(~isnan(mcMin), 2); end
if ~any(gecerli)
    error('Hicbir MC dosyasindan gecerli sonuc alinamadi.');
end
dzA = dzA(gecerli);  mcMin = mcMin(gecerli,:);  mcAd = mc(gecerli);
nGec = numel(dzA);

% --- Dosya bazinda uc kontrol ---
okStrokA  = abs(dzA) + ayar.termal <= ayar.strok;
okPencere = dzA >= pencere(1) & dzA <= pencere(2);   % Parca 1 ile ayni isaret kurali
sartVar   = ayar.A.mtfHesapla && ~isempty(ayar.A.mtfSart);
if sartVar
    okMtfA = all(mcMin >= ayar.A.mtfSart(:)', 2);
else
    okMtfA = true(nGec,1);
end

st.dzOrt   = mean(dzA);
st.dzSigma = std(dzA);
st.dzP95   = yuzdelik(abs(dzA), 95);
st.dzMax   = max(abs(dzA));
st.strok   = mean(okStrokA) * 100;
st.pencere = mean(okPencere) * 100;
st.mtf     = mean(okMtfA) * 100;
st.hepsi   = mean(okStrokA & okPencere & okMtfA) * 100;

fprintf('\n===== KOMPANZATOR HAREKETI (%d MC kosusu) =====\n', nGec);
fprintf('Ortalama          : %+.4f mm\n', st.dzOrt);
fprintf('Std sapma (1s)    : %.4f mm\n', st.dzSigma);
fprintf('%%95 |dz|          : %.4f mm\n', st.dzP95);
fprintf('En buyuk |dz|     : %.4f mm (+ termal %.2f = %.4f / strok %.2f)\n', ...
        st.dzMax, ayar.termal, st.dzMax + ayar.termal, ayar.strok);
fprintf('Kabul penceresi   : %+.4f ... %+.4f mm\n', pencere(1), pencere(2));

if ayar.A.mtfHesapla
    fprintf('\n===== KOMPANZASYON SONRASI MTF (en kotu alan, T/S) =====\n');
    fprintf('%-18s %10s %10s %10s %10s %10s\n', 'Frekans', 'Nominal', 'Ortalama', '%90 verim', 'En kotu', 'Sart');
    for q = 1:nFrek
        if sartVar, sartMetin = sprintf('%.3f', ayar.A.mtfSart(q)); else, sartMetin = '-'; end
        fprintf('%-18s %10.3f %10.3f %10.3f %10.3f %10s\n', char(frekAd(q)), nomMin(q), ...
            mean(mcMin(:,q)), verimDegeri(mcMin(:,q), 90), min(mcMin(:,q)), sartMetin);
    end
    fprintf('(%%90 verim: dosyalarin en az %%90''inin sagladigi MTF)\n');
end

%% ===================== KARAR =====================
fprintf('\n===== KARAR (hedef: dosyalarin en az %%%g''i) =====\n', ayar.A.verimHedef);
k1ok = st.strok >= ayar.A.verimHedef;
k2ok = st.pencere >= ayar.A.verimHedef;
fprintf('1) STROK    : dosyalarin %%%.1f''inde hareket + termal pay stroka sigiyor            -> %s\n', ...
    st.strok, evetHayir(k1ok));
fprintf('2) KAPASITE : dosyalarin %%%.1f''inde hareket Parca 1 penceresinin icinde           -> %s\n', ...
    st.pencere, evetHayir(k2ok));
if pencereAlt && ~k2ok
    fprintf('              (pencere tarama sinirina dayanmisti: ayar.B.guc araligini genisletip tekrar dene)\n');
end
if sartVar
    k3ok = st.mtf >= ayar.A.verimHedef;
    fprintf('3) MTF      : dosyalarin %%%.1f''inde kompanzasyon sonrasi MTF sarti sagliyor         -> %s\n', ...
        st.mtf, evetHayir(k3ok));
    fprintf('   Uc kontrolu birden saglayan dosya orani: %%%.1f\n', st.hepsi);
    fprintf('\nSONUC: ');
    if k1ok && k2ok && k3ok
        fprintf('EVET. Kompanzator hareketi yeterli ve kompanzasyon sonrasi MTF sarti sagliyor.\n');
    elseif k1ok && k2ok
        fprintf(['KISMEN. Kompanzator odak hatasini duzeltmeye yetiyor, ama kompanzasyon sonrasi MTF sarti\n' ...
                 '       saglamiyor. Kayip, kompanzatorun duzeltemedigi aberasyonlardan (koma, astigmatizma,\n' ...
                 '       kuresel) geliyor: toleranslari sikilastir ya da baska bir ayar (yanal ayar, clocking) dusun.\n']);
    else
        fprintf(['HAYIR. Gereken kompanzator hareketi stroku ya da kompanzatorun kapasitesini asiyor:\n' ...
                 '       en buyuk kompanzator hareketini veren toleranslari sikilastir ya da stroku artir.\n']);
    end
else
    fprintf('3) MTF      : sart girilmedi (ayar.A.mtfSart bos)');
    if ayar.A.mtfHesapla
        fprintf('; yukaridaki MTF tablosunu kendi sartinla karsilastir.\n');
    else
        fprintf(' ve MTF hesaplanmadi (ayar.A.mtfHesapla = false).\n');
    end
    fprintf('\nSONUC: ');
    if k1ok && k2ok
        fprintf(['Kompanzator hareketi yeterli (kontrol 1 ve 2). Bu GEREKLI kosuldur, YETERLI degildir:\n' ...
                 '       kompanzasyonun basarili oldugunu soylemek icin MTF sartina (kontrol 3) bakilmali.\n']);
    else
        fprintf(['HAYIR. Gereken kompanzator hareketi stroku ya da kompanzatorun kapasitesini asiyor:\n' ...
                 '       en buyuk kompanzator hareketini veren toleranslari sikilastir ya da stroku artir.\n']);
    end
end

T2 = table(string(mcAd(:)), dzA, okStrokA, okPencere, ...
    'VariableNames', {'Dosya', 'dz_komp_mm', 'Strok_OK', 'Pencere_Icinde'});
if ayar.A.mtfHesapla
    for q = 1:nFrek
        T2.(sprintf('MTFmin_%s', frekEtiket(ayar.mtfFrek(q)))) = mcMin(:,q);
    end
    if sartVar, T2.MTF_OK = okMtfA; end
end
writetable(T2, [cikti '_parca2.csv']);
fprintf('\nCiktilar yazildi: %s\n', ciktiKlasor);

%% ===================== GRAFIKLER (PARCA 2) =====================
figure('Name', 'Parca 2 - Kompanzator hareketi dagilimi');
histogram(dzA, max(10, round(sqrt(nGec)))); hold on; grid on;
xline(pencere(1), 'g-', 'LineWidth', 2); xline(pencere(2), 'g-', 'kabul penceresi', 'LineWidth', 2);
xline(-izinStrok, 'r--'); xline(izinStrok, 'r--', 'strok (termal dusulmus)');
xlabel('Kompanzator hareketi \Deltaz_{komp} [mm]'); ylabel('Kosu sayisi');
title(sprintf('Ihtiyac vs kapasite: %%%.1f pencere icinde', st.pencere));
resimKaydet(gcf, [cikti '_parca2_hareket_dagilim.png']);

if ayar.A.mtfHesapla
    figure('Name', 'Parca 2 - Kompanzasyon sonrasi MTF dagilimi', 'Position', [100 100 450 * nFrek 420]);
    for q = 1:nFrek
        subplot(1, nFrek, q);
        histogram(mcMin(:,q), max(10, round(sqrt(nGec)))); hold on; grid on;
        xline(nomMin(q), 'k-', 'nominal', 'LineWidth', 1.5);
        if sartVar, xline(ayar.A.mtfSart(q), 'r--', 'sart', 'LineWidth', 1.5); end
        xlabel('MTF (en kotu alan, T/S)'); ylabel('Kosu sayisi');
        title(char(frekAd(q)));
    end
    resimKaydet(gcf, [cikti '_parca2_MTF_dagilim.png']);

    % MTF kaybi kompanzator hareketiyle birlikte mi artiyor?
    % Evetse kayip odak turundendir; dagilim rastgeleyse baska aberasyonlardan gelir.
    figure('Name', 'Parca 2 - MTF ve kompanzator hareketi', 'Position', [100 100 450 * nFrek 420]);
    for q = 1:nFrek
        subplot(1, nFrek, q);
        plot(dzA, mcMin(:,q), '.', 'MarkerSize', 10); hold on; grid on;
        yline(nomMin(q), 'k-', 'nominal');
        if sartVar, yline(ayar.A.mtfSart(q), 'r--', 'sart'); end
        xline(pencere(1), 'g-'); xline(pencere(2), 'g-', 'pencere');
        xlabel('Kompanzator hareketi \Deltaz_{komp} [mm]'); ylabel('MTF (en kotu alan, T/S)');
        title(char(frekAd(q)));
    end
    resimKaydet(gcf, [cikti '_parca2_MTF_hareket.png']);
end

%% ===================== YARDIMCI FONKSIYONLAR =====================
function [T, S] = mtfOkuCok(TheSystem, frek, nAlan)
    % FFT MTF'i verilen frekanslarda okur. T, S: (frekans sayisi) x (alan sayisi)
    a = TheSystem.Analyses.New_FftMtf();
    s = a.GetSettings();
    s.MaximumFrequency = 1.5 * max(frek);
    s.SampleSize = ZOSAPI.Analysis.SampleSizes.S_64x64;
    try, s.ShowDiffractionLimit = false; catch, end
    a.ApplyAndWaitForCompletion();
    r = a.GetResults();
    T = nan(numel(frek), nAlan); S = T;
    for i = 0:min(r.NumberOfDataSeries, nAlan) - 1
        ds = r.GetDataSeries(i);
        x = ds.XData.Data.double;
        y = ds.YData.Data.double;     % sutun 1: tanjansiyel, sutun 2: sagital
        T(:, i+1) = interp1(x, y(:,1), frek(:));
        S(:, i+1) = interp1(x, y(:,2), frek(:));
    end
    a.Close();
end

function ad = frekEtiket(kat)
    % Sutun adi icin frekans etiketi: 0.5 -> 'N050', 1 -> 'N100'
    ad = sprintf('N%03d', round(100 * kat));
end

function v = verimDegeri(x, yuzde)
    % Degerlerin en az %yuzde'sinin sagladigi (esit veya ustunde kaldigi) deger
    x = sort(x(:), 'descend');
    v = x(max(1, ceil(yuzde / 100 * numel(x))));
end

function m = evetHayir(d)
    if d, m = 'EVET'; else, m = 'HAYIR'; end
end

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
