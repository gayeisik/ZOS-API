%% MONTE CARLO IRR DAĞILIMI - YÜZEY YÜZEY (TEZI / TIRR), ÖLÇÜMLE KARŞILAŞTIRMA
% NE YAPAR
%   OpticStudio tolerans sonuç penceresinden kaydedilen .txt dosyalarını
%   (Monte Carlo tablosu) okur. Her yüzey için TEZI ve TIRR operandlarının
%   200 MC dosyasında ürettiği düzensizliği interferometre büyüklüklerine
%   (IRR P-V, IRR RMS, ASTI; fringe) çevirir, dağılımlarını çizer ve ölçüm
%   Excel'indeki değerlerle karşılaştırır. Farklı istatistiklerle koşulmuş
%   birden fazla .txt aynı anda seçilebilir; her biri ayrı seri olarak çizilir.
%   Amaç, ölçümü en iyi taklit eden operand + istatistik kombinasyonunu bulmak.
%
% NASIL ÇALIŞIR
%   1. Bir veya birden fazla sonuç .txt dosyası seçilir (UTF-16 / UTF-8,
%      sekmeyle ayrılmış, ondalık virgül veya nokta).
%   2. Sütun adları satırında TEZI / TIRR sütunları bulunur. TEZI sütununu
%      izleyen "Zernike n" sütunları o TEZI'nin MC katsayılarıdır. Yüzey
%      numarası, Min ve Max başlık bilgi satırlarından ("Surf: 4") okunur.
%      Satırların çoğu sayı olanlar MC veri satırıdır (MC_1, MC_2, ...).
%   3. Ölçüm Excel'i ve sayfası seçilir (IRR_Olcum_Oku.m). Değerler nm
%      cinsindendir. Her yüzeyin MC dağılımı yalnızca "SURFACE (ZEMAX)"
%      sütununda aynı numarayı taşıyan ölçümlerle karşılaştırılır. Excel'de
%      olmayan yüzeyler analize alınmaz.
%   4. Her yüzey için şekil, tüm yüzeyler için birleştirilmiş şekil, Excel ve
%      .mat çıktısı üretilir.
%
% YÜZEY NUMARALARI
%   .txt'deki "Surf:" numaraları Tolerance Data Editor'daki numaralardır, yani
%   NOMİNAL dosyanın yüzey numaralarıdır. MC_T*.zmx dosyalarında koordinat
%   kırılmaları eklendiği için numaralar kayar, ama .txt bundan etkilenmez.
%   Bu yüzden MC dosyalarını tek tek açmaya gerek yoktur; Excel'deki yüzey
%   numarası nominal dosyadaki numara olmalıdır.
%
% BİRİMLER  (her şey fringe'e çevrilir; 1 fringe = lambda/2 = 316,5 nm, lambda = 633 nm)
%   Excel (IRR P-V, IRR RMS, ASTI) : nm   -> / (lambda_nm/2)
%   TEZI ve Zernike katsayıları    : mm   -> / (lambda_nm*1e-6/2)
%   TIRR                           : fringe (olduğu gibi)
%
% HESAPLAR
%   TEZI (Zernike Standard, Noll sırası):
%     W(rho,theta) = sum c_j Z_j      (excludeTerms hariç: piston, tilt, power)
%     IRR P-V  = max(W) - min(W)      birim diskte gridSize örneklemle
%     IRR RMS  = sqrt(sum c_j^2)      (Zernike Standard ortonormal)
%     ASTI     = 2*sqrt(6)*sqrt(c5^2 + c6^2)   astigmatizma teriminin P-V'si
%     Kontrol  : sqrt(sum c_j^2) ile TEZI sütunundaki değer karşılaştırılır.
%   TIRR: sütundaki değer fringe cinsinden düzensizliktir -> IRR P-V = |değer|.
%     TIRR katsayıları dosyaya yazılmadığı için RMS ve ASTI hesaplanmaz.
%   Normalizasyon yarıçapı Zemax'taki yarı çaptır; ölçüm temiz açıklıkta
%   yapılıyorsa ikisinin aynı olduğundan emin olun.
%   Karşılaştırma: iki örneklemli Kolmogorov-Smirnov D (0 = aynı dağılım,
%   1 = hiç örtüşmüyor). Ölçümle en küçük D'yi veren seri en yakın modeldir.
%
% AYARLAR (User settings)
%   lambda_nm, excludeTerms, gridSize, nBins, measSheet
%
% ÇIKTILAR (ilk .txt dosyasının yanında "MC IRR DAĞILIMI" klasörüne)
%   IRR_S<yüzey>.png        : yüzey başına; IRR P-V yoğunluk ve CDF, TEZI varsa
%                             IRR RMS ve ASTI/IRR yoğunluğu (ölçüm gri)
%   IRR_Tum_Yuzeyler.png    : her seri için tüm yüzeyler birleştirilmiş
%   MC_IRR_Sonuclar.xlsx    : "Per File" (MC dosyası x yüzey), "Statistics"
%                             (ortalama, medyan, std, P5, P95, min, max, KS D),
%                             "Measurements", "Settings"
%   MC_IRR_Data.mat         : tüm seriler ve ölçüm
%   Command Window          : yüzey başına IRR P-V özeti ve KS D sıralaması

clc;
clear;
close all;

%% User settings
lambda_nm    = 633;          % test wavelength (nm); 1 fringe = lambda/2 of surface
excludeTerms = 1:4;          % Zernike terms not counted as irregularity (piston, tilt, power)
gridSize     = [101 360];    % radial x angular samples on the unit disk for P-V
nBins        = 25;           % histogram bins
measSheet    = '';           % measurement sheet, '' = choose in a dialog (e.g. 'Sheet1')

%% Select result files
[fNames, fPath] = uigetfile({'*.txt', 'Zemax tolerance results (*.txt)'}, ...
    'Select Monte Carlo result file(s)', 'MultiSelect', 'on');

if isequal(fNames, 0)
    return;
end

fNames = cellstr(fNames);
frMM = lambda_nm * 1e-6 / 2;   % one fringe in mm (TEZI / Zernike coefficients are in mm)

%% Read and evaluate
series = struct('run', {}, 'op', {}, 'surf', {}, 'label', {}, 'mc', {}, 'value', {}, ...
    'pv', {}, 'rms', {}, 'asti', {}, 'astiRatio', {}, 'pvRms', {});

for f = 1:numel(fNames)
    [~, run] = fileparts(fNames{f});
    R = readMcTable(fullfile(fPath, fNames{f}));
    fprintf('%s: %d MC rows, %d irregularity operands\n', run, numel(R.mc), numel(R.ops));

    for k = 1:numel(R.ops)
        op = R.ops(k);
        val = R.data(:, op.col);
        nMC = numel(val);
        rms = NaN(nMC, 1);
        asti = NaN(nMC, 1);

        switch op.type
            case 'TEZI'
                use = ~ismember(op.zNum, excludeTerms);
                if any(use)
                    Cfr = R.data(:, op.zCol(use)) / frMM;
                    [pv, rms, asti] = zernikeIrr(op.zNum(use), Cfr, gridSize);

                    % Check: RSS of all saved coefficients vs the TEZI value
                    rss = sqrt(sum(R.data(:, op.zCol(op.zNum > 1)).^2, 2));
                    relDiff = max(abs(rss - abs(val)) ./ abs(val));
                    fprintf('   TEZI S%-3d terms %d-%d: RSS(coef)/TEZI max diff %.1e%s\n', op.surf, ...
                        min(op.zNum), max(op.zNum), relDiff, ...
                        repmat(' (TEZI value = surface RMS)', 1, relDiff < 1e-3));
                else
                    warning('TEZI on surface %d: no Zernike coefficient columns found, RMS from TEZI value only.', op.surf);
                    pv = NaN(nMC, 1);
                    rms = abs(val) / frMM;
                end
            case 'TIRR'
                pv = abs(val);
            otherwise
                continue;
        end

        series(end+1) = struct('run', run, 'op', op.type, 'surf', op.surf, ...
            'label', sprintf('%s - %s', run, op.type), 'mc', R.mc, 'value', val, ...
            'pv', pv, 'rms', rms, 'asti', asti, 'astiRatio', asti ./ pv, 'pvRms', pv ./ rms); %#ok<SAGROW>
    end
end

if isempty(series)
    error('No TEZI or TIRR data found in the selected files.');
end

%% Measurement (nm)
[mName, mPath] = uigetfile({'*.xlsx;*.xlsm;*.xls', 'Excel files'}, 'Select measurement Excel', fPath);

if isequal(mName, 0)
    return;
end

mFile = fullfile(mPath, mName);
sh = measSheet;

if isempty(sh)
    sheets = cellstr(sheetnames(mFile));
    [sel, ok] = listdlg('ListString', sheets, 'SelectionMode', 'single', ...
        'PromptString', 'Select measurement sheet:', 'ListSize', [300 150]);
    if ~ok
        return;
    end
    sh = sheets{sel};
end

[M, mInfo] = IRR_Olcum_Oku(mFile, sh, 'nm', lambda_nm);
meas = struct('label', sprintf('Measured - %s', sh), 'file', mFile, 'sheet', sh, ...
    'table', M, 'pv', M.IRR_PV_fr, 'rms', M.IRR_RMS_fr, 'asti', M.ASTI_fr, ...
    'astiRatio', M.ASTI_fr ./ M.IRR_PV_fr, 'pvRms', M.IRR_PV_fr ./ M.IRR_RMS_fr, ...
    'surf', M.Surface);
fprintf('Measurement: %s, sheet %s, n = %d (%s)\n', mName, sh, height(M), mInfo.columns);

%% Keep only the surfaces given in the Excel
if any(isfinite(M.Surface))
    us = unique(M.Surface(isfinite(M.Surface)))';
    fprintf('Measured surfaces: %s\n', strjoin(arrayfun(@(u) sprintf('S%d (n = %d)', u, ...
        sum(M.Surface == u)), us, 'UniformOutput', false), ', '));

    skipped = setdiff(unique([series.surf]), us);
    if ~isempty(skipped)
        fprintf('Not in the Excel, skipped: %s\n', strjoin(arrayfun(@(u) sprintf('S%d', u), skipped, ...
            'UniformOutput', false), ', '));
    end
    missing = setdiff(us, [series.surf]);
    if ~isempty(missing)
        warning('Measured surface(s) %s have no TEZI/TIRR column in the MC results.', mat2str(missing));
    end

    series = series(ismember([series.surf], us));
    if isempty(series)
        error('None of the measured surfaces has a TEZI/TIRR column in the MC results.');
    end
else
    warning('No "SURFACE" column in the measurement sheet: every surface is compared with all measurements.');
end

%% Groups: each surface, then all surfaces pooled per series
surfList = unique([series.surf]);
groups = struct('name', {}, 'file', {}, 'items', {}, 'surfs', {});

for s = surfList
    groups(end+1) = struct('name', sprintf('Surface %d', s), 'file', sprintf('IRR_S%d.png', s), ...
        'items', series([series.surf] == s), 'surfs', s); %#ok<SAGROW>
end

[keys, ~, kIdx] = unique({series.label}, 'stable');
pooled = series([]);
for k = 1:numel(keys)
    sk = series(kIdx == k);
    p = sk(1);
    p.surf = NaN;
    for fld = {'mc', 'value', 'pv', 'rms', 'asti', 'astiRatio', 'pvRms'}
        p.(fld{1}) = vertcat(sk.(fld{1}));
    end
    p.label = sprintf('%s (%d surfaces)', p.label, numel(sk));
    pooled(end+1) = p; %#ok<SAGROW>
end
groups(end+1) = struct('name', 'All surfaces', 'file', 'IRR_Tum_Yuzeyler.png', 'items', pooled, 'surfs', surfList);

%% Statistics, command window summary, plots
outFolder = fullfile(fPath, ['MC IRR DA' char(286) 'ILIMI']);   % "MC IRR DAĞILIMI"
if ~isfolder(outFolder)
    mkdir(outFolder);
end

metrics = {'pv', 'IRR_PV_fr'; 'rms', 'IRR_RMS_fr'; 'asti', 'ASTI_fr'; 'astiRatio', 'ASTI_over_IRR'; 'pvRms', 'PV_over_RMS'};
statRows = {};

for g = 1:numel(groups)
    items = groups(g).items;
    gm = measFor(meas, groups(g).surfs);   % measurements of this group's surface(s)

    fprintf('\n===== %s : IRR P-V [fr] =====\n', groups(g).name);
    if ~isempty(gm)
        q = pct(gm.pv, [0.5 0.95]);
        fprintf('%-40s n %4d  mean %.3f  median %.3f  P95 %.3f  max %.3f\n', gm.label, ...
            sum(isfinite(gm.pv)), mean(gm.pv, 'omitnan'), q(1), q(2), max(gm.pv));
    end

    D = NaN(numel(items), 1);
    for i = 1:numel(items)
        v = items(i).pv;
        q = pct(v, [0.5 0.95]);
        if ~isempty(gm)
            D(i) = ks2(v, gm.pv);
        end
        fprintf('%-40s n %4d  mean %.3f  median %.3f  P95 %.3f  max %.3f  KS D %.3f\n', items(i).label, ...
            sum(isfinite(v)), mean(v, 'omitnan'), q(1), q(2), max(v), D(i));
    end
    if ~isempty(gm) && numel(items) > 1 && any(isfinite(D))
        [~, iBest] = min(D);
        fprintf('Closest to measurement: %s\n', items(iBest).label);
    end

    % Statistics rows
    for m = 1:size(metrics, 1)
        fld = metrics{m, 1};
        if ~isempty(gm) && any(isfinite(gm.(fld)))
            statRows(end+1, :) = statRow(groups(g).name, gm.label, '', 'Measured', metrics{m, 2}, gm.(fld), []); %#ok<SAGROW>
        end
        for i = 1:numel(items)
            v = items(i).(fld);
            if ~any(isfinite(v))
                continue;
            end
            ref = [];
            if ~isempty(gm)
                ref = gm.(fld);
            end
            statRows(end+1, :) = statRow(groups(g).name, items(i).label, items(i).run, items(i).op, ...
                metrics{m, 2}, v, ref); %#ok<SAGROW>
        end
    end

    plotGroup(groups(g).name, items, gm, nBins, fullfile(outFolder, groups(g).file));
end

%% Save results
outputFile = fullfile(outFolder, 'MC_IRR_Sonuclar.xlsx');
if isfile(outputFile)
    delete(outputFile);
end

perFile = {};
for i = 1:numel(series)
    s = series(i);
    n = numel(s.pv);
    perFile = [perFile; cellstr(repmat(string(s.run), n, 1)), cellstr(s.mc), num2cell(repmat(s.surf, n, 1)), ...
        cellstr(repmat(string(s.op), n, 1)), num2cell(s.value), num2cell(s.pv), num2cell(s.rms), ...
        num2cell(s.asti), num2cell(s.astiRatio), num2cell(s.pvRms)]; %#ok<AGROW>
end
PerFile = cell2table(perFile, 'VariableNames', {'Run', 'MC', 'Surface', 'Operand', 'OperandValue', ...
    'IRR_PV_fr', 'IRR_RMS_fr', 'ASTI_fr', 'ASTI_over_IRR', 'PV_over_RMS'});
writetable(PerFile, outputFile, 'Sheet', 'Per File');

Stats = cell2table(statRows, 'VariableNames', {'Group', 'Series', 'Run', 'Operand', 'Metric', 'N', ...
    'Mean', 'Median', 'Std', 'P5', 'P95', 'Min', 'Max', 'KS_D_vs_Measured', 'KS_p'});
writetable(Stats, outputFile, 'Sheet', 'Statistics');

if ~isempty(meas)
    writetable(meas.table, outputFile, 'Sheet', 'Measurements');
end

measText = 'none';
if ~isempty(meas)
    measText = sprintf('%s | %s', meas.file, meas.sheet);
end
Settings = table( ...
    {'Result files'; 'Test wavelength (nm)'; 'Fringe (nm)'; 'Units'; 'Excluded Zernike terms'; ...
     'P-V grid (radial x angular)'; 'Measurement'}, ...
    {strjoin(fNames, ', '); num2str(lambda_nm); num2str(lambda_nm / 2); 'Excel nm, TEZI mm, TIRR fringe -> fringe'; ...
     mat2str(excludeTerms); mat2str(gridSize); measText}, ...
    'VariableNames', {'Setting', 'Value'});
writetable(Settings, outputFile, 'Sheet', 'Settings');

save(fullfile(outFolder, 'MC_IRR_Data.mat'), 'series', 'meas', 'lambda_nm', 'excludeTerms', 'fNames');
fprintf('\nResults saved to: %s\n', outFolder);

%% Local functions
function R = readMcTable(file)
    txt = readTextAuto(file);
    lines = regexp(txt, '\r\n|\n|\r', 'split');
    cells = cell(size(lines));
    for i = 1:numel(lines)
        cells{i} = strtrim(strsplit(lines{i}, char(9), 'CollapseDelimiters', false));
    end

    iName = find(cellfun(@(c) any(strcmpi(c, 'TEZI') | strcmpi(c, 'TIRR')), cells), 1);
    if isempty(iName)
        error('No TEZI / TIRR column found in %s.', file);
    end
    names = cells{iName};
    nCols = numel(names);

    % Data rows: label in the first cell, mostly numbers after it
    isData = false(size(cells));
    for i = iName+1:numel(cells)
        c = cells{i};
        if numel(c) < 2 || isempty(c{1})
            continue;
        end
        v = toNum(c(2:end));
        isData(i) = sum(isfinite(v)) >= 0.5 * (nCols - 1);
    end
    iData = find(isData);
    if isempty(iData)
        error('No Monte Carlo data rows found in %s.', file);
    end

    % Header row without the label column -> shift by one
    nData = median(cellfun(@numel, cells(iData)));
    shift = 0;
    if nData == nCols + 1 && ~isempty(names{1})
        shift = 1;
        names = [{''}, names];
        nCols = nCols + 1;
    end

    % Column info rows ("Surf: 4", "Min: -2", ...)
    attr = repmat({''}, 1, nCols);
    for r = iName+1:iData(1)-1
        c = [repmat({''}, 1, shift), cells{r}];
        for j = 1:min(numel(c), nCols)
            if ~isempty(c{j})
                attr{j} = [attr{j} ' ' c{j}];
            end
        end
    end

    nMC = numel(iData);
    data = NaN(nMC, nCols);
    mc = strings(nMC, 1);
    for k = 1:nMC
        c = cells{iData(k)};
        mc(k) = string(c{1});
        v = toNum(c);
        n = min(numel(v), nCols);
        data(k, 1:n) = v(1:n);
    end

    ops = struct('type', {}, 'surf', {}, 'col', {}, 'zNum', {}, 'zCol', {}, 'minV', {}, 'maxV', {});
    cur = 0;
    for j = 1:nCols
        nm = upper(strtrim(names{j}));
        if any(strcmp(nm, {'TEZI', 'TIRR', 'TEXI'}))
            if strcmp(nm, 'TEXI')
                warning('TEXI (Zernike Fringe) columns are not supported and are skipped.');
            end
            ops(end+1) = struct('type', nm, 'surf', attrNum(attr{j}, 'Surf'), 'col', j, ...
                'zNum', [], 'zCol', [], 'minV', attrNum(attr{j}, 'Min'), 'maxV', attrNum(attr{j}, 'Max')); %#ok<AGROW>
            cur = numel(ops);
        elseif cur > 0 && strcmp(ops(cur).type, 'TEZI') && ~isempty(regexp(nm, '^ZERNIKE\s*\d+$', 'once'))
            ops(cur).zNum(end+1) = str2double(regexp(nm, '\d+', 'match', 'once'));
            ops(cur).zCol(end+1) = j;
        else
            cur = 0;
        end
    end
    ops = ops(~strcmp({ops.type}, 'TEXI'));

    R.mc = mc;
    R.data = data;
    R.ops = ops;
    R.names = names;
end

function txt = readTextAuto(file)
    fid = fopen(file, 'r');
    if fid < 0
        error('Could not open %s.', file);
    end
    b = fread(fid, inf, '*uint8')';
    fclose(fid);

    if numel(b) >= 2 && b(1) == 255 && b(2) == 254            % UTF-16 LE BOM
        txt = utf16(b(3:end), true);
    elseif numel(b) >= 2 && b(1) == 254 && b(2) == 255        % UTF-16 BE BOM
        txt = utf16(b(3:end), false);
    elseif numel(b) >= 3 && b(1) == 239 && b(2) == 187 && b(3) == 191   % UTF-8 BOM
        txt = native2unicode(b(4:end), 'UTF-8');
    elseif numel(b) >= 4 && mean(b(2:2:end) == 0) > 0.3       % UTF-16 LE without BOM
        txt = utf16(b, true);
    else
        txt = native2unicode(b, 'UTF-8');
    end
end

function t = utf16(b, little)
    b = double(b(1:2 * floor(numel(b) / 2)));
    if little
        t = char(b(1:2:end) + 256 * b(2:2:end));
    else
        t = char(256 * b(1:2:end) + b(2:2:end));
    end
end

function v = toNum(c)
    v = str2double(strrep(c, ',', '.'));
end

function v = attrNum(a, key)
    tok = regexp(a, ['(?<![A-Za-z#])' key ':\s*([-+0-9.,Ee]+)'], 'tokens', 'once');
    v = NaN;
    if ~isempty(tok)
        v = str2double(strrep(tok{1}, ',', '.'));
    end
end

function [pv, rms, asti] = zernikeIrr(zNum, Cfr, gridSize)
    % Cfr: nMC x nZ coefficients in fringes, Zernike Standard (Noll) terms zNum
    rho = linspace(0, 1, gridSize(1));
    th = (0:gridSize(2)-1) * 2*pi / gridSize(2);
    [RHO, TH] = ndgrid(rho, th);
    B = zeros(numel(RHO), numel(zNum));
    for i = 1:numel(zNum)
        B(:, i) = zernikeNoll(zNum(i), RHO(:), TH(:));
    end

    bad = any(~isfinite(Cfr), 2);
    Cfr(bad, :) = 0;
    W = B * Cfr.';                                   % nPoints x nMC
    pv = (max(W, [], 1) - min(W, [], 1)).';
    rms = sqrt(sum(Cfr(:, zNum > 1).^2, 2));
    asti = NaN(size(pv));
    i5 = find(zNum == 5, 1);
    i6 = find(zNum == 6, 1);
    if ~isempty(i5) && ~isempty(i6)
        asti = 2 * sqrt(6) * hypot(Cfr(:, i5), Cfr(:, i6));
    end
    pv(bad) = NaN;
    rms(bad) = NaN;
    asti(bad) = NaN;
end

function Z = zernikeNoll(j, rho, theta)
    % OpticStudio Zernike Standard = Noll: Z5 = sqrt(6) rho^2 sin2t, Z6 = sqrt(6) rho^2 cos2t, ...
    n = 0;
    j1 = j - 1;
    while j1 > n
        n = n + 1;
        j1 = j1 - n;
    end
    m = (-1)^j * (mod(n, 2) + 2 * floor((j1 + mod(n + 1, 2)) / 2));
    am = abs(m);

    R = zeros(size(rho));
    for k = 0:(n - am) / 2
        R = R + (-1)^k * factorial(n - k) / (factorial(k) * factorial((n + am)/2 - k) * ...
            factorial((n - am)/2 - k)) * rho.^(n - 2*k);
    end

    if m == 0
        Z = sqrt(n + 1) * R;
    elseif m > 0
        Z = sqrt(2 * (n + 1)) * R .* cos(am * theta);
    else
        Z = sqrt(2 * (n + 1)) * R .* sin(am * theta);
    end
end

function q = pct(x, p)
    x = sort(x(isfinite(x)));
    n = numel(x);
    if n == 0
        q = NaN(size(p));
        return;
    end
    pos = 1 + (n - 1) * p;
    lo = floor(pos);
    hi = ceil(pos);
    q = x(lo).' + (pos - lo) .* (x(hi).' - x(lo).');
end

function [D, p] = ks2(a, b)
    a = a(isfinite(a));
    b = b(isfinite(b));
    if isempty(a) || isempty(b)
        D = NaN;
        p = NaN;
        return;
    end
    x = unique([a(:); b(:)]).';
    Fa = mean(a(:) <= x, 1);
    Fb = mean(b(:) <= x, 1);
    D = max(abs(Fa - Fb));
    ne = numel(a) * numel(b) / (numel(a) + numel(b));
    lam = (sqrt(ne) + 0.12 + 0.11 / sqrt(ne)) * D;
    j = (1:100)';
    p = min(max(2 * sum((-1).^(j - 1) .* exp(-2 * j.^2 * lam^2)), 0), 1);
end

function m = measFor(meas, surfs)
    % Measurement rows of the given surface(s); all rows if the sheet has no surface column
    m = meas;
    if isempty(meas) || ~any(isfinite(meas.surf))
        return;
    end
    sel = ismember(meas.surf, surfs);
    if ~any(sel)
        m = [];
        return;
    end
    for fld = {'pv', 'rms', 'asti', 'astiRatio', 'pvRms', 'surf'}
        m.(fld{1}) = meas.(fld{1})(sel);
    end
    if isscalar(surfs)
        m.label = sprintf('%s S%d', meas.label, surfs);
    end
end

function row = statRow(group, label, run, op, metric, v, ref)
    v = v(isfinite(v));
    q = pct(v, [0.05 0.5 0.95]);
    D = NaN;
    p = NaN;
    if ~isempty(ref)
        [D, p] = ks2(v, ref);
    end
    row = {group, label, run, op, metric, numel(v), mean(v), q(2), std(v), q(1), q(3), ...
        min(v), max(v), D, p};
end

function plotGroup(name, items, meas, nBins, imageFile)
    hasRms = any(arrayfun(@(s) any(isfinite(s.rms)), items));
    if hasRms
        fig = figure('Color', 'w', 'Position', [60 60 1250 860]);
        tl = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    else
        fig = figure('Color', 'w', 'Position', [60 60 1250 450]);
        tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    end
    title(tl, sprintf('%s - Monte Carlo irregularity vs measurement', name), 'Interpreter', 'none');

    pdfTile(nexttile(tl), items, meas, 'pv', 'IRR P-V [fringe]', nBins);
    cdfTile(nexttile(tl), items, meas, 'pv', 'IRR P-V [fringe]');
    if hasRms
        pdfTile(nexttile(tl), items, meas, 'rms', 'IRR RMS [fringe]', nBins);
        pdfTile(nexttile(tl), items, meas, 'astiRatio', 'ASTI / IRR P-V', nBins);
    end

    try
        exportgraphics(fig, imageFile, 'Resolution', 200);   % R2020a+
    catch
        saveas(fig, imageFile);
    end
    close(fig);
end

function pdfTile(ax, items, meas, fld, xl, nBins)
    hold(ax, 'on');
    allv = vertcat(items.(fld));
    hasMeas = ~isempty(meas) && any(isfinite(meas.(fld)));
    if hasMeas
        allv = [allv; meas.(fld)];
    end
    allv = allv(isfinite(allv));
    if isempty(allv)
        axis(ax, 'off');
        return;
    end
    edges = linspace(0, max(max(allv) * 1.02, eps), nBins + 1);
    co = lines(max(numel(items), 1));

    if hasMeas
        v = meas.(fld)(isfinite(meas.(fld)));
        histogram(ax, v, edges, 'Normalization', 'pdf', 'FaceColor', [0.75 0.75 0.75], ...
            'EdgeColor', [0.55 0.55 0.55], 'FaceAlpha', 0.6, ...
            'DisplayName', sprintf('%s (n = %d)', meas.label, numel(v)));
    end
    for i = 1:numel(items)
        v = items(i).(fld)(isfinite(items(i).(fld)));
        if isempty(v)
            continue;
        end
        lbl = sprintf('%s (n = %d)', items(i).label, numel(v));
        if hasMeas
            lbl = sprintf('%s, KS D %.2f', lbl, ks2(v, meas.(fld)));
        end
        histogram(ax, v, edges, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', ...
            'EdgeColor', co(i, :), 'LineWidth', 1.6, 'DisplayName', lbl);
    end
    xlabel(ax, xl);
    ylabel(ax, 'Probability density');
    legend(ax, 'Location', 'best', 'Box', 'off', 'Interpreter', 'none', 'FontSize', 8);
    grid(ax, 'on');
    box(ax, 'on');
end

function cdfTile(ax, items, meas, fld, xl)
    hold(ax, 'on');
    co = lines(max(numel(items), 1));
    if ~isempty(meas) && any(isfinite(meas.(fld)))
        v = sort(meas.(fld)(isfinite(meas.(fld))));
        stairs(ax, [0; v], [0; (1:numel(v))' / numel(v)], 'Color', [0.35 0.35 0.35], ...
            'LineWidth', 2.5, 'DisplayName', meas.label);
    end
    for i = 1:numel(items)
        v = sort(items(i).(fld)(isfinite(items(i).(fld))));
        if isempty(v)
            continue;
        end
        stairs(ax, [0; v], [0; (1:numel(v))' / numel(v)], 'Color', co(i, :), 'LineWidth', 1.5, ...
            'DisplayName', items(i).label);
    end
    ylim(ax, [0 1]);
    xlabel(ax, xl);
    ylabel(ax, 'Cumulative probability');
    legend(ax, 'Location', 'southeast', 'Box', 'off', 'Interpreter', 'none', 'FontSize', 8);
    grid(ax, 'on');
    box(ax, 'on');
end
