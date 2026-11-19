%% IRR ÖLÇÜM DAĞILIMI - FORMÜL ÇIKARMA (MATLAB)
% NE YAPAR
%   Ölçüm Excel'indeki IRR P-V değerlerine (her sayfa = bir malzeme, ör. Silikon,
%   Germanyum, ZnS) olasılık dağılımları uydurur, en iyi uyan dağılımın
%   formülünü ve parametrelerini verir. Bu parametrelerle MATLAB'da yeni IRR
%   değerleri üretmek için IRR_Uret.m fonksiyonu kullanılır (ör. Monte Carlo
%   dosyalarına işlemek üzere).
%
% NASIL ÇALIŞIR
%   1. Excel dosyası ve sayfa(lar) seçilir.
%   2. Her sayfada irrHeader ("IRR P-V", "IRR RMS" veya "ASTI") başlıklı sütun
%      IRR_Olcum_Oku.m ile bulunur. Sayfada aynı
%      başlık iki blokta varsa ("[birimler nm]" ve "[birimler fr]"), sourceUnit
%      ile seçilen bloktaki sütun kullanılır. Başlığın altındaki tüm sayısal
%      hücreler veri kabul edilir (sıfır ve boş hücreler atlanır). "Lens no"
%      sütunu varsa birleştirilmiş hücreler aşağı doğru doldurularak her ölçüme
%      lens numarası atanır.
%   3. nm verisi fringe'e çevrilir:  x[fr] = IRR[nm] / (lambda_nm / 2)
%   4. Lognormal, Gamma ve Weibull dağılımları en çok olabilirlik (MLE) ile
%      uydurulur (Statistics Toolbox gerekmez). En iyi model AIC ile seçilir.
%   5. Seçilen modelden nCheck değer üretilir ve ölçümle karşılaştırılır
%      (üretecin doğru çalıştığının kontrolü).
%
% HESAPLAR / FORMÜLLER  (x = IRR P-V, fringe)
%   Lognormal : p(x) = exp(-(ln x - mu)^2 / (2 sigma^2)) / (x sigma sqrt(2 pi))
%               mu = mean(ln x), sigma = std(ln x) (1/n)
%               Üretim: x = exp(mu + sigma*randn(N,1))
%   Gamma     : p(x) = x^(k-1) exp(-x/theta) / (Gamma(k) theta^k)
%               Üretim: IRR_Uret (Marsaglia-Tsang) veya gamrnd(k,theta,N,1)
%   Weibull   : p(x) = (k/lambda) (x/lambda)^(k-1) exp(-(x/lambda)^k)
%               Üretim: x = lambda*(-log(rand(N,1))).^(1/k)
%   AIC = 2*2 - 2*logL   (küçük olan daha iyi)
%   KS  = Kolmogorov-Smirnov uyum testi. Parametreler aynı veriden
%         hesaplandığı için p değeri iyimserdir, sadece kaba kontrol içindir.
%
% AYARLAR (User settings)
%   lambda_nm, irrHeader, sourceUnit, sheetNames, models, truncMax_fr,
%   nCheck, rngSeed
%
% ÇIKTILAR (Excel dosyasının yanında "IRR DAĞILIM ANALİZİ" klasörüne)
%   IRR_Dagilim_<sayfa>.png      : histogram + model yoğunlukları, ölçülen ve
%                                  model CDF'leri, üretilen değerlerin CDF'i
%   IRR_Dagilim_Sonuclar.xlsx    : "Parameters" (tüm modeller, AIC, KS, formül,
%                                  MATLAB üretim kodu), "Measurements", "Settings"
%   IRR_Dagilim_Modeller.mat     : irrModel (sayfa başına en iyi model) ve
%                                  allFits; IRR_Uret ile kullanılır
%   Command Window               : her sayfa için formül ve parametreler
%
% ÜRETİM ÖRNEĞİ
%   S = load('IRR_Dagilim_Modeller.mat');
%   m = S.irrModel(strcmp({S.irrModel.sheet}, 'Silikon'));
%   rng(1);  irr_fr = IRR_Uret(m, 500);          % fringe
%   irr_mm = IRR_Uret(m, 500, 'mm');             % yüzey P-V, mm

clc;
clear;
close all;

%% User settings
lambda_nm   = 632.8;      % interferometer test wavelength (nm); 1 fringe = lambda/2 of surface
irrHeader   = 'IRR P-V';  % column to fit: 'IRR P-V', 'IRR RMS' or 'ASTI'
sourceUnit  = 'nm';       % 'nm': use the "[birimler nm]" block and convert; 'fr': use the fringe block
sheetNames  = {};         % {} = choose in a dialog, e.g. {'Silikon'}
models      = {'Lognormal', 'Gamma', 'Weibull'};
truncMax_fr = [];         % optional upper limit for generated values (fr), [] = no limit
nCheck      = 1e5;        % values generated from the best model to check the generator
rngSeed     = 1;          % random seed for the check

%% Select Excel file and sheets
[fName, fPath] = uigetfile({'*.xlsx;*.xlsm;*.xls', 'Excel files'}, 'Select measurement Excel file');

if isequal(fName, 0)
    return;
end

excelFile = fullfile(fPath, fName);
allSheets = cellstr(sheetnames(excelFile));   % R2019b+

if isempty(sheetNames)
    [sel, ok] = listdlg('ListString', allSheets, 'SelectionMode', 'multiple', ...
        'PromptString', 'Select sheets (materials):', 'ListSize', [300 150]);
    if ~ok
        return;
    end
    sheetNames = allSheets(sel);
end

outFolder = fullfile(fPath, ['IRR DA' char(286) 'ILIM ANAL' char(304) 'Z' char(304)]);   % "IRR DAĞILIM ANALİZİ"

if ~isfolder(outFolder)
    mkdir(outFolder);
end

%% Fit each sheet
paramRows = {};
dataRows = {};
irrModel = struct('sheet', {}, 'header', {}, 'type', {}, 'p', {}, 'pNames', {}, ...
    'lambda_nm', {}, 'truncMax_fr', {}, 'n', {}, 'formula', {}, 'matlabCode', {});
allFits = struct('sheet', {}, 'fits', {});

for s = 1:numel(sheetNames)
    sheet = sheetNames{s};
    [M, info] = IRR_Olcum_Oku(excelFile, sheet, sourceUnit, lambda_nm);
    switch upper(regexprep(irrHeader, '[\s\-_]', ''))
        case 'IRRPV'
            x = M.IRR_PV_fr;
        case 'IRRRMS'
            x = M.IRR_RMS_fr;
        case 'ASTI'
            x = M.ASTI_fr;
        otherwise
            error('irrHeader must be ''IRR P-V'', ''IRR RMS'' or ''ASTI''.');
    end
    keep = isfinite(x) & x > 0;
    x = x(keep);
    lensNo = M.LensNo(keep);

    if numel(x) < 5
        warning('Sheet %s: only %d values found, skipped.', sheet, numel(x));
        continue;
    end

    fits = fitModels(x, models);
    [~, iBest] = min([fits.aic]);
    best = fits(iBest);

    % Model quantiles (numerical inverse CDF, no toolbox)
    q = [0.5 0.9 0.95 0.99];
    xq = arrayfun(@(qq) modelInv(best.name, best.p, qq, max(x)), q);

    fprintf('\n===== %s  (%s, unit %s, n = %d) =====\n%s\n', sheet, irrHeader, info.unit, numel(x), info.columns);
    fprintf('Measured  : mean %.3f  median %.3f  std %.3f  min %.3f  max %.3f fr\n', ...
        mean(x), median(x), std(x), min(x), max(x));
    for m = 1:numel(fits)
        mark = '';
        if m == iBest
            mark = '   <-- best';
        end
        fprintf('%-10s: %s = %.4f, %s = %.4f   AIC %.2f   KS D %.3f (p %.2f)%s\n', fits(m).name, ...
            fits(m).pNames{1}, fits(m).p(1), fits(m).pNames{2}, fits(m).p(2), fits(m).aic, fits(m).ksD, fits(m).ksP, mark);
    end
    fprintf('Formula   : %s\n', best.formula);
    fprintf('MATLAB    : %s\n', best.code);
    fprintf('Model     : median %.3f  P90 %.3f  P95 %.3f  P99 %.3f fr\n', xq);

    % Generator check
    rng(rngSeed);
    mdl = struct('type', best.name, 'p', best.p, 'lambda_nm', lambda_nm, 'truncMax_fr', truncMax_fr);
    xGen = IRR_Uret(mdl, nCheck, 'fr');
    fprintf('Generated : mean %.3f  median %.3f  std %.3f fr  (N = %d)\n', mean(xGen), median(xGen), std(xGen), nCheck);

    % Store
    irrModel(end+1) = struct('sheet', sheet, 'header', irrHeader, 'type', best.name, 'p', best.p, ...
        'pNames', {best.pNames}, 'lambda_nm', lambda_nm, 'truncMax_fr', truncMax_fr, 'n', numel(x), ...
        'formula', best.formula, 'matlabCode', best.code); %#ok<SAGROW>
    allFits(end+1) = struct('sheet', sheet, 'fits', fits); %#ok<SAGROW>

    for m = 1:numel(fits)
        paramRows(end+1, :) = {sheet, fits(m).name, m == iBest, numel(x), fits(m).pNames{1}, fits(m).p(1), ...
            fits(m).pNames{2}, fits(m).p(2), fits(m).logL, fits(m).aic, fits(m).ksD, fits(m).ksP, ...
            fits(m).formula, fits(m).code}; %#ok<SAGROW>
    end
    for i = 1:numel(x)
        dataRows(end+1, :) = {sheet, lensNo(i), x(i) * lambda_nm / 2, x(i)}; %#ok<SAGROW>
    end

    plotSheet(x, fits, iBest, xGen, sheet, irrHeader, lambda_nm, ...
        fullfile(outFolder, ['IRR_Dagilim_' matlab.lang.makeValidName(sheet) '.png']));
end

if isempty(irrModel)
    error('No sheet with enough data.');
end

%% Save results
outputFile = fullfile(outFolder, 'IRR_Dagilim_Sonuclar.xlsx');

if isfile(outputFile)
    delete(outputFile);
end

Params = cell2table(paramRows, 'VariableNames', {'Sheet', 'Model', 'Best', 'N', 'Param1', 'Value1', ...
    'Param2', 'Value2', 'LogL', 'AIC', 'KS_D', 'KS_p', 'Formula', 'MatlabCode'});
writetable(Params, outputFile, 'Sheet', 'Parameters');

Data = cell2table(dataRows, 'VariableNames', {'Sheet', 'LensNo', 'IRR_PV_nm', 'IRR_PV_fr'});
writetable(Data, outputFile, 'Sheet', 'Measurements');

Settings = table( ...
    {'Excel file'; 'Header'; 'Source unit'; 'Test wavelength (nm)'; 'Fringe (nm)'; 'Truncation max (fr)'; 'Check N'; 'Seed'}, ...
    {excelFile; irrHeader; sourceUnit; num2str(lambda_nm); num2str(lambda_nm / 2); num2str(truncMax_fr); ...
    num2str(nCheck); num2str(rngSeed)}, ...
    'VariableNames', {'Setting', 'Value'});
writetable(Settings, outputFile, 'Sheet', 'Settings');

save(fullfile(outFolder, 'IRR_Dagilim_Modeller.mat'), 'irrModel', 'allFits', 'excelFile', 'lambda_nm');

fprintf('\nResults saved to: %s\n', outFolder);

%% Local functions
function fits = fitModels(x, models)
    x = x(:);
    n = numel(x);
    lx = log(x);
    fits = struct('name', {}, 'p', {}, 'pNames', {}, 'logL', {}, 'aic', {}, 'ksD', {}, 'ksP', {}, ...
        'formula', {}, 'code', {});

    for m = 1:numel(models)
        name = models{m};
        switch lower(name)
            case 'lognormal'
                name = 'Lognormal';
                mu = mean(lx);
                sg = sqrt(mean((lx - mu).^2));
                p = [mu sg];
                pNames = {'mu', 'sigma'};
                formula = sprintf(['p(x) = exp(-(ln x - mu)^2 / (2 sigma^2)) / (x sigma sqrt(2 pi)),  ' ...
                    'mu = %.4f, sigma = %.4f  (median = exp(mu) = %.4f fr)'], mu, sg, exp(mu));
                code = sprintf('x = exp(%.4f + %.4f*randn(N,1));', mu, sg);
            case 'gamma'
                name = 'Gamma';
                sl = log(mean(x)) - mean(lx);
                k = (3 - sl + sqrt((sl - 3)^2 + 24*sl)) / (12*sl);
                for it = 1:100
                    dk = (log(k) - psi(k) - sl) / (1/k - psi(1, k));
                    k = k - dk;
                    if abs(dk) < 1e-12 * k
                        break;
                    end
                end
                p = [k mean(x)/k];
                pNames = {'k', 'theta'};
                formula = sprintf('p(x) = x^(k-1) exp(-x/theta) / (Gamma(k) theta^k),  k = %.4f, theta = %.4f', p);
                code = sprintf('x = IRR_Uret(struct(''type'',''Gamma'',''p'',[%.4f %.4f]), N);  %% or gamrnd(%.4f, %.4f, N, 1)', p, p);
            case 'weibull'
                name = 'Weibull';
                f = @(k) sum(x.^k .* lx) / sum(x.^k) - 1/k - mean(lx);
                k = fzero(f, [0.05 50]);
                p = [k mean(x.^k)^(1/k)];
                pNames = {'k', 'lambda'};
                formula = sprintf('p(x) = (k/lambda) (x/lambda)^(k-1) exp(-(x/lambda)^k),  k = %.4f, lambda = %.4f', p);
                code = sprintf('x = %.4f*(-log(rand(N,1))).^(1/%.4f);', p(2), p(1));
            otherwise
                error('Unknown model: %s', name);
        end

        logL = sum(log(modelPdf(name, p, x)));

        % Kolmogorov-Smirnov
        xs = sort(x);
        F = modelCdf(name, p, xs);
        i = (1:n)';
        D = max(max(i/n - F, F - (i-1)/n));
        lam = (sqrt(n) + 0.12 + 0.11/sqrt(n)) * D;
        j = (1:100)';
        pKS = min(max(2 * sum((-1).^(j-1) .* exp(-2 * j.^2 * lam^2)), 0), 1);

        fits(end+1) = struct('name', name, 'p', p, 'pNames', {pNames}, 'logL', logL, ...
            'aic', 2*2 - 2*logL, 'ksD', D, 'ksP', pKS, 'formula', formula, 'code', code); %#ok<AGROW>
    end
end

function y = modelPdf(name, p, x)
    y = zeros(size(x));
    ok = x > 0;
    xo = x(ok);
    switch lower(name)
        case 'lognormal'
            y(ok) = exp(-(log(xo) - p(1)).^2 / (2*p(2)^2)) ./ (xo * p(2) * sqrt(2*pi));
        case 'gamma'
            y(ok) = exp((p(1) - 1)*log(xo) - xo/p(2) - p(1)*log(p(2)) - gammaln(p(1)));
        case 'weibull'
            y(ok) = (p(1)/p(2)) * (xo/p(2)).^(p(1) - 1) .* exp(-(xo/p(2)).^p(1));
    end
end

function F = modelCdf(name, p, x)
    F = zeros(size(x));
    ok = x > 0;
    xo = x(ok);
    switch lower(name)
        case 'lognormal'
            F(ok) = 0.5 * erfc(-(log(xo) - p(1)) / (p(2) * sqrt(2)));
        case 'gamma'
            F(ok) = gammainc(xo / p(2), p(1));
        case 'weibull'
            F(ok) = 1 - exp(-(xo/p(2)).^p(1));
    end
end

function xq = modelInv(name, p, q, xMax)
    hi = xMax;
    while modelCdf(name, p, hi) < q
        hi = 2 * hi;
    end
    xq = fzero(@(t) modelCdf(name, p, t) - q, [1e-12 hi]);
end

function plotSheet(x, fits, iBest, xGen, sheet, header, lambda_nm, imageFile)
    n = numel(x);
    xMaxPlot = max(1.15 * max(x), modelInv(fits(iBest).name, fits(iBest).p, 0.995, max(x)));
    xg = linspace(xMaxPlot * 1e-3, xMaxPlot, 500)';
    styles = {'-', '--', ':', '-.'};
    colors = [0.12 0.47 0.71; 0.85 0.37 0.01; 0.17 0.63 0.17; 0.58 0.40 0.74];

    fig = figure('Color', 'w', 'Position', [100 100 1150 460]);
    tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('%s - %s, n = %d   (fringe = %s/2, %s = %.1f nm)', sheet, header, n, ...
        char(955), char(955), lambda_nm), 'Interpreter', 'none');

    % PDF
    ax1 = nexttile(tl);
    hold(ax1, 'on');
    histogram(ax1, x, 'Normalization', 'pdf', 'FaceColor', [0.79 0.84 0.91], ...
        'EdgeColor', [0.36 0.48 0.65], 'DisplayName', 'Measured');
    for m = 1:numel(fits)
        lw = 1.3;
        lbl = sprintf('%s (AIC %.1f)', fits(m).name, fits(m).aic);
        if m == iBest
            lw = 2.5;
            lbl = [lbl ' - best']; %#ok<AGROW>
        end
        plot(ax1, xg, modelPdf(fits(m).name, fits(m).p, xg), styles{m}, 'Color', colors(m, :), ...
            'LineWidth', lw, 'DisplayName', lbl);
    end
    yl = ylim(ax1);
    plot(ax1, x, -0.03 * yl(2) * ones(size(x)), '|', 'Color', [0.2 0.2 0.2], 'MarkerSize', 10, ...
        'HandleVisibility', 'off');
    ylim(ax1, [-0.06 * yl(2) yl(2)]);
    xlim(ax1, [0 xMaxPlot]);
    xlabel(ax1, 'IRR P-V [fringe]');
    ylabel(ax1, 'Probability density');
    legend(ax1, 'Location', 'northeast', 'Box', 'off', 'Interpreter', 'none');
    text(ax1, 0.98, 0.45, sprintf('mean %.2f  median %.2f\nstd %.2f  max %.2f fr', ...
        mean(x), median(x), std(x), max(x)), 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'FontSize', 8, 'Color', [0.3 0.3 0.3]);
    grid(ax1, 'on');
    box(ax1, 'on');

    % CDF
    ax2 = nexttile(tl);
    hold(ax2, 'on');
    xs = sort(x);
    stairs(ax2, [0; xs; xMaxPlot], [0; (1:n)'/n; 1], 'k', 'LineWidth', 1.5, 'DisplayName', 'Measured');
    for m = 1:numel(fits)
        lw = 1.3;
        if m == iBest
            lw = 2.5;
        end
        plot(ax2, xg, modelCdf(fits(m).name, fits(m).p, xg), styles{m}, 'Color', colors(m, :), ...
            'LineWidth', lw, 'DisplayName', fits(m).name);
    end
    gs = sort(xGen);
    idx = unique(round(linspace(1, numel(gs), 400)));
    plot(ax2, gs(idx), idx / numel(gs), 'o', 'Color', colors(iBest, :), 'MarkerSize', 3, ...
        'DisplayName', sprintf('Generated (IRR\\_Uret, N = %d)', numel(gs)));
    xlim(ax2, [0 xMaxPlot]);
    ylim(ax2, [0 1]);
    xlabel(ax2, 'IRR P-V [fringe]');
    ylabel(ax2, 'Cumulative probability');
    legend(ax2, 'Location', 'southeast', 'Box', 'off');
    grid(ax2, 'on');
    box(ax2, 'on');

    try
        exportgraphics(fig, imageFile, 'Resolution', 300);   % R2020a+
    catch
        saveas(fig, imageFile);
    end
end
