function x = IRR_Uret(model, N, unit)
%IRR_URET  Ölçümden uydurulan dağılıma göre rastgele IRR P-V değerleri üretir.
%   x = IRR_Uret(model, N)         N x 1 değer, fringe
%   x = IRR_Uret(model, N, unit)   unit: 'fr'    (varsayılan) fringe
%                                        'nm'    yüzey P-V, nm      (fr * lambda/2)
%                                        'mm'    yüzey P-V, mm      (Zernike katsayıları için)
%                                        'dalga' yüzey P-V, test dalga boyu cinsinden (fr / 2)
%
%   model: IRR_Dagilim_Analizi.m'nin kaydettiği irrModel'in bir elemanı
%          (IRR_Dagilim_Modeller.mat) veya şu alanları olan bir struct:
%            type        'Lognormal' | 'Gamma' | 'Weibull'
%            p           [p1 p2]  Lognormal: [mu sigma], Gamma: [k theta], Weibull: [k lambda]
%            lambda_nm   test dalga boyu (nm), 'nm'/'mm' için gerekli (varsayılan 632.8)
%            truncMax_fr (isteğe bağlı) bu değerin üstü yeniden çekilir (fr)
%
%   Statistics Toolbox gerekmez. Tekrarlanabilirlik için önce rng(seed) çağırın.
%
%   Örnek:
%     S = load('IRR_Dagilim_Modeller.mat');
%     m = S.irrModel(strcmp({S.irrModel.sheet}, 'Silikon'));
%     rng(1);  irr_fr = IRR_Uret(m, 500);
%
%     m = struct('type', 'Lognormal', 'p', [-0.8085 0.5994]);   % elle
%     irr_fr = IRR_Uret(m, 1000);

    if nargin < 3 || isempty(unit)
        unit = 'fr';
    end

    x = drawValues(model.type, model.p, N);

    if isfield(model, 'truncMax_fr') && ~isempty(model.truncMax_fr)
        bad = x > model.truncMax_fr;
        it = 0;
        while any(bad)
            x(bad) = drawValues(model.type, model.p, nnz(bad));
            bad = x > model.truncMax_fr;
            it = it + 1;
            if it > 10000
                error('truncMax_fr = %g is too small for this distribution.', model.truncMax_fr);
            end
        end
    end

    lambda_nm = 632.8;
    if isfield(model, 'lambda_nm') && ~isempty(model.lambda_nm)
        lambda_nm = model.lambda_nm;
    end

    switch lower(unit)
        case 'fr'
        case 'nm'
            x = x * lambda_nm / 2;
        case 'mm'
            x = x * lambda_nm / 2 * 1e-6;
        case 'dalga'
            x = x / 2;
        otherwise
            error('Unknown unit: %s', unit);
    end
end

function x = drawValues(type, p, N)
    switch lower(type)
        case 'lognormal'
            x = exp(p(1) + p(2) * randn(N, 1));
        case 'gamma'
            x = p(2) * gammaDraw(p(1), N);
        case 'weibull'
            x = p(2) * (-log(rand(N, 1))).^(1 / p(1));
        otherwise
            error('Unknown distribution type: %s', type);
    end
end

function g = gammaDraw(k, N)
    % Marsaglia & Tsang (2000), scale 1
    boost = ones(N, 1);
    if k < 1
        boost = rand(N, 1).^(1 / k);
        k = k + 1;
    end
    d = k - 1/3;
    c = 1 / sqrt(9 * d);
    g = zeros(N, 1);
    todo = true(N, 1);
    while any(todo)
        m = nnz(todo);
        z = randn(m, 1);
        v = (1 + c * z).^3;
        v(v <= 0) = NaN;
        u = rand(m, 1);
        ok = log(u) < 0.5 * z.^2 + d - d * v + d * log(v);
        idx = find(todo);
        g(idx(ok)) = d * v(ok);
        todo(idx(ok)) = false;
    end
    g = g .* boost;
end
