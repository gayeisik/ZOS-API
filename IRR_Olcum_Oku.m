function [M, info] = IRR_Olcum_Oku(file, sheet, unit, lambda_nm)
%IRR_OLCUM_OKU  Ölçüm Excel'inden IRR P-V, IRR RMS ve ASTI değerlerini fringe olarak okur.
%   [M, info] = IRR_Olcum_Oku(file, sheet)
%   [M, info] = IRR_Olcum_Oku(file, sheet, unit, lambda_nm)
%
%   unit       'nm' (varsayılan): "[birimler nm]" bloğu okunur ve fringe'e çevrilir
%              'fr'             : "[birimler fr]" bloğu doğrudan okunur
%   lambda_nm  test dalga boyu (varsayılan 632.8), 1 fringe = lambda/2
%
%   M     tablo: LensNo, Surface, IRR_PV_fr, IRR_RMS_fr, ASTI_fr
%         IRR P-V > 0 olan her ölçüm bir satır; sayfada olmayan sütunlar NaN.
%   info  bulunan sütunlar, başlık satırı ve birim
%
%   Sütunlar başlıklarından bulunur ("IRR P-V", "IRR RMS", "ASTI"). Aynı başlık
%   hem nm hem fr bloğunda varsa, üstteki blok başlığına ("... [birimler nm]")
%   bakılarak istenen blok seçilir; blok başlığı yoksa birim değerlerden
%   tahmin edilir (medyan > 10 ise nm). "Lens no" ve "SURFACE (ZEMAX)"
%   sütunlarındaki birleştirilmiş hücreler aşağı doğru doldurulur
%   (sütun yoksa LensNo = "", Surface = NaN).

    if nargin < 3 || isempty(unit)
        unit = 'nm';
    end
    if nargin < 4 || isempty(lambda_nm)
        lambda_nm = 632.8;
    end

    C = readcell(file, 'Sheet', sheet);
    raw = cellfun(@cellText, C, 'UniformOutput', false);
    key = cellfun(@normKey, raw, 'UniformOutput', false);

    [r0, cPV, colUnit] = findHeader(raw, key, 'IRR P-V', unit, [], []);
    if isempty(r0)
        error('Header "IRR P-V" not found on sheet %s.', sheet);
    end
    [~, cRMS] = findHeader(raw, key, 'IRR RMS', colUnit, r0, cPV);
    [~, cAST] = findHeader(raw, key, 'ASTI', colUnit, r0, cPV);

    pv = colValues(C, r0, cPV);
    rms = NaN(size(pv));
    ast = NaN(size(pv));
    if ~isempty(cRMS)
        rms = colValues(C, r0, cRMS);
    end
    if ~isempty(cAST)
        ast = colValues(C, r0, cAST);
    end
    keep = isfinite(pv) & pv > 0;

    if isempty(colUnit)
        if median(pv(keep)) > 10
            colUnit = 'nm';
        else
            colUnit = 'fr';
        end
        warning('Sheet %s: unit of column %s not found in the block title, assumed "%s".', ...
            sheet, colLetter(cPV), colUnit);
    end

    scale = 1;
    if strcmp(colUnit, 'nm')
        scale = 1 / (lambda_nm / 2);
    end

    % Lens number (merged cells -> forward fill)
    lensNo = strings(numel(pv), 1);
    lc = find(strcmp(key(r0, :), 'LENSNO'), 1);
    if ~isempty(lc)
        last = "";
        for i = 1:numel(pv)
            t = C{r0 + i, lc};
            if (isnumeric(t) && isscalar(t) && isfinite(t)) || ischar(t) || (isstring(t) && ~ismissing(t))
                last = string(t);
            end
            lensNo(i) = last;
        end
    end

    % Zemax surface number ("SURFACE (ZEMAX)", merged cells -> forward fill)
    surfNo = NaN(numel(pv), 1);
    sc = find(startsWith(key(r0, :), 'SURF'), 1);
    if ~isempty(sc)
        last = NaN;
        for i = 1:numel(pv)
            v = cellNumber(C{r0 + i, sc});
            if isfinite(v)
                last = v;
            end
            surfNo(i) = last;
        end
    end

    M = table(lensNo(keep), surfNo(keep), pv(keep) * scale, rms(keep) * scale, ast(keep) * scale, ...
        'VariableNames', {'LensNo', 'Surface', 'IRR_PV_fr', 'IRR_RMS_fr', 'ASTI_fr'});

    info.headerRow = r0;
    info.unit = colUnit;
    info.columns = sprintf('IRR P-V: %s, IRR RMS: %s, ASTI: %s (row %d)', colLetter(cPV), ...
        colLetter(cRMS), colLetter(cAST), r0);
end

function [r, c, u] = findHeader(raw, key, header, unit, row, nearCol)
    r = [];
    c = [];
    u = '';
    [rr, cc] = find(strcmp(key, normKey(header)));
    if ~isempty(row)
        sel = rr == row;
        rr = rr(sel);
        cc = cc(sel);
    end
    if isempty(rr)
        return;
    end

    units = arrayfun(@(k) blockUnit(raw, rr(k), cc(k)), 1:numel(rr), 'UniformOutput', false);
    pick = strcmp(units, unit);
    if ~any(pick)
        pick = cellfun(@isempty, units);
    end
    if ~any(pick)
        return;
    end
    rr = rr(pick);
    cc = cc(pick);
    units = units(pick);

    if isempty(nearCol)
        i = 1;
    else
        [~, i] = min(abs(cc - nearCol));
    end
    r = rr(i);
    c = cc(i);
    u = units{i};
end

function v = colValues(C, r0, c)
    v = cellfun(@cellNumber, C(r0+1:end, c));
end

function t = cellText(v)
    if ischar(v) || (isstring(v) && isscalar(v) && ~ismissing(v))
        t = lower(char(v));
    else
        t = '';
    end
end

function k = normKey(t)
    k = upper(t);
    k = strrep(k, char(304), 'I');            % 'İ'
    k = regexprep(k, '[\s\-_\.]', '');
end

function v = cellNumber(c)
    if isnumeric(c) && isscalar(c)
        v = double(c);
    elseif ischar(c) || (isstring(c) && isscalar(c) && ~ismissing(c))
        v = str2double(strrep(char(c), ',', '.'));
    else
        v = NaN;
    end
end

function u = blockUnit(raw, r, c)
    % First text found to the left on the rows above, e.g. "silikon [birimler nm]"
    u = '';
    for rr = r-1:-1:1
        for cc = c:-1:1
            t = raw{rr, cc};
            if ~isempty(t)
                tok = regexp(t, '\<(nm|fr)\>', 'tokens', 'once');
                if ~isempty(tok)
                    u = tok{1};
                    return;
                end
                break;   % text without a unit: try the row above
            end
        end
    end
end

function s = colLetter(c)
    s = '';
    if isempty(c)
        s = '-';
        return;
    end
    while c > 0
        m = mod(c - 1, 26);
        s = [char(65 + m) s]; %#ok<AGROW>
        c = floor((c - 1) / 26);
    end
end
