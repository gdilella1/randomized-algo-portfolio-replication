function plot_sp100_range(S, tStart, tEnd)
%PLOT_SP100_RANGE  Plotta la serie S&P100 tra due date (giornaliero, forward-fill).
%   plot_sp100_range(S, tStart, tEnd)
%   - S: struct con campi date (p.es. S.dates / S.Date / S.time) e valori
%        (p.es. S.marketcapUSD / S.marketcap / S.close / S.value / S.price / S.index)
%   - tStart, tEnd: stringhe 'yyyy-MM-dd' o datetime
%
% Preferenze:
% - Allineamento daily con forward-fill (valore del giorno precedente)
% - Gestione time zone coerente (rimozione TZ)
% - Clip all’intervallo dati, warning se fuori range
% - Asse X con sigle dei mesi (xtickformat('MMM'))
% - Scala Y “intelligente” (USD / miliardi / bilioni) se l’ordine di grandezza lo richiede

    % ---------- Estrai date ----------
    t = pickDateVector(S);
    if isempty(t)
        error('Non trovo un campo date nella struct (attesi: dates, Date, time, Time, ...).');
    end
    % datetime robusto
    if iscellstr(t) || isstring(t) || ischar(t)
        t = datetime(t,'InputFormat','yyyy-MM-dd');
    elseif isnumeric(t)
        t = datetime(t,'ConvertFrom','datenum');
    elseif ~isdatetime(t)
        error('Il campo date deve essere string/cellstr/char/numeric o datetime.');
    end
    % normalizza a inizio giorno e rimuovi TZ
    t = dateshift(t(:), 'start','day');
    try, t.TimeZone = ''; catch, end

    % ---------- Estrai valori ----------
    y = pickValueVector(S);
    if isempty(y)
        error(['Non trovo un campo valori nella struct (attesi: marketcapUSD, marketcap, ', ...
               'close, Close, value, Value, price, index, ...).']);
    end
    y = y(:);
    if numel(y) ~= numel(t)
        error('Lunghezze non coerenti: vettore date e vettore valori.');
    end

    % ---------- Ordina per data ----------
    [t, idx] = sort(t);
    y = y(idx);

    % ---------- Timetable e retime giornaliero (forward-fill) ----------
    TT  = timetable(t, y, 'VariableNames', {'Series'});
    TTd = retime(TT, 'daily', 'previous');          % usa valore precedente
    TTd.Series = fillmissing(TTd.Series, 'next');   % fallback per eventuale buco iniziale

    % ---------- Input date ----------
    if ischar(tStart) || isstring(tStart), tStart = datetime(tStart,'InputFormat','yyyy-MM-dd'); end
    if ischar(tEnd)   || isstring(tEnd),   tEnd   = datetime(tEnd,  'InputFormat','yyyy-MM-dd'); end
    tStart = dateshift(tStart,'start','day');
    tEnd   = dateshift(tEnd,  'start','day');
    try, tStart.TimeZone = ''; tEnd.TimeZone = ''; catch, end
    if tEnd < tStart, error('tEnd deve essere >= tStart.'); end

    % ---------- Clip all’intervallo disponibile ----------
    tMin = TTd.t(1); tMax = TTd.t(end);
    if tStart < tMin
        warning('tStart (%s) antecede i dati (%s). Si parte da %s.', string(tStart), string(tMin), string(tMin));
        tStart = tMin;
    end
    if tEnd > tMax
        warning('tEnd (%s) oltre i dati (%s). Si usa %s.', string(tEnd), string(tMax), string(tMax));
        tEnd = tMax;
    end

    mask = (TTd.t >= tStart) & (TTd.t <= tEnd);
    if ~any(mask), error('Nessun dato nell’intervallo richiesto.'); end
    tplot = TTd.t(mask);
    yplot = TTd.Series(mask);

    % ---------- Scala Y intelligente ----------
    maxv = max(yplot);
    if maxv >= 1e12
        factor = 1e12; unit = 'bilioni USD';
    elseif maxv >= 1e9
        factor = 1e9;  unit = 'miliardi USD';
    else
        factor = 1;    unit = 'unità';
    end

    % ---------- Plot ----------
    figure; %#ok<FIGURE>
    plot(tplot, yplot./factor, 'LineWidth', 1.8);
    grid on;
    xlabel('Data');
    ylabel(['Valore (' unit ')']);
    title(sprintf('S&P 100 — dal %s al %s', datestr(tStart,'dd-mm-yyyy'), datestr(tEnd,'dd-mm-yyyy')));
    xtickformat('MMM');   % sigle dei mesi (in italiano se locale impostata)
end

% ======= Helper: trova campo date più probabile =======
function t = pickDateVector(S)
    cand = {'dates','date','date_iso','Date','Dates','time','Time','timestamp','Timestamp'};
    t = [];
    for k = 1:numel(cand)
        if isfield(S, cand{k})
            t = S.(cand{k});
            return;
        end
    end
end

% ======= Helper: trova campo valori più probabile =======
function y = pickValueVector(S)
    cand = {'marketcapUSD','marketcap','MarketCapUSD','MarketCap', ...
            'close','Close','value','Value','price','Price','index','Index'};
    y = [];
    for k = 1:numel(cand)
        if isfield(S, cand{k})
            y = S.(cand{k});
            return;
        end
    end
end
