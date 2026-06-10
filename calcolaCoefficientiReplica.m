function coeff = calcolaCoefficientiReplica(sp100_data, datiStruct, dataInizio, dataFine, dataFutura, algoritmo, k, sketch)
    if nargin < 8
        sketch = randn(30,100)/sqrt(30); %la rigenero ogni volta
    end
    if nargin < 7
        k = 20;
    end

    P = calcolaCorrelazioneStocks(datiStruct, dataInizio, dataFine);
    [I, ~] = scegliAlgoritmo(algoritmo, P,k, sketch);
    indiciSelezionati = I(1:k);

    nomiCampi = fieldnames(datiStruct);
    nomiScelti = nomiCampi(indiciSelezionati);
    numScelti = length(indiciSelezionati);
    
    % Definiamo i tre punti temporali
    tInizio  = dateshift(datetime(dataInizio), 'start', 'day');
    tFine    = dateshift(datetime(dataFine), 'end', 'day');
    tFuturo  = dateshift(datetime(dataFutura), 'end', 'day');
    
    % Range per l'ottimizzazione (In-Sample) e per il grafico totale
    tr_opt   = timerange(tInizio, tFine, 'closed');
    tr_total = timerange(tInizio, tFuturo, 'closed');

    % --- 1. COSTRUZIONE TARGET (S&P100) TOTALE ---
    if isfield(sp100_data, 'date_iso')
        d_idx = datetime(sp100_data.date_iso);
    else
        d_idx = datetime(sp100_data.Date);
    end
    
    prezzi_target = isfield(sp100_data, 'adjclose');
    if prezzi_target, p_vals = sp100_data.adjclose(:); else, p_vals = sp100_data.close(:); end
    
    targetTT_total = timetable(d_idx(:), p_vals, 'VariableNames', {'TargetIndex'});

    % --- 2. COSTRUZIONE TITOLI TOTALE (Intraday -> Daily) ---
    elencoTT = cell(numScelti, 1);
    for i = 1:numScelti
        stk = datiStruct.(nomiScelti{i});
        ts_s = datetime(strcat(stk.Date, {' '}, stk.Time));
        tempTT = timetable(ts_s(:), stk.Close(:), 'VariableNames', {nomiScelti{i}});
        elencoTT{i} = retime(tempTT, 'daily', 'lastvalue');
    end

    % --- 3. SINCRONIZZAZIONE SUL RANGE TOTALE ---
    masterTable_total = synchronize(targetTT_total, elencoTT{:}, 'intersection');
    masterTable_total = masterTable_total(tr_total, :);
    
    % --- 4. OTTIMIZZAZIONE SOLO SUL PRIMO PERIODO (In-Sample) ---
    masterTable_opt = masterTable_total(tr_opt, :);
    masterTable_opt = masterTable_opt(all(~isnan(table2array(masterTable_opt)), 2), :);
    
    dataMat_opt = table2array(masterTable_opt);
    y_opt = dataMat_opt(:, 1);
    X_opt = dataMat_opt(:, 2:end);
    
    coeff = lsqnonneg(X_opt, y_opt);
    
    % --- 5. CALCOLO REPLICA SUL RANGE TOTALE ---
    % Applichiamo i coeff calcolati prima a tutta la matrice dei prezzi
    X_total = table2array(masterTable_total(:, 2:end));
    y_real_total = masterTable_total.TargetIndex;
    replica_total = X_total * coeff;

    % --- 6. GRAFICO E RISULTATI ---
    figure;
    plot(masterTable_total.Time, y_real_total, 'k', 'LineWidth', 2); hold on;
    plot(masterTable_total.Time, replica_total, 'r', 'LineWidth', 1.5);
    
    % Aggiungiamo una linea verticale per separare i periodi
    yl = ylim;
    line([tFine tFine], yl, 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 2);

    legend('S&P100 Reale', 'Portafoglio Replica');
    title(['Tracking S&P100']);
    xlabel('Data'); ylabel('Prezzo'); grid on;

    % Statistiche separate
    replica_opt = replica_total(masterTable_total.Time <= tFine);
    y_opt_check = y_real_total(masterTable_total.Time <= tFine);
    corr_in = corr(replica_opt, y_opt_check, 'Rows', 'complete');
    
    fprintf('Correlazione In-Sample (fino al %s): %.4f\n', dataFine, corr_in);
end