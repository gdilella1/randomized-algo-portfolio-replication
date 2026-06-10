function [erroriL1, tabellaLatex] = graficaRepliche(sp100_data, datiStruct, dataInizio, dataFine, dataFutura, k, tipoSketch)
    if nargin < 7, tipoSketch = 'GAUSS'; end
    if nargin < 6, k = 20; end
     
    l = k + 7; 
    P = calcolaCorrelazioneStocks(datiStruct, dataInizio, dataFine);
    sketch = scegliSketch(tipoSketch, l, size(P,1)); 
    
    % Inizializzazione matrice indici e vettore errori
    I = zeros(6, size(P,1));
    erroriL1 = zeros(1, 6);
    alg_nomi = {'LUPP', 'LUPP1', 'CPQR', 'CPQR1', 'DEIM', 'LS'};
    
    [I(1,:), ~] = scegliAlgoritmo('LUPP', P,k, sketch);
    [I(2,:), ~] = scegliAlgoritmo('LUPP1', P,k, sketch);
    [I(3,:), ~] = scegliAlgoritmo('CPQR', P,k, sketch);
    [I(4,:), ~] = scegliAlgoritmo('CPQR1', P,k, sketch);
    [I(5,1:k), ~] = scegliAlgoritmo('DEIM', P,k, sketch);
    [I(6,1:k), ~] = scegliAlgoritmo('LS', P,k, sketch);

    % --- PREPARAZIONE RANGE TEMPORALI ---
    tInizio  = dateshift(datetime(dataInizio), 'start', 'day');
    tFine    = dateshift(datetime(dataFine), 'end', 'day');
    tFuturo  = dateshift(datetime(dataFutura), 'end', 'day');
    
    tr_opt   = timerange(tInizio, tFine, 'closed');
    tr_total = timerange(tInizio, tFuturo, 'closed');
    tr_error = timerange(tFine, tFuturo, 'closed');
    
    if isfield(sp100_data, 'date_iso')
        d_idx = datetime(sp100_data.date_iso);
    else
        d_idx = datetime(sp100_data.Date);
    end
    
    prezzi_target = isfield(sp100_data, 'adjclose');
    if prezzi_target, p_vals = sp100_data.adjclose(:); else, p_vals = sp100_data.close(:); end
    
    targetTT_total = timetable(d_idx(:), p_vals, 'VariableNames', {'TargetIndex'});
    targetTT_total = unique(targetTT_total);
    targetTT_total = targetTT_total(tr_total, :); 

    % --- INIZIO GRAFICO ---
    figure('Color', 'w'); hold on;
    plot(targetTT_total.Time, targetTT_total.TargetIndex, 'k', 'LineWidth', 2, 'DisplayName', 'S&P100 Reale');
    
    yl = [min(targetTT_total.TargetIndex)*0.95, max(targetTT_total.TargetIndex)*1.05];
    line([tFine tFine], yl, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Fine In-Sample');

    % --- CICLO PORTAFOGLI ---
    nomiCampi = fieldnames(datiStruct);
    colori = [1 0 0; 0.85 0.33 0.1; 0 0 1; 0 1 1; 1 0 1; 0 0.5 0];
    
    for alg = 1:6
        indiciSelezionati = I(alg, 1:k);
        nomiScelti = nomiCampi(indiciSelezionati);
        numScelti = length(indiciSelezionati);
        
        elencoTT = cell(numScelti, 1);
        for i = 1:numScelti
            stk = datiStruct.(nomiScelti{i});
            ts_s = datetime(strcat(stk.Date, {' '}, stk.Time));
            tempTT = timetable(ts_s(:), stk.Close(:), 'VariableNames', {nomiScelti{i}});
            elencoTT{i} = retime(unique(tempTT), 'daily', 'lastvalue');
        end

        masterTable_total = synchronize(targetTT_total, elencoTT{:}, 'union', 'linear');
        masterTable_total = masterTable_total(tr_total, :);
        masterTable_total = fillmissing(masterTable_total, 'nearest');
        
        masterTable_opt = masterTable_total(tr_opt, :);
        dataMat_opt = table2array(masterTable_opt);
        y_opt = dataMat_opt(:, 1);
        X_opt = dataMat_opt(:, 2:end);
        
        coeff = lsqnonneg(X_opt, y_opt);
        
        % Calcolo Replica
        X_total = table2array(masterTable_total(:, 2:end));
        replica_total = X_total * coeff;
        
        % --- CALCOLO ERRORE NORMA 1 ---
        masterTable_error_period = masterTable_total(tr_error, :);
        target_vals_error = masterTable_error_period.TargetIndex;
        X_error_period = table2array(masterTable_error_period(:, 2:end));
        replica_vals_error = X_error_period * coeff;
        
        erroriL1(alg) = sum(abs(target_vals_error - replica_vals_error));
        
        plot(masterTable_total.Time, replica_total, 'Color', colori(alg,:), 'LineWidth', 1.2, 'DisplayName', alg_nomi{alg});
    end

    title(['Tracking S&P100 (Errore calcolato da ', dataFine, ' a ', dataFutura, ')']);
    xlabel('Data'); ylabel('Prezzo ($)'); 
    grid on; xlim([tInizio tFuturo]);
    legend('Location', 'bestoutside');
    hold off;

    % --- GENERAZIONE TABELLA LATEX ---
    % nl è il carattere di a capo reale, non la stringa "\n"
    nl = newline; 
    
    % Costruzione dell'intestazione
    % Nota: per LaTeX \hline serve un solo backslash. In sprintf ne scriviamo due.
    % Nota: per il fine riga LaTeX serve \\. In sprintf ne scriviamo quattro.
    header = [sprintf('\\begin{tabular}{|l|r|}'), nl, ...
              sprintf('\\hline'), nl, ...
              sprintf('Algoritmo & Errore $L_1$ \\\\'), nl, ...
              sprintf('\\hline'), nl];
    
    corpo = '';
    for j = 1:6
        % Costruiamo ogni riga: Nome & Valore \\
        riga = sprintf('%s & %.4f \\\\', alg_nomi{j}, erroriL1(j));
        corpo = [corpo, riga, nl];
    end
    
    footer = [sprintf('\\hline'), nl, sprintf('\\end{tabular}')];
    
    tabellaLatex = [header, corpo, footer];
end