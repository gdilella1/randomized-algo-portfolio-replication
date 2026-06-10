function [erroriMedi, tabellaLatex] = calcolaErroriMedi(sp100_data, datiStruct, dataInizio, dataFine, dataFutura, k, tipoSketch, numProve)
    % Se non specificato, esegue 10 prove di default
    if nargin < 8, numProve = 25; end
    if nargin < 7, tipoSketch = 'GAUSS'; end
    if nargin < 6, k = 20; end

    fprintf('Inizio calcolo medi su %d prove...\n', numProve);

    % --- 1. PREPARAZIONE DATI COMUNI (Eseguita una sola volta) ---
    P = calcolaCorrelazioneStocks(datiStruct, dataInizio, dataFine);
    alg_nomi = {'LUPP', 'LUPP1', 'CPQR', 'CPQR1', 'DEIM', 'LS'};
    accumulatoreErrori = zeros(1, 6);
    
    % Range temporali
    tInizio  = dateshift(datetime(dataInizio), 'start', 'day');
    tFine    = dateshift(datetime(dataFine), 'end', 'day');
    tFuturo  = dateshift(datetime(dataFutura), 'end', 'day');
    tr_opt   = timerange(tInizio, tFine, 'closed');
    tr_total = timerange(tInizio, tFuturo, 'closed');
    tr_error = timerange(tFine, tFuturo, 'closed');

    % Preparazione Target (S&P100)
    d_idx = datetime(sp100_data.date_iso); % Assumiamo date_iso come standard
    p_vals = if_then_else(isfield(sp100_data, 'adjclose'), sp100_data.adjclose(:), sp100_data.close(:));
    
    targetTT_total = unique(timetable(d_idx(:), p_vals, 'VariableNames', {'TargetIndex'}));
    targetTT_total = targetTT_total(tr_total, :);
    nomiCampi = fieldnames(datiStruct);

    % --- 2. CICLO DI MONTE CARLO (Prove non deterministiche) ---
    for p = 1:numProve
        fprintf('  Prova %d/%d...\n', p, numProve);
        l = k + 7; 
        sketch = scegliSketch(tipoSketch, l, size(P,1));
        
        % Matrice temporanea per gli indici di questa prova
        I_current = zeros(6, size(P,1));
        [I_current(1,:), ~] = scegliAlgoritmo('LUPP', P, k, sketch);
        [I_current(2,:), ~] = scegliAlgoritmo('LUPP1', P, k, sketch);
        [I_current(3,:), ~] = scegliAlgoritmo('CPQR', P, k, sketch);
        [I_current(4,:), ~] = scegliAlgoritmo('CPQR1', P, k, sketch);
        [I_current(5,1:k), ~] = scegliAlgoritmo('DEIM', P, k, sketch);
        [I_current(6,1:k), ~] = scegliAlgoritmo('LS', P, k, sketch);

        % Calcolo errore per ogni algoritmo in questa prova
        for alg = 1:6
            indiciSelezionati = I_current(alg, 1:k);
            nomiScelti = nomiCampi(indiciSelezionati);
            
            % Sincronizzazione rapida per i titoli scelti
            elencoTT = cell(length(nomiScelti), 1);
            for i = 1:length(nomiScelti)
                stk = datiStruct.(nomiScelti{i});
                ts_s = datetime(strcat(stk.Date, {' '}, stk.Time));
                tempTT = timetable(ts_s(:), stk.Close(:), 'VariableNames', {nomiScelti{i}});
                elencoTT{i} = retime(unique(tempTT), 'daily', 'lastvalue');
            end
            
            masterTable = synchronize(targetTT_total, elencoTT{:}, 'union', 'linear');
            masterTable = fillmissing(masterTable(tr_total, :), 'nearest');
            
            % Ottimizzazione pesi (In-Sample)
            mOpt = masterTable(tr_opt, :);
            coeff = lsqnonneg(table2array(mOpt(:, 2:end)), table2array(mOpt(:, 1)));
            
            % Calcolo Errore (Out-of-Sample)
            mError = masterTable(tr_error, :);
            target_vals = mError.TargetIndex;
            replica_vals = table2array(mError(:, 2:end)) * coeff;
            
            accumulatoreErrori(alg) = accumulatoreErrori(alg) + sum(abs(target_vals - replica_vals));
        end
    end

    % --- 3. CALCOLO MEDIE E GENERAZIONE OUTPUT ---
    erroriMedi = accumulatoreErrori / numProve;

    % Generazione Latex
    nl = newline; 
    header = [sprintf('\\begin{tabular}{|l|r|}'), nl, ...
              sprintf('\\hline'), nl, ...
              sprintf('Algoritmo & Errore Medio $L_1$ (%d prove) \\\\', numProve), nl, ...
              sprintf('\\hline'), nl];
    
    corpo = '';
    for j = 1:6
        corpo = [corpo, sprintf('%s & %.4f \\\\', alg_nomi{j}, erroriMedi(j)), nl];
    end
    
    tabellaLatex = [header, corpo, sprintf('\\hline'), nl, sprintf('\\end{tabular}')];
    fprintf('Calcolo completato.\n');
end

% Funzione di utilità interna per pulizia codice
% (Se non disponibile, puoi sostituirla con un IF classico)
function out = if_then_else(cond, valTrue, valFalse)
    if cond, out = valTrue; else, out = valFalse; end
end