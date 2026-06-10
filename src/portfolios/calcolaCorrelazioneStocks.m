function corrMatrix = calcolaCorrelazioneStocks(datiStruct, dataInizioStr, dataFineStr)
    % calcolaCorrelazioneStocks: Matrice di correlazione 100x100.
    % Risolve l'errore 'Unrecognized table variable name ts'.

    nomiTitoli = fieldnames(datiStruct);
    numTitoli = length(nomiTitoli);
    
    % Definizione dei limiti temporali (giornata intera)
    tInizio = dateshift(datetime(dataInizioStr), 'start', 'day');
    tFine = dateshift(datetime(dataFineStr), 'end', 'day');
    
    elencoTimetable = cell(numTitoli, 1);

    for i = 1:numTitoli
        titolo = nomiTitoli{i};
        dataStock = datiStruct.(titolo);
        
        % Conversione corretta dei timestamp 
        % Usiamo (:) per assicurarci che siano vettori colonna
        ts = datetime(strcat(dataStock.Date, {' '}, dataStock.Time));
        prezzi = dataStock.Close(:);
        
        % Creazione della timetable
        % ts(:) diventa l'indice temporale (proprietà .Time)
        tt = timetable(ts(:), prezzi, 'VariableNames', {titolo});
        
        % CORREZIONE: Per filtrare una timetable si usa la funzione timerange
        % oppure si accede alla proprietà .Time invece di .ts
        tr = timerange(tInizio, tFine, 'closed'); 
        tt = tt(tr, :);
        
        elencoTimetable{i} = tt;
    end

    % Sincronizzazione dei dati su un asse temporale comune
    % 'union' mantiene tutti gli orari, 'linear' interpola i valori mancanti
    masterTable = synchronize(elencoTimetable{:}, 'union', 'linear');
    
    % Estrazione matrice prezzi e calcolo rendimenti logaritmici
    prezziMatrix = table2array(masterTable);
    rendimenti = diff(log(prezziMatrix));
    
    % Pulizia righe con NaN (es. mercati chiusi o dati assenti in certi intervalli)
    rendimenti = rendimenti(all(~isnan(rendimenti), 2), :);

    % Calcolo matrice di correlazione
    corrMatrix = corr(rendimenti);
    
    fprintf('Matrice calcolata con successo su %d campionamenti temporali.\n', size(rendimenti, 1));
end