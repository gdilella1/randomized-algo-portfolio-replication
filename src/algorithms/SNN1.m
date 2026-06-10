function A = SNN1(sz)
    % sz: [m, n] oppure [] per il default 1e6 x 1e6
    
    %% 1. Definizione Dimensioni
    if isempty(sz)
        m = 1e6; n = m; 
    elseif length(sz) < 2
        m = sz(1); n = m;
    else 
        m = sz(1); n = sz(2);
    end
    
    %% 2. Parametri di Generazione (da codice originale)
    r = 400;            % Rango interno della fattorizzazione
    s = 2/r;            % Densità dello sprand (sparsity)
    r0 = 100; a = 2; b = 1;
    
    % Generazione componenti sparse non-negative
    X = sprand(m, r, s); 
    Y = sprand(r, n, s); 
    
    %% 3. Decadimento dei Pesi (Valori Singolari)
    % Questa è la formula esponenziale che permette errori molto bassi
    d = 2.^(-(1:r) .* (16/r));
    
    % Modifica dei primi r0 componenti
    d(1:r0) = d(1:r0) * a; 
    d(r0+1:r) = d(r0+1:r) * b;
    
    %% 4. Costruzione Output
    % Per gestire matrici 1e6 x 1e6, il codice restituisce 
    % obbligatoriamente la forma fattorizzata {AL, AR}
    AL = (X .* d); 
    AR = Y;
    
    if (m * n) > 1e7
        % Se la matrice è troppo grande per la RAM, restituisce i fattori
        A = {AL, AR};
    else
        % Per taglie piccole (es. 1000x1000) restituisce la matrice piena
        A = AL * AR;
    end
end