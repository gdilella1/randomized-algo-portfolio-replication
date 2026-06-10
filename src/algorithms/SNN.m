function A = SNN(type)
    % Genera le matrici SNN basate sulle specifiche del paper
    
    switch type
        case 'SNN1e3'
            m = 1000; n = 1000; r = 1000; % Dimensioni SNN1e3
            density = 0.015; 
        case 'SNN1e6'
            m = 1e6; n = 1e6; r = 400;    % Dimensioni SNN1e6
            density = 0.00002; % Densità bassissima necessaria per scale 1e6
        otherwise
            error('Usa SNN1e3 o SNN1e6.');
    end

    % 1. Definizione dei coefficienti s_i
    s = (1:r)';
    s(1:100) = 2 ./ s(1:100);    % s_i = 2/i per i=1..100
    if r > 100
        s(101:r) = 1 ./ s(101:r); % s_i = 1/i per i > 100
    end

    % 2. Generazione vettori sparsi non-negativi x_i e y_i
    X = sprand(m, r, density);
    Y = sprand(n, r, density);

    % 3. Costruzione della matrice A = sum(s_i * x_i * y_i')
    if strcmp(type, 'SNN1e6')
        warning('Attenzione: Generare A intera per 1e6 richiede troppa memoria. Restituisco i fattori.');
        A = {X .* sqrt(s'), (Y .* sqrt(s'))'};
    else
        % Per SNN1e3 possiamo calcolare il prodotto completo
        A = (X .* s') * Y'; 
    end
end