function Gamma = scegliSketch(type, l, m, zeta)
    % GET_SUBSPACE_EMBEDDING Genera matrici di embedding per RandNLA.
    %
    % Input:
    %   type : stringa ('GAUSS', 'SRTT', 'SSM')
    %   l    : dimensione target (numero di righe)
    %   m    : dimensione originale (numero di colonne)
    %   zeta : parametro di sparsità (richiesto solo per 'SSM')
    if nargin < 4
        zeta = 8;
    end

    switch upper(type)
        
        case 'GAUSS'
            % (1) Gaussian embeddings: S_ij ~ N(0, 1/l)
            % Dalla foto: Gamma in R^{l x m} con i.i.d. Gaussian entries.
            % Usiamo mu = 0 come da standard.
            Gamma = (1/sqrt(l)) * randn(l, m);

        case 'SRTT'
            % (2) Subsampled Randomized Trigonometric Transforms
            % Formula: Gamma = sqrt(m/l) * Pi_{m->l} * T * Phi * Pi_{m->m}
            
            % a. Pi_{m->m}: Matrice di permutazione casuale
            P = eye(m);
            P = P(randperm(m), :);
            
            % b. Phi: Matrice diagonale di segni casuali (Rademacher)
            phi = sign(randn(m, 1));
            Phi = diag(phi);
            
            % c. T: Matrice unitaria di Hartley (per R)
            % La trasformata di Hartley si ottiene da FFT: H = real(FFT) - imag(FFT)
            % Per renderla unitaria, dividiamo per sqrt(m)
            H = real(fft(eye(m))) - imag(fft(eye(m)));
            T = H / sqrt(m);
            
            % d. Pi_{m->l}: Selezione casuale di l righe su m
            indices = randperm(m, l);
            Pi_select = eye(m);
            Pi_select = Pi_select(indices, :);
            
            % Assemblaggio finale
            Gamma = sqrt(m/l) * (Pi_select * T * Phi * P);

        case 'SSM'
            % (3) Sparse Sign Matrices: Gamma = sqrt(m/zeta) * [s1, ..., sm]
            % Ogni colonna s_j ha zeta entrate Rademacher i.i.d. in posizioni casuali.
            if nargin < 4 || isempty(zeta)
                error('Per il tipo SSM è necessario fornire il parametro zeta (2 <= zeta <= l).');
            end
            
            rows = zeros(zeta * m, 1);
            cols = zeros(zeta * m, 1);
            vals = zeros(zeta * m, 1);
            
            for j = 1:m
                % Seleziona zeta coordinate casuali uniformi tra 1 e l
                pos = randperm(l, zeta);
                
                % Indici per la costruzione della matrice sparsa
                idx_start = (j-1)*zeta + 1;
                idx_end = j*zeta;
                
                rows(idx_start:idx_end) = pos;
                cols(idx_start:idx_end) = j;
                % Valori Rademacher (+1 o -1)
                vals(idx_start:idx_end) = sign(randn(zeta, 1));
            end
            
            % Creazione matrice sparsa per efficienza e applicazione scala
            Gamma = sqrt(m/zeta) * sparse(rows, cols, vals, l, m);
            % Converti in full se necessario (opzionale)
            % Gamma = full(Gamma);

        otherwise
            error('Tipo di embedding non riconosciuto. Usa GAUSS, SRTT o SSM.');
    end
end