function benchmarkSketching(A, m_vec, l_vec)
    % A     : Matrice sorgente (m_max x n)
    % m_vec : Vettore delle dimensioni m da testare (es. logspace(2, 5, 10))
    % l_vec : Vettore dei valori di k (target dimension l)
    
    % Configurazione stili per matchare la foto
    tipi = {'GAUSS', 'SRTT', 'SSM'};
    nomi_legenda = {'Gaussian', 'SRFT', 'Sparse Sign \zeta = 8'};
    colori = {'r', 'g', 'b'};
    marker = {'o', 'd', 's'}; % Cerchio, Diamante, Quadrato
    
    zeta_fisso = 8; % Valore preso dal grafico di riferimento
    n_subplot = length(l_vec);
    n_cols = size(A, 2);
    
    figure('Color', 'w', 'Position', [100, 100, 1200, 400]);
    
    for i = 1:n_subplot
        l_curr = l_vec(i);
        subplot(1, n_subplot, i);
        hold on;
        
        for t = 1:length(tipi)
            tempi = zeros(size(m_vec));
            
            for j = 1:length(m_vec)
                m_curr = round(m_vec(j));
                A_sub = A(1:m_curr, :);
                
                % Misurazione del tempo totale (generazione + prodotto)
                Omega = scegliSketch(tipi{t}, l_curr, m_curr, zeta_fisso);
                % Calcolo dello sketch
                tic;
                Y = Omega * A_sub; 
                tempi(j) = toc;
            end
            
            % Plot in scala logaritmica
            semilogy(m_vec, tempi, 'Color', colori{t}, 'Marker', marker{t}, ...
                   'LineStyle', '-', 'MarkerSize', 4, 'LineWidth', 1);
        end
        
        % Formattazione grafico
        grid on;
        xlabel('m');
        ylabel('Time (s)');
        title(['k = ', num2str(l_curr), ', n = ', num2str(n_cols)]);
        
        if i == 1
            legend(nomi_legenda, 'Location', 'northwest', 'FontSize', 9);
        end
        
        % Opzionale: imposta i limiti per coerenza visiva
        % ylim([1e-6, 1e2]);
    end
end