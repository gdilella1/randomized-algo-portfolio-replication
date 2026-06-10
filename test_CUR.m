function test_CUR(A, k_vec, algos, varargin)
    %% 1. Gestione Input e Parametri Opzionali
    doRF = false;
    tag = 'Analisi_CUR';
    if nargin > 3
        for i = 1:length(varargin)
            if strcmp(varargin{i}, 'rangefinder')
                doRF = true;
            elseif ischar(varargin{i})
                tag = varargin{i}; 
            end
        end
    end

    [m, n] = size(A);
    [U, S, V] = svd(full(A), 'econ');
    AL = U * S; AR = V';

    %% 2. Calcolo spettro e norme di riferimento
    [~, TL] = qr(full(AL), 0);
    [~, TR] = qr(full(AR'), 0);
    spec = svd(TL * TR'); 
    norm_A_fro = norm(spec);    
    norm_A_2 = max(spec);       
    
    if ~exist('algos','var') || isempty(algos)
        algos = {'SRCUR','CPQR','LUPP','RSVDLS'}; 
    end
    
    time = struct(); errfro = struct(); err2 = struct();
    n_k = length(k_vec);

    %% 3. Esecuzione Algoritmi
    if doRF
        fprintf('Calcolo in corso: rangefinder \n')
        time.rangefinder = zeros(size(k_vec));
        errfro.rangefinder = zeros(size(k_vec));
        err2.rangefinder = zeros(size(k_vec));
        for t = 1:n_k
            k = k_vec(t);
            tic;
            A_approx = rangefinder(A, k, 'gauss');
            time.rangefinder(t) = toc;
            errfro.rangefinder(t) = norm(A - A_approx, 'fro') / norm_A_fro;
            err2.rangefinder(t) = norm(A - A_approx) / norm_A_2;
        end
    end

    for idx = 1:length(algos)
        algo = algos{idx};
        fprintf('Calcolo in corso: %s \n', algo)
        time.(algo) = zeros(size(k_vec));
        errfro.(algo) = zeros(size(k_vec));
        err2.(algo) = zeros(size(k_vec));
        
        for t = 1:n_k
            k = k_vec(t);
            tic;
            %[i, j] = scegliAlgoritmo(algo, A, k, []);
            [i, j] = scegliAlgoritmo(algo, A, k, 'gauss');
            time.(algo)(t) = toc; 
            
            CR = AR(:, j(1:k)); [Qcr, ~] = qr(full(CR), 0);
            RL = AL(i(1:k), :); [Qrl, ~] = qr(full(RL'), 0);
            
            Ecore = eye(size(AL,2)) - (Qcr * (Qcr' * Qrl)) * Qrl';
            Euinv = (TL * Ecore) * TR';
            
            errfro.(algo)(t) = norm(Euinv, 'fro') / norm_A_fro;
            err2.(algo)(t) = norm(Euinv) / norm_A_2;
        end
    end

    %% 4. Sezione Plotting
    figure('Name', ['Replication: ' tag], 'Position', [100, 100, 1400, 480], 'Color', 'w');
    colors = [
        1.00, 0.00, 0.00; % Rosso
        0.85, 0.33, 0.10; % Arancione
        0.00, 0.00, 1.00; % Blu
        0.00, 1.00, 1.00; % Ciano
        1.00, 0.00, 1.00; % Magenta
        0.00, 0.50, 0.00; % Verde
        0.00, 0.00, 0.00; % Nero (Rangefinder)
    ];

    titles = {'Frobenius norm error', 'Spectral norm error', 'Runtime'};
    ylabels = {'||A - CUR||_F / ||A||_F', '||A - CUR||_2 / ||A||_2', 'Time (s)'};
    data_fields = {errfro, err2, time};

    for p = 1:3
        subplot(1,3,p); hold on;
        current_data = data_fields{p};
        for idx = 1:length(algos)
            semilogy(k_vec, current_data.(algos{idx}), '-', 'Color', colors(idx,:), ...
                'LineWidth', 1.5, 'DisplayName', algos{idx});
        end
        if doRF
            semilogy(k_vec, current_data.rangefinder, '-', 'Color', colors(7,:), ...
                'LineWidth', 1.5, 'DisplayName', 'rangefinder');
        end
        
        % Elementi essenziali del plot inseriti direttamente
        grid on;
        xlabel('k');
        ylabel(ylabels{p});
        title(titles{p});
    end
    legend('Location', 'northeastoutside', 'Interpreter', 'none');
end