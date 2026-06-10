function [time, errfro, err2] = test_CUR_large1(sz, ranks, algos, tag)
    %% 1. Configurazione e Generazione Matrice (Logica Autori)
    if isempty(sz)
        m = 1e6; n = m; 
        tag = 'snn-1e6-1e6-a2b1-k100-r400-s2or';
    elseif length(sz) < 2
        m = sz(1); n = m;
    else 
        m = sz(1); n = sz(2);
    end
    
    r = 400; s = 2/r;
    r0 = 100; a = 2; b = 1;
    X = sprand(m,r,s); Y = sprand(r,n,s); 
    d = 2.^(-(1:r).*(16/r));
    d(1:r0) = d(1:r0)*a; d(r0+1:r) = d(r0+1:r)*b;
    AL = (X.*d); AR = Y;
    A = {AL, AR};
    [~,TL] = qr(full(AL),0);
    [~,TR] = qr(full(AR'),0);
    spec = svd(TL*TR');
    
    % --- NECESSARIO PER ERRORE RELATIVO ---
    norm_A_fro = norm(spec);
    norm_A_2 = max(spec);
    % --------------------------------------
    
    disp('Target generated')
    
    if ~exist('algos','var') || isempty(algos)
        algos = {'SRCUR','CPQR','CPQR2pass','LUPP','LUPP2pass','RSVDDEIM','RSVDLS'};
    end
    if ~exist('tag','var') || isempty(tag), tag = 'aux'; end
    
    k = ranks;
    time = struct(); errfro = struct(); err2 = struct();
    
    %% 2. Esecuzione Algoritmi (Runtime e Accuratezza)
    for idx = 1:length(algos)
        algo = algos{idx};
        fprintf('%s \n', algo)
        time.(algo) = zeros(size(k));
        errfro.(algo) = zeros(size(k));
        err2.(algo) = zeros(size(k));
        
        % Misurazione Tempo
        for t = 1:length(k)
            tic;
            [i,j] = cur_algos(algo, A, k(t));
            time.(algo)(t) = toc;
            fprintf('k = %d: %.4f\n', k(t), time.(algo)(t))
        end
        
        % Calcolo Errore
        for t = 1:length(k)
            CR = AR(:,j(1:k(t)));
            [Qcr,~] = qr(full(CR),0);
            RL = AL(i(1:k(t)),:);
            [Qrl,~] = qr(full(RL'),0);
            
            Ecore = eye(r) - (Qcr * (Qcr' * Qrl)) * Qrl';
            Euinv = (TL * Ecore) * TR';
            
            % --- ERRORE RELATIVO ---
            errfro.(algo)(t) = norm(Euinv, 'fro') / norm_A_fro;
            err2.(algo)(t) = norm(Euinv) / norm_A_2;
            % -----------------------
            fprintf('%d / %d\t', t, length(k));
        end
        fprintf('\n')
    end

    %% 3. Sezione Plotting (Sostituisce i salvataggi)
    % Usiamo una chiamata robusta a figure per evitare l'errore degli "input arguments"
    hFig = figure; 
    set(hFig, 'Name', ['Results: ' char(tag)], 'Position', [100, 100, 1400, 480], 'Color', 'w');
    
    colors = lines(length(algos));
    
    % (A) Frobenius norm error
    subplot(1,3,1); hold on;
    for idx = 1:length(algos)
        semilogy(k, errfro.(algos{idx}), '-', 'Color', colors(idx,:), 'LineWidth', 1.5, 'DisplayName', algos{idx});
    end
    configura_asse_log('||A - CUR||_F / ||A||_F', 'Frobenius norm error');

    % (B) Spectral norm error
    subplot(1,3,2); hold on;
    for idx = 1:length(algos)
        semilogy(k, err2.(algos{idx}), '-', 'Color', colors(idx,:), 'LineWidth', 1.5, 'DisplayName', algos{idx});
    end 
    configura_asse_log('||A - CUR||_2 / ||A||_2', 'Spectral norm error');

    % (C) Runtime
    subplot(1,3,3); hold on;
    for idx = 1:length(algos)
        semilogy(k, time.(algos{idx}), '-', 'Color', colors(idx,:), 'LineWidth', 1.5, 'DisplayName', algos{idx});
    end
    configura_asse_log('Time (s)', 'Runtime');
    
    legend('Location', 'northeastoutside', 'Interpreter', 'none');
end

% --- Sottofunzione per la formattazione logaritmica richiesta ---
function configura_asse_log(y_label_text, title_text)
    grid on;
    set(gca, 'YScale', 'log', 'YMinorGrid', 'on', 'TickDir', 'in');
    xlabel('k'); ylabel(y_label_text); title(title_text);
    
    % Forza i limiti per vedere le tacche 10^x
    yl = ylim;
    l_min = floor(log10(yl(1)));
    l_max = ceil(log10(yl(2)));
    if l_min == l_max, l_min = l_min - 1; l_max = l_max + 1; end
    ylim([10^l_min, 10^l_max]);
    
    % Genera etichette pulite (es. 10^0 invece di 10^-0)
    yticks_vec = 10.^(l_min:l_max);
    set(gca, 'YTick', yticks_vec);
    set(gca, 'YTickLabel', arrayfun(@(x) sprintf('10^{%d}', round(log10(x))), yticks_vec, 'UniformOutput', false));
end