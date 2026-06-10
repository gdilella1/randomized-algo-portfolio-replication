function [time, errfro, err2] = test_CUR_large(sz, ranks, algos, tag)
%%
%   A: (m,n) matrix
%   k: true rank of decomposition (i.e., denoted as l = target rank + O(1) in the manuscript)

    %%
    if isempty(sz)
        % default target: (1e6, 1e6) sparse non-negative matrix of rank r = 1000
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
    disp('Target generated')
    
    if ~exist('algos','var') || isempty(algos)
        algos = {'SRCUR',...
                 'CPQR',...
                 'CPQR2pass',...
                 ...
                 'LUPP',...
                 'LUPP2pass',...
                 ...
                 'RSVDDEIM',...
                 ...
                 'RSVDLS'};
    end
    
    if ~exist('tag','var') || isempty(tag)
        tag = 'aux';
    end
    
    save_target(spec, tag);
    %% Initialization
    rank_file = sprintf('rank_%s.mat',tag);
    if isfile(rank_file)
        k = load(rank_file); k = k.k;
        warning('Rank file exists. Using existing ranks.')
    else
        k = ranks;
    end
    save(rank_file,'k')
    fprintf('write out: %s \n', rank_file)
    
    err_file = sprintf('errfro_%s.mat', tag);
    if isfile(err_file)
        errfro = load(err_file);
    else
        errfro = struct();
    end
    
    err2_file = sprintf('err2_%s.mat', tag);
    if isfile(err2_file)
        err2 = load(err2_file);
    else
        err2 = struct();
    end
    
    time_file = sprintf('time_%s.mat',tag);
    if isfile(time_file)
        time = load(sprintf('time_%s.mat',tag));
    else
        time = struct();
    end
    
    %% Experiments
    for idx = 1:length(algos)
        algo = algos{idx};
        fprintf('%s \n', algo)
        time.(algo) = zeros(size(k));
        errfro.(algo) = zeros(size(k));
        err2.(algo) = zeros(size(k));
        
        for t = 1:length(k)
            tic;
            [i,j] = cur_algos(algo, A, k(t));
            time.(algo)(t) = toc;
            
            fprintf('k = %d: %.4f\n', k(t), time.(algo)(t))
        end
        

        for t = 1:length(k)
            % A = AL * Proj(CR) * Proj(RL) * AR
            CR = AR(:,j(1:k(t))); % (r,l)
            [Qcr,~] = qr(full(CR),0); % (r,l)
            RL = AL(i(1:k(t)),:); % (l,r)
            [Qrl,~] = qr(full(RL'),0); % (r,l)
            
            Ecore = eye(r) - (Qcr * (Qcr' * Qrl)) * Qrl'; % (r,r)
            Euinv = (TL * Ecore) * TR'; % (r,r)
            errfro.(algo)(t) = norm(Euinv, 'fro');
            err2.(algo)(t) = norm(Euinv);
            fprintf('%d / %d\t', t, length(k));
        end
        fprintf('\n')

        
        
        save_time(time, tag);
        save_err(errfro, tag,'fro')
        save_err(err2, tag,'2')
    end

end


function save_target(spec, tag)
    file = sprintf('target-spec_%s.mat',tag);
    out = struct('tag',tag,...
                 'spec',spec);
    save(file,'-struct','out')
    fprintf('write out: %s \n', file)
end


function save_time(time, tag)
    time_file = sprintf('time_%s.mat',tag);
    save(time_file,'-struct','time')
    fprintf('write out: %s \n', time_file)
end

function save_err(err, tag, errmeas)
    err_file = sprintf('err%s_%s.mat',errmeas,tag);
    save(err_file,'-struct','err')
    fprintf('write out: %s \n', err_file)
end

