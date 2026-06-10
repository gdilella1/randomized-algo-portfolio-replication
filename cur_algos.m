function [i,j] = cur_algos(algo, A, l)
%   {'SRCUR',...
%    'LUPP', 'LUPP2pass', 'CSLUPP', 'LUPPstream',...
%    'CPQR', 'CPQR2pass', 'CSCPQR', 'CPQRstream',...
%    'RSVDDEIM',...
%    'RSVDLS'}
    if iscell(A)
        AL = A{1}; % (m,r)
        AR = A{2}; % (r,n)
        m = size(AL,1); n = size(AR,2);
        oracle = 1;
    else
        m = size(A,1); n = size(A,2);
        oracle = 0;
    end
%%
    switch algo
        % srcur
        case 'SRCUR'
            if oracle
                [i,j] = srcur(AL*AR, l);
            else
                [i,j] = srcur(A, l);
            end
        % lupp
        case 'LUPP'
            sketch = 'gauss';
            Smpr = embed(m, l, sketch);
            if oracle
                XL = Smpr(AL); % (l,r)
                X = XL * AR; % (l,n)
                [~,~,j] = lu(full(X'),'vector'); 
                j = j(1:l);
                C = AL * AR(:,j); %(m,l)
                [~,~,i] = lu(full(C),'vector');
                i = i(1:l);
            else
                X = Smpr(A); % (l,n)
                [~,~,j] = lu(full(X'),'vector'); 
                j = j(1:l);
                C = A(:,j); %(m,l)
                [~,~,i] = lu(full(C),'vector');
                i = i(1:l);
            end
        case 'LUPP2pass'
            sketch = 'gauss';
            Smpr = embed(n, l, sketch);
            if oracle
                Y = AL * Smpr(AR')'; % (m,l)
                X = AR' * (AL' * Y); %(n,l)
                [~,~,j] = lu(full(X),'vector');
                j = j(1:l);
                C = AL * AR(:,j); %(m,l)
                [~,~,i] = lu(full(C),'vector');
                i = i(1:l);
            else
                Yt = Smpr(A'); % (l,m)
                X = Yt*A; %(l,n)
                [~,~,j] = lu(full(X'),'vector');
                j = j(1:l);
                C = A(:,j); %(m,l)
                [~,~,i] = lu(full(C),'vector');
                i = i(1:l);
            end
        case 'CSLUPP'
            sketch = 'sparse8';
            Smpr = embed(m, l, sketch);
            if oracle
                XL = Smpr(AL); % (l,r)
                X = XL * AR; % (l,n)
                [~,~,j] = lu(full(X'),'vector'); 
                j = j(1:l);
                C = AL * AR(:,j); %(m,l)
                [~,~,i] = lu(full(C),'vector');
                i = i(1:l);
            else
                X = Smpr(A); % (l,n)
                [~,~,j] = lu(full(X'),'vector'); 
                j = j(1:l);
                C = A(:,j); %(m,l)
                [~,~,i] = lu(full(C),'vector');
                i = i(1:l);
            end
        case 'LUPPstream'
            sketch = 'gauss';
            Smpr_m = embed(m, l, sketch);
            Smpr_n = embed(n, l, sketch);
            if oracle
                X = Smpr_m(AL) * AR; % (l,n)
                [~,~,j] = lu(full(X'),'vector'); 
                j = j(1:l);
                Y = AL * Smpr_n(AR')'; %(m,l)
                [~,~,i] = lu(full(Y),'vector');
                i = i(1:l);
            else
                X = Smpr_m(A); % (l,n)
                Yt = Smpr_n(A'); % (l,m)
                [~,~,j] = lu(full(X'),'vector');
                [~,~,i] = lu(full(Yt'),'vector');
                j = j(1:l);
                i = i(1:l);
            end
        % cpqr
        case 'CPQR'
            sketch = 'gauss';
            Smpr = embed(m, l, sketch);
            if oracle
                XL = Smpr(AL); % (l,r)
                X = XL * AR; % (l,n)
                [~,~,j] = qr(full(X),0);
                j = j(1:l);
                C = AL * AR(:,j); %(m,l)
                [~,~,i] = qr(full(C'),0);
                i = i(1:l);
            else
                X = Smpr(A); % (l,n)
                [~,~,j] = qr(full(X),0);
                j = j(1:l);
                C = A(:,j); %(m,l)
                [~,~,i] = qr(full(C'),0);
                i = i(1:l);
            end          
        case 'CPQR2pass'
            sketch = 'gauss';
            Smpr = embed(n, l, sketch);
            if oracle
                Y = AL * Smpr(AR')'; % (m,l)
                X = AR' * (AL' * Y); %(n,l)
                [~,~,j] = qr(full(X'),0);
                j = j(1:l);
                C = AL * AR(:,j); %(m,l)
                [~,~,i] = qr(full(C'),0);
                i = i(1:l);
            else
                Yt = Smpr(A'); % (l,m)
                X = Yt*A; %(l,n)
                [~,~,j] = qr(full(X'),0);
                j = j(1:l);
                C = A(:,j); %(m,l)
                [~,~,i] = qr(full(C'),0);
                i = i(1:l);
            end
        case 'CSCPQR'
            sketch = 'sparse8';
            Smpr = embed(m, l, sketch);
            if oracle
                XL = Smpr(AL); % (l,r)
                X = XL * AR; % (l,n)
                [~,~,j] = qr(full(X),0);
                j = j(1:l);
                C = AL * AR(:,j); %(m,l)
                [~,~,i] = qr(full(C'),0);
                i = i(1:l);
            else
                X = Smpr(A); % (l,n)
                [~,~,j] = qr(full(X),0);
                j = j(1:l);
                C = A(:,j); %(m,l)
                [~,~,i] = qr(full(C'),0);
                i = i(1:l);
            end        
        case 'CPQRstream'
            sketch = 'gauss';
            Smpr_m = embed(m, l, sketch);
            Smpr_n = embed(n, l, sketch);
            if oracle
                X = Smpr_m(AL) * AR; % (l,n)
                [~,~,j] = qr(full(X),0);
                j = j(1:l);
                Y = AL * Smpr_n(AR')'; %(m,l)
                [~,~,i] = qr(full(Y'),0);
                i = i(1:l);
            else
                X = Smpr_m(A); % (l,n)
                Yt = Smpr_n(A'); % (l,m)
                [~,~,j] = qr(full(X),0);
                [~,~,i] = qr(full(Yt),0);
                j = j(1:l);
                i = i(1:l);
            end
        % deim
        case 'RSVDDEIM'
            if oracle
                [V,~,W] = rsvd_oracle(AL,AR,l);
            else
                [V,~,W] = RSVD(A,l);
            end
            [~,~,i] = lu(V,'vector');
            [~,~,j] = lu(W,'vector');
        % leverage score sampling
        case 'RSVDLS'
            if oracle
                [V,~,W] = rsvd_oracle(AL,AR,l);
            else
                [V,~,W] = RSVD(A,l);
            end
            i = LeverageScoreSampling(V); 
            j = LeverageScoreSampling(W); 
        otherwise 
            warning('The CUR algorithm %s is not found!',algo);
            i = 1:l; j = 1:l;
    end    
end

function [U,S,V] = rsvd_oracle(AL,AR,l)
    n = size(AR,2);
    sketch = 'gauss';
    Smpr = embed(n, l, sketch);
    Y = AL * Smpr(AR')'; % (m,l)
    [Q, ~] = qr(full(Y),0); %aggiustato

    B = (Q' * AL) * AR;
    [Uh,S,V] = svd(B,'econ');
    U = Q * Uh;
end

function p = LeverageScoreSampling(V)
%   V = [m*k matrix] top k left / right singular vectors
%   p = [1*k] row skeleton indices selected by LS
    ls = sum(V.^2,2)./size(V,2);
%     p = sampleProb((1:size(V,1))', ls, size(V,2));
    [~,p] = maxk(ls, size(V,2));
end