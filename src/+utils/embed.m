function S = embed(d, k, randmat)
% Input
%   d: [int] ambience dimension
%   k: [int] subspace dimension, k << d
%   randmat: [str] {'gauss','srft','sparse(n)'}
% Output
%   S: [anonymous func] A (d,n) -> S(A) (k,n) map of subspace embedding
%%
    if nargin < 3, randmat = 'gauss'; end
    
    if strcmpi(randmat,'srft')
        scale = exp(1i.*rand(d,1).*2.*pi);
        select = randperm(d,k);
        S = @(A) srft(A, select, scale) / sqrt(k);
        
    elseif length(randmat)>=length('sparse') && strcmpi(randmat(1:length('sparse')),'sparse')
        if length(randmat)==length('sparse')
           zeta = 3;
        else
           zeta = str2double(randmat(length('sparse')+1:length(randmat)));
        end
        Saux = spsign(d,k,zeta);
        S = @(A) Saux * A(1:d,:);
        
    else % if strcmpi(randmat,'gauss')
        G = randn(k,d)/sqrt(k);
        S = @(A) G * A(1:d,:);
    end
end

%%
function Y = srft(A, select, scale)
    aux = scale.*A;
    if issparse(aux), aux = full(aux); end
    aux = fft(aux);
    Y = aux(select,:);
end

%%
function S = spsign(d, k, zeta)
% Input
%   A: [d by n] -> S*A: [k by n]
%   k: [int] embedding dim, k << d
%   zeta: sparsity parameter, zeta >= 2
%       number of nonzeros in each of the d columns of S / 
%       number of hashes of each of the d coordinates / 
%       number of throws of each of the d balls into k bins
% Output
%   S: [(k,d) spmatrix] sparse sign matrix
%       nonzero = 1 or -1 w.p. 1/2 each 
    j = repelem(1:d, zeta);
    imat = zeros(zeta,d);
    for c = 1:d
        imat(:,c) = randperm(k,min(zeta,k))';
    end
    i = reshape(imat,1,[]);
    vals = (-1).^randi(2, 1, d*zeta);
    S = sparse(i, j, vals, k, d) / sqrt(zeta);
end
