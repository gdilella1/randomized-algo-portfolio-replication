function [U,s,V] = RSVD(A,k,randmat,l)
%   A = [m*n matrix]
%   k << min(m,n) sample rank
%   randmat: [str] {'gauss','srft','sparse(n)'}
%   U = [m*k]
%   s = [k*1]
%   V = [n*k]

    [m,n] = size(A);
    if nargin<3, randmat = 'gauss'; end
    if nargin<4, l = min([k+5,m,n]); end
    l = max(l,k);
    
    S = embed(n, l, randmat);
    Y = S(A')'; % (m,l)
    [Q,~] = qr(full(Y),0); % (m,l)
    [Ul,Sl,Vl] = svd(Q'*A,'econ'); % Ul:(l,l), V:(n,l)
    U = Q*Ul(:,1:k); % (m,k)
    sl = diag(Sl); s = sl(1:k);
    V = Vl(:,1:k);
end

   