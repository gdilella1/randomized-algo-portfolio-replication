function U = CUR_U_eval(A,C,R)
% Input
%   A: (m,n) matrix
%   C: (m,k) permutation vector
%   R: (l,n) permutation vector
% Output
%   U: (k,k) matrix in CUR
    if issparse(A)
        [B, T] = qr(C, A, 0); % B:(k,n), T:(k,k), T\B = least_sq(C\A):(k,n)
        [G, L] = qr(R', B', 0); % G:(l,m), L:(l,l), L\G = least_sq(R'\B'):(l,k)
        U = T\(L\G)'; % (k,l)
    else
        X = C\A; % X = [k*n]
        U = (R'\X')'; % U = [k*l]
    end
end