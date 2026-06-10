function [i,j,U] = srcur(A,l,f)
    if ~exist('f','var') || isempty(f)
        f = 10;
    end
    
    [i,j,LU] = fptlu(A,l);
    [i,j,~] = fsrp(i,j,LU,l,f);
    
    if nargout==3
        C = A(:,j(1:l));
        R = A(i(1:l),:);
        U = CUR_U_eval(A,C,R);
    end
end

function [i,j,LU] = fsrp(i,j,LU,l,f)
    m = size(LU,1); n = size(LU,2);
    maxiter = ceil(l*log(m*n));
    p2 = min(m-l, ceil(log(m-l))+5);
    p3 = min(l+1, ceil(log(l+1))+5);
    S1 = randn(p2, m-l);
    S2 = randn(p3, l+1);
    
    S = LU(l+1:m,l+1:n);
    R = S1*S;
    [sr,sc,alpha] = eelm(S,R);
    [LU,i,j] = swap(LU,i,j,l+1,sc,sr);
    L11aug = tril(LU(1:l+1,1:l+1),-1) + eye(l+1);
    U11aug = triu(LU(1:l+1,1:l+1));
    IA11augT = (U11aug\(L11aug\eye(l+1)))';
    W = S2*IA11augT;
    [ar,ac,beta] = eelm(IA11augT, W);
    
    for iter = 1:maxiter
        if alpha*beta > f
            break;
        end
        
        [LU,i,j] = swap(LU,i,j,l+1,ac,ar);
        S = LU(l+1:m,l+1:n);
        R = S1*S;
        [sr,sc,alpha] = eelm(S,R);
        [LU,i,j] = swap(LU,i,j,l+1,sc,sr);
        L11aug = tril(LU(1:l+1,1:l+1),-1) + eye(l+1);
        U11aug = triu(LU(1:l+1,1:l+1));
        IA11augT = (U11aug\(L11aug\eye(l+1)))';
        W = S2*IA11augT;
        [ar,ac,beta] = eelm(IA11augT, W);
    end
end


function [i,j,LU] = fptlu(A,l)
    m = size(A,1); n = size(A,2);
    p1 = min(n, ceil(log(n))+5);
    S1 = randn(p1,m);
    B = S1*A;
    i = 1:m; j = 1:n;
    LU = A;
    
    [r,c,~] = eelm(A,B); 
    [LU,i,j] = swap(LU,i,j,1,c,r);
    B = sketch_update(B,S1,LU,1);
    for k = 1:l
        [r,c,~] = eelm(A,B); 
        [LU,i,j] = swap(LU,i,j,k,c,r);
        LU(k:m,k) = LU(k:m,k) - LU(k:m,1:k-1)*LU(1:k-1,k);
        LU(k,k+1:n) = LU(k,k+1:n) - LU(k,1:k-1)*LU(1:k-1,k+1:n);
        B = sketch_update(B,S1,LU,k);
    end
end

function [r,c,x] = eelm(A,R)
% Inputs
%   A: (m,n), R: (p,n)
%   r: row pivot
%   c: col pivot
%   x: A(r,c)

    [~,c] = max(sum(R.^2,1));
    [x,r] = max(A(:,c));
end

function B = sketch_update(B,S,LU,k)
%   B:(p1,n); S:(p1,m); LU:(m,n); 
    n = size(B,2); m = size(S,2);
    B(:,k+1:n) = B(:,k+1:n) - ( S(:,k)+S(:,k+1:m)*LU(k+1:m,k) ) * LU(k,k+1:n);
end

function [A,i,j] = swap(A,i,j,k,c,r)
if exist('c','var') && ~isempty(c)
% swap row and update i
    aux = A(k,:);
    A(k,:) = A(r,:);
    A(r,:) = aux;
    aux = i(k); i(k) = i(r); i(r) = aux;  
end

if exist('r','var') && ~isempty(r)
% sway col and update j
    aux = A(:,k);
    A(:,k) = A(:,c);
    A(:,c) = aux;
    aux = j(k); j(k) = j(c); j(c) = aux;
end
end