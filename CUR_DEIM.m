function [i,j,U] = CUR_DEIM(A,k,randmat,stream,V,W)
%   Inputs:
%       A = [m*n matrix]
%       k = [int] sample size = O(rank(A))
%       V, W = left [m*n matrix] / right [n*n matrix] singular vectors
%       randmat = {'gauss' (default),'srft','sparse(n)','stream'}
%       stream = {0 (default),1}
%   Outputs:
%       i = [k row idx]; j = [k col idx]
%       C = A(:,j); R = A(i,:)
%       U = [(k,k) matrix]
    if nargin < 3 
        randmat = 0; stream = 0;
    end
    if nargin < 4
        stream = 0; 
    end
    if nargin < 6
        if randmat
            if stream
                l = 2*k; s = 2*l;
                [V,~,W] = stream_svd(A,k,randmat,l,s);
            else
                l = k+5;
                [V,~,W] = RSVD(A,k,randmat,l);
            end
        else
            [V,~,W] = svd(full(A),'econ');
        end
    end
    
%%
    j = LUselect(W(:,1:k)); i = LUselect(V(:,1:k));
    
    i = i(1:k); j = j(1:k);
    if nargout==3
        C = A(:,j(1:k)); R = A(i(1:k),:);
        U = CUR_U_eval(A,C,R);
    end
end

function p = LUselect(V)
%   V = [m*k matrix] top k left / right singular vectors
%   p = [1*m] perm vec with leading k row skeleton indices selected by DEIM
    [~,~,p] = lu(V,'vector');
%     p = p(1:size(V,2));
end

function p = DEIM(V)
%   V = [m*k matrix] top k left / right singular vectors
%   p = [1*m] perm vec with leading k row skeleton indices selected by DEIM
    [m,k] = size(V);
    v = V(:,1); 
    [~,p1] = max(abs(v));
    p = 1:m;
    p(1) = p1; p(p1) = 1;
    p1 = p(1); p2 = p(2:end);
    
    for j = 2:k
        c = V(p1,1:j-1)\V(p1,j);
        r = V(p2,j) - V(p2,1:j-1)*c;
        [~,pj] = max(abs(r));
        p(pj+j-1) = p(j);
        p(j) = p2(pj);
        
        p1 = p(1:j); 
        p2 = p(min(j+1,k):end);
    end
end