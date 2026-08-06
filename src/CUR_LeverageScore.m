function [i,j,U] = CUR_LeverageScore(A,k,ls,stream,V,W)
%   Inputs:
%       A = [m*n matrix]
%       k = [int] sample size = O(rank(A))
%       V, W = left [m*n matrix] / right [n*n matrix] singular vectors
%       ls = 'full'    -> using SVD to find V and W
%          = {'gauss','srft','sparse(n)'} -> using RSVD to find V and W
%          = 'lsa' -> using LeverageScoreApproximation.m [Drineas2012]
%       stream = {0 (default),1}
%   Outputs:
%       U = [k*k matrix]
%       ik, jk = [1*k row / column indices]
%       i = [1*m row idx]; j = [1*n col idx]
%       C = A(:,j(1:k)); R = A(i(1:k),:) 
    if nargin < 3
        ls = 'gauss'; stream = 0;
    end
    if nargin < 4
        stream = 0;
    end
    
    if nargin >= 6
        i = LeverageScoreSampling(V(:,1:k)); 
        j = LeverageScoreSampling(W(:,1:k)); 
    elseif strcmpi(ls, 'lsa')
        eps = 1e-1;
        lsi = LeverageScoreApproximation(A, eps);
        lsj = LeverageScoreApproximation(A', eps);
        [~,i] = maxk(lsi, k);
        [~,j] = maxk(lsj, k);
    elseif strcmpi(ls, 'full')
        [V,~,W] = svd(full(A),'econ');
        i = LeverageScoreSampling(V(:,1:k)); 
        j = LeverageScoreSampling(W(:,1:k));
    else
        if stream
            l = 2*k; s = 2*l;
            [V,~,W] = stream_svd(A,k,ls,l,s);
        else
            l = k+5;
            [V,~,W] = RSVD(A,k,ls,l);
        end
        i = LeverageScoreSampling(V(:,1:k)); 
        j = LeverageScoreSampling(W(:,1:k)); 
    
    end
%%
    i = i(1:k); j = j(1:k);
    if nargout==3
        C = A(:,j); R = A(i,:);
        U = CUR_U_eval(A,C,R);
    end
end

function p = LeverageScoreSampling(V)
%   V = [m*k matrix] top k left / right singular vectors
%   p = [1*k] row skeleton indices selected by LS
    ls = sum(V.^2,2)./size(V,2);
%     p = sampleProb((1:size(V,1))', ls, size(V,2));
    [~,p] = maxk(ls, size(V,2));
end