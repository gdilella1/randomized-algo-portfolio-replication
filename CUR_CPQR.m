function [i,j,U] = CUR_CPQR(A,k,sketch,piter,ortho)
%   Inputs:
%       A = [m*n matrix]
%       k = [int] sample size = O(rank(A))
%       randmat = {'gauss','fsrft','srft','gsrft','sparse(n)'}
%       q = [int] number of power iterations
%       ortho = [bool] 1 = orthogonalize
%   Outputs:
%       U = [k*k matrix]
%       ik, jk = [1*k row / column indices]
%       i = [1*m row idx]; j = [1*n col idx]
%       C = A(:,j(1:k)); R = A(i(1:k),:) 

    if ~exist('sketch','var')
        sketch = 'gauss';
    end
    if ~exist('piter','var') || isempty(piter)
        piter = 0;
    end
    if ~exist('ortho','var') || isempty(ortho)
        ortho = 0;
    end
    [m,n] = size(A);
    
    if isempty(sketch)
        [~,~,j] = qr(full(A),0);
        C = A(:,j(1:k)); %(m,k)
        [~,~,i] = qr(full(C'),0);
        
    else
        l = min(k+10, min(m,n));

        S_col = embed(m, l, sketch);
        X = S_col(A); % (l,n)
        if piter==0
            [~,~,j] = qr(full(X),0);
            C = A(:,j(1:k)); %(m,k)
            [~,~,i] = qr(full(C'),0);
        else
            if ortho==0
                Y = X*A'; %(l,m)
                [~,~,i] = qr(full(Y),0);
            else
                [Q,~] = qr(X',0); %(n,l)
                Y = (A*Q)'; %(l,m)
                [~,~,i] = qr(full(Y),0);
            end
            R = A(i(1:k),:); %(k,n)
            [~,~,j] = qr(full(R),0);
        end
    end
    
    if nargout == 3
        C = A(:,j(1:k));
        R = A(i(1:k),:);
        U = CUR_U_eval(A,C,R);
    end
    
end


