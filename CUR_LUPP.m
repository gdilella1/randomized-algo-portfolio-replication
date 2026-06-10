function [I,J] = CUR_LUPP(A,k,sketch,stream,power,ortho)
    if ~exist('sketch','var')
        sketch = [];
    end
    if isempty(sketch) || ~exist('stream','var') || isempty(stream)
        stream = [];
    end
    if isempty(sketch) || ~isempty(stream) || ~exist('power','var') || isempty(power)
        power = 0;
    end
    if ~exist('ortho','var')
        ortho = 0;
    end
    
    m = size(A,1); n = size(A,2);
    os = 10;
    l = min([m,n,k+os]);
    
    if isempty(sketch)
        %% deterministic
        [~,~,I] = lu(A, 'vector');
        [~,~,J] = lu(A(I(1:k),:)', 'vector');
    
    elseif isempty(stream)
        %% col-sampling
        Sy = embed(n,l,sketch);
        Y = Sy(A')'; %(m,l)
        if power==0
            [~,~,I] = lu(Y, 'vector');
            [~,~,J] = lu(A(I(1:k),:)','vector');
        else
            if ortho==0
                Z = A'*Y; %(n,l)
            else
                [Q,~] = qr(Y,0);
                Z = A'*Q; %(n,l)
            end
            [~,~,J] = lu(Z,'vector');
            [~,~,I] = lu(A(:,J(1:k)),'vector');
        end
        
    else
        %% row+col sampling
        Sy = embed(n,l,sketch);
        Sx = embed(m,l,sketch);
        Y = Sy(A')'; %(m,l)
        X = Sx(A)'; % (n,l);
        [~,~,I] = lu(Y, 'vector');
        [~,~,J] = lu(X, 'vector');
    end
        
end