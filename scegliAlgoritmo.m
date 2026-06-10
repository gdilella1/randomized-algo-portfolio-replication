function [I,J] = scegliAlgoritmo(algoritmo,A,k,sketch,stream,ortho)
if nargin < 5
    stream = [];
end
if nargin < 6
    ortho = 0;
end
ls = 'gauss';
    switch algoritmo
    case 'LUPP'
        [I,J] = CUR_LUPP(A,k,sketch,stream,0,ortho);
    case 'CPQR'
        [I,J, ~] = CUR_CPQR(A,k,sketch,0,ortho);
    case 'LUPP1'
        [I,J] = CUR_LUPP(A,k,sketch,stream,1,ortho);
    case 'CPQR1'
        [I,J, ~] = CUR_CPQR(A,k,sketch,1,ortho);
    case 'DEIM'
        [I,J, ~] = CUR_DEIM(A,k,sketch,stream);
    case 'LS'
        [I,J, ~] = CUR_LeverageScore(A,k,ls,stream);
    otherwise
        error("Algoritmo non riconosciuto.");
    end
end