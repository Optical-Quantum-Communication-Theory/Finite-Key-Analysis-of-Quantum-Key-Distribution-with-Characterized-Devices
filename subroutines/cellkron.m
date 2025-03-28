%% 
% Given two cells of matrices, 
% compute tensor product of each pair of combinations 
function C = cellkron(A,B)
    C = {};
    for i = 1:length(A)
        Asym = (A{i}+A{i}')/2;
        for j = 1:length(B)
            Bsym = (B{j}+B{j}')/2;
            C{end+1} = kron( Asym,Bsym );
        end
    end
end