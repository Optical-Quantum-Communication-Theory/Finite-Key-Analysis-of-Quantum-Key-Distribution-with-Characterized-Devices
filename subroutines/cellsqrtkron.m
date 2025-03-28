%% 
% Given two cells of matrices, 
% compute tensor product of square roots of each pair of combinations 
function C = cellsqrtkron( A,B )
    C = {};
    for i = 1:length(A)
        sqrtmA = sqrtm(A{i});
        sqrtmA = (sqrtmA+sqrtmA')/2;
        for j = 1:length(B)
            
            sqrtmB = sqrtm(B{j});
            sqrtmB = (sqrtmB+sqrtmB')/2;
            C{end+1} = kron( sqrtmA,sqrtmB );
        end
    end
end