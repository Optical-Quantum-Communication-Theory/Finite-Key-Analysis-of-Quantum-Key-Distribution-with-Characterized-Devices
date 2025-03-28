% find the inner product of two vectors y, A, assuming they have the same length. 
function z = inner_prod(y,A)
    z = 0;
    for i = 1:length(A)
        z = z + y(i)*A{i};
    end
end