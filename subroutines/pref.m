% add a prefactor to each element in a cell of vectors/matrices
function C = pref(a, B)
    for i = 1:length(B)
        C{i} = a*B{i};
    end
end