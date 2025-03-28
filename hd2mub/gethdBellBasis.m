%% This function generates the Bell basis in dimension d: {\ket{\Phi_{j,k}} = 1/sqrt{d} \sum_k \omega^{sk} \ket{s}\ket{s+j} \forall j,k}
%
function basis = gethdBellBasis(d)
    basis = cell(d,d);
    for j =1:d
        for k = 1:d
            basis{j,k} =  hdBellState(d,j-1,k-1);
        end
    end
end