%% This function generates the corresponding projectors to the Bell basis in dimension d: {\ket{U_{r,s}} = 1/sqrt{d} \sum_k \omega^{ks} \ket{k+r}\ket{k} \forall r, s}
%
function projectors = gethdBellStateProjectors(d)
    basis = gethdBellBasis(d);
    projectors = cell(d,d);
    for r =1:d
        for s = 1:d
            projectors{r,s} =  basis{r,s}*basis{r,s}';
        end
    end
end