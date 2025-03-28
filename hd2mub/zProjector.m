%% FUNCTION NAME: zProjector
% It outpus a projection operator \ket{j}\bra{j} in the computational basis (z basis)
% when the index = j, counting from 1 to dimension = dim.
function outputProjector = zProjector(dim, index )
    idMatrix=eye(dim);
    jColumn=idMatrix(:,index);
    outputProjector= jColumn * jColumn';
end