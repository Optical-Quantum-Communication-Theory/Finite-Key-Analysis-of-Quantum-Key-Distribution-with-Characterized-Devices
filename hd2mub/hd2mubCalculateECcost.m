%% 
% This function is to calculate EC cost for high-dimensional 2-MUB protocols.
% Note that the ECcost returned by this function 
% does not include fEC.
function ECcost = hd2mubCalculateECcost(dim, simulated_state, AlicePOVM,BobPOVM, options) 
  
    nElms = numel(AlicePOVM);
    mElms = numel(BobPOVM);
    key_generation_freqs= zeros(nElms,mElms);
    for i =1:nElms
        for j =1:mElms
             key_generation_freqs(i,j) = trace(kron(AlicePOVM{i},BobPOVM{j})*simulated_state);
        end
    end
    key_generation_freqs = real(key_generation_freqs);
    ECcost = twoMUBECcost(key_generation_freqs, dim, options);
end
