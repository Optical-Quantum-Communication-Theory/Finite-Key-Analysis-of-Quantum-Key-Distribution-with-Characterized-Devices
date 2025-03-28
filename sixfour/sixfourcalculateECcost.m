%% 
% This function is to calculate EC cost for the six-four protocol.
% Note that the ECcost returned by this function 
% does not include fEC.

function ECcost = sixfourcalculateECcost(simulated_state, AlicePOVM,BobPOVM, options)
    dim = 2;
    nElms = numel(AlicePOVM);
    mElms = numel(BobPOVM);
    key_generation_freqs= zeros(nElms,mElms);
    for i =1:nElms
        for j =1:mElms
             key_generation_freqs(i,j) = trace(kron(AlicePOVM{i},BobPOVM{j})*simulated_state);
        end
    end
    key_generation_freqs = real(key_generation_freqs);
    
    if options.ZbasisOnly
            
            key_generation_freqs = key_generation_freqs(1:dim,1:dim);
            pSift = sum(sum(key_generation_freqs));
            key_generation_freqs = key_generation_freqs ./ pSift;
            ECcost = pSift * calculateEC(key_generation_freqs);
    else
        % both Z basis and X basis
            key_generation_freqs(1:dim,dim+1:2*dim)=zeros(dim,dim);
            key_generation_freqs(dim+1:2*dim,1:dim)=zeros(dim,dim);
            key_generation_freqs(dim*2+1:3*dim,1:dim)=zeros(dim,dim);
            key_generation_freqs(dim*2+1:3*dim,dim+1:2*dim)=zeros(dim,dim);
            pSift = sum(sum(key_generation_freqs));
            key_generation_freqs = key_generation_freqs ./ pSift;
            ECcost = pSift * calculateEC(key_generation_freqs);
    end
end
  