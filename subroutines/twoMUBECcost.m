%% 
% This function is to calculate EC cost for protocols with two MUBs.
% Note that the ECcost returned by this function 
% does not include fEC.
function ECcost = twoMUBECcost(key_generation_freqs, dim, options) 
   
    if options.ZbasisOnly
         % use only Z basis to generate keys        
        key_generation_freqs = key_generation_freqs(1:dim,1:dim);
        % renormalize
        pSift = sum(sum(key_generation_freqs));
        key_generation_freqs = key_generation_freqs ./ pSift;
        ECcost =  pSift*calculateEC(key_generation_freqs);
    else
        % use both bases
        % remove mismatched basis choices
        key_generation_freqs(1:dim,dim+1:2*dim)=zeros(dim,dim);
        key_generation_freqs(dim+1:2*dim,1:dim)=zeros(dim,dim);
        % renormalize
        pSift = sum(sum(key_generation_freqs));
        key_generation_freqs = key_generation_freqs ./ pSift;
        ECcost =  pSift*calculateEC(key_generation_freqs);
    end
end
