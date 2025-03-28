%% 
% This function is to calculate EC cost for BB84-type protocols.
% Note that the ECcost returned by this function 
% does not include fEC.
function ECcost = opticalBB84CalculateECcost(key_generation_freqs, options) 
   
    if options.ZbasisOnly
            
        key_generation_freqs = key_generation_freqs(1:2,1:2);
        pSift = sum(sum(key_generation_freqs));
        key_generation_freqs = key_generation_freqs ./ pSift;
        ECcost =  pSift*calculateEC(key_generation_freqs);
    else
        % use both bases
        % remove mismatched basis choices
        key_generation_freqs(1:2,3:4)=zeros(2,2);
        key_generation_freqs(3:4,1:2)=zeros(2,2);
        pSift = sum(sum(key_generation_freqs));
        key_generation_freqs = key_generation_freqs ./ pSift;
        ECcost =  pSift*calculateEC(key_generation_freqs);
    end
end
