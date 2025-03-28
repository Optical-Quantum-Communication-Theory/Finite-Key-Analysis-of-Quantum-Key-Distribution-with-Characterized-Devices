%% This function generates Z basis and X basis projections for high-dimensional QKD.
% Make sure everything is Hermitian
function projectors = gethdONBProjectors(d,isAlice)
    Id = eye(d);
    
    zBasisProjectors = cell(d,1);
    xBasisProjectors = cell(d,1);
    F = fourier_operator(d);
    if isAlice
        for r = 1:d
            zBasisProjectors{r} = Id(:,r)*Id(r,:);
            xBasisProjectors{r} = F(:,r)*F(:,r)';
            xBasisProjectors{r} = (xBasisProjectors{r}+ xBasisProjectors{r}')/2;
        end
        
    else
         for r = 1:d
            zBasisProjectors{r} = Id(:,r)*Id(r,:);
            xBasisProjectors{r} =  conj(F(:,r))*F(:,r).';
            xBasisProjectors{r} = (xBasisProjectors{r}+ xBasisProjectors{r}')/2;
        end
        
    end
    
    
    projectors = [zBasisProjectors;xBasisProjectors];
    
end