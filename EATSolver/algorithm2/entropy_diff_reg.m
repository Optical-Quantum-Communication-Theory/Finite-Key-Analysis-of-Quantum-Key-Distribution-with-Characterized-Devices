function dfval= entropy_diff_reg(rho, Kp, Kpk)
    % completely regularized formula for the derivative
    
    nP = numel(Kp);
    dfval = 0;
    
    for p = 1:nP
        dfval = dfval + dVNent( rho,Kp{p} );
        for k = 1:numel(Kpk{p})
            dfval = dfval - dVNent( rho,Kpk{p}{k} );
        end
    end
end

function grad = dVNent( rho,K )

    f = @(x)( x.*log2(x+(x==0)) );
    rhos = sqrtm(rho);
    [U,D] = eig( rhos*K*K'*rhos );
    grad = rhos\ ((U*diag(f(diag(D)))*U')/rhos );
    grad = full((grad + grad')/2);
end