function fval = entropy_fun(rho, Kp, Kpk)
    % objective function using von Neumann entropy 
    nP = numel(Kp);
    
    fval = 0;
    
    for p = 1:nP
        fval = fval - VNent( Kp{p}*rho*Kp{p}' );
        for k = 1:numel(Kpk{p})
            fval = fval + VNent( Kpk{p}{k} * rho * Kpk{p}{k}' );
        end
    end
    fval = real(fval);
end

function y = VNent( rho )
     % von Neumann entropy function of rho
    e = double(eig(rho));
    y = -e'*log2(e+(e==0));
end