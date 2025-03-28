%% Function name: EAT_step2
% Implementation of Eq. (47) of the paper with the replacement t by t/c_0


function [minTF,fval] = EAT_step2(protocol,problem,rho1)
    
    n_pe = length(protocol.POVM_pe);
    POVM_pe = protocol.POVM_pe;
    p_accept = problem.p_accept;
    Kp = protocol.Kp;
    Kpk = protocol.Kpk;
    dim = protocol.dim;
    n_pe = length(POVM_pe);
    pT = protocol.pT;
    
    % parameters of error term
    c0 = problem.c0;
    c1 = problem.c1;
    % optimize dual
    
    epsilon = 1e-12;
    
    cvx_begin sdp quiet
        cvx_precision best
        variable minTF(n_pe)
        variables t u v
        obj_fun = p_accept*minTF - c0*t;
        maximize obj_fun
            % constraint: t >= sqrt( 1 + c1^2*(u - v)^2)
            [t-1,       c1*(u-v);...
             c1*(u-v),  t+1] >= 0;
            for i = 1:n_pe
                u >= minTF(i);
                v <= minTF(i);
            end
            entropy_diff_reg(rho1,Kp,Kpk) - inner_prod(minTF, POVM_pe) - epsilon*eye(dim) == hermitian_semidefinite(dim);
    cvx_end
    
    % verifying if minTF is valid
    
    if  ((min(eig( entropy_diff_reg(rho1,Kp,Kpk) - inner_prod(minTF, POVM_pe) ))>= 0)...
        && (min(eig(rho1) >= 0)))
        fval = obj_fun;
    else
        fval = -inf;
    end
end