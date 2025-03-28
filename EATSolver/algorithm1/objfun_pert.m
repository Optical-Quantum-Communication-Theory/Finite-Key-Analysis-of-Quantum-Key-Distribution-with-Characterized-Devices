%% Function name: objfun_pert
%% Date of last modification: August 28, 2023
% Perturbed objective function

function [fval,epsilon_max] = objfun_pert(rho, Kp, Kpk)
  
    nP = numel(Kp);
    fval = 0;
    epsilon_max= 0;
    
    for iElm = 1: nP
        nPK = numel(Kpk{iElm});
        for jElm = 1: nPK
            rho2= Kpk{iElm}{jElm}*rho*Kpk{iElm}{jElm}';
            [rhoPrime,realEpsilon] = perturbation_channel(rho2);   
            epsilon_max = max(epsilon_max,realEpsilon);   
            fval = fval - real(trace(rhoPrime* logm(rhoPrime)));
        end
        rho3=Kp{iElm}*rho*Kp{iElm}';
        [rhoPrime,realEpsilon] = perturbation_channel(rho3);
        epsilon_max = max(epsilon_max,realEpsilon);
        fval = fval + real(trace(rhoPrime* logm(rhoPrime)));
    end

   fval = real(fval)/log(2);
end

