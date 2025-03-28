%% Function name: gradf_pert
%% Date of last modification: August 28, 2023

% Note: make everything in log base 2 so that we can the correct dual variables
% from linearized SDP.

function [Dfval,epsilon_max] = gradf_pert(rho, Kp, Kpk)
    nP = numel(Kp); 
    Dfval = 0;
    epsilon_max= 0;
    dim = size(rho,1);
    
    for iElm = 1: nP
        nPK = numel(Kpk{iElm});
        for jElm = 1: nPK
            rho2= Kpk{iElm}{jElm}*rho*Kpk{iElm}{jElm}';
           
            [rhoPrime,realEpsilon] = perturbation_channel(rho2);  
            epsilon_max = max(epsilon_max,realEpsilon);
            Dfval = Dfval - ((1-realEpsilon)*Kpk{iElm}{jElm}*logm(rhoPrime)*Kpk{iElm}{jElm}+realEpsilon/dim*trace(rhoPrime)*eye(dim));
        end
        
        rho3=Kp{iElm}*rho*Kp{iElm}';
       
        [rhoPrime,realEpsilon] = perturbation_channel(rho3);
        epsilon_max = max(epsilon_max,realEpsilon);
        Dfval = Dfval + ((1-realEpsilon)*Kp{iElm}*logm(rhoPrime)*Kp{iElm}+realEpsilon/dim*trace(rhoPrime)*eye(dim));
    end
    Dfval = Dfval/log(2);

end