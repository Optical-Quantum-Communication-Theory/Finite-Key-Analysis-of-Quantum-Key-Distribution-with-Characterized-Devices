%% FUNCTION NAME:    find_min_tradeoff_coeff
% Date of last modification: August 28, 2023

% Similar to previous step 2 solver, except no constraint violation is considered so that no need to duplicate gamma

function [dualY, zetaEp] = find_min_tradeoff_coeff(rho, protocol_setup, problem_setup)
    
    warning('off','MATLAB:logm:nonPosRealEig');

    Kpk = protocol_setup.Kpk;
    Kp =  protocol_setup.Kp;
    Gamma_eq = problem_setup.eq_observables;
    gamma_eq = problem_setup.eq_expectations;
    [~, epsilon1] = objfun_pert(rho, Kp, Kpk);
 
    [gradf, epsilon2] = gradf_pert(rho,Kp, Kpk); 
    gradf = (gradf + gradf')/2;
    dprime = size(rho, 1);
    
    epsilon = max(epsilon1, epsilon2);
 
    if epsilon> 1/(exp(1)*(dprime-1))
      ME = MException('dualSDP:epsilon too large','Theorem cannot be applied. Please have a better rho to start with');
      throw(ME);
    end
  
    % remove linear dependence if necessary 
    %[Gamma_eq, independentCols] = removeLinearDependence(Gamma_eq);
    %Gamma_eq = Gamma_eq(independentCols);
    
    dualY = submaxproblem(Gamma_eq,gamma_eq, gradf);

    if epsilon == 0
        zetaEp = 0;
    else 
        nTerms = 0;
        for i= 1:numel(Kp)
           nTerms = nTerms + numel(Kpk{i}); 
        end
        % this is the correction term due to epsilon perturbation of rho.
       	%zetaEp = (numel(Kpk)+1)*(numel(Kp))* epsilon * (dprime-1) * log2(dprime/(epsilon*(dprime-1)));
        zetaEp = nTerms* epsilon * (dprime-1) * log2(dprime/(epsilon*(dprime-1)));
 
    end
   
 
end


function dualY = submaxproblem(Gamma_eq, gamma_eq,  gradf)
    nEqConstraints = length(gamma_eq);
    totalDim = size(gradf, 1);
   
    cvx_begin sdp quiet
        cvx_precision best
        variable dualY(nEqConstraints) 
        maximize  sum(gamma_eq  .* dualY)
        gradf - sdpCondition(dualY, Gamma_eq)  == hermitian_semidefinite(totalDim)  
 
    cvx_end
    
%     if ~strcmp(cvx_status, 'Solved')
%         try
%             checkValue = lambda_min(gradf -  sdpCondition(dualY, Gamma_eq));
%         catch 
%            checkValue = inf; 
%         end
%         flag = strcat(cvx_status, num2str(checkValue));
%     else 
%         flag=lambda_min(gradf -  sdpCondition(dualY, Gamma_eq));
%     end
end


function result = sdpCondition(dualY, GammaVector)
   	result =0;
    for iConstraint = 1 : length(dualY)
        result = result + dualY(iConstraint) * GammaVector{iConstraint};
    end
        
end