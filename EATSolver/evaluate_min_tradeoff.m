%% evaluate_min_tradeoff
% sdp for optimizing over the acceptance set Q
% define the acceptance set Q as t-ball around the preferred frequency
% distribution acceptance_F.
% We need to define Q by the statistics without perp symbol.

function [lowerBound, Fvector] = evaluate_min_tradeoff(coeff, acceptance_F, acceptance_t, perturbation_correction)
    len=length(coeff);
    
    cvx_begin sdp quiet
        cvx_precision best
        variable Fvector(len) 
        minimize  sum(coeff .*Fvector)
       
        norm(acceptance_F-Fvector,1)<= acceptance_t
        Fvector>=zeros(len,1)
        Fvector<=ones(len,1)
        sum(Fvector) == 1
       
    cvx_end
  

    % note that g(q') = f(q). 
    lowerBound =  sum(coeff .*  Fvector) - perturbation_correction;
   
    
end