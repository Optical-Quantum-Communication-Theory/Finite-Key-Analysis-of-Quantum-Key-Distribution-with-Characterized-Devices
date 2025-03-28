%% 
% outputs statistics.
function prob_dist = opticalebBB84Sim(protocol_parameters, exp_parameters, mu, fine_grained) 
   LA =  exp_parameters.LA; % Alice - source distance
   etadA = exp_parameters.etadA; % Alice's detector efficiency
   LB =  exp_parameters.LB;% Bob - source distance
   etadB = exp_parameters.etadB; % Bob's detector efficiency
   ed = exp_parameters.ed; % intrinsic detector error
   lambda = mu/2; % mean photon number of each half signal
   etaA = 10^(-0.02*LA)*etadA; % Alice's total transmittance 
   etaB = 10^(-0.02*LB)*etadB; % Bob's total transmittance
   pdA = exp_parameters.pdA; % Alice's detector dark count
   pdB = exp_parameters.pdB; % Bob's detector dark count
   
   pAz = protocol_parameters.pAz;
   pBz = protocol_parameters.pBz;
   pAx = 1 - pAz;
   pBx = 1 - pBz;
   
   paramA = (1-pdA)*(1-pdB)/(1+(etaA+etaB-etaA*etaB)*lambda)^2;
   paramB =  (1-pdB)/(1+etaB*lambda)^2-(1-pdA)*(1-pdB)/(1+(etaA+etaB-etaA*etaB)*lambda)^2;
   paramC =  (1-pdA)/(1+etaA*lambda)^2-(1-pdA)*(1-pdB)/(1+(etaA+etaB-etaA*etaB)*lambda)^2;
   paramQ = 1-  (1-pdA)/(1+etaA*lambda)^2 - (1-pdB)/(1+etaB*lambda)^2 + (1-pdA)*(1-pdB)/(1+(etaA+etaB-etaA*etaB)*lambda)^2;
   paramEQ = paramQ/2 - 2*(1/2-ed)*etaA*etaB*lambda*(1+lambda)/((1+etaA*lambda)*(1+etaB*lambda)*(1+etaA*lambda+etaB*lambda-etaA*etaB*lambda));
   if fine_grained
    prob_dist = [paramEQ*pAz*pBz/2,(paramQ-paramEQ)*pAz*pBz/2,paramQ*pAz*pBx/4, paramQ*pAz*pBx/4, paramB*pAz/2,;...
     	(paramQ-paramEQ)*pAz*pBz/2, paramEQ*pAz*pBz/2,paramQ*pAz*pBx/4, paramQ*pAz*pBx/4,paramB*pAz/2;...
      paramQ*pAx*pBz/4, paramQ*pAx*pBz/4,paramEQ*pAx*pBx/2,(paramQ-paramEQ)*pAx*pBx/2,paramB*pAx/2;...
        paramQ*pAx*pBz/4, paramQ*pAx*pBz/4,(paramQ-paramEQ)*pAx*pBx/2,paramEQ*pAx*pBx/2,paramB*pAx/2;...
        paramC*pBz/2,paramC*pBz/2,paramC*pBx/2,paramC*pBx/2,paramA];
  
   else
      % error rate, 1- error rate, no detection
      prob_dist = [paramEQ,(1-paramEQ/paramQ)*paramQ,1-paramQ]; 
   end
 
        
        
   
    
   
end