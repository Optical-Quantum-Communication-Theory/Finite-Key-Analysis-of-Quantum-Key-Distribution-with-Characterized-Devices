function [krausOp, keyMap, observables, expectations] = hd2mubDataGen(dim, pZ, errorRate)
    pAz = pZ;
    pAx = 1-pZ;
    pBz = pZ;
    pBx = 1-pBz;

    unscaledAlicePOVM = gethdONBProjectors(dim,true);
    unscaledBobPOVM = gethdONBProjectors(dim,false);
 
    scaledAlicePOVM = cell(dim*2,1);
    scaledBobPOVM = cell(dim*2,1);
    
   
    for i =1:dim
       
        scaledAlicePOVM{i} = pAz*unscaledAlicePOVM{i};
     	scaledBobPOVM{i} = pBz*unscaledBobPOVM{i};
        
       
        scaledAlicePOVM{i+dim} = pAx*unscaledAlicePOVM{i+dim};
        scaledBobPOVM{i+dim} = pBx*unscaledBobPOVM{i+dim};
       
    end
    
    
    krausOp = {};
    
    keyMap = {};
    for i =1:dim
        krausOp =[krausOp, kron(kron(zket(dim,i),sqrtm(scaledAlicePOVM{i})),sqrt(pBz)*eye(dim))];
        keyMap = [keyMap,kron(zProjector(dim,i),eye(dim^2))];
    end
    
    % coarse-grained POVM;
  	phaseErrorOp = double(eye(dim^2)- xcorrop(dim));    
 	observables = {phaseErrorOp,eye(dim^2)-phaseErrorOp};
    
    % fine-grained POVM 
%     observables = {};
%     for i = 1:dim^2
%         for j = 1:dim^2
%             observables = [observables,kron(scaledAlicePOVM{i},scaledBobPOVM{j})];
%         end
%     end
    
    %
    
    
    
    % simulation
    basis = gethdBellBasis(dim);
    Q = errorRate;
    qjk = Q/(dim-1)*ones(1,dim);
    qjk(1) = 1-Q;
    lambda = zeros(dim,dim);
    dpRho = 0;
    for j = 1:dim
        for k = 1:dim
            tmp = 0;
            for s = 1:dim
                tmp = tmp + qjk(mod((s-1)*(j-1)-(k-1),dim)+1);
            end
            lambda(j,k) = (tmp + qjk(j)-1)/dim;
            dpRho = dpRho + lambda(j,k)*(basis{j,k}*basis{j,k}');
        end
    end
    rho_simul = (dpRho + dpRho')/2;
    expectations = zeros(length(observables),1);
    % calculate expectations
    for i = 1:length(observables)
        expectations(i) = trace(rho_simul * observables{i});
    end
    
end