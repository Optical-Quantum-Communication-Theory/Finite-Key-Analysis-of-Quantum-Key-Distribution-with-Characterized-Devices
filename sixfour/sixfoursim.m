%% Generate POVM operators as well as simulated statistics.

function outputRho = sixfoursim(errorRate, rotation_axis, theta) 
    phiPlus = 1/sqrt(2) * [1;0;0;1];
    phiMinus = 1/sqrt(2)*[1;0;0;-1];
    psiPlus = 1/sqrt(2)*[0;1;1;0];
    psiMinus = 1/sqrt(2)*[0;1;-1;0];
    
    
    %We depolarize the state by just writing down the output state
    dpRho = (1-3/2*errorRate)*(phiPlus*phiPlus') + (errorRate/2)*(phiMinus*phiMinus' + psiMinus*psiMinus' +psiPlus*psiPlus');
    
    
    sigmaZ = [1 0 ; 0 -1];
    sigmaX = [0 1 ; 1 0];
    sigmaY = [0 -1i;1i 0];
    %Now we rotate the state to get the output state
    switch lower(rotation_axis)
        case 'z'
            rot = expm(1i*theta*sigmaZ);
        case 'x'
            rot = expm(1i*theta*sigmaX);
        case 'y'
            rot = expm(1i*theta*sigmaY);
        otherwise
            rot = eye(2);
    end
    
    outputRho = kron(eye(2),rot)*dpRho* (kron(eye(2),rot)');
    outputRho = (outputRho + outputRho')/2;
    
  


   
    
    
   
end