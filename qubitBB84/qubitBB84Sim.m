%% Generate simulated statistics.
% just to fit in the generic simulation method interface
function [gamma, outputRho] = qubitBB84Sim(jointPOVM, errorRate, other_parameters) 
    phiPlus = 1/sqrt(2) * [1;0;0;1];
    phiMinus = 1/sqrt(2)*[1;0;0;-1];
    psiPlus = 1/sqrt(2)*[0;1;1;0];
    psiMinus = 1/sqrt(2)*[0;1;-1;0];

    %We depolarize the state by just writing down the output state
    dpRho = (1-3/2*errorRate)*(phiPlus*phiPlus') + (errorRate/2)*(phiMinus*phiMinus' + psiMinus*psiMinus' +psiPlus*psiPlus');

	outputRho = (dpRho + dpRho')/2;
	if ~isfield(other_parameters,'table')
        other_parameters.table = false;
    end
    
    if other_parameters.table 
        % output gamma as a probability distribution table
        [nRows,mCols] = size(jointPOVM);
        gamma = zeros(nRows,mCols);
        for i =1:nRows
            for j = 1:mCols
                gamma(i,j) =  trace(dpRho*jointPOVM{i,j});
            end
        end
            
    else
        % output gamma as a vector
        nElements = numel(jointPOVM);
        gamma = zeros(nElements,1);
  
        for i=1:numel(jointPOVM)

            gamma(i) = trace(dpRho*jointPOVM{i});
        end
    end
        
        

    
   
end