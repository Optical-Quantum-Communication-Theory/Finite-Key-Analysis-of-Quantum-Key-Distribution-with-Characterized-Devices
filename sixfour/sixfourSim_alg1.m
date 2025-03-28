%% Generate simulated statistics and format it for algorithm 1 solvers.

function [gamma, outputRho] = sixfourSim_alg1(jointPOVM, errorRate, other_parameters)

    rotation_axis = other_parameters.rotation_axis;
    theta = other_parameters.rotation_angle;


    phiPlus = 1/sqrt(2) * [1;0;0;1];
    phiMinus = 1/sqrt(2)*[1;0;0;-1];
    psiPlus = 1/sqrt(2)*[0;1;1;0];
    psiMinus = 1/sqrt(2)*[0;1;-1;0];


    % We depolarize the state by just writing down the output state
    dpRho = (1-3/2*errorRate)*(phiPlus*phiPlus') + (errorRate/2)*(phiMinus*phiMinus' + psiMinus*psiMinus' +psiPlus*psiPlus');


    sigmaZ = [1 0 ; 0 -1];
    sigmaX = [0 1 ; 1 0];
    sigmaY = [0 -1i;1i 0];
    % Now we rotate the state to get the output state
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


	% output gamma as a vector


    if other_parameters.table
        % output gamma as a probability distribution table
        [nRows,mCols] = size(jointPOVM);
        gamma = zeros(nRows,mCols);
        for i =1:nRows
            for j = 1:mCols
                gamma(i,j) =  real(trace(outputRho*jointPOVM{i,j}));
            end
        end


    else
        % output gamma as a vector
        nElements = numel(jointPOVM);
        gamma = zeros(nElements,1);

        for i=1:numel(jointPOVM)


            gamma(i) = real(trace(outputRho*jointPOVM{i}));
        end
    end





end
