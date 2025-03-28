%% Function: evaluateKeyRate
% date of last modification: August 31, 2023
% This is a function that assembles all pieces of information for key rate
% calculation together. This is a helper function.

function [keyrate,  Out] = evaluateKeyRate(test_prob,  signal_size, coeff, termH, EC_leakage, protocol_parameters, options)
    % test_prob: the testing probability
    % signal_size: the number of signals n
    % coeff: the coeffients of crossover min-tradeoff function (note that a rescaling of (1-test_prob) is needed if using M^data map)
    % termH: the result of evaluation g(q) for q' in the acceptance set Q.
    % EC_leakage: cost of error correction (this includes f_EC)
    % protocol_parameters: a struct that stores protocol-related parameters
    %               for calculating second-order correction terms
    % options: a struct that stores EAT-related options

    if ~isfield(options, 'EATversion')
       options.version = 'v3';
    end

	  if ~isfield(options, 'crossover')
       options.crossover = false;
    end


    variable_parameters = struct;
    variable_parameters.test_prob = test_prob;
    variable_parameters.signal_size =signal_size;


    if options.crossover
        % per Eqs. (54-57) for crossover min-tradeoff functions
        maxh = max(coeff);
        minh = min(coeff);
        variable_parameters.maxf = maxh;
        variable_parameters.minf = (1-1/test_prob)*maxh + 1/test_prob*minh;
        variable_parameters.minf_sigma = minh;
        variable_parameters.varf = 1/test_prob*(maxh-minh)^2; % upper bound
    else
        % per Eq. (14) and an upper bound on the variance
        variable_parameters.maxf = max(coeff);
        variable_parameters.minf = min(coeff);
        variable_parameters.minf_sigma = min(coeff); % this is not optimal
        variable_parameters.varf = (max(coeff)-min(coeff))^2/4; % this is an upper bound
    end


    [FOT, SOT, CT] = calculate_correction(protocol_parameters, variable_parameters, options);

    correction_term = FOT*signal_size + SOT * sqrt(signal_size) + CT;
    Out.FOT = FOT;
    Out.SOT = SOT;
    Out.CT = CT;
    Out.correction_term = correction_term;
    keyrate = termH -  EC_leakage - correction_term/signal_size;

end
