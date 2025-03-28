%% Correction terms for QKD via EAT
% This code allows several options for EAT
% It allows the possibility to use EAT version 2 - [Dupuis and Fawzi (2019)], version 3 [Dupuis (2021)].
% Note that there is no error correction leakage term,
% which is a separate thing from correction terms calculated here.
% All terms are to be subtracted from the leading order term. (All signs are flipped.)

%% Output:
% FOT: first-order term, which needs to be multipled by N, signal size;
% SOT: second-order term, which needs to be multipled by sqrt(N);
% CT: constant term
% Out: a struct that stores extra information if needed, e.g. for debugging
% purpose.
function [FOT,SOT,CT,Out] = calculate_correction(protocol_parameters, variable_parameters,options)
    if ~isfield(options, 'EATversion')
       options.version = 'v3';
    end


    Out = struct;
    dA = protocol_parameters.dimAbar;
    dB = protocol_parameters.dimBbar;
    dS = protocol_parameters.dimS;

    sizeA = protocol_parameters.sizeA;
    sizeB = protocol_parameters.sizeB;

    epsilon_bar = protocol_parameters.epsilon_bar;
    epsilon_EC = protocol_parameters.epsilon_EC;
    epsilon_PA = protocol_parameters.epsilon_PA;

    test_prob = variable_parameters.test_prob;

    correction1  = test_prob*log2(sizeA*sizeB);

    FOT = correction1;
	  SOT = 0;
    CT = 0;

 	  maxf= variable_parameters.maxf;
	  minf_sigma = variable_parameters.minf_sigma;
 	  varf = double(variable_parameters.varf);
    signal_size = variable_parameters.signal_size;

    K_alpha_fun = @(alpha) 1/(6*(2-alpha)^3*log(2))*2^((alpha-1)*(log2(dS)+maxf-minf_sigma))*(log(2^(log2(dS)+maxf-minf_sigma)+exp(2)))^3;

    % note that the K_alpha function already uses the fact that S is classical
    V = sqrt(varf+2)+log2(2*dS^2+1);

    optimoptions = optimset('TolX',1e-14,'MaxIter',1000,'MaxFunEvals',1000);

    if strcmp(options.EATversion,'v2')
        % Theorem 15
        epsilon_acc = epsilon_bar + epsilon_PA + epsilon_EC;

        total_alpha_fun = @(alpha)  (alpha-1)*log(2)/2*V^2 + 1/(alpha-1)*double(log2(32/(epsilon_bar^2*epsilon_acc^2)))/signal_size +  (alpha-1)^2 * K_alpha_fun(alpha);

        % use fminbnd to optimize alpha
        [op_alpha,correction2] = fminbnd(total_alpha_fun,1,2,optimoptions);
        FOT = FOT+correction2;

        CT = 2*log2(1-sqrt(1-epsilon_bar^2/16))+2*double(log2(2/epsilon_PA)); % constant term

        SOT = 2*log2(1+2*dS)*sqrt(1-2*log2(epsilon_bar/4*epsilon_acc)); %second-order term, i.e., \sqrt{N} term

        Out.op_alpha= op_alpha;


    elseif strcmp(options.EATversion,'v3')
        % Theorem 5
        epsilon_sec = epsilon_bar + epsilon_PA;
        epsilon_acc = epsilon_sec + epsilon_EC;

        total_beta_fun = @(beta)  (beta-1)*log(2)/2*V^2 + beta/(beta-1)*double(log2(1/epsilon_acc))/signal_size +  (beta-1)^2 * K_alpha_fun(beta);

        total_delta_fun = @(delta) -delta/(1-delta)*double(log2(epsilon_acc))/signal_size;

        total_alpha_fun = @(alpha) (alpha)/(alpha-1)*double(log2(1/epsilon_sec))/signal_size;

        total_func = @(x) total_alpha_fun((-x(1)+x(2))/(-1+2*x(2)-x(1)*x(2))) + total_beta_fun(x(1))+total_delta_fun(x(2));

        % use fmincon to optimize alpha, beta, delta
        [x,correction2]=fmincon(total_func, [1.000000001,0.5001],[],[],[],[],[1+eps,1/2+eps],[2-eps,1-eps], [], optimoptions);


        FOT = FOT + correction2;
        CT = -1;
        SOT = 0;
        Out.x = x;

    end
end
