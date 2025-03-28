% This script is to run six-four protocol with coarse-grained constraints

% before running this script,
% remember to add CVX and EATSolver, subroutines folders to the PATH

filename = 'sixfouralg1.mat';




    % experimental parameters
    fEC = 1.16; % error correction efficiency
    errorRate = 0.01; % simulation

    % choice of protocol security parameters

    epsilons = {};
    epsilons.EC = vpa(2.5e-9);
    epsilons.bar = vpa(2.5e-9);
    epsilons.PA = vpa(2.5e-9);
    epsilons.EA = vpa(2.5e-9);

    % define the size of t-ball for acceptance set Q.
    acceptance_t = 0.005;

    EAT_option = struct;

     % choices for  for EAT variations
    EAT_option.EATversion = 'v3'; % possible options:  v2, v3;
    EAT_option.crossover = true; % can decide whether to use crossover min-tradeoff function
    EAT_option.ZbasisOnly = true;

    % solver options
    solver_option = struct;
    solver_option.verbose = 1;
    solver_option.linearconstrainttolerance = 1e-10;
    solver_option.maxgap = 1e-8;
    solver_option.initmethod = 2;
    solver_option.maxiter = 200;


    protocol_parameters = struct;
    protocol_parameters.sizeA = 2;
    protocol_parameters.sizeB = 2;
    protocol_parameters.dimS = 9;
    protocol_parameters.dimAbar = 2;
    protocol_parameters.dimBbar = 2;
    protocol_parameters.epsilon_EC = epsilons.EC;
    protocol_parameters.epsilon_EA = epsilons.EA;
    protocol_parameters.epsilon_bar = epsilons.bar;
    protocol_parameters.epsilon_PA = epsilons.PA;


    simulation_function = @sixfourSim_alg1;


	  preferred_error_rate = 0.01; % choose a preferred parameters to define the acceptance_set

    search_parameter_range = 0.005:0.005:0.07;


    expRange = 5:0.5:14;

    n_range = round(10.^expRange);

    nFunctions = numel(search_parameter_range);
    pZ_range =1-sqrt(10.^( -[2:0.1:4]));
    pZ_range2 = 1-sqrt(10.^( -[3:0.2:7]));

    nPzRange = length(pZ_range);
    nBlockSize = length(n_range);


    ECList = zeros(nBlockSize,nPzRange);
    termHList = cell(nBlockSize,nPzRange);
    fullKeyList = cell(nBlockSize,nPzRange);
    all_correction_term_list= cell(nBlockSize,nPzRange);
    all_h_min_tradeoff_list = cell(nBlockSize,nPzRange);
    maxKeyList = zeros(1,nBlockSize);

    for i = 1: nBlockSize

        signal_size =  n_range(i);

        keyList = zeros(1,nPzRange);
        for j = 1:nPzRange
            all_correction_term_list{i,j} = cell(1, nFunctions);
            termHList{i,j} = zeros(1, nFunctions);
            fullKeyList{i,j} = zeros(1, nFunctions);

            if signal_size < 1e10
                pZ = pZ_range(j);
            else
                pZ = pZ_range2(j);
            end
            pAList = [pZ, (1-pZ)/2, (1-pZ)/2];
            pBList = [pZ, 1-pZ];
            test_prob = (1-pZ)^2;
            protocol_setup = gen_protocol_sixfour( pAList, pBList, test_prob, EAT_option );

            POVM_pe = protocol_setup.POVM_pe;



            AlicePOVM = protocol_setup.AlicePOVM;
            BobPOVM = protocol_setup.BobPOVM;

            other_parameters = struct;
            other_parameters.rotation_angle = 11/180*pi;
            other_parameters.rotation_axis = 'z';
            other_parameters.table = false;

            [preferred_freqs, simulated_state] = simulation_function(POVM_pe, preferred_error_rate, other_parameters);


            [h_min_tradeoff_list,zeta_list] = optimize_min_tradeoff(protocol_setup, simulation_function,  search_parameter_range, other_parameters, solver_option);
            all_h_min_tradeoff_list{i,j}= h_min_tradeoff_list;


            % generate statistics for error correction purpose

            ECcost =  fEC * sixfourcalculateECcost(simulated_state, AlicePOVM,BobPOVM, EAT_option);
            ECcost = ECcost + double(log2(2/epsilons.EC)/signal_size);
            ECList(i,j)= ECcost;



            for k = 1:nFunctions
                coeff = h_min_tradeoff_list{k};
                zeta = zeta_list(k);
                [termH, ~] = evaluate_min_tradeoff(coeff,preferred_freqs, acceptance_t, zeta);
                termHList{i,j}(k) = termH;

                [keyrate,  Out] = evaluateKeyRate(test_prob,  signal_size, coeff, termH, ECcost, protocol_parameters, EAT_option);
                fullKeyList{i,j}(k) = keyrate;
                all_correction_term_list{i,j}{k} = Out.correction_term;

            end

         keyList(j) = max(max(fullKeyList{i,j}),0);
        end

        [maxKeyList(i),index2] = max(keyList);
    end

save(filename);
