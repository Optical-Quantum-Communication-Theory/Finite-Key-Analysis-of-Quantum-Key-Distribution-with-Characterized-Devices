% This script is to generate the algorithm 1 curve of Fig. 2
% Date of last modification: Dec. 19, 2024


% before running this script,
% remember to add CVX and EATSolver, subroutines folders to the PATH

filename = 'BB84fig2curve1.mat';



    % experimental parameters
    fEC = 1.16; % error correction efficiency
    errorRate = 0.01; % simulation

    % choice of protocol security parameters
    % currently randomly pick epsilon's.
    epsilons = {};
    epsilons.EC = vpa(1/3*1e-8);
    epsilons.bar = vpa(1/3*1e-8);
    epsilons.PA = vpa(1/3*1e-8);
    epsilons.EA = vpa(1);

    % define the size of t-ball for acceptance set Q.
    acceptance_t = 0.005;

    EAT_option = struct;
     % choices for  for EAT variations
    EAT_option.EATversion = 'v3'; % possible options: v2, v3;
    EAT_option.crossover = true; % can decide whether to use crossover min-tradeoff function
    EAT_option.ZbasisOnly = true;




    % solver options
    solver_option = struct;
    solver_option.verbose = 1;
    solver_option.linearconstrainttolerance = 1e-10;
    solver_option.maxgap = 1e-8;
    solver_option.initmethod = 2;
    solver_option.maxiter = 200;
    solver_option.epsilon = 0; % 0<=epsilon<=1/(e(d'-1)), epsilon is related to perturbed channel
    solver_option.perturbation = 1e-15;


    epsilons = {};
    epsilons.EC = vpa(0.25e-8);
    epsilons.bar = vpa(0.25e-8);
    epsilons.PA = vpa(0.25e-8);
    epsilons.EA = vpa(0.25e-8);


    protocol_parameters = struct;
    protocol_parameters.sizeA = 2; %
    protocol_parameters.sizeB = 2; %
    protocol_parameters.dimS = 9;
    protocol_parameters.dimAbar = 2;
    protocol_parameters.dimBbar = 2;
    protocol_parameters.epsilon_EC = epsilons.EC;
    protocol_parameters.epsilon_EA = epsilons.EA;
    protocol_parameters.epsilon_bar = epsilons.bar;
    protocol_parameters.epsilon_PA = epsilons.PA;


    other_parameters = {};

    simulation_function = @qubitBB84Sim;


	  preferred_error_rate = 0.01; % choose a preferred parameters to define the acceptance_set

    search_parameter_range = 0.005:0.005:0.07;

    % specify data range
    expRange = 5:0.5:14;
    n_range = 10.^expRange;
    nFunctions = numel(search_parameter_range);


    pZ_range =  1-sqrt(10.^( -[2:0.1:4]));
    pZ_range2 =  1-sqrt(10.^( -[3:0.2:7]));


    nBlockSize = length(n_range);

    nPzRange = length(pZ_range);
    keyList = zeros(1,nBlockSize);
    ECList = zeros(1,nBlockSize);
    termHList = cell(1,nBlockSize);
    fullKeyList = cell(1,nBlockSize);
    all_correction_term_list= cell(1,nBlockSize);
    all_h_min_tradeoff_list = cell(1,nBlockSize);




parfor i = 1: nBlockSize
    all_correction_term_list{i} = cell(nPzRange, nFunctions);
   	termHList{i} = zeros(nPzRange, nFunctions);
   	fullKeyList{i} = zeros(nPzRange, nFunctions);
   	signal_size =  n_range(i);

	all_h_min_tradeoff_list{i}=cell(nPzRange,1);

	for k = 1:nPzRange

        if signal_size <1e10
            pZ = pZ_range(k);
        else
            pZ = pZ_range2(k);
       	end

       	pT= (1-pZ)^2;


      	protocol = protocol_setup_BB84(pZ,pT, EAT_option);
       	POVM_pe = protocol.POVM_pe;
       	AlicePOVM = protocol.AlicePOVM;
     	  BobPOVM = protocol.BobPOVM;

       	[preferred_freqs, simulated_state]= simulation_function(POVM_pe, preferred_error_rate, other_parameters);

      	%fprintf(['pz = ', num2str(pz),'; ']);
      	%fprintf(['signal size = ', num2str(signal_size),'\n']);
      	%fprintf(['block size = ',num2str(signal_size),' ; ']);
      	%fprintf(['test_prob = ', num2str(test_prob), '\n']);
     	  [h_min_tradeoff_list,zeta_list] = optimize_min_tradeoff(protocol, simulation_function, search_parameter_range, other_parameters, solver_option);
      	all_h_min_tradeoff_list{i}{k}= h_min_tradeoff_list;

        % generate statistics for error correction purpose

    	   ECcost=  fEC * qubitBB84calculateECcost(simulated_state, AlicePOVM, BobPOVM, EAT_option);
         ECcost = ECcost + double(log2(2/epsilons.EC)/signal_size);


         for j = 1:nFunctions
            coeff = h_min_tradeoff_list{j};
            zeta = zeta_list(j);
         	  [termH, ~] = evaluate_min_tradeoff(coeff, preferred_freqs, acceptance_t, zeta);
            termHList{i}(k,j) = termH;

        	%fprintf(['current min-tradeoff function value = ', num2str(termH),'\n']);
         	[keyrate,  Out] =  evaluateKeyRate(pT,  signal_size, coeff, termH, ECcost, protocol_parameters, EAT_option);
         	fullKeyList{i}(k,j) = keyrate;
          all_correction_term_list{i}{k,j} = Out.correction_term;


        end
	end
  [maxtmpList,index1] = max(fullKeyList{i});
	[keyList(i),index2] =max(maxtmpList);
end


maxKeyList = max(0,keyList);


save(filename);
