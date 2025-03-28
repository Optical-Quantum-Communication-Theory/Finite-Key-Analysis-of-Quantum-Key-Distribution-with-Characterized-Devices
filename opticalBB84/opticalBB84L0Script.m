%% protocol: optical BB84
% For key generation, we use ZZ basis.
% number of signals n = 10^7 - 10^15;
% noise = 0.01;
% fix epsilons 0.25e-8 for each one: epsilon_EC, epsilon_PA, epsilon_EA,
% and bar{\epsilon}.
% calculate the key rate that includes all correction terms

%%

% before running this script,
% remember to add CVX and EATSolver, subroutines folders to the PATH

filename = 'opBB84v3L0t.mat';



    epsilons = {};
    epsilons.EC = vpa(0.25e-8);
    epsilons.bar = vpa(0.25e-8);
    epsilons.PA = vpa(0.25e-8);
    epsilons.EA = vpa(0.25e-8);
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

    exp_parameters = struct;
    total_distance = 0.0;
    exp_parameters.LA = total_distance/2; % Alice - source distance
    exp_parameters.etadA=0.8; % Alice's detector efficiency
    exp_parameters.LB = total_distance/2;% Bob - source distance
    exp_parameters.etadB =0.8; % Bob's detector efficiency
    exp_parameters.ed = 0.01; % intrinsic detector error
    %exp_parameters.mu = 0.3; % mean photon number of each half signal
    exp_parameters.pdA= 1e-7;
    exp_parameters.pdB = 1e-7;

    % list of fixed parameters
    n_range = 10.^[7:0.5:15];
    mu_range = 0.2:0.02:0.4;


    dim = 3; % Alice's dimension = Bob's dimension = qubit + vacuum
    acceptance_t = 0.005;
    fEC = 1.16;
    EAT_option = struct;
    EAT_option.EATversion = 'v3'; % possible options: v2, v3;
  	EAT_option.crossover = true;
    EAT_option.ZbasisOnly = true;
    % list of free parameters

    pZ_range = 1-sqrt(10.^( -[3:0.2:7]));




tic


maxKeyList = zeros(1,length(n_range));
minTFcell = cell(1,length(n_range));

exitList=  cell(1,length(n_range));
upperList = cell(1,length(n_range));
keyList = cell(1,length(n_range));
parfor i = 1:length(n_range)


   n = n_range(i);
   tmpkeyList =zeros(length(pZ_range),length(mu_range));
   tmpupList  = zeros(length(pZ_range),length(mu_range));
   tmpexitList = zeros(length(pZ_range),length(mu_range));
   tmpminCell = cell(length(pZ_range),length(mu_range));

  for j = 1:length(pZ_range)
        pZ = pZ_range(j);

      	protocol_setup = struct;
      	protocol_setup.pAz = pZ;
     		protocol_setup.pBz = pZ;
        for k = 1:length(mu_range)
            mu = mu_range(k);
            prob_dist = opticalebBB84Sim(protocol_setup, exp_parameters, mu, false);
            if prob_dist(1)/(1- prob_dist(3)) > 0.11
                 continue
        		end
        		key_generation_freqs = opticalebBB84Sim(protocol_setup, exp_parameters, mu, true);
        		key_generation_freqs= key_generation_freqs(1:2,1:2);



        		pT = (1-pZ)^2;

        		% generate protocol
        		protocol = gen_protocol_opticalBB84(dim, pZ,pT, EAT_option);

            % generate problem information
            problem = struct;

            problem.p_accept = prob_dist;

            problem.rho0 = eye(dim^2);
            problem.c0 = double(2*sqrt(log(2))*sqrt( 1 - 2*log2(epsilons.bar*epsilons.EA) )/sqrt(n));
            problem.c1 = 1/sqrt(2*pT);

            try

            	[rho1,fval_primal,gap1] = EAT_step1( protocol,problem);
            	exitflag = 1;
            	[minTF,fval_dual] = EAT_step2( protocol,problem,rho1 );

            	termH = evaluate_min_tradeoff(minTF,problem.p_accept.'/sum(problem.p_accept), acceptance_t,0);


            	ECcost=  fEC * opticalBB84CalculateECcost(key_generation_freqs, EAT_option);
            	ECcost = ECcost + double(log2(2/epsilons.EC)/n);
            	[keyrate,  Out] = evaluateKeyRate(pT,  n, minTF, termH, ECcost, protocol_parameters, EAT_option);



          	catch ME
                keyrate = 0;
                fval_primal = -inf;
                fval_dual = inf;
                exitflag = 0;
                minTF = zeros(1,length(problem.p_accept));
            end
            % post analysis
            gap = fval_primal - fval_dual;


            % Save EVERYTHING
            tmpupList(j,k) = fval_primal;
            tmpkeyList(j,k) = keyrate;
            tmpminCell{j,k} = minTF;
            tmpexitList(j,k) = exitflag;
        end
    end
    minTFcell{i} = tmpminCell;
    exitList{i} =  tmpexitList;
    upperList{i} = tmpupList;
    keyList{i}  = tmpkeyList;
    maxKeyList(i) = max(max(max(tmpkeyList)),0);
end
toc
save(filename);
