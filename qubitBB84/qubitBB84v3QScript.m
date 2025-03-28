%% protocol: qubit BB84
% For key generation, we use ZZ basis.
% number of signals n = 10^5 - 10^14;
% noise = 0.01;
% fix epsilons 0.25e-8 for each one: epsilon_EC, epsilon_PA, epsilon_EA,
% and bar{\epsilon}.
% calculate the key rate that includes all correction terms
%%


% before running this script,
% remember to add CVX and EATSolver, subroutines folders to the PATH

filename = 'bb84fig2curve2.mat';

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


n_range = 10.^[5:0.5:14];

errorRate = 0.01;

acceptance_t = 0.005;
fEC = 1.16; % Error correction inefficiency
EAT_option = struct;
EAT_option.EATversion = 'v3'; % possible options: v2, v3;
EAT_option.crossover = true;
EAT_option.ZbasisOnly = true;

% list of free parameters
pZ_range =  1-sqrt(10.^( -[2:0.1:4]));
pZ_range2 =  1-sqrt(10.^( -[3:0.2:7]));



tic
phiPlus = 1/sqrt(2) * [1;0;0;1];
phiMinus = 1/sqrt(2)*[1;0;0;-1];
psiPlus = 1/sqrt(2)*[0;1;1;0];
psiMinus = 1/sqrt(2)*[0;1;-1;0];


% We depolarize the state by just writing down the output state
rho_simul = (1-3/2*errorRate)*(phiPlus*phiPlus') + (errorRate/2)*(phiMinus*phiMinus' + psiMinus*psiMinus' +psiPlus*psiPlus');


maxKeyList = zeros(1,length(n_range));
minTFcell = cell(1,length(n_range));

exitList=  cell(1,length(n_range));
upperList = cell(1,length(n_range));
keyList =cell(1,length(n_range));

parfor i = 1:length(n_range)

   tmpkeyList =zeros(length(pZ_range),1);
   tmpupList  = zeros(length(pZ_range),1);
   tmpexitList = zeros(length(pZ_range),1);
   tmpminCell = cell(length(pZ_range),1);
   n = n_range(i);

    for j = 1:length(pZ_range)
        % choose the pZ (and thus testing probability) differently for
        % different block sizes
        % heuristic choices
        if n<1e10
            pZ = pZ_range(j);
        else
            pZ = pZ_range2(j);
       	end


            % testing part: XX block
            pT = (1-pZ)^2;

            % obtain the protocol setup
            protocol = protocol_setup_BB84(pZ,pT, EAT_option);

            % generate problem information
            problem = struct;

            problem.p_accept = eval_operators(rho_simul,protocol.POVM_pe);
            problem.rho0 = rho_simul;
            problem.c0 = double(2*sqrt(log(2))*sqrt( 1 - 2*log2(epsilons.bar*epsilons.EA) )/sqrt(n));
            problem.c1 = 1/sqrt(2*pT);
            % solve primal and linerarised dual problems

            try

                %exitflag = 1;
                [rho1,fval_primal,~] = EAT_step1( protocol,problem);
                [minTF,fval_dual] = EAT_step2( protocol,problem,rho1 );

                termH = evaluate_min_tradeoff(minTF,problem.p_accept.'/sum(problem.p_accept), acceptance_t,0);

                ECcost=  fEC * qubitBB84calculateECcost(rho_simul, protocol.AlicePOVM, protocol.BobPOVM, EAT_option);
                ECcost = ECcost + double(log2(2/epsilons.EC)/n);
                [keyrate,  Out] = evaluateKeyRate(pT,  n, minTF, termH, ECcost, protocol_parameters, EAT_option);


            catch ME
                keyrate = 0;
                fval_primal = -inf;
                fval_dual = inf;
                exitflag = 0;
                minTF = 0;
                fprintf("error occurs")
            end
            % post analysis
            gap = fval_primal - fval_dual;


            % Save all the data
            tmpupList(j) = fval_primal;
            tmpkeyList(j) = keyrate;
            tmpminCell{j} = minTF;
            %tmpexitList(j) = exitflag;

    end
    minTFcell{i} = tmpminCell;
    exitList{i} =  tmpexitList;
    upperList{i} = tmpupList;
    keyList{i}  = tmpkeyList;
    maxKeyList(i) = max(max(tmpkeyList),0);
end
toc
save(filename);
