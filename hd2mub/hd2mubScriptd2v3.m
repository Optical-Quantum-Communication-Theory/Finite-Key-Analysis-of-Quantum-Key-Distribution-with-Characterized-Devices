%% protocol: high-dimensional 2-mub
% noise = 0.01; t = 0, dim = 2
% calculate the key rate that includes all correction terms
%%


% before running this script,
% remember to add CVX and EATSolver, subroutines folders to the PATH

filename = 'hd2dv3.mat';


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

% list of fixed parameters
n_range = 10.^[4.5:0.5:16];

errorRate = 0.01;
dim = 2; % prime number only due to data simulation

acceptance_t = 0;
fEC = 1.16;
EAT_option = struct;
EAT_option.EATversion = 'v3'; % possible options: v2, v3;
EAT_option.crossover = true;
EAT_option.ZbasisOnly = true;
% list of free parameters
pZ_range = 1-sqrt(10.^( -[2:0.1:4]));
pZ_range2 = 1-sqrt(10.^( -[3:0.2:7]));


tic

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
maxKeyList = zeros(1,length(n_range));
minTFcell = cell(1,length(n_range));

exitList=  cell(1,length(n_range));
upperList = cell(1,length(n_range));
keyList =cell(1,length(n_range));
parfor i = 1:length(n_range)

	n = n_range(i);
	tmpkeyList =zeros(length(pZ_range),1);
	tmpupList  = zeros(length(pZ_range),1);
	tmpexitList = zeros(length(pZ_range),1);
	tmpminCell = cell(length(pZ_range),1);

	for j = 1:length(pZ_range)
    	if n<1e10
            pZ = pZ_range(j);
      else
            pZ = pZ_range2(j);
      end

      pT = (1-pZ)^2;

      % generate protocol
      protocol = gen_protocol_hd2mub(dim, pZ,pT, EAT_option );

      % generate problem information
      problem = struct;
      problem.p_accept = eval_operators( rho_simul,protocol.POVM_pe);
      problem.rho0 = rho_simul;

      problem.c0 = double(2*sqrt(log(2))*sqrt( 1 - 2*log2(epsilons.bar*epsilons.EA) )/sqrt(n));
      problem.c1 = 1/sqrt(2*pT);

      try
        	[rho1,fval_primal,gap1] = EAT_step1( protocol,problem);
          [minTF,fval_dual] = EAT_step2( protocol,problem,rho1 );

          termH = evaluate_min_tradeoff(minTF,problem.p_accept.'/sum(problem.p_accept), acceptance_t,0);

          ECcost=  fEC * hd2mubCalculateECcost(dim,rho_simul, protocol.AlicePOVM, protocol.BobPOVM, EAT_option);
          ECcost = ECcost + double(log2(2/epsilons.EC)/n);
          [keyrate,  Out] = evaluateKeyRate(pT,  n, minTF, termH, ECcost, protocol_parameters, EAT_option);
          exitflag = 1;

     	catch ME
          keyrate = 0;
         	fval_primal = -inf;
          fval_dual = inf;
          exitflag = 0;
        	minTF = zeros(1,length(problem.p_accept));
        end

      gap = fval_primal - fval_dual;
     	tmpupList(j) = fval_primal;
       tmpkeyList(j) = keyrate;
    	tmpminCell{j} = minTF;
      tmpexitList(j) = exitflag;

    end
    minTFcell{i} = tmpminCell;
    exitList{i} =  tmpexitList;
    upperList{i} = tmpupList;
    keyList{i}  = tmpkeyList;
    maxKeyList(i) = max(max(tmpkeyList));
end
toc
save(filename);
