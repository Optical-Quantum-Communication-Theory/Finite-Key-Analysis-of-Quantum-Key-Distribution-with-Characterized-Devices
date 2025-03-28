%% FUNCTION NAME: findOptimalRho
% new formulation, still Frank-Wolfe
% mostly unchanged. objective function and gradient were swapped.
%%
function [rho,fval,gap] = find_optimal_rho(protocol_setup, problem_setup,options)

    % options
    defaultOptions.maxiter = 300; % Maximum number of iterations
    defaultOptions.maxgap = 1e-5; % Maximum gap as specified by the Frank-Wolfe algorithm
    defaultOptions.linesearchprecision = 1e-20;
    defaultOptions.linesearchminstep = 1e-3;
    defaultOptions.linearconstrainttolerance = 1e-10;
    defaultOptions.initmethod = 2; % 1 for closest to rho0, 2 for maximize minimum eigenvalue
    defaultOptions.verbose = 0;
   
    
    if nargin == 2
        if ~isfield(options,'maxiter')
            options.maxiter = defaultOptions.maxiter;
        end
        if ~isfield(options,'maxgap')
            options.maxgap = defaultOptions.maxgap;
        end
        if ~isfield(options,'linesearchprecision')
            options.linesearchprecision = defaultOptions.linesearchprecision;
        end
        if ~isfield(options,'linesearchminstep')
            options.linesearchminstep = defaultOptions.linesearchminstep;
        end
        if ~isfield(options,'linearconstrainttolerance')
            options.linearconstrainttolerance = defaultOptions.linearconstrainttolerance;
        end
        if ~isfield(options,'initmethod')
            options.initmethod = defaultOptions.initmethod;
        end
        if ~isfield(options,'verbose')
            options.verbose = defaultOptions.verbose;
        end
        if ~isfield(options,'epsilon')
            options.epsilon = defaultOptions.epsilon;
        end
        if ~isfield(options,'perturbation')
            options.perturbation = defaultOptions.perturbation;
        end
    else
        options = defaultOptions;
    end
    
   
    rho0 = problem_setup.rho0;
    Kpk = protocol_setup.Kpk;
    Kp =  protocol_setup.Kp;
    eq_observables = problem_setup.eq_observables;
    eq_expectations = problem_setup.eq_expectations;
    
    
    fval = 0;
    gap = Inf;
    %flag = 0; % success flag

    %[observables,independentCols] = removeLinearDependence(observables);
    %expectations = expectations(independentCols);

    % project rho0 onto the set of density matrices consistent with observations
    rho = closestDensityMatrix(rho0,eq_observables,eq_expectations,options);
    rho = full(rho);

    if lambda_min(rho) < 0
        %flag = 1;
        if options.verbose
            display('Error: minimium eigenvalue less than 0.');
            display('Try increasing the constraint tolerance.');
        end
        return;
    end

 
    % Main optimization loop
    for i = 1:options.maxiter
        gradf = gradf_pert(rho,Kp,Kpk);
        deltaRho = subproblem(rho,eq_observables,eq_expectations, gradf,options);

        % perform an exact line search
        optimoptions = optimset('TolX',options.linesearchprecision);
        stepSize = fminbnd(@(t)objfun_pert(rho+t*deltaRho,Kp,Kpk),options.linesearchminstep,1,optimoptions);
        stepSize = real(stepSize);
        gap = abs(trace(rho*gradf)-trace((rho+deltaRho)*gradf));

        if gap < options.maxgap
            
            rho = rho + stepSize*deltaRho;
            rho = (rho+rho')/2;
            break;
        end
        
        rho = rho + stepSize*deltaRho;
        rho = (rho+rho')/2;
        if i == options.maxiter
            %flag = 1;
            display('Warning: Maximum iteration limit reached.');
        end
       
    end
    
    if options.verbose
        display(['Current gap: ',num2str(gap),'  Num iters: ',num2str(i)]);
    end
    
    if lambda_min(rho) < 0
        %flag = 1;
        if options.verbose
            display('Warning: minimium eigenvalue less than 0.');
        end
    end
    
    fval = objfun_pert(rho,Kp,Kpk);
end

function deltaRho = subproblem(rho,observables,expectations, gradf, options)
    n = size(rho,1);

    cvx_begin sdp quiet
        cvx_precision best
        variable deltaRho(n,n) hermitian
        minimize real(trace(gradf*deltaRho))
        for i = 1:numel(observables)
            abs(trace(observables{i}*(rho + deltaRho)) - expectations(i)) <= options.linearconstrainttolerance
        end
        rho + deltaRho == hermitian_semidefinite(n)
    cvx_end
end

function rho = closestDensityMatrix(rho0,observables,expectations,options)
    dim = size(rho0,1);
    
    cvx_begin sdp quiet
        cvx_precision best
        variable rho(dim,dim) hermitian semidefinite
        if options.initmethod == 1
            minimize norm(rho0-rho)
        elseif options.initmethod == 2
            minimize -lambda_min(rho)
        end
        for i = 1:numel(observables)
            abs(trace(observables{i}*rho) - expectations(i)) <= options.linearconstrainttolerance
        end
    cvx_end
end