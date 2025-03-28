%% FUNCTION NAME: EAT_step1
% EAT algorithm 2 based on Frank-Wolfe algorithm
% 
%%
function [rho,fval,gap] = EAT_step1(protocol_setup, problem_setup,options)

    % options
    defaultOptions.maxiter = 300; % Maximum number of iterations
    defaultOptions.maxgap = 1e-6; % Maximum gap as specified by the Frank-Wolfe algorithm
    defaultOptions.linesearchprecision = 1e-20;
    defaultOptions.linesearchminstep = 1e-3;
    defaultOptions.linearconstrainttolerance = 1e-10;
    defaultOptions.initmethod = 2; % 1 for closest to rho0, 2 for maximize minimum eigenvalue
    defaultOptions.verbose = 0;
    defaultOptions.epsilon = 0; % 0<=epsilon<=1/(e(d'-1)), epsilon is related to perturbed channel
    defaultOptions.perturbation = 1e-15;
    
    if nargin == 3
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
    
   
    % access Kraus operators and the parameter estimation POVM
    Kpk = protocol_setup.Kpk;
    Kp =  protocol_setup.Kp;
    POVM_pe = protocol_setup.POVM_pe;
    
    rho0 = problem_setup.rho0;
    expectations = problem_setup.p_accept;

    c0 = problem_setup.c0;
    c1 = problem_setup.c1;
 
    lx = 1; % scaling function for error vector
    % we always set lx = 1; not testing other values

    fval = 0;
    gap = Inf;
   

    %[POVM_pe,independentCols] = removeLinearDependence(POVM_pe);
    %expectations = POVM_pe(independentCols);

    n_pe = length(POVM_pe); % # of parameter-estimation POVM elements
   
    % initial point for x
    x = ones(n_pe,1)/n_pe *(c0*c1)/lx / 100; 
   
    % project rho0 onto the set of density matrices consistent with observations
    rho = closestDensityMatrix(rho0, x, POVM_pe, expectations, options);
    rho = full(rho);
    
   
    best_upper = entropy_fun(rho,Kp,Kpk) - error_fun(x,c0,c1,lx);
    best_lower= 0;
    % Main optimization loop
    for i = 1:options.maxiter
        gradf = entropy_diff_reg(rho,Kp,Kpk);
        
        graderror = error_diff(x,c0,c1,lx);

        fval =  real(entropy_fun(rho,Kp,Kpk) - error_fun(x,c0,c1,lx));
        [deltaRho,deltaX] = subproblem(rho, x, POVM_pe, expectations, gradf, graderror,c0*c1, options);

        % perform an exact line search
        optimoptions = optimset('TolX',options.linesearchprecision);
        stepSize = fminbnd(@(t)entropy_fun(rho+t*deltaRho,Kp,Kpk) - error_fun(x+t*deltaX,c0,c1,lx),options.linesearchminstep,1,optimoptions);
        stepSize = real(stepSize);
        
        rho = rho + stepSize*deltaRho;
        rho = (rho+rho')/2;
        x = x + stepSize*deltaX;
        
        best_lower= max(best_lower, fval + real(trace(gradf*(stepSize*deltaRho)) - sum(graderror .* (stepSize*deltaX))));
        
        gap = abs(trace(gradf*(stepSize*deltaRho)) - sum(graderror .* (stepSize*deltaX)));
        if gap < options.maxgap

            break;
        end
       
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
   
    fval =  entropy_fun(rho,Kp,Kpk) - error_fun(x,c0,c1,lx);
end

function [deltaRho,deltaX] = subproblem(rho,x,observables,expectations, grad_f, grad_error, upper_limit, options)
    n = size(rho,1);
    n_pe = numel(observables);
    cvx_begin sdp quiet
        cvx_precision best
        variable deltaRho(n,n) hermitian
        variable deltaX(n_pe)
        minimize real(trace(grad_f*deltaRho) - sum(grad_error .* deltaX))
        for i = 1:n_pe
            real(trace(observables{i}*(rho + deltaRho)) - expectations(i)) <= options.linearconstrainttolerance + (x(i)+deltaX(i))
            real(trace(observables{i}*(rho + deltaRho)) - expectations(i)) >= - options.linearconstrainttolerance - (x(i)+deltaX(i))
        end
        sum(x+deltaX) <= 2*upper_limit
        trace(rho + deltaRho) == 1
        rho + deltaRho == hermitian_semidefinite(n)
    cvx_end
end

function rho = closestDensityMatrix(rho0,x, observables,expectations,options)
    % find a feasible point in the feasible set
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
            real(trace(observables{i}*rho) - expectations(i)) <= options.linearconstrainttolerance + x(i)
            real(trace(observables{i}*rho) - expectations(i)) >= - options.linearconstrainttolerance - x(i)
        end
        trace(rho) == 1
        
    cvx_end
end

function y = error_fun(x,c0,c1,lx)
    % error function (second-order)
    if sum(x)*lx > 2*(c0*c1) 
        y = +inf;
    else
        y = c0*sqrt( 1 - ( sum(x)*lx/(2*c0*c1) )^2 );
    end
end

function y = error_diff(x,c0,c1,lx)
    % derivative of the error function
    if sum(x)*lx >2*(c0*c1)
        y = +inf*ones(length(x),1);
    else                                       
        y = -ones(length(x),1) *sum(x)*c0*lx^2/(2*c0*c1)^2 /sqrt( 1 - ( sum(x)*lx/(2*c0*c1) )^2 );
    end
end