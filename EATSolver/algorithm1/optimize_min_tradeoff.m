%% Name: optimize_min_tradeoff.m
% Date of last modification: August 31, 2023
% Description: optimize the choice of min-tradeoff function heuristically.
% It assumes one-dim search.
% It returns a cell of min-tradeoff function coefficients as well as a list of
% zeta correction term for the perturbation.
%% Explanations 
% search_parameter_range is a list of values for the searchable parameter.
% other_parameters is a struct that contains all other parameters not used
% for searching, i.e., they are needed for the construction of q* but
% different construction of q* use the same parameters in other_parameters
% and only search_parameter changes to give different q*'s.

function [min_tradeoff_list,zeta_list] = optimize_min_tradeoff(protocol_setup, simulation_function, search_parameter_range, other_parameters, solver_option)
    nPoints = length(search_parameter_range);
    min_tradeoff_list = cell(nPoints,1);
    zeta_list = zeros(nPoints,1);
    POVM_pe = protocol_setup.POVM_pe;
    for iPoint = 1:nPoints
    
        gamma_eq = simulation_function(POVM_pe,search_parameter_range(iPoint), other_parameters);
        
        problem_setup = {};
        problem_setup.eq_observables = protocol_setup.POVM_pe;
        problem_setup.eq_expectations = gamma_eq;
    
        problem_setup.rho0 = eye(protocol_setup.dim);
        
        % the usual two-step algorithm
        try
            rho = find_optimal_rho(protocol_setup, problem_setup, solver_option);
            [min_tradeoff_list{iPoint}, zeta_list(iPoint)] = find_min_tradeoff_coeff(rho,protocol_setup, problem_setup);
            if any(isnan(min_tradeoff_list{iPoint}))
                min_tradeoff_list{iPoint} = zeros(size(min_tradeoff_list{iPoint}));
            end
        catch ME
            
            min_tradeoff_list{iPoint} = zeros(length(gamma_eq),1);
            zeta_list(iPoint) = 0;
        end
    end
    
end