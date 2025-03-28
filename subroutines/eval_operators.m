function gamma = eval_operators(rho,Gamma)
    gamma = zeros(1,length(Gamma));
    for i = 1:length(Gamma)
        gamma(i) = real(trace(Gamma{i}*rho));
    end
end