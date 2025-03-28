%% Function name: fourier_operator
% This function constructs the Fourier transform operator for constructing
% high-dimensional analog of X basis.
function  Y  = fourier_operator(d)
    Y=zeros(d,d);
    for j=1:d
        for k=1:d
   
            Y(j,k) = vpa(exp(-2 * pi * 1i * (j-1)*(k-1)/d)/sqrt(d));
   
        end
    end

end