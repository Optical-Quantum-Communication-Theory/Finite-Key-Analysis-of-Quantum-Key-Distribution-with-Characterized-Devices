function  Y  = xcorrop(d)
    Y = 0;
    Fop = fourier_operator(d);
    for j = 1:d

        Y = Y + vpa(kron(Fop*zProjector(d,j)*Fop' , conj(Fop)*zProjector(d,j)*(Fop.')));
   
    end
    Y  = (Y+Y')/2; % make sure it is Hermitian
end