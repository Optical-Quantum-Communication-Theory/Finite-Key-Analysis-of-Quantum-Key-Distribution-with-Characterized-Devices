%% This function generates the Bell state in dimension d: \ket{U_{j,k}} = 1/sqrt{d} \sum_s \omega^{sk} \ket{s}\ket{s+j}.
% d: dimension of the Hilbert space of a single system
% r, s: parameters for the Bell states
function Ujk = hdBellState(d,j,k)
    Id= eye(d);
    omega = exp(2i*pi/d);
    Ujk = 0;
    for s =1:d
        Ujk = Ujk + omega^((s-1)*k)*kron(Id(:,s),Id(:,mod(s-1+j,d)+1));
    end
    Ujk = 1/sqrt(d)*Ujk;
end