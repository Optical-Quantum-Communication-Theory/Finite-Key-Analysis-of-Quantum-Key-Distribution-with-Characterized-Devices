function Y = zcorrop( d )

    Y = 0;
    for j=1:d

        Y = Y + kron(zProjector(d,j),zProjector(d,j));
   
    end

end