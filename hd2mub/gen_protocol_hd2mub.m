function protocol = gen_protocol_hd2mub( dim, pZ,pT, options )
    pAz = pZ;
    pAx = 1-pZ;
    pBz = pZ;
    pBx = 1-pBz;

    unscaledAlicePOVM = gethdONBProjectors(dim,true);
    unscaledBobPOVM = gethdONBProjectors(dim,false);
 
    scaledAlicePOVM = cell(dim*2,1);
    scaledBobPOVM = cell(dim*2,1);
    
   
    for i =1:dim
       
        scaledAlicePOVM{i} = pAz*unscaledAlicePOVM{i};
     	scaledBobPOVM{i} = pBz*unscaledBobPOVM{i};
       
        scaledAlicePOVM{i+dim} = (1-pAz)*unscaledAlicePOVM{i+dim};
        scaledBobPOVM{i+dim} = (1-pBz)*unscaledBobPOVM{i+dim};
       
    end
    
    Kpk = {};
    Kp = {};

    % cross-over min-tradeoff only
 
	
    phaseErrorOp = double(eye(dim^2)- xcorrop(dim)); 
    POVM_pe = {phaseErrorOp,eye(dim^2)-phaseErrorOp};
       
	% XX test
    Kp{end+1}  =       sqrt( pAx*pBx )* eye(dim^2);
    Kpk{end+1} = pref( sqrt( pAx*pBx ), cellsqrtkron( unscaledAlicePOVM(1+dim:2*dim),unscaledBobPOVM(1+dim:2*dim) ));  
     	   
	% ZZ generation
	Kp{end+1}  =       sqrt( pAz*pBz )* eye(dim^2);
	Kpk{end+1} = pref( sqrt( pAz*pBz ), cellsqrtkron( unscaledAlicePOVM(1:dim),{eye(dim)} ));       
   

    protocol = struct;    
    protocol.Kpk = Kpk;
    protocol.Kp = Kp;
    protocol.POVM_pe  = POVM_pe;
    protocol.dim = dim^2;
    protocol.pT = pT;
    protocol.AlicePOVM = scaledAlicePOVM;
    protocol.BobPOVM = scaledBobPOVM;
    
end
