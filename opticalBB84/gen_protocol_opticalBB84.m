function protocol = gen_protocol_opticalBB84(dim, pZ,pT, options)

    
    pAz = pZ;
    pAx = 1-pZ;
    pBz = pZ;
    pBx = 1-pBz;
    pX = 1-pZ;
 
    localPOVM = { pZ*diag([1,0,0]), pZ*diag([0,1,0]), pX*[1,1,0;1,1,0;0,0,0]/2, pX*[1,-1,0;-1,1,0;0,0,0]/2,diag([0,0,1])};
    sigmaX = [0,1,0;1,0,0;0,0,0];

    subID = diag([1,1,0]);
    subID2 = kron(subID,subID);
    phaseError = (subID2-kron(sigmaX,sigmaX))/2;
    XCorr =(subID2+kron(sigmaX,sigmaX))/2;
    unscaledPhaseErrorPOVM = {phaseError,XCorr,eye(dim^2)-subID2};
    
    Kpk = {};
    Kp = {};
  

        
        % cross-over min-tradeoff only
     
        %corrZZOp = kron(diag([1 0]),diag([1 0]))+ kron(diag([0 1]),diag([0 1]));
       
     	POVM_pe = unscaledPhaseErrorPOVM;
        %POVM_pe = {phaseErrorOp, eye(4)-phaseErrorOp};
        % XX test
    	Kp{end+1}  =       sqrt( pAx*pBx )*subID2;
     	Kpk{end+1} = cellsqrtkron( localPOVM(3:4),localPOVM(3:4) );  
     	% Alice no detection
      	Kp{end+1}  =       sqrt( pAx*pBx )*kron(diag([0,0,1]),subID);
       	Kpk{end+1} = pref( sqrt( pAx),cellsqrtkron( localPOVM(5),localPOVM(3:4))); 
        % Bob no detection
      	Kp{end+1}  =        sqrt( pAx*pBx )*kron(subID,diag([0,0,1]));
       	Kpk{end+1} = pref(sqrt( pBx),cellsqrtkron( localPOVM(3:4),localPOVM(5))); 
        % Both no detection
      	%Kp{end+1}  =   sqrt( pAx*pBx )*kron(diag([0,0,1]),diag([0,0,1]));
       	%Kpk{end+1} = pref(sqrt( pAx*pBx),cellsqrtkron( localPOVM(5),localPOVM(5)));     
        
        % ZZ generation
       	Kp{end+1}  =       sqrt( pAz*pBz )* subID2;
      	Kpk{end+1} =  cellsqrtkron( localPOVM(1:2),{pBz*subID} );       
   

    protocol = struct;    
    protocol.Kpk = Kpk;
    protocol.Kp = Kp;
    protocol.POVM_pe  = POVM_pe;
    protocol.dim = dim^2;
    protocol.pT = pT;
    protocol.AlicePOVM = localPOVM;
    protocol.BobPOVM = localPOVM;
end






