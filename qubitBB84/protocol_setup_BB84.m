% set up qubit-based BB84 
function protocol = protocol_setup_BB84(pZ,pT, options) 
    % parameters: pZ  -- probability of choosing the Z basis (key-generation basis)
    %             pT  -- probability of testing (not used if sifting and testing is combined)
    %             options -- additional protocol choices
    
 
    pX = 1-pZ;

    ZKraus = { diag([1,0]),diag([0,1]) };
    XKraus = { [1;1]*[1,1]/2,[1;-1]*[1,-1]/2 };
    
    Kpk = {};
    Kp = {};
    sigmaX = [0,1;1,0];
    localPOVM = { pZ*diag([1 0]), pZ*diag([0 1]), pX*[1 1;1 1]/2, pX*[1 -1;-1 1]/2};

    phaseErrorOp = (eye(4) - kron(sigmaX,sigmaX))/2;
    POVM_pe = {phaseErrorOp, eye(4)-phaseErrorOp};
    % XX test -  testing probability is determined by basis choice
    %               probability
    Kp{end+1}  =       sqrt( pX^2 )* eye(4);
    Kpk{end+1} = pref( sqrt( pX^2 ), cellkron( XKraus,XKraus ));        
    % ZZ generation
    Kp{end+1}  =       sqrt( pZ^2 )* eye(4);
    Kpk{end+1} = pref( sqrt( pZ^2 ), cellkron( ZKraus,{eye(2)} ));       
  
       
    % construct the protocol struct
    protocol = struct;
    protocol.Kpk = Kpk;
    protocol.Kp = Kp;
    protocol.POVM_pe  = POVM_pe;
    protocol.dim = 4; % dimension
    protocol.pT = pT;
    protocol.AlicePOVM = localPOVM;
    protocol.BobPOVM = localPOVM;
end