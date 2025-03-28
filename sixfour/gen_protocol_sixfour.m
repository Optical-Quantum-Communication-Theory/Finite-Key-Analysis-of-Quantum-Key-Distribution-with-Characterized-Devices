%%
% Protocol setup for six-state four-state protocol
% This is only for cross-over min-tradeoff functions
function protocol = gen_protocol_sixfour( pAList, pBList,pT, options )
   
    pAz = pAList(1);
    pAx = pAList(2);
    pAy = 1-pAz - pAx;
    pBz = pBList(1);
    pBx = 1-pBz;
    
    ZKraus = { diag([1,0]),diag([0,1]) };
    XKraus = { [1;1]*[1,1]/2,[1;-1]*[1,-1]/2 };
    YKraus = { [1;1i]*[1,-1i]/2,[1;-1i]*[1,1i]/2 };
    
    Kpk = {};
    Kp = {};
    sigmaX = [0,1;1,0];
    AlicePOVM = { pAz*diag([1 0]), pAz*diag([0 1]), pAx*[1 1;1 1]/2, pAx*[1 -1;-1 1]/2, pAy*[1,-1i;1i,1]/2, pAy*[1,1i;-1i,1]/2} ;
    BobPOVM = { pBz*diag([1 0]), pBz*diag([0 1]), pBx*[1 1;1 1]/2, pBx*[1 -1;-1 1]/2} ;
 
        
    % cross-over min-tradeoff only
    corrYXOp = kron([1,-1i;1i,1]/2,[1,1;1,1]/2)+ kron([1,1i;-1i,1]/2,[1 -1;-1 1]/2);
    %corrZZOp = kron(diag([1 0]),diag([1 0]))+ kron(diag([0 1]),diag([0 1]));
	
    phaseErrorXXOp = (eye(4) - kron(sigmaX,sigmaX))/2;
        
    prob_all = pAList(2)*pBList(2)+pAList(3)*pBList(2);
    POVM_pe = {pAList(2)*pBList(2)/prob_all*phaseErrorXXOp, pAList(2)*pBList(2)/prob_all*(eye(4)-phaseErrorXXOp),...
                pAList(3)*pBList(2)/prob_all*corrYXOp, pAList(3)*pBList(2)/prob_all*(eye(4)-corrYXOp)};
    %POVM_pe = {phaseErrorOp, eye(4)-phaseErrorOp};
    
    % XX test
    Kp{end+1}  =       sqrt( pAx*pBx )* eye(4);
    Kpk{end+1} = pref( sqrt( pAx*pBx ), cellsqrtkron( XKraus,XKraus ));  
    % YX test
    Kp{end+1}  =       sqrt( pAy*pBx )* eye(4);
    Kpk{end+1} = pref( sqrt( pAy*pBx ), cellsqrtkron( YKraus,XKraus ));        
    % ZZ generation
    Kp{end+1}  =       sqrt( pAz*pBz )* eye(4);
    Kpk{end+1} = pref( sqrt( pAz*pBz ), cellsqrtkron( ZKraus,{eye(2)} ));       
 
    
    
    protocol = struct;    
    protocol.Kpk = Kpk;
    protocol.Kp = Kp;   
    protocol.POVM_pe  = POVM_pe;
    protocol.dim = 4;
    protocol.pT = pT;
    protocol.AlicePOVM =AlicePOVM;
    protocol.BobPOVM = BobPOVM;
end
