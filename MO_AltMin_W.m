function [Wfull, WRF, WBB ] = MO_AltMin_W( Wopt,Fopt,H, NRF, noisevar, epsilon )






% addpath('./manopt');
% D:\myResearch_101\MIMODeepUnfolding\TWC\AltMin\Narrowband\MO-AltMin\manopt
[Nr, Ns] = size(Wopt);
y = [];
WRF = exp( 1i*unifrnd(0,2*pi,Nr,NRF) );
% FRF = PerEntryQuantizationNoZeros(Fopt,Nt);
while(isempty(y) || abs(y(1)-y(2))>epsilon)
    
    Wmmse = ((Fopt'*(H'*H)*Fopt+noisevar*Ns*eye(Ns))\Fopt'*H')';
    Ess = 1/Ns*eye(Ns);
    Eyy = H*Fopt*Ess*Fopt'*H' + noisevar*eye(Nr);
    WBB = (WRF'*Eyy*WRF)\(WRF'*Eyy*Wmmse);
%     WRF = conj(WRF);
%     WBB = conj(WBB);
%     WBB = pinv(WRF) * Wopt;
    %     FBB = sqrt(Ns)*FBB/norm(FRF*FBB,'fro');
    y(1) = norm(Wopt - WRF * WBB,'fro')^2;
    [WRF, y(2)] = sig_manif(Wopt, WRF, WBB);
%     WRF = conj(WRF);
end
Wfull = WRF*WBB;
Wfull = Wfull./sqrt(sum(diag(Wfull'*Wfull))/NRF);
end