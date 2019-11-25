function [Ffull, FRF,FBB ] = MO_AltMin( Fopt, NRF, epsilon )

% addpath('./manopt');
% D:\myResearch_101\MIMODeepUnfolding\TWC\AltMin\Narrowband\MO-AltMin\manopt
[Nt, Ns] = size(Fopt);
y = [];
FRF = exp( 1i*unifrnd(0,2*pi,Nt,NRF) );
% FRF = PerEntryQuantizationNoZeros(Fopt,Nt);
while(isempty(y) || abs(y(1)-y(2))>epsilon)
    FBB = pinv(FRF) * Fopt;
%     FBB = sqrt(Ns)*FBB/norm(FRF*FBB,'fro');
    y(1) = norm(Fopt - FRF * FBB,'fro')^2;
    [FRF, y(2)] = sig_manif(Fopt, FRF, FBB);
end
Ffull = FRF*FBB;
Ffull = Ffull./sqrt(sum(diag(Ffull'*Ffull))/NRF);
end