function [Wrf,Wbb] = findWrfWbb(Hin,Ns,NrRF,Ar,noisevar)
% use paper convention, convert from Comm convention
% H = Hin.';
H = Hin;
[~,~,v] = svd(H);
Fopt = v(:,1:Ns);
[Nr,Nt] = size(H);
% Fbb = Fbb.';
% Frf = Frf.';



% 
% Wmmse = ((Fbb'*Frf'*(H'*H)*Frf*Fbb+noisevar*Ns*eye(Ns))\Fbb'*Frf'*H')';
% Ess = 1/Ns*eye(Ns);
% Eyy = H*Frf*Fbb*Ess*Fbb'*Frf'*H'+noisevar*eye(Nr);
% Wrf = Ar;
% Wbb = (Wrf'*Eyy*Wrf)\(Wrf'*Eyy*Wmmse);

Wmmse = ((Fopt'*(H'*H)*Fopt+noisevar*Ns*eye(Ns))\Fopt'*H')';
Ess = 1/Ns*eye(Ns);
Eyy = H*Fopt*Ess*Fopt'*H'+noisevar*eye(Nr);
Wrf = Ar;
Wbb = (Wrf'*Eyy*Wrf)\(Wrf'*Eyy*Wmmse);

Wrf = conj(Wrf);
Wbb = conj(Wbb);

