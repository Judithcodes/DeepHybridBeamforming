function [Fbb,Frf,Wbb,Wrf] = helperOMPHybridWeights(Hin,NtRF,NrRF,Ns,At,Ar,noisevar)
% This function helperOMPHybridWeights is only in support of
% HybridPrecodingExample. It may change in a future release.

% Copyright 2017 The MathWorks, Inc.

% use paper convention, convert from Comm convention
% H = Hin.';
H = Hin;
[~,~,v] = svd(H);
Fopt = v(:,1:Ns);
[Nr,Nt] = size(H);

Frf = complex(zeros(Nt,NtRF));
Fres = Fopt;
for m = 1:NtRF
    Psi = At'*Fres;
    [~,k] = max(diag(Psi*Psi'));
    Frf(:,m) = At(:,k);
    Fbb = (Frf(:,1:m)'*Frf(:,1:m))\Frf(:,1:m)'*Fopt;
    temp = Fopt-Frf(:,1:m)*Fbb;
    Fres = temp/norm(temp,'fro');
end
Fbb = sqrt(Ns)*Fbb/norm(Frf*Fbb,'fro');

Wmmse = ((Fbb'*Frf'*(H'*H)*Frf*Fbb+noisevar*Ns*eye(Ns))\Fbb'*Frf'*H')';
%Wmmse = (1/Ns*(Fbb'*Frf'*H')/(1/Ns*H*(Frf*(Fbb*Fbb')*Frf')*H'+noisevar*eye(Nr)))';
Wrf = complex(zeros(Nr,NrRF));
Wres = Wmmse;
Ess = 1/Ns*eye(Ns);
Eyy = H*Frf*Fbb*Ess*Fbb'*Frf'*H'+noisevar*eye(Nr);
for m = 1:NrRF
    Psi = Ar'*Eyy*Wres;
    [~,k] = max(diag(Psi*Psi'));
    Wrf(:,m) = Ar(:,k);
    Wbb = (Wrf(:,1:m)'*Eyy*Wrf(:,1:m))\(Wrf(:,1:m)'*Eyy*Wmmse);
    temp = Wmmse-Wrf(:,1:m)*Wbb;
    Wres = temp/norm(temp,'fro');
end

% % sort based on diagonal term
% ChanEff = Wbb'*Wrf'*H*Frf*Fbb;
% [~,idx] = sort(diag(ChanEff),'descend');
% Fbb = Fbb(:,idx);
% Frf = Frf(:,idx);
% Wrf = Wrf(:,idx);
% Wbb = Wbb(:,idx);


% convert back to comm convention
% match Wbb'*Wrf'*H.'*Frf*Fbb*X.' to X*Fbb.'*Frf.'*H*conj(Wrf)*conj(Wbb) to
% X*Fbb*Frf*H*Wrf*Wbb

% Fbb = Fbb.';
% Frf = Frf.';
Wrf = conj(Wrf);
Wbb = conj(Wbb);

