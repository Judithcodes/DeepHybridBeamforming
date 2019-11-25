function [Fbb,Frf] = helperOMPHybridPrecodingWeights(H,NtRF,Ns,At)
% This function helperOMPHybridPrecodingWeights is only in support of
% HybridPrecodingExample. It may change in a future release.

% Copyright 2017 The MathWorks, Inc.

% use Comms convention
[~,~,v] = svd(H.');
Fopt = v(:,1:Ns);
Nt = size(H,1);

Frf = complex(zeros(Nt,NtRF));
Fres = Fopt;
for m = 1:NtRF
    Psi = At'*Fres;
    [~,k] = max(diag(Psi*Psi'));
    Frf(:,m) = At(:,k);
    Fbb = (Frf(:,1:m)'*Frf(:,1:m))\(Frf(:,1:m)'*Fopt);
    temp = Fopt-Frf(:,1:m)*Fbb;
    Fres = temp/norm(temp,'fro');
end
Fbb = sqrt(Ns)*Fbb/norm(Frf*Fbb,'fro');

% match H.'*Frf*Fbb*X.' to X*Fbb*Frf*H

Fbb = Fbb.';
Frf = Frf.';

