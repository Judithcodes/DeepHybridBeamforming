function [Frf,Fbb] = findFrfFbb(Hin,Ns,NtRF,At)

% H = Hin.';
H = Hin;
[~,Sigma,v] = svd(H);
Fopt = v(:,1:Ns);
[Nr,Nt] = size(H);

Frf = At;
Fbb = (Frf'*Frf)\Frf'*Fopt;
Fbb = sqrt(Ns)*Fbb/norm(Frf*Fbb,'fro');



% 
% 
% Frf = complex(zeros(Nt,NtRF));
% Fres = Fopt;
% for m = 1:NtRF
%     Psi = At'*Fres;
%     [~,k] = max(diag(Psi*Psi'));
%     Frf(:,m) = At(:,k);
%     Fbb = (Frf(:,1:m)'*Frf(:,1:m))\Frf(:,1:m)'*Fopt;
%     temp = Fopt-Frf(:,1:m)*Fbb;
%     Fres = temp/norm(temp,'fro');
% end
% Fbb = sqrt(Ns)*Fbb/norm(Frf*Fbb,'fro');
% 
% 
% 
% 


% 
% 
% Fbb = Fbb.';
% Frf = Frf.';




