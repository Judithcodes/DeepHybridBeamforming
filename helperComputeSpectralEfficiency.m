function R = helperComputeSpectralEfficiency(H,F,W,Ns,snr)
% This function helperComputeSpectralEfficiency is only in support of
% HybridPrecodingExample. It may change in a future release.

% Copyright 2017 The MathWorks, Inc.
H = H.';
F = F.';

% Heff = (F*H*W).';
% Weff = W.';
% R = log2(det(eye(Ns)+snr/Ns*(real(conj(Weff)*Weff.')\real(conj(Heff)*Heff.'))));

temp = F(1:Ns,:)*H*W(:,1:Ns);
R = log2(det(eye(Ns)+snr/Ns*(real(W(:,1:Ns)'*W(:,1:Ns))\real(temp'*temp))));
