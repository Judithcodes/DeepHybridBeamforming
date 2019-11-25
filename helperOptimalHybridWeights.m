function [Fopt,Wopt] = helperOptimalHybridWeights(Hin,Ns,noisevar)
% This function helperOptimalgWeights is only in support of
% HybridPrecodingExample. It may change in a future release.

% Copyright 2017 The MathWorks, Inc.

% use paper convention, convert from Comm convention
% H = Hin.'; % H is Nr x Nt in paper.
H = Hin;
[~,~,v] = svd(H);
Fopt = v(:,1:Ns);

Wopt = ((Fopt'*(H'*H)*Fopt+noisevar*Ns*eye(Ns))\(Fopt'*H'))';

% convert back to comm convention
% match Wopt'*H.'*Fopt*X.' to X*Fopt.'*H*conj(Wopt) to
% X*Fopt*H*Wopt

% Fopt = Fopt.';
Wopt = conj(Wopt);

