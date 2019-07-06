function [gx] = g_SFC_simple_bias(x,Phi,~,~)

% [gx] = g_TwoFreqDiscrim_Physio(x,Phi,u,in)
%
% Computes the probability of the subject choosing 's1 > s2' and
% derives the action of the subject following a softmax rule.
% IN:
%   - x: the state posterior sufficient statistics.
%   - in: options set in options.inG
% OUT:
%   - gx: the predicted subject's output.

% up-dated conditional sufficient statistics on states

m1p = x(1);
v1p = x(2);
m2p = x(3);
v2p = x(4);
m1  = x(5);
v1  = x(6);
m2  = x(7);
v2  = x(8);

bias = Phi(1);

%% Behavioral response
%-----------------------

Dm = m1p - m2p;     % difference between posterior means
Sv = v1p + v2p;     % sum of posterior variances


gx = 0.5*(1 + erf(Dm/(sqrt(2)*sqrt(Sv)))) + bias; % Prob F1>F2


end










