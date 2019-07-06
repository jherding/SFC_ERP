function [fx] = f_SFC_simple_null(x,P,u,in)

% Computes VB update rules for hidden states sufficient statistics
%
% [fx] = f_SFC_simple_null(x,P,u,in)
%
% This is the update rule for the posterior sufficient
% statistics of a 2 interval AFC task (tactile frequency discrimination)
%
% IN:
%   - x: states
%       m1 = x(1);
%       v1 = x(2);
%       m2 = x(3);
%       v2 = x(4);
%       m1 = x(5);
%       v1 = 
%   - P: the perceptual model parameter vector, ie. P = [log(sigma),log(gamma),log(mu_global),log(v_global)]
%   - u: the two trial inputs (frequency of stim 1 and stim 2)
%   - in: options set in options.inF
% OUT:
%   - fx: the updated posterior sufficient statistics (having observed u),
%   as well as book keeping of the past believes.
%
% Jeremie - 06/03/2013
% Jan     - 02/16/2015

% stimulus vals
f1 = u(1);          % f1 (freq on log scale)
f2 = u(2);          % f2 (freq on log scale)

% model parameters
%--------------------
sigma     = P(1);      % precision of encoding (variance of stimulus likelihood for f1 and f2)
gamma     = P(2);      % decay parameter (multiplied with variance of f1 likelihood or posterior to model increase in memory uncertainty due to retention)
mu_global = P(3);      % mean of global prior (on log scale)
v_global  = P(4);      % variance of global prior

% model states
% --------------
m1  = mu_global;
v1  = v_global;
m2 = m1;
v2 = v1;

%% Up-dating the hidden states pertaining to the first stimulation


v1p = sigma; % variance of posterior for f1
m1p = f1;    % mean of posterior for f1 (freq on log scale)

m2p  = f2;                           % mean of likelihood of f2 serves as mean of posterior for f2
v2p  = sigma;                        % variance likelihood of f2 serves as variance of posterior for f2

%% Conditional posteriors on states
fx = zeros(8,1);
fx(1) = m1p;
fx(2) = v1p;
fx(3) = m2p;
fx(4) = v2p;
fx(5) = m1;
fx(6) = v1;
fx(7) = m2;
fx(8) = v2;

disp = in.Disp;
if disp
%% Display
%------------

xfig = log(10:0.1:50);

Yprior_f1 = normpdf(xfig,m1,sqrt(v1));           % prior of f1
Yprior_f1 = Yprior_f1/sum(Yprior_f1);   % ... normalized
Ylike_f1  = normpdf(xfig,f1,sqrt(v_lhood1));         % Likelihood of f1
Ylike_f1  = Ylike_f1/sum(Ylike_f1);     % ... normalized
Ypost_f1  = normpdf(xfig,m1p,sqrt(v1p));         % posterior of f1
Ypost_f1  = Ypost_f1/sum(Ypost_f1);     % ... normalized

% ... same procedure for f2
Yprior_f2 = normpdf(xfig,m2,sqrt(v2)); 
Yprior_f2 = Yprior_f2/sum(Yprior_f2);
Ylike_f2  = normpdf(xfig,f2,sqrt(v_lhood2));
Ylike_f2  = Ylike_f2/sum(Ylike_f2);
Ypost_f2  = normpdf(xfig,m2p,sqrt(v2p));
Ypost_f2  = Ypost_f2/sum(Ypost_f2);
% ypost_tmp = Yprior_f2.*Ylike_f2;
% ypost_tmp = ypost_tmp/(sum(ypost_tmp));


%% Behavioral response
%-----------------------

Dm = m1p - m2p;                             % difference of posterior means
Sv = v1p + v2p;                             % sum of posterior variance
chose_f1 = 0.5*(1 + erf(Dm/(sqrt(2)*sqrt(Sv))));  % Prob F1>F2

subj_diff = exp(m2p)-exp(m1p); % in Hz

figure(in.Hdisp);
set(gcf,'color','white');

% stim1
subplot(1,3,1);
plot(xfig,Yprior_f1,'b-','LineWidth',3); hold on
plot(xfig,Ylike_f1,'k-','LineWidth',3);
plot(xfig,Ypost_f1,'r--','LineWidth',3);
title(sprintf('f1 = %2.1f Hz', exp(f1)),'FontSize',10,'FontWeight','bold');
legend('Prior','Likelihood','Posterior');
set(gca,'FontSize',16,'FontWeight','bold');
axis([xfig(1) xfig(end) 0 0.03]);
hold off

% stim 2
subplot(1,3,2);
plot(xfig,Yprior_f2,'b-','LineWidth',3); hold on
plot(xfig,Ylike_f2,'k-','LineWidth',3);
plot(xfig,Ypost_f2,'r--','LineWidth',3);
% plot(xfig,ypost_tmp,'g:','Linewidth', 4);
title(sprintf('f2 = %2.1f Hz', exp(f2)),'FontSize',10,'FontWeight','bold');
legend('Prior','Likelihood','Posterior');
set(gca,'FontSize',16,'FontWeight','bold');
axis([xfig(1) xfig(end) 0 0.03]);
hold off

% Compare Posterior
subplot(1,3,3);
plot(xfig,Ypost_f1,'b-','LineWidth',3);
hold on
plot(xfig,Ypost_f2,'r-','LineWidth',3);
plot([f1 f1],[0 max(Ypost_f1)],'b--','LineWidth',2);
plot([f2 f2],[0 max(Ypost_f2)],'r--','LineWidth',2);
legend('Post. f1','Post. f2','f1','f2');
text(log(35), 0.03, sprintf('diff: %1.2f\nprob: %1.2f', m2p - m1p, chose_f1))
title(sprintf('f2 - f1 = %2.1f Hz', subj_diff),'FontSize',10,'FontWeight','bold');
set(gca,'FontSize',16,'FontWeight','bold');
axis([xfig(1) xfig(end) 0 0.03]);
hold off
drawnow

pause;
%close(h1);
else
end


