% Model accounts for time order effect in sequential comparison task

clear all

data_ID = '_SFC_DD_NM';

% load data
%-------------------

project_dir_SFC_DD = '/home/jan/data1/EEG/SFC_DD';
log_dir_SFC_DD = fullfile(project_dir_SFC_DD, 'log_files');

project_dir_SFC_DD_BP = '/home/jan/data1/EEG/SFC_DD_BP';
log_dir_SFC_DD_BP = fullfile(project_dir_SFC_DD_BP, 'log_files');

project_dir_SFC_DD_BP_NM = '/home/jan/data1/EEG/SFC_DD_BP_NM';
log_dir_SFC_DD_BP_NM = fullfile(project_dir_SFC_DD_BP_NM, 'log_files');

project_dir_SFC_DD_NM = '/home/jan/data1/EEG/SFC_DD_NM';
log_dir_SFC_DD_NM = fullfile(project_dir_SFC_DD_NM, 'log_files');


file_tag = '_log.mat'; %'_SFCwE_log.mat'; %'_dmts_dec_log.mat';
filenames = [];

subjects = {};

log_dir = eval(['log_dir' data_ID]);

if isempty(subjects)
    filenames = dir(fullfile(log_dir, ['*' file_tag]));
    filenames = {filenames.name};
else
    for i=1:length(subjects)
        filenames{i} = [subjects{i} file_tag];
    end
end


%%
% stimulus params
%------------

stim_set = [16 16 16 16 16 20 20 20 20 20 24 24 24 24 24 28 28 28 28 28;...
            12 14 16 18 20 16 18 20 22 24 20 22 24 26 28 24 26 28 30 32];
f_mean   = mean(stim_set(:));

plot_results    = 0;
invert_model    = 1;
plot_singTrials = 0;
n_rep           = 2000;

% param for simulation & prior means
%--------
sigma = 0.05;    % variance parameter of likelihood (default: 0.05)
lambda_stim = 1/sigma;
gamma = 1;      % relative scaling of likelihood/posterior variances [if inF.Gamma = 'posterior': post_var1 *= gamma, if inF.Gamma = 'likelihood': lhood_var1 *= gamma]

mu_prior = f_mean;

% variance of global prior: FWHM of the range of all frequencies
fwhm         = max(log(stim_set(2,:)))-min(log(stim_set(2,:)));
v_prior      = (fwhm/(2*sqrt(2*log(2))))^2;
% v_prior      = 0.05;
lambda_prior = 1/v_prior;

% bias term - phi parameter
bias = 0;

param2est     = [1 0 0 1];          % sigma, gamma, mu_prior, v_prior
phi_param2est = [1];

est_var       = 1;                  % variance of to-be-estimated posterior distributions
param2est     = est_var.*param2est;
phi_param2est = est_var.*phi_param2est;


%% Evolution and Observation functions
f_fname = @f_SFC_simple;              % evolution function perception process
g_fname = @g_SFC_simple_bias;         % whole decision process

%% loop over all subjects & invert model & save parameters

for n = 1:length(filenames)
    
    try
        disp(['==== Trying to load: ' fullfile(log_dir,filenames{n}) ' ===='])

        load(fullfile(log_dir,filenames{n}))

        stim_pairs                  = reshape(permute(mylog.flutter, [3 2 1]), 2, []);
        behav_resp                  = reshape(mylog.choice_content', 1, []);
        behav_resp(behav_resp == 0) = -1;           % set unanswered trials to -1
        behav_resp(behav_resp == 2) = 0;            % recode data from 1=f1>f2 & 2=f2>f1 to 1=f1>f2 & 0=f2>f1

        valid_idx = all([stim_pairs(1,:) ~= 0 ; behav_resp ~= -1]); 
    catch
        stim_pairs = [];
        behav_resp = [];
        disp(['>>>>>>>> NO DATA LOADED: ' filenames{n}])
        continue
    end
    
    if ~isempty(stim_pairs) && ~isempty(behav_resp)
        u = stim_pairs(:,valid_idx);    
        y = behav_resp(valid_idx);       % binary vector with ones for trials with choice f1 > f2 and zeros elsewhere
    else
        u = [];
        y = [];
    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % prior mean set
    %-----------------
    m1 = log(f_mean);    % prior mean for the first stim
    m2 = m1;             % prior mean for the second stim
    v1 = v_prior;        % prior variance for the first stim
    v2 = v1;             % prior variance for the second stim

    % Option for g and f functions
    %-------------------------------------
    inF.PriorF2   = 'none'; % 'none' (use likelihood of f2), 'PostF1' (use posterior of f1) or 'Same' (use same as for f1)
    inF.Gamma     = 'Likelihood'; % 'Posterior' or 'Likelihood'

    inF.Disp = plot_singTrials;% plot distribution
    if inF.Disp
        inF.Hdisp = figure;
    end

    % Theta, Phi, x0
    %----------------
    theta = [sigma; gamma; log(mu_prior); v_prior]; % f parameters
    phi   = bias; % g parameters
    fprintf('sigma:\t%1.3f\ngamma:\t%1.3f\nmu_p:\t%1.3f\nv_p:\t%1.3f\nbias:\t%1.3f\n', theta, phi)
    
    x0 = repmat([m1;v1;m2;v2],2,1); % hidden states priors
 
    %------------------------------------
    % choose initial conditions (priors)
    %------------------------------------

    dim = struct('n',length(x0),...     % number of hidden states in the probability learning model
        'n_theta',length(theta),...     % evolution parameters
        'n_phi',numel(phi));                     % observation parameters

    options = [];
    options.binomial = 1;
    
    priors.muPhi      = phi;
    priors.muTheta    = theta;
    priors.muX0       = x0;
    priors.SigmaX0    = 0*eye(dim.n);
    priors.SigmaPhi   = 0*eye(dim.n_phi);
    priors.SigmaTheta = 0*eye(dim.n_theta);
    priors.a_alpha    = Inf;
    priors.b_alpha    = 0;

    options.priors   = priors;
    options.inF      = inF;

    %% Simulate model and format data for further inversion
    %-------------------------------------------------------

    % create stimulus sequence
    if isempty(u)
        if plot_singTrials
            u = stim_set;        
        else
            Nt = n_rep;
            NStim   = size(stim_set,2);
            NRepet  = Nt/NStim;
            u = [];
            for i = 1:NRepet
                Irnd = randperm(NStim);
                u = [u stim_set(:,Irnd)];
            end
        end
    end

    Nt = length(u); % number of trials

    % simulate behavioral responses
    if isempty(y)

        alpha_simu = Inf; % precision of the stochastic innovations
        sigma_simu = Inf; % precision of the measurement error -> does not play a role in binomial model

        [y,x,x0,eta,e] = simulateNLSS(Nt,f_fname,g_fname,theta,phi,log(u),alpha_simu,sigma_simu,options,x0);
    end



    %% Results
    %---------------------------
 

    [cTOE(n,:),PerfGlobal] = Plot_TOE(u,y,plot_results);
 
    


    if invert_model
        %----------------------------------------------------------------------
        %% MODEL INVERSION
        %----------------------------------------------------------------------
        options.priors.a_alpha  = Inf;
        options.priors.b_alpha  = 0;
        options.priors.muX0     = x0;
        options.priors.SigmaX0  = diag([0 5 0 5 0 5 0 5]);
        
        options.priors.muPhi    = phi;
        options.priors.SigmaPhi = diag(phi_param2est);

        options.priors.muTheta    = theta;
        options.priors.SigmaTheta = diag(param2est);

        options.isYout = zeros(size(y,1),size(y,2));
        
        options.gradF       = 1;       % optimize free energy (1) or variational energy (0)
        options.backwardLag = 1; % size of short sighted backward passes (default: 1)
        
        % controlling convergence of VBA
        options.MaxIter = 10;   % default: 10
        options.TolFun = 1e-4;  % default: 1e-4
        
        % controlling convergence of inner loop Gauss-Newton
        options.GnMaxIter = 10;  % default: 10
        options.GnTolFun = 1e-4; % default: 1e-4

        [posterior,out] = VBA_NLStateSpaceModel(y,log(u),f_fname,g_fname,dim,options);

        if plot_results
            % display posterior means
            xtick = 1:out.dim.n_theta;
            if ~out.options.OnLine
                V = VBA_getVar(posterior.SigmaTheta);
                muTheta = posterior.muTheta;
            else
                V = VBA_getVar(posterior.SigmaTheta{end});
                muTheta = posterior.muTheta(:,end);
            end
            plotUncertainTimeSeries(muTheta,V,[]);
            set(gca, 'Xtick', xtick)
            title('theta')
            plot(theta,'go')
        end

        subj_posterior(n) = posterior;
        fitting_stats(n)  = out;
        subj_ID(n)        = str2num(subjects{n});
        
    end
end

save('SFC_DD_NM_posteriors_simple_sigma_v_prior_bias_std_erf_direct_last_4', 'subj_posterior', 'fitting_stats', 'subj_ID')

