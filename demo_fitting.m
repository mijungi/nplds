%% a demo to demonstrate how to use code for PLDS with nonstationary firing rates

% mijung wrote on Oct 2, 2015

clear all;
close all;
clc;

addpath core_functions/

% essential quantities (fixing pinp==0)
k = 2; % dimensionality of hidden states
p = 10; % output dimension (e.g., # of neurons)
d = 4; % input dimension
T = 100; % time 1:T (for training)
r = 10; % # of recordings

%% defining hyper-parameters and draw samples for parameters.

Model =  'NSFR';

% generate parameters 
ind_train = 1:r;
tau_init = 10*rand; 
params = generate_params_from_prior(k, d, p, r, Model, ind_train, tau_init);

% generate inputs 
params.inpn = 0.5*randn(d, T);

params.ind_train = ind_train; % all of them are used for training 
params.tau_init = tau_init; 
params.m_h_init = zeros(k,1);
params.sig_init = 10*abs(rand); 

% generate data
xyzinpn = generate_data_PLDS_multiple_recordings(T, params, r, Model);

meanfiringrate = zeros(p,r);
for i=1:r
    meanfiringrate(:,i) = mean(xyzinpn{i}.y,2);
end

figure(1);
subplot(2,4,[1 2 3 4]); plot(meanfiringrate, 'o-'); title('mean firing rate');
subplot(2,4,5); hinton(params.A,['true A1  ' num2str(max(max(params.A)),'%.4f')],'standard');
subplot(2,4,6); hinton(params.B,['true B  ' num2str(max(max(params.B)),'%.4f')],'standard');
subplot(2,4,7); hinton(params.C,['true C  ' num2str(max(max(params.C)),'%.4f')],'standard');
subplot(2,4,8); plot(params.h'); title(' true h');

%% VBEM

params.maxIter = 5; % define how many EM iterations you want to run 
datastruct = VBEM_PLDSnonstationary(xyzinpn, r, params, Model); 

%% generate samples for Y and Z (log firing rate) from the model with estimated parameters, and true parameters
% for sanity check

nsamps = 500;

% assuming T is same for all recordings
Zn = zeros(p, r*T, nsamps);
Yn = zeros(p, r*T, nsamps);
Yn_estParam = zeros(p, r*T, nsamps);
Zn_estParam = zeros(p, r*T, nsamps);

datastruct.Mstep = datastruct.Mstep{end}; % take the results from the last iteration
datastruct.Mstep.inpn = params.inpn;

for i=1:nsamps
    fprintf(['Sample generation: ' num2str(i) 'th sample, out of ' num2str(nsamps) 'samples' '\n'])
    
    sanitycheck = generate_data_PLDS_multiple_recordings(T,datastruct.Mstep,r,Model);
    sanitychecktrue = generate_data_PLDS_multiple_recordings(T,params,r,Model);
    
    for rr = 1: r
        Zn_estParam(:, (rr-1)*T+1:rr*T, i) = sanitycheck{rr}.z;
        Yn_estParam(:, (rr-1)*T+1:rr*T, i) = sanitycheck{rr}.y;
        Zn(:,(rr-1)*T+1:rr*T,i) = sanitychecktrue{rr}.z;
        Yn(:,(rr-1)*T+1:rr*T,i) = sanitychecktrue{rr}.y;
    end
end

%  then, check the mean firing rates, covariances, etc. 

% mean vector
meanvec = zeros(p,r);
meanvec_est = zeros(p,r);

covmat = zeros(p, p, r);
covmat_est = zeros(p, p, r);

for rr = 1: r

    %% test mean
    meanvec(:,rr) = mean(mean(Zn(:,(rr-1)*T+1:rr*T,:),3),2);
    meanvec_est(:,rr) = mean(mean(Zn_estParam(:, (rr-1)*T+1:rr*T, :),3),2);
    
    %% test covariance
    
    Znrshp = reshape(Zn(:,(rr-1)*T+1:rr*T,:), p, []);
    Znestrshp = reshape(Zn_estParam(:,(rr-1)*T+1:rr*T,:), p, []);
    
    covmat(:,:,rr) = cov(Znrshp');
    covmat_est(:,:,rr) = cov(Znestrshp');
    
    % onestep delayed covariance
    Znrshp2 = Znrshp(:,[2:end,1]);
    Znestrshp2 = Znestrshp(:, [2:end,1]);
        
    %% figure
    figure;
    subplot(2,2,[1 2]); plot([meanvec(:,rr) meanvec_est(:,rr)]); legend('mean of Z using true params', 'mean of Z using estimated params'); xlabel('neurons');
    subplot(2,2,3); hinton(covmat(:,:,rr),['cov z' num2str(max(max(covmat(:,:,rr))),'%.4f')],'standard');
    subplot(2,2,4); hinton(covmat_est(:,:,rr),['cov z est  ' num2str(max(max(covmat_est(:,:,rr))),'%.4f')],'standard');
        
end
