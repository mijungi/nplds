%% to fit NPLDS
% mijung  edited on Oct 3, 2015

clear all;
close all;
clc;

seed = 1;
oldRng = rng();
rng(seed);


addpath ../core_functions/

%%
% load data
load mat_files/all_NSFR.mat

r = 100;

% put all the true params
params.A = A;
params.B = B;
params.C = C;
params.d = d;
params.h = h;

params.x0 = x0;
params.V0 = V0;
params.Q = eye(k);

newdataset = cell(r,1);

for i=1:r
    newdataset{i} = xyzinpn{i};
end

xyzinpn = newdataset;

params.ind_train = 1:r; % all of them are used for training 
params.tau_init = 10*rand; 
params.m_h_init = zeros(k,1);
params.sig_init = 10*abs(rand); 

%% start fitting here

Model = 'NSFR';
params.maxIter = 10; 
datastruct = VBEM_PLDSnonstationary(xyzinpn, r, params, Model); 

save mat_files/fitting_NPLDS.mat

%% visualising results

load mat_files/fitting_NPLDS.mat

datastruct.Mstep = datastruct.Mstep{end};
datastruct.Estep = datastruct.Estep{end};

% (1) log mean firing rates

zz = zeros(p, T*size(y,3));
for i=1:size(y,3)
    zz(:, 1+(i-1)*T:i*T) = z(:,:,i);
end

z1_est = zeros(r, T);
z2_est = zeros(r, T);

z1 = zeros(r, T);
z2 = zeros(r, T);

corr_z_est = zeros(r,1);
corr_z = zeros(r, 1);

zestmat = zeros(p, T, r);

for trial_to_check = 1:r
    
    CC = datastruct.Mstep.C;
    hh = datastruct.Mstep.h(:, trial_to_check);
    covhh = datastruct.Mstep.covh(:,:,trial_to_check);
    mu = datastruct.Estep{trial_to_check}.mumarg;
    
    Cmud = zeros(p, T);
    for t=1:T
        Cmud(:,t) = CC*(mu(:,t)+hh) + datastruct.Mstep.d;
    end
    
    z1_est(trial_to_check,:) = sum(Cmud(1:p/2,:));
    z2_est(trial_to_check,:) = sum(Cmud(p/2+1:p,:));
    
    zestmat(:,:,trial_to_check) = Cmud;
     
    corr_z_est(trial_to_check) = corr(z1_est(trial_to_check,:)', z2_est(trial_to_check,:)');
    
    invsig = datastruct.Estep{trial_to_check}.inv_sigmarg;
    covz_errbar = zeros(p, p, T);
    errorbar = zeros(p, T);
    for t=1:T
        covz_errbar(:,:,t) = CC*(inv(invsig(:,:,t))+covhh)*CC';
        errorbar(:,t) =  diag(covz_errbar(:,:,t));
    end
    
    z1(trial_to_check,:) = sum(xyzinpn{trial_to_check}.z(1:p/2,:,1));
    z2(trial_to_check,:) = sum(xyzinpn{trial_to_check}.z(p/2+1:p,:,1));
    
    corr_z(trial_to_check) = corr(z1(trial_to_check,:)', z2(trial_to_check,:)'); 

end

figure
subplot(211);
plot(1:r, mean(z1,2)/(p/2),'r', 1:r, mean(z2,2)/(p/2), 'b', ...
    1:r, mean(z1_est,2)/(p/2), 'r--', 1:r, mean(z2_est,2)/(p/2), 'b--')
set(gca, 'ylim', [-3.0 -0.5]); legend('true z (grp1)', 'estimated z(grp1)', 'true z (grp2)', 'estimated z(grp2)');

% (2) covariances 

numtimebins = T;
autocorr = zeros(p, numtimebins+1);

for whichcell = 1:p
    autocorr(whichcell, :) = xcov(zz(whichcell,:), numtimebins/2, 'unbiased');
end

avgautocorr_acrcells = mean(autocorr);

autocorr_condi = zeros(p, numtimebins+1, r);

for whichcell = 1:p
    
    for whichtrial = 1:r
        autocorr_condi(whichcell, :, whichtrial) = xcov(z(whichcell,:,whichtrial), numtimebins/2, 'unbiased');
    end
    
end

autocorr_per_eachtrial = squeeze(mean(autocorr_condi));
avgautocorr_condi = mean(autocorr_per_eachtrial,2);

figure
subplot(211); 
plot(1:numtimebins+1, avgautocorr_acrcells/max(avgautocorr_acrcells), 'k', 1:numtimebins+1, avgautocorr_condi/max(avgautocorr_condi), 'r')
legend('total covariance', 'conditional covariance');

zzz = [];
for i=1:r
    zzz = [zz zestmat(:,:,i)];
end

autocorr = zeros(p, numtimebins+1);

for whichcell = 1:p
    autocorr(whichcell, :) = xcov(zzz(whichcell,:), numtimebins/2, 'unbiased');
end

avgautocorr_acrcells = mean(autocorr);

autocorr_condi = zeros(p, numtimebins+1, r);

for whichcell = 1:p
    
    for whichtrial = 1:r
        autocorr_condi(whichcell, :, whichtrial) = xcov(zestmat(whichcell,:,whichtrial), numtimebins/2, 'unbiased');
    end
    
end

autocorr_per_eachtrial = squeeze(mean(autocorr_condi));
avgautocorr_condi = mean(autocorr_per_eachtrial,2);

subplot(212);
plot(1:numtimebins+1, avgautocorr_acrcells/max(avgautocorr_acrcells), 'k--', 1:numtimebins+1, avgautocorr_condi/max(avgautocorr_condi), 'r--')
legend('estimate of total covariance', 'estimate of conditional covariance');

