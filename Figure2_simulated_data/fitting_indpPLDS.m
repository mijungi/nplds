% fitting independent PLDS
% mijung  edited on Oct 3, 2015

clear all;
close all;
clc;

seed = 1;
oldRng = rng();
rng(seed);


addpath ../core_functions/

%%
load mat_files/all_NSFR.mat


params = cell(r,1);
xyzinpn = cell(r,1);

for i=1:r
    % put all the true params
    params{i}.A = A;
    params{i}.B = B;
    params{i}.x0 = x0;
    params{i}.V0 = V0;
    params{i}.Q = eye(k);
    params{i}.d = d;
    params{i}.C = C;
    params{i}.h = h(:,i);
    params{i}.inpn = inpn;
    
    params{i}.ind_train = 1; % all of them are used for training
    params{i}.tau_init = 10*rand;
    params{i}.m_h_init = zeros(k,1);
    params{i}.sig_init = 10*abs(rand);
      
    %     xyzinpn{i}.x = x(:,:,i);
    xyzinpn{i}.y = y(:,:,i);
    xyzinpn{i}.z = z(:,:,i);
    xyzinpn{i}.T = T;
    xyzinpn{i}.inpn = inpn ;
end

%% run VBEM code ...

datastruct = cell(r,1);
for i = 1:r
    Model = 'NSFR';
    params{i}.maxIter = 10;
    datastruct{i} = VBEM_PLDSnonstationary(xyzinpn{i}, 1, params{i}, Model);
end

save mat_files/fitting_indpPLDS.mat

%% visualising the results 

load mat_files/fitting_indpPLDS.mat

% (1) log firing rate

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
    
    CC = datastruct{trial_to_check}.Mstep{end}.C;
    hh = datastruct{trial_to_check}.Mstep{end}.h;
    covhh = datastruct{trial_to_check}.Mstep{end}.covh;
    mu = datastruct{trial_to_check}.Estep{end}.mumarg;
    
    Cmud = zeros(p, T);
    for t=1:T
        Cmud(:,t) = CC*(mu(:,t)+hh) + datastruct{trial_to_check}.Mstep{end}.d;
    end
    
    zestmat(:,:,trial_to_check) = Cmud;
    
    z1_est(trial_to_check,:) = sum(Cmud(1:p/2,:));
    z2_est(trial_to_check,:) = sum(Cmud(p/2+1:p,:));
    
    corr_z_est(trial_to_check) = corr(z1_est(trial_to_check,:)', z2_est(trial_to_check,:)');
    
    invsig = datastruct{trial_to_check}.Estep{end}.inv_sigmarg;
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

% numtimebins = T/4;
autocorr = zeros(p, numtimebins+1);

for whichcell = 1:p
    %    autocorr(whichcell, :) = xcov(yy(whichcell,:), numtimebins/2, 'unbiased');
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

