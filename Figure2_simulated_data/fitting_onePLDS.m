%% to fit one PLDS
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

r = 1;

% put all the true params
params.A = A;
params.B = B;
params.C = C;
params.d = d;
params.h = h;

params.x0 = x0;
params.V0 = V0;
params.Q = eye(k);

xyzinpn.inpn = repmat(inpn, [1 100]);

yy = zeros(p, T*size(y,3));
for i=1:size(y,3)
    yy(:, 1+(i-1)*T:i*T) = y(:,:,i);
end

xyzinpn.y = yy;
% xyzinpn.z = zz; 
params.inpn = xyzinpn.inpn;

params.ind_train = 1:r; % all of them are used for training 
params.tau_init = 10*rand; 
params.m_h_init = zeros(k,1);
params.sig_init = 10*abs(rand); 

%% start fitting here


Model = 'NSFR';
params.maxIter = 10; 
datastruct = VBEM_PLDSnonstationary(xyzinpn, r, params, Model); 

save mat_files/fitting_onePLDS.mat


%% visualising results 

load mat_files/fitting_onePLDS.mat

fromMstep = datastruct.Mstep{end};
fromEstep = datastruct.Estep{end};

% (1) log firing rates

zz = zeros(p, T*size(y,3));
for i=1:size(y,3)
    zz(:, 1+(i-1)*T:i*T) = z(:,:,i);
end

xyzinpn.z = zz;

CC = fromMstep.C;
hh = fromMstep.h;
covhh = fromMstep.covh;
mu = fromEstep.mumarg;

numsubsets = 100; 

Cmud = zeros(p, numsubsets*T);
for t=1:numsubsets*T
    Cmud(:,t) = CC*(mu(:,t)+hh) + fromMstep.d;
end


z1_est = zeros(numsubsets, T);
z2_est = zeros(numsubsets, T);

z1 = zeros(numsubsets, T);
z2 = zeros(numsubsets, T);

corr_z_est = zeros(numsubsets,1);
corr_z = zeros(numsubsets, 1);

zestmat = zeros(p, T, numsubsets); 

invsigmat = fromEstep.inv_sigmarg;

for trial_to_check = 1: numsubsets

    zestmat(:,:,trial_to_check) = Cmud(:, (trial_to_check-1)*T + 1: trial_to_check*T);
    
    z1_est(trial_to_check,:) = sum(Cmud(1:p/2,(trial_to_check-1)*T + 1: trial_to_check*T));
    z2_est(trial_to_check,:) = sum(Cmud(p/2+1:p,(trial_to_check-1)*T + 1: trial_to_check*T));
    
    corr_z_est(trial_to_check) = corr(z1_est(trial_to_check,:)', z2_est(trial_to_check,:)');
    
    invsig = invsigmat(:,:,(trial_to_check-1)*T + 1: trial_to_check*T);
    covz_errbar = zeros(p, p, T);
    errorbar = zeros(p, T);
    for t=1:T
        covz_errbar(:,:,t) = CC*(inv(invsig(:,:,t))+covhh)*CC';
        errorbar(:,t) =  diag(covz_errbar(:,:,t));
    end
    
    z1(trial_to_check,:) = sum(xyzinpn.z(1:p/2,(trial_to_check-1)*T + 1: trial_to_check*T));
    z2(trial_to_check,:) = sum(xyzinpn.z(p/2+1:p,(trial_to_check-1)*T + 1: trial_to_check*T));
    
    corr_z(trial_to_check) = corr(z1(trial_to_check,:)', z2(trial_to_check,:)'); 

end

figure
subplot(211);
plot(1:numsubsets, mean(z1,2)/(p/2),'r', 1:numsubsets, mean(z2,2)/(p/2), 'b', ...
    1:numsubsets, mean(z1_est,2)/(p/2), 'r--', 1:numsubsets, mean(z2_est,2)/(p/2), 'b--')
 legend('true z (grp1)', 'estimated z(grp1)', 'true z (grp2)', 'estimated z(grp2)');

% (2) covariances

numtimebins = T;
autocorr = zeros(p, numtimebins+1);

for whichcell = 1:p
     autocorr(whichcell, :) = xcov(zz(whichcell,:), numtimebins/2, 'unbiased');
end

avgautocorr_acrcells = mean(autocorr);

autocorr_condi = zeros(p, numtimebins+1, numsubsets);

for whichcell = 1:p
    
    for whichtrial = 1:numsubsets
        autocorr_condi(whichcell, :, whichtrial) = xcov(z(whichcell,:,whichtrial), numtimebins/2, 'unbiased');
    end
    
end

autocorr_per_eachtrial = squeeze(mean(autocorr_condi));
avgautocorr_condi = mean(autocorr_per_eachtrial,2);

figure;
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

autocorr_condi = zeros(p, numtimebins+1, numsubsets);

for whichcell = 1:p
    
    for whichtrial = 1:numsubsets
        autocorr_condi(whichcell, :, whichtrial) = xcov(zestmat(whichcell,:,whichtrial), numtimebins/2, 'unbiased');
    end
    
end

autocorr_per_eachtrial = squeeze(mean(autocorr_condi));
avgautocorr_condi = mean(autocorr_per_eachtrial,2);

subplot(212);
plot(1:numtimebins+1, avgautocorr_acrcells/max(avgautocorr_acrcells), 'k--', 1:numtimebins+1, avgautocorr_condi/max(avgautocorr_condi), 'r--')
legend('estimated total covariance', 'estimated conditional covariance');
