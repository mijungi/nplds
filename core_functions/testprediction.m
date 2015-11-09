% test prediction


clear all;
close all;
clc;

% add path to use 'hinton' figure to show estimates
addpath ../mattBealsCode_v3_4_1/
addpath ../standardEM/
addpath ../testcode_nonstationaryPLDS_A_multipleRecordings/

%%
% load data
load all_NSFR.mat

% r = 100;
% leave every 5th trial for testing
ind_test = 10:10:r;
rtest = length(ind_test);
ind_train = find(rem([1:r],10)); 
rtrain = length(ind_train);

params.ind_train = ind_train;

% put all the true params
params.A = A;
params.B = B;
params.C = C;
params.d = d;
params.h = h(:,ind_train);

params.x0 = x0;
params.V0 = V0;
params.Q = eye(k);

newdataset = cell(rtrain,1);

for i=1:rtrain
    newdataset{i}.x = xyzinpn{ind_train(i)}.x;
    newdataset{i}.y = xyzinpn{ind_train(i)}.y;
    newdataset{i}.z = xyzinpn{ind_train(i)}.z;
    newdataset{i}.T = T;
    newdataset{i}.inpn = xyzinpn{ind_train(i)}.inp;
end

% xyzinpn = newdataset;

%% start fitting here

% Model = 'NSFR';
% datastruct = VBEM_PLDSnonstationary(xyzinpn, rtrain, params, Model); 

% save datastruct_testprediction_1to100_10percentTesting datastruct
% save datastruct_testprediction_11to30 datastruct
% 

% load  datastruct_testprediction_11to30

%% sanity check

% load datastruct_testprediction_1to100;
load datastruct_testprediction_1to100_10percentTesting;

z1_est = zeros(rtrain, T);
z2_est = zeros(rtrain, T);

z1 = zeros(rtrain, T);
z2 = zeros(rtrain, T);

corr_z_est = zeros(rtrain,1);
corr_z = zeros(rtrain, 1);

fromMstep = datastruct.Mstep{20};
fromEstep = datastruct.Estep{20};

for trial_to_check = 1:rtrain
    
    CC = fromMstep.C;
    hh = fromMstep.h(:, trial_to_check);
    covhh = fromMstep.covh(:,:,trial_to_check);
    mu = fromEstep{trial_to_check}.mumarg;
    
    Cmud = zeros(p, T);
    for t=1:T
        Cmud(:,t) = CC*(mu(:,t)+hh) + fromMstep.d;
    end
    
    z1_est(trial_to_check,:) = sum(Cmud(1:p/2,:));
    z2_est(trial_to_check,:) = sum(Cmud(p/2+1:p,:));
    
    corr_z_est(trial_to_check) = corr(z1_est(trial_to_check,:)', z2_est(trial_to_check,:)');
    
    invsig = fromEstep{trial_to_check}.inv_sigmarg;
    covz_errbar = zeros(p, p, T);
    errorbar = zeros(p, T);
    for t=1:T
        covz_errbar(:,:,t) = CC*(inv(invsig(:,:,t))+covhh)*CC';
        errorbar(:,t) =  diag(covz_errbar(:,:,t));
    end
    
    
    % plot(1:T, k1k2mat(2,:,1)', 'k', 1:T, k1k2mat_est(2,:,1), 'r', 1:T, k1k2mat_est(2,:,1)-1.64*sqrt(sum(errorbar(1+p/2:p,:))), 'r--', 1:T, k1k2mat_est(2,:,1)+1.64*sqrt(sum(errorbar(1+p/2:p,:))), 'r--');
    % set(gca, 'ylim', [-45 -35])
    
    z1(trial_to_check,:) = sum(xyzinpn{trial_to_check}.z(1:p/2,:,1));
    z2(trial_to_check,:) = sum(xyzinpn{trial_to_check}.z(p/2+1:p,:,1));
    
    corr_z(trial_to_check) = corr(z1(trial_to_check,:)', z2(trial_to_check,:)'); 

    
end

subplot(311); plot(1:rtrain, corr_z, 'o-', 1:rtrain, corr_z_est, 'o-'); 
% set(gca, 'ylim', [0 1]);
subplot(312); plot(1:T, sum(z1), 'k', 1:T, sum(z1_est), 'r')
subplot(313); plot(1:T, sum(z2), 'k', 1:T, sum(z2_est), 'r')

%% predict h for held-out datasets

epsilon = 1e-3;

Kstar = makeK_test(fromMstep.sig2, epsilon, fromMstep.tau2, rtrain, k, rtest, ind_test, ind_train);
K = makeK_test(fromMstep.sig2, epsilon, fromMstep.tau2, rtrain, k, rtrain, ind_train, ind_train);

hmu = reshape(fromMstep.h,[],1);
% hmu = zeros(k2*rtrain,1);
% for i=1:rtrain
%     hmu(k2*(i-1)+1:i*k2) = reshape(fromMstep.h(:,i)', [],1);
% end

hpred = (Kstar/K)*hmu;

hpred = reshape(hpred, k, rtest);

imagesc(K)


hall = zeros(k, r);
hall(:,ind_train) = fromMstep.h;
hall(:,ind_test) = hpred;

%% compute z on test data

z1_est_tst = zeros(rtest, T);
z2_est_tst = zeros(rtest, T);

z1_tst = zeros(rtest, T);
z2_tst = zeros(rtest, T);

Cest = fromMstep.C;
Aest = fromMstep.A;
Best = fromMstep.B;
x0est = fromMstep.x0;
V0est = fromMstep.V0;
inp = xyzinpn{1}.inp;

for trial_to_check = 1:rtest
    
    hest = hpred(:, trial_to_check);
    
    Cmud = zeros(p, T);
    x = zeros(k,T);
    for t=1:T
        if t==1
            x(:,t) = x0est + sqQ*randn(k,1);
        else           
            x(:,t) = Aest*x(:,t-1)+ Best*inp(:,t)+ sqQ*randn(k,1);
            Cmud(:,t) = Cest*(x(:,t)+hest) + fromMstep.d;
        end
    end
    
    z1_est_tst(trial_to_check,:) = sum(Cmud(1:p/2,:));
    z2_est_tst(trial_to_check,:) = sum(Cmud(p/2+1:p,:));
    
    z1_tst(trial_to_check,:) = sum(z(1:p/2,:,ind_test(trial_to_check)));
    z2_tst(trial_to_check,:) = sum(z(p/2+1:p,:,ind_test(trial_to_check)));
    
%     corr_z(trial_to_check) = corr(z1(trial_to_check,:)', z2(trial_to_check,:)'); 

    
end

subplot(211); plot(1:T, sum(z1_tst), 'k', 1:T, sum(z1_est_tst), 'r')
subplot(212); plot(1:T, sum(z2_tst), 'k', 1:T, sum(z2_est_tst), 'r')


%% generate data form hpred

Model = 'NSFR';
params = fromMstep;

hall = zeros(k, r);
hall(:,ind_train) = params.h;
hall(:,ind_test) = hpred;

params.h = hall;
params.inpn = xyzinpn{1}.inp;

nsamps = 1000;
xyzinpn_tst = cell(nsamps,1);
yytst = zeros(p, r*T, nsamps);
zztst = zeros(p, r*T, nsamps);

for i=1:nsamps
    [i nsamps]
    xyzinpn_tst{i} = generate_data_PLDS_multiple_recordings(T,params,r,Model);
    
    for wt = 1:r
        yytst(:,(wt-1)*T+1:wt*T,i) = xyzinpn_tst{i}{wt}.y;
        zztst(:,(wt-1)*T+1:wt*T,i) = xyzinpn_tst{i}{wt}.z;
    end
%     yytst(:,:,i) = [xyzinpn_tst{i}{1}.y xyzinpn_tst{i}{2}.y xyzinpn_tst{i}{3}.y xyzinpn_tst{i}{4}.y];
%     zztst(:,:,i) = [xyzinpn_tst{i}{1}.z xyzinpn_tst{i}{2}.z xyzinpn_tst{i}{3}.z xyzinpn_tst{i}{4}.z];
end

%%
% samplemean = mean(yytst, 3);
% samplemeanzz = mean(zztst, 3);
% 
% whichtrial  = 10; 
% % mean([mean(samplemean(1:p/2,(whichtrial-1)*T+1:whichtrial*T))' mean(y(1:p/2,:,ind_test(whichtrial)))'])
% mean([mean(samplemeanzz(1:p/2,(whichtrial-1)*T+1:whichtrial*T))' mean(z(1:p/2,:,ind_test(whichtrial)))'])
% 
% % mean([mean(samplemean(p/2+1:p,(whichtrial-1)*T+1:whichtrial*T))' mean(y(p/2+1:p,:,ind_test(whichtrial)))'])
% mean([mean(samplemeanzz(p/2+1:p,(whichtrial-1)*T+1:whichtrial*T))' mean(z(p/2+1:p,:,ind_test(whichtrial)))'])


%%

samplemean = mean(zztst, 3);
samplemeanyy = mean(yytst, 3);
% 
% ytruetst = z(:,:,ind_test);

mfr_truezz = zeros(p, r);
mfr_estzz = zeros(p, r);
mfr_trueyy = zeros(p, r);
mfr_estyy = zeros(p, r);

for whichtrial  = 1: r
    mfr_truezz(:,whichtrial) = mean(z(:,:,whichtrial),2);
    mfr_estzz(:, whichtrial) = mean(samplemean(:,(whichtrial-1)*T+1:whichtrial*T),2); 
    mfr_trueyy(:,whichtrial) = mean(y(:,:,whichtrial),2);
    mfr_estyy(:, whichtrial) = mean(samplemeanyy(:,(whichtrial-1)*T+1:whichtrial*T),2);
end

figure;
subplot(211); plot(1:r, mfr_truezz', 'k', 1:r, mfr_estzz', 'r');
subplot(212); plot(1:r, mfr_trueyy', 'k', 1:r, mfr_estyy', 'r');

%%

subplot(211); plot(1:rtest, mfr_truezz(:,ind_test)', 'k', 1:rtest, mfr_estzz(:,ind_test)', 'r');
subplot(212); plot(1:rtest, mfr_trueyy(:,ind_test)', 'k', 1:rtest, mfr_estyy(:,ind_test)', 'r');
