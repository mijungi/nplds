ls
%% testscript for PLDS with nonstationary firing rates
% testscriptPLDS_NSFR
% May 8, 2014
% wrote by Mijung 

clear all;
close all;
clc;

% add path to use 'hinton' figure to show estimates
addpath ../mattBealsCode_v3_4_1/
addpath ../standardEM/
addpath ../testcode_nonstationaryPLDS_A_multipleRecordings/

% essential quantities (fixing pinp==0)
k = 2; % dimensionality of hidden states
p = 10; % output dimension (e.g., # of neurons)
d = 4; % input dimension
T = 100; % time 1:T (for training)
r = 10; % # of recordings

%% defining hyper-parameters and draw samples for parameters.

Model =  'NSFR';
params = generate_params_from_prior(k, d, p, r, Model);

% params.h = zeros(size(params.h));
% params.B = zeros(size(params.B));
% params.AB = zeros(size(params.AB));

% generate inputs and observations

if strcmp(Model, 'NSFR') || strcmp(Model, 'PLDS')
    params.inpn = 0.5*randn(d, T);
else % NSLS
    params.inpn = [];
end

% we assume T is same for each recording for now
xyzinpn = generate_data_PLDS_multiple_recordings(T, params, r, Model);

meanfiringrate = zeros(p,r);
for i=1:r
    meanfiringrate(:,i) = mean(xyzinpn{i}.y,2);
end

figure(1);
subplot(2,4,[1 2 3 4]); plot(meanfiringrate, 'o-'); title('mean firing rate');
% subplot(2,4,[1 2]); plot([xyzinpn{1}.y xyzinpn{2}.y xyzinpn{3}.y]', 'o'); set(gca, 'xlim', [0 r*T]);  title('spike counts');
% subplot(2,4,[3 4]); plot([xyzinpn{1}.x xyzinpn{2}.x xyzinpn{3}.x]'); set(gca, 'xlim', [0 r*T]); title('true x');
subplot(2,4,5); hinton(params.A,['true A1  ' num2str(max(max(params.A)),'%.4f')],'standard');
subplot(2,4,6); hinton(params.B,['true B  ' num2str(max(max(params.B)),'%.4f')],'standard');
subplot(2,4,7); hinton(params.C,['true C  ' num2str(max(max(params.C)),'%.4f')],'standard');
subplot(2,4,8); plot(params.h');

%% VBEM

datastruct = VBEM_PLDSnonstationary(xyzinpn, r, params, Model); 


%% sanity check

% nsamps = 500;
% 
% % assuming T is same for all recordings
% Zn = zeros(p, r*T, nsamps);
% Yn = zeros(p, r*T, nsamps);
% Yn_estParam = zeros(p, r*T, nsamps);
% Zn_estParam = zeros(p, r*T, nsamps);
% 
% datastruct.Mstep.inpn = params.inpn;
% 
% for i=1:nsamps
%     [i nsamps]
%     sanitycheck = generate_data_PLDS_multiple_recordings(T,datastruct.Mstep,r,Model);
%     sanitychecktrue = generate_data_PLDS_multiple_recordings(T,params,r,Model);
%     
%     for rr = 1: r
%         Zn_estParam(:, (rr-1)*T+1:rr*T, i) = sanitycheck{rr}.z;
%         Yn_estParam(:, (rr-1)*T+1:rr*T, i) = sanitycheck{rr}.y;
%         Zn(:,(rr-1)*T+1:rr*T,i) = sanitychecktrue{rr}.z;
%         Yn(:,(rr-1)*T+1:rr*T,i) = sanitychecktrue{rr}.y;
%     end
% end

%%

% % mean vector
% meanvec = zeros(p,r);
% meanvec_est = zeros(p,r);
% 
% covmat = zeros(p, p, r);
% covmat_est = zeros(p, p, r);
% 
% for rr = 1: r
% 
%     %% test mean
%     meanvec(:,rr) = mean(mean(Zn(:,(rr-1)*T+1:rr*T,:),3),2);
%     meanvec_est(:,rr) = mean(mean(Zn_estParam(:, (rr-1)*T+1:rr*T, :),3),2);
%     
%     %% test covariance
%     
%     Znrshp = reshape(Zn(:,(rr-1)*T+1:rr*T,:), p, []);
%     Znestrshp = reshape(Zn_estParam(:,(rr-1)*T+1:rr*T,:), p, []);
%     
%     covmat(:,:,rr) = cov(Znrshp');
%     covmat_est(:,:,rr) = cov(Znestrshp');
%     
%     % onestep delayed covariance
%     Znrshp2 = Znrshp(:,[2:end,1]);
%     Znestrshp2 = Znestrshp(:, [2:end,1]);
%         
%     %% figure
%     figure;
%     
%     subplot(2,2,[1 2]); plot([meanvec(:,rr) meanvec_est(:,rr)]); legend('mean of Z using true params', 'mean of Z using estimated params'); xlabel('neurons');
% %     set(gca,'ylim', [params.d(1)-0.2, params.d(1)+0.2]);
%     
%     subplot(2,2,3); hinton(covmat(:,:,rr),['cov z' num2str(max(max(covmat(:,:,rr))),'%.4f')],'standard');
%     subplot(2,2,4); hinton(covmat_est(:,:,rr),['cov z est  ' num2str(max(max(covmat_est(:,:,rr))),'%.4f')],'standard');
%         
% end

%%

z1_est = zeros(r, T);
z2_est = zeros(r, T);

z1 = zeros(r, T);
z2 = zeros(r, T);

corr_z_est = zeros(r,1);
corr_z = zeros(r, 1);

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
    
    corr_z_est(trial_to_check) = corr(z1_est(trial_to_check,:)', z2_est(trial_to_check,:)');
    
    invsig = datastruct.Estep{trial_to_check}.inv_sigmarg;
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

subplot(311); plot(1:r, corr_z, 'o-', 1:r, corr_z_est, 'o-'); set(gca, 'ylim', [0 1]);
subplot(312); plot(1:T, sum(z1), 'k', 1:T, sum(z1_est), 'r')
subplot(313); plot(1:T, sum(z2), 'k', 1:T, sum(z2_est), 'r')
% subplot(212); plot(1:r, squeeze(Q(1,2,:))./squeeze(Q(1,1,:)), 'k', 1:r, corr_z_est, 'ro-')
%%
% plot([z1' z2']);

% figure;
alpha = 1.96;
% alpha =10;
subplot(211); plot(1:T, sum(z1), 'k', 1:T, sum(z1_est), 'r', 1:T, sum(z1_est)-alpha*sum(sqrt(errorbar(1:p/2,:))), 'r--', 1:T, sum(z1_est)+alpha*sum(sqrt(errorbar(1:p/2,:))), 'r--')
subplot(212); plot(1:T, sum(z2), 'k', 1:T, sum(z2_est), 'r', 1:T, sum(z2_est)-alpha*sum(sqrt(errorbar(1+p/2:p,:))), 'r--', 1:T, sum(z2_est)+alpha*sum(sqrt(errorbar(1+p/2:p,:))), 'r--')
