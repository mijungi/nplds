%% generate data for Figure 2
% mijung wrote on April 22, 2014
% edited on Oct 3, 2015


clear all;
close all;
clc;

seed = 1;
oldRng = rng();
rng(seed);

addpath ../core_functions

% essential quantities (fixing pinp==0)
k = 4; % dimensionality of hidden states
p = 40; % output dimension (e.g., # of neurons)
pinp = 2; % input dimension
T = 200; % time 1:T (for training)
r = 100; % # of recordings

%% generate params (A, B, C, h, d)

% (1) A
A = zeros(k); A(3,3) = 0.9; A(4,4) = 0.9;

% (2) B
B = zeros(k, pinp); B(1:2,1:2) = [5 0; 0 6]; 

% (3) C
C = zeros(p,k); 
maxC1 = 0.1;
maxC2 = 0.1;
C(1:p/2,1) = maxC1*ones(p/2,1); C(1:p/2,3) = maxC1*ones(p/2,1);
C(p/2+1:p,2) = maxC2*ones(p/2,1); C(p/2+1:p,4) = maxC2*ones(p/2,1);
 
C(:,3) = C(:,3)/3;
C(:,4) = C(:,4)/3;

% (4) d
d = [-2.2*ones(p/2,1); -1.2*ones(p/2,1)];

% (5) h
h = zeros(k,r);
h13 = normpdf(1:r, r/4, r/8);
h13 = 12*h13/max(h13);
h24 = - fliplr(h13);
h(1,:) = h13; h(3,:) = h13;
h(2,:) = h24; h(4,:) = h24; 

% (6) priors on x and noise variance for x
x0 = zeros(k, 1); V0 = eye(k); sqQ = sqrtm(eye(k));

% generate inputs
twoperiods = 100;
% inp - u by Tn input sequence
ti = [1:T]/twoperiods*4*pi;
window = [zeros(1,T/4) ones(1,T/2) zeros(1,T/4)];
inpn = 0.4*[cos(ti).*window; sin(ti).*window]; % assumes u = 2;

% putting parameters to params structure
params.A = A;
params.B = B;
params.C = C;
params.d = d;
params.Q = sqQ;
params.h = h;
params.inpn = inpn;
params.x0 = x0;
params.V0 = V0;

Model = 'NSFR';
xyzinpn = generate_data_PLDS_multiple_recordings(T,params,r,Model);

meanspikecount = zeros(p, r);
meanz = zeros(p,r);
meanx = zeros(k,r);
y = zeros(p, T, r);
z = zeros(p, T, r);

for i=1:r
    y(:,:,i) = xyzinpn{i}.y;  
    x = xyzinpn{i}.x;
    meanspikecount(:,i) = mean(xyzinpn{i}.y,2);
    meanz(:,i) = mean(xyzinpn{i}.z,2);
    meanx(:,i) = mean(x,2);
    z(:,:,i) = xyzinpn{i}.z;
end

save mat_files/all_NSFR.mat;

%% visualise log mean firing rate of generated data (Figure2: A and B)

mean_zz_g1 = zeros(r,1);
mean_zz_g2 = zeros(r,1);
for i=1:r
    mean_zz_g1(i,:) = mean(mean(z(1:p/2,:,i)));
    mean_zz_g2(i,:) = mean(mean(z(p/2+1:p, :,i)));
end

figure; 
subplot(211);
plot(1:r, mean_zz_g1,'r')
subplot(212);
plot(1:r, mean_zz_g2, 'b')

%% raster plot (Figure 2: C and D) 

binsize = 50; % 50ms bin, which makes each recording 10s long.

rasterplot = zeros(p, binsize*T, r); 
    
for trial = 1:r   
    for whichneuron = 1: p      
        for t=1:T
 
            if y(whichneuron,t,trial)>0
                howmanyspikes = y(whichneuron,t,trial);
                randomspikelocations = rand(howmanyspikes,1);
                spikelocations_ineachbin = 1 + (t-1)*binsize + floor(binsize.*randomspikelocations);               
                % how many spikes?
                rasterplot(whichneuron, spikelocations_ineachbin, trial) = 1;                
            end
        end
    end
end

figure; subplot(2,3,[ 4 5 6]);

hold on;
% Plot spikes
height = 1 / p;
trial_to_plot = 25; 

for i=1:p
    for t = 1:binsize*T
        if rasterplot(i, t, trial_to_plot)
            plot([t, t], [(i - 1) * height, i * height], 'k', 'LineWidth', 1);  ylabel('neurons'); xlabel('time');
        end
    end
end

hold off

%% visualise covariances (total, and conditional. i.e., Figure 2 E)

% total covariance
yy = [];
zz = [];
for i=1:r
    yy = [yy y(:,:,i)];
    zz = [zz z(:,:,i)];
end

% conditional covariance
condicov = zeros(p,p,r);
for i=1:r
    condicov(:,:,i) = cov(y(:,:,i)',1);
end

meanmat = zeros(p,r);
for i=1:r
    meanmat(:,i) = mean(y(:,:,i),2);
end
condicov2 = cov(meanmat',1);

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
plot(1:numtimebins+1, avgautocorr_acrcells/max(avgautocorr_acrcells), 'k', 1:numtimebins+1, avgautocorr_condi/max(avgautocorr_condi), 'r')
legend('total covariance', 'conditional covariance');
