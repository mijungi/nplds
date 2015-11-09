%This demo shows the whole procedure of fitting our model to training data,
%then predicting the mean firing rates of neurons on left-out trials

clear all;
close all;
clc;

addpath core_functions/
addpath gpml-matlab/gpml
addpath Figure3_and_4/functions

% essential quantities (fixing pinp==0)
k = 2; % dimensionality of hidden states
p = 10; % output dimension (e.g., # of neurons)
d = 4; % input dimension
T = 100; % time within each recording 1:T
r = 10; % # of recordings

%% defining hyper-parameters and draw samples for parameters.

Model =  'NSFR';

% generate parameters 
tau_init = 10*rand; 
params = generate_params_from_prior(k, d, p, r, Model, 1:r, tau_init);

% generate inputs 
params.inpn = 0.5*randn(d, T);

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


%% Split data into training and validation (this shall be done multiple times for a proper analysis)
  ind_test = 6; %out of 10 trials;
  ind_train = 1:r; ind_train(6) = [];
  rtrain = length(ind_train);

  params.ind_train = ind_train;
  params.ind_test = ind_test;
  
 %% Fit the same data with 3 different latent dimensionalities
 pred_errors = zeros(3,1);
 h_all = cell(3,1);
 for k = 1:3
 % Initialize the fitting parameters based on the training data.
  fprintf('\n\nFitting k = %d latent dimensionality\n', k);
  params.maxIter = 3; %FITPARAM Maximum number of VBEM iterations
  
  params.A = 0.9*eye(k); %FITPARAM Inialially contractive dynamics matrix (0.9)
  params.C = randn(p,k); %neuron num by latent dim
  params.B = randn(k,d); %latent dim by input dim
  %params.inpn  is already set, in experimental settings we can assume we know the true input. 
  

  params.x0 = zeros(k,1);
  params.V0 = eye(k);
  params.Q = eye(k);

  newdataset = cell(rtrain,1);
  yy =[];

  for i=1:rtrain
      newdataset{i} = xyzinpn{ind_train(i)};
      yy = [yy newdataset{i}.y]; %concatante all training trials to get initial estimates for parameters
  end
  
  fprintf('Fitting GP to initialize tau.\n');
  params.tau_init = GPfittingTrainingData(yy);

  meanyy = mean(yy,2);
  params.d = log(meanyy);
  params.yy = yy; 
  
  tmp = reshape(params.yy,size(params.yy,1),length(params.ind_train),[]);
  tmp = mean(tmp,3);
  params.sig_init = 1.5*max(max(tmp,[],2) -min(tmp,[],2)); %FITPARAM Initial noise parameter in the kernel. 1.5 * maximum change in firing rate of a neuron over observed period
  params.m_h_init = zeros(k,1); % Initialize the expected mean of the latent h
  
  %% Run the fitting on the training data
    datastruct = VBEM_PLDSnonstationary(newdataset, rtrain, params, Model); 
    
    
  %% Make predictions for the left-out test trial
    fprintf('\nPredicting left-out trial for k = %d latent dimensionality\n', k);
    [fr_pred_error, hpred] = make_predictions_demo(xyzinpn, datastruct, params, k, Model);
   
    pred_errors(k) = fr_pred_error;
    h_all{k} = zeros(k,r);
    h_all{k}(:,params.ind_train) = datastruct.Mstep{end}.h;
    h_all{k}(:,params.ind_test) = hpred;
    
 end
    
 
 %% Analyze results
 
 %Plot fitted and predicted h
 for k=1:3
    figure(k+1);
    h1= plot(1:r, params.h,'black');
    hold on;
    h2 = plot(params.ind_train, h_all{k}(:, params.ind_train),'r');
    h3 = scatter(repmat(params.ind_test,k,1), h_all{k}(:, params.ind_test), 55, 'gx');
    legend([h1(1); h2(1); h3(1)], {'True h', 'Fitted h', 'Predicted h'})
    xlabel('Trial number');
    ylabel('Estimated h');
    title(['Estimating h with k = ' num2str(k) ' latent dimensionality']);
 end
 
 %Output RMSE values 
 fprintf('\n\n The estimated RMSE in predicting the mean firing rates on the left-out trial for different latent dimensionalities are: \n   k     Error  \n   1   %f4 \n   2   %f4 \n   3   %f4 \n',pred_errors(1), pred_errors(2), pred_errors(3))
 %Note that because all the data generation etc is random, this will change
 %from run to run. Also, due to the data provided here is insufficient to
 %properly regularize the model, as well as the validation is only for a
 %single trial instead of cross-validation, in many cases the true latent
 %dimensionality will not provide the least prediction error.