function [fr_pred_error, hpred] = make_predictions_demo( data, datastruct, params, k, Model)
%MAKE_PREDICTIONS Used just for the demo (not using the file structure as input).
% After the fitting has been done, use the converged model
%parameters to make predictions on left-out trials, and evaluate the
%performance.


addpath core_functions
addpath gpml-matlab/gpml
addpath Figure3_and_4/functions

% 
%% predict h for held-out datasets


Estepresults = datastruct.Estep{end};
Mstepresults = datastruct.Mstep{end};

if strcmp(Model, 'NSFR')
  epsilon = 1e-3;

  Kstar = makeK_test(Mstepresults.sig2, epsilon, Mstepresults.tau2, length(params.ind_train), k, length(params.ind_test), params.ind_test, params.ind_train);
  K = makeK_test(Mstepresults.sig2, epsilon, Mstepresults.tau2, length(params.ind_train), k, length(params.ind_train), params.ind_train, params.ind_train);

  hmu = reshape(Mstepresults.h,[],1);

  hpred = (Kstar/K)*hmu;

  hpred = reshape(hpred, k, length(params.ind_test));
elseif strcmp(Model, 'PLDS')
  hpred = zeros(k,length(params.ind_test));
end


nsamps = 1000; %for each left out trial, generate nsamps samples and then take the mean of the firing rates as the predicted mean from h
xyzinpn_tst = cell(nsamps,1);
fr_mean_pred_all = zeros(size(params.C,1), length(params.ind_test),nsamps); % num_cells x test indices x samples

%For each left out trial
for tst1 = 1:numel(params.ind_test)
  %generate the samples, compute mean
  Mstepresults.h = hpred(:,tst1);
  Mstepresults.inpn = params.inpn;
  for i1=1:nsamps
    if mod(i1,100)==0, fprintf('Drawing samples for prediction, progress %d per cent\n', floor(i1/10)); end;
    xyzinpn_tst{i1} = generate_data_PLDS_NSFR(size(data{1}.x,2), Mstepresults);
    fr_mean_pred_all(:,tst1,i1) = mean(exp(xyzinpn_tst{i1}.z),2);
  end
end

%Compute predicted mean firing rate for all trials and all cells 
fr_mean_pred = median(fr_mean_pred_all,3); % num_cells x test indices


%% Get true mean firing rates from left out trials, then compute error
fr_mean_true = [];
for tst1 = 1:numel(params.ind_test)
    fr_mean_true(:,:,tst1) = data{tst1}.y; %neuron by time by trials
end
fr_mean_true = squeeze(mean(fr_mean_true,2)); %mean firing rate, neurons by trials


fr_pred_error = sqrt(sum(sum((fr_mean_pred - fr_mean_true).^2))); %RMSE value over all neurons and test trials


end










