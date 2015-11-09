function fr_pred_error = make_predictions_all( orig_data_file, saved_data_file, saved_params_file, code_dir, k, Model, redo)
%MAKE_PREDICTIONS Simulates data not only for left-out trials but also for
%trials we fitted over, purely for visualization purposes, this is not used to
%assess performance


if nargin < 7
  redo = 0;
end



%Check if it was already done
if exist([saved_data_file(1:end-4) '_predfr_all.mat'],'file') && redo==0
  load([saved_data_file(1:end-4) '_predfr_all.mat'], 'fr_pred_error');
  return;
end

[status, seed] = system('od /dev/urandom --read-bytes=4 -tu | awk ''{print $2}''');
seed=str2double(seed);
rng(mod(seed, intmax));

%% draw samples for y, compute mean(y)

load(orig_data_file, 'resps', 'stims');
load(saved_data_file, 'datastruct')
load(saved_params_file, 'params')

cd(code_dir);
addpath core_functions
addpath gpml-matlab/gpml
addpath Figure3_and_4/functions


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

hall = zeros(k, size(resps,3));
hall(:,params.ind_train) = Mstepresults.h;
hall(:,params.ind_test) = hpred;


all_data = cell(1,size(resps,3));
nsamps = 1000; %for each left out trial, generate nsamps samples and then take the mean of the firing rates as the predicted mean from h
xyzinpn_tst = cell(nsamps,1);
fr_mean_pred_all = zeros(size(params.C,1), size(resps,3),nsamps); % num_cells x test indices x samples

%For each left out trial
for tst1 = 1:size(resps,3)
  tst1
  %generate the samples, compute mean
  Mstepresults.h = hall(:,tst1);
  Mstepresults.inpn = stims(:,:,tst1);
  for i1=1:nsamps
    xyzinpn_tst{i1} = generate_data_PLDS_NSFR(size(resps,2), Mstepresults);
    fr_mean_pred_all(:,tst1,i1) = mean(exp(xyzinpn_tst{i1}.z),2);
  end
  all_data{tst1} = xyzinpn_tst;
end

%Compute predicted mean firing rate for all trials and all cells 
fr_mean_pred = median(fr_mean_pred_all,3); % num_cells x trial indices


%% Get true mean firing rates from left out trials, then compute error

fr_mean_true = resps(:,:,:);
fr_mean_true = squeeze(mean(fr_mean_true,2));


fr_pred_error = sqrt(sum(sum((fr_mean_pred - fr_mean_true).^2))); %RMSE value

Mstepresults = datastruct.Mstep{end};
save([saved_data_file(1:end-4) '_predfr_all.mat'], 'fr_pred_error', 'fr_mean_pred_all', 'fr_mean_true', 'Mstepresults', 'hpred', 'hall','params', 'all_data','-v7.3');

end

