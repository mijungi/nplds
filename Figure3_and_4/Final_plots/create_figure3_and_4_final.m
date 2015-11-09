%% Global settings
%HARDCODED_DIRECTORY_PATHS
code_dir = '/nfs/nhome/live/gbohner/Git_main/nonstat_plds/nplds_code_NIPS2015'; % Github main code directory
data_file = '/nfs/data3/gergo/Mijung/Figure4/Input_data/alexdata_session2_org0_bin50.mat'; % Data file used
NSFR_fit_file = '/nfs/data3/gergo/Mijung/Figure4/Example_run/k7_run5/final_NSFR_data.mat'; % Fitted data for example run
NSFR_sim_file = '/nfs/data3/gergo/Mijung/Figure4/Example_run/k7_run5/final_NSFR_data_predfr_all.mat'; % Predicted data for example run
PLDS_sim_file = '/nfs/data3/gergo/Mijung/Figure4/Example_run/k7_run5/final_PLDS_data_predfr_all.mat'; % Predicted data for example run

cd(code_dir);

DATA_COLOR = [0 0 0];
NSFR_COLOR = [1 0 0];
PLDS_COLOR = [.4 .4 .4];

%% Panel 3A - Example data, plus NSFR and PLDS fits

% % Extract the data (commented out)
% load(data_file, 'resps');
% resps = double(resps);
% resps = squeeze(mean(resps,2));
% pop_mean = mean(mean(mean(resps)))*20; %Average spike / bin * 20 (due to 50 ms bins) = average firing rate;
% resps_orig = resps;
% % 
% %Do smoothing then least square to find non-stationarity
% for i1 = 1:size(resps,1)
%   resps(i1,:) = smooth(resps(i1,:),'moving',15);
% end
% resps_smooth = resps;
% resps = bsxfun(@minus, resps, mean(resps,2));
% resps = sqrt(mean(resps.^2,2));
% [~, ix] = sort(resps,'descend');
% 
% % %Do linear regression to find non-stationarity
% % lincoeffs = zeros(size(resps,1),1);
% % for i1 = 1:size(resps,1)
% %   p = polyfit(1:100,resps(i1,:),1);
% %   lincoeffs(i1) = p(1);
% % end
% % [~,ix] = sort(abs(lincoeffs),'descend');
%  
% %Choose which neurons to show
% neurons_to_show_l = ix(1:5); %5 most nonstationary
% neurons_to_show_r = ix(end-6:end-2); %5 most stationary (but still has spikes, last 2 does not)
% 
% %Load simulated data
% load(NSFR_sim_file, 'fr_mean_pred_all');
% fr_mean_pred_NSFR = median(fr_mean_pred_all,3);
% load(PLDS_sim_file, 'fr_mean_pred_all');
% fr_mean_pred_PLDS = median(fr_mean_pred_all,3);
% 
% save('Figure3_and_4/Final_plots/Panel_3A','resps_orig','fr_mean_pred_NSFR','fr_mean_pred_PLDS','neurons_to_show_l', 'neurons_to_show_r')

%Make the figure;
load('Figure3_and_4/Final_plots/Panel_3A');

% Show most nonstationary neurons
fig1 = figure(1); clf
hold on;
for i1 = 1:length(neurons_to_show_l)
  plot(resps_orig(neurons_to_show_l(i1),:)'*20+(i1-1)*20, 'Color', DATA_COLOR);
  line([0,100],[(i1-1)*20,(i1-1)*20],'Color','k');
end

for i1 = 1:length(neurons_to_show_l)
  plot(fr_mean_pred_NSFR(neurons_to_show_l(i1),:)'*20+(i1-1)*20,'Color', NSFR_COLOR);
end
 

for i1 = 1:length(neurons_to_show_l)
  plot(fr_mean_pred_PLDS(neurons_to_show_l(i1),:)'*20+(i1-1)*20,'Color',PLDS_COLOR);
end

set(gca,'XTick',0:25:100);
ylim([0,20*length(neurons_to_show_l)]);
set(gca,'YTick',0:5:(20*length(neurons_to_show_l)-1));
set(gca,'YTickLabel',repmat(0:5:15,1,length(neurons_to_show_l)));
xlabel('Trial')
ylabel('Firing rate (Hz)');

set(fig1,'Units','centimeters');
set(fig1,'Position',[1,1,12,24]);

print('-depsc2', 'Figure3_and_4/Final_plots/Panel_3Al.eps');

% Show most stationary neurons
fig11 = figure(11); clf
hold on;
for i1 = 1:length(neurons_to_show_r)
  plot(resps_orig(neurons_to_show_r(i1),:)'*20+(i1-1)*1, 'Color', DATA_COLOR);
  line([0,100],[(i1-1)*20,(i1-1)*20],'Color','k');
end

for i1 = 1:length(neurons_to_show_r)
  plot(fr_mean_pred_NSFR(neurons_to_show_r(i1),:)'*20+(i1-1)*1,'Color', NSFR_COLOR);
end
 

for i1 = 1:length(neurons_to_show_r)
  plot(fr_mean_pred_PLDS(neurons_to_show_r(i1),:)'*20+(i1-1)*1,'Color',PLDS_COLOR);
end

set(gca,'XTick',0:25:100);
ylim([0,1*length(neurons_to_show_r)]);
set(gca,'YTick',0:0.5:(1*length(neurons_to_show_r)-1));
set(gca,'YTickLabel',repmat(0:0.5:0.5,1,length(neurons_to_show_r)));
xlabel('Trial')
ylabel('Firing rate (Hz)');

set(fig11,'Units','centimeters');
set(fig11,'Position',[13,1,12,24]);
print('-depsc2', 'Figure3_and_4/Final_plots/Panel_3Ar.eps');

%% Panel 4B - Example modulation over trials

% Extract the data (commented out)
%See h-predictions (they're on the line for h);
% load(NSFR_fit_file, 'datastruct');
% load(NSFR_sim_file, 'hall'); % To compute h for all trials you need the predictions for the mean of h via the GP as well as the fitted ones.
% fig2 = figure(2);
% 
% C = datastruct.Mstep{end}.C;

% save('Figure3_and_4/Final_plots/Panel_4B','hall','C')

% Make the figure;
load('Figure3_and_4/Final_plots/Panel_4B');

fig2 = figure(2);
plot(1:100, hall); % Quite meaningless in terms of neural firing rates, see other options below
ylabel('Latent modulation');

% plot(1:100, exp(C*hall)*20); % AFter projecting through C matrix and exponentiating (individual neurons), *20 due to 50 ms bins to get Hz
% ylabel('Firing rate modulation (Hz)');

% plot(1:100, (C*hall) + log(20)) % AFter projecting through C matrix, but before exponentiating (individual neurons)
% ylabel('Log firing rate modulation (log Hz)');

xlabel('Trial')
set(fig2,'Units','centimeters');
set(fig2,'Position',[1,1,20,10]);
print('-depsc2', 'Figure3_and_4/Final_plots/Panel_4B.eps');

%% Panel 4A - Tau histogram
load('Figure3_and_4/Final_plots/Panel_4A','tau2_end'); % k x runs tau^2 after last iteration
fig3 = figure(3);
tau2_end = tau2_end(:);
hist(sqrt(tau2_end),0:1:20); %Create histogram
xlim([0,20])
ylim([0,30]);
hold on;
mean_tau = mean(sqrt(tau2_end));
line([mean_tau mean_tau], [0,25],'Color','r'); %Add line to the mean

xlabel('tau (trials)');
ylabel('Counts');
title('Tau is conserved over different latent dimensionalities and training sets');
set(fig3,'Units','centimeters');
set(fig3,'Position',[1,1,20,10]);
print('-depsc2', 'Figure3_and_4/Final_plots/Panel_4A.eps');

%% Panel 4C - Total Covariance
load('Figure3_and_4/Final_plots/Panel_4C');

%All XCov*_total are cell by timeshifts
%All XCov*_condi are trial by cell by timeshifts

%Compute normalized covariances
to_plot(1,:) = squeeze(mean(XCov_data_total,1))/max(mean(XCov_data_total,1));
to_plot(2,:) = squeeze(mean(XCov_NSFR_total,1))/max(mean(XCov_NSFR_total,1));
to_plot(3,:) = squeeze(mean(XCov_PLDS_total,1))/max(mean(XCov_PLDS_total,1));
to_plot(4,:) = squeeze(mean(mean(XCov_data_condi,1),2))/max(mean(mean(XCov_data_condi,1),2));
to_plot(5,:) = squeeze(mean(mean(XCov_NSFR_condi,1),2))/max(mean(mean(XCov_NSFR_condi,1),2));
to_plot(6,:) = squeeze(mean(mean(XCov_PLDS_condi,1),2))/max(mean(mean(XCov_PLDS_condi,1),2));


fig4 = figure(4); hold off;
%Total covariances
plot(-500:50:500,to_plot(1,:),'Color',DATA_COLOR);
hold on;
plot(-500:50:500,to_plot(2,:),'Color',NSFR_COLOR)
plot(-500:50:500,to_plot(3,:),'Color',PLDS_COLOR)

% %Conditional covariances
% plot(-500:50:500,to_plot(4,:),'Color',DATA_COLOR,'LineStyle','--'); 
% plot(-500:50:500,to_plot(5,:),'Color',NSFR_COLOR,'LineStyle','--')
% plot(-500:50:500,to_plot(6,:),'Color',PLDS_COLOR,'LineStyle','--')

legend({'Data','N-PLDS (k==7)','PLDS (k==7)'})
xlabel('Time lag (ms)')
ylabel('Normalized autocovariance')
set(fig4,'Units','centimeters');
set(fig4,'Position',[1,1,15,10]);
print('-depsc2', 'Figure3_and_4/Final_plots/Panel_4C.eps');

%% Panel 3Br - Error values
load('Figure3_and_4/Final_plots/Panel_3Br.mat')
fig5 = figure(5); 
errorbar((1:8)-.05, mean(NSFR_error,2)/640*20, std(NSFR_error,0,2)/640*20, 'rx'); %/640 due to 10 left out trials and 64 neurons/trials to get error / single bin, *20 to get back to Hz from 50 ms bins.
hold on; errorbar((1:8)+.05, mean(PLDS_error,2)/640*20, std(PLDS_error,0,2)/640*20, 'bx');
ylim([0.002 0.005]*20)
xlabel('k')
title('RMSE in firing rates on left out trials (10%) after 30 iter')
set(gca,'XTick',1:8);
ylabel('RMSE'); %Population mean firing rate 4.5412 over everything;
legend('N-PLDS','PLDS')
set(fig5,'Units','centimeters');
set(fig5,'Position',[1,1,20,10]);
print('-depsc2', 'Figure3_and_4/Final_plots/Panel_3Br.eps');

%% Panel F - Error values on just certain neurons, based on non-stationarity
load('Figure3_and_4/Final_plots/Panel_3Blm','NSFR_error_each', 'PLDS_error_each','neuron_stationarity_sort')

%Most non-stationary 5 neurons
neurons_to_include = neuron_stationarity_sort(1:5); %1 - most nonstationary, end - most stationary
NSFR_error = mean(NSFR_error_each(:,:,neurons_to_include),3);
PLDS_error = mean(PLDS_error_each(:,:,neurons_to_include),3);

fig6 = figure(6); 
errorbar((1:8)-.05, mean(NSFR_error,2)/(10*length(neurons_to_include))*20, std(NSFR_error,0,2)/(10*length(neurons_to_include))*20, 'rx'); %/(10*length(neurons_to_include)) due to 10 left out trials and certain number of chosen neurons/trials to get error / single bin, *20 to get back to Hz from 50 ms bins.
hold on; errorbar((1:8)+.05, mean(PLDS_error,2)/(10*length(neurons_to_include))*20, std(PLDS_error,0,2)/(10*length(neurons_to_include))*20, 'bx');
xlabel('k')
title('RMSE in firing rates of most non-stationary 5 neurons')
set(gca,'XTick',1:8);
ylabel('RMSE'); %Population mean firing rate 4.5412 over everything;
legend('N-PLDS','PLDS')
set(fig6,'Units','centimeters');
set(fig6,'Position',[1,1,20,10]);
print('-depsc2', 'Figure3_and_4/Final_plots/Panel_3Bl.eps');

%Most stationary 5 neurons
neurons_to_include = neuron_stationarity_sort(end-4:end); %1 - most nonstationary, end - most stationary
NSFR_error = mean(NSFR_error_each(:,:,neurons_to_include),3);
PLDS_error = mean(PLDS_error_each(:,:,neurons_to_include),3);

fig7 = figure(7); 
errorbar((1:8)-.05, mean(NSFR_error,2)/(10*length(neurons_to_include))*20, std(NSFR_error,0,2)/(10*length(neurons_to_include))*20, 'rx'); %/(10*length(neurons_to_include)) due to 10 left out trials and certain number of chosen neurons/trials to get error / single bin, *20 to get back to Hz from 50 ms bins.
hold on; errorbar((1:8)+.05, mean(PLDS_error,2)/(10*length(neurons_to_include))*20, std(PLDS_error,0,2)/(10*length(neurons_to_include))*20, 'bx');
xlabel('k')
title('RMSE in firing rates of most non-stationary 5 neurons')
set(gca,'XTick',1:8);
ylabel('RMSE'); %Population mean firing rate 4.5412 over everything;
legend('N-PLDS','PLDS')
set(fig7,'Units','centimeters');
set(fig7,'Position',[1,1,20,10]);
print('-depsc2', 'Figure3_and_4/Final_plots/Panel_3Bm.eps');


