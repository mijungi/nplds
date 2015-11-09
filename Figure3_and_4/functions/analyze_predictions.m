data_dir = '/nfs/data3/gergo/Mijung/Figure4/Output_data_retesting';
data_file = '/nfs/data3/gergo/Mijung/Figure4/Input_data/alexdata_session2_org0_bin50.mat';
code_dir = '/nfs/nhome/live/gbohner/Git_main/nonstat_plds_code/Gergo';
output_data_file_temp = 'final_MODEL_data.mat';
output_params_file = 'init_data.mat';

Models = {'NSFR','PLDS'};

% sig_hist = zeros(10,10,31); % k x run x [initial iter]
% tau_hist = zeros(10,10,31); % k x run x [initial iter]
% mh_hist = zeros(10,10,10,31); % size(mh) x k x run x [initial iter]
tau_vals = zeros(1,51);
tau_p_vals = zeros(9,51);

% load('hists')

for mod = 1:2
Model = Models{mod};
output_data_file = strrep(output_data_file_temp, 'MODEL', Model);
fr_pred_error = zeros(9,10); % # of ks time # of runs
for k = 1:2
    for runnum = 0:3
      

      [mod k runnum]
      output_data_dir = [data_dir filesep 'k' num2str(k) '_run' num2str(runnum)];
      saved_data_file = [output_data_dir filesep output_data_file];
      saved_params_file = [output_data_dir filesep output_params_file];
      fr_pred_error(k, runnum+1) = make_predictions( data_file, saved_data_file, saved_params_file, code_dir, k, Model);
%       if fr_pred_error(k,runnum+1) > 2.5
%         fr_pred_error(k, runnum+1) = make_predictions( data_file, saved_data_file, saved_params_file, code_dir, k, Model, 1);
%       end
%       continue;
%       if mod == 1
%         load(saved_params_file);
%         sig_hist(k,runnum+1,1) = params.sig_init;
%         tau_hist(k,runnum+1,1) = params.tau_init;
%         mh_hist(1:k, k,runnum+1,1) = params.m_h_init;
%         load(saved_data_file);
%         for i1=1:30
%           sig_hist(k,runnum+1,i1+1) = datastruct.Mstep{i1}.sig2;
%           tau_hist(k,runnum+1,i1+1) = datastruct.Mstep{i1}.tau2;
%           mh_hist(1:k, k,runnum+1,i1+1) = datastruct.Mstep{i1}.m_h;
%         end
%       end
    end
%     if mod ==1
%       for runnum =0:9
%         output_data_dir = [data_dir filesep 'k' num2str(k) '_run' num2str(runnum)];
%         saved_data_file = [output_data_dir filesep output_data_file];
%         saved_params_file = [output_data_dir filesep output_params_file];
%         load([saved_data_file(1:end-4) '_predfr.mat'], 'Mstepresults', 'hpred', 'params','fr_mean_pred_all', 'fr_mean_true');
%         fr_mean_pred = median(fr_mean_pred_all,3);
%         fr_pred_error_trial = sqrt(mean((fr_mean_pred - fr_mean_true).^2,1)); %RMSE value
%         fr_pred_error(k, runnum+1) = sqrt(sum(sum((fr_mean_pred - fr_mean_true).^2)));
% %   %       figure(3); scatter(runnum*10 + (1:10), fr_pred_error_trial); hold on;
% %   %       figure(3); scatter(params.ind_test, fr_pred_error_trial); hold on;
% %         load([saved_data_file])
% %         tmp = fr_pred_error(k,:);
% %         tmp = ceil(127*mat2gray(tmp,[min(tmp),max(tmp)])+1);
% %         clr = colormap(jet(128));
% %         tmp = clr(tmp,:);
% % %         figure(2); hold on; plot(1:30, squeeze(tau_hist(k,runnum+1,2:31))', 'Color', tmp(runnum+1,:));
% %         figure(3); hold on; plot(params.ind_train, datastruct.Mstep{end}.h)%, 'Color', tmp(runnum+1,:)); colorbar('YTickLabel', {sort(fr_pred_error(k,:))}); 
% %         figure(4); hold on; scatter(repmat(k,k,1), mean(datastruct.Mstep{end}.h,2));%, 'Color', tmp(runnum+1,:));
% % %         title('h changes and predictions (k == 4, color ~ RMSE)'); xlabel('Trial'); ylabel('<h>'); 
% % %         scatter(params.ind_test, mean(hpred,1)); 
% %         figure(3); hold on; scatter(reshape(repmat(params.ind_test,size(hpred,1),1),1,[])', hpred(:));
%       end
%     end
%     figure(13);
%     subplot(2,5,k)
%     [tau_vals, tmp] = show_p_value_contour(data_dir, output_data_file, k, 30);
%     tau_p_vals(k,:) = tmp;
end

if mod == 1
NSFR_error = fr_pred_error;
else
PLDS_error = fr_pred_error;
end

end

% figure(13);
% print_to = ['/nfs/nhome/live/gbohner/Dropbox/Gatsby/Random/Mijung_NPLDS/Code/Figure4/Plots/contours' num2str(k) '.eps'];
% print('-depsc2', print_to);
% 
% % figure('visible','off');
% figure(9);
% scatter(reshape(repmat((1:10)-.1,10,1),100,1),reshape(NSFR_error',100,1),505,'b.')
% hold on;
% scatter(reshape(repmat((1:10)+.1,10,1),100,1),reshape(PLDS_error',100,1),505,'r.')
% xlabel('k')
% title('RMSE in firing rates on left out trials (10%) after 30 iter')
% ylabel('RMSE')
% ylim([1,4])
% legend('NSFR','PLDS')
% print('-depsc2', 'Figure4/Plots/error_scatter.eps');
% 
% % figure('visible', 'off'); 
% figure(10);
% errorbar((1:10)-.1, mean(NSFR_error,2), std(NSFR_error,0,2), 'bx');
% hold on; errorbar((1:10)+.1, mean(PLDS_error,2), std(PLDS_error,0,2), 'rx');
% ylim([1 4])
% xlabel('k')
% title('RMSE in firing rates on left out trials (10%) after 30 iter')
% ylabel('RMSE')
% legend('NSFR','PLDS')
% print('-depsc2', 'Figure4/Plots/error_std.eps');

% sig_min = hp_pos(find(hp_val_mean==min(hp_val_mean)),1);
% 
% inds = find(hp_pos(:,1)==sig_min);
% 
% % errorbar(hp_pos(inds,2),hp_val_mean(inds),hp_val_std(inds));
% % plot(hp_pos(inds,2),hp_val(inds,:));
% % surf(tau_vals, sig_vals, log(log(hp_val_mean_2d)));
% % plot3_errorbars_surf(tau_vals, sig_vals, hp_val_mean_2d, hp_val_std_2d);
% 
% plot3_errorbars_surf(hp_pos(:,1),hp_pos(:,2),hp_val_mean,hp_val_std);

%  for k = 1:10
% figure(5); hold on; scatter(repmat(k,1,10), fr_pred_error(k,:),55,'x')
% end
% xlim([0,5])
% xlabel('k')
% title('RMSE in firing rates on left out trials (10%)')
% ylabel('RMSE')
% title('nPLDS RMSE in firing rates on left out trials (10%)')


% %Plot norm of m_h vs k
% figure(71);
% scatter(reshape((randn(10,1)/225+1)*(1:10),1,[]),reshape(squeeze(sqrt(sum((squeeze(mh_hist(:,:,:,end)).^2),1)))',1,[]),300,'.')
% xlabel('k');
% ylabel('norm(m_h)');
% print('-depsc2','/nfs/nhome/live/gbohner/Dropbox/Gatsby/Random/Mijung_NPLDS/Code/Figure4/Plots/mh_final_val.eps')
% 
% %Plot sig vs k
% figure(72);
% scatter(reshape((randn(10,1)/225+1)*(1:10),1,[]),reshape(squeeze(sig_hist(:,:,end))',1,[]),300,'.')
% xlabel('k');
% ylabel('sig^2');
% print('-depsc2','/nfs/nhome/live/gbohner/Dropbox/Gatsby/Random/Mijung_NPLDS/Code/Figure4/Plots/sig_final_val.eps')
% 
% %Plot tau vs k
% figure(73);
% scatter(reshape((randn(10,1)/225+1)*(1:10),1,[]),reshape(squeeze(tau_hist(:,:,end))',1,[]),300,'.')
% xlabel('k');
% ylabel('tau^2');
% print('-depsc2','/nfs/nhome/live/gbohner/Dropbox/Gatsby/Random/Mijung_NPLDS/Code/Figure4/Plots/tau_final_val.eps')
% 
% %Convergence - changes in tau
% figure(74);
% plot(1:30,sum(diff(reshape(tau_hist,100,31),1,2)~=0,1))
% xlabel('Iteration')
% ylim([0,100])
% ylabel('% of runs where tau changed')
% print('-depsc2','/nfs/nhome/live/gbohner/Dropbox/Gatsby/Random/Mijung_NPLDS/Code/Figure4/Plots/tau_change_perc.eps')
% 
% %Convergence - changes in sig
% figure(74);
% plot(1:30,sum(diff(reshape(sig_hist,100,31),1,2)~=0,1))
% xlabel('Iteration')
% ylim([0,100])
% ylabel('% of runs where \sigma changed')
% print('-depsc2','/nfs/nhome/live/gbohner/Dropbox/Gatsby/Random/Mijung_NPLDS/Code/Figure4/Plots/sig_change_perc.eps')
