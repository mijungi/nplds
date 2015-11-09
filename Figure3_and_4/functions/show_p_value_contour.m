function [tau_out, p_out] = show_p_value_contour(data_dir, output_data_file, k, varargin)



for runnum = 0:9
  output_data_dir = [data_dir filesep 'k' num2str(k) '_run' num2str(runnum)];
  try
	load([output_data_dir filesep output_data_file], 'datastruct');
  if runnum == 0, hp_pos = datastruct.Mstep{end}.logHyperLandscape(:,1:2); end;
  hp_val(:,runnum+1) = datastruct.Mstep{end}.logHyperLandscape(:,3);
  end
end


hp_val_mean = mean(hp_val,2);
hp_val_std = std(hp_val,[],2);

sig_vals = unique(hp_pos(:,1));
tau_vals = unique(hp_pos(:,2));
% sig_vals = [-4:0.5:4]';
% tau_vals = [-6, -2:1:2, 2,2.5,3:0.05:5,5.5,6,7]';

hp_val_mean_2d = reshape(hp_val_mean,length(sig_vals), length(tau_vals));
hp_val_std_2d = reshape(hp_val_std,length(sig_vals), length(tau_vals));

[min_val, min_ind] = min(hp_val_mean);
min_std = hp_val_std(min_ind) + eps;

t_vals = (hp_val_mean - min_val)./sqrt(hp_val_std.^2/size(hp_val,2) +  min_std.^2/size(hp_val,2));

p_vals = 1-tcdf(t_vals,size(hp_val,2)-1);

%Collapse into tau (by choosing maximum wrt sigma)

tau_out = exp(tau_vals);
p_out = max(reshape(p_vals,length(sig_vals), length(tau_vals)),[],1);

% surf(tau_vals, sig_vals, reshape(p_vals,length(sig_vals), length(tau_vals)));
% try
% C = contour(exp(tau_vals),exp(sig_vals), reshape(p_vals,length(sig_vals), length(tau_vals)));
% clabel(C);
% % title(['P-value contours from 10-fold bootstrap on the hyperparameter likelihood surface (k=' num2str(k) ')'])
% % xlabel('tau (trials)');
% % ylabel('sigma (noise to signal ratio)');
% xlim([0, 200])
% ylim([0,exp(max(sig_vals(:)))])
% % print_to = ['/nfs/nhome/live/gbohner/Dropbox/Gatsby/Random/Mijung_NPLDS/Code/Figure4/Plots/contour_k' num2str(k) '.eps'];
% % print('-depsc2', print_to);
% end
end