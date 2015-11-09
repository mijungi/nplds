function tau_init = GPfittingTrainingData(resps)

% fit data to GP for initialization of tau (length scale)

n_neurons = size(resps, 1);
rtrain = size(resps, 3);

%%
rates=squeeze(mean(resps,2));

for i=1:n_neurons
    %for each neuron and each orientation, we do four things
    
    stats.mean(i)=mean(rates(i,:));
    stats.var(i)=var(rates(i,:),0);
    stats.rates(i,:)=rates(i,:);
    
    %now, lets use GP regression to find the best smoothing kernel for each
    %neuron. This is probably a bit of an overkill but maybe important and
    %useful later on:
    x=[[1:rtrain]'];
    y=rates(i,:)' - mean(rates(i,:)');
    covfunc = {'covSum', {'covSEiso','covNoise'}};
    loghyper = [log(1.0); log(1.0); log(0.1)];
    loghyper = minimize(loghyper, 'gpr', -100, covfunc, x, y);
    stats.gp.tau(i)=exp(loghyper(1)); % this is exp(-(xi-xj)^2/tau), length scale
    stats.gp.var(i)=exp(loghyper(2));
    stats.gp.noisevar(i)=exp(loghyper(3));
    [mu] = gpr(loghyper, covfunc, x, y,x);
    stats.gp.smoothrate(i,:)= mu + mean(rates(i,:)');
    %keyboard
    
end

tau_init = median(stats.gp.tau); 