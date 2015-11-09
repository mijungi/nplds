
function hyperparams = generate_hyperparams(k, Model)


if strcmp(Model, 'NSFR') || strcmp(Model, 'PLDS')
   %% PLDS with nonstationary firing rates
      
    % (1) hyperparameters for h
    hyperparams.sig2 = 10*abs(rand); % marginal variance on h
    hyperparams.tau2 = 300; % lengthscale
    hyperparams.epsilon = 1e-3;
        
    % (2) hyperparameters on A
    hyperparams.Anrm = 0.7;
    hyperparams.alpha = 10*abs(rand);
    
    % (3) hyperparameters on B
    hyperparams.beta = 10*abs(rand); 
    
    % (2) hyperparameters for C
%     hyperparams.gamma = 10*abs(rand);
%     hyperparams.gamma = 1;
    
else
    
    %% for less messy optimization
    % 1/sig2 and gamma shouldn't be too large...
    
    % (1) hyperparameters for A
    hyperparams.sig2 = 10*abs(rand); % marginal variance on A
    hyperparams.tau2 = 10*abs(rand); % lengthscale
    hyperparams.epsilon = 1e-3;
    hyperparams.Anrm = 0.9;
    
    Amean = randn(k);
    Amean = Amean/max(abs(eigs(Amean)))*hyperparams.Anrm;
    
    hyperparams.Amean = reshape(Amean', [], 1);
    
    % (2) hyperparameters for C
    % hyperparams.gamma = 10*abs(rand);
    hyperparams.gamma = 1;
    
end