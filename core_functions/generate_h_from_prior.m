function [h, K] = generate_h_from_prior(hyperparams, k, r, ind_train, tau_init)
   
    %% unpack hyperparams
    
    sig2 = hyperparams.sig2;
%     tau2 = hyperparams.tau2;
    tau2 = tau_init; 
    epsilon = hyperparams.epsilon;

    %% generate covariance K
    
    K = makeK(sig2, epsilon, tau2, r, k, ind_train);
    
%     figure(1); subplot(2,4,8); imagesc(K); axis('image'); title('K');
    
    % draw samples of h
    meanh = zeros(1, k*r);
    h = mvnrnd(meanh, K);
    
    % reshape a to matrix h:
    h = reshape(h, k, []);
    
