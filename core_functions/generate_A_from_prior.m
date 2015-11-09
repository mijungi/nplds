function A = generate_A_from_prior(hyperparams, k, r, Model)

%%
if strcmp(Model, 'NSFR') || strcmp(Model, 'PLDS')
    
    A = 0.1*rand(k);
    A = A/max(abs(eigs(A)))*hyperparams.Anrm;
    
else
    
    %% generate A from prior of A
    
    % amean: mean of A
    % k: size(A) = k by k
    % sig2: marginal variance of A
    % tau2: length scale
    % epsilon: adding some small value to diag of K\
    % r: number of recordings
    
    % note: both sig2 and tau2 should be non-negative
    
    %% generate covariance K
    
    k2 = k^2;
    K = makeK(sig2, epsilon, tau2, r, k2);
    
%     figure(1); subplot(2,4,8); imagesc(K); axis('image'); title('K');
    
    % draw samples of a
    a = mvnrnd(repmat(amean, r, 1)', K);
    
    % reshape a to matrix A:
    A = permute(reshape(a, k, k, r), [2 1 3]);
    
    for i=1:r
        A(:,:,i) = A(:,:,i)/max(abs(eigs(A(:,:,i))))*Anrm;
    end
    
end