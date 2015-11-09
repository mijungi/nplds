function B = generate_B_from_prior(hyperparams, d, k, Model)


if strcmp(Model, 'NSFR') || strcmp(Model, 'PLDS')
    
    %% generate B from prior of B
    
    covb = eye(d*k);
    b = 0.1*mvnrnd(zeros(1, d*k), 1/hyperparams.beta*covb);
    
    B = reshape(b, d, k)';
    
else
    
    B = [];
    
end
