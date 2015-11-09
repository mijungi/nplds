function params = generate_params_from_prior(k, d, p, r, Model, ind_train, tau_init)

%% generate hyperparameters 

params = generate_hyperparams(k, Model);

% verycertainprecision = 1e6;
abitcertainprecision = 10;

%% (1) generate A from prior of A

params.A = generate_A_from_prior(params, k, r, Model);

covA = 1/abitcertainprecision*eye(k);
params.covA = repmat(covA, [ 1 1 r]);

if Model=='NSLS'
    AA = zeros(k, k, r);
    for i=1:r
        AA(:,:,i) = params.A(:,:,i)'*params.A(:,:,i) + k*params.covA(:,:,i);
    end
    
else % Single A
    AA = params.A'*params.A + k*covA;
end

params.AA = AA;

%% (2) generate B

params.B = generate_B_from_prior(params, d, k, Model);

covB = 1/abitcertainprecision*eye(d);
params.covB = repmat(covB, [ 1 1 k]);

if strcmp(Model, 'NSFR') || strcmp(Model, 'PLDS')
    params.AB = params.A'*params.B;
end

%% (4) generate C from prior of C

if strcmp(Model, 'NSFR') || strcmp(Model, 'PLDS')
    % there is no prior on C
    C = 0.2*randn(p,k);
    params.C = C;
    
else
    
    Cgamma = params.gamma;
    params.C = generate_C_from_prior(p, k, Cgamma);
    covC = 1/verycertainprecision*eye(k);
    params.covC = repmat(covC, [1 1 p]);
end

%% (4) the rest

params.d = - 2*ones(p,1);

params.V0 = eye(k); params.x0 = zeros(k,1); params.Q = eye(k);

%% (5) h

if strcmp(Model, 'NSFR') || strcmp(Model, 'PLDS')
    [params.h, params.K] = generate_h_from_prior(params, k, r, ind_train, tau_init);
end

covh = 1/abitcertainprecision*eye(k);
params.covh = repmat(covh, [ 1 1 r]);


