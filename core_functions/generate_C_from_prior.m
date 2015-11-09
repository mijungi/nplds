function C = generate_C_from_prior(p, k, Cgamma)
%% generate C from prior of C

ccov = eye(p*k);
c = mvnrnd(zeros(1, p*k), 1/Cgamma*ccov);

C = reshape(c, k, p)';

