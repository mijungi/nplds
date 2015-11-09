
function [A, B, AA, AB, covA, covB, alpha, beta] = VBcomputeAB_NSFR(r, fromVBEstep, initparams)


%% (1) estimate joint mean and covariance of A and B

% step0: unpack latent sufficient statistics
k = size(initparams.A,1);
p = size(initparams.C,1);
pinp = size(initparams.B,2);

WAmat = zeros(k,k,r);
SAmat = zeros(k,k,r);
GAmat = zeros(k,pinp,r);
Mtilmat = zeros(pinp, k, r); 

if r==1
    
    WAmat = fromVBEstep.suffstat.WA;
    SAmat = fromVBEstep.suffstat.SA;
    GAmat = fromVBEstep.suffstat.GA;
    Mtilmat = fromVBEstep.suffstat.Mtil;
    
else
    
    for i=1:r
        
        WAmat(:,:,i) = fromVBEstep{i}.suffstat.WA;
        SAmat(:,:,i) = fromVBEstep{i}.suffstat.SA;
        GAmat(:,:,i) = fromVBEstep{i}.suffstat.GA;
        Mtilmat(:,:,i) = fromVBEstep{i}.suffstat.Mtil;
        
    end
    
end

Ik = eye(k);

W = kron(Ik,sum(WAmat,3)); 
G = kron(Ik,sum(GAmat,3));

sumSA = sum(SAmat,3);
s = sumSA(:);

sumM = sum(Mtilmat,3);
m = sumM(:);

% Udot (same from every trial)

if r==1
    U = kron(Ik, fromVBEstep.suffstat.Udot);
else
    U = kron(Ik, fromVBEstep{1}.suffstat.Udot);
end

%% update mu_a, mu_b, cov_a, cov_b, cov_ab given initial hyperparams

alpha = initparams.alpha;
beta = initparams.beta;

% Hessians

Ik2 = eye(k^2);
Ha = alpha*Ik2 + W;
Hab = G;
Irpinp = eye(k*pinp);
Hb = beta*Irpinp + r*U;

% covariances

covA = inv(Ha - (Hab/Hb)*Hab');
covB = inv(Hb- (Hab'/Ha)*Hab);
covAB = - covA*(Hab/Hb);

% means

mu_a = covA*(s- (Hab/Hb)*m);
mu_b = covB*(m-(Hab'/Ha)*s);

%% optimize hyperparams given mu_a, mu_b, cov_a, cov_b, cov_ab

logalpha = [-6:0.2:6];
logbeta = [-6:0.2:6];

[st1, st2] = ndgrid(logalpha, logbeta);

funval = zeros(length(st1(:)),1);
st12 = [st1(:) st2(:)];

fun = @(prs) (KLAB(prs, k, pinp, mu_a, mu_b, covA, covB, covAB));

for i=1:length(st1(:))
    
    funval(i) = fun(st12(i,:));
    
end

% figure(101);
% imagesc(exp(logalpha), exp(logbeta), reshape(funval, length(logalpha), []))
[minvalu, minloca] = min(funval);
% 
% % dost this match 
%[exp(st12(minloca,:))]

alpha = exp(st12(minloca,1));
beta = exp(st12(minloca,2));

%% update mu_a, mu_b, cov_a, cov_b, cov_ab given a new set of hyperparams

Ik2 = eye(k^2);
Ha = alpha*Ik2 + W;
Hab = G;
Irpinp = eye(k*pinp);
Hb = beta*Irpinp + r*U;

% covariances

covA = inv(Ha - (Hab/Hb)*Hab');
covB = inv(Hb- (Hab'/Ha)*Hab);
covAB = - covA*(Hab/Hb);

% means

mu_a = covA*(s- (Hab/Hb)*m);
mu_b = covB*(m-(Hab'/Ha)*s);

% return these

A = reshape(mu_a, k, k)';
B = reshape(mu_b, pinp, k)';
covA = covA(1:k, 1:k);
covB = covB(1:pinp, 1:pinp);
covAB = covAB(1:k, 1:pinp);

AA = A'*A + k*covA;
BB = B'*B + k*covB;
AB = A'*B - (covAB/covB)*BB;

%% =====================================
function l = KLAB(prs, k, pinp, mu_a, mu_b, covA, covB, covAB)

%% unpack hyperparams

alpha = exp(prs(1));
beta = exp(prs(2));

%% computing KL(a,b)

Ik2 = eye(k^2);
Irpinp = eye(k*pinp);
Lambda_prior = [1/alpha*Ik2 zeros(k^2, k*pinp); zeros(k*pinp, k^2) 1/beta*Irpinp];
Lambda_post = [covA covAB; covAB' covB];
mu_ab = [mu_a; mu_b];

l = -0.5*logdetns(Lambda_prior\Lambda_post)+ 0.5*trace(Lambda_prior\(Lambda_post + mu_ab*mu_ab'));





