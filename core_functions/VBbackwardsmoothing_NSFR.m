function bs = VBbackwardsmoothing_NSFR(Yn, initparams, opts)

%% unpack params

A = initparams.A;
AA = initparams.AA;

B = initparams.B;

AB = initparams.AB;

C = initparams.C;

h = initparams.h;
covh = initparams.covh;

d = initparams.d; 

inpn = initparams.inpn;

%% backward smoothing 

T = size(Yn,2);
k = size(A,1);
p = size(C,1);

maxexp = 1e2; 
I = eye(k);

etatilde = zeros(k, T);
inv_phitilde = zeros(k, k, T);

inv_phistar = zeros(k, k, T);

mubwrd = zeros(k, T);
inv_sigbwrd = zeros(k, k, T);

% for t = T, inv(sigbwrd) = 0 to ensure beta(x_T) = 1
t = T;

% compute mean
a0 = zeros(k,1);
fun = @(a)(obj_bwrd0_NSFR(a, Yn(:,t), C, d, h, covh, maxexp));
etatilde(:,t) = fsolve(fun, a0, opts);

% compute cov
quad = diag(C*covh*C');
g = exp(d + C*(etatilde(:,t)+h) + 0.5*quad);
% g = exp(min(d + C*(etatilde(:,t)+h) + 0.5*quad, maxexp*ones(p,1)));

cctrp = zeros(k,k,p);
for i=1:p
    cstrp = C(i,:);
    cctrp(:,:,i) = cstrp'*cstrp;
end

frstTrm = bsxfun(@times, cctrp, reshape(g,[1 1 p]));

inv_phitilde(:,:,t) = sum(frstTrm,3);

inv_phistar(:,:,t) = I + inv_phitilde(:,:,t);

tmp = AA - A'*(inv_phistar(:,:,t)\A);
% diagtmp = diag(tmp);
% diagtmp(diagtmp==0) = 1./maxexp;
% tmp(logical(eye(k))) = diagtmp;

inv_sigbwrd(:,:,t-1) = tmp; 

% mubwrd(:,t-1) = inv_sigbwrd(:,:,t-1)\( (A'/inv_phistar(:,:,t))*(inv_phitilde(:,:,t)*etatilde(:,t)));
mubwrd(:,t-1) = inv_sigbwrd(:,:,t-1)\( (A'/inv_phistar(:,:,t))*(B*inpn(:,t) + inv_phitilde(:,:,t)*etatilde(:,t)) - AB*inpn(:,t) );

%%

for t = T-1:-1:2

%     fun = @(a) C'*(Yn(:,t) - exp(min(C*a+d, maxexp*ones(p,1)))) - inv_sigbwrd(:,:,t)*(a - mubwrd(:,t));
    a0 = etatilde(:,t+1);
    fun = @(a)(obj_bwrd_NSFR(a, mubwrd(:,t), inv_sigbwrd(:,:,t), Yn(:,t), C, d, h, covh, maxexp));
    etatilde(:,t) = fsolve(fun, a0, opts);
    
    % compute cov
    g = exp(d + C*(etatilde(:,t)+h) + 0.5*quad);
    frstTrm = bsxfun(@times, cctrp, reshape(g,[1 1 p]));
    inv_phitilde(:,:,t) = sum(frstTrm,3) + inv_sigbwrd(:,:,t);
    
    inv_phistar(:,:,t) = I + inv_phitilde(:,:,t);
    
    % to make sure inv_sigbwrd is always dividable (to avoid getting NaN in
    % mubwrd)
    tmp = AA - A'*(inv_phistar(:,:,t)\A);
%     diagtmp = diag(tmp);
%     diagtmp(diagtmp==0) = 1/maxexp;
%     tmp(logical(eye(k))) = diagtmp;
    
    inv_sigbwrd(:,:,t-1) = tmp;
    
%     mubwrd(:,t-1) = inv_sigbwrd(:,:,t-1)\( (A'/inv_phistar(:,:,t))*( inv_phitilde(:,:,t)*etatilde(:,t)) );
    mubwrd(:,t-1) = inv_sigbwrd(:,:,t-1)\( (A'/inv_phistar(:,:,t))*(B*inpn(:,t) + inv_phitilde(:,:,t)*etatilde(:,t)) - AB*inpn(:,t) );
    
end

%%
t = 1;

a0 = etatilde(:,t+1);
fun = @(a)(obj_bwrd_NSFR(a, mubwrd(:,t), inv_sigbwrd(:,:,t), Yn(:,t), C, d, h, covh, maxexp));
etatilde(:,t) = fsolve(fun, a0, opts);

% compute cov
g = exp(d + C*(etatilde(:,t)+h) + 0.5*quad);
frstTrm = bsxfun(@times, cctrp, reshape(g,[1 1 p]));
inv_phitilde(:,:,t) = sum(frstTrm,3) + inv_sigbwrd(:,:,t);

inv_phistar(:,:,t) = I + inv_phitilde(:,:,t);

% to make sure inv_sigbwrd is always dividable (to avoid getting NaN in
% mubwrd)
tmp = AA - A'*(inv_phistar(:,:,t)\A);
% diagtmp = diag(tmp);
% diagtmp(diagtmp < 1/maxexp) = 1/maxexp;
% tmp(logical(eye(k))) = diagtmp;

inv_sigbwrd0 = tmp;

% mubwrd0 = inv_sigbwrd0\( (A'/inv_phistar(:,:,t))*( inv_phitilde(:,:,t)*etatilde(:,t)) );
mubwrd0 = inv_sigbwrd0\( (A'/inv_phistar(:,:,t))*(B*inpn(:,t) + inv_phitilde(:,:,t)*etatilde(:,t)) - AB*inpn(:,t) );

%% return these

bs.mubwrd = mubwrd;
bs.inv_sigbwrd = inv_sigbwrd;
bs.mubwrd0 = mubwrd0;
bs.inv_sigbwrd0 = inv_sigbwrd0;


%% ===========================================================
function obj = obj_bwrd0_NSFR(a, y, C, d, h, covh, maxexp)

% [quad, lin] = xSigcx(a, covC);
p = size(C,1);
% obj = C'*y - (C' + lin)*exp(min(d +  C*a + 0.5*quad, maxexp*ones(p,1)));

quad = diag(C*covh*C');
secondTrm_in_mean = exp(min(d + C*(a+h) + 0.5*quad, maxexp*ones(p,1)));

obj = C'*y - C'*secondTrm_in_mean;

%% ===========================================================
function obj = obj_bwrd_NSFR(a, eta, inv_phi, y, C, d, h, covh, maxexp)

% [quad, lin] = xSigcx(a, covC);
p = size(C,1);
% obj = C'*y - (C' + lin)*exp(min(d + C*a + 0.5*quad, maxexp*ones(p,1))) - inv_phi*(a - eta);
quad = diag(C*covh*C');
secondTrm_in_mean = exp(min(d + C*(a+h) + 0.5*quad, maxexp*ones(p,1)));

obj = C'*y - C'*secondTrm_in_mean - inv_phi*(a - eta);


