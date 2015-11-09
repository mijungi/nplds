function [h, covh, m_h, sig2, tau2, logHyperLandscape] = VBcomputeh_C_NSFR(xyzinpn, r, fromVBEstep, params, opts2, ind_train)

%% unpack init params

C = params.C;
d = params.d;

%% (1) given initial hyperparams, update mu_h and cov_h

m_h = params.m_h;
sig2 = params.sig2;
tau2 = params.tau2;
epsilon = params.epsilon;

[p, k] = size(C);
% generate covariance K
K = makeK(sig2, epsilon, tau2, r, k, ind_train);

prs0 = randn(k*r,1)+repmat(m_h,r,1); % h^{(1:r)}
fun = @(prs) (loglikeliFunc_h(prs, C, d, fromVBEstep, r, xyzinpn, K, m_h));

%fprintf(['M-step VBcomputeh_C_NFSR mu_h cov_h started ' datestr(clock, 'yy-mm-dd-HH:MM:SS') '\n']);
h = fminunc(fun, prs0, opts2);
[~,~,invcov] = fun(h);
covh = inv(invcov);
%fprintf(['M-step VBcomputeh_C_NFSR mu_h cov_h done ' datestr(clock, 'yy-mm-dd-HH:MM:SS') '\n']);

%% (2) given mu_h and cov_h, find hyperparams that minimize KL(h)

  m_h = mean(reshape(h,k,r),2); % Update the prior mean of h
  fun = @(prs)(KLA(prs, epsilon, r, k, h, covh, m_h, ind_train));

  %FITPARAM - Set the grid for hyperparam update
  %Define the grid over which we compute noise values and temporal correlation values
  %in the kernel hyperparams
  logsig2 = [-3:0.2:2];
  logtau2 = [-6, -2:1:2, 2.5,3:0.05:5,5.5,6,7]; % first one for no trial-to-trial interaction at all, then to model all reasonable interactions from timescale of single-trial correlation to full experiment correlation

  [st1, st2] = ndgrid(logsig2, logtau2);

  funval = zeros(length(st1(:)),1);
  st12 = [st1(:) st2(:)];

  for i=1:length(st1(:))

      funval(i) = fun(st12(i,:));
      %fprintf(['M-step VBcomputeh_C_NFSR hyperparam iter ' num2str(i) '/' num2str(length(st12(i,:))) ' done ' datestr(clock, 'yy-mm-dd-HH:MM:SS') '\n']);

  end
  % 
  % % plot(logsig2, funval, 'o-')
  % 
  % figure(101);
  % imagesc(exp(logsig2), exp(logtau2), reshape(funval, length(logsig2), []))
  [minvalu, minloca] = min(funval);
  % 
  % % dost this match 
  % [exp(st12(minloca,:))]

  logHyperLandscape = [st12, funval];

  sig2 = exp(st12(minloca,1));
  tau2 = exp(st12(minloca,2));
  

%% (3) once setting hyperparams, update mu_h and cov_h again

K = makeK(sig2, epsilon, tau2, r, k, ind_train);

prs0 = randn(k*r,1)+repmat(m_h,r,1); % h^{(1:r)}
fun = @(prs) (loglikeliFunc_h(prs, C, d, fromVBEstep, r, xyzinpn, K, m_h));

h = fminunc(fun, prs0, opts2);
[~,~,invcov] = fun(h);

% return these
% (1) mean:
h = reshape(h, k, r);
% (2) cov:
covhtot = inv(invcov);
covh = zeros(k,k,r);
for i=1:r
    covh(:,:,i) = covhtot(k*(i-1)+1:i*k,k*(i-1)+1:i*k);
end

% ===========================================================
function [l, dl, ddl] = loglikeliFunc_h(prs, C, d, fromVBEstep, r, xyzinpn, K, m_h)

Ir = eye(r);
Cbd = kron(Ir,C);

p = size(C,1);
ymat = zeros(p,r);
gmat = zeros(p,r);

for i=1:r
    
    %% unpack mu/cov of latent variables
    if r==1
        mumarg = fromVBEstep.mumarg;
        inv_sigmarg = fromVBEstep.inv_sigmarg;
        ymat(:,i) = sum(xyzinpn.y,2);
    else
        
        mumarg = fromVBEstep{i}.mumarg;
        inv_sigmarg = fromVBEstep{i}.inv_sigmarg;
        ymat(:,i) = sum(xyzinpn{i}.y,2);
    end
        
    T = size(mumarg,2);
    obj_scndTrm = zeros(p,T);
    
    for t=1:T
        
        CUpsilon = C/inv_sigmarg(:,:,t);
        f = diag(CUpsilon*C');
        g = C*mumarg(:,t) + 0.5*f + d;
        
        tmp = exp(g);
        
        obj_scndTrm(:,t) = tmp;
        
    end
    
    gmat(:,i) = sum(obj_scndTrm,2);
end

yy = reshape(ymat,[],1);
gg = reshape(gmat,[],1);
obj = prs'*(Cbd'*yy) - sum(exp(Cbd*prs).*gg) - 0.5*((prs-repmat(m_h,r,1))'/K)*(prs-repmat(m_h,r,1));

l = -obj;

if nargout >= 2
  dobj = Cbd'*yy - Cbd'*(exp(Cbd*prs).*gg) - K\(prs-repmat(m_h,r,1));
  dl = -dobj;
end

if nargout >=3
  ddobj = - Cbd'*diag(exp(Cbd*prs).*gg)*Cbd - inv(K);
  ddl = -ddobj;
end

%% =============================================================================
function [l] = KLA(prs, epsilon, r, k, jointMeanA, jointCovA, m_h, ind_train)

% unpack params
sig2 = exp(prs(1));
tau2 = exp(prs(2));

jointMeanA = jointMeanA -repmat(m_h,r,1); %Correct for the nonzero prior mean
% jointCovA = jointCovA - (mean(jointMeanA,2)*mean(jointMeanA,2)');

K = makeK(sig2, epsilon, tau2, r, k, ind_train);

KinvJointAs = (K\(jointCovA + (jointMeanA)*(jointMeanA)'));
l = - 0.5*logdetns(K\jointCovA) + 0.5*trace(KinvJointAs);



