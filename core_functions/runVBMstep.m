function fromVBMstep = runVBMstep(xyzinpn, fromVBEstep, initparams, r, ind_train,Model)
% fromVBMstep = runVBMstep(Yn, fromVBEstep, initparams, trueparams, r);

%% (1) update x0 and V0

x0 = initparams.x0;
V0 = initparams.V0;

% x0 = fromEstep.mumarg0;
% V0 = inv(fromEstep.inv_sigmarg0);

%% (2) update A and B

% A = initparams.A;
% B = initparams.B;
% AA = initparams.AA;
% AB = initparams.AB;
% covA = initparams.covA;
% covB = initparams.covB;

%fprintf(['M-step started ' datestr(clock, 'yy-mm-dd-HH:MM:SS') '\n']);

% tic;
[A, B, AA, AB, covA, covB, alpha, beta] = VBcomputeAB_NSFR(r, fromVBEstep, initparams);

fprintf(['M-step VBcomputeAB_NSFR done ' datestr(clock, 'yy-mm-dd-HH:MM:SS') '\n']);
% toc;

%% (4) update C and h

% C = initparams.C;
% h = initparams.h;
% covh = initparams.covh;

% sig2 = initparams.sig2;
% tau2 = initparams.tau2;

opts1 = optimset('Display', 'off', 'GradObj' , 'on', 'Hessian', 'off', 'maxIter', 1e3);
% tic;
C = VBcomputeC_h_NSFR(xyzinpn, r, fromVBEstep, initparams, opts1);

fprintf(['M-step VBcomputeC_h_NSFR done ' datestr(clock, 'yy-mm-dd-HH:MM:SS') '\n']);
% toc;

% update C from optimization above
initparams.C = C;

opts2 = optimset('Display', 'off', 'GradObj' , 'on', 'Hessian', 'on', 'maxIter', 1e3);
% tic; 
if strcmp(Model,'NSFR')
  [h, covh, m_h, sig2, tau2, logHyperLandscape] = VBcomputeh_C_NSFR(xyzinpn, r, fromVBEstep, initparams, opts2, ind_train);
elseif strcmp(Model, 'PLDS')
  h = zeros(size(initparams.h));
  covh = repmat(zeros(size(h,1)),[1,1,size(h,2)]);
  sig2 = 1;
  tau2 = 1;
  m_h = zeros(size(h,1),1);
  logHyperLandscape = [1,1,1];
end

fprintf(['M-step VBcomputeh_C_NSFR done ' datestr(clock, 'yy-mm-dd-HH:MM:SS') '\n']);
% toc;

%%
% update h from optimization above
% initparams.h = h;
% initparams.sig2 = sig2;
% initparams.tau2 = tau2;
% initparams.covh = covh;
% 
% % then optimize C given h again
% C = VBcomputeC_h_NSFR(xyzinpn, r, fromVBEstep, initparams, opts1);
% 
% initparams.C = C;
% [h, covh, sig2, tau2] = VBcomputeh_C_NSFR(xyzinpn, r, fromVBEstep, initparams, opts2);

%% return these:

fromVBMstep.A = A;
fromVBMstep.AA = AA;
fromVBMstep.covA = covA;
fromVBMstep.alpha = alpha;

fromVBMstep.B = B;
fromVBMstep.covB = covB;
fromVBMstep.AB = AB;
fromVBMstep.beta = beta;

fromVBMstep.C = C;

fromVBMstep.h = h;
fromVBMstep.covh = covh;
fromVBMstep.m_h = m_h;
fromVBMstep.sig2 = sig2;
fromVBMstep.tau2 = tau2;

fromVBMstep.d = initparams.d; 

fromVBMstep.x0 = x0;
fromVBMstep.V0 = V0;
fromVBMstep.Q = initparams.Q;

fromVBMstep.logHyperLandscape = logHyperLandscape;

