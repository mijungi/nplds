function fromVBEstep = VBcomputingmarginals_NSFR(Yn, ff, bs, initparams)

% fromVBEstep = VBcomputingmarginals_A(inpn, Yn, ff, bs, params);

%% unpack msgs

mufwrd = ff.mufwrd;
sigfwrd = ff.sigfwrd;
sigstar = ff.sigstar;
sigstar0 = ff.sigstar0;

mubwrd0 = bs.mubwrd0;
inv_sigbwrd0 = bs.inv_sigbwrd0;
mubwrd = bs.mubwrd;
inv_sigbwrd = bs.inv_sigbwrd;

%% unpack params

A = initparams.A;
C = initparams.C;
d = initparams.d;
h = initparams.h;
covh = initparams.covh;

x0 = initparams.x0;
V0 = initparams.V0;

inpn = initparams.inpn;
pinp = size(inpn,1);

%% computing marignals and cross-covariances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = size(A,1);
T = size(mufwrd,2);
p = size(Yn,1);
mumarg = zeros(k, T);
inv_sigmarg = zeros(k, k, T);
crsscov = zeros(k, k, T-1);

% second deriv of log-likelihood
maxexp = 1e2; 
% W = @(a) C'*diag(exp(min(C*a+d, maxexp*ones(p,1))))*C;

for t=1:T
    
    tmp = inv(sigfwrd(:,:,t)) + inv_sigbwrd(:,:,t);
%     diagtmp = diag(tmp);
%     diagtmp(diagtmp==0) = 1./maxexp;
%     tmp(logical(eye(k))) = diagtmp;
    
    inv_sigmarg(:,:,t) = tmp;

    mumarg(:,t) = inv_sigmarg(:,:,t)\(sigfwrd(:,:,t)\mufwrd(:,t) + inv_sigbwrd(:,:,t)*mubwrd(:,t));
end

inv_sigmarg0 = inv(V0) + inv_sigbwrd0;

tmp = inv_sigmarg0;
% diagtmp = diag(tmp);
% diagtmp(diagtmp==0) = 1./maxexp;
% tmp(logical(eye(k))) = diagtmp;

inv_sigmarg0 = tmp;

mumarg0 = inv_sigmarg0\(inv_sigbwrd0*mubwrd0 + V0\x0);

%% this is different from the standard EM

% t=0,1
% W =  computeHessian(mumarg(:,1), covC, C, d, maxexp);
W =  computeHessian_NSFR(mumarg(:,1), C, d, h, covh, maxexp);
I = eye(k);
RRterm = W  + inv_sigbwrd(:,:,1) + I;
crosscov0 = (inv(sigstar0) - (A'/RRterm)*A)\(A'/RRterm);

% crossvar

for t=1:T-1
%     t
%     W =  computeHessian(mumarg(:,t+1), covC, C, d, maxexp);
    W =  computeHessian_NSFR(mumarg(:,t+1), C, d, h, covh, maxexp);

    RRterm = W + inv_sigbwrd(:,:,t+1) + I;  
    crsscov(:, :, t) = (inv(sigstar(:,:,t)) - (A'/RRterm)*A)\(A'/RRterm);
end



%% computing sufficient statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WA = zeros(k, k, T);
SA = zeros(k, k, T);
SC = zeros(k, p, T); 
WC = zeros(k, k, T);
GA = zeros(k, pinp, T); 
Mtil = zeros(pinp, k, T);
Udot = zeros(pinp, pinp, T);

for t = 1: T
   
    SC(:,:,t) = mumarg(:,t)*Yn(:,t)';
    WC(:,:,t) =  inv(inv_sigmarg(:,:,t)) + mumarg(:,t)*mumarg(:,t)';
    Mtil(:,:,t) = inpn(:,t)*mumarg(:,t)';
    Udot(:,:,t) = inpn(:,t)*inpn(:,t)';
    
    if t ==1
        WA(:,:,t) = inv(inv_sigmarg0) + mumarg0*mumarg0'; 
        SA(:,:,t) = crosscov0 + mumarg0*mumarg(:,t)';
        GA(:,:,t) = mumarg0*inpn(:,t)';
    else
        WA(:,:,t) = inv(inv_sigmarg(:,:,t-1)) + mumarg(:,t-1)*mumarg(:,t-1)';
        SA(:,:,t) = crsscov(:,:,t-1) + mumarg(:,t-1)*mumarg(:,t)';
        GA(:,:,t) = mumarg(:,t-1)*inpn(:,t)';
    end
end
        
        
suffstat.WA = sum(WA,3);
suffstat.SA = sum(SA,3);
suffstat.SC = sum(SC,3);
suffstat.WC = sum(WC,3);
suffstat.GA = sum(GA,3);
suffstat.Mtil = sum(Mtil,3);
suffstat.Udot = sum(Udot,3);

%% return these

fromVBEstep.mumarg = mumarg;
fromVBEstep.inv_sigmarg = inv_sigmarg;
fromVBEstep.mumarg0 = mumarg0;
fromVBEstep.inv_sigmarg0 = inv_sigmarg0;
fromVBEstep.crosscov0 = crosscov0;
fromVBEstep.crosscov = crsscov;
fromVBEstep.suffstat = suffstat;

