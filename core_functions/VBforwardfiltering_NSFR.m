function ff = VBforwardfiltering_NSFR(Yn, initparams, opts)

% ff = VBforwardfiltering_NSFR(Yn, params, opts);

%% unpack params

A = initparams.A;
AA = initparams.AA;

B = initparams.B;
covB = initparams.covB;

AB = initparams.AB;

C = initparams.C;

h = initparams.h;
covh = initparams.covh;

d = initparams.d; 

x0 = initparams.x0;
V0 = initparams.V0;

inpn = initparams.inpn;

%% forward filtering

% matrices for storing mean/cov of forward messages
T = size(Yn,2);
k = length(x0);

sigstar = zeros(k, k, T);

mutilde = zeros(k, T);
sigtilde = zeros(k, k, T);

mufwrd = zeros(k, T);
sigfwrd = zeros(k, k, T);

mufwrd(:,1) = x0;
sigfwrd(:,:,1) = V0;

I = eye(k);

for t=1:T
    
    if t==1
        
        sigstar0 = inv(inv(V0) + AA);
        
        sigtilde(:,:,t) = inv(I - A*sigstar0*A');
        
%         mutilde(:,t) = sigtilde(:,:,t)*(A*sigstar0*(V0\x0));
        mutilde(:,t) = sigtilde(:,:,t)*(B*inpn(:,t) + A*(sigstar0/V0)*x0 - A*sigstar0*AB*inpn(:,t));
        
        [mufwrd(:,t), sigfwrd(:,:,t)] = VBfindModeCov_NSFR(mutilde(:,t), sigtilde(:,:,t), C, d, h, covh, Yn(:,t), opts);
        
    else
        
        
        %(1) compute sigstar
        sigstar(:,:,t-1) = inv(inv(sigfwrd(:,:,t-1)) + AA);
        
        %(2) compute mutilde and sigtilde
        sigtilde(:,:,t) = inv(I - A*sigstar(:,:,t-1)*A');
        
%         mutilde(:,t) = sigtilde(:,:,t)*A*sigstar(:,:,t-1)*(sigfwrd(:,:,t-1)\mufwrd(:,t-1));
        mutilde(:,t) = sigtilde(:,:,t)*(B*inpn(:,t) + A*(sigstar(:,:,t-1)/sigfwrd(:,:,t-1))*mufwrd(:,t-1) - A*sigstar(:,:,t-1)*AB*inpn(:,t));
        
        [mufwrd(:,t), sigfwrd(:,:,t)] = VBfindModeCov_NSFR(mutilde(:,t), sigtilde(:,:,t), C, d, h, covh, Yn(:,t), opts);
        
    end
    
end

%% return these

ff.mufwrd = mufwrd;
ff.sigfwrd = sigfwrd;
ff.sigstar = sigstar;
ff.sigstar0 = sigstar0;
