function [xyzinpn, ddotU, dotY] = generate_data_PLDS_NSFR(T, params)
%% generate data from PLDS

%% (1) unpack params

% prior on x
V0 = params.V0; x0 = params.x0;

% latent noise covariance
Q = params.Q;

A = params.A; B = params.B; C = params.C; d = params.d;

h = params.h;

inp = params.inpn;

%% extracting essential quantities

[p, k] = size(C);
pinp = size(B,2);

y = zeros(p,T);
x = zeros(k,T);
z = zeros(p,T);

sqQ = sqrtm(Q);

% generate hidden states
x(:,1)=  sqrtm(V0)*randn(k,1)+ x0+ B*inp(:,1);
z(:,1) = C*(x(:,1)+h) + d;
y(:,1) = poissrnd(exp(z(:,1)));

Ax = zeros(k, T-1);
Binp = zeros(k, T-1);

for t=2:T
    
    Ax(:,t-1) = A*x(:,t-1);
    Binp(:,t-1) = B*inp(:,t);
    
    x(:,t) = A*x(:,t-1)+ B*inp(:,t)+ sqQ*randn(k,1);
    
    z(:,t) = C*(x(:,t)+h) + d;
    y(:,t) = poissrnd(exp(z(:,t)));
end

% some of statistics that are required in M step
ddotUmat = zeros(pinp, pinp,T);
for i = 1: T
    ddotUmat(:,:,i) = inp(:,i)*inp(:,i)';
end

ddotU = sum(ddotUmat,3);
dotY = sum(y, 2);

% you want to make sure Ax and Binp have higher pwr than noise
% [cov(Ax') cov(Binp') sqQ]

% generate input_tst and x_tst for computing prediction error
inpn_tst = generate_inputs(1, pinp);
x_tst = A*x(:,end) + B*inpn_tst + sqQ*randn(k,1);
z_tst = C*(x_tst+h) + d;


%% return these

xyzinpn.x = x;
xyzinpn.y = y;
xyzinpn.z = z;
xyzinpn.inp = inp;

xyzinpn.x_tst = x_tst;
xyzinpn.z_tst = z_tst;
xyzinpn.inp_tst = inpn_tst;




