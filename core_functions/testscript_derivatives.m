% test derivatives w.r.t. h

p = 3; % # neurons
k = 2; % latent dimension
T = 10; % # of time steps

C = rand(p,k);
d = -1*ones(p,1);
h = rand(k,1);

omega = randn(k, T); 
Upsilon = repmat(0.5*eye(k), [1 1 T]);

%% (1) expression using summation (i-th trial alone)

dfdh_1 = zeros(k,p);
sumovert = zeros(p,T);

for s=1:p
    cstrp = C(s,:);
    ds = d(s);
    
    for t=1:T
        sumovert(s,t) = exp(cstrp*omega(:,t)+ 0.5*cstrp*Upsilon(:,:,t)*cstrp');
    end
    
    dfdh_1(:,s) = cstrp*exp(ds+cstrp*h)*sum(sumovert(s,:));
end
    
sum(dfdh_1,2)

%% (2) expression using matrix notation (i-th trial alone)

gC = sum(sumovert,2);

C'*diag(exp(C*h + d))*gC

%% compare now for several trials

clear all;
clc;

p = 4; % # neurons
k = 3; % latent dimension
T = 10; % # of time steps
r = 5; 

C = rand(p,k);
d = -1*ones(p,1);
h = rand(r*k,1);

omega = randn(k, T, r); 
Upsilon = repmat(0.5*eye(k), [1 1 T r]);

%% summation

dfdh_all  = zeros(k, r);
gC = zeros(p*r,1);

for i=1:r
    
    dfdh_1 = zeros(k,p);
    sumovert = zeros(p,T);
    hi = h(k*(i-1)+1:k*i);
    
    for s=1:p
        cstrp = C(s,:);
        ds = d(s);
        
        for t=1:T
            sumovert(s,t) = exp(cstrp*omega(:,t,i)+ 0.5*cstrp*Upsilon(:,:,t,i)*cstrp');
        end
        
        dfdh_1(:,s) = cstrp*exp(ds+cstrp*hi)*sum(sumovert(s,:));
    end
    
    gC(p*(i-1)+1:i*p) = sum(sumovert,2);
    
    dfdh_all(:,i) = sum(dfdh_1,2);
    
end

dfdh_all

%% matrix

Ones_r = ones(r,1);
Cr = kron(Ones_r, C);
Cbd = kron(eye(r), C);
dr = repmat(d, [r 1]);

reshape(Cbd'*diag(exp(Cbd*h+dr))*gC, k, [])



