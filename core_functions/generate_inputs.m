function inpn = generate_inputs(T, pinp)

twoperiods = 50; ti = [1:T]/twoperiods*pi;
if pinp==0
    inpn = zeros(pinp, T);
else
    inpn = randn(pinp, T);
%     inpn = [cos(ti); sin(ti); rand(pinp-2,T)-1]; 
end