%% Call to the network
load('net_trained_in_time.mat')
mu0 = linspace(1e-4, 1e-2,20);
mu = normalize(mu0,'range');
mu0 = mu0(2:end-1);
mu = mu(2:end-1);

cellmu = cell(18,1);
for i = 1:18
    s = mu(i);
    s = num2cell(s,1);
    cellmu(i) = s;
end

taupred = predict(net_trained_in_time,cellmu);

%% Data used for the comparison
mu = mu0(3);
tau = double(abs(taupred{3}));

%% Without stabilization - only CG
[errors,solutions,femregion,Dati]=C_main2D('Test1',3,0,mu);
errors

%% Theoretical Estimate for tau
[errors,solutions,femregion,Dati]=C_main2D_teorico('Test1',3,mu);
errors

%% ANN Estimate for tau
[errors,solutions,femregion,Dati]=C_main2D('Test1',3,tau,mu);
errors
