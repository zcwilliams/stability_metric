clear all; close all;

% (0) enter directory path containing files in 'dirloc'
dirloc = '~/Desktop/williams_mcnamara_natureReports/github/use'; % change to your folder location here
cd(dirloc)

% (1) set up Lorenz model params to generate time series
r          = 45;    % model param
s          = 30;    % model param 
b          = 8/3;   % model param
stochastic = 1;     % multiplicative noise = 1/ no noise = 0,
sigma      = 0.2;   % multiplicative noise variance
dt         = 0.005; % model time step
time       = 200;   % total time for Lorenz (points time/dt)
transient  = 1000;  % discard the first 'transient'-points of Lorenz time series output
n_lib      = 2000;  % number phase space points in library attractor
n_test     = 2000;  % number phase space points in test attractor

% phase space reconstruction - delay-embedding parameters
E          = 3;     % embedding dimension
tau        = 20;    % embedding time lag

% stability metric parameters
L          = [dt:dt:dt*500]'; % local time measured from a given point
npt        = 1000;            % total number of points used in test attractor

if stochastic==1
    [XYZ] = lorenz_stochastic(r, s, b, sigma, time, dt); % Lorenz with multiplicative noise time series
else
    [XYZ] = lorenz_RK(r, s, b, time, dt); % Lorenz-83 model time series
end


% % (2) delay-embed X_lib & X_test into E-dimensional space, using embedding lag tau
index_lib   = transient+1:transient+n_lib+(E-1)*tau;
index_test  = transient+n_lib+(E-1)*tau+1:transient+n_lib+n_test+2*(E-1)*tau;
[embd_lib]  = delay_embed(XYZ(index_lib,1),E,tau);  % reconstruct library attractor 
[embd_test] = delay_embed(XYZ(index_test,1),E,tau); % reconstruct test attractor

% % % or switch to full attractor using x,y, and vars. 
% embd_lib=XYZ(transient+1:transient+n_lib,:);
% embd_test=XYZ(transient+n_lib+1:transient+n_lib+n_test,:);

% determine nearest neighbor (in library attractor) for each point in test attractor
ns            = createns(embd_lib);
[xn,dstore]   = knnsearch(ns,embd_test,'k',1);

% (4) calculate lambda^- and lambda^+
[L,lamF,lamB] = getlambda(E,L,dt,dstore,xn,embd_lib,embd_test,npt);

% take the difference between backward and forward time lambda
lambda_diff   = lamB-lamF;
[dLam,index]  = max(lambda_diff);

% plot
plot(L,lambda_diff,'k','linewidth',2)
hold on
xlabel('L','fontsize',14)
ylabel('\lambda^--\lambda^+','fontsize',14)
plot(L(index),dLam,'rx','markersize',10,'linewidth',3)
text(L(index),dLam+(1/12)*max(lambda_diff),'\Delta \lambda','fontsize',20,'color',[1 0 0])
ylim([4/3*min(lambda_diff) 4/3*max(lambda_diff)])
